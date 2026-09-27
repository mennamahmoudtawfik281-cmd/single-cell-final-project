# ==========================================
# 1. Load Required Libraries
# ==========================================

library(Seurat)
library(tidyverse)
# ==========================================
# 2. Load Final Object
# ==========================================
obj <- readRDS(file.choose())

obj

colnames(obj@meta.data)
head(obj@meta.data)
table(obj$condition)
table(obj$region)
table(obj$orig.ident)
table(obj$broad_cell_type)
table(obj$condition, obj$region)
# ==========================================
# 3. Check cell counts per donor x cell type x region
# ==========================================

cell_counts <- obj@meta.data %>%
  group_by(region, orig.ident, condition, broad_cell_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(region, broad_cell_type, orig.ident)

print(cell_counts, n = 200)

low_count_groups <- cell_counts %>%
  filter(n_cells < 10)

print(low_count_groups, n = 100)
# ==========================================
# 4. Load DESeq2 
# ==========================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
# ==========================================
# 5. Subset object to one region (CB) pilot
# ==========================================
obj_CB <- subset(obj, subset = region == "CB")

table(obj_CB$orig.ident, obj_CB$condition)
table(obj_CB$broad_cell_type)
# ==========================================
# STEP 6: Identify valid groups (donor x cell_type) with n_cells >= 10
# ==========================================
valid_groups_CB <- cell_counts %>%
  filter(region == "CB", n_cells >= 10) %>%
  mutate(group_id = paste(orig.ident, broad_cell_type, sep = "_"))

# Check how many groups passed the threshold
nrow(valid_groups_CB)
valid_groups_CB
# ==========================================
# STEP 7: Pseudobulk aggregation using Seurat
# ==========================================
pseudobulk_CB <- AggregateExpression(
  obj_CB,
  group.by = c("orig.ident", "broad_cell_type"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE)

# Check the resulting matrix shape and sample names
dim(pseudobulk_CB$RNA)
colnames(pseudobulk_CB$RNA)[1:10]
# ==========================================
# STEP 8: Build metadata matching pseudobulk columns
# ==========================================

# Get the pseudobulk sample names (column names)
pb_colnames <- colnames(pseudobulk_CB$RNA)

# Reconstruct a matching data frame from valid_groups_CB
# (fixing underscores to dashes in orig.ident, same way Seurat did)
valid_groups_CB <- valid_groups_CB %>%
  mutate(
    orig_ident_fixed = gsub("_", "-", orig.ident),
    pb_colname = paste(orig_ident_fixed, broad_cell_type, sep = "_"))

# Check that these constructed names actually match the real column names
sum(valid_groups_CB$pb_colname %in% pb_colnames)
nrow(valid_groups_CB)

# ==========================================
# STEP 9: Filter pseudobulk matrix to valid columns only
# ==========================================

# Keep only the valid (>=10 cells) pseudobulk samples
pb_counts_CB <- pseudobulk_CB$RNA[, valid_groups_CB$pb_colname]

# Check dimensions after filtering
dim(pb_counts_CB)
# ==========================================
# STEP 10: Build sample metadata for DESeq2
# ==========================================

# Build metadata in the SAME order as the columns of pb_counts_CB
sample_metadata_CB <- valid_groups_CB %>%
  select(pb_colname, orig.ident, condition, broad_cell_type, n_cells) %>%
  column_to_rownames("pb_colname")

# Make sure the row order matches the column order of the count matrix EXACTLY
sample_metadata_CB <- sample_metadata_CB[colnames(pb_counts_CB), ]

# Double check alignment
identical(rownames(sample_metadata_CB), colnames(pb_counts_CB))

# Make condition a factor with Control as the reference level
sample_metadata_CB$condition <- factor(sample_metadata_CB$condition, levels = c("Control", "HD"))

head(sample_metadata_CB)
# ==========================================
# STEP 11: Run DESeq2 separately for each cell type in CB
# ==========================================

# Get list of cell types present in CB (with enough samples)
cell_types_CB <- unique(sample_metadata_CB$broad_cell_type)
cell_types_CB

# Empty list to store DESeq2 results per cell type
deseq_results_CB <- list()
# ==========================================
# STEP 12: Check sample counts per condition for each cell type
# ==========================================

sample_check_CB <- sample_metadata_CB %>%
  group_by(broad_cell_type, condition) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n_samples, values_fill = 0)

sample_check_CB
# ==========================================
# STEP 13: Exclude cell types with insufficient replicates
# ==========================================

# Cell types with at least 2 samples in each condition
valid_cell_types_CB <- sample_check_CB %>%
  filter(Control >= 2, HD >= 2) %>%
  pull(broad_cell_type)

valid_cell_types_CB
length(valid_cell_types_CB)
# ==========================================
# STEP 14: Run DESeq2 for each valid cell type (loop)
# ==========================================

for (ct in valid_cell_types_CB) {
  
  message("Running DESeq2 for: ", ct)
  
  # Subset metadata to this cell type only
  meta_ct <- sample_metadata_CB %>% filter(broad_cell_type == ct)
  
  # Subset counts matrix to matching samples (same order)
  counts_ct <- pb_counts_CB[, rownames(meta_ct)]
  
  # Build DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_ct,
    colData = meta_ct,
    design = ~ condition)
  
  # Pre-filter low-count genes (keep genes with at least 10 reads total across samples)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results: HD vs Control
  res <- results(dds, contrast = c("condition", "HD", "Control"))
  
  # Store as a clean data frame, sorted by adjusted p-value
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  
  # Save into the results list
  deseq_results_CB[[ct]] <- res_df}

# Confirm all cell types ran
names(deseq_results_CB)
# ==========================================
# STEP 15: Summarize DEG counts per cell type
# ==========================================

deg_summary_CB <- map_dfr(names(deseq_results_CB), function(ct) {
  df <- deseq_results_CB[[ct]]
  
  data.frame(
    cell_type = ct,
    total_genes_tested = nrow(df),
    total_DEGs = sum(df$padj < 0.05, na.rm = TRUE),
    upregulated = sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE))})

deg_summary_CB
# ==========================================
# STEP 16: Save DEG tables to CSV files
# ==========================================

dir.create("results/DEG_tables/CB", recursive = TRUE, showWarnings = FALSE)

for (ct in names(deseq_results_CB)) {
  file_name <- paste0("results/DEG_tables/CB/", gsub(" ", "_", ct), "_HD_vs_Control.csv")
  write.csv(deseq_results_CB[[ct]], file_name, row.names = FALSE)
}

list.files("results/DEG_tables/CB")
# ==========================================
# STEP 17: Volcano plot for Cerebellar Granule Cells
# ==========================================

library(ggplot2)

# Pick the cell type to plot
ct_to_plot <- "Cerebellar Granule Cells"
res_plot <- deseq_results_CB[[ct_to_plot]]

# Remove genes with NA padj (DESeq2 assigns NA to genes it excluded from testing)
res_plot <- res_plot %>% filter(!is.na(padj))

# Add a significance category for coloring
res_plot <- res_plot %>%
  mutate(sig_status = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Up in HD",
    padj < 0.05 & log2FoldChange < 0 ~ "Down in HD",
    TRUE ~ "Not significant"
  ))

# Volcano plot
volcano_plot <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig_status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up in HD" = "red", "Down in HD" = "blue", "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = paste0("Volcano Plot: ", ct_to_plot, " (HD vs Control) - CB region"),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Status"
  ) +
  theme_minimal()

volcano_plot
# ==========================================
# STEP 18: Save volcano plot to file
# ==========================================

dir.create("results/figures/CB", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "results/figures/CB/Volcano_Cerebellar_Granule_Cells_HD_vs_Control.png",
  plot = volcano_plot,
  width = 8, height = 6, dpi = 300
)

list.files("results/figures/CB")
# ==========================================
# STEP 19: Install and load pheatmap (if needed)
# ==========================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)
# ==========================================
# STEP 20: Identify cell types that have at least 1 significant DEG
# ==========================================

cell_types_with_DEGs <- deg_summary_CB %>%
  filter(total_DEGs > 0) %>%
  pull(cell_type)

cell_types_with_DEGs
# ==========================================
# STEP 21: Generate heatmap for each cell type with DEGs
# ==========================================

dir.create("results/figures/CB", recursive = TRUE, showWarnings = FALSE)

for (ct in cell_types_with_DEGs) {
  
  message("Building heatmap for: ", ct)
  
  # Re-subset counts and metadata for this cell type
  meta_ct <- sample_metadata_CB %>% filter(broad_cell_type == ct)
  counts_ct <- pb_counts_CB[, rownames(meta_ct)]
  
  # Rebuild DESeq2 object and apply VST normalization
  dds_ct <- DESeqDataSetFromMatrix(
    countData = counts_ct,
    colData = meta_ct,
    design = ~ condition)
  dds_ct <- dds_ct[rowSums(counts(dds_ct)) >= 10, ]
  vst_ct <- vst(dds_ct, blind = FALSE)
  vst_mat <- assay(vst_ct)
  
  # Get top significant DEGs for this cell type (up to 30, sorted by padj)
  top_genes <- deseq_results_CB[[ct]] %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(gene)
  
  # Subset the VST matrix to those genes
  heatmap_mat <- vst_mat[top_genes, , drop = FALSE]
  
  # Scale each gene (row) for better visual contrast (z-score across samples)
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
  
  # Build column annotation (condition)
  annotation_col <- data.frame(Condition = meta_ct$condition)
  rownames(annotation_col) <- rownames(meta_ct)
  
  # Draw and save the heatmap
  file_name <- paste0("results/figures/CB/Heatmap_", gsub(" ", "_", ct), "_HD_vs_Control.png")
  
  pheatmap(
    heatmap_mat_scaled,
    annotation_col = annotation_col,
    main = paste0("Top DEGs: ", ct, " (HD vs Control) - CB region"),
    fontsize_row = 7,
    filename = file_name,
    width = 8, height = 8)}

list.files("results/figures/CB")
# ==========================================
# STEP 22: Build a reusable function for the full pipeline
# ==========================================
run_pseudobulk_DE_pipeline <- function(obj_full, region_name, cell_counts_df, min_cells = 10, min_replicates = 2) {
  
message("=== Starting pipeline for region: ", region_name, " ===")
                          
obj_region <- subset(obj_full, subset = region == region_name)
                          
valid_groups <- cell_counts_df %>%
filter(region == region_name, n_cells >= min_cells) %>%
mutate(
orig_ident_fixed = gsub("_", "-", orig.ident),
pb_colname = paste(orig_ident_fixed, broad_cell_type, sep = "_"))

pseudobulk <- AggregateExpression(
obj_region,
group.by = c("orig.ident", "broad_cell_type"),
assays = "RNA",
slot = "counts",
return.seurat = FALSE)
                          
pb_counts <- pseudobulk$RNA[, valid_groups$pb_colname]
                          
sample_metadata <- valid_groups %>%
select(pb_colname, orig.ident, condition, broad_cell_type, n_cells) %>%
column_to_rownames("pb_colname")
sample_metadata <- sample_metadata[colnames(pb_counts), ]
sample_metadata$condition <- factor(sample_metadata$condition, levels = c("Control", "HD"))
                          
stopifnot(identical(rownames(sample_metadata), colnames(pb_counts)))
                          
sample_check <- sample_metadata %>%
group_by(broad_cell_type, condition) %>%
summarise(n_samples = n(), .groups = "drop") %>%
pivot_wider(names_from = condition, values_from = n_samples, values_fill = 0)
                          
valid_cell_types <- sample_check %>%
filter(Control >= min_replicates, HD >= min_replicates) %>%
pull(broad_cell_type)
                          
message("Valid cell types for ", region_name, ": ", paste(valid_cell_types, collapse = ", "))
                          
deseq_results <- list()
                          
for (ct in valid_cell_types) {
message("Running DESeq2 for: ", ct)
                            
meta_ct <- sample_metadata %>% filter(broad_cell_type == ct)
counts_ct <- pb_counts[, rownames(meta_ct)]
                            
dds <- DESeqDataSetFromMatrix(
countData = counts_ct,
colData = meta_ct,
design = ~ condition
                            )
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
                            
res <- results(dds, contrast = c("condition", "HD", "Control"))
res_df <- as.data.frame(res) %>%
rownames_to_column("gene") %>%
arrange(padj)
                            
deseq_results[[ct]] <- res_df
                          }
                          
deg_summary <- map_dfr(names(deseq_results), function(ct) {
df <- deseq_results[[ct]]
data.frame(
cell_type = ct,
total_genes_tested = nrow(df),
total_DEGs = sum(df$padj < 0.05, na.rm = TRUE),
upregulated = sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE),
downregulated = sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)
                            )
                          })
                          
deg_dir <- paste0("results/DEG_tables/", region_name)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
                          
for (ct in names(deseq_results)) {
file_name <- paste0(deg_dir, "/", gsub(" ", "_", ct), "_HD_vs_Control.csv")
write.csv(deseq_results[[ct]], file_name, row.names = FALSE)
}
                          
fig_dir <- paste0("results/figures/", region_name)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
                          
cell_types_with_DEGs <- deg_summary %>% filter(total_DEGs > 0) %>% pull(cell_type)
                          
for (ct in cell_types_with_DEGs) {
res_plot <- deseq_results[[ct]] %>%
filter(!is.na(padj)) %>%
mutate(sig_status = case_when(
padj < 0.05 & log2FoldChange > 0 ~ "Up in HD",
padj < 0.05 & log2FoldChange < 0 ~ "Down in HD",
TRUE ~ "Not significant"
))
                            
volcano_plot <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig_status)) +
geom_point(alpha = 0.6, size = 1.5) +
scale_color_manual(values = c("Up in HD" = "red", "Down in HD" = "blue", "Not significant" = "grey70")) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
labs(
title = paste0("Volcano Plot: ", ct, " (HD vs Control) - ", region_name, " region"),
x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Status") +
theme_minimal()
                            
ggsave(
filename = paste0(fig_dir, "/Volcano_", gsub(" ", "_", ct), "_HD_vs_Control.png"),
plot = volcano_plot, width = 8, height = 6, dpi = 300)
                            
dds_ct <- DESeqDataSetFromMatrix(
countData = pb_counts[, rownames(sample_metadata %>% filter(broad_cell_type == ct))],
colData = sample_metadata %>% filter(broad_cell_type == ct),
design = ~ condition)
dds_ct <- dds_ct[rowSums(counts(dds_ct)) >= 10, ]
vst_ct <- vst(dds_ct, blind = FALSE)
vst_mat <- assay(vst_ct)
                            
top_genes <- deseq_results[[ct]] %>%
filter(padj < 0.05) %>%
arrange(padj) %>%
slice_head(n = 30) %>%
pull(gene)
                            
heatmap_mat <- vst_mat[top_genes, , drop = FALSE]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
                            
meta_ct <- sample_metadata %>% filter(broad_cell_type == ct)
annotation_col <- data.frame(Condition = meta_ct$condition)
rownames(annotation_col) <- rownames(meta_ct)
                            
if (nrow(heatmap_mat_scaled) >= 2) {
  pheatmap(
    heatmap_mat_scaled,
    annotation_col = annotation_col,
    main = paste0("Top DEGs: ", ct, " (HD vs Control) - ", region_name, " region"),
    fontsize_row = 7,
    filename = paste0(fig_dir, "/Heatmap_", gsub(" ", "_", ct), "_HD_vs_Control.png"),
    width = 8, height = 8)} else {
  message("Skipping heatmap for ", ct, " - not enough significant genes (n < 2)")}}
                          
message("=== Finished pipeline for region: ", region_name, " ===")
                          
list(
pb_counts = pb_counts,
sample_metadata = sample_metadata,
deseq_results = deseq_results,
deg_summary = deg_summary)}  
###############################################                        
 exists("run_pseudobulk_DE_pipeline")
# ==========================================

# ==========================================

results_HIP <- run_pseudobulk_DE_pipeline(
  obj_full = obj,
  region_name = "HIP",
  cell_counts_df = cell_counts)

results_HIP$deg_summary
#=========================================
results_IFG <- run_pseudobulk_DE_pipeline(
  obj_full = obj,
  region_name = "IFG",
  cell_counts_df = cell_counts)

results_IFG$deg_summary
# ==========================================

# Create a filtered version of cell_counts for CN that excludes the
# unreliable donor CN_HD3 (has conflicting Control/HD labels - see team note)
cell_counts_CN_clean <- cell_counts %>%
  filter(!(region == "CN" & orig.ident == "CN_HD3_scRNA"))

# Confirm what's left for CN
cell_counts_CN_clean %>% filter(region == "CN") %>% distinct(orig.ident, condition)
# ==========================================
#  Run the pipeline on CN region (excluding CN_HD3)
# ==========================================

results_CN <- run_pseudobulk_DE_pipeline(
  obj_full = obj,
  region_name = "CN",
  cell_counts_df = cell_counts_CN_clean)

results_CN$deg_summary
# ==========================================
# Combine DEG summaries across all 4 regions
# ==========================================

deg_summary_CB$region <- "CB"
results_HIP$deg_summary$region <- "HIP"
results_IFG$deg_summary$region <- "IFG"
results_CN$deg_summary$region <- "CN"

deg_summary_all_regions <- bind_rows(
  deg_summary_CB,
  results_HIP$deg_summary,
  results_IFG$deg_summary,
  results_CN$deg_summary
) %>%
  select(region, cell_type, total_genes_tested, total_DEGs, upregulated, downregulated) %>%
  arrange(region, desc(total_DEGs))

# View the full summary
deg_summary_all_regions

# Save it as a master summary table
write.csv(deg_summary_all_regions, "results/DEG_tables/Summary_all_regions.csv", row.names = FALSE)
# ==========================================
#  Prepare upregulated/downregulated gene lists for Pathway analysis
# ==========================================

dir.create("results/gene_lists_for_pathway", recursive = TRUE, showWarnings = FALSE)

# Collect all DESeq2 result tables from all regions into one named list
all_deseq_results <- list(
  CB  = deseq_results_CB,
  HIP = results_HIP$deseq_results,
  IFG = results_IFG$deseq_results,
  CN  = results_CN$deseq_results)

# Loop through every region and cell type, extract up/down gene lists
for (region in names(all_deseq_results)) {
  for (ct in names(all_deseq_results[[region]])) {
    
    df <- all_deseq_results[[region]][[ct]]
    
    up_genes <- df %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene)
    down_genes <- df %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene)
    
    # Only save files if there's at least one gene (avoid empty files)
    ct_clean <- gsub(" ", "_", ct)
    
    if (length(up_genes) > 0) {
      writeLines(up_genes, paste0("results/gene_lists_for_pathway/", region, "_", ct_clean, "_UP.txt"))
    }
    if (length(down_genes) > 0) {
      writeLines(down_genes, paste0("results/gene_lists_for_pathway/", region, "_", ct_clean, "_DOWN.txt"))
    }}}

# Confirm what got created
list.files("results/gene_lists_for_pathway")
