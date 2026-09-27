# =============================================================
#cell_communication.R
# Cell-Cell Communication analysis: HD vs Control
#
# Question this script addresses:
#   Do ligand-receptor interactions between the HD-enriched glial
#   subpopulations change between HD and Control, and is there
#   evidence they coordinate rather than respond independently?
#
# Steps 1-6:
# (createCellChat -> DB -> subsetData -> overexpressed genes ->
#  computeCommunProb -> filterCommunication -> pathway level ->
#  aggregateNet -> pathway list -> circle plot -> LR table),
# run once per condition. Steps 7+ extend this to a direct
# HD vs Control comparison.
#
# Input : the group's finalized, LOCKED annotated object
#         (must contain a subcluster/cell-type column and a
#          "condition" column with values "HD" / "Control")
# Output: results/cellchat_HD.rds, results/cellchat_Control.rds,
#         results/cellchat_merged.rds, figures/*.png, results/*.csv
# =============================================================

set.seed(42)   # reproducibility for any stochastic steps

# ---- Step 0: Libraries (matches tutorial) ----
library(CellChat)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("single-cell-final-project")
getwd()   # confirm it's now inside single-cell-final-project
file.exists("data/merged_subclustered_final.rds")   # should now be TRUE

# ---- Step 1: Load the group's LOCKED annotated object ----
# Use a path RELATIVE to the project root, not an absolute
# machine-specific path.
merged <- readRDS("data/merged_subclustered_final.rds")
# NOTE: use the version saved AFTER the subcluster-assembly fix
# (Excluded_contaminant cells removed / properly labeled), not an
# earlier save from before that correction.

# ---- Known upstream issue: some cells retain their old broad-category
# label instead of a fine subcluster label (leftover from contamination
# removal during subclustering upstream). We filter these out here so they don't get treated as
# legitimate populations in the communication analysis. Document this
# as a known limitation in the README.
leftover_broad_labels <- c("Astrocytes", "Cerebellar Granule Cells",
                            "Excitatory Neurons", "Inhibitory Neurons",
                            "Microglia", "Oligodendrocytes",
                            "OPCs", "Vascular Cells")
# NOTE: "Peripheral Immune Cells" is NOT in this list - that population
# was intentionally never subclustered further, so its broad label is
# the correct, final label, not leftover contamination.

n_before <- ncol(merged)
merged <- subset(merged, subset = !(subcluster %in% leftover_broad_labels))
merged$subcluster <- droplevels(as.factor(merged$subcluster))
cat("Removed", n_before - ncol(merged),
    "cells with leftover broad-category labels out of", n_before, "total\n")

# Confirm the default assay is RNA (log-normalized "data" slot) -
# CellChat pulls the data slot of the DEFAULT assay when you pass
# a Seurat object directly, so this must be RNA, not an
# integration-only assay.
DefaultAssay(merged)          # should print "RNA"
# If not, uncomment:
# DefaultAssay(merged) <- "RNA"

# Sanity checks before doing anything else
stopifnot("condition" %in% colnames(merged@meta.data))
stopifnot(all(c("HD", "Control") %in% unique(merged$condition)))

# Set this to your group's actual finalized annotation column
# (e.g. "subcluster" from the subclustering script).
annotation_col <- "subcluster"
stopifnot(annotation_col %in% colnames(merged@meta.data))

table(merged@meta.data[[annotation_col]], merged$condition)

# ---- Step 2: Split into HD and Control ----
# CellChat compares networks BETWEEN two separately-run objects,
# not within one object with a condition column - so we build
# one CellChat object per condition.
seurat_HD      <- subset(merged, subset = condition == "HD")
seurat_Control <- subset(merged, subset = condition == "Control")

# ---- Step 3: Build + run the CellChat pipeline per condition 

run_cellchat_pipeline <- function(seurat_obj, group_col) {

  # 1. Initialization & Database Setup
  cellchat <- createCellChat(object = seurat_obj, group.by = group_col)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)

  # 2. Identify Overexpressed Genes & Interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 3. Infer Cell-Cell Communication 
  cellchat <- computeCommunProb(cellchat)   # default raw.use, as taught
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # 4. Build the Communication Network 
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}

cellchat_HD      <- run_cellchat_pipeline(seurat_HD, annotation_col)
cellchat_Control <- run_cellchat_pipeline(seurat_Control, annotation_col)

# ---- Step 4: Save intermediate objects ----
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

saveRDS(cellchat_HD,      "results/cellchat_HD.rds")
saveRDS(cellchat_Control, "results/cellchat_Control.rds")

# ---- Step 5: Circle plots per condition 
png("figures/circle_overall_Control.png", width = 700, height = 600)
netVisual_circle(cellchat_Control@net$weight,
                  vertex.weight = as.numeric(table(cellchat_Control@idents)),
                  weight.scale = TRUE,
                  label.edge = FALSE,
                  title.name = "Control: Overall Communication Network")
dev.off()

png("figures/circle_overall_HD.png", width = 700, height = 600)
netVisual_circle(cellchat_HD@net$weight,
                  vertex.weight = as.numeric(table(cellchat_HD@idents)),
                  weight.scale = TRUE,
                  label.edge = FALSE,
                  title.name = "HD: Overall Communication Network")
dev.off()

# ---- Step 6: Extract + view L-R results per condition 
lr_HD      <- subsetCommunication(cellchat_HD)
lr_Control <- subsetCommunication(cellchat_Control)

head(lr_HD)
head(lr_Control)

write.csv(lr_HD,      "results/interactions_HD.csv", row.names = FALSE)
write.csv(lr_Control, "results/interactions_Control.csv", row.names = FALSE)

# ---- Step 7: List + inspect active signaling pathways-
# Look through these lists for pathways relevant to HD pathology -
# e.g. inflammatory (TNF, IL1, CSF), complement, or stress-related
# pathways - before deciding which ones to visualize individually.
cellchat_HD@netP$pathways
cellchat_Control@netP$pathways



setdiff(cellchat_HD@netP$pathways,      cellchat_Control@netP$pathways)   # only in HD
setdiff(cellchat_Control@netP$pathways, cellchat_HD@netP$pathways)       # only in Control
netVisual_aggregate(cellchat_Control, signaling = "CXCL", layout = "circle")
# (no HD version needed for CXCL - it's absent, that absence IS the finding)

# Edit "PATHWAY_NAME" to one that actually appears above
# and is relevant to your question, then repeat for each condition):
netVisual_aggregate(cellchat_HD,      signaling = "COMPLEMENT", layout = "circle")
netVisual_aggregate(cellchat_Control, signaling = "COMPLEMENT", layout = "circle")
png("figures/complement_HD.png", width = 1400, height = 1200, res = 150)
netVisual_aggregate(cellchat_HD, signaling = "COMPLEMENT", layout = "circle",
                    vertex.label.cex = 0.6)
dev.off()

png("figures/complement_Control.png", width = 1400, height = 1200, res = 150)
netVisual_aggregate(cellchat_Control, signaling = "COMPLEMENT", layout = "circle",
                    vertex.label.cex = 0.6)
dev.off()

# =============================================================

# ---- Step 8: Merge for direct HD vs Control comparison ----
object.list <- list(Control = cellchat_Control, HD = cellchat_HD)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat_merged, "results/cellchat_merged.rds")

# Overall information flow comparison across all pathways
p_rankNet <- rankNet(cellchat_merged, mode = "comparison",
                      stacked = TRUE, do.stat = TRUE)
ggsave("figures/rankNet_pathways.png", p_rankNet, width = 8, height = 10)

# Differential number/strength of interactions (HD vs Control)
png("figures/diff_interactions.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
dev.off()

# ---- Step 9: Focus on the HD-enriched / stress-response populations ----


populations_of_interest <- c("Astro_Thyme", "Oligo_Oak", "Mglia_Violet", "Mglia_Daisy")

#Bubble plot: all significant L-R pairs between these populations,
# split by condition, for direct comparison
p_bubble <- netVisual_bubble(
  cellchat_merged,
  sources.use = populations_of_interest,
  targets.use = populations_of_interest,
  comparison = c(1, 2),        # Control = 1, HD = 2
  angle.x = 45)


ggsave("figures/bubble_enriched_populations.png", p_bubble,
      width = 10, height = 8)

# ---- Trimmed bubble plot for presentation ----
# The full bubble plot (all significant L-R pairs among our 4
# populations of interest) has too many rows to read on a slide.
# This keeps only the top 25 strongest interactions (by
# communication probability, pooled across HD and Control) so the
# figure is legible for the presentation. The complete, untrimmed
# table is still saved separately (interactions_HD/Control_focus.csv)
# for the repo and for anyone who wants the full picture.
top_pairs <- dplyr::bind_rows(lr_HD_focus, lr_Control_focus) %>%
  dplyr::distinct(interaction_name, .keep_all = TRUE) %>%
  dplyr::arrange(dplyr::desc(prob)) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(interaction_name)

p_bubble_top <- netVisual_bubble(
  cellchat_merged,
  sources.use = populations_of_interest,
  targets.use = populations_of_interest,
  comparison = c(1, 2),
  angle.x = 45,
  pairLR.use = data.frame(interaction_name = top_pairs)
)

ggsave("figures/bubble_top25.png", p_bubble_top, width = 10, height = 10)

# ---- Step 10: Export interaction tables restricted to populations of interest ----
lr_HD_focus <- lr_HD %>%
 filter(source %in% populations_of_interest,
        target %in% populations_of_interest)
lr_Control_focus <- lr_Control %>%
  filter(source %in% populations_of_interest,
         target %in% populations_of_interest)

write.csv(lr_HD_focus,
          "results/interactions_HD_enriched_populations.csv", row.names = FALSE)
write.csv(lr_Control_focus,
          "results/interactions_Control_enriched_populations.csv", row.names = FALSE)



# ---- Step 11: Session info for reproducibility ----
writeLines(capture.output(sessionInfo()), "results/sessionInfo_cellchat.txt")



# =============================================================
# =============================================================
# See README.md for full interpretation notes, limitations,
# and the final write-up of results and conclusions.
# =============================================================
# =============================================================
