# ============================================
# Integration
# HD scRNA-seq Analysis
# Following Bøstrand et al. 2024 methods
#
# Input:  24 QC-filtered RDS files (from Script 01)
# Output: merged_filtered.rds
#
# Pipeline:
# 1.  Load 24 RDS files
# 2.  Add metadata (condition + region)
# 3.  Merge all 24 samples
# 4.  JoinLayers + FindVariableFeatures + ScaleData + PCA
# 5.  ElbowPlot → choose dims
# 6.  UMAP before integration
# 7.  Harmony integration (by orig.ident)
# 8.  UMAP after integration
# 9.  Clustering — compare resolutions, choose 0.25
# 10. Post-clustering QC — remove donor-dominated clusters
# 11. Save figures and final object
# ============================================

set.seed(42)

# Load libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)

# SET YOUR WORKING DIRECTORY BEFORE RUNNING
# Change this to wherever you saved the project
# on your computer. Examples:
# Windows: setwd("C:/Users/YourName/Documents/HD_Project_local")
# Mac:     setwd("/Users/YourName/Documents/HD_Project_local")
# ============================================
setwd("YOUR/PATH/TO/HD_Project_local")

# Create output folders
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ============================================
# Step 1 - Load all 24 QC-filtered RDS files
# Files are QC-filtered and normalized by upstream Script 01
# ============================================

# CB samples
CB_HD1   <- readRDS("CB_HD1.rds")
CB_HD2   <- readRDS("CB_HD2.rds")
CB_HD3   <- readRDS("CB_HD3.rds")
CB_Ctrl1 <- readRDS("CB_Ctrl1.rds")
CB_Ctrl2 <- readRDS("CB_Ctrl2.rds")
CB_Ctrl3 <- readRDS("CB_Ctrl3.rds")

# IFG samples
IFG_HD1   <- readRDS("IFG_HD1.rds")
IFG_HD2   <- readRDS("IFG_HD2.rds")
IFG_HD3   <- readRDS("IFG_HD3.rds")
IFG_Ctrl1 <- readRDS("IFG_Ctrl1.rds")
IFG_Ctrl2 <- readRDS("IFG_Ctrl2.rds")
IFG_Ctrl3 <- readRDS("IFG_Ctrl3.rds")

# CN samples
CN_HD1   <- readRDS("CN_HD1.rds")
CN_HD2   <- readRDS("CN_HD2.rds")
CN_HD3   <- readRDS("CN_HD3.rds")
CN_ctrl1 <- readRDS("CN_ctrl1.rds")
CN_ctrl2 <- readRDS("CN_ctrl2.rds")
CN_ctrl3 <- readRDS("CN_ctrl3.rds")

# HIP samples
HIP_HD1   <- readRDS("HIP_HD1.rds")
HIP_HD2   <- readRDS("HIP_HD2.rds")
HIP_HD3   <- readRDS("HIP_HD3.rds")
HIP_Ctrl1 <- readRDS("HIP_Ctrl1.rds")
HIP_Ctrl2 <- readRDS("HIP_Ctrl2.rds")
HIP_Ctrl3 <- readRDS("HIP_Ctrl3.rds")

# ============================================
# Step 2 - Add condition and region metadata
# Done BEFORE merging so each cell carries
# the correct label into the merged object
# ============================================

# CB samples
CB_HD1$condition <- "HD";        CB_HD1$region <- "CB"
CB_HD2$condition <- "HD";        CB_HD2$region <- "CB"
CB_HD3$condition <- "HD";        CB_HD3$region <- "CB"
CB_Ctrl1$condition <- "Control"; CB_Ctrl1$region <- "CB"
CB_Ctrl2$condition <- "Control"; CB_Ctrl2$region <- "CB"
CB_Ctrl3$condition <- "Control"; CB_Ctrl3$region <- "CB"

# IFG samples
IFG_HD1$condition <- "HD";        IFG_HD1$region <- "IFG"
IFG_HD2$condition <- "HD";        IFG_HD2$region <- "IFG"
IFG_HD3$condition <- "HD";        IFG_HD3$region <- "IFG"
IFG_Ctrl1$condition <- "Control"; IFG_Ctrl1$region <- "IFG"
IFG_Ctrl2$condition <- "Control"; IFG_Ctrl2$region <- "IFG"
IFG_Ctrl3$condition <- "Control"; IFG_Ctrl3$region <- "IFG"

# CN samples
CN_HD1$condition <- "HD";        CN_HD1$region <- "CN"
CN_HD2$condition <- "HD";        CN_HD2$region <- "CN"
CN_HD3$condition <- "HD";        CN_HD3$region <- "CN"
CN_ctrl1$condition <- "Control"; CN_ctrl1$region <- "CN"
CN_ctrl2$condition <- "Control"; CN_ctrl2$region <- "CN"
CN_ctrl3$condition <- "Control"; CN_ctrl3$region <- "CN"

# HIP samples
HIP_HD1$condition <- "HD";        HIP_HD1$region <- "HIP"
HIP_HD2$condition <- "HD";        HIP_HD2$region <- "HIP"
HIP_HD3$condition <- "HD";        HIP_HD3$region <- "HIP"
HIP_Ctrl1$condition <- "Control"; HIP_Ctrl1$region <- "HIP"
HIP_Ctrl2$condition <- "Control"; HIP_Ctrl2$region <- "HIP"
HIP_Ctrl3$condition <- "Control"; HIP_Ctrl3$region <- "HIP"

# ============================================
# Step 3 - Merge all 24 samples
# ============================================

merged <- merge(
  CB_HD1,
  y = list(
    CB_HD2, CB_HD3, CB_Ctrl1, CB_Ctrl2, CB_Ctrl3,
    IFG_HD1, IFG_HD2, IFG_HD3, IFG_Ctrl1, IFG_Ctrl2, IFG_Ctrl3,
    HIP_HD1, HIP_HD2, HIP_HD3, HIP_Ctrl1, HIP_Ctrl2, HIP_Ctrl3,
    CN_HD1, CN_HD2, CN_HD3, CN_ctrl1, CN_ctrl2, CN_ctrl3
  ),
  add.cell.ids = c(
    "CB_HD1", "CB_HD2", "CB_HD3", "CB_Ctrl1", "CB_Ctrl2", "CB_Ctrl3",
    "IFG_HD1", "IFG_HD2", "IFG_HD3", "IFG_Ctrl1", "IFG_Ctrl2", "IFG_Ctrl3",
    "HIP_HD1", "HIP_HD2", "HIP_HD3", "HIP_Ctrl1", "HIP_Ctrl2", "HIP_Ctrl3",
    "CN_HD1", "CN_HD2", "CN_HD3", "CN_ctrl1", "CN_ctrl2", "CN_ctrl3"
  ),
  project = "HD_Integration"
)

# Verify merge
merged
table(merged$condition, merged$region)

# ============================================
# Step 4 - Preprocessing
# Normalization already done in upstream Script 01
# JoinLayers required for Seurat v5
# ============================================

merged <- JoinLayers(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)

# Save checkpoint after PCA — expensive step
saveRDS(merged, "merged_pca.rds")

# ============================================
# Step 5 - ElbowPlot
# dims = 1:24 chosen — curve flattens after PC20-25
# Consistent across UMAP before, UMAP after,
# FindNeighbors and FindClusters
# ============================================

ElbowPlot(merged, ndims = 50)

ggsave("figures/ElbowPlot.png",
       width = 8, height = 6, dpi = 300)

# ============================================
# Step 6 - UMAP Before Integration
# Shows batch effects before Harmony correction
# Expected: cells separate by sample/region not cell type
# ============================================

merged <- RunUMAP(merged,
                  reduction = "pca",
                  dims = 1:24,
                  reduction.name = "umap_before",
                  reduction.key = "UMAPbefore_")

# Create before integration plots
p_before_condition <- DimPlot(merged,
                              reduction = "umap_before",
                              group.by = "condition",
                              pt.size = 0.1) +
  ggtitle("Before Integration - by Condition")

p_before_region <- DimPlot(merged,
                           reduction = "umap_before",
                           group.by = "region",
                           pt.size = 0.1) +
  ggtitle("Before Integration - by Region")

p_before_sample <- DimPlot(merged,
                           reduction = "umap_before",
                           group.by = "orig.ident",
                           label = TRUE,
                           label.size = 2,
                           pt.size = 0.1) +
  NoLegend() +
  ggtitle("Before Integration - by Sample")

# Save before integration plots
ggsave("figures/UMAP_Before_by_Condition.png",
       plot = p_before_condition, width = 8, height = 6, dpi = 300)

ggsave("figures/UMAP_Before_by_Region.png",
       plot = p_before_region, width = 8, height = 6, dpi = 300)

ggsave("figures/UMAP_Before_by_Sample.png",
       plot = p_before_sample, width = 10, height = 8, dpi = 300)

# ============================================
# Step 7 - Harmony Integration
# Corrects for donor/sample batch effects
# group.by.vars = "orig.ident" — corrects for sample/donor
# NEVER use "condition" — that removes the biology you want to study
# dims = 1:24 consistent with ElbowPlot choice
# ============================================

merged <- RunHarmony(merged,
                     group.by.vars = "orig.ident",
                     plot_convergence = TRUE,
                     nclust = 50,
                     max_iter = 10,
                     early_stop = TRUE,
                     reduction = "pca",
                     reduction.save = "harmony")

ggsave("figures/Harmony_convergence.png",
       width = 8, height = 6, dpi = 300)

# Save after Harmony
saveRDS(merged, "merged_harmony.rds")

# ============================================
# Step 8 - UMAP After Integration
# Expected: cells now cluster by cell type not by sample
# ============================================

merged <- RunUMAP(merged,
                  reduction = "harmony",
                  dims = 1:24,
                  reduction.name = "umap_after",
                  reduction.key = "UMAPafter_")

# Create after integration plots
p_after_condition <- DimPlot(merged,
                             reduction = "umap_after",
                             group.by = "condition",
                             pt.size = 0.1) +
  ggtitle("After Integration - by Condition")

p_after_region <- DimPlot(merged,
                          reduction = "umap_after",
                          group.by = "region",
                          pt.size = 0.1) +
  ggtitle("After Integration - by Region")

p_after_sample <- DimPlot(merged,
                          reduction = "umap_after",
                          group.by = "orig.ident",
                          label = TRUE,
                          label.size = 2,
                          pt.size = 0.1) +
  NoLegend() +
  ggtitle("After Integration - by Sample")

# Save individual after integration plots
ggsave("figures/UMAP_After_by_Condition.png",
       plot = p_after_condition, width = 8, height = 6, dpi = 300)

ggsave("figures/UMAP_After_by_Region.png",
       plot = p_after_region, width = 8, height = 6, dpi = 300)

ggsave("figures/UMAP_After_by_Sample.png",
       plot = p_after_sample, width = 10, height = 8, dpi = 300)

# Save combined before vs after plots
ggsave("figures/UMAP_Before_After_Condition.png",
       plot = p_before_condition | p_after_condition,
       width = 16, height = 6, dpi = 300)

ggsave("figures/UMAP_Before_After_Region.png",
       plot = p_before_region | p_after_region,
       width = 16, height = 6, dpi = 300)

ggsave("figures/UMAP_Before_After_Sample.png",
       plot = p_before_sample | p_after_sample,
       width = 16, height = 6, dpi = 300)

# ============================================
# Step 9 - Clustering
# Compared resolutions 0.3, 0.5, 0.8 before choosing
# Resolution 0.25 chosen — gives biologically meaningful
# clusters consistent with Bøstrand et al. 2024
# dims = 1:24 consistent with UMAP after integration
# ============================================

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:24)

# Compare resolutions before choosing
merged <- FindClusters(merged, resolution = 0.3)
p_res03 <- DimPlot(merged, reduction = "umap_after",
                   label = TRUE) + ggtitle("Resolution 0.3")

merged <- FindClusters(merged, resolution = 0.5)
p_res05 <- DimPlot(merged, reduction = "umap_after",
                   label = TRUE) + ggtitle("Resolution 0.5")

merged <- FindClusters(merged, resolution = 0.8)
p_res08 <- DimPlot(merged, reduction = "umap_after",
                   label = TRUE) + ggtitle("Resolution 0.8")

p_res03 | p_res05 | p_res08

ggsave("figures/Resolution_comparison_0.3_0.5_0.8.png",
       width = 24, height = 8, dpi = 300)

# Apply chosen resolution 0.25
merged <- FindClusters(merged, resolution = 0.25)

DimPlot(merged,
        reduction = "umap_after",
        label = TRUE,
        pt.size = 0.1) +
  ggtitle("Clusters after Integration - Resolution 0.25")

ggsave("figures/Clusters_Resolution_0.25.png",
       width = 10, height = 8, dpi = 300)

# ============================================
# Step 10 - Post-Clustering QC
# Remove donor-dominated and small clusters
# Following Bøstrand et al. 2024 criteria:
# - One donor >50% of cluster → remove
# - Fewer than 5 donors contributing ≥2% → remove
# - Fewer than 100 cells → remove
# ============================================

# Check cluster sizes
table(merged$RNA_snn_res.0.25)

# Check donor contribution per cluster
cluster_donor <- table(merged$RNA_snn_res.0.25, merged$orig.ident)
cluster_donor_pct <- prop.table(cluster_donor, margin = 1) * 100

# Maximum donor contribution per cluster
round(apply(cluster_donor_pct, 1, max), 1)

# Count donors contributing ≥2% per cluster
donors_per_cluster <- apply(cluster_donor_pct, 1,
                            function(x) sum(x >= 2))
donors_per_cluster

# Remove cluster 17 — 62 cells, 100% from one donor (CN_HD3)
merged <- subset(merged, subset = RNA_snn_res.0.25 != "17")
merged$RNA_snn_res.0.25 <- droplevels(merged$RNA_snn_res.0.25)

# Remove cluster 9 — 99.9% from one donor (CB_Ctrl1)
# Remove cluster 11 — 92.3% from one donor (CB_Ctrl1)
# Remove cluster 12 — 50.5% from one donor + fewer than 5 donors ≥2%
merged <- subset(merged,
                 subset = RNA_snn_res.0.25 != "9" &
                   RNA_snn_res.0.25 != "11" &
                   RNA_snn_res.0.25 != "12")

merged$RNA_snn_res.0.25 <- droplevels(merged$RNA_snn_res.0.25)

# Set active identity to chosen resolution
Idents(merged) <- "RNA_snn_res.0.25"

# Verify final cell count and cluster distribution
cat("Final cell count:", ncol(merged), "\n")
table(merged$RNA_snn_res.0.25)

# Visualize final clusters after QC
DimPlot(merged,
        reduction = "umap_after",
        label = TRUE,
        pt.size = 0.1) +
  ggtitle("Final Clusters after QC - Resolution 0.25")

ggsave("figures/Clusters_final_after_QC.png",
       width = 10, height = 8, dpi = 300)

# ============================================
# Step 11 - Save final object and session info
# ============================================

saveRDS(merged, "merged_filtered.rds")

sink("results/sessionInfo_integration.txt")
sessionInfo()
sink()