#Sat Sep. 12, after updating R last time and repeating every step:

# Reproducibility
set.seed(42)

  # Step 1 - Load libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

# Step 2 - Set working directory
setwd("C:/Users/eman5/OneDrive/Documents/HD_Project")

# Step 3 - Load all 24 processed RDS files from upstream analysis
# Each file contains QC-filtered, normalized cells from one brain region sample

# CB samples
CB_HD1   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_HD1.rds")
CB_HD2   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_HD2.rds")
CB_HD3   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_HD3.rds")
CB_Ctrl1 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_Ctrl1.rds")
CB_Ctrl2 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_Ctrl2.rds")
CB_Ctrl3 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CB_Ctrl3.rds")

# IFG samples
IFG_HD1   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_HD1.rds")
IFG_HD2   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_HD2.rds")
IFG_HD3   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_HD3.rds")
IFG_Ctrl1 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_Ctrl1.rds")
IFG_Ctrl2 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_Ctrl2.rds")
IFG_Ctrl3 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/IFG_Ctrl3.rds")

# CN samples
CN_HD1   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_HD1.rds")
CN_HD2   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_HD2.rds")
CN_HD3   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_HD3.rds")
CN_ctrl1 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_ctrl1.rds")
CN_ctrl2 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_ctrl2.rds")
CN_ctrl3 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/CN_ctrl3.rds")

# HIP samples
HIP_HD1   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_HD1.rds")
HIP_HD2   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_HD2.rds")
HIP_HD3   <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_HD3.rds")
HIP_Ctrl1 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_Ctrl1.rds")
HIP_Ctrl2 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_Ctrl2.rds")
HIP_Ctrl3 <- readRDS("C:/Users/eman5/OneDrive/Documents/HD_Project/HIP_Ctrl3.rds")

# Step 4: Add condition and region metadata to each sample

# CB samples
CB_HD1$condition <- "HD";    CB_HD1$region <- "CB"
CB_HD2$condition <- "HD";    CB_HD2$region <- "CB"
CB_HD3$condition <- "HD";    CB_HD3$region <- "CB"
CB_Ctrl1$condition <- "Control"; CB_Ctrl1$region <- "CB"
CB_Ctrl2$condition <- "Control"; CB_Ctrl2$region <- "CB"
CB_Ctrl3$condition <- "Control"; CB_Ctrl3$region <- "CB"

# IFG samples
IFG_HD1$condition <- "HD";    IFG_HD1$region <- "IFG"
IFG_HD2$condition <- "HD";    IFG_HD2$region <- "IFG"
IFG_HD3$condition <- "HD";    IFG_HD3$region <- "IFG"
IFG_Ctrl1$condition <- "Control"; IFG_Ctrl1$region <- "IFG"
IFG_Ctrl2$condition <- "Control"; IFG_Ctrl2$region <- "IFG"
IFG_Ctrl3$condition <- "Control"; IFG_Ctrl3$region <- "IFG"

# HIP samples
HIP_HD1$condition <- "HD";    HIP_HD1$region <- "HIP"
HIP_HD2$condition <- "HD";    HIP_HD2$region <- "HIP"
HIP_HD3$condition <- "HD";    HIP_HD3$region <- "HIP"
HIP_Ctrl1$condition <- "Control"; HIP_Ctrl1$region <- "HIP"
HIP_Ctrl2$condition <- "Control"; HIP_Ctrl2$region <- "HIP"
HIP_Ctrl3$condition <- "Control"; HIP_Ctrl3$region <- "HIP"

# CN samples
CN_HD1$condition <- "HD";    CN_HD1$region <- "CN"
CN_HD2$condition <- "HD";    CN_HD2$region <- "CN"
CN_HD3$condition <- "HD";    CN_HD3$region <- "CN"
CN_ctrl1$condition <- "Control"; CN_ctrl1$region <- "CN"
CN_ctrl2$condition <- "Control"; CN_ctrl2$region <- "CN"
CN_ctrl3$condition <- "Control"; CN_ctrl3$region <- "CN"

# Step 5: Merge
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

# Step 6: Verify merge
merged
table(merged$condition, merged$region)

# Step 7: Join layers
merged <- JoinLayers(merged)

# Step 8: Variable features
merged <- FindVariableFeatures(merged, nfeatures = 2000)

# Step 9: Check layers after joining
merged

# Step 10: Scale data
merged <- ScaleData(merged)

# Step 11: Principal Component Analysis (PCA)
merged <- RunPCA(merged, npcs = 50)

# Save to LOCAL folder immediately
saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged_pca.rds")

# Check dimensions of PCA
dim(merged@reductions$pca)

# Step 12: Elbow Plot
ElbowPlot(
  merged,
  ndims = 30
)

merged <- merged %>%
  RunUMAP(
    reduction = "pca",
    dims = 1:24,
    reduction.name = "umap_before"
  )

# Step 13: UMAP before integration - colored by condition
p_before <- DimPlot(
  merged,
  reduction = "umap_before",
  group.by = "condition",
  pt.size = 0.1
) +
  ggtitle("Before Integration - by Condition")

# Show the plot
p_before

# Step 14: Complete before-integration visualization (by region and sample)
# Colored by region
DimPlot(merged, reduction = "umap_before",
        group.by = "region",
        pt.size = 0.1) +
  ggtitle("Before Integration - by Region")

# Colored by sample
DimPlot(merged, reduction = "umap_before",
        group.by = "orig.ident",
        pt.size = 0.1) +
  ggtitle("Before Integration - by Sample")

DimPlot(merged, reduction = "umap_before",
        group.by = "orig.ident",
        pt.size = 0.1) +
  ggtitle("Before Integration - by Sample") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# Step 15: Harmony Integration
merged <- RunHarmony(
  merged,
  group.by.vars = "orig.ident",
  plot_convergence = TRUE,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  reduction = "pca",
  reduction.save = "harmony"
)

# Step 16 - UMAP after integration
merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:24,
  reduction.name = "umap_after",
  reduction.key = "UMAPafter_"
)

#Save
saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged_harmony.rds")
# Step 17 - Clustering
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:24)
merged <- FindClusters(merged, resolution = 0.25)
DimPlot(merged, reduction = "umap_after", label = TRUE) + 
  ggtitle("Resolution 0.25")


# Step 18 - Visualize
DimPlot(merged, reduction = "umap_after",
        group.by = "condition",
        pt.size = 0.1) +
  ggtitle("UMAP After Integration - by Condition")

DimPlot(merged, reduction = "umap_after",
        group.by = "region",
        pt.size = 0.1) +
  ggtitle("UMAP After Integration - by Region")

# By sample - after integration
DimPlot(merged,
        reduction = "umap_after",
        group.by = "orig.ident",
        pt.size = 0.1,
        label = TRUE,
        label.size = 2) +
  NoLegend() +
  ggtitle("UMAP After Integration - by Sample")

#Before vs After Integration
# By condition
p1 <- DimPlot(merged,
              reduction = "umap_before",
              group.by = "condition",
              pt.size = 0.1) +
  ggtitle("Before Integration - by Condition")

p2 <- DimPlot(merged,
              reduction = "umap_after",
              group.by = "condition",
              pt.size = 0.1) +
  ggtitle("After Integration - by Condition")

# Combine side by side
p1 | p2

saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged_final.rds")

list.files("C:/Users/eman5/Documents/HD_Project_local/")

setwd("C:/Users/eman5/Documents/HD_Project_local")

merged
Reductions(merged)

merged <- readRDS("merged_final.rds")


# By region
p3 <- DimPlot(merged,
              reduction = "umap_before",
              group.by = "region",
              pt.size = 0.1) +
  ggtitle("Before Integration - by Region")

p4 <- DimPlot(merged,
              reduction = "umap_after",
              group.by = "region",
              pt.size = 0.1) +
  ggtitle("After Integration - by Region")

# Combine side by side
p3 | p4

p5 <- DimPlot(merged,
              reduction = "umap_before",
              group.by = "orig.ident",
              pt.size = 0.1,
              label = TRUE,
              label.size = 2) +
  NoLegend() +
  ggtitle("Before Integration - by Sample")

p6 <- DimPlot(merged,
              reduction = "umap_after",
              group.by = "orig.ident",
              pt.size = 0.1,
              label = TRUE,
              label.size = 2) +
  NoLegend() +
  ggtitle("After Integration - by Sample")

# Combine side by side
p5 | p6

#Visualize Clusters
DimPlot(merged, reduction = "umap_after",
        label = TRUE,
        pt.size = 0.1) +
  ggtitle("Clusters after Integration")


# See all resolution columns stored
colnames(merged@meta.data)
# Step 19 - Save
saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged.rds")

#Check 1 — Clusters with fewer than 100 cells:

table(merged$RNA_snn_res.0.25)

# Check donor contribution for cluster 17
cluster_donor <- table(merged$RNA_snn_res.0.25, merged$orig.ident)
cluster_donor_pct <- prop.table(cluster_donor, margin = 1) * 100
round(cluster_donor_pct["17",], 1)

# Remove cluster 17
merged <- subset(merged,
                 subset = RNA_snn_res.0.25 != "17")

# Verify it's gone
table(merged$RNA_snn_res.0.25)

# Drop empty factor levels
merged$RNA_snn_res.0.25 <- droplevels(merged$RNA_snn_res.0.25)

# Verify cluster 17 is completely gone
table(merged$RNA_snn_res.0.25)

# Update seurat clusters and identity
Idents(merged) <- "RNA_snn_res.0.25"

# Check total cell count
ncol(merged)

# Save
saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged_filtered.rds")
# Check 1 - cluster sizes
table(merged$seurat_clusters)

# Check 2 - donor contribution per cluster
cluster_donor <- table(merged$seurat_clusters, merged$orig.ident)
cluster_donor_pct <- prop.table(cluster_donor, margin = 1) * 100
round(cluster_donor_pct, 1)

# Count donors contributing ≥2% for each cluster
donors_per_cluster <- apply(cluster_donor_pct, 1, 
                            function(x) sum(x >= 2))
donors_per_cluster

# Remove only clearly donor-dominated clusters
merged <- subset(merged,
                 subset = RNA_snn_res.0.25 != "9" &
                   RNA_snn_res.0.25 != "11" &
                   RNA_snn_res.0.25 != "12")

# Clean up factor levels
merged$RNA_snn_res.0.25 <- droplevels(merged$RNA_snn_res.0.25)

# Verify
table(merged$RNA_snn_res.0.25)
ncol(merged)

# Save
saveRDS(merged, "C:/Users/eman5/Documents/HD_Project_local/merged_filtered.rds")

# Cluster 1 has only 5 donors contributing ≥2%
# and CB_HD3 contributes 38.2% — borderline
round(cluster_donor_pct["1",], 1)

# Visualize final clusters
DimPlot(merged,
        reduction = "umap_after",
        label = TRUE,
        pt.size = 0.1) +
  ggtitle("Final Clusters after QC - Resolution 0.25")

# Create results folder
dir.create("results")

# Save session Info to results folder
sink("C:/Users/eman5/Documents/HD_Project_local/results/sessionInfo.txt")
sessionInfo()
sink()

# Verify it saved
file.exists("C:/Users/eman5/Documents/HD_Project_local/results/sessionInfo.txt")

# Create all needed folders
dir.create("C:/Users/eman5/Documents/HD_Project_local/figures", recursive = TRUE)
dir.create("C:/Users/eman5/Documents/HD_Project_local/scripts", recursive = TRUE)
dir.create("C:/Users/eman5/Documents/HD_Project_local/data", recursive = TRUE)

# Verify all folders exist
dir.exists("C:/Users/eman5/Documents/HD_Project_local/figures")
dir.exists("C:/Users/eman5/Documents/HD_Project_local/scripts")
dir.exists("C:/Users/eman5/Documents/HD_Project_local/data")

# Save all your UMAP plots to figures folder

# Before vs After by condition
ggsave("C:/Users/eman5/Documents/HD_Project_local/figures/UMAP_condition_before_after.png",
       plot = p1 | p2, width = 14, height = 6, dpi = 300)

# Before vs After by region
ggsave("C:/Users/eman5/Documents/HD_Project_local/figures/UMAP_region_before_after.png",
       plot = p3 | p4, width = 14, height = 6, dpi = 300)

# Before vs After by sample
ggsave("C:/Users/eman5/Documents/HD_Project_local/figures/UMAP_sample_before_after.png",
       plot = p5 | p6, width = 14, height = 6, dpi = 300)

# Final clusters
ggsave("C:/Users/eman5/Documents/HD_Project_local/figures/UMAP_clusters_final.png",
       width = 10, height = 8, dpi = 300)