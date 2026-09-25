# ============================================
# Cell Type Annotation
# HD scRNA-seq Analysis
# Following Bøstrand et al. 2024 methods
#
# Input:  merged_filtered.rds (from Script 02)
# Output: merged_annotated.rds
#
# Annotation evidence hierarchy:
# 1. DotPlot — PRIMARY (exact expression per cluster)
# 2. FeaturePlot — SECONDARY (spatial confirmation)
# 3. VlnPlot — ADDITIONAL (distribution confirmation)
# 4. FindAllMarkers — CLUSTER DE (supports annotation)
# ============================================

set.seed(42)

# Load libraries
# Note: Install MAST once with BiocManager::install("MAST") if needed
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(MAST)

# Set working directory
setwd("C:/Users/eman5/Documents/HD_Project_local")

# Create output folders if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Load filtered object
merged <- readRDS("merged_filtered.rds")

# Set resolution 0.25 as active identity — used throughout this script
Idents(merged) <- "RNA_snn_res.0.25"

# ============================================
# Step 1 - Cluster numbers UMAP
# Visualize cluster locations before annotation
# ============================================

DimPlot(merged,
        reduction = "umap_after",
        label = TRUE,
        label.size = 5,
        pt.size = 0.1) +
  ggtitle("Clusters - Resolution 0.25") +
  NoLegend()

ggsave("figures/UMAP_cluster_numbers.png",
       width = 10, height = 8, dpi = 300)

# ============================================
# Step 2 - DotPlot (PRIMARY annotation evidence)
# Most reliable method — shows exact expression
# level and percentage per cluster
# Used as primary evidence before FeaturePlot
# ============================================

dp <- DotPlot(merged,
              features = c(
                "ALDH1L1", "GFAP",      # Astrocytes
                "CX3CR1",  "ITGAM",     # Microglia
                "MOG",     "MBP",       # Oligodendrocytes
                "PDGFRA",  "GPR17",     # OPCs/COPs
                "VWF",     "CLDN5",     # Endothelial
                "PDGFRB",  "NOTCH3",    # Pericytes
                "SLC17A7", "GAD2",      # Excitatory/Inhibitory neurons
                "PPP1R1B", "RELN",      # Medium spiny/Cerebellar neurons
                "CD2",     "CD8A",      # Peripheral immune
                "SNAP25"                # Pan-neuronal
              ),
              cols = c("lightblue", "darkred")) +
  RotatedAxis() +
  ggtitle("Canonical Marker Expression per Cluster")

dp

ggsave("figures/DotPlot_canonical_markers.png",
       width = 16, height = 8, dpi = 300)

# Pull exact DotPlot values for precise annotation
# Used to confirm cell type assignments per cluster
write.csv(dp$data,
          "results/dotplot_data.csv",
          row.names = FALSE)

# ============================================
# Step 3 - FeaturePlot (SECONDARY evidence)
# Confirms DotPlot findings spatially
# Canonical markers from Bøstrand et al. 2024
# ============================================

FeaturePlot(merged,
            reduction = "umap_after",
            features = c(
              "ALDH1L1",  # Astrocytes
              "GFAP",     # Astrocytes
              "CX3CR1",   # Microglia
              "ITGAM",    # Microglia
              "MOG",      # Oligodendrocytes
              "MBP",      # Oligodendrocytes
              "PDGFRA",   # OPCs
              "GPR17",    # COPs
              "VWF",      # Endothelial
              "CLDN5",    # Endothelial
              "PDGFRB",   # Pericytes
              "NOTCH3",   # Pericytes
              "SLC17A7",  # Excitatory neurons
              "GAD2",     # Inhibitory neurons
              "PPP1R1B",  # Medium spiny neurons (CN)
              "RELN"      # Cerebellar granule cells (CB)
            ),
            ncol = 4)

ggsave("figures/FeaturePlot_canonical_markers.png",
       width = 20, height = 20, dpi = 300)

# ============================================
# Step 4 - VlnPlot (ADDITIONAL confirmation)
# Distribution of marker expression per cluster
# ============================================

VlnPlot(merged,
        features = c("SNAP25", "GFAP", "MBP",
                     "CX3CR1", "PDGFRA", "CLDN5"),
        ncol = 3,
        pt.size = 0)

ggsave("figures/VlnPlot_key_markers.png",
       width = 16, height = 8, dpi = 300)

# ============================================
# Step 5 - Check ambiguous clusters
# Cluster 6 and 7: VlnPlot to confirm identity
# Cluster 7: OLIG lineage markers confirm OPCs
# Cluster 10: FeaturePlot confirms vascular cells
# ============================================

# Reset identity before checking ambiguous clusters
Idents(merged) <- "RNA_snn_res.0.25"

# Check clusters 6 and 7
VlnPlot(merged,
        features = c("MBP", "SNAP25", "GAD2",
                     "SLC17A7", "GFAP", "CX3CR1"),
        ncol = 3,
        pt.size = 0,
        idents = c("6", "7"))

ggsave("figures/VlnPlot_clusters_6_7.png",
       width = 14, height = 6, dpi = 300)

# Confirm cluster 7 as OPCs using oligodendrocyte lineage markers
VlnPlot(merged,
        features = c("OLIG1", "OLIG2", "SOX10",
                     "NRGN",  "SYT1",  "STMN2"),
        ncol = 3,
        pt.size = 0,
        idents = "7")

ggsave("figures/VlnPlot_cluster_7_OPC_markers.png",
       width = 14, height = 6, dpi = 300)

# Confirm cluster 10 as vascular cells
FeaturePlot(merged,
            reduction = "umap_after",
            features = c("GFAP", "VWF",
                         "CLDN5", "PDGFRB"),
            ncol = 2)

ggsave("figures/FeaturePlot_cluster_10_vascular.png",
       width = 10, height = 8, dpi = 300)

# ============================================
# Step 6 - FindAllMarkers using MAST
# Parameters from Bøstrand et al. 2024:
# min.pct = 0.25 (expressed in >= 25% of cluster)
# logfc.threshold = 0.25
# top 5 markers ranked by log2FC
# ============================================

all_markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500,
  test.use = "MAST"
)

# Top 5 markers per cluster ranked by log2FC
# as per Bøstrand et al. 2024
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Print results — read carefully before assigning cell type names
print(top5_markers, n = 70)

# Add cell type labels to marker results for reference
# Note: this labels the CSV files only
# Actual annotation is applied to Seurat object in Step 7
cell_type_labels <- c(
  "0"  = "Oligodendrocytes",
  "1"  = "Cerebellar Granule Cells",
  "2"  = "Excitatory Neurons",
  "3"  = "Astrocytes",
  "4"  = "Microglia",
  "5"  = "Medium Spiny Neurons",
  "6"  = "Inhibitory Neurons",
  "7"  = "OPCs",
  "8"  = "Inhibitory Neurons",
  "10" = "Vascular Cells",
  "13" = "Excitatory Neurons",
  "14" = "Inhibitory Neurons",
  "15" = "Peripheral Immune Cells",
  "16" = "Excitatory Neurons"
)

all_markers$cell_type <- cell_type_labels[as.character(all_markers$cluster)]
top5_markers$cell_type <- cell_type_labels[as.character(top5_markers$cluster)]

# Save marker results
write.csv(all_markers,
          "results/all_cluster_markers_annotated.csv",
          row.names = FALSE)

write.csv(top5_markers,
          "results/top5_cluster_markers_annotated.csv",
          row.names = FALSE)

# ============================================
# Step 7 - Assign final cell type labels
# Based on DotPlot + FeaturePlot + cluster DE
# evidence from Steps 2-6
# Labels follow Bøstrand et al. 2024 cell types
# ============================================

# Set to cluster numbers before RenameIdents
Idents(merged) <- "RNA_snn_res.0.25"

# Rename clusters with confirmed cell type labels
merged <- RenameIdents(merged,
                       "0"  = "Oligodendrocytes",
                       "1"  = "Cerebellar Granule Cells",
                       "2"  = "Excitatory Neurons",
                       "3"  = "Astrocytes",
                       "4"  = "Microglia",
                       "5"  = "Medium Spiny Neurons",
                       "6"  = "Inhibitory Neurons",
                       "7"  = "OPCs",
                       "8"  = "Inhibitory Neurons",
                       "10" = "Vascular Cells",
                       "13" = "Excitatory Neurons",
                       "14" = "Inhibitory Neurons",
                       "15" = "Peripheral Immune Cells",
                       "16" = "Excitatory Neurons"
)

# Store annotation in metadata
merged$broad_cell_type <- Idents(merged)

# Verify cell type counts
table(merged$broad_cell_type)

# ============================================
# Step 8 - Visualize final annotation
# ============================================

DimPlot(merged,
        reduction = "umap_after",
        group.by = "broad_cell_type",
        label = TRUE,
        label.size = 4,
        pt.size = 0.1,
        repel = TRUE) +
  ggtitle("Broad Cell Type Annotation") +
  NoLegend()

ggsave("figures/UMAP_broad_annotation_final.png",
       width = 10, height = 8, dpi = 300)

# ============================================
# Step 9 - Save annotated object and session info
# ============================================

saveRDS(merged, "merged_annotated.rds")

sink("results/sessionInfo_annotation.txt")
sessionInfo()
sink()