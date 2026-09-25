# ============================================
# Subclustering
# HD scRNA-seq Analysis
# Following Bøstrand et al. 2024 methods:
# - Microglia       → flower names
# - Astrocytes      → herb names
# - Oligodendroglia → tree names
# - Neurons         → descriptive + regional specificity
# - Vascular cells  → lineage and function
#
# IMPORTANT INSTRUCTIONS BEFORE RUNNING:
# 1. Run ONE cell type at a time
# 2. After each ElbowPlot → update dims_use before continuing
# 3. After each resolution comparison → choose best resolution
# 4. After FindAllMarkers → read top5 markers BEFORE naming
# 5. Run table(Idents(x)) BEFORE RenameIdents to confirm cluster count
# 6. NEVER copy paper names without verifying your own marker evidence
# ============================================

set.seed(42)

# Load libraries
# Note: Install MAST once with BiocManager::install("MAST") if needed
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(MAST)

# Set working directory
setwd("C:/Users/eman5/Documents/HD_Project_local")

# Create output folders if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Load annotated object
merged <- readRDS("merged_annotated.rds")

# Confirm broad cell types and cell counts
table(merged$broad_cell_type)

# Confirm reduction names — use these in all DimPlot calls
Reductions(merged)

# ============================================
# HELPER FUNCTION
# Applies the same clustering pipeline as
# integration step for each cell type subset
# as per Bøstrand et al. 2024
#
# PARAMETERS:
# seurat_obj     → full merged object
# cell_type_name → exact name from broad_cell_type column
# resolutions    → resolutions to test (compare all before choosing)
# npcs           → how many PCs to compute (30 is safe default)
# dims_use       → PCs to use for UMAP and clustering
#                  MUST be updated after reading each ElbowPlot
# ============================================

subcluster_celltype <- function(seurat_obj,
                                cell_type_name,
                                resolutions = c(0.2, 0.3, 0.5),
                                npcs = 30,
                                dims_use = 15) {
  
  cat("\n=============================\n")
  cat("Subclustering:", cell_type_name, "\n")
  cat("Cells:", sum(seurat_obj$broad_cell_type == cell_type_name), "\n")
  cat("Using dims = 1:", dims_use, "\n")
  cat("=============================\n")
  
  subset_obj <- subset(seurat_obj,
                       subset = broad_cell_type == cell_type_name)
  
  subset_obj <- FindVariableFeatures(subset_obj, nfeatures = 2000)
  subset_obj <- ScaleData(subset_obj)
  subset_obj <- RunPCA(subset_obj, npcs = npcs)
  
  # Read ElbowPlot — find where curve flattens, use that number as dims_use
  print(ElbowPlot(subset_obj, ndims = npcs) +
          ggtitle(paste("Elbow Plot -", cell_type_name,
                        "| Current dims_use =", dims_use)))
  
  # Harmony integration by donor
  subset_obj <- RunHarmony(subset_obj,
                           group.by.vars = "orig.ident",
                           reduction = "pca",
                           reduction.save = "harmony")
  
  # UMAP using dims_use chosen from ElbowPlot
  subset_obj <- RunUMAP(subset_obj,
                        reduction = "harmony",
                        dims = 1:dims_use,
                        reduction.name = "umap.sub")
  
  subset_obj <- FindNeighbors(subset_obj,
                              reduction = "harmony",
                              dims = 1:dims_use)
  
  # Test all resolutions — compare before choosing
  for(res in resolutions) {
    subset_obj <- FindClusters(subset_obj, resolution = res)
    cat("Resolution", res, "→",
        length(unique(subset_obj$seurat_clusters)), "clusters\n")
  }
  
  return(subset_obj)
}

# ============================================
# 1. MICROGLIA SUBCLUSTERING
# Named with flower names (Bøstrand et al. 2024)
#
# ElbowPlot result → dims_use = 12
# Resolution chosen → 0.2 (10 clusters, cleanest separation)
# Contaminated clusters removed: 4 (neuronal), 5 (oligo),
#                                7 (astrocyte), 9 (T cell)
# Final subclusters: 0, 1, 2, 3, 6, 8
# Known markers:
#   Homeostatic: P2RY12, CX3CR1, TMEM119
#   DAM: TREM2, APOE, LPL
#   Proliferating: MKI67, TOP2A
#   BAMs: MARCO, SIGLEC1, LILRB4
# ============================================

# Run subclustering — read ElbowPlot before continuing
microglia <- subcluster_celltype(merged, "Microglia",
                                 resolutions = c(0.2, 0.3, 0.5),
                                 npcs = 30,
                                 dims_use = 15)

ggsave("figures/ElbowPlot_microglia.png", width = 8, height = 6, dpi = 300)

# Rerun with correct dims_use = 12 (from ElbowPlot)
microglia <- RunHarmony(microglia,
                        group.by.vars = "orig.ident",
                        reduction = "pca",
                        reduction.save = "harmony")

microglia <- RunUMAP(microglia,
                     reduction = "harmony",
                     dims = 1:12,
                     reduction.name = "umap.sub")

microglia <- FindNeighbors(microglia,
                           reduction = "harmony",
                           dims = 1:12)

microglia <- FindClusters(microglia, resolution = 0.2)
microglia <- FindClusters(microglia, resolution = 0.3)
microglia <- FindClusters(microglia, resolution = 0.5)

# Compare resolutions — choose cleanest separation
p1 <- DimPlot(microglia, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.2",
              label = TRUE) + ggtitle("Microglia - Res 0.2")
p2 <- DimPlot(microglia, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.3",
              label = TRUE) + ggtitle("Microglia - Res 0.3")
p3 <- DimPlot(microglia, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.5",
              label = TRUE) + ggtitle("Microglia - Res 0.5")

p1 | p2 | p3

ggsave("figures/Microglia_resolution_comparison.png",
       width = 18, height = 6, dpi = 300)

# Set chosen resolution (0.2)
Idents(microglia) <- "RNA_snn_res.0.2"

# Check cluster counts BEFORE RenameIdents
table(Idents(microglia))

# Check donor contribution for small clusters
cluster_donor <- table(Idents(microglia), microglia$orig.ident)
cluster_donor_pct <- prop.table(cluster_donor, margin = 1) * 100
round(cluster_donor_pct[c("6", "7", "8", "9"),], 1)
# Result: all clusters pass donor check — kept all

# Find marker genes
microglia_markers <- FindAllMarkers(
  microglia,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)

top5_microglia <- microglia_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Read top5 markers carefully before assigning names
print(top5_microglia, n = 50)

# Remove contaminated clusters
microglia <- subset(microglia,
                    idents = c("0", "1", "2", "3", "6", "8"),
                    invert = FALSE)

cat("Clean microglia:", ncol(microglia), "\n")
table(Idents(microglia))

# Confirm with paper markers on UMAP
FeaturePlot(microglia,
            reduction = "umap.sub",
            features = c("CX3CR1", "P2RY12", "TMEM119",  # homeostatic
                         "APOE", "TREM2",                 # DAM-like
                         "LPL", "ITGAX", "CLEC7A",        # disease-associated
                         "LILRB4",                        # BAMs
                         "MKI67"),                        # proliferating
            ncol = 3)

ggsave("figures/FeaturePlot_microglia_paper_markers.png",
       width = 14, height = 12, dpi = 300)

# Assign flower names based on marker evidence
microglia <- RenameIdents(microglia,
                          "0" = "Mglia_Rose",    # homeostatic — CX3CR1+ P2RY12+ TMEM119+
                          "1" = "Mglia_Violet",  # homeostatic variant — CX3CR1+ P2RY12 moderate
                          "2" = "Mglia_Lily",    # DAM-like — APOE+ TREM2+
                          "3" = "Mglia_Daisy",   # disease-associated — LPL+ CLEC7A+ TREM2+
                          "6" = "Mglia_BAMs",    # border-associated — LILRB4+
                          "8" = "Mglia_Tulip"    # proliferating — MKI67+
)

microglia$subcluster <- Idents(microglia)
table(microglia$subcluster)

DimPlot(microglia,
        reduction = "umap.sub",
        label = TRUE, label.size = 4,
        repel = TRUE, pt.size = 0.5) +
  ggtitle("Microglia Subclusters") + NoLegend()

ggsave("figures/UMAP_microglia_subclusters_final.png",
       width = 8, height = 6, dpi = 300)

saveRDS(microglia, "microglia_subclustered.rds")

write.csv(microglia_markers,
          "results/microglia_markers.csv",
          row.names = FALSE)

write.csv(top5_microglia,
          "results/top5_microglia_markers.csv",
          row.names = FALSE)

# ============================================
# 2. ASTROCYTE SUBCLUSTERING
# Named with herb names (Bøstrand et al. 2024)
#
# ElbowPlot result → dims_use = 15
# Resolution chosen → 0.2 (7 clusters, matches paper's 7 subtypes)
# Contaminated clusters removed: 1 (neuronal), 5 (neuronal),
#                                6 (microglia)
# Final subclusters: 0, 2, 3, 4
# Known markers:
#   Reactive: GFAP, VIM, CD44
#   Homeostatic: AQP4, GJA1, SLC1A2, DIO2
#   Stress-response: HSPH1, HSPA4L, BAG3
# ============================================

astrocytes <- subcluster_celltype(merged, "Astrocytes",
                                  resolutions = c(0.2, 0.3, 0.5),
                                  npcs = 30,
                                  dims_use = 15)

ggsave("figures/ElbowPlot_astrocytes.png", width = 8, height = 6, dpi = 300)

# Rerun with correct dims_use = 15 (confirmed from ElbowPlot)
astrocytes <- RunHarmony(astrocytes,
                         group.by.vars = "orig.ident",
                         reduction = "pca",
                         reduction.save = "harmony")

astrocytes <- RunUMAP(astrocytes,
                      reduction = "harmony",
                      dims = 1:15,
                      reduction.name = "umap.sub")

astrocytes <- FindNeighbors(astrocytes,
                            reduction = "harmony",
                            dims = 1:15)

astrocytes <- FindClusters(astrocytes, resolution = 0.2)
astrocytes <- FindClusters(astrocytes, resolution = 0.3)
astrocytes <- FindClusters(astrocytes, resolution = 0.5)

# Compare resolutions
p1 <- DimPlot(astrocytes, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.2",
              label = TRUE) + ggtitle("Astrocytes - Res 0.2")
p2 <- DimPlot(astrocytes, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.3",
              label = TRUE) + ggtitle("Astrocytes - Res 0.3")
p3 <- DimPlot(astrocytes, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.5",
              label = TRUE) + ggtitle("Astrocytes - Res 0.5")

p1 | p2 | p3

ggsave("figures/Astrocytes_resolution_comparison.png",
       width = 18, height = 6, dpi = 300)

# Set chosen resolution (0.2)
Idents(astrocytes) <- "RNA_snn_res.0.2"

# Check cluster counts BEFORE RenameIdents
table(Idents(astrocytes))

# Check donor contribution
cluster_donor_ast <- table(Idents(astrocytes), astrocytes$orig.ident)
cluster_donor_pct_ast <- prop.table(cluster_donor_ast, margin = 1) * 100
round(apply(cluster_donor_pct_ast, 1, max), 1)

# Check cluster 6 donors contributing ≥2%
round(cluster_donor_pct_ast["6",], 1)
# Result: good donor mixing — keep all clusters

# Find marker genes
astrocyte_markers <- FindAllMarkers(
  astrocytes,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)

top5_astrocytes <- astrocyte_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Read top5 markers carefully before assigning names
print(top5_astrocytes, n = 50)

write.csv(astrocyte_markers,
          "results/astrocyte_markers.csv",
          row.names = FALSE)

write.csv(top5_astrocytes,
          "results/top5_astrocyte_markers.csv",
          row.names = FALSE)

# Remove contaminated clusters (1=neuronal, 5=neuronal, 6=microglia)
astrocytes <- subset(astrocytes,
                     idents = c("0", "2", "3", "4"),
                     invert = FALSE)

cat("Clean astrocytes:", ncol(astrocytes), "\n")
table(Idents(astrocytes))

# Assign herb names based on marker evidence
astrocytes <- RenameIdents(astrocytes,
                           "0" = "Astro_Basil",  # Reactive — CD44+ APLNR+
                           "2" = "Astro_Mint",   # Homeostatic — DIO2+ GJB6+
                           "3" = "Astro_Thyme",  # Stress-response — HSPH1+ HSPA4L+
                           "4" = "Astro_Sage"    # Specialized/regional — DAB2+ PAX3+
)

astrocytes$subcluster <- Idents(astrocytes)
table(astrocytes$subcluster)

DimPlot(astrocytes,
        reduction = "umap.sub",
        label = TRUE, label.size = 4,
        repel = TRUE, pt.size = 0.5) +
  ggtitle("Astrocyte Subclusters") + NoLegend()

ggsave("figures/UMAP_astrocyte_subclusters_final.png",
       width = 8, height = 6, dpi = 300)

saveRDS(astrocytes, "astrocytes_subclustered.rds")

# ============================================
# 3. OLIGODENDROGLIA SUBCLUSTERING
# Named with tree names (Bøstrand et al. 2024)
# Includes BOTH Oligodendrocytes AND OPCs combined
#
# ElbowPlot result → dims_use = 15
# Resolution chosen → 0.3 (6 clusters, good balance)
# Contaminated cluster removed: 5 (astrocyte + neuronal)
# Final subclusters: 0, 1, 2, 3, 4
# Known markers:
#   Mature oligos: MBP, PLP1, MOG, MOBP, OPALIN
#   OPCs: PDGFRA, VCAN, CSPG4
#   COPs: GPR17, ENPP6
#
# NOTE FOR DE TEAM: Bøstrand et al. EXCLUDED cerebellum
# oligodendroglia from DE due to very low CB cell counts.
# ============================================

# Extract BOTH oligodendrocyte populations together
oligo_cells <- subset(merged,
                      subset = broad_cell_type %in%
                        c("Oligodendrocytes", "OPCs"))

cat("Total oligodendroglia cells:", ncol(oligo_cells), "\n")

oligo_cells <- FindVariableFeatures(oligo_cells, nfeatures = 2000)
oligo_cells <- ScaleData(oligo_cells)
oligo_cells <- RunPCA(oligo_cells, npcs = 30)

# Read ElbowPlot before proceeding
ElbowPlot(oligo_cells, ndims = 30) +
  ggtitle("Elbow Plot - Oligodendroglia")

ggsave("figures/ElbowPlot_oligodendroglia.png",
       width = 8, height = 6, dpi = 300)

# dims_use = 15 based on ElbowPlot
oligo_cells <- RunHarmony(oligo_cells,
                          group.by.vars = "orig.ident",
                          reduction = "pca",
                          reduction.save = "harmony")

oligo_cells <- RunUMAP(oligo_cells,
                       reduction = "harmony",
                       dims = 1:15,
                       reduction.name = "umap.sub")

oligo_cells <- FindNeighbors(oligo_cells,
                             reduction = "harmony",
                             dims = 1:15)

for(res in c(0.2, 0.3, 0.5)) {
  oligo_cells <- FindClusters(oligo_cells, resolution = res)
  cat("Resolution", res, "→",
      length(unique(oligo_cells$seurat_clusters)), "clusters\n")
}

# Compare resolutions
p1 <- DimPlot(oligo_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.2",
              label = TRUE) + ggtitle("Oligodendroglia - Res 0.2")
p2 <- DimPlot(oligo_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.3",
              label = TRUE) + ggtitle("Oligodendroglia - Res 0.3")
p3 <- DimPlot(oligo_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.5",
              label = TRUE) + ggtitle("Oligodendroglia - Res 0.5")

p1 | p2 | p3

ggsave("figures/Oligodendroglia_resolution_comparison.png",
       width = 18, height = 6, dpi = 300)

# Set chosen resolution (0.3)
Idents(oligo_cells) <- "RNA_snn_res.0.3"

# Check cluster counts BEFORE RenameIdents
table(Idents(oligo_cells))

# Check donor contribution
cluster_donor_oligo <- table(Idents(oligo_cells), oligo_cells$orig.ident)
cluster_donor_pct_oligo <- prop.table(cluster_donor_oligo, margin = 1) * 100
round(apply(cluster_donor_pct_oligo, 1, max), 1)
# Result: all clusters pass — excellent donor mixing

# Find marker genes
oligo_markers <- FindAllMarkers(
  oligo_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)

top5_oligo <- oligo_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Read top5 markers carefully before assigning names
print(top5_oligo, n = 70)

# Confirm marker specificity with DotPlot
DotPlot(oligo_cells,
        features = c("MBP", "PLP1", "MOBP",
                     "OPALIN", "ENPP6", "PDGFRA",
                     "CSPG4", "GPR17"),
        group.by = "RNA_snn_res.0.3") +
  RotatedAxis() +
  ggtitle("Oligodendroglia Marker Confirmation")

ggsave("figures/DotPlot_oligodendroglia_markers.png",
       width = 10, height = 6, dpi = 300)

write.csv(oligo_markers,
          "results/oligodendroglia_markers.csv",
          row.names = FALSE)

write.csv(top5_oligo,
          "results/top5_oligodendroglia_markers.csv",
          row.names = FALSE)

# Remove contaminated cluster 5 (astrocyte + neuronal: GFAP+ AQP4+ SYT1+)
oligo_cells <- subset(oligo_cells,
                      idents = c("0", "1", "2", "3", "4"),
                      invert = FALSE)

cat("Clean oligodendroglia:", ncol(oligo_cells), "\n")
table(Idents(oligo_cells))

# Assign tree names based on marker evidence
oligo_cells <- RenameIdents(oligo_cells,
                            "0" = "Oligo_Birch",  # Most mature — OPALIN+ MBP+ PLP1+
                            "1" = "Oligo_Oak",    # Stress-response — HSPA1A+ BAG3+
                            "2" = "Oligo_Cedar",  # COPs — ENPP6+ DHCR24+
                            "3" = "Oligo_Elder",  # Metabolically active — NEAT1+ MT-ND1+
                            "4" = "OPC_Ash"       # OPCs — PDGFRA+ CSPG4+
)

oligo_cells$subcluster <- Idents(oligo_cells)
table(oligo_cells$subcluster)

DimPlot(oligo_cells,
        reduction = "umap.sub",
        label = TRUE, label.size = 4,
        repel = TRUE, pt.size = 0.3) +
  ggtitle("Oligodendroglia Subclusters") + NoLegend()

ggsave("figures/UMAP_oligodendroglia_subclusters_final.png",
       width = 8, height = 6, dpi = 300)

saveRDS(oligo_cells, "oligodendroglia_subclustered.rds")

# ============================================
# 4. NEURON SUBCLUSTERING
# Named descriptively with regional specificity
# Regional suffix added ONLY if cluster >50% from one region
#
# ElbowPlot result → dims_use = 20 (neurons most diverse)
# Resolution chosen → 0.3 (15 clusters)
# Removed clusters: 5 (69.7% one donor), 12 (53.7% one donor),
#                   14 (73.1% one donor)
# Final: 12 subclusters across all neuronal populations
# Known markers:
#   Excitatory: SLC17A7, CAMK2A, NEFL, NEFM
#   Inhibitory: GAD2, GAD1, SST, VIP
#   Medium spiny: PPP1R1B, DRD1, DRD2
#   Cerebellar granule: RELN, PAX6, CBLN1
# ============================================

# Extract all neuronal populations together
neuron_cells <- subset(merged,
                       subset = broad_cell_type %in%
                         c("Excitatory Neurons",
                           "Inhibitory Neurons",
                           "Cerebellar Granule Cells",
                           "Medium Spiny Neurons"))

cat("Total neuron cells:", ncol(neuron_cells), "\n")

neuron_cells <- FindVariableFeatures(neuron_cells, nfeatures = 2000)
neuron_cells <- ScaleData(neuron_cells)
neuron_cells <- RunPCA(neuron_cells, npcs = 30)

# Read ElbowPlot — neurons are diverse, may need more dims
ElbowPlot(neuron_cells, ndims = 30) +
  ggtitle("Elbow Plot - Neurons")

ggsave("figures/ElbowPlot_neurons.png", width = 8, height = 6, dpi = 300)

# dims_use = 20 based on ElbowPlot
neuron_cells <- RunHarmony(neuron_cells,
                           group.by.vars = "orig.ident",
                           reduction = "pca",
                           reduction.save = "harmony")

neuron_cells <- RunUMAP(neuron_cells,
                        reduction = "harmony",
                        dims = 1:20,
                        reduction.name = "umap.sub")

neuron_cells <- FindNeighbors(neuron_cells,
                              reduction = "harmony",
                              dims = 1:20)

for(res in c(0.2, 0.3, 0.5)) {
  neuron_cells <- FindClusters(neuron_cells, resolution = res)
  cat("Resolution", res, "→",
      length(unique(neuron_cells$seurat_clusters)), "clusters\n")
}

# Compare resolutions
p1 <- DimPlot(neuron_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.2",
              label = TRUE) + ggtitle("Neurons - Res 0.2")
p2 <- DimPlot(neuron_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.3",
              label = TRUE) + ggtitle("Neurons - Res 0.3")
p3 <- DimPlot(neuron_cells, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.5",
              label = TRUE) + ggtitle("Neurons - Res 0.5")

p1 | p2 | p3

ggsave("figures/Neurons_resolution_comparison.png",
       width = 18, height = 6, dpi = 300)

# Set chosen resolution (0.3)
Idents(neuron_cells) <- "RNA_snn_res.0.3"

# Check cluster counts
table(Idents(neuron_cells))

# Check donor contribution — remove donor-dominated clusters
cluster_donor_neu <- table(Idents(neuron_cells), neuron_cells$orig.ident)
cluster_donor_pct_neu <- prop.table(cluster_donor_neu, margin = 1) * 100
round(apply(cluster_donor_pct_neu, 1, max), 1)

# Check small clusters specifically
round(cluster_donor_pct_neu[c("13", "14"),], 1)

# Check cluster 13 — 44.4% but 7 donors ≥2% → passes threshold
sum(cluster_donor_pct_neu["13",] >= 2)

# Check region distribution BEFORE removing clusters
# Regional suffix added only if cluster >50% from one region
DimPlot(neuron_cells, reduction = "umap.sub",
        group.by = "region", pt.size = 0.3) +
  ggtitle("Neurons - by Region")

ggsave("figures/UMAP_neurons_by_region.png",
       width = 8, height = 6, dpi = 300)

# Proportion table to justify regional names
region_proportions <- prop.table(
  table(Idents(neuron_cells), neuron_cells$region),
  margin = 1) * 100

round(region_proportions, 1)

# Remove donor-dominated clusters 5, 12, 14
neuron_cells <- subset(neuron_cells,
                       idents = c("0", "1", "2", "3", "4",
                                  "6", "7", "8", "9", "10",
                                  "11", "13"),
                       invert = FALSE)

cat("Clean neurons:", ncol(neuron_cells), "\n")
table(Idents(neuron_cells))

# Find marker genes
neuron_markers <- FindAllMarkers(
  neuron_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)

# Filter out removed clusters when getting top5
top5_neurons <- neuron_markers %>%
  filter(!cluster %in% c("5", "12", "14")) %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top5_neurons, n = 60)

# Check ambiguous cluster 10
FeaturePlot(neuron_cells,
            reduction = "umap.sub",
            features = c("PDGFRA", "EYA4", "SNAP25", "SLC17A7"),
            ncol = 2)

ggsave("figures/FeaturePlot_cluster10_check.png",
       width = 10, height = 8, dpi = 300)
# Result: SNAP25+ SLC17A7+ confirms cluster 10 = Excitatory Neurons

# Save markers (filtered to clean clusters only)
write.csv(neuron_markers %>% filter(!cluster %in% c("5", "12", "14")),
          "results/neuron_markers.csv",
          row.names = FALSE)

write.csv(top5_neurons,
          "results/top5_neuron_markers.csv",
          row.names = FALSE)

# Assign descriptive names based on marker + regional evidence
# Regional suffix justified by proportion table results
Idents(neuron_cells) <- "RNA_snn_res.0.3"

neuron_cells <- RenameIdents(neuron_cells,
                             "0"  = "Granule_CB",         # PAX6+ CRTAM+ — CB 96.5%
                             "1"  = "Excitatory_1",       # NEFL+ NEFM+ — HIP+IFG mixed
                             "2"  = "Excitatory_2_IFG",   # LINC00507+ — IFG 59.3%
                             "3"  = "Inhibitory_SST",     # SST+ LHX6+ SOX6+
                             "4"  = "Inhibitory_CB",      # GABRA6+ CBLN3+ — CB 99.5%
                             "6"  = "Inhibitory_VIP",     # VIP+ CHRNA2+
                             "7"  = "Excitatory_3",       # RORB+
                             "8"  = "Excitatory_4",       # FOXP2+
                             "9"  = "Cerebellar_CB",      # TFAP2B+ TFAP2A+ — CB 97.4%
                             "10" = "Excitatory_5",       # EYA4+ SLC17A7+
                             "11" = "Excitatory_HIP",     # FEZF2+ — HIP 64.4%
                             "13" = "Cerebellar_2_CB"     # PPP1R17+ — CB 86.4%
)

neuron_cells$subcluster <- Idents(neuron_cells)
table(neuron_cells$subcluster)

DimPlot(neuron_cells,
        reduction = "umap.sub",
        label = TRUE, label.size = 3,
        repel = TRUE, pt.size = 0.3) +
  ggtitle("Neuron Subclusters") + NoLegend()

ggsave("figures/UMAP_neuron_subclusters_final.png",
       width = 10, height = 8, dpi = 300)

saveRDS(neuron_cells, "neurons_subclustered.rds")

# ============================================
# 5. VASCULAR CELL SUBCLUSTERING
# Named by lineage and function (Bøstrand et al. 2024)
#
# ElbowPlot result → dims_use = 14
# Initial resolutions 0.2/0.3/0.5 too high for small population
# Lower resolutions tested → 0.1 chosen (cleanest)
# Removed: clusters 5 and 6 (<100 cells)
#          then cluster 3 (microglia contamination: C1QA+ CSF1R+)
# Final subclusters: Pericyte, Endothelial, Smooth_Muscle, Endothelial_2
# Known markers:
#   Endothelial: VWF, CLDN5, ABCB1, FLT1, TIE1
#   Pericytes: PDGFRB, NOTCH3, RGS5, ABCC9
# ============================================

vascular <- subset(merged,
                   subset = broad_cell_type == "Vascular Cells")

cat("Vascular cells:", ncol(vascular), "\n")

vascular <- FindVariableFeatures(vascular, nfeatures = 2000)
vascular <- ScaleData(vascular)
vascular <- RunPCA(vascular, npcs = 20)

ElbowPlot(vascular, ndims = 20) +
  ggtitle("Elbow Plot - Vascular Cells")

ggsave("figures/ElbowPlot_vascular.png",
       width = 8, height = 6, dpi = 300)

# dims_use = 14 based on ElbowPlot
vascular <- RunHarmony(vascular,
                       group.by.vars = "orig.ident",
                       reduction = "pca",
                       reduction.save = "harmony")

vascular <- RunUMAP(vascular,
                    reduction = "harmony",
                    dims = 1:14,
                    reduction.name = "umap.sub")

vascular <- FindNeighbors(vascular,
                          reduction = "harmony",
                          dims = 1:14)

# Test standard resolutions first
vascular <- FindClusters(vascular, resolution = 0.2)
vascular <- FindClusters(vascular, resolution = 0.3)
vascular <- FindClusters(vascular, resolution = 0.5)

p1 <- DimPlot(vascular, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.2",
              label = TRUE) + ggtitle("Vascular - Res 0.2")
p2 <- DimPlot(vascular, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.3",
              label = TRUE) + ggtitle("Vascular - Res 0.3")
p3 <- DimPlot(vascular, reduction = "umap.sub",
              group.by = "RNA_snn_res.0.5",
              label = TRUE) + ggtitle("Vascular - Res 0.5")

p1 | p2 | p3

ggsave("figures/Vascular_resolution_comparison.png",
       width = 18, height = 6, dpi = 300)

# Standard resolutions too high for small population — test lower
vascular <- FindClusters(vascular, resolution = 0.1)
vascular <- FindClusters(vascular, resolution = 0.15)

p_01  <- DimPlot(vascular, reduction = "umap.sub",
                 group.by = "RNA_snn_res.0.1",
                 label = TRUE) + ggtitle("Vascular - Res 0.1")
p_015 <- DimPlot(vascular, reduction = "umap.sub",
                 group.by = "RNA_snn_res.0.15",
                 label = TRUE) + ggtitle("Vascular - Res 0.15")

p_01 | p_015

ggsave("figures/Vascular_low_resolution_comparison.png",
       width = 12, height = 6, dpi = 300)

# Set chosen resolution (0.1)
Idents(vascular) <- "RNA_snn_res.0.1"

# Check cluster sizes and donor contribution
table(Idents(vascular))

cluster_donor_vas <- table(Idents(vascular), vascular$orig.ident)
cluster_donor_pct_vas <- prop.table(cluster_donor_vas, margin = 1) * 100
round(apply(cluster_donor_pct_vas, 1, max), 1)

# Check cluster 0 (51% borderline) and small clusters 5, 6
round(cluster_donor_pct_vas["0",], 1)
round(cluster_donor_pct_vas[c("5", "6"),], 1)
# Cluster 0: 51% borderline but multiple donors — keep
# Clusters 5, 6: below 100 cells — remove

# Remove small clusters 5 and 6
vascular <- subset(vascular,
                   idents = c("0", "1", "2", "3", "4"),
                   invert = FALSE)

cat("After removing small clusters:", ncol(vascular), "\n")
table(Idents(vascular))

# Confirm endothelial vs pericyte identity with FeaturePlot
FeaturePlot(vascular, reduction = "umap.sub",
            features = c("VWF", "CLDN5",      # endothelial
                         "ABCB1", "FLT1",      # BBB endothelial
                         "PDGFRB", "NOTCH3",   # pericytes
                         "RGS5", "ABCC9"),     # pericytes
            ncol = 4)

ggsave("figures/FeaturePlot_vascular_markers.png",
       width = 16, height = 8, dpi = 300)

# Find marker genes
vascular_markers <- FindAllMarkers(
  vascular,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)

top5_vascular <- vascular_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Read top5 markers carefully before assigning names
print(top5_vascular, n = 25)
# Cluster 3: C1QA, C1QB, CSF1R → microglia contamination → remove

# Remove microglia contamination (cluster 3)
vascular <- subset(vascular,
                   idents = c("0", "1", "2", "4"),
                   invert = FALSE)

cat("Clean vascular:", ncol(vascular), "\n")
table(Idents(vascular))

# Assign names by lineage and function based on marker evidence
vascular <- RenameIdents(vascular,
                         "0" = "Pericyte",       # TNC+ CRB2+ PDGFRB+
                         "1" = "Endothelial",    # TIE1+ VWF+ CLDN5+
                         "2" = "Smooth_Muscle",  # NOTCH3+ CSPG4+
                         "4" = "Endothelial_2"   # CEMIP+ C7+
)

vascular$subcluster <- Idents(vascular)
table(vascular$subcluster)

DimPlot(vascular,
        reduction = "umap.sub",
        label = TRUE, label.size = 4,
        repel = TRUE, pt.size = 0.5) +
  ggtitle("Vascular Cell Subclusters") + NoLegend()

ggsave("figures/UMAP_vascular_subclusters_final.png",
       width = 8, height = 6, dpi = 300)

# Save markers excluding contamination cluster 3
vascular_markers %>%
  filter(!cluster %in% c("3")) %>%
  write.csv("results/vascular_markers.csv", row.names = FALSE)

vascular_markers %>%
  filter(!cluster %in% c("3")) %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  write.csv("results/top5_vascular_markers.csv", row.names = FALSE)

saveRDS(vascular, "vascular_subclustered.rds")

# ============================================
# 6. ADD ALL SUBCLUSTERS BACK TO MERGED OBJECT
# Combines all subcluster labels into one column
# Cells removed during subclustering retain their
# broad cell type label — this is expected
# ============================================

# Start with broad cell type labels as default
merged$subcluster <- as.character(merged$broad_cell_type)

# Add each cell type's subclusters by matching barcodes
microglia_labels <- as.character(microglia$subcluster)
names(microglia_labels) <- colnames(microglia)
merged$subcluster[names(microglia_labels)] <- microglia_labels

astrocyte_labels <- as.character(astrocytes$subcluster)
names(astrocyte_labels) <- colnames(astrocytes)
merged$subcluster[names(astrocyte_labels)] <- astrocyte_labels

oligo_labels <- as.character(oligo_cells$subcluster)
names(oligo_labels) <- colnames(oligo_cells)
merged$subcluster[names(oligo_labels)] <- oligo_labels

neuron_labels <- as.character(neuron_cells$subcluster)
names(neuron_labels) <- colnames(neuron_cells)
merged$subcluster[names(neuron_labels)] <- neuron_labels

vascular_labels <- as.character(vascular$subcluster)
names(vascular_labels) <- colnames(vascular)
merged$subcluster[names(vascular_labels)] <- vascular_labels

# Verify all subclusters assigned correctly
table(merged$subcluster)

# After building merged$subcluster — identify contaminated cells
# These are cells that kept their broad label despite being
# in a subclustered cell type

# Get barcodes of all clean subclustered cells
clean_barcodes <- c(
  colnames(microglia),
  colnames(astrocytes),
  colnames(oligo_cells),
  colnames(neuron_cells),
  colnames(vascular)
)

# Get barcodes of cells that WERE subclustered (broad type)
# but are NOT in clean objects
subclustered_types <- c("Microglia", "Astrocytes",
                        "Oligodendrocytes", "OPCs",
                        "Excitatory Neurons", "Inhibitory Neurons",
                        "Cerebellar Granule Cells",
                        "Medium Spiny Neurons", "Vascular Cells")

# 7. Find contaminated cells — in subclustered types but not in clean objects
contaminated <- rownames(merged@meta.data)[
  merged$broad_cell_type %in% subclustered_types &
    !rownames(merged@meta.data) %in% clean_barcodes
]

cat("Total contaminated/removed cells:", length(contaminated), "\n")

# Label them explicitly
merged$subcluster[contaminated] <- "Removed_contamination"

# Verify
table(merged$subcluster)

# Visualize all subclusters on main UMAP
# umap_after confirmed from Reductions(merged)
DimPlot(merged,
        reduction = "umap_after",
        group.by = "subcluster",
        label = TRUE,
        label.size = 3,
        pt.size = 0.1,
        repel = TRUE) +
  ggtitle("All Subclusters") +
  NoLegend()

ggsave("figures/UMAP_all_subclusters.png",
       width = 12, height = 10, dpi = 300)

# ============================================
# 8. SAVE FINAL OBJECT AND SESSION INFO
# ============================================

saveRDS(merged, "merged_subclustered_final.rds")

sink("results/sessionInfo_subclustering.txt")
sessionInfo()
sink()