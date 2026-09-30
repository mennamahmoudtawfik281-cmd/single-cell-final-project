# ==========================================
# 1. Load Required Libraries
# ==========================================

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)


# ==========================================
# 2. Load HIP H5 Files
# ==========================================

hip_dir <- "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/"

if (!dir.exists(hip_dir)) {
  stop("The folder data/hippocampus/ was not found.")
}

# Check available H5 files
list.files(
  hip_dir,
  pattern = "\\.h5$",
  full.names = TRUE
)


# Function to find the correct H5 file by GSM
get_file <- function(gsm) {
  
  files <- list.files(
    hip_dir,
    pattern = paste0(gsm, ".*\\.h5$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    stop(
      paste(
        "No H5 file was found for:",
        gsm
      )
    )
  }
  
  files[1]
}


# ==========================================
# 3. Read HIP Samples
# ==========================================

HIP_Ctrl1 <- Read10X_h5(
  get_file("GSM8610265")
)

HIP_Ctrl2 <- Read10X_h5(
  get_file("GSM8610269")
)

HIP_Ctrl3 <- Read10X_h5(
  get_file("GSM8610270")
)

HIP_HD1 <- Read10X_h5(
  get_file("GSM8610264")
)

HIP_HD2 <- Read10X_h5(
  get_file("GSM8610268")
)

HIP_HD3 <- Read10X_h5(
  get_file("GSM8610276")
)


# ==========================================
# 4. Create Seurat Objects
# ==========================================

HIP_Ctrl1 <- CreateSeuratObject(
  counts = HIP_Ctrl1,
  project = "HIP_Ctrl1"
)

HIP_Ctrl2 <- CreateSeuratObject(
  counts = HIP_Ctrl2,
  project = "HIP_Ctrl2"
)

HIP_Ctrl3 <- CreateSeuratObject(
  counts = HIP_Ctrl3,
  project = "HIP_Ctrl3"
)

HIP_HD1 <- CreateSeuratObject(
  counts = HIP_HD1,
  project = "HIP_HD1"
)

HIP_HD2 <- CreateSeuratObject(
  counts = HIP_HD2,
  project = "HIP_HD2"
)

HIP_HD3 <- CreateSeuratObject(
  counts = HIP_HD3,
  project = "HIP_HD3"
)


# ==========================================
# 5. Raw Cell Counts
# ==========================================

raw_counts <- c(
  HIP_HD1   = ncol(HIP_HD1),
  HIP_HD2   = ncol(HIP_HD2),
  HIP_HD3   = ncol(HIP_HD3),
  HIP_Ctrl1 = ncol(HIP_Ctrl1),
  HIP_Ctrl2 = ncol(HIP_Ctrl2),
  HIP_Ctrl3 = ncol(HIP_Ctrl3)
)

raw_counts


# ==========================================
# 6. Calculate Mitochondrial Percentage
# ==========================================

HIP_Ctrl1[["percent.mt"]] <- PercentageFeatureSet(
  HIP_Ctrl1,
  pattern = "^MT-"
)

HIP_Ctrl2[["percent.mt"]] <- PercentageFeatureSet(
  HIP_Ctrl2,
  pattern = "^MT-"
)

HIP_Ctrl3[["percent.mt"]] <- PercentageFeatureSet(
  HIP_Ctrl3,
  pattern = "^MT-"
)

HIP_HD1[["percent.mt"]] <- PercentageFeatureSet(
  HIP_HD1,
  pattern = "^MT-"
)

HIP_HD2[["percent.mt"]] <- PercentageFeatureSet(
  HIP_HD2,
  pattern = "^MT-"
)

HIP_HD3[["percent.mt"]] <- PercentageFeatureSet(
  HIP_HD3,
  pattern = "^MT-"
)


# ==========================================
# 7. QC Plots Before Filtering
# ==========================================

VlnPlot(
  HIP_Ctrl1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

VlnPlot(
  HIP_Ctrl2,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

VlnPlot(
  HIP_Ctrl3,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

VlnPlot(
  HIP_HD1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

VlnPlot(
  HIP_HD2,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

VlnPlot(
  HIP_HD3,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)


# ==========================================
# 8. QC Filtering
# ==========================================

HIP_Ctrl1 <- subset(
  HIP_Ctrl1,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

HIP_Ctrl2 <- subset(
  HIP_Ctrl2,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

HIP_Ctrl3 <- subset(
  HIP_Ctrl3,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

HIP_HD1 <- subset(
  HIP_HD1,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

HIP_HD2 <- subset(
  HIP_HD2,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

HIP_HD3 <- subset(
  HIP_HD3,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)


# ==========================================
# 9. Cell Counts After QC
# ==========================================

qc_counts <- c(
  HIP_HD1   = ncol(HIP_HD1),
  HIP_HD2   = ncol(HIP_HD2),
  HIP_HD3   = ncol(HIP_HD3),
  HIP_Ctrl1 = ncol(HIP_Ctrl1),
  HIP_Ctrl2 = ncol(HIP_Ctrl2),
  HIP_Ctrl3 = ncol(HIP_Ctrl3)
)

qc_counts


# ==========================================
# 10. Convert to SingleCellExperiment
# ==========================================

sce_Ctrl1 <- as.SingleCellExperiment(HIP_Ctrl1)
sce_Ctrl2 <- as.SingleCellExperiment(HIP_Ctrl2)
sce_Ctrl3 <- as.SingleCellExperiment(HIP_Ctrl3)

sce_HD1 <- as.SingleCellExperiment(HIP_HD1)
sce_HD2 <- as.SingleCellExperiment(HIP_HD2)
sce_HD3 <- as.SingleCellExperiment(HIP_HD3)


# ==========================================
# 11. Detect Doublets
# ==========================================

set.seed(100)

sce_Ctrl1 <- scDblFinder(sce_Ctrl1)

set.seed(100)

sce_Ctrl2 <- scDblFinder(sce_Ctrl2)

set.seed(100)

sce_Ctrl3 <- scDblFinder(sce_Ctrl3)

set.seed(100)

sce_HD1 <- scDblFinder(sce_HD1)

set.seed(100)

sce_HD2 <- scDblFinder(sce_HD2)

set.seed(100)

sce_HD3 <- scDblFinder(sce_HD3)


# ==========================================
# 12. Add Doublet Score and Class
# ==========================================

HIP_Ctrl1$doublet_score <- colData(sce_Ctrl1)$scDblFinder.score

HIP_Ctrl1$doublet_class <- colData(sce_Ctrl1)$scDblFinder.class

HIP_Ctrl2$doublet_score <- colData(sce_Ctrl2)$scDblFinder.score

HIP_Ctrl2$doublet_class <- colData(sce_Ctrl2)$scDblFinder.class

HIP_Ctrl3$doublet_score <- colData(sce_Ctrl3)$scDblFinder.score

HIP_Ctrl3$doublet_class <- colData(sce_Ctrl3)$scDblFinder.class

HIP_HD1$doublet_score <- colData(sce_HD1)$scDblFinder.score

HIP_HD1$doublet_class <- colData(sce_HD1)$scDblFinder.class

HIP_HD2$doublet_score <- colData(sce_HD2)$scDblFinder.score

HIP_HD2$doublet_class <- colData(sce_HD2)$scDblFinder.class

HIP_HD3$doublet_score <- colData(sce_HD3)$scDblFinder.score

HIP_HD3$doublet_class <- colData(sce_HD3)$scDblFinder.class


# ==========================================
# 13. Show Doublet Results
# ==========================================

table(HIP_Ctrl1$doublet_class)

table(HIP_Ctrl2$doublet_class)

table(HIP_Ctrl3$doublet_class)

table(HIP_HD1$doublet_class)

table(HIP_HD2$doublet_class)

table(HIP_HD3$doublet_class)


# ==========================================
# 14. Remove Doublets
# ==========================================

HIP_Ctrl1 <- subset(
  HIP_Ctrl1,
  subset = doublet_class == "singlet"
)

HIP_Ctrl2 <- subset(
  HIP_Ctrl2,
  subset = doublet_class == "singlet"
)

HIP_Ctrl3 <- subset(
  HIP_Ctrl3,
  subset = doublet_class == "singlet"
)

HIP_HD1 <- subset(
  HIP_HD1,
  subset = doublet_class == "singlet"
)

HIP_HD2 <- subset(
  HIP_HD2,
  subset = doublet_class == "singlet"
)

HIP_HD3 <- subset(
  HIP_HD3,
  subset = doublet_class == "singlet"
)


# ==========================================
# 15. Singlet Counts
# ==========================================

singlet_counts <- c(
  HIP_HD1   = ncol(HIP_HD1),
  HIP_HD2   = ncol(HIP_HD2),
  HIP_HD3   = ncol(HIP_HD3),
  HIP_Ctrl1 = ncol(HIP_Ctrl1),
  HIP_Ctrl2 = ncol(HIP_Ctrl2),
  HIP_Ctrl3 = ncol(HIP_Ctrl3)
)

singlet_counts


# ==========================================
# 16. Create Cell Count Summary
# ==========================================

cell_counts <- data.frame(
  Sample = names(raw_counts),
  Raw = raw_counts,
  After_QC = qc_counts[
    names(raw_counts)
  ],
  Doublets_removed =
    qc_counts[names(raw_counts)] -
    singlet_counts[names(raw_counts)],
  After_doublets =
    singlet_counts[names(raw_counts)],
  Percent_kept =
    round(
      singlet_counts[names(raw_counts)] /
        raw_counts *
        100,
      1
    )
)

cell_counts <- rbind(
  cell_counts,
  data.frame(
    Sample = "Total",
    Raw = sum(raw_counts),
    After_QC = sum(qc_counts),
    Doublets_removed =
      sum(qc_counts) -
      sum(singlet_counts),
    After_doublets =
      sum(singlet_counts),
    Percent_kept =
      round(
        sum(singlet_counts) /
          sum(raw_counts) *
          100,
        1
      )
  )
)

print(
  cell_counts,
  row.names = FALSE
)


# ==========================================
# 17. Normalize Data
# ==========================================

HIP_Ctrl1 <- NormalizeData(
  HIP_Ctrl1,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

HIP_Ctrl2 <- NormalizeData(
  HIP_Ctrl2,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

HIP_Ctrl3 <- NormalizeData(
  HIP_Ctrl3,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

HIP_HD1 <- NormalizeData(
  HIP_HD1,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

HIP_HD2 <- NormalizeData(
  HIP_HD2,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

HIP_HD3 <- NormalizeData(
  HIP_HD3,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ==========================================
# 18. Find Variable Features
# ==========================================

HIP_Ctrl1 <- FindVariableFeatures(
  HIP_Ctrl1,
  selection.method = "vst",
  nfeatures = 2000
)

HIP_Ctrl2 <- FindVariableFeatures(
  HIP_Ctrl2,
  selection.method = "vst",
  nfeatures = 2000
)

HIP_Ctrl3 <- FindVariableFeatures(
  HIP_Ctrl3,
  selection.method = "vst",
  nfeatures = 2000
)

HIP_HD1 <- FindVariableFeatures(
  HIP_HD1,
  selection.method = "vst",
  nfeatures = 2000
)

HIP_HD2 <- FindVariableFeatures(
  HIP_HD2,
  selection.method = "vst",
  nfeatures = 2000
)

HIP_HD3 <- FindVariableFeatures(
  HIP_HD3,
  selection.method = "vst",
  nfeatures = 2000
)


# ==========================================
# 19. Variable Feature Plots
# ==========================================

VariableFeaturePlot(HIP_Ctrl1)

VariableFeaturePlot(HIP_Ctrl2)

VariableFeaturePlot(HIP_Ctrl3)

VariableFeaturePlot(HIP_HD1)

VariableFeaturePlot(HIP_HD2)

VariableFeaturePlot(HIP_HD3)


# ==========================================
# 20. Save Cell Count Table
# ==========================================

dir.create(
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/results",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  cell_counts,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/results/HIP_cell_counts.csv",
  row.names = FALSE
)


# ==========================================
# 21. Save Processed HIP Objects
# ==========================================

dir.create(
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  HIP_Ctrl1,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_Ctrl1_processed.rds"
)

saveRDS(
  HIP_Ctrl2,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_Ctrl2_processed.rds"
)

saveRDS(
  HIP_Ctrl3,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_Ctrl3_processed.rds"
)

saveRDS(
  HIP_HD1,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_HD1_processed.rds"
)

saveRDS(
  HIP_HD2,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_HD2_processed.rds"
)

saveRDS(
  HIP_HD3,
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed/HIP_HD3_processed.rds"
)


# ==========================================
# 22. Final Check
# ==========================================

print(cell_counts)

list.files(
  "C:/Users/DELL/Documents/Projects/single-cell-final-project/data/hippocampus/processed",
  full.names = TRUE
)