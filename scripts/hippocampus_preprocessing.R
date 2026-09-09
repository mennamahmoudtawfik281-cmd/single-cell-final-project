# 1.Library required packages
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(hdf5r)

# 2. Load Data (Auto-detect exact filenames)
hip_dir <- "data/hippocampus/"

get_file <- function(gsm) {
  files <- list.files(hip_dir, pattern = paste0("^", gsm, ".*\\.h5$"), full.names = TRUE)
  if (length(files) == 0) stop(paste("File not found for GSM:", gsm))
  return(files[1])
}

HIP_Ctrl1 <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610265")), project = "HIP_Ctrl1_scRNA", min.cells = 3)
HIP_Ctrl2 <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610269")), project = "HIP_Ctrl2_scRNA", min.cells = 3)
HIP_Ctrl3 <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610270")), project = "HIP_Ctrl3_scRNA", min.cells = 3)

HIP_HD1   <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610264")), project = "HIP_HD1_scRNA",   min.cells = 3)
HIP_HD2   <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610268")), project = "HIP_HD2_scRNA",   min.cells = 3)
HIP_HD3   <- CreateSeuratObject(counts = Read10X_h5(get_file("GSM8610276")), project = "HIP_HD3_scRNA",   min.cells = 3)


# 3. calculate perecnt.mt
HIP_Ctrl1[["percent.mt"]] <- PercentageFeatureSet(HIP_Ctrl1, pattern = "^MT-")
HIP_Ctrl2[["percent.mt"]] <- PercentageFeatureSet(HIP_Ctrl2, pattern = "^MT-")
HIP_Ctrl3[["percent.mt"]] <- PercentageFeatureSet(HIP_Ctrl3, pattern = "^MT-")
HIP_HD1[["percent.mt"]] <- PercentageFeatureSet(HIP_HD1, pattern = "^MT-")
HIP_HD2[["percent.mt"]] <- PercentageFeatureSet(HIP_HD2, pattern = "^MT-")
HIP_HD3[["percent.mt"]] <- PercentageFeatureSet(HIP_HD3, pattern = "^MT-")


# 4. Visualization before Filtering
VlnPlot(HIP_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_Ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(HIP_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_Ctrl2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(HIP_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_Ctrl3, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(HIP_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_HD1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(HIP_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_HD2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(HIP_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(HIP_HD3, feature1 = "nCount_RNA", feature2 = "percent.mt")


# 5. Cell Filtering
HIP_Ctrl1 <- subset(HIP_Ctrl1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HIP_Ctrl2 <- subset(HIP_Ctrl2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HIP_Ctrl3 <- subset(HIP_Ctrl3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HIP_HD1 <- subset(HIP_HD1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HIP_HD2 <- subset(HIP_HD2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HIP_HD3 <- subset(HIP_HD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(HIP_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# 6.convert into SingleCellExperiment
HIP_Ctrl1_d <- as.SingleCellExperiment(HIP_Ctrl1)
HIP_Ctrl2_d <- as.SingleCellExperiment(HIP_Ctrl2)
HIP_Ctrl3_d <- as.SingleCellExperiment(HIP_Ctrl3)
HIP_HD1_d <- as.SingleCellExperiment(HIP_HD1)
HIP_HD2_d <- as.SingleCellExperiment(HIP_HD2)
HIP_HD3_d <- as.SingleCellExperiment(HIP_HD3)


# 7. create doublets
set.seed(100)

HIP_Ctrl1_d <- scDblFinder(HIP_Ctrl1_d)
HIP_Ctrl2_d <- scDblFinder(HIP_Ctrl2_d)
HIP_Ctrl3_d <- scDblFinder(HIP_Ctrl3_d)
HIP_HD1_d <- scDblFinder(HIP_HD1_d)
HIP_HD2_d <- scDblFinder(HIP_HD2_d)
HIP_HD3_d <- scDblFinder(HIP_HD3_d)


# 8.Add calssification and score to data
HIP_Ctrl1$doublet_score <- colData(HIP_Ctrl1_d)$scDblFinder.score
HIP_Ctrl1$doublet_class <- colData(HIP_Ctrl1_d)$scDblFinder.class

HIP_Ctrl2$doublet_score <- colData(HIP_Ctrl2_d)$scDblFinder.score
HIP_Ctrl2$doublet_class <- colData(HIP_Ctrl2_d)$scDblFinder.class

HIP_Ctrl3$doublet_score <- colData(HIP_Ctrl3_d)$scDblFinder.score
HIP_Ctrl3$doublet_class <- colData(HIP_Ctrl3_d)$scDblFinder.class

HIP_HD1$doublet_score <- colData(HIP_HD1_d)$scDblFinder.score
HIP_HD1$doublet_class <- colData(HIP_HD1_d)$scDblFinder.class

HIP_HD2$doublet_score <- colData(HIP_HD2_d)$scDblFinder.score
HIP_HD2$doublet_class <- colData(HIP_HD2_d)$scDblFinder.class

HIP_HD3$doublet_score <- colData(HIP_HD3_d)$scDblFinder.score
HIP_HD3$doublet_class <- colData(HIP_HD3_d)$scDblFinder.class


# 9.show results
table(HIP_Ctrl1$doublet_class)
VlnPlot(HIP_Ctrl1, features = "doublet_score", group.by = "doublet_class")

table(HIP_Ctrl2$doublet_class)
VlnPlot(HIP_Ctrl2, features = "doublet_score", group.by = "doublet_class")

table(HIP_Ctrl3$doublet_class)
VlnPlot(HIP_Ctrl3, features = "doublet_score", group.by = "doublet_class")

table(HIP_HD1$doublet_class)
VlnPlot(HIP_HD1, features = "doublet_score", group.by = "doublet_class")

table(HIP_HD2$doublet_class)
VlnPlot(HIP_HD2, features = "doublet_score", group.by = "doublet_class")

table(HIP_HD3$doublet_class)
VlnPlot(HIP_HD3, features = "doublet_score", group.by = "doublet_class")


# 10.remove doublets
HIP_Ctrl1 <- subset(HIP_Ctrl1, subset = doublet_class == "singlet")
VlnPlot(HIP_Ctrl1, features = "doublet_score", group.by = "doublet_class")

HIP_Ctrl2 <- subset(HIP_Ctrl2, subset = doublet_class == "singlet")
VlnPlot(HIP_Ctrl2, features = "doublet_score", group.by = "doublet_class")

HIP_Ctrl3 <- subset(HIP_Ctrl3, subset = doublet_class == "singlet")
VlnPlot(HIP_Ctrl3, features = "doublet_score", group.by = "doublet_class")

HIP_HD1 <- subset(HIP_HD1, subset = doublet_class == "singlet")
VlnPlot(HIP_HD1, features = "doublet_score", group.by = "doublet_class")

HIP_HD2 <- subset(HIP_HD2, subset = doublet_class == "singlet")
VlnPlot(HIP_HD2, features = "doublet_score", group.by = "doublet_class")

HIP_HD3 <- subset(HIP_HD3, subset = doublet_class == "singlet")
VlnPlot(HIP_HD3, features = "doublet_score", group.by = "doublet_class")


# 11.Normalization
HIP_Ctrl1 <- NormalizeData(HIP_Ctrl1, normalization.method = "LogNormalize", scale.factor = 10000)
HIP_Ctrl2 <- NormalizeData(HIP_Ctrl2, normalization.method = "LogNormalize", scale.factor = 10000)
HIP_Ctrl3 <- NormalizeData(HIP_Ctrl3, normalization.method = "LogNormalize", scale.factor = 10000)
HIP_HD1 <- NormalizeData(HIP_HD1, normalization.method = "LogNormalize", scale.factor = 10000)
HIP_HD2 <- NormalizeData(HIP_HD2, normalization.method = "LogNormalize", scale.factor = 10000)
HIP_HD3 <- NormalizeData(HIP_HD3, normalization.method = "LogNormalize", scale.factor = 10000)


# 12. Highly Variable Features
HIP_Ctrl1 <- FindVariableFeatures(HIP_Ctrl1, selection.method = "vst", nfeatures = 2000)
HIP_Ctrl2 <- FindVariableFeatures(HIP_Ctrl2, selection.method = "vst", nfeatures = 2000)
HIP_Ctrl3 <- FindVariableFeatures(HIP_Ctrl3, selection.method = "vst", nfeatures = 2000)
HIP_HD1 <- FindVariableFeatures(HIP_HD1, selection.method = "vst", nfeatures = 2000)
HIP_HD2 <- FindVariableFeatures(HIP_HD2, selection.method = "vst", nfeatures = 2000)
HIP_HD3 <- FindVariableFeatures(HIP_HD3, selection.method = "vst", nfeatures = 2000)


# 13. Visualize Highly Variable Features
VariableFeaturePlot(HIP_Ctrl1)
VariableFeaturePlot(HIP_Ctrl2)
VariableFeaturePlot(HIP_Ctrl3)
VariableFeaturePlot(HIP_HD1)
VariableFeaturePlot(HIP_HD2)
VariableFeaturePlot(HIP_HD3)


# 14. Save Processed Objects to exact target directory
if (!dir.exists("data/hippocampus/processed")) {
  dir.create("data/hippocampus/processed", recursive = TRUE)
}

saveRDS(HIP_HD1, "data/hippocampus/processed/HIP_HD1.rds")
saveRDS(HIP_HD2, "data/hippocampus/processed/HIP_HD2.rds")
saveRDS(HIP_HD3, "data/hippocampus/processed/HIP_HD3.rds")
saveRDS(HIP_Ctrl1, "data/hippocampus/processed/HIP_Ctrl1.rds")
saveRDS(HIP_Ctrl2, "data/hippocampus/processed/HIP_Ctrl2.rds")
saveRDS(HIP_Ctrl3, "data/hippocampus/processed/HIP_Ctrl3.rds")