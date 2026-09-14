# 1.Library required packages

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(hdf5r)

# 2. Load Data

IFG_HD1 = Read10X_h5("E:/SC_RNA/GSM8610271_14931PJpool03-S_HD_31_NBB17-060_IFG_filtered_feature_bc_matrix.h5")
IFG_HD1 =CreateSeuratObject(
  counts = IFG_HD1,
  project = "IFG_HD1_scRNA")

IFG_HD2 = Read10X_h5("E:/SC_RNA/GSM8610277_14931PJpool03-S_HD_6_L-HD240_FrCx_filtered_feature_bc_matrix.h5")
IFG_HD2 = CreateSeuratObject(
  counts = IFG_HD2,
  project = "IFG_HD2_scRNA")

IFG_HD3 = Read10X_h5("E:/SC_RNA/GSM8610279_14931PJpool03-S_HD_35_L-HD241_FrCx_filtered_feature_bc_matrix.h5")
IFG_HD3 = CreateSeuratObject(
  counts = IFG_HD3,
  project = "IFG_HD3_scRNA")


IFG_Ctrl1= Read10X_h5("E:/SC_RNA/GSM8610255_14931PJpool01-S_HD_27_NBB17-005_IFG_filtered_feature_bc_matrix.h5")
IFG_Ctrl1= CreateSeuratObject(
  counts = IFG_Ctrl1,
  project = "IFG_Ctrl1_scRNA")

IFG_Ctrl2= Read10X_h5("E:/SC_RNA/GSM8610262_14931PJpool02-S_HD_42_NBB95-310_IFG_filtered_feature_bc_matrix.h5")
IFG_Ctrl2= CreateSeuratObject(
  counts = IFG_Ctrl2,
  project = "IFG_Ctrl2_scRNA")

IFG_Ctrl3= Read10X_h5("E:/SC_RNA/GSM8610263_14931PJpool02-S_HD_1_NBB16-056_IFG_filtered_feature_bc_matrix.h5")
IFG_Ctrl3= CreateSeuratObject(
  counts = IFG_Ctrl3,
  project = "IFG_Ctrl3_scRNA")


# 3. calculate perecnt.mt
IFG_HD1[["percent.mt"]] = PercentageFeatureSet(IFG_HD1,pattern = "^MT-")
IFG_HD2[["percent.mt"]] = PercentageFeatureSet(IFG_HD2,pattern = "^MT-")
IFG_HD3[["percent.mt"]] = PercentageFeatureSet(IFG_HD3,pattern = "^MT-")
IFG_Ctrl1[["percent.mt"]] = PercentageFeatureSet(IFG_Ctrl1,pattern = "^MT-")
IFG_Ctrl2[["percent.mt"]] = PercentageFeatureSet(IFG_Ctrl2,pattern = "^MT-")
IFG_Ctrl3[["percent.mt"]] = PercentageFeatureSet(IFG_Ctrl3,pattern = "^MT-")

# 4. Visualization before Filtering
VlnPlot(IFG_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_HD1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(IFG_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_HD2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(IFG_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_HD3, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(IFG_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_Ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(IFG_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_Ctrl2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(IFG_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(IFG_Ctrl3, feature1 = "nCount_RNA", feature2 = "percent.mt")

# 5. Cell Filtering
IFG_HD1 = subset(IFG_HD1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

IFG_HD2 = subset(IFG_HD2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

IFG_HD3 = subset(IFG_HD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

IFG_Ctrl1 = subset(IFG_Ctrl1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

IFG_Ctrl2 = subset(IFG_Ctrl2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

IFG_Ctrl3 = subset(IFG_Ctrl3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(IFG_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

# 6.convert into SingleCellExperiment
IFG_HD1_d = as.SingleCellExperiment(IFG_HD1)
IFG_HD2_d = as.SingleCellExperiment(IFG_HD2)
IFG_HD3_d = as.SingleCellExperiment(IFG_HD3)
IFG_Ctrl1_d = as.SingleCellExperiment(IFG_Ctrl1)
IFG_Ctrl2_d = as.SingleCellExperiment(IFG_Ctrl2)
IFG_Ctrl3_d = as.SingleCellExperiment(IFG_Ctrl3)

# 7. create doublets
set.seed(100)
IFG_HD1_d = scDblFinder(IFG_HD1_d)
IFG_HD2_d = scDblFinder(IFG_HD2_d)
IFG_HD3_d = scDblFinder(IFG_HD3_d) 
IFG_Ctrl1_d = scDblFinder(IFG_Ctrl1_d)
IFG_Ctrl2_d = scDblFinder(IFG_Ctrl2_d)
IFG_Ctrl3_d = scDblFinder(IFG_Ctrl3_d)

# 8.Add calssification and score to data
IFG_HD1$doublet_score = colData(IFG_HD1_d)$scDblFinder.score
IFG_HD1$doublet_class = colData(IFG_HD1_d)$scDblFinder.class

IFG_HD2$doublet_score = colData(IFG_HD2_d)$scDblFinder.score
IFG_HD2$doublet_class = colData(IFG_HD2_d)$scDblFinder.class

IFG_HD3$doublet_score = colData(IFG_HD3_d)$scDblFinder.score
IFG_HD3$doublet_class = colData(IFG_HD3_d)$scDblFinder.class

IFG_Ctrl1$doublet_score = colData(IFG_Ctrl1_d)$scDblFinder.score
IFG_Ctrl1$doublet_class = colData(IFG_Ctrl1_d)$scDblFinder.class

IFG_Ctrl2$doublet_score = colData(IFG_Ctrl2_d)$scDblFinder.score
IFG_Ctrl2$doublet_class = colData(IFG_Ctrl2_d)$scDblFinder.class

IFG_Ctrl3$doublet_score = colData(IFG_Ctrl3_d)$scDblFinder.score
IFG_Ctrl3$doublet_class = colData(IFG_Ctrl3_d)$scDblFinder.class

# 9.show results
table(IFG_HD1$doublet_class)
VlnPlot(IFG_HD1, features = "doublet_score", group.by = "doublet_class")

table(IFG_HD2$doublet_class)
VlnPlot(IFG_HD2, features = "doublet_score", group.by = "doublet_class")

table(IFG_HD3$doublet_class)
VlnPlot(IFG_HD3, features = "doublet_score", group.by = "doublet_class")

table(IFG_Ctrl1$doublet_class)
VlnPlot(IFG_Ctrl1, features = "doublet_score", group.by = "doublet_class")

table(IFG_Ctrl2$doublet_class)
VlnPlot(IFG_Ctrl2, features = "doublet_score", group.by = "doublet_class")

table(IFG_Ctrl3$doublet_class)
VlnPlot(IFG_Ctrl3, features = "doublet_score", group.by = "doublet_class")


# 10.remove doublets
IFG_HD1 = subset(IFG_HD1, subset = doublet_class == "singlet")
VlnPlot(IFG_HD1, features = "doublet_score", group.by = "doublet_class")

IFG_HD2 = subset(IFG_HD2, subset = doublet_class == "singlet")
VlnPlot(IFG_HD2, features = "doublet_score", group.by = "doublet_class")

IFG_HD3 = subset(IFG_HD3, subset = doublet_class == "singlet")
VlnPlot(IFG_HD3, features = "doublet_score", group.by = "doublet_class")

IFG_Ctrl1 = subset(IFG_Ctrl1, subset = doublet_class == "singlet")
VlnPlot(IFG_Ctrl1, features = "doublet_score", group.by = "doublet_class")

IFG_Ctrl2 = subset(IFG_Ctrl2, subset = doublet_class == "singlet")
VlnPlot(IFG_Ctrl2, features = "doublet_score", group.by = "doublet_class")

IFG_Ctrl3 = subset(IFG_Ctrl3, subset = doublet_class == "singlet")
VlnPlot(IFG_Ctrl3, features = "doublet_score", group.by = "doublet_class")

# 11.Normalization
IFG_HD1 = NormalizeData(IFG_HD1, normalization.method = "LogNormalize", scale.factor = 10000)
IFG_HD2 = NormalizeData(IFG_HD2, normalization.method = "LogNormalize", scale.factor = 10000)
IFG_HD3 = NormalizeData(IFG_HD3, normalization.method = "LogNormalize", scale.factor = 10000)
IFG_Ctrl1 = NormalizeData(IFG_Ctrl1, normalization.method = "LogNormalize", scale.factor = 10000)
IFG_Ctrl2 = NormalizeData(IFG_Ctrl2, normalization.method = "LogNormalize", scale.factor = 10000)
IFG_Ctrl3 = NormalizeData(IFG_Ctrl3, normalization.method = "LogNormalize", scale.factor = 10000)

# 12. Highly Variable Features
IFG_HD1 = FindVariableFeatures(IFG_HD1, selection.method = "vst", nFeatures = 2000)
IFG_HD2 = FindVariableFeatures(IFG_HD2, selection.method = "vst", nfeatures = 2000)
IFG_HD3 = FindVariableFeatures(IFG_HD3, selection.method = "vst", nfeatures = 2000)
IFG_Ctrl1 = FindVariableFeatures(IFG_Ctrl1, selection.method = "vst", nfeatures = 2000)
IFG_Ctrl2 = FindVariableFeatures(IFG_Ctrl2, selection.method = "vst", nfeatures = 2000)
IFG_Ctrl3 = FindVariableFeatures(IFG_Ctrl3, selection.method = "vst", nfeatures = 2000)

# 13. Visualize Highly Variable Features
VariableFeaturePlot(IFG_HD1)
VariableFeaturePlot(IFG_HD2)
VariableFeaturePlot(IFG_HD3)
VariableFeaturePlot(IFG_Ctrl1)
VariableFeaturePlot(IFG_Ctrl2)
VariableFeaturePlot(IFG_Ctrl3)


saveRDS(IFG_HD1, "E:/Project 1/IFG_HD1.rds")
saveRDS(IFG_HD2, "E:/Project 1/IFG_HD2.rds")
saveRDS(IFG_HD3, "E:/Project 1/IFG_HD3.rds")
saveRDS(IFG_Ctrl1, "E:/Project 1/IFG_Ctrl1.rds")
saveRDS(IFG_Ctrl2, "E:/Project 1/IFG_Ctrl2.rds")
saveRDS(IFG_Ctrl3, "E:/Project 1/IFG_Ctrl3.rds")

dir.create("E:/Project 1", recursive = TRUE)
save.image("E:/Project 1/full_workspace.RData")
