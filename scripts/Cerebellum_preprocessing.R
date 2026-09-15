# 1.Library required packages
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(hdf5r)

# 2. Load Data
CB_HD1 = Read10X_h5("E:/GSM8610257_14931PJpool01-S_HD_11_NBB17-081_CB_filtered_feature_bc_matrix.h5")
CB_HD1 =CreateSeuratObject(counts = CB_HD1, project = "CB_HD1_scRNA")

CB_HD2 = Read10X_h5("E:/GSM8610258_14931PJpool01-S_HD_15_L-HD241_CB_filtered_feature_bc_matrix.h5")
CB_HD2 = CreateSeuratObject(counts = CB_HD2, project = "CB_HD2_scRNA")

CB_HD3 = Read10X_h5("E:/GSM8610266_14931PJpool02-S_HD_8_L-HD240_CB_filtered_feature_bc_matrix.h5")
CB_HD3 = CreateSeuratObject(counts = CB_HD3, project = "CB_HD3_scRNA")

CB_Ctrl1 = Read10X_h5("E:/GSM8610267_14931PJpool02-S_HD_47_NBB95-310_CB_filtered_feature_bc_matrix.h5")
CB_Ctrl1 = CreateSeuratObject(counts = CB_Ctrl1, project = "CB_Ctrl1_scRNA")

CB_Ctrl2 = Read10X_h5("E:/GSM8610283_14931PJpool04-S_HD_41_NBB17-005_CB_filtered_feature_bc_matrix.h5")
CB_Ctrl2 = CreateSeuratObject(counts = CB_Ctrl2, project = "CB_Ctrl2_scRNA")

CB_Ctrl3 = Read10X_h5("E:/GSM8610288_14931PJpool04-S_HD_24_NBB16-056_CB_filtered_feature_bc_matrix.h5")
CB_Ctrl3 = CreateSeuratObject(counts = CB_Ctrl3, project = "CB_Ctrl3_scRNA")

# 3. calculate perecnt.mt
CB_HD1[["percent.mt"]] = PercentageFeatureSet(CB_HD1,pattern = "^MT-")
CB_HD2[["percent.mt"]] = PercentageFeatureSet(CB_HD2,pattern = "^MT-")
CB_HD3[["percent.mt"]] = PercentageFeatureSet(CB_HD3,pattern = "^MT-")
CB_Ctrl1[["percent.mt"]] = PercentageFeatureSet(CB_Ctrl1,pattern = "^MT-")
CB_Ctrl2[["percent.mt"]] = PercentageFeatureSet(CB_Ctrl2,pattern = "^MT-")
CB_Ctrl3[["percent.mt"]] = PercentageFeatureSet(CB_Ctrl3,pattern = "^MT-")

# 4. Visualization before Filtering
VlnPlot(CB_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_HD1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CB_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_HD2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CB_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_HD3, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CB_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_Ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CB_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_Ctrl2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CB_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))
FeatureScatter(CB_Ctrl3, feature1 = "nCount_RNA", feature2 = "percent.mt")

# 5. Cell Filtering
CB_HD1 = subset(CB_HD1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

CB_HD2 = subset(CB_HD2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

CB_HD3 = subset(CB_HD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

CB_Ctrl1 = subset(CB_Ctrl1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

CB_Ctrl2 = subset(CB_Ctrl2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

CB_Ctrl3 = subset(CB_Ctrl3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CB_Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol= 3))

# 6.convert into SingleCellExperiment
CB_HD1_d = as.SingleCellExperiment(CB_HD1)
CB_HD2_d = as.SingleCellExperiment(CB_HD2)
CB_HD3_d = as.SingleCellExperiment(CB_HD3)
CB_Ctrl1_d = as.SingleCellExperiment(CB_Ctrl1)
CB_Ctrl2_d = as.SingleCellExperiment(CB_Ctrl2)
CB_Ctrl3_d = as.SingleCellExperiment(CB_Ctrl3)

# 7. create doublets
set.seed(100)
CB_HD1_d = scDblFinder(CB_HD1_d)
CB_HD2_d = scDblFinder(CB_HD2_d)
CB_HD3_d = scDblFinder(CB_HD3_d)
CB_Ctrl1_d = scDblFinder(CB_Ctrl1_d)
CB_Ctrl2_d = scDblFinder(CB_Ctrl2_d)
CB_Ctrl3_d = scDblFinder(CB_Ctrl3_d)


# 8.Add calssification and score to data
CB_HD1$doublet_score = colData(CB_HD1_d)$scDblFinder.score
CB_HD1$doublet_class = colData(CB_HD1_d)$scDblFinder.class

CB_HD2$doublet_score = colData(CB_HD2_d)$scDblFinder.score
CB_HD2$doublet_class = colData(CB_HD2_d)$scDblFinder.class

CB_HD3$doublet_score = colData(CB_HD3_d)$scDblFinder.score
CB_HD3$doublet_class = colData(CB_HD3_d)$scDblFinder.class

CB_Ctrl1$doublet_score = colData(CB_Ctrl1_d)$scDblFinder.score
CB_Ctrl1$doublet_class = colData(CB_Ctrl1_d)$scDblFinder.class

CB_Ctrl2$doublet_score = colData(CB_Ctrl2_d)$scDblFinder.score
CB_Ctrl2$doublet_class = colData(CB_Ctrl2_d)$scDblFinder.class

CB_Ctrl3$doublet_score = colData(CB_Ctrl3_d)$scDblFinder.score
CB_Ctrl3$doublet_class = colData(CB_Ctrl3_d)$scDblFinder.class

# 9.show results
table(CB_HD1$doublet_class)
VlnPlot(CB_HD1, features = "doublet_score", group.by = "doublet_class")

table(CB_HD2$doublet_class)
VlnPlot(CB_HD2, features = "doublet_score", group.by = "doublet_class")

table(CB_HD3$doublet_class)
VlnPlot(CB_HD3, features = "doublet_score", group.by = "doublet_class")

table(CB_Ctrl1$doublet_class)
VlnPlot(CB_Ctrl1, features = "doublet_score", group.by = "doublet_class")

table(CB_Ctrl2$doublet_class)
VlnPlot(CB_Ctrl2, features = "doublet_score", group.by = "doublet_class")

table(CB_Ctrl3$doublet_class)
VlnPlot(CB_Ctrl3, features = "doublet_score", group.by = "doublet_class")

# 10.remove doublets
CB_HD1 = subset(CB_HD1, subset = doublet_class == "singlet")
VlnPlot(CB_HD1, features = "doublet_score", group.by = "doublet_class")

CB_HD2 = subset(CB_HD2, subset = doublet_class == "singlet")
VlnPlot(CB_HD2, features = "doublet_score", group.by = "doublet_class")

CB_HD3 = subset(CB_HD3, subset = doublet_class == "singlet")
VlnPlot(CB_HD3, features = "doublet_score", group.by = "doublet_class")

CB_Ctrl1 = subset(CB_Ctrl1, subset = doublet_class == "singlet")
VlnPlot(CB_Ctrl1, features = "doublet_score", group.by = "doublet_class")

CB_Ctrl2 = subset(CB_Ctrl2, subset = doublet_class == "singlet")
VlnPlot(CB_Ctrl2, features = "doublet_score", group.by = "doublet_class")

CB_Ctrl3 = subset(CB_Ctrl3, subset = doublet_class == "singlet")
VlnPlot(CB_Ctrl3, features = "doublet_score", group.by = "doublet_class")

# 11.Normalization
CB_HD1 = NormalizeData(CB_HD1, normalization.method = "LogNormalize", scale.factor = 10000)
CB_HD2 = NormalizeData(CB_HD2, normalization.method = "LogNormalize", scale.factor = 10000)
CB_HD3 = NormalizeData(CB_HD3, normalization.method = "LogNormalize", scale.factor = 10000)
CB_Ctrl1 = NormalizeData(CB_Ctrl1, normalization.method = "LogNormalize", scale.factor = 10000)
CB_Ctrl2 = NormalizeData(CB_Ctrl2, normalization.method = "LogNormalize", scale.factor = 10000)
CB_Ctrl3 = NormalizeData(CB_Ctrl3, normalization.method = "LogNormalize", scale.factor = 10000)

# 12. Highly Variable Features
CB_HD1 = FindVariableFeatures(CB_HD1, selection.method = "vst", nFeatures = 2000)
CB_HD2 = FindVariableFeatures(CB_HD2, selection.method = "vst", nfeatures = 2000)
CB_HD3 = FindVariableFeatures(CB_HD3, selection.method = "vst", nfeatures = 2000)
CB_Ctrl1 = FindVariableFeatures(CB_Ctrl1, selection.method = "vst", nfeatures = 2000)
CB_Ctrl2 = FindVariableFeatures(CB_Ctrl2, selection.method = "vst", nfeatures = 2000)
CB_Ctrl3 = FindVariableFeatures(CB_Ctrl3, selection.method = "vst", nfeatures = 2000)

# 13. Visualize Highly Variable Features
VariableFeaturePlot(CB_HD1)
VariableFeaturePlot(CB_HD2)
VariableFeaturePlot(CB_HD3)
VariableFeaturePlot(CB_Ctrl1)
VariableFeaturePlot(CB_Ctrl2)
VariableFeaturePlot(CB_Ctrl3)
