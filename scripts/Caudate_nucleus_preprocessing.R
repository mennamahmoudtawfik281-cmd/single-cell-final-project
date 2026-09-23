
# 1.Library required packages
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(hdf5r)

# 2. Load Data

CN_HD1 = Read10X_h5("data/CN/raw/GSM8610254_14931PJpool01-S_HD_7_L-HD240_CN_filtered_feature_bc_matrix.h5")
CN_HD1 = CreateSeuratObject(counts = CN_HD1, project = "CN_HD1_scRNA")

CN_HD2 = Read10X_h5("data/CN/raw/GSM8610272_14931PJpool03-S_HD_30_L_HD241-CN_CN_filtered_feature_bc_matrix.h5")
CN_HD2 = CreateSeuratObject(counts = CN_HD2, project = "CN_HD2_scRNA")

CN_HD3 = Read10X_h5("data/CN/raw/GSM8610275_14931PJpool03-S_HD_28_NBB17-060_CN_filtered_feature_bc_matrix.h5")
CN_HD3 = CreateSeuratObject(counts = CN_HD3, project = "CN_HD3_scRNA")

CN_ctrl1 = Read10X_h5("data/CN/raw/GSM8610253_14931PJpool01-S_HD_43_EBBSD025_13_CN_filtered_feature_bc_matrix.h5")
CN_ctrl1 = CreateSeuratObject(counts = CN_ctrl1, project = "CN_ctrl1_scRNA")

CN_ctrl2 = Read10X_h5("data/CN/raw/GSM8610259_14931PJpool01-S_HD_29_NBB17-005_CN_filtered_feature_bc_matrix.h5")
CN_ctrl2 = CreateSeuratObject(counts = CN_ctrl2, project = "CN_ctrl2_scRNA")

CN_ctrl3 = Read10X_h5("data/CN/raw/GSM8610261_14931PJpool01-S_HD_44_NBB16-056_CN_filtered_feature_bc_matrix.h5")
CN_ctrl3 = CreateSeuratObject(counts = CN_ctrl3, project = "CN_ctrl3_scRNA")




# 3. calculate perecnt.mt
CN_HD1[["percent.mt"]] = PercentageFeatureSet(CN_HD1,pattern  = "^MT-")
CN_HD2[["percent.mt"]] = PercentageFeatureSet(CN_HD2,pattern  = "^MT-")
CN_HD3[["percent.mt"]] = PercentageFeatureSet(CN_HD3,pattern  = "^MT-")

CN_ctrl1[["percent.mt"]] = PercentageFeatureSet(CN_ctrl1,pattern  = "^MT-")
CN_ctrl2[["percent.mt"]] = PercentageFeatureSet(CN_ctrl2,pattern  = "^MT-")
CN_ctrl3[["percent.mt"]] = PercentageFeatureSet(CN_ctrl3,pattern  = "^MT-")


# 4. Visualization before Filtering
VlnPlot(CN_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_HD1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CN_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_HD2, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CN_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_HD3, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CN_ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")

VlnPlot(CN_ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_ctrl2, feature1 = "nCount_RNA", feature2 = "percent.mt")


VlnPlot(CN_ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)
FeatureScatter(CN_ctrl3, feature1 = "nCount_RNA", feature2 = "percent.mt")



# 5. Cell Filtering
CN_HD1 = subset(CN_HD1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_HD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

CN_HD2 = subset(CN_HD2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_HD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

CN_HD3 = subset(CN_HD3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_HD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

CN_ctrl1 = subset(CN_ctrl1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

CN_ctrl2 = subset(CN_ctrl2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)

CN_ctrl3 = subset(CN_ctrl3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(CN_ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol= 3)




# 6.convert into SingleCellExperiment
CN_HD1_d = as.SingleCellExperiment(CN_HD1)
CN_HD2_d = as.SingleCellExperiment(CN_HD2)
CN_HD3_d = as.SingleCellExperiment(CN_HD3)

CN_ctrl1_d = as.SingleCellExperiment(CN_ctrl1)
CN_ctrl2_d = as.SingleCellExperiment(CN_ctrl2)
CN_ctrl3_d = as.SingleCellExperiment(CN_ctrl3)

# 7. create doublets
set.seed(100)
CN_HD1_d = scDblFinder(CN_HD1_d)
CN_HD2_d = scDblFinder(CN_HD2_d)
CN_HD3_d = scDblFinder(CN_HD3_d)

CN_ctrl1_d = scDblFinder(CN_ctrl1_d)
CN_ctrl2_d = scDblFinder(CN_ctrl2_d)
CN_ctrl3_d = scDblFinder(CN_ctrl3_d)



# 8.Add calssification and score to data
CN_HD1$doublet_score = colData(CN_HD1_d)$scDblFinder.score
CN_HD1$doublet_class = colData(CN_HD1_d)$scDblFinder.class


CN_HD2$doublet_score = colData(CN_HD2_d)$scDblFinder.score
CN_HD2$doublet_class = colData(CN_HD2_d)$scDblFinder.class

CN_HD3$doublet_score = colData(CN_HD3_d)$scDblFinder.score
CN_HD3$doublet_class = colData(CN_HD3_d)$scDblFinder.class


CN_ctrl1$doublet_score = colData(CN_ctrl1_d)$scDblFinder.score
CN_ctrl1$doublet_class = colData(CN_ctrl1_d)$scDblFinder.class

CN_ctrl2$doublet_score = colData(CN_ctrl2_d)$scDblFinder.score
CN_ctrl2$doublet_class = colData(CN_ctrl2_d)$scDblFinder.class

CN_ctrl3$doublet_score = colData(CN_ctrl3_d)$scDblFinder.score
CN_ctrl3$doublet_class = colData(CN_ctrl3_d)$scDblFinder.class


# 9.show results
table(CN_HD1$doublet_class)
VlnPlot(CN_HD1, features = "doublet_score", group.by = "doublet_class")

table(CN_HD2$doublet_class)
VlnPlot(CN_HD2, features = "doublet_score", group.by = "doublet_class")

table(CN_HD3$doublet_class)
VlnPlot(CN_HD3, features = "doublet_score", group.by = "doublet_class")


table(CN_ctrl1$doublet_class)
VlnPlot(CN_ctrl1, features = "doublet_score", group.by = "doublet_class")

table(CN_ctrl2$doublet_class)
VlnPlot(CN_ctrl2, features = "doublet_score", group.by = "doublet_class")

table(CN_ctrl3$doublet_class)
VlnPlot(CN_ctrl3, features = "doublet_score", group.by = "doublet_class")


# 10.remove doublets
CN_HD1 = subset(CN_HD1, subset = doublet_class == "singlet")
VlnPlot(CN_HD1, features = "doublet_score", group.by = "doublet_class")

CN_HD2 = subset(CN_HD2, subset = doublet_class == "singlet")
VlnPlot(CN_HD2, features = "doublet_score", group.by = "doublet_class")

CN_HD3 = subset(CN_HD3, subset = doublet_class == "singlet")
VlnPlot(CN_HD3, features = "doublet_score", group.by = "doublet_class")

CN_ctrl1 = subset(CN_ctrl1, subset = doublet_class == "singlet")
VlnPlot(CN_ctrl1, features = "doublet_score", group.by = "doublet_class")

CN_ctrl2 = subset(CN_ctrl2, subset = doublet_class == "singlet")
VlnPlot(CN_ctrl2, features = "doublet_score", group.by = "doublet_class")

CN_ctrl3 = subset(CN_ctrl3, subset = doublet_class == "singlet")
VlnPlot(CN_ctrl3, features = "doublet_score", group.by = "doublet_class")


# 11.Normalization
CN_HD1 = NormalizeData(CN_HD1, normalization.method = "LogNormalize", scale.factor = 10000)
CN_HD2 = NormalizeData(CN_HD2, normalization.method = "LogNormalize", scale.factor = 10000)
CN_HD3 = NormalizeData(CN_HD3, normalization.method = "LogNormalize", scale.factor = 10000)

CN_ctrl1 = NormalizeData(CN_ctrl1, normalization.method = "LogNormalize", scale.factor = 10000)
CN_ctrl2 = NormalizeData(CN_ctrl2, normalization.method = "LogNormalize", scale.factor = 10000)
CN_ctrl3 = NormalizeData(CN_ctrl3, normalization.method = "LogNormalize", scale.factor = 10000)


# 12. Highly Variable Features
CN_HD1 = FindVariableFeatures(CN_HD1, selection.method = "vst", nfeatures = 2000)
CN_HD2 = FindVariableFeatures(CN_HD2, selection.method = "vst", nfeatures = 2000)
CN_HD3 = FindVariableFeatures(CN_HD3, selection.method = "vst", nfeatures = 2000)

CN_ctrl1 = FindVariableFeatures(CN_ctrl1, selection.method = "vst", nfeatures = 2000)
CN_ctrl2 = FindVariableFeatures(CN_ctrl2, selection.method = "vst", nfeatures = 2000)
CN_ctrl3 = FindVariableFeatures(CN_ctrl3, selection.method = "vst", nfeatures = 2000)


# 13. Visualize Highly Variable Features
VariableFeaturePlot(CN_HD1)
VariableFeaturePlot(CN_HD2)
VariableFeaturePlot(CN_HD3)

VariableFeaturePlot(CN_ctrl1)
VariableFeaturePlot(CN_ctrl2)
VariableFeaturePlot(CN_ctrl3)




# Save your processed Seurat object
dir.create("data/CN/processed", recursive = TRUE)

saveRDS(CN_HD1, "data/CN/processed/CN_HD1.rds")
saveRDS(CN_HD2, "data/CN/processed/CN_HD2.rds")
saveRDS(CN_HD3, "data/CN/processed/CN_HD3.rds")
saveRDS(CN_ctrl1, "data/CN/processed/CN_ctrl1.rds")
saveRDS(CN_ctrl2, "data/CN/processed/CN_ctrl2.rds")
saveRDS(CN_ctrl3, "data/CN/processed/CN_ctrl3.rds")

