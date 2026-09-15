# 1. Library required Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)


# 2. Load Data

CB_Ctrl1 = readRDS("E:/project scRNA/CB_Ctrl1.rds")
CB_Ctrl2 = readRDS("E:/project scRNA/CB_Ctrl2.rds")
CB_Ctrl3 = readRDS("E:/project scRNA/CB_Ctrl3.rds")
CB_HD1 = readRDS("E:/project scRNA/CB_HD1.rds")
CB_HD2 = readRDS("E:/project scRNA/CB_HD2.rds")
CB_HD3 = readRDS("E:/project scRNA/CB_HD3.rds")
CN_Ctrl1 = readRDS("E:/project scRNA/CN_ctrl1.rds")
CN_Ctrl2 = readRDS("E:/project scRNA/CN_ctrl2.rds")
CN_Ctrl3 = readRDS("E:/project scRNA/CN_ctrl3.rds")
CN_HD1 = readRDS("E:/project scRNA/CN_HD1.rds")
CN_HD2 = readRDS("E:/project scRNA/CN_HD2.rds")
CN_HD3 = readRDS("E:/project scRNA/CN_HD3.rds")
HIP_Ctrl1 = readRDS("E:/project scRNA/HIP_Ctrl1.rds")
HIP_Ctrl2 = readRDS("E:/project scRNA/HIP_Ctrl2.rds")
HIP_Ctrl3 = readRDS("E:/project scRNA/HIP_Ctrl3.rds")
HIP_HD1 = readRDS("E:/project scRNA/HIP_HD1.rds")
HIP_HD2 = readRDS("E:/project scRNA/HIP_HD2.rds")
HIP_HD3 = readRDS("E:/project scRNA/HIP_HD3.rds")
IFG_Ctrl1 = readRDS("E:/project scRNA/IFG_Ctrl1.rds")
IFG_Ctrl2 = readRDS("E:/project scRNA/IFG_Ctrl2.rds")
IFG_Ctrl3 = readRDS("E:/project scRNA/IFG_Ctrl3.rds")
IFG_HD1 = readRDS("E:/project scRNA/IFG_HD1.rds")
IFG_HD2 = readRDS("E:/project scRNA/IFG_HD2.rds")
IFG_HD3 = readRDS("E:/project scRNA/IFG_HD3.rds")


# 3. merge all Seurat objects

HD_pro = merge(CB_Ctrl1, y = list ( CB_Ctrl2, CB_Ctrl3, CB_HD1, CB_HD2, CB_HD3,
                                    CN_Ctrl1, CN_Ctrl2, CN_Ctrl3, CN_HD1, CN_HD2, 
                                    CN_HD3 ,HIP_Ctrl1, HIP_Ctrl2, HIP_Ctrl3, HIP_HD1,
                                    HIP_HD2, HIP_HD3, IFG_Ctrl1, IFG_Ctrl2, IFG_Ctrl3,
                                    IFG_HD1, IFG_HD2, IFG_HD3),  project = "HD_scRNA" )


table(HD_pro$orig.ident)


# 4. Create metadat 

HD_pro$Condition <- c(
  rep("Control", ncol(CB_Ctrl1)),  rep("Control", ncol(CB_Ctrl2)), rep("Control", ncol(CB_Ctrl3)),
  rep("HD", ncol(CB_HD1)),        rep("HD", ncol(CB_HD2)),        rep("HD", ncol(CB_HD3)),
  rep("Control", ncol(CN_Ctrl1)), rep("Control", ncol(CN_Ctrl2)), rep("Control", ncol(CN_Ctrl3)),
  rep("HD", ncol(CN_HD1)),        rep("HD", ncol(CN_HD2)),        rep("HD", ncol(CN_HD3)),
  rep("Control", ncol(HIP_Ctrl1)),rep("Control", ncol(HIP_Ctrl2)),rep("Control", ncol(HIP_Ctrl3)),
  rep("HD", ncol(HIP_HD1)),       rep("HD", ncol(HIP_HD2)),       rep("HD", ncol(HIP_HD3)),
  rep("Control", ncol(IFG_Ctrl1)),rep("Control", ncol(IFG_Ctrl2)),rep("Control", ncol(IFG_Ctrl3)),
  rep("HD", ncol(IFG_HD1)),       rep("HD", ncol(IFG_HD2)),       rep("HD", ncol(IFG_HD3))
)
HD_pro$Region = c(
  rep("Cerebellum", ncol(CB_Ctrl1)), rep("Cerebellum", ncol(CB_Ctrl2)), rep("Cerebellum", ncol(CB_Ctrl3)),
  rep("Cerebellum", ncol(CB_HD1)), rep("Cerebellum", ncol(CB_HD2)), rep("Cerebellum", ncol(CB_HD3)),
  rep("Caudata Nucleus", ncol(CN_Ctrl1)), rep("Caudata Nucleus", ncol(CN_Ctrl2)), rep("Caudata Nucleus", ncol(CN_Ctrl3)),
  rep("Caudata Nucleus", ncol(CN_HD1)), rep("Caudata Nucleus", ncol(CN_HD2)), rep("Caudata Nucleus", ncol(CN_HD3)),
  rep("Hippocampus", ncol(HIP_Ctrl1)), rep("Hippocampus", ncol(HIP_Ctrl2)), rep("Hippocampus", ncol(HIP_Ctrl3)),
  rep("Hippocampus", ncol(HIP_HD1)), rep("Hippocampus", ncol(HIP_HD2)), rep("Hippocampus", ncol(HIP_HD3)),
  rep("Inferior Frontal", ncol(IFG_Ctrl1)), rep("Inferior Frontal", ncol(IFG_Ctrl2)), rep("Inferior Frontal", ncol(IFG_Ctrl3)),
  rep("Inferior Frontal", ncol(IFG_HD1)), rep("Inferior Frontal", ncol(IFG_HD2)), rep("Inferior Frontal", ncol(IFG_HD3))
)


head(HD_pro@meta.data)
table(HD_pro@meta.data$Region)

# 5.Scaling
HD_pro = HD_pro %>% ScaleData(verbose = FALSE)

# 6. PCA
HD_pro = HD_pro %>% RunPCA(features = VariableFeatures(HD_pro), npcs = 30, verbose = FALSE)

# 7.ElbowPlot
ElbowPlot(HD_pro, ndims = 30)

# 8. Run Umap before Integration
HD_pro = HD_pro %>% RunUMAP(reduction = "pca", dims = 1:20, reduction.name = "Umap_before")

# 9. Visualization according to difference in Condition
b_before1= DimPlot(HD_pro, reduction = "Umap_before", group.by = "Condition", pt.size = 0.1) +
  ggtitle("Before Integration")

b_before1

# 10. Visualization according to diffenence in Regions
b_before2= DimPlot(HD_pro, reduction = "Umap_before", group.by = "Region", pt.size = 0.1) +
  ggtitle("Before Integration")
b_before2  

# 11. Integration
HD_pro = RunHarmony(HD_pro, group.by.vars = "orig.ident", plot_convergence = TRUE, max_iter = 10)

# 12. Run Umap after Integration
HD_pro = HD_pro %>% RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "Umap_after")

# 13. Visualization according to difference in Condition
b_after1= DimPlot(HD_pro, reduction = "Umap_after", group.by = "Condition", pt.size = 0.1) +
  ggtitle("After Integration")
b_after1

# 14. Visualization according to diffenence in Regions
b_after2= DimPlot(HD_pro, reduction = "Umap_after", group.by = "Region", pt.size = 0.1) +
  ggtitle("After Integration")
b_after2

# 15. Visualization before and after Integration
b_before1 + b_after1
b_before2 + b_after2

# 16. Clustering
HD_pro = FindNeighbors(HD_pro, reduction = "harmony", dims = 1:20)
HD_pro = FindClusters(HD_pro, resolution = 0.5)

# 17. Visualization after Integration
DimPlot(HD_pro, reduction = "Umap_after", group.by = "seurat_clusters",
        label = TRUE, pt.size = 0.1) + ggtitle("Clusters after Harmony Integration")

