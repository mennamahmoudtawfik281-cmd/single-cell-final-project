# 1. Load Packages

library(Seurat)
library(decoupleR)
library(dorothea)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Load Data
HD_project = readRDS("E:/merged_subclustered_final.rds")

# 3.Load DoRothEA TF-target Network
data(
  dorothea_hs,
  package = "dorothea"
)

net <- dorothea_hs %>%
  +   dplyr::filter(
    +     confidence %in% c("A", "B", "C")
    +   ) %>%
  +   dplyr::select(
    +     source = tf,
    +     target = target,
    +     mor = mor
    +   )

# 4.Get Expression Matrix
  mat <- GetAssayData(
    HD_project,
    assay = "RNA",
    layer = "data"
  )

# 5.Infer TF Activity Using ULM
mat_substance = mat[VariableFeatures(HD_project),]

acts <- run_ulm(
  mat = mat_substance,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)  

# 6.Create TF Activity Matrix
  tf_activity <- acts %>%
  dplyr::select(
    source,
    condition,
    score
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = score
  ) %>%
  tibble::column_to_rownames("source")
  
# 7.Add TF Activity to Seurat Object
  HD_project[["TF_activity"]] <- CreateAssay5Object(
    data = tf_activity
  )

DefaultAssay(HD_project) <- "TF_activity"

# 8.Visualize HSF1 TF Activity
  p_activity <- FeaturePlot(
    HD_project,
    features = "HSF1",
    reduction = "umap_after",
    order = TRUE
  ) +
  scale_color_gradient2(low = "blue",
                        mid = "lightgrey",
                        high = "red",
                        midpoint = 0) +
  ggtitle("HSF1 TF Activity")

p_activity

# 9. Visualize HSF1 TF Gene Expression
  DefaultAssay(HD_project) <- "RNA"

p_expression <- FeaturePlot(
  HD_project,
  features = "HSF1",
  reduction = "umap_after",
  order = TRUE
) +
  scale_color_gradient2(low = "blue",
                        mid = "lightgrey",
                        high = "red",
                        midpoint = 0) +  
  ggtitle("HSF1 Gene Expression")


p_expression

# 10. Calculate Mean TF Activity per Cell Type
  DefaultAssay(HD_project) <- "TF_activity"

tf_mean <- as.data.frame(
  t(
    GetAssayData(
      HD_project,
      assay = "TF_activity"
    )
  )
)

tf_mean$subcluster <- Idents(HD_project)

tf_mean <- tf_mean %>%
  dplyr::group_by(subcluster) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      mean,
      na.rm = TRUE
    )
  ) %>%
  tibble::column_to_rownames("subcluster")
  
# 11.Select Top 10 Variable TFs
  top_tfs <- apply(
    tf_mean,
    2,
    sd,
    na.rm = TRUE
  ) %>%
  sort(
    decreasing = TRUE
  ) %>%
  head(10) %>%
  names()


heatmap_mat <- tf_mean[
  ,
  top_tfs,
  drop = FALSE
]

top_tfs  

# 12.TF Activity Heatmap
  pheatmap(
    heatmap_mat,
    scale = "none",
    border_color = NA,
    main = "Top 10 TF Activities"
  )
  
  

