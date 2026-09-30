# 1. Load Required Libraries

library(Seurat)
library(dplyr)
library(ggplot2)


# 2. Load Final Object

obj <- readRDS(file.choose())

obj

colnames(obj@meta.data)

table(obj$condition)

table(obj$region)

table(obj$orig.ident)

table(obj$broad_cell_type)


# 3. Check CN sample labels

table(
  obj$orig.ident,
  obj$condition
)

sum(
  grepl(
    "^CN_ctrl1_",
    rownames(obj@meta.data)
  )
)

sum(
  grepl(
    "^CN_HD3_",
    rownames(obj@meta.data)
  )
)


# 4. Fix CN_ctrl1 sample label

obj$orig.ident[
  grepl(
    "^CN_ctrl1_",
    rownames(obj@meta.data)
  )
] <- "CN_Ctrl1_scRNA"

obj$condition[
  grepl(
    "^CN_ctrl1_",
    rownames(obj@meta.data)
  )
] <- "Control"


table(
  obj$orig.ident,
  obj$condition
)

table(
  obj$region,
  obj$condition
)


# 5. Calculate cell counts

da_counts <- obj@meta.data %>%
  count(
    orig.ident,
    condition,
    region,
    broad_cell_type,
    name = "n_cells"
  )

print(
  da_counts,
  n = 200
)


# 6. Calculate total cells per sample

sample_totals <- obj@meta.data %>%
  count(
    orig.ident,
    name = "total_cells"
  )

da_counts$total_cells <-
  sample_totals$total_cells[
    match(
      da_counts$orig.ident,
      sample_totals$orig.ident
    )
  ]


# 7. Calculate cell proportions

da_counts$proportion <-
  da_counts$n_cells /
  da_counts$total_cells


# Check proportions

check_proportions <- da_counts %>%
  group_by(orig.ident) %>%
  summarise(
    total_proportion = sum(proportion)
  )

range(
  check_proportions$total_proportion
)


# 8. Select glial cell types

glial_da <- da_counts %>%
  filter(
    broad_cell_type %in% c(
      "Oligodendrocytes",
      "Astrocytes",
      "Microglia",
      "OPCs"
    )
  )

table(
  glial_da$region,
  glial_da$broad_cell_type,
  glial_da$condition
)


# 9. Prepare sample-level table

glial_sample_table <- glial_da %>%
  select(
    orig.ident,
    condition,
    region,
    broad_cell_type,
    n_cells,
    total_cells,
    proportion
  ) %>%
  arrange(
    region,
    broad_cell_type,
    condition,
    orig.ident
  )

glial_sample_table


# 10. Differential abundance analysis

da_results <- glial_sample_table %>%
  group_by(
    region,
    broad_cell_type
  ) %>%
  summarise(
    p_value = wilcox.test(
      proportion ~ condition,
      exact = TRUE
    )$p.value,
    
    control_median =
      median(
        proportion[
          condition == "Control"
        ]
      ),
    
    HD_median =
      median(
        proportion[
          condition == "HD"
        ]
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

da_results


# 11. Add difference and direction

da_final <- da_results %>%
  mutate(
    difference =
      HD_median - control_median,
    
    direction = ifelse(
      difference > 0,
      "Higher in HD",
      "Lower in HD"
    )
  ) %>%
  arrange(
    region,
    broad_cell_type
  )

da_final


# 12. Save DA results

dir.create(
  "results/DA",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  da_final,
  "results/DA/DA_final_results.csv",
  row.names = FALSE
)


# 13. Plot glial cell proportions

da_plot <- ggplot(
  glial_sample_table,
  aes(
    x = condition,
    y = proportion,
    fill = condition
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2
  ) +
  facet_grid(
    broad_cell_type ~ region,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    title = "Glial Cell Proportions",
    x = NULL,
    y = "Cell proportion"
  ) +
  theme_classic()

da_plot


# 14. Save figure

dir.create(
  "results/figures/DA",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename =
    "results/figures/DA/DA_glial_proportions.png",
  plot = da_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# 15. Final check

list.files(
  "results/DA"
)

list.files(
  "results/figures/DA"
)
