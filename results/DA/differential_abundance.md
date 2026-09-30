# Differential Abundance Analysis

## Overview

This analysis compares the relative proportions of glial cell types between Huntington's disease (HD) and Control samples across the four brain regions in the dataset.

The analysis was performed using the processed Seurat object after checking the sample and condition labels. Cell proportions were calculated at the sample level and used for the statistical comparisons.

## Cell Types and Regions

The analysis included four glial cell types:

* Oligodendrocytes
* Astrocytes
* Microglia
* OPCs

The analysis was performed separately for:

* CB
* CN
* HIP
* IFG

## Analysis Steps

The workflow was:

1. Load the processed Seurat object and inspect the metadata.
2. Check the sample and condition labels.
3. Correct the CN sample label before the analysis.
4. Count cells for each sample and broad cell type.
5. Calculate the total number of cells in each sample.
6. Calculate the proportion of each cell type within each sample.
7. Keep the four selected glial cell types.
8. Compare Control and HD proportions for each region and cell type.
9. Apply Benjamini-Hochberg FDR correction.
10. Calculate the difference between the HD and Control median proportions.
11. Generate a plot showing the glial cell proportions.

## CN Label Correction

During the metadata check, the cells with the `CN_ctrl1_` prefix were found under the `CN_HD3` sample label.

The cell names were used to correct the sample and condition labels before continuing with the analysis.

After correction, the sample and condition tables were checked again to confirm the updated labels.

## Cell Proportion Calculation

For each sample, cell-type proportions were calculated as:

```text
cell proportion = number of cells of the cell type / total cells in the sample
```

The proportions were checked to make sure that the cell-type proportions for each sample summed to 1.

## Statistical Analysis

The comparisons were performed using the sample-level proportions.

For each combination of region and glial cell type, the Control and HD samples were compared using a two-sided Wilcoxon rank-sum test.

There were 3 Control and 3 HD samples per region, giving a total of 16 comparisons:

```text
4 glial cell types × 4 regions = 16 comparisons
```

P-values were adjusted using the Benjamini-Hochberg method.

## Results

None of the 16 comparisons reached statistical significance after FDR correction.

The smallest observed p-value was **0.10**, and the smallest FDR was **0.533**.

The main descriptive differences were:

| Region | Cell type        | Control median | HD median |
| ------ | ---------------- | -------------: | --------: |
| CB     | Oligodendrocytes |          6.66% |     1.98% |
| CB     | Astrocytes       |          10.4% |     3.42% |
| CB     | Microglia        |          2.01% |    0.628% |
| CB     | OPCs             |          2.31% |    0.664% |
| CN     | Oligodendrocytes |          39.9% |     53.2% |
| HIP    | Microglia        |          11.0% |     4.96% |
| IFG    | Oligodendrocytes |          24.5% |     31.1% |

These values describe the observed sample-level differences and were not statistically significant after FDR correction.

## Visualization

The glial cell proportions were visualized for Control and HD samples using boxplots with individual sample points.

The plot was separated by brain region and glial cell type to make the regional differences easier to compare.

## Files

The analysis produces:

* `differential_abundance.R`
* `results/DA/DA_final_results.csv`
* `results/figures/DA/DA_glial_proportions.png`

## Notes

This analysis was performed as a sample-level exploratory differential abundance analysis. Because only three Control and three HD samples were available per region, the statistical power is limited, and the observed differences should be interpreted as descriptive trends rather than significant changes in glial abundance.
