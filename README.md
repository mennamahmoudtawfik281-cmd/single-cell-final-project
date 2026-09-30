# Glial Changes Across Four Brain Regions in Huntington's Disease: a snRNA-seq Re-analysis

> Group project re-analysing Bøstrand et al. (2024), *Mapping the glial transcriptome in Huntington's disease using snRNAseq*. We re-process the authors' public data for four brain regions: **caudate nucleus (CN)**, **cerebellum (CB)**, **inferior frontal gyrus / frontal cortex (IFG)** and **hippocampus (HIP)**. We integrate, annotate and subcluster the data, then test how glial cell–cell communication changes in HD.

**Our question:** Do ligand–receptor interactions between stress-response and disease-associated glial subpopulations change between HD and control, and is there evidence that these populations coordinate rather than respond independently?

### Key findings

- **Glial communication is rewired rather than uniformly increased or decreased in HD.** Some pathways are detected only in control (CXCL, MHC-I, WNT, CD200) and others only in HD (SOMATOSTATIN, CCK, PACAP, KLK).
- **Links between the focal glial populations are kept in HD**. Microglial C3 complement autocrine signalling and APP–SORL1 signalling from Oligo_Oak to microglia (Mglia_Violet, Mglia_Daisy) are present at similar strength in both conditions (APP–SORL1 probability 0.194–0.195 in control, 0.207 in HD).
- **Supporting analyses:** pseudobulk DE shows the strongest transcriptional change in CN medium spiny neurons (162 DEGs). Glial proportions show trends but no significant HD vs control difference with 3 donors per group.
- **Astrocytic GJA1 gap-junction signalling** in Astro_Thyme is about 15% lower in HD (0.261 → 0.222), the only change between the focal populations.
---

## Pipeline at a glance

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1. QC and preprocessing | `01_preprocessing_{CN,CB,IFG,HIP}.R` | 24 GEO `.h5` files | 24 filtered `.rds` files |
| 2. Integration | `02_integration.R` | Filtered `.rds` files | `merged_filtered.rds` |
| 3. Broad annotation | `03_annotation.R` | `merged_filtered.rds` | `merged_annotated.rds` |
| 4. Subclustering | `04_subclustering.R` | `merged_annotated.rds` | `merged_subclustered_final.rds` |
| 5. Differential abundance | `differential_abundance.R` | `merged_subclustered_final.rds` | `results/DA/` |
| 6. Differential expression | `Differential_Expression_Pseudobulk.R` | `merged_subclustered_final.rds` | `results/DEG_tables/` |
| 7. TF activity (HSF1) | TF_Activity.R | `merged_subclustered_final.rds` |
| 8. Cell–cell communication | `06_cell_communication.R` | `merged_subclustered_final.rds` | CellChat objects, interaction tables |

**Steps 5–8 are independent analyses. Each starts from merged_subclustered_final.rds and none uses the output of another, so they can be run in any order.**

---

## The data

**Paper:** Bøstrand SMK, Seeker LA, Bestard-Cuche N, et al. *Mapping the glial transcriptome in Huntington's disease using snRNAseq: selective disruption of glial signatures across brain regions.* Acta Neuropathologica Communications 12, 165 (2024). [doi:10.1186/s40478-024-01871-3](https://doi.org/10.1186/s40478-024-01871-3) · [PubMed 39428482](https://pubmed.ncbi.nlm.nih.gov/39428482/)

**GEO series:** [GSE281069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281069)

**Technology:** Single-nucleus RNA-seq, 10x Genomics Chromium 3′ v3, sequenced on Illumina NovaSeq and aligned to GRCh38 with CellRanger v3.0.2 by the authors. We start from the authors' `filtered_feature_bc_matrix.h5` files.

**Samples: 3 HD and 3 control donors per region (24 samples)**

Donor IDs are taken from the GEO file names. Brain bank: NBB = Netherlands Brain Bank, Leiden = Leiden University Medical Center, EBB = Edinburgh Brain Bank.

### Caudate nucleus (CN)

| Sample   | Condition | GEO accession | Donor ID | Brain bank |
|----------|-----------|---------------|----------|------------|
| CN_HD1   | HD        | [GSM8610254](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610254) | HD240       | Leiden |
| CN_HD2   | HD        | [GSM8610272](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610272) | HD241       | Leiden |
| CN_HD3   | HD        | [GSM8610275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610275) | NBB17-060   | NBB |
| CN_ctrl1 | Control   | [GSM8610253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610253) | EBBSD025_13 | EBB |
| CN_ctrl2 | Control   | [GSM8610259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610259) | NBB17-005   | NBB |
| CN_ctrl3 | Control   | [GSM8610261](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610261) | NBB16-056   | NBB |

### Cerebellum (CB)

| Sample   | Condition | GEO accession | Donor ID | Brain bank |
|----------|-----------|---------------|----------|------------|
| CB_HD1   | HD        | [GSM8610257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610257) | NBB17-081 | NBB |
| CB_HD2   | HD        | [GSM8610258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610258) | HD241     | Leiden |
| CB_HD3   | HD        | [GSM8610266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610266) | HD240     | Leiden |
| CB_Ctrl1 | Control   | [GSM8610267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610267) | NBB95-310 | NBB |
| CB_Ctrl2 | Control   | [GSM8610283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610283) | NBB17-005 | NBB |
| CB_Ctrl3 | Control   | [GSM8610288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610288) | NBB16-056 | NBB |

### Inferior frontal gyrus / frontal cortex (IFG)

GEO labels some of these files `IFG` and others `FrCx`; the paper calls this region frontal cortex (FrCx).

| Sample    | Condition | GEO accession | Donor ID | Brain bank |
|-----------|-----------|---------------|----------|------------|
| IFG_HD1   | HD        | [GSM8610271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610271) | NBB17-060 | NBB |
| IFG_HD2   | HD        | [GSM8610277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610277) | HD240     | Leiden |
| IFG_HD3   | HD        | [GSM8610279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610279) | HD241     | Leiden |
| IFG_Ctrl1 | Control   | [GSM8610255](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610255) | NBB17-005 | NBB |
| IFG_Ctrl2 | Control   | [GSM8610262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610262) | NBB95-310 | NBB |
| IFG_Ctrl3 | Control   | [GSM8610263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610263) | NBB16-056 | NBB |

### Hippocampus (HIP)

| Sample    | Condition | GEO accession | Donor ID | Brain bank |
|-----------|-----------|---------------|----------|------------|
| HIP_HD1   | HD        | [GSM8610264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610264) | NBB17-060   | NBB |
| HIP_HD2   | HD        | [GSM8610268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610268) | HD231       | Leiden |
| HIP_HD3   | HD        | [GSM8610276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610276) | NBB17-081   | NBB |
| HIP_Ctrl1 | Control   | [GSM8610265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610265) | NBB16-056   | NBB |
| HIP_Ctrl2 | Control   | [GSM8610269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610269) | NBB17-005   | NBB |
| HIP_Ctrl3 | Control   | [GSM8610270](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8610270) | EBBSD025_13 | EBB |

> **CN sample labels.** In the final merged object, the CN_ctrl1 nuclei carry the sample label `CN_HD3`. The `CN_HD3` label covers 2,336 nuclei: 1,168 with condition `Control` and 1,168 with condition `HD`. Their cell barcodes still start with `CN_ctrl1_`, so the true sample can be recovered. Each downstream step handled this differently: the DA step restored the labels from the barcodes, and the DE step excluded `CN_HD3` entirely. See Limitations.

### How to download

1. Open each GSM link above and download the `*_filtered_feature_bc_matrix.h5` file from **Supplementary file(s)**.
2. Put the files in the matching region folder under `data/`, keeping the original GEO file names:

```
data/
├── CN/raw/    6 × CN  .h5 files
├── CB/raw/    6 × CB  .h5 files
├── IFG/raw/   6 × IFG/FrCx .h5 files
└── HIP/raw/   6 × HIP .h5 files
```

The `data/` folder and all `.rds` objects are git-ignored, so create the folder yourself. No data files are committed to this repository.

---

## How to run it

Open `single-cell-final-project.Rproj` so all paths resolve from the repository root, then run the scripts in order:

```r
# Step 1: one script per region (independent, any order)
source("scripts/01_preprocessing_CN.R")
source("scripts/01_preprocessing_CB.R")
source("scripts/01_preprocessing_IFG.R")
source("scripts/01_preprocessing_HIP.R")

# Steps 2 to 4: must run in this order
source("scripts/02_integration.R")
source("scripts/03_annotation.R")
source("scripts/04_subclustering.R")

# Steps 5 to 8: all start from merged_subclustered_final.rds (any order)
source("scripts/differential_abundance.R")
source("scripts/Differential_Expression_Pseudobulk.R")
source("scripts/TF_Activity.R")
source("scripts/06_cell_communication.R")
```

---

## Step 1: QC and preprocessing

Each of the 24 samples is processed on its own, with the same steps and thresholds in every region:

1. **Load:** read each `.h5` with `Read10X_h5()` and create a Seurat object. Compute `percent.mt` from genes starting with `MT-`.
2. **Filter nuclei:** keep nuclei with `200 < nFeature_RNA < 6000` and `percent.mt < 10`.
3. **Remove doublets:** convert to `SingleCellExperiment`, run `scDblFinder` with `set.seed(100)`, and keep nuclei classed as `singlet`.
4. **Normalize:** `NormalizeData()`, LogNormalize, scale factor 10,000.
5. **Highly variable genes:** `FindVariableFeatures()`, `vst`, 2,000 genes.
6. **Save:** one `.rds` per sample in `data/<REGION>/processed/`.

**Why these thresholds**

- **More than 200 genes:** removes empty droplets and debris.
- **Fewer than 6,000 genes:** removes the most extreme high-complexity barcodes, which are likely multiplets. scDblFinder handles the remaining doublets.
- **Mitochondrial reads below 10%:** isolating nuclei strips most mitochondria, so a high mito fraction points to cytoplasmic contamination or damaged nuclei. Most nuclei here are well below 5%.
- **One set of thresholds for all regions:** keeps the four regions directly comparable when they are merged. The QC violin plots of every sample were checked before filtering to confirm the cut-offs sit outside the main distribution.

| Script | Input | Output |
|--------|-------|--------|
| `01_preprocessing_CN.R`  | `data/CN/raw/*.h5`  | `data/CN/processed/CN_{HD1,HD2,HD3,ctrl1,ctrl2,ctrl3}.rds` |
| `01_preprocessing_CB.R`  | `data/CB/raw/*.h5`  | `data/CB/processed/CB_{HD1,HD2,HD3,Ctrl1,Ctrl2,Ctrl3}.rds` |
| `01_preprocessing_IFG.R` | `data/IFG/raw/*.h5` | `data/IFG/processed/IFG_{HD1,HD2,HD3,Ctrl1,Ctrl2,Ctrl3}.rds` |
| `01_preprocessing_HIP.R` | `data/HIP/raw/*.h5` | `data/HIP/processed/HIP_{HD1,HD2,HD3,Ctrl1,Ctrl2,Ctrl3}.rds` |


### Nuclei kept at each step

| Sample   | Raw  | After QC filter | Doublets removed | After doublet removal | % kept |
|----------|------|-----------------|------------------|-----------------------|--------|
| CN_HD1   | 3,440 | 2,529 | 138 | 2,391 | 69.5 |
| CN_HD2   | 1,337 | 1,098 | 42  | 1,056 | 79.0 |
| CN_HD3   | 1,528 | 1,438 | 65  | 1,373 | 89.9 |
| CN_ctrl1 | 2,431 | 2,194 | 104 | 2,090 | 86.0 |
| CN_ctrl2 | 1,093 | 1,007 | 42  | 965   | 88.3 |
| CN_ctrl3 | 5,621 | 4,343 | 290 | 4,053 | 72.1 |
| **CN total** | **15,450** | **12,609** | **681** | **11,928** | **77.2** |
| CB_HD1   | 2,980 | 2,914 | 159 | 2,755 | 92.4 |
| CB_HD2   | 1,985 | 1,838 | 87  | 1,751 | 88.2 |
| CB_HD3   | 7,425 | 6,747 | 688 | 6,059 | 81.6 |
| CB_Ctrl1 | 6,466 | 6,160 | 594 | 5,566 | 86.1 |
| CB_Ctrl2 | 4,533 | 4,443 | 349 | 4,094 | 90.3 |
| CB_Ctrl3 | 4,103 | 3,421 | 219 | 3,202 | 78.0 |
| **CB total** | **27,492** | **25,523** | **2,096** | **23,427** | **85.2** |
| IFG_HD1   | 1,328 | 1,279 | 54  | 1,225 | 92.2 |
| IFG_HD2   | 3,164 | 2,198 | 126 | 2,072 | 65.5 |
| IFG_HD3   | 8,377 | 7,390 | 732 | 6,658 | 79.5 |
| IFG_Ctrl1 | 2,593 | 2,340 | 120 | 2,220 | 85.6 |
| IFG_Ctrl2 | 2,280 | 862   | 30  | 832   | 36.5 |
| IFG_Ctrl3 | 3,846 | 3,354 | 216 | 3,138 | 81.6 |
| **IFG total** | **21,588** | **17,423** | **1,278** | **16,145** | **74.8** |
| HIP_HD1   | 1,161 | 995   | 52  | 943   | 81.2 |
| HIP_HD2   | 2,304 | 1,702 | 86  | 1,616 | 70.1 |
| HIP_HD3   | 4,070 | 3,889 | 250 | 3,639 | 89.4 |
| HIP_Ctrl1 | 4,430 | 3,966 | 258 | 3,708 | 83.7 |
| HIP_Ctrl2 | 4,929 | 4,578 | 338 | 4,240 | 86.0 |
| HIP_Ctrl3 | 810   | 761   | 30  | 731   | 90.2 |
| **HIP total** | **17,704** | **15,891** | **1,014** | **14,877** | **84.0** |
| **All regions** | **82,234** | **71,446** | **5,069** | **66,377** | **80.7** |

In total, 66,377 of 82,234 nuclei (80.7%) passed QC and doublet removal across the 24 samples; 61,558 remain after integration and cluster-level QC. The paper reports 127,205 nuclei from 44 samples after its own QC.

### QC metrics after filtering

All retained nuclei fall within 200–6,000 genes and below 10% mitochondrial reads. One representative sample per region is shown.

| **CN (CN_HD2)** | **CB (CB_HD1)** |
|---|---|
| ![QC, CN_HD2](<figures/Caudate nucleus/QC/CN_HD2-QC.png>) | ![QC, CB_HD1](<figures/Cerebellum/QC/CB_HD1_QC.png>) |
| **IFG (IFG_Ctrl1)** | **HIP (HIP_Ctrl1)** |
| ![QC, IFG_Ctrl1](<figures/Inferior frontal gyrus/QC/IFG_CRTL1_QC.png>) | ![QC, HIP_Ctrl1](<figures/hippocampus/QC/HIP_Ctrl1_QC.png>) |

**Doublet scores of retained singlets (CN_ctrl2 shown).** Most scores are close to zero after removal.

![Doublet scores, CN_ctrl2](<figures/Caudate nucleus/Doublets/CN-Ctrl2-Doublet.png>)

**Highly variable genes (CN_HD2 shown).** 2,000 variable genes (red) out of 33,538.

![Variable features, CN_HD2](<figures/Caudate nucleus/VariableFeatures/CN-HD2-HVG.png>)

Plots for all 24 samples are in [`figures/`](figures/), in one subfolder per region.

---

## Step 2: Integration

The filtered samples are merged, batch-corrected with Harmony, and clustered at resolution 0.25, the same resolution the paper used.

1. Load the filtered `.rds` files and add `condition` (HD / Control) and `region` (CB / CN / HIP / IFG) metadata.
2. Merge into one Seurat object, `JoinLayers()`, then `FindVariableFeatures()`, `ScaleData()` and `RunPCA()`.
3. Plot the UMAP **before** integration, coloured by sample, condition and region, to confirm the batch effect.
4. Run Harmony on `orig.ident` (one value per sample). We **never** correct for `condition`, which would remove the HD signal we want to study.
5. Plot the UMAP **after** integration to confirm samples now mix.
6. Cluster at resolutions 0.25, 0.3, 0.5 and 0.8, and keep 0.25.
7. Post-clustering QC: remove donor-dominated or very small clusters (criteria below), then save `merged_filtered.rds`.

### Key parameters

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `nfeatures` | 2,000 | Standard number of variable features |
| `npcs` | 50 computed, 24 used | ElbowPlot flattens after PC 20–25 |
| `dims` | 1:24 | From the ElbowPlot |
| Harmony `group.by.vars` | `orig.ident` | Corrects sample/donor batch, not condition |
| `resolution` | 0.25 | Matches the paper; chosen after comparing 0.3, 0.5 and 0.8 |

### Post-clustering QC

Following the paper, a cluster was removed if **one donor contributed more than 50% of it**, if **fewer than 5 donors each contributed at least 2%**, or if it had **fewer than 100 nuclei**.

| Cluster | Nuclei | Reason for removal |
|---------|--------|--------------------|
| 17 | 62          | 100% from CN_HD3, and under 100 nuclei |
| 9  | -           | 99.9% from CB_Ctrl1 |
| 11 | -           | 92.3% from CB_Ctrl1 |
| 12 | -           | 50.5% from one donor|

**Final dataset after QC:** 61,558 cells across 14 clusters (0–8, 10, 13–16)

### Key figures

**Before vs after integration, by condition.** Before Harmony, nuclei separate by sample. After Harmony, they group by cell type, and HD and control nuclei overlap.

![Before vs after integration, condition](<figures/Integration/UMAP Before After Integration -by Condition.jpg>)

**Before vs after integration, by region.** Some regional separation remains after integration. This is expected, because the four regions contain different cell types (for example, cerebellar granule cells).

![Before vs after integration, region](<figures/Integration/UMAP Before and After Intergation -by Region.jpg>)

**Final clusters at resolution 0.25, after post-clustering QC.**

![Final clusters](<figures/Integration/Clusters after Integration 0.25.jpg>)

<details>
<summary>All outputs of 02_integration.R</summary>

```
merged_filtered.rds
results/sessionInfo_integration.txt
figures/ElbowPlot.png
figures/Harmony_convergence.png
figures/UMAP_Before_by_{Condition,Region,Sample}.png
figures/UMAP_After_by_{Condition,Region,Sample}.png
figures/UMAP_Before_After_{Condition,Region,Sample}.png
figures/Resolution_comparison_0.3_0.5_0.8.png
figures/Clusters_Resolution_0.25.png
figures/Clusters_final_after_QC.png
```
</details>

---

## Step 3: Broad cell type annotation

Each of the 14 clusters is given a broad cell type label using the canonical marker panel from Bøstrand et al. 2024.

1. **DotPlot** of canonical markers across clusters, showing expression level and percentage. This is the primary evidence.
2. **FeaturePlot and VlnPlot** to confirm where each marker is expressed on the UMAP.
3. **Lineage-specific markers** for ambiguous clusters (clusters 6 and 7, and the vascular cluster 10).
4. **`FindAllMarkers()` (MAST)** as supporting evidence. These p-values are used only to rank markers for naming clusters, not to claim significance, because the same genes were used to build the clusters.
5. Save `merged_annotated.rds`.

### Annotation summary

| Cell type | Clusters | Nuclei | Key markers |
|-----------|----------|--------|-------------|
| Oligodendrocytes | 0 | 16,794 | MOG, MBP |
| Cerebellar granule cells | 1 | 14,038 | RELN, SLC17A7 |
| Excitatory neurons | 2, 13, 16 | 9,862 | SLC17A7 |
| Astrocytes | 3 | 4,736 | ALDH1L1, GFAP |
| Microglia | 4 | 3,781 | CX3CR1, ITGAM |
| Medium spiny neurons | 5 | 3,615 | PPP1R1B |
| Inhibitory neurons | 6, 8, 14 | 4,590 | GAD2 |
| OPCs | 7 | 2,048 | PDGFRA, GPR17 |
| Vascular cells | 10 | 1,803 | VWF, CLDN5 |
| Peripheral immune cells | 15 | 291 | CD2, CD8A |
| **Total** | **14 clusters** | **61,558** | |

### Key figures

**Canonical marker DotPlot**, the primary annotation evidence.

![DotPlot of canonical markers](<figures/Cell Annotation/DotPlot_canonical_markers.png>)

**Broad annotation UMAP**, with 10 cell types across 61,558 nuclei.

![Broad annotation UMAP](<figures/Cell Annotation/UMAP_broad_annotation_final.png>)

**Canonical marker FeaturePlots**, the secondary annotation evidence.

![FeaturePlot of canonical markers](<figures/Cell Annotation/FeaturePlot_canonical_markers.png>)

<details>
<summary>All outputs of 03_annotation.R</summary>

```
merged_annotated.rds
results/dotplot_data.csv
results/all_cluster_markers_annotated.csv
results/top5_cluster_markers_annotated.csv
results/sessionInfo_annotation.txt
figures/UMAP_cluster_numbers.png
figures/DotPlot_canonical_markers.png
figures/FeaturePlot_canonical_markers.png
figures/VlnPlot_key_markers.png
figures/VlnPlot_clusters_6_7.png
figures/VlnPlot_cluster_7_OPC_markers.png
figures/FeaturePlot_cluster_10_vascular.png
figures/UMAP_broad_annotation_final.png
```
</details>

---

## Step 4: Subclustering

Each broad population (microglia, astrocytes, oligodendroglia, neurons, vascular cells) is extracted and re-clustered on its own.

1. Re-run `FindVariableFeatures()`, `ScaleData()`, PCA, Harmony, UMAP and `FindClusters()` on each subset.
2. Choose the number of PCs per cell type from its ElbowPlot, and the resolution with the cleanest biological separation.
3. Remove contaminated clusters (neuronal doublets, cross-lineage contamination) and donor-dominated clusters.
4. Name subclusters from `FindAllMarkers()` (MAST) marker evidence. Glia follow the paper's naming convention: flowers for microglia, herbs for astrocytes, trees for oligodendroglia.
5. Write all labels back to the merged object. Removed nuclei are labelled `"Removed_contamination"`.
6. Save `merged_subclustered_final.rds` for downstream analysis.

> **Naming note:** the names follow the paper's convention, but a shared name does not always mean the same population as in the paper. Our Mglia_Rose, Mglia_Daisy, Astro_Thyme, Astro_Sage and Oligo_Birch have markers that match the paper's clusters of the same name. Our Mglia_Violet, Oligo_Oak and Oligo_Cedar have different marker profiles from the paper's clusters of those names.

### Microglia (flower names)

| Subcluster | Nuclei | Identity | Key markers |
|------------|--------|----------|-------------|
| Mglia_Rose   | 1,087 | Homeostatic | CX3CR1, P2RY12, TMEM119 |
| Mglia_Violet | 1,022 | Homeostatic variant | CX3CR1, moderate P2RY12 |
| Mglia_Lily   | 591   | DAM-like | APOE, TREM2 |
| Mglia_Daisy  | 527   | Disease-associated | LPL, CLEC7A |
| Mglia_BAMs   | 91    | Border-associated macrophages | LILRB4 |
| Mglia_Tulip  | 59    | Proliferating | MKI67 |

### Astrocytes (herb names)

| Subcluster | Nuclei | Identity | Key markers |
|------------|--------|----------|-------------|
| Astro_Basil | 1,409 | Reactive | CD44, APLNR |
| Astro_Mint  | 828   | Homeostatic | DIO2, GJB6 |
| Astro_Thyme | 765   | Stress response | HSPH1, HSPA4L |
| Astro_Sage  | 467   | Regional (cerebellar) | DAB2, PAX3 |

### Oligodendroglia (tree names)

| Subcluster | Nuclei | Identity | Key markers |
|------------|--------|----------|-------------|
| Oligo_Birch | 4,914 | Mature, myelinating | OPALIN, MBP, PLP1 |
| Oligo_Oak   | 4,749 | Stress response | HSPA1A, BAG3 |
| Oligo_Elder | 3,251 | Metabolically active | NEAT1, MT-ND1 |
| Oligo_Cedar | 3,413 | Committed OPCs (COPs) | ENPP6, DHCR24 |
| OPC_Ash     | 2,051 | OPCs | PDGFRA, CSPG4 |

### Neurons (descriptive and regional names)

| Subcluster | Nuclei | Identity | Main region |
|------------|--------|----------|-------------|
| Granule_CB       | 12,274 | Cerebellar granule cells | CB 96.5% |
| Excitatory_1     | 6,225  | Excitatory neurons | HIP + IFG |
| Excitatory_2_IFG | 2,634  | Excitatory neurons | IFG 59.3% |
| Inhibitory_SST   | 2,277  | SST+ inhibitory neurons | Mixed |
| Inhibitory_CB    | 1,707  | Cerebellar inhibitory neurons | CB 99.5% |
| Inhibitory_VIP   | 1,454  | VIP+ inhibitory neurons | Mixed |
| Excitatory_3     | 1,263  | Excitatory neurons | HIP + IFG |
| Excitatory_4     | 975    | Excitatory neurons | HIP + IFG |
| Cerebellar_CB    | 460    | Cerebellar neurons | CB 97.4% |
| Excitatory_5     | 460    | Excitatory neurons | HIP + IFG |
| Excitatory_HIP   | 421    | Excitatory neurons | HIP 64.4% |
| Cerebellar_2_CB  | 81     | Cerebellar neurons | CB 86.4% |



### Vascular cells (lineage and function)

| Subcluster | Nuclei | Identity | Key markers |
|------------|--------|----------|-------------|
| Pericyte      | 720 | Pericytes | TNC, CRB2, PDGFRB |
| Endothelial   | 406 | Endothelial cells | TIE1, VWF, CLDN5 |
| Smooth_Muscle | 360 | Smooth muscle cells | NOTCH3, CSPG4 |
| Endothelial_2 | 88  | Endothelial variant | CEMIP, C7 |

### Key figures

**All subclusters on the main integrated UMAP.** Nuclei removed during subclustering are shown as `Removed_contamination`.

![All subclusters](<figures/Subclustering/UMAP_all_subclusters.png>)

**Microglia subclusters**: homeostatic, disease-associated, border-associated and proliferating states.

![Microglia subclusters](<figures/Subclustering/UMAP_microglia_subclusters_final.png>)

**Neuron subclusters**: excitatory, inhibitory and cerebellar populations, with regional labels where one region dominates.

![Neuron subclusters](<figures/Subclustering/UMAP_neuron_subclusters.png>)

<details>
<summary>All outputs of 04_subclustering.R</summary>

```
merged_subclustered_final.rds
{microglia,astrocytes,oligodendroglia,neurons,vascular}_subclustered.rds
results/{microglia,astrocyte,oligodendroglia,neuron,vascular}_markers.csv
results/top5_{microglia,astrocyte,oligodendroglia,neuron,vascular}_markers.csv
results/sessionInfo_subclustering.txt
figures/ElbowPlot_{microglia,astrocytes,oligodendroglia,neurons,vascular}.png
figures/{Microglia,Astrocytes,Oligodendroglia,Neurons,Vascular}_resolution_comparison.png
figures/UMAP_{microglia,astrocyte,oligodendroglia,neuron,vascular}_subclusters_final.png
figures/UMAP_all_subclusters.png
```
</details>

### Notes for downstream analysis

- **Drop removed nuclei first.** Exclude them before any downstream analysis:
  ```r
  merged_clean <- subset(merged, subset = subcluster != "Removed_contamination")
  ```
- **Check the CN sample labels** (see "CN sample labels" under The data). Barcodes starting with `CN_ctrl1_` belong to CN_ctrl1, even where `orig.ident` says `CN_HD3`.
- **Cerebellar oligodendroglia:** as in the paper, exclude CB oligodendroglia from subcluster DE because there are very few of them.
- **Compare donors, not nuclei.** HD vs control tests should use sample-level values (proportions or pseudobulk), not individual nuclei.

---

## Step 5: Differential abundance

**Question for this step:** do the proportions of the four main glial types change between HD and control in each region?

1. Restore the CN_ctrl1 sample and condition labels from the `CN_ctrl1_` barcode prefix, and re-check the sample × condition tables.
2. Count nuclei per sample and broad cell type, and divide by each sample's total to get proportions. Each sample's proportions sum to 1.
3. Keep oligodendrocytes, astrocytes, microglia and OPCs.
4. For each region × cell type (4 × 4 = 16 tests), compare the 3 control and 3 HD sample proportions with a two-sided Wilcoxon rank-sum test.
5. Correct for multiple testing with Benjamini–Hochberg FDR, and report the HD minus control difference in medians.

**Result:** None of the 16 region × cell type comparisons is significant after FDR correction (smallest p = 0.10, smallest FDR = 0.533). With 3 samples per group, a Wilcoxon test cannot give a p-value below 0.10, so these are descriptive trends only.

- **Cerebellum:** all four glial types are lower in HD, e.g. astrocytes 10.4% → 3.4% and oligodendrocytes 6.7% → 2.0%. The oligodendrocyte gap rests largely on one control sample with about 50% oligodendrocytes.
- **Caudate and frontal cortex:** oligodendrocytes are higher in HD (CN 39.9% → 53.2%; IFG 24.5% → 31.1%).
- **Hippocampus and frontal cortex:** microglia are lower in HD (HIP 11.0% → 5.0%; IFG 6.8% → 3.9%).

![Glial proportions, HD vs control](<results/figures/DA/DA_glial_proportions.png>)

**Conclusion:** Glial abundance does not clearly change with HD at this sample size. Because proportions are shares of all nuclei, a higher oligodendrocyte share in HD caudate may partly reflect neuronal loss rather than more glia. The HD effects in this dataset show up more in gene expression and signalling than in cell numbers.
**Outputs:** `results/DA/DA_final_results.csv`, `results/figures/DA/DA_glial_proportions.png`

---

## Step 6: Differential expression (pseudobulk)

**Question for this step:** which genes change between HD and control within each cell type and region, and is the signal strongest in the CN and its medium spiny neurons?

1. **Pseudobulk:** sum raw counts per sample × broad cell type with `AggregateExpression()`. The sample, not the nucleus, is the replicate.
2. **Filtering:** drop sample × cell type groups with fewer than 10 nuclei. Test a cell type only if it has at least 2 pseudobulk samples per condition. Drop genes with fewer than 10 total reads.
3. **DESeq2**, design `~ condition`, contrast HD vs control, run separately for every region × cell type.
4. **Significance:** `padj < 0.05`, with no log2 fold-change cut-off.

**Replicates per region:** CB, HIP and IFG use 3 vs 3. CN uses 2 vs 2 (CN_ctrl2, CN_ctrl3 vs CN_HD1, CN_HD2), because `CN_HD3` was excluded for its mixed labels.

**Result:** the DE signal is concentrated in CN neurons. The other regions show few or no DEGs at this sample size.

| Region | Cell type | DEGs | Up | Down |
|--------|-----------|-----:|---:|-----:|
| CN  | Medium spiny neurons     | 162 | 146 | 16 |
| CN  | Excitatory neurons       | 96  | 75  | 21 |
| CB  | Cerebellar granule cells | 29  | 10  | 19 |
| CN  | Oligodendrocytes         | 22  | 4   | 18 |
| IFG | Microglia                | 14  | 5   | 9  |
| CN  | Astrocytes               | 12  | 8   | 4  |
| CN  | Microglia                | 8   | 2   | 6  |
| IFG | Medium spiny neurons     | 3   | 1   | 2  |
| IFG | Vascular cells           | 3   | 0   | 3  |
| IFG | Astrocytes               | 2   | 1   | 1  |
| HIP | Excitatory neurons       | 1   | 1   | 0  |
| CN  | OPCs                     | 1   | 0   | 1  |

All other region × cell type pairs had 0 significant DEGs. In total 353 genes are significant, and 301 of them (85%) are in the caudate. The hippocampus has a single DEG.

![Volcano plot, CN medium spiny neurons](<figures/DE_Pseudobulk/CN/Volcano_Medium_Spiny_Neurons_HD_vs_Control.png>)

In caudate medium spiny neurons, 146 of 162 significant genes are up in HD, many with log2 fold changes of 2–10. Such large changes suggest genes that are near zero in control, and with only 2 samples per group in CN they need confirming.

**Conclusion:** Transcriptional change is concentrated in the caudate, the region HD damages first, and mainly in its neurons (medium spiny and excitatory). Glia in the other regions change little at this sample size. Because the CN comparison uses only 2 vs 2 samples, treat it as provisional; a missing DEG elsewhere may reflect low power rather than no biological change.


**Outputs:** `results/DEG_tables/<REGION>/` (one CSV per cell type: `gene`, `log2FoldChange`, `pvalue`, `padj`), `results/DEG_tables/Summary_all_regions.csv`, volcano plots and heatmaps in `results/figures/<REGION>/`, and gene lists for pathway analysis in `results/gene_lists_for_pathway/`.

---

## Step 7: Transcription factor activity (HSF1)

**Question for this step:** is HSF1, the master regulator of the heat-shock response, more active in the stress-response glial subclusters?

1. Load the human **DoRothEA** TF–target network and keep interactions with confidence A, B or C.
2. Infer TF activity per nucleus from the normalized RNA data with the **ULM** method in `decoupleR` (minimum 5 targets per TF), using the 2,000 variable genes.
3. Store the scores as a new `TF_activity` assay.
4. Plot HSF1 **activity** and HSF1 **gene expression** on the same UMAP. Activity reflects HSF1's target genes, which can differ from the level of HSF1 mRNA itself.
5. Average TF activity per broad cell type, and show the 10 most variable TFs in a heatmap.

**Result:** Across broad cell types, inferred HSF1 activity is highest in peripheral immune cells and moderately positive in microglia, oligodendrocytes and cerebellar granule cells. It is close to zero or below in the other neuronal types. The other top variable TFs follow each lineage as expected: SOX10 is strongly active in oligodendrocytes and OPCs, and immune regulators (SPI1, RUNX1, RFX5, LYL1) are most active in microglia. This shows the inference behaves sensibly.
**Conclusion** HSF1, the master regulator of the heat-shock response, is active in glia, especially microglia and oligodendrocytes, more than in most neurons. This fits the glial stress-response and chaperone signature described in the paper, and the stress-response subclusters (Astro_Thyme, Oligo_Oak) we focus on. Activity was averaged by cell type, not split by condition, so this analysis does not show whether HSF1 activity is higher in HD than in control. Comparing HSF1 activity between HD and control within the stress-response subclusters is the next step.

| HSF1 TF activity | HSF1 gene expression |
|---|---|
| ![HSF1 activity](<figures/TF Activity and Expression/TF_Activity.png>) | ![HSF1 expression](<figures/TF Activity and Expression/TF_Gene Expression.png>) |

![Top 10 variable TFs across subclusters](<figures/TF Activity and Expression/Heatmap_most expressed TF.png>)

**Outputs:** `figures/TF Activity and Expression/TF_Activity.png`, `figures/TF Activity and Expression/TF_Gene Expression.png`, `figures/TF Activity and Expression/Heatmap_most expressed TF.png`

---

## Step 8: Cell–cell communication (main question)

**Populations analysed:**

| Subcluster | Why it was chosen |
|------------|-------------------|
| Astro_Thyme  | Stress-response astrocytes (HSPH1, HSPA4L); matches the paper's HD-enriched Astro Thyme |
| Oligo_Oak    | Stress-response oligodendrocytes (HSPA1A, BAG3) |
| Mglia_Violet | Microglia; strongest source of autocrine C3 signalling |
| Mglia_Daisy  | Disease-associated microglia (TREM2, CLEC7A, LPL) |

We did not test whether these subclusters are enriched in HD ourselves (Step 5 tested broad cell types only). The microglial populations were chosen from marker evidence, without a separate check against heat-shock genes.

**Method:**

1. Start from `merged_subclustered_final.rds`. Drop nuclei still carrying a broad label from the contamination-removal step (for example `Astrocytes` or `Microglia`), but keep `Peripheral Immune Cells`, which were never subclustered.
2. Build separate **CellChat** objects for HD and control, grouped by `subcluster`.
3. Run the standard pipeline: `subsetData` → `identifyOverExpressedGenes` → `identifyOverExpressedInteractions` → `computeCommunProb` → `filterCommunication(min.cells = 10)` → `computeCommunProbPathway` → `aggregateNet`.
4. Merge the two objects with `mergeCellChat()` and compare HD with control directly.

### Results

**1. Pathways are gained and lost in HD.**
- Detected only in control: CXCL, MHC-I, WNT, CD200 and others.
- Detected only in HD: SOMATOSTATIN, CCK, PACAP, KLK.

![Pathways ranked by information flow, HD vs control](<figures/Cell to Cell comunication/rankNet_pathways.png>)

**2. Microglial C3 complement signalling is present in both conditions.**
Mglia_Violet and Mglia_Daisy signal to themselves through `C3 → C3AR1` and `C3 → ITGAM+ITGB2 / ITGAX+ITGB2`, at similar strength in HD and control. This relates to the paper's discussion of complement-mediated synapse engulfment by microglia, but we see no HD-specific increase.

| HD | Control |
|---|---|
| ![Complement signalling, HD](<figures/Cell to Cell comunication/complement_HD.png>) | ![Complement signalling, control](<figures/Cell to Cell comunication/complement_Control.png>) |

**3. APP–SORL1 signalling from oligodendrocytes to microglia is present in both conditions.**
Oligo_Oak → Mglia_Daisy and Oligo_Oak → Mglia_Violet via APP–SORL1 have a communication probability of 0.194–0.195 in control and 0.207 in HD. This is a stable link between the focal populations, not an HD-specific change.

**4. Astrocyte gap-junction signalling is reduced in HD.**
GJA1–GJA1 autocrine signalling in Astro_Thyme is present in both conditions but lower in HD (communication probability 0.261 in control, 0.222 in HD, about 15% lower). It is reduced, not lost, and CellChat gives no HD vs control test for it, so it is a lead to validate rather than a confirmed result.

![Top 25 interactions, HD vs control](<figures/Cell to Cell comunication/bubble_top25.png>)

### Answer to our question

Communication is rewired in HD rather than uniformly increased or shut down: some pathways are detected only in control (CXCL, MHC-I, WNT, CD200) and others only in HD (SOMATOSTATIN, CCK, PACAP, KLK). Between our four focal populations, however, the main links are kept: both microglial subclusters share the same C3 autocrine programme, and Oligo_Oak signals to both through APP–SORL1, in HD and control alike. This is partial evidence of a shared, coordinated glial network, but because it is equally present in control it does not show coordination caused by HD. The only change we see between the focal populations is a modest drop in GJA1–GJA1 gap-junction signalling within Astro_Thyme (0.261 → 0.222, about 15%).

**Outputs:** CellChat objects (`.rds`) and `session info` in `results/`; `results/interactions_HD_enriched_populations.csv`, `results/interactions_Control_enriched_populations.csv`; rankNet, complement circle plots, overall circle plots, differential-interaction plot and bubble plot in `figures/`. The bubble plot shows only the 25 strongest interactions; the full tables are in `results/`.

---

## Requirements

Full `sessionInfo()` outputs are in `results/`.

**Step 1: QC and preprocessing**

| Package              | CN     | CB     | IFG    | HIP    |
|----------------------|--------|--------|--------|--------|
| R                    | 4.6.1  | 4.6.1 | 4.6.1  | 4.6.1  |
| Seurat               | 5.5.1  | 5.5.1 | 5.5.1  | 5.5.1  |
| SeuratObject         | 5.4.0  | 5.4.0 | 5.4.0  | 5.4.0  |
| SingleCellExperiment | 1.34.0 | 1.34.0 | 1.34.0 | 1.34.0 |
| scDblFinder          | 1.26.7 | 1.26.7 | 1.26.7 | 1.26.7 |
| hdf5r                | 1.3.12 | 1.3.12 | 1.3.15 | 1.3.12 |

The random seed is `set.seed(100)`, set once before the six `scDblFinder()` calls in each script. Keep the calls in the same order to reproduce the doublet calls exactly.

**Steps 2–4: integration, annotation, subclustering**

| Package      | Integration | Annotation | Subclustering |
|--------------|-------------|------------|---------------|
| R            | 4.6.1       | 4.6.1      | 4.6.1  |
| Seurat       | 5.5.1       | 5.5.1      | 5.5.1  |
| SeuratObject | 5.4.0       | 5.4.0      | 5.4.0  |
| harmony      | 2.0.5       | –          | 2.0.5  |
| MAST         | –           | 1.38.0     | 1.38.0 |
| ggplot2      | 4.0.3       | 4.0.3      | 4.0.3  |
| dplyr        | 1.2.1       | 1.2.1      | 1.2.1  |
| patchwork    | 1.3.2       | 1.3.2      | 1.3.2  |

**Steps 5–8: downstream analysis**

The cell–cell communication step was run with R 4.6.1, Seurat 5.5.1, SeuratObject 5.4.0, ggplot2 4.0.3, dplyr 1.2.1, patchwork 1.3.2 and igraph 2.3.3.

Install MAST, DESeq2, decoupleR and dorothea with `BiocManager::install()`, harmony with `install.packages("harmony")`, and CellChat with `devtools::install_github("jinworks/CellChat")`.

---

## Differences from the paper

| Step              | Bøstrand et al. 2024                                | This analysis |
|-------------------|-----------------------------------------------------|---------------|
| Input counts      | CellRanger + Velocyto (spliced and unspliced reads) | CellRanger filtered matrix from GEO only |
| Samples           | 44 samples from 12 donors (one control excluded)    | 24 samples; CN_ctrl1 merged under the CN_HD3 label |
| QC tool           | Scater                                              | Seurat |
| Gene filter       | Genes expressed in ≤200 nuclei removed              | All genes kept |
| Nuclei filter     | Per-sample UMI, gene and mito limits (Table S2)     | Fixed: 200–6,000 genes, mito < 10% |
| Doublets          | scDblFinder, removed if score ≥ 0.94                | scDblFinder default class (stricter) |
| Normalization     | scran, per brain region                             | LogNormalize, per sample |
| Batch correction  | Seurat v4 CCA between 10x chips, 37 PCs             | Harmony by sample, 24 PCs |
| Clustering        | Louvain, resolution 0.25                            | Louvain, resolution 0.25 (same) |
| Cluster QC        | Donor-contribution and size filters                 | Same criteria |
| Marker genes      | MAST                                                | MAST |
| Subcluster names  | Flowers / herbs / trees                             | Same convention, but not always the same populations |
| Differential abundance | MiloR on KNN neighbourhoods, subcluster level, all regions combined | Sample-level proportions of broad glial types, Wilcoxon per region |
| Differential expression | MAST on single nuclei, \|log2FC\| ≥ 0.8, padj < 0.05 | DESeq2 pseudobulk per sample, padj < 0.05 |
| TF activity       | Not done                                            | DoRothEA + decoupleR ULM (HSF1) |
| Cell–cell communication | Not done                                      | CellChat, HD vs control |

Because of these differences, our nucleus counts, subclusters and statistics will not match the paper's exactly.

---

## Limitations

- **CN sample labels.** In the merged object, the CN_ctrl1 nuclei are labelled `CN_HD3`, so `CN_HD3` holds 2,336 nuclei split exactly 1,168 Control / 1,168 HD. The barcode prefix still identifies the true sample:
  - Harmony and the donor-share checks for cluster removal treated the two samples as one (for example cluster 17).
  - DA restored the labels from the barcodes and uses 3 vs 3 in CN.
  - DE excluded `CN_HD3` and uses 2 vs 2 in CN.
  - Cell–cell communication groups nuclei by condition only.

  CN results should be read with caution.
- One fixed set of QC thresholds is used for every sample and region, whereas the paper tuned them per sample. The cerebellum, dominated by small granule cells, may need a lower gene-count floor.
- Only 3 donors per condition per region, and not always the same donors across regions. With 3 vs 3, a Wilcoxon test cannot reach p < 0.10, so DA results are descriptive only.
- Two of the removed clusters (9 and 11) came almost entirely from CB_Ctrl1, so that sample contributes less to the final data than the others.
- The paper added intronic (unspliced) reads with Velocyto; the GEO matrices we use may not include them, so genes and UMIs per nucleus can be lower than in the paper.
- No ambient RNA correction (for example SoupX or CellBender) was applied.
- DE uses `padj < 0.05` with no fold-change cut-off, so some significant genes have small effects. The absence of DEGs in a region may reflect low power rather than no biological change.
- TF activity is inferred from a generic, non-brain-specific network, using only the 2,000 variable genes.
- CellChat probabilities reflect ligand–receptor co-expression, not measured interactions. "Not detected" means below CellChat's threshold, not absent. The communication findings are hypotheses to validate.

---

## Team

| Member | Role | GitHub |
|--------|------|--------|
| Shayma Alsweis | Caudate nucleus (CN) preprocessing, Cell–cell communication| [@shaymasweis-94](https://github.com/shaymasweis-94) |
| Ali | Cerebellum (CB) preprocessing, Transcription factor activity | [@alinawar99](https://github.com/alinawar99) |
| Menna | Inferior frontal gyrus (IFG) preprocessing, Differential expression (pseudobulk) | [@Mennae4](https://github.com/Mennae4) |
| Menna Tawfik | Hippocampus (HIP) preprocessing, Differential abundance| [@mennamahmoudtawfik281-cmd](https://github.com/mennamahmoudtawfik281-cmd) |
| Eman Abadi | Integration, annotation, subclustering,  | [@ejabadi](https://github.com/ejabadi) |

**Note on commit authorship:** Commits `397db93` ("Add CN preprocessing script and processed data") and `57e9b90` ("Add cell-cell communication script, figures, and result summaries") appear under the GitHub username *shaymahaji*. This happened because the local Git email was mistyped as shayma@gmail.com instead of shayma.sweis@gmail.com, and GitHub links commits to accounts by email. These commits were authored by Shayma Alsweis ([@shaymasweis-94](https://github.com/shaymasweis-94)). The history was left unchanged to preserve the original commit times; the configuration was corrected on 29 Sep 2026.

Repository: [mennamahmoudtawfik281-cmd/single-cell-final-project](https://github.com/mennamahmoudtawfik281-cmd/single-cell-final-project)

