# panAD-MONET

Analysis code accompanying **"An Integrated Single-Nucleus Atlas Resolves Cell-Type-Specific Programs and Molecular Subtypes in Alzheimer's Disease"** (Rahimzadeh et al.) -- a 13-cohort, >3,000,000-nucleus snRNA-seq atlas (791 donors, Control/MCI/AD) used to resolve four reproducible AD molecular subtypes (the "MONET" framework) and their cell-type-specific transcriptional programs, gene regulatory networks, and cell-cell communication changes.

The processed atlas, metadata, DE results, and GRN edges are browsable at the [Stellar Atlas](https://swaruplab.bio.uci.edu/panad_atlas/) (installable package: [swaruplab/stellar-atlas](https://github.com/swaruplab/stellar-atlas) -- a separate repository from this one).

---

## Repository Structure

Modules are numbered in the order their results first appear in the manuscript (Fig. 1 -> Fig. 8; Fig. 9 lives in a separate repository).

```
.
├── analysis/
│   ├── 00_preprocessing_qc/     # scVI/scANVI (Allen reference) -> Leiden -> annotation -> ENDO rescue (Fig. 1)
│   ├── 01_deg_analysis/         # Differential expression: MAST, Wilcoxon, DESeq2
│   ├── 02_deg_summary/          # DEG visualization and cross-method comparison
│   ├── 03_hdwgcna_modules/      # Weighted gene co-expression networks (hdWGCNA) (Fig. 2-3)
│   ├── 04_cellchat/             # Cell-cell communication (CellChat) (Fig. 4, 8b)
│   ├── 05_scenic/               # Gene regulatory networks: pySCENIC + NetREm (Fig. 5-6)
│   ├── 06_sex_differences/      # Sex-stratified DEG x hdWGCNA module enrichment (Fig. 6a-b)
│   ├── 07_monet_subtyping/      # Masked VAE + GMM ensemble subtyping, CMap, clinical stats (Fig. 7-8)
│   └── utils/                   # Shared helpers (ORA.R)
├── results/                     # Output directory for generated figures and tables
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies (pip)
└── README.md
```

## Analysis Modules

### 0. Preprocessing, Integration & Annotation (`00_preprocessing_qc/`) -- Fig. 1
Per-cohort processing (CellBender, Scanpy QC, Scrublet) is not included. This module picks up after each of the 13 cohorts has been independently processed. It integrates them by scVI/scANVI reference mapping onto the Allen Human Multiple Cortical Areas SMART-seq reference, with each sample as a batch, and clusters the scANVI latent space with Leiden (resolution 0.3; 17 clusters). The saved models from the original run can be loaded to reproduce the published embedding exactly; a second mode reruns the mapping following the current scvi-tools tutorials. Major cell types are assigned where the scANVI label transfer and marker-based annotation of the Leiden clusters agree (script 00b). A small endothelial population missed by the initial annotation is then added as its own cluster. The substate subclustering was done interactively and is described, not scripted, in the module README. Also includes the Fig. 1b UMAP and pie-chart scripts and the Suppl. Fig. 1 QC figure.

### 1. Differential Expression (`01_deg_analysis/`)
Differential gene expression across AD subtypes (1--4), MCI, and controls using multiple frameworks:
- **MAST** (R/Seurat) -- the method reported in Methods ("Differential expression analysis") (`min.pct=0, logfc.threshold=0`), for the 3-way disease-stage comparison (`run_mast.R`) and the 35-test subtype battery (`subtype_deg*.R`, `subcelltype_run_findmarkers.R`)
- **Wilcoxon rank-sum** (scanpy) and **DESeq2** (via rpy2) -- earlier/cross-check approaches, not the manuscript's reported method


### 2. DEG Summary (`02_deg_summary/`)
Cross-method comparison (DESeq2 vs MAST scatter plots), multi-cell-type DEG heatmaps, overlap with hdWGCNA hub genes, and the subtype DEG UpSet plots (Fig. 7d).

### 3. hdWGCNA Modules (`03_hdwgcna_modules/`) -- Fig. 2-3
Per-cell-type network construction (`00_hdwgcna_construct_network.r`: metacells, signed networks, module eigengenes, hub genes, GO enrichment), reported for five cell types (EXC, INH, MIC, ODC, ASC). Module-trait (Braak, CERAD, MMSE) Spearman correlations with BH-FDR correction are computed in **ROSMAP samples only (n=176)**, as disclosed in the manuscript. DME lollipop plots (Figs. 2a, 3a). Module numbers (e.g. MIC has M1,M2,M4,M5,M10 but no M3/M6-M9) come from a shared WGCNA color-to-number lookup applied identically across cell types -- MIC/ODC simply produced fewer distinct modules than ASC/EXC/INH, nothing is filtered or hidden.

### 4. CellChat (`04_cellchat/`) -- Fig. 4, 8b
Sex- and diagnosis-stratified cell-cell communication inference using CellChat v2 on scANVI-predicted cell types, plus a subtype-stratified panel (Control vs. AD1-4) covering the ANGPTL/RA/TGF-beta signaling loss described in Fig. 8b.

### 5. SCENIC + NetREm (`05_scenic/`) -- Fig. 5-6
Gene regulatory network inference using pySCENIC (GRNBoost2 -> cisTarget -> AUCell; per cell type and sex-stratified run scripts included), refined with [NetREm](https://github.com/SaniyaKhullar/NetREm)'s TF-TF cooperative/antagonistic classification (Fig. 5), plus sex-stratified regulon activity analysis and TF-target networks (Fig. 6c-d).

> **NetREm's regression implementation is not included in this repository.** It was performed by a collaborator (S. Khullar) using the separately-published NetREm method; this repo's Fig. 5 scripts consume NetREm's output TF-to-TG tables as external input data (see that module's README).

### 6. Sex Differences (`06_sex_differences/`) -- Fig. 6a-b
Module-level enrichment of AD-only sex-differential MAST DEGs (male vs. female) within hdWGCNA co-expression modules, with a statistical annotation gate (FDR<0.05, overlap q>=4, module overlap>=2%).

### 7. MONET Subtyping (`07_monet_subtyping/`) -- Fig. 7-8
Masked variational autoencoder + Gaussian mixture model ensemble framework ("MONET") that resolves four reproducible AD molecular subtypes (AD1=Metal-Ion Stress, AD2=Neuroinflammatory, AD3=Synaptic Integrity, AD4=Tissue Remodeling) from covariate-adjusted, multi-cell-type pseudobulk expression. Trains across 20 random seeds and derives a consensus clustering from the resulting pairwise sample co-assignment matrix. Also includes the per-seed silhouette diagnostic, CLUE/CMap gene lists, result parsing and heatmaps, Braak/CERAD/MMSE clinical-trait comparisons, and subtype module-eigengene plots (Fig. 8a). An earlier exploratory approach (1D convolutional autoencoder + Harmony + UMAP + DBSCAN), not used for any results in the manuscript, is not included.

**Limitations stated in the manuscript:** silhouette analysis (k=2-10) favored k=2-3; k=4 was chosen to resolve finer-grained heterogeneity. Subtyping is performed at the **library/sample level, not the donor level**: of the 391 subtyped samples, 21 (5.4%) belong to donors whose duplicate samples received different subtype calls (13 of 36 SEA-AD donors with replicate libraries and 8 of 12 cross-dataset donors).

## Manuscript Figure Map

| Figure | Analysis | Module |
|--------|----------|--------|
| Fig. 1 | Atlas construction, QC, integration, annotation; UMAP and pie charts | `00_preprocessing_qc/` |
| Fig. 2-3 | hdWGCNA co-expression modules (EXC/INH, MIC/ODC/ASC) | `03_hdwgcna_modules/` |
| Fig. 4 | CellChat: Control/MCI/AD differential signaling | `04_cellchat/` |
| Fig. 5 | pySCENIC + NetREm gene regulatory networks | `05_scenic/` |
| Fig. 6 | Sex-stratified DE, regulons, CellChat | `06_sex_differences/`, `05_scenic/`, `04_cellchat/` |
| Fig. 7 | MONET subtyping, clinical-trait comparison; subtype DEG UpSet plots (7d) | `07_monet_subtyping/`, `02_deg_summary/` |
| Fig. 8 | Subtype module eigengenes (8a), subtype CellChat (8b), CMap/CLUE (8c) | `07_monet_subtyping/`, `04_cellchat/` |
| Fig. 9 | Stellar Atlas | separate repository ([swaruplab/stellar-atlas](https://github.com/swaruplab/stellar-atlas)) |
| Suppl. Fig. 1 | QC metrics by study | `00_preprocessing_qc/supp_fig1_qc_violins.py` |
| Suppl. Table 1 | Dataset/cohort metadata | -- (see manuscript) |
| Suppl. Table 2-3 | Module characterization, hub genes | `03_hdwgcna_modules/` |
| Suppl. Table 4 | CMap candidate compounds | `07_monet_subtyping/` (`04b_parse_clue_gct.R`, `05_cmap_heatmap.py`) |

## Setup

### Option A: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate panad-monet
```

### Option B: pip (Python dependencies only)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

R packages (CellChat, Seurat, MAST, hdWGCNA, WGCNA) must be installed separately. See the [R package installation notes](#r-dependencies) below.

### R Dependencies

Install from an R session:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "cowplot", "forcats",
                    "tidyr", "tidyverse", "glue", "scales", "Matrix", "future",
                    "reticulate", "optparse", "patchwork", "tibble", "rlang",
                    "boot", "viridis", "RColorBrewer", "circlize", "ggpubr",
                    "ggsignif", "tidytext", "hdf5r", "harmony", "anndata", "enrichR"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("MAST", "WGCNA", "ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"))

# GitHub-only packages
devtools::install_github("jinworks/CellChat")
devtools::install_github("smorabit/hdWGCNA", ref = "dev")
devtools::install_github("satijalab/seurat-data")     # SeuratData
devtools::install_github("mojaveazure/seurat-disk")   # SeuratDisk
devtools::install_github("carmonalab/UCell")          # UCell
devtools::install_github("cmap/cmapR")                # cmapR (04b_parse_clue_gct.R)
```

## Data

Input data files (`.h5ad`, `.rds`, `.csv`) are not included due to size and data-use restrictions.

- **Raw snRNA-seq data:** all 13 datasets are publicly available from GEO and Synapse under their respective data use agreements. Supplementary Table 1 of the manuscript lists each study with its PubMed ID and repository accession.
- **Integrated atlas, processed data, metadata and differential expression results:** https://swaruplab.bio.uci.edu/panad_atlas/

Scripts that read from a `data/` directory expect this layout (create it at the repository root):

```
data/
├── h5ad/
│   ├── 13AD_matched_label_final.h5ad        # Main annotated snRNA-seq object
│   ├── 13adata_ct_female.h5ad                # Control female subset
│   ├── 13adata_ct_male.h5ad                  # Control male subset
│   ├── 13adata_mci_female.h5ad               # MCI female subset
│   ├── 13adata_mci_male.h5ad                 # MCI male subset
│   ├── 13adata_dem_female.h5ad               # Dementia female subset
│   └── 13adata_dem_male.h5ad                 # Dementia male subset
├── csv/
│   ├── 20_seeds_final_subgroups_with_meta.csv
│   └── rosmap_pathology_cluster_ctrl_mci_ad.csv
├── CellTypes/
│   ├── Microglia/hdWGCNA/...                 # hdWGCNA Seurat objects (.rds)
│   ├── Inhibitory/...
│   ├── Excitatory/...
│   ├── ODC/...
│   ├── Astrocyte/...
│   ├── opc/...
│   └── Endo/...
└── scenic/
    ├── hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    ├── motifs-v9-nr.hgnc-m0.001-o0.0.tbl
    ├── gwas_genes_ad_ftd.xlsx
    └── netrem_tables/
        ├── NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_Ast_Nov_1_2025_bySaniya_Final.csv
        ├── ..._ExN_..., ..._InN_..., ..._Mic_..., ..._Oli_...  (one per cell type)
        └── TF_TF_BRAIN_{AST,EXN,INH,MIC,OLI}.csv  (NetREm TF-TF cooperative/antagonistic tables)
```

The NetREm TF-TG tables under `data/scenic/netrem_tables/` are NetREm outputs produced externally by a collaborator and are not generated by any script here (see `analysis/05_scenic/README.md`).

## Reproducing the Analyses

Each module under `analysis/` is self-contained. Scripts read their inputs by file name, either from `data/` or from the working directory; each module README lists the expected inputs and the run order.

```bash
# Example: run the all-cells DEG analysis
python analysis/01_deg_analysis/deg_all.py

# Example: run CellChat for control females
Rscript analysis/04_cellchat/cellchat_ct_fe.R

# Example: run MAST DEG (takes h5ad path as argument)
Rscript analysis/01_deg_analysis/run_mast.R data/h5ad/celltype_file.h5ad
```

> **Note:** Some scripts assume HPC execution (SLURM). See inline comments in `subtype_deg_array.R` and `subcelltype_run_findmarkers.R`.

## Known Limitations & Gaps

This repository does not provide full raw-to-figure reproducibility. Not included or not scripted:

1. **Per-cohort processing (Fig. 1).** CellBender, Scanpy QC and Scrublet were run per cohort before integration; those per-cohort runs are not included.
2. **Atlas annotation (Fig. 1).** Integration (scVI/scANVI reference mapping, Leiden) and the major cell-type annotation (label transfer + Leiden markers, concordant cells kept) are scripted in `00_preprocessing_qc/`. The substate subclustering was done interactively and is described, not scripted.
3. **hdWGCNA settings for EXC, INH, OPC and ENDO.** These networks were built with the same template as MIC/ODC/ASC, but their original soft-power and module-detection settings were not preserved; `00_hdwgcna_construct_network.r` uses hdWGCNA defaults for them.
4. **NetREm (Fig. 5).** Run by S. Khullar with the separately published [NetREm](https://github.com/SaniyaKhullar/NetREm) method; the Fig. 5 scripts here read NetREm's output tables as input data.
5. **Per-sample regulon activity vs. disease traits.** The code correlating sample-averaged AUCell scores with disease traits is not included.
6. **Neff et al. MSBB subtype concordance (Fig. 7 / Discussion).** The code for this comparison is not included.

Notes on interpretation:

- **pySCENIC networks are not seeded.** GRNBoost2 is stochastic; reruns reproduce regulons approximately.
- **Sex-stratified adjacency filter.** `05_scenic/01_filter_adjacency_noncoding.py` also removes some protein-coding target genes (symbols ending in P+digits); see that script's docstring.
- **k=4 was not the silhouette-optimal solution** (silhouette favored k=2-3 across seeds; recorded by `03_vae_ensemble_pipeline.py`); k=4 was chosen to resolve finer-grained heterogeneity.
- **Two different clinical sample sets.** hdWGCNA module-trait correlations (Figs. 2-3) use ROSMAP samples only (n=176); the MONET clinical-trait comparison (Fig. 7c) pools ROSMAP + SEA-AD (n=217) with CERAD rescaled between the two cohorts' scales.

## Citation

If you use this code, please cite:

> Negin Rahimzadeh, Samuel Morabito, Saniya Khullar, Zechuan Shi, Zhenkun Cao, Vivek Swarup. **An Integrated Single-Nucleus Atlas Resolves Cell-Type-Specific Programs and Molecular Subtypes in Alzheimer's Disease.** *bioRxiv* 2026.09.03.747935; doi: [https://doi.org/10.64898/2026.09.03.747935](https://doi.org/10.64898/2026.09.03.747935)

```bibtex
@article{Rahimzadeh2026panAD,
  title   = {An Integrated Single-Nucleus Atlas Resolves Cell-Type-Specific Programs and Molecular Subtypes in Alzheimer's Disease},
  author  = {Rahimzadeh, Negin and Morabito, Samuel and Khullar, Saniya and Shi, Zechuan and Cao, Zhenkun and Swarup, Vivek},
  journal = {bioRxiv},
  year    = {2026},
  note    = {2026.09.03.747935},
  doi     = {10.64898/2026.09.03.747935}
}
```

## License

Released under the [MIT License](LICENSE).
