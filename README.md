# Single-Nucleus Transcriptomic Analysis and Molecular Subtyping of Alzheimer’s Disease

Analysis workflows and scripts supporting the panAD single-nucleus transcriptomic atlas and MONET molecular subtyping framework, including differential expression, co-expression networks, cell–cell communication, and gene regulatory network analyses.

> **Note:** This repository will be made public upon publication.

---

## Repository Structure

```
.
├── analysis/
│   ├── 01_monet_subtyping/     # Masked VAE + GMM ensemble subtyping (MONET)
│   ├── 02_deg_analysis/        # Differential expression (Wilcoxon, MAST, DESeq2)
│   ├── 03_deg_summary/         # DEG visualization and cross-method comparison
│   ├── 04_cellchat/            # Cell-cell communication (CellChat)
│   ├── 05_hdwgcna_modules/     # Weighted gene co-expression networks (hdWGCNA)
│   └── 06_scenic/              # Gene regulatory network inference (pySCENIC)
├── data/                       # Input data (not tracked; see data/README.md)
├── results/                    # Output directory for generated figures and tables
├── environment.yml             # Conda environment specification
├── requirements.txt            # Python dependencies (pip)
└── README.md
```

## Analysis Modules

### 1. MONET Subtyping (`01_monet_subtyping/`)
Masked variational autoencoder + Gaussian mixture model ensemble framework ("MONET") that resolves four reproducible AD molecular subtypes from covariate-adjusted, multi-cell-type pseudobulk expression. Trains across 20 random seeds and derives a consensus clustering from the resulting pairwise sample co-assignment matrix. 

### 2. Differential Expression (`02_deg_analysis/`)
Differential gene expression across AD subtypes (1--4), MCI, and controls using multiple frameworks:
- **Wilcoxon rank-sum** (scanpy) -- all-cells and per-cell-type analyses
- **MAST** (R/Seurat) -- single-cell-aware hurdle model, SLURM-compatible
- **DESeq2** (via rpy2) -- pseudobulk approach aggregating counts per sample

### 3. DEG Summary (`03_deg_summary/`)
Cross-method comparison (DESeq2 vs MAST scatter plots), multi-cell-type DEG heatmaps, and overlap with hdWGCNA hub genes.

### 4. CellChat (`04_cellchat/`)
Sex- and diagnosis-stratified cell-cell communication inference using CellChat on scANVI-predicted cell types. Separate analyses for Control, MCI, and Dementia groups by sex.

### 5. hdWGCNA Modules (`05_hdwgcna_modules/`)
Donor-level weighted gene co-expression network analysis across seven major cell types. Includes module-trait (Braak, CERAD, MMSE) Spearman correlations with FDR correction, differentially expressed module (DME) lollipop plots, and publication-ready multi-panel figures.

### 6. SCENIC (`06_scenic/`)
Gene regulatory network inference using pySCENIC. Includes sex-stratified regulon activity analysis, TF-target gene networks (igraph), and differential regulon heatmaps between male and female AD samples.

## Setup

### Option A: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate negin-ad
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
                    "tidyr", "glue", "scales", "Matrix", "future", "reticulate"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("MAST", "WGCNA"))

# CellChat
devtools::install_github("jinworks/CellChat")

# hdWGCNA
devtools::install_github("smorabit/hdWGCNA", ref = "dev")
```

## Data

Input data files (`.h5ad`, `.rds`, `.csv`) are not included due to size and access restrictions. See [`data/README.md`](data/README.md) for the expected directory layout and file descriptions.

## Reproducing the Analyses

Each module under `analysis/` is self-contained. Scripts are designed to be run from the **repository root** with data files placed in `data/` as described above.

```bash
# Example: run the all-cells DEG analysis
python analysis/02_deg_analysis/deg_all.py

# Example: run CellChat for control females
Rscript analysis/04_cellchat/cellchat_ct_fe.R

# Example: run MAST DEG (takes h5ad path as argument)
Rscript analysis/02_deg_analysis/run_mast.R data/h5ad/celltype_file.h5ad
```

> **Note:** Some scripts assume HPC execution (SLURM). See inline comments in `subtype_deg_array.R` and `subcelltype_run_findmarkers.R`.

## Citation

If you use this code, please cite the associated publication (details to be added upon publication).

## License

All rights reserved. This code is shared for reproducibility of the associated publication.
