# Data

This directory is a placeholder for input data files required by the analysis scripts.

Due to size and privacy constraints, data files are **not** included in this repository. To reproduce the analyses, place the following files in the appropriate subdirectories:

## Expected structure

```
data/
├── h5ad/
│   ├── 12AD_matched_label_final.h5ad        # Main annotated snRNA-seq object
│   ├── 12adata_ct_female.h5ad                # Control female subset
│   ├── 12adata_ct_male.h5ad                  # Control male subset
│   ├── 12adata_mci_female.h5ad               # MCI female subset
│   ├── 12adata_mci_male.h5ad                 # MCI male subset
│   ├── 12adata_dem_female.h5ad               # Dementia female subset
│   └── 12adata_dem_male.h5ad                 # Dementia male subset
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
    └── gwas_genes_ad_ftd.xlsx
```

## Data flow for `analysis/01_monet_subtyping/`

This module has its own multi-stage data flow, run from the repository root
with per-cell-type raw h5ad files (e.g. `12ad_Microglia.h5ad`, one per major
cell type) as the starting input:

```
12ad_{CellType}.h5ad                                  (raw, per cell type -- input to 00a)
    -> 00a_h5ad_to_seurat.r  -> 12ad_{celltype}_hdWGCNA_object.rds
    -> 00_pseudobulk.r       -> pseudobulk_{CellType}.csv
    -> 01_regress_covariates.r -> Regressed_data_DEM_{CT}_thisexp~Study+Age+Sex.csv
    -> 02_build_feature_vectors.py -> input_dnn_celltypes_pseudobulked_mapped_VectorsofFloats.csv
    -> 03_vae_ensemble_pipeline.py  -> final_assignments.csv, per-seed cluster CSVs
```

Place the raw per-cell-type h5ad files in `data/h5ad/` (or update the paths
in `00a_h5ad_to_seurat.r`); intermediate CSVs are produced by each script in
turn and consumed by the next.

## Data access

Contact the corresponding author for data availability.
