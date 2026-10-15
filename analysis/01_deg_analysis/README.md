# Differential Expression Analysis

Multiple DE frameworks applied across AD subtypes, MCI, and controls at both single-cell and pseudobulk resolution.

## Scripts

### Python (Wilcoxon / DESeq2)

| Script | Description |
|--------|-------------|
| `deg_all.py` | Wilcoxon one-vs-REST across all condition groups (all cells). |
| `deg_all_dx.py` | Pairwise and one-vs-REST Wilcoxon across all diagnosis groups (1--4, MCI, Control). |
| `deg_ad1-4.py` | Pairwise and one-vs-REST Wilcoxon restricted to AD subtypes 1--4 only. |
| `per_celltype_deg.py` | Per-cell-type Wilcoxon (all groups: 1--4, MCI, Control). |
| `per_celltype_deg_1-4.py` | Per-cell-type Wilcoxon restricted to AD subtypes 1--4. |
| `pseudobulk_allcells_DEG.py` | Pseudobulk DESeq2 analysis aggregating all cells per sample. |

### R (MAST)

| Script | Description |
|--------|-------------|
| `run_mast.R` | MAST DE for a single cell-type h5ad (AD vs Control, MCI vs Control, AD vs MCI). |
| `subtype_deg.R` | Loops over all h5ad files in a directory, running MAST per cell subtype. |
| `subtype_deg_array.R` | SLURM array-compatible version of subtype MAST DE. |
| `subcelltype_run_findmarkers.R` | SLURM array-compatible MAST for sub-cell-type level analysis. |

## Inputs

- `data/h5ad/13AD_matched_label_final.h5ad` -- main annotated snRNA-seq object
- `data/csv/20_seeds_final_subgroups_with_meta.csv` -- subtype metadata
