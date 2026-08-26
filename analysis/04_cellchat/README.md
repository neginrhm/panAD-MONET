# CellChat Analysis

Sex- and diagnosis-stratified cell-cell communication inference using [CellChat](https://github.com/jinworks/CellChat) on scANVI-predicted cell types, plus differential interaction analysis across disease stages (Control, MCI, Dementia) matching the manuscript's Fig. 4.

## Scripts

### Per-condition/sex CellChat objects

| Script | Condition | Sex |
|--------|-----------|-----|
| `cellchat_ct_fe.R` | Control | Female |
| `cellchat_ct_male.R` | Control | Male |
| `cellchat_mci_fe.R` | MCI | Female |
| `cellchat_mci_male.R` | MCI | Male |
| `cellchat_dem_fe.R` | Dementia | Female |

All scripts follow the same pipeline: load h5ad via reticulate, build CellChat object, identify overexpressed genes/interactions, compute communication probabilities, and generate circle plots.

> **Note:** there is no Dementia/Male script in this set — only Female is present for the Dementia condition. Worth confirming whether that combination was intentionally excluded or is simply missing from this copy.

`cellchat_function.r` factors the shared pipeline logic above into a single reusable function (h5ad in, CellChat RDS + standard plots out) and is the basis the per-condition/sex scripts follow.

### Differential analysis (pairwise disease-stage comparisons)

| Script | Comparison |
|--------|-----------|
| `cellchat_differential.r` | Dementia vs. Control |
| `cellchat_differential_ct_mci.r` | MCI vs. Control |
| `cellchat_differential_mci_ad.r` | Dementia vs. MCI |
| `cellchat_differential_sexes.r` | Male vs. Female |

These load the pre-computed per-condition CellChat RDS objects (not sourced as code — `cellchat_function.r` is run first to produce them), merge them, and generate differential interaction plots (circle, heatmap, chord, bubble), signaling role analysis, and per-pathway visualizations. This is the analysis underlying the manuscript's Fig. 4 differential interaction heatmaps and chord diagrams.

## Inputs

- `data/h5ad/12adata_{condition}_{sex}.h5ad` -- condition- and sex-specific AnnData subsets
- Differential scripts additionally require the `.rds` CellChat objects produced by the per-condition/sex scripts above

## Outputs

- `cellchat_by_{condition}_{sex}.rds` -- serialized CellChat object
- `celltype_global_circles{condition}_{sex}.pdf` -- communication circle plot
- Differential scripts: comparison-specific PDF plots and `.RData` merged-object files (see each script's header for exact filenames)
