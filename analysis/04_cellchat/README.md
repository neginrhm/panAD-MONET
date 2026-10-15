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
| `cellchat_dem_male.R` | Dementia | Male |

All scripts follow the same pipeline: load h5ad via reticulate, build CellChat object, identify overexpressed genes/interactions, compute communication probabilities, and generate circle plots.

`cellchat_function.r` factors the shared pipeline logic above into a single reusable function (h5ad in, CellChat RDS + standard plots out) and is the basis the per-condition/sex scripts follow.

### Differential analysis (pairwise disease-stage comparisons)

| Script | Comparison |
|--------|-----------|
| `cellchat_differential.r` | Dementia vs. Control |
| `cellchat_differential_ct_mci.r` | MCI vs. Control |
| `cellchat_differential_mci_ad.r` | Dementia vs. MCI |
| `cellchat_differential_ad4_vs_ct.r` | AD Subtype 4 vs. Control |

These load the pre-computed per-condition CellChat RDS objects (not sourced as code — `cellchat_function.r` is run first to produce them), merge them, and generate differential interaction plots (circle, heatmap, chord, bubble), signaling role analysis, and per-pathway visualizations. This is the analysis underlying the manuscript's Fig. 4 differential interaction heatmaps and chord diagrams.

### AD-subtype comparison (Control vs. AD1-4)

| Script | Role |
|--------|------|
| `run_cellchat_pair_subtype.R` | Generic pairwise comparison script (same structure as `run_cellchat_pair_sex.R`) -- run once per subtype to produce the AD{n}-vs-Control comparison suite, including a full per-pathway sweep covering ANGPTL, RA, and TGF-beta. |
| `cellchat_subtype_panel_fig8b.R` | Combines the 5 per-group CellChat objects (Control, AD1-4) into the 4 Control-vs-AD{n} merged objects and produces the unified rankNet panel with one consistent Control+AD1-4 palette -- the multi-subtype comparison behind Fig. 8b / Suppl. Fig. 6a. |

This is the analysis underlying the manuscript's subtype-stratified CellChat findings (Fig. 8b: ANGPTL/RA/TGF-beta present in Control, AD3 = near-complete loss of all three, AD2 = selective ANGPTL loss, AD1/AD4 approx. Control).

`cellchat_differential_ad4_vs_ct.r` above covers the AD4-vs-Control contrast specifically with the same full-plot-suite treatment as the other differential scripts; `run_cellchat_pair_subtype.R` + `cellchat_subtype_panel_fig8b.R` cover all four subtypes together via the generic pathway/rankNet route.

### Sex comparison (male vs. female, per diagnosis group)

| Script | Comparison |
|--------|-----------|
| `run_cellchat_pair_sex.R` | Male vs. Female (run once each for Control, MCI, Dementia) |

A generic, parameterized script (`--female_rds`, `--male_rds`, `--outdir`, etc.) that loads a female/male CellChat RDS pair for one diagnosis group, rebuilds inference if needed, merges them (`object.list <- list(female = ..., male = ...)`), and generates the full comparison suite: `compareInteractions`, differential interaction circle and heatmap plots, per-dataset circles, signaling-role scatter/heatmaps, rankNet info-flow bar charts, bubble plots, functional/structural similarity embeddings, and per-pathway circle/heatmap/chord plots. This is the analysis underlying the manuscript's sex-stratified CellChat supplementary figures.

Run once per diagnosis group, e.g.:
```
Rscript run_cellchat_pair_sex.R \
  --female_rds cellchat_by_ct_fe.rds --male_rds cellchat_by_ct_male.rds \
  --outdir cellchat_diff_results_CT_female_vs_male
```

**Color convention:** because `object.list` is built as `list(female = ..., male = ...)`, CellChat's standard diverging red/blue outputs (`netVisual_diffInteraction`, `netVisual_heatmap`) follow its dataset-order convention: **red = higher in male** (dataset 2), **blue = higher in female** (dataset 1). `rankNet` instead uses an explicit categorical palette (`female = "#CC79A7"` pink, `male = "#0072B2"` blue) that identifies which sex a bar belongs to, not direction of change.

## Inputs

- `data/h5ad/13adata_{condition}_{sex}.h5ad` -- condition- and sex-specific AnnData subsets
- Differential scripts additionally require the `.rds` CellChat objects produced by the per-condition/sex scripts above

## Outputs

- `cellchat_by_{condition}_{sex}.rds` -- serialized CellChat object
- `celltype_global_circles{condition}_{sex}.pdf` -- communication circle plot
- Differential scripts: comparison-specific PDF plots and `.RData` merged-object files (see each script's header for exact filenames)
