# Sex-Difference Analysis (MAST DE x hdWGCNA modules)

Module-level enrichment of AD-only sex-differential genes (male vs. female, MAST) within hdWGCNA co-expression modules, matching the manuscript's Fig. 6a-b.

## Scripts

| Script | Role |
|--------|------|
| `01_sex_deg_pipeline.r` | Steps 1-3: removes X/Y-chromosome genes from the per-cell-type male-vs-female MAST DEG tables, filters to `p_val_adj < 0.05`, Venn-compares Dementia vs. Control to isolate AD-only sex DEGs, and runs a Module-Eigengene (ME) Wilcoxon test by sex. Also produces a first-draft ORA enrichment heatmap (its own Step 4) and a compact ME-difference bar chart (Step 5) — superseded for the published heatmap by the script below, but kept as this pipeline's own diagnostic output. |
| `02_sex_deg_ora_heatmap_v4_gated.r` | The publication ORA heatmap. Reads the AD-only sex DEGs written by `01_sex_deg_pipeline.r` and re-runs the module enrichment against the **full** hdWGCNA module gene sets (hundreds of genes/module, not just hub genes), with a `MIN_DEG >= 10` power filter and a hard annotation gate (stars/OR/% labels render only for `FDR < 0.05 AND overlap q >= 4 AND module overlap >= 2%`); failing cells keep their color but no text. Produces the combined supplementary table (`sex_ORA_supplementary_table_v4.csv`) with a `passes_heatmap_gate` flag per (cell type, direction, module). |

Run `01_sex_deg_pipeline.r` first; `02_sex_deg_ora_heatmap_v4_gated.r` depends on its `ad_only_degs/*.csv` output.

The module-number relabeling in `02_sex_deg_ora_heatmap_v4_gated.r` (`module_rename`, e.g. `brown -> ASC-M1`) is the same WGCNA color-to-number lookup used throughout the hdWGCNA modules -- see [`../03_hdwgcna_modules/README.md`](../03_hdwgcna_modules/README.md) for why MIC/ODC module numbers aren't contiguous (they simply produced fewer distinct modules than ASC/EXC/INH, so only the numbers matching colors that actually appeared are used).

## Inputs

- `sex_DEG_analysis/MAST_results/{Dementia,Control}/{celltype}_male_vs_female_normalized.csv` -- per-cell-type, per-diagnosis MAST male-vs-female DE results
- `ipynb/not_ipynb/13ad_metadata.csv` -- sample-level Sex/Diagnosis metadata
- `CellTypes/{CellType}/[hdWGCNA/]{celltype}_MEs_by_sample_manual.csv` -- module eigengenes by sample
- `CellTypes/{CellType}/[hdWGCNA/]13ad_{celltype}_modules_table.csv` -- full hdWGCNA module gene-to-color assignments (used by `02_...v4_gated.r`)
- `CellTypes/{CellType}/[hdWGCNA/]13ad_{celltype}_modules_hub_df_table.csv` -- hub-gene-only module tables (used by `01_sex_deg_pipeline.r`'s own first-draft ORA)
- `../utils/ORA.R` -- shared Fisher's-exact over-representation helper (both scripts `source()` this)

## Outputs

- `sex_DEG_analysis/pipeline_v2/` -- cleaned DEGs, Venn diagrams, AD-only/shared/control-only DEG CSVs, ME-by-sex boxplots and Wilcoxon results, first-draft ORA heatmaps, bar charts
- `sex_DEG_analysis/module_enrichment_v4_gated/` -- publication ORA heatmaps (per cell type + combined), supplementary table, color-scale summary
