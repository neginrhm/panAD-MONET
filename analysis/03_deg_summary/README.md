# DEG Summary and Visualization

Post-processing, cross-method comparison, and publication-quality visualization of differential expression results.

## Scripts

| Script | Description |
|--------|-------------|
| `extract_deg_5celltypes.py` | Extracts DEGs across five cell types (EXC, ODC, INH, AST, MIC), overlaps with hdWGCNA hub genes, and produces summary heatmaps. |
| `scatter_with_gene_labels.py` | Annotated scatter plots comparing DESeq2 vs MAST log2 fold-change estimates per contrast and cell type. CLI-based with `--deseq_dir`, `--mast_dir`, `--outdir` arguments. |
| `markers.py` | Canonical and subtype-specific marker gene sets for major brain cell types; generates dot plots. |
