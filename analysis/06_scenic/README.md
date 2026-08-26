# SCENIC -- Gene Regulatory Network Inference

Gene regulatory network analysis using [pySCENIC](https://github.com/aertslab/pySCENIC), including sex-stratified regulon activity comparison in Alzheimer's disease.

## Scripts

| Script | Description |
|--------|-------------|
| `scenic_ctx_aucell_viz.py` | Full pySCENIC pipeline: ctx (regulon prediction), AUCell (cellular enrichment), UMAP, heatmaps, and Z-score clustermaps for six sex-diagnosis groups. |
| `scenic_regulon.py` | Regulon Z-score computation and UMAP visualization for microglia. |
| `sex_scenic_heatmaps_networks.py` | Male vs Female comparison: per-cell-type regulon heatmaps, differential regulon analysis (M-F), and TF-TG networkx plots. |
| `sex_scenic_igraph_networks.py` | TF-TG heatmaps and igraph network visualizations with sex-specific edge identification and optional GWAS gene highlighting. |

## Inputs

- pySCENIC adjacency CSVs and loom files
- Ranking database: `hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- Motif annotations: `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`
- GWAS gene list: `gwas_genes_ad_ftd.xlsx` (optional)
- Sex-stratified h5ad files

## Outputs

- Regulon CSVs and AUCell matrices per group
- Regulon Z-score heatmaps and UMAP plots
- TF-TG network plots (per-sex and comparison)
- Sex-specific edge tables
