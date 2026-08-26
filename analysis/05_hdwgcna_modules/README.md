# hdWGCNA Module Analysis

Weighted gene co-expression network analysis using [hdWGCNA](https://github.com/smorabit/hdWGCNA) across seven major brain cell types, with module-trait correlation and differentially expressed module (DME) analyses.

## Scripts

| Script | Description |
|--------|-------------|
| `hdwgcna_module_trait_correlation.r` | Module-trait correlation heatmap for Inhibitory neurons (Braak, CERAD, MMSE). |
| `signed_cor_hdwgcna_modules_celltypes.r` | Donor-level Spearman module-trait correlations across all seven cell types with BH-FDR correction. |
| `signed_cor_hdwgcna_modules_LAYOUT3.r` | Publication-ready grouped panels (Neurons, Glia, VascOPC) with mean and median pseudobulking. |
| `hdwgcna_DME_lollipop.r` | Differentially expressed modules between diagnosis groups; lollipop plots of log2FC across cell types. |

## Inputs

- `data/CellTypes/{celltype}/hdWGCNA/*.rds` -- hdWGCNA Seurat objects per cell type
- `data/csv/rosmap_pathology_cluster_ctrl_mci_ad.csv` -- donor clinical metadata

## Cell types analyzed

Excitatory, Inhibitory, Astrocyte, Microglia, ODC, OPC, Endothelial
