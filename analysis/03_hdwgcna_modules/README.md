# hdWGCNA Module Analysis

Weighted gene co-expression network analysis using [hdWGCNA](https://github.com/smorabit/hdWGCNA) across seven major brain cell types, with module-trait correlation and differentially expressed module (DME) analyses. Figs. 2-3 report five of them (EXC, INH, MIC, ODC, ASC); the OPC and endothelial networks feed the MONET pseudobulk input (`07_monet_subtyping/00_pseudobulk.r`).

## Scripts

| Script | Description |
|--------|-------------|
| `00_hdwgcna_construct_network.r` | Network construction per cell type: Seurat object from `13ad_{CellType}.h5ad`, Harmony (Study), metacells (k=25, max_shared=10, per Celltype_markers x Sample), signed network, module eigengenes (harmonized by Sample), kME, top-10 hub genes (Suppl. Table 3), enrichR GO enrichment. Run once per cell type: `Rscript 00_hdwgcna_construct_network.r Microglia`. Network settings recorded per cell type in `CT_SETTINGS` (MIC/ODC/AST from the original runs; EXC/INH/OPC/ENDO use hdWGCNA defaults because their original settings were not preserved). |
| `hdwgcna_module_trait_correlation.r` | Module-trait correlation heatmap for Inhibitory neurons (Braak, CERAD, MMSE). |
| `signed_cor_hdwgcna_modules_celltypes.r` | Donor-level (ROSMAP: one sample per donor) Spearman module-trait correlations across all seven cell types with BH-FDR correction. |
| `signed_cor_hdwgcna_modules_LAYOUT3.r` | Publication-ready grouped panels (Neurons, Glia, VascOPC) with mean and median pseudobulking. |
| `hdwgcna_DME_lollipop.r` | Differentially expressed modules between diagnosis groups (hdWGCNA `FindDMEs`, `test.use = 'wilcox'`); lollipop plots of log2FC across cell types (Figs. 2a, 3a). |

## Inputs

- `data/CellTypes/{celltype}/hdWGCNA/*.rds` -- hdWGCNA Seurat objects per cell type
- `data/csv/rosmap_pathology_cluster_ctrl_mci_ad.csv` -- donor clinical metadata

## Cell types analyzed

Excitatory, Inhibitory, Astrocyte, Microglia, ODC, OPC, Endothelial
