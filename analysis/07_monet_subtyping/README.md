# AD Molecular Subtyping (MONET)

Masked variational autoencoder (VAE) + Gaussian mixture model (GMM) ensemble
framework ("MONET": Multi-seed Optimization of Neural Embeddings for
subTyping) that resolves four reproducible AD molecular subtypes (AD1-AD4)
from covariate-adjusted, multi-cell-type pseudobulk expression.

Full pipeline, in order:

| Script | Description |
|--------|-------------|
| `00_pseudobulk.r` | Aggregates raw counts into sample-level pseudobulk profiles per cell type, from the per-cell-type hdWGCNA Seurat objects built by `../03_hdwgcna_modules/00_hdwgcna_construct_network.r`. |
| `01_regress_covariates.r` | Bootstrap linear regression (R=1000, `thisexp ~ Study + Age + Sex`) to remove technical/demographic covariates; median coefficient across bootstraps used for adjustment. |
| `02_build_feature_vectors.py` | Assembles a 7-element `[AST, MIC, ENDO, ODC, OPC, INH, EXC]` covariate-adjusted expression vector per gene per sample. |
| `03_vae_ensemble_pipeline.py` | Core MONET pipeline: masked VAE (128-dim latent) trained across 20 random seeds -> GMM (k=4) clustering per seed -> pairwise co-assignment matrix -> agglomerative consensus clustering -> Sankey diagrams. Also records GMM silhouette scores for k=2-10 per seed (`silhouette_by_seed_k.csv`, `silhouette_by_k.pdf`). Includes `load_data()` for parsing/imputing/padding the raw feature-vector CSV. |
| `04_qc_consensus.py` | Co-assignment heatmap QC, ordered by final consensus cluster. |
| `04a_make_clue_gene_lists.py` | Writes the top-150 up/down gene lists per subtype (from Wilcoxon subtype-vs-rest DEG tables, adjusted p < 0.05, ranked by log fold change) that were uploaded to the CLUE L1000 query tool. |
| `04b_parse_clue_gct.R` | Parses a raw CLUE/CMap query result (.gct) for one subtype signature into the filtered top-compound CSV consumed by `05_cmap_heatmap.py`. Run once per subtype (S1-S4). |
| `05_cmap_heatmap.py` | CLUE/CMap drug-repurposing connectivity heatmap per subtype (Fig. 8c / Suppl. Table 4), hierarchically clustered, vector (editable-text) output, full top-30-per-signature compound set. **The published figure panel was further hand-curated in Adobe Illustrator afterward** (a subset of compound labels removed for clarity) -- that step is editorial, not scripted, and is not reproduced here. |
| `06_clinical_boxplots.py` | Braak stage / CERAD score / MMSE boxplots by subtype only (AD1-4). |
| `07_clinical_stats.py` | Braak/CERAD/MMSE comparison across **Control, MCI, and AD1-4** (Fig. 7c / Suppl. Fig. 6b) -- sample-level Kruskal-Wallis + pairwise Mann-Whitney U with BH-FDR. Produces the fuller sample-level table `06_clinical_boxplots.py` does not include Control/MCI in. |
| `08_subtype_module_eigengene_boxplots.py` | hdWGCNA module eigengenes per sample across Control, MCI and AD1-4, all modules per cell type (Fig. 8a shows a curated subset of modules). |

## Method summary

For each sample, cell-type-specific pseudobulk profiles are covariate-adjusted
and concatenated into a single multi-cell-type feature vector, then encoded
by a masked VAE into a 128-dimensional latent space. This is repeated across
20 independently trained random seeds; a Gaussian mixture model (k=4) is
fitted to each seed's latent embedding, and a consensus clustering is derived
from the resulting pairwise sample co-assignment matrix (agglomerative
clustering, average linkage). See the manuscript Methods for the full
mathematical description (masked reconstruction + KL loss, co-assignment and
consensus-distance formulas).

Note: silhouette analysis across the 20 seeds favored k=2 or k=3
(`03_vae_ensemble_pipeline.py` records it); k = 4 was chosen to resolve
finer-grained molecular heterogeneity, as stated in the manuscript.

Subtyping is performed per sample (library), not per donor: donors with
several libraries, or present in more than one cohort, contribute several
samples (see the manuscript Methods).

## Earlier exploratory approach (not included)

An earlier deep-learning trial (1D convolutional autoencoder + Harmony + UMAP
+ DBSCAN) was tried before settling on the masked-VAE + GMM ensemble approach
above. It was not used for any results in the manuscript and is not included
in this repository.

## Data flow

This module has its own multi-stage data flow, run from the repository root
with per-cell-type raw h5ad files (e.g. `13ad_Microglia.h5ad`, one per major
cell type) as the starting input:

```
13ad_{CellType}.h5ad                                  (raw, per cell type)
    -> 03_hdwgcna_modules/00_hdwgcna_construct_network.r -> 13ad_{celltype}_hdWGCNA_object.rds
    -> 00_pseudobulk.r       -> pseudobulk_{CellType}.csv
    -> 01_regress_covariates.r -> Regressed_data_DEM_{CT}_thisexp~Study+Age+Sex.csv
    -> 02_build_feature_vectors.py -> input_dnn_celltypes_pseudobulked_mapped_VectorsofFloats.csv
    -> 03_vae_ensemble_pipeline.py  -> final_assignments.csv, per-seed cluster CSVs
```

Place the raw per-cell-type h5ad files in the working directory of each
script; intermediate CSVs are produced by each script in
turn and consumed by the next.

### Drug repurposing (`04b_parse_clue_gct.R` -> `05_cmap_heatmap.py`)

```
DEG_*_vs_REST.csv subtype DEG tables
    -> 04a_make_clue_gene_lists.py -> top-150 up/down gene lists per subtype
    -> submitted to CLUE (https://clue.io) as an L1000 query, one per subtype
    -> CLUE sig_queryl1k_tool result directory (arfs/TAG/query_result.gct), one per subtype
    -> 04b_parse_clue_gct.R (run once per subtype S1-S4) -> data/cmap/s{N}_cmap_top_compounds.csv
    -> 05_cmap_heatmap.py -> figures/cmap_heatmap_clean_colors.png/.pdf/.svg
    -> (published Fig. 8c) hand-curated in Adobe Illustrator: a subset of
       compound labels removed from the full heatmap above for clarity --
       not a scripted step
```

`data/cmap/` is expected to contain:
```
data/cmap/
├── s1_nothreshold_cmap_top_compounds.csv   # AD1: no compounds passed the significance filter;
│                                             unfiltered top hits shown for comparison (Supplementary Table 4)
├── s2_cmap_top_compounds.csv               # AD2
├── s3_cmap_top_compounds.csv               # AD3
└── s4_cmap_top_compounds.csv               # AD4
```

Raw CLUE `.gct` query results are not included (regenerate via the CLUE
website/API from the subtype DEG signatures); the parsed/filtered CSVs above
are also not committed as they are generated outputs.
