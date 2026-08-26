# AD Molecular Subtyping (MONET)

Masked variational autoencoder (VAE) + Gaussian mixture model (GMM) ensemble
framework ("MONET": Multi-seed Optimization of Neural Embeddings for
subTyping) that resolves four reproducible AD molecular subtypes (AD1-AD4)
from covariate-adjusted, multi-cell-type pseudobulk expression.

Full pipeline, in order:

| Script | Description |
|--------|-------------|
| `00a_h5ad_to_seurat.r` | Reads per-cell-type h5ad (scVI/scANVI-integrated), builds a Seurat object, runs Harmony + hdWGCNA network construction. Shown for microglia; other cell types follow the same pattern. |
| `00_pseudobulk.r` | Aggregates raw counts into sample-level pseudobulk profiles per cell type. |
| `01_regress_covariates.r` | Bootstrap linear regression (R=1000, `thisexp ~ Study + Age + Sex`) to remove technical/demographic covariates; median coefficient across bootstraps used for adjustment. |
| `02_build_feature_vectors.py` | Assembles a 7-element `[AST, MIC, ENDO, ODC, OPC, INH, EXC]` covariate-adjusted expression vector per gene per sample. |
| `03_vae_ensemble_pipeline.py` | Core MONET pipeline: masked VAE (128-dim latent) trained across 20 random seeds -> GMM (k=4) clustering per seed -> pairwise co-assignment matrix -> agglomerative consensus clustering -> Sankey diagrams. Includes `load_data()` for parsing/imputing/padding the raw feature-vector CSV. |
| `04_qc_consensus.py` | Co-assignment heatmap QC, ordered by final consensus cluster. |
| `05_cmap_heatmap.py` | CLUE/CMap drug-repurposing connectivity heatmap per subtype. |
| `06_clinical_boxplots.py` | Braak stage / CERAD score / MMSE boxplots by subtype (ROSMAP + SEA-AD, CERAD harmonized across cohorts). |

## Method summary

For each donor, cell-type-specific pseudobulk profiles are covariate-adjusted
and concatenated into a single multi-cell-type feature vector, then encoded
by a masked VAE into a 128-dimensional latent space. This is repeated across
20 independently trained random seeds; a Gaussian mixture model (k=4) is
fitted to each seed's latent embedding, and a consensus clustering is derived
from the resulting pairwise sample co-assignment matrix (agglomerative
clustering, average linkage). See the manuscript Methods for the full
mathematical description (masked reconstruction + KL loss, co-assignment and
consensus-distance formulas).

Note: `k = 4` was fixed a priori to resolve finer-grained molecular
heterogeneity than the silhouette-optimal solution (which favored k=2/3
across seeds); silhouette analysis was computed as a diagnostic, not as the
basis for selecting k.

## Note on `00a_h5ad_to_seurat.r`

Update the `reticulate::use_condaenv(...)` call (or set the `PYR_CONDAENV`
environment variable) to point at your own conda environment before running.

## Earlier exploratory approach

An earlier deep-learning trial (1D convolutional autoencoder + Harmony +
UMAP + DBSCAN) is preserved for reference in
[`legacy/01_cnn_subtyping_early_trial`](../../legacy/01_cnn_subtyping_early_trial).
It was not used for any results in the manuscript.
