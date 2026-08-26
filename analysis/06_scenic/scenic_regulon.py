"""
scenic_regulon.py
=================
Compute regulon Z-scores and UMAP visualization for microglia from pySCENIC
AUCell output.

Inputs:
    - Microglia motif/regulon CSV: <SCENIC_DIR>/results_raw_count/microglia_dem_reg.csv
    - Microglia AUCell matrix:     <SCENIC_DIR>/results_raw_count/microglia_auc_matrix.csv
    - Pre-loaded AnnData object `all` with obs column 'Celltype_markers'

Outputs:
    - microglia_regulons.dat              (pickled regulon objects)
    - microglia_scenic_zscore_df_heatmap.csv (regulon Z-score table)
    - UMAP PDF colored by Sex, Diagnosis, Tissue

Usage:
    python scenic_regulon.py
"""

import pickle

import pandas as pd
import scanpy as sc

from pyscenic.utils import load_motifs

# --- Configuration -----------------------------------------------------------

BASE_DIR = "."
SCENIC_DIR = BASE_DIR

MOTIF_CSV = f"{SCENIC_DIR}/results_raw_count/microglia_dem_reg.csv"
AUC_CSV = f"{SCENIC_DIR}/results_raw_count/microglia_auc_matrix.csv"
REGULONS_DAT_FNAME = "microglia_regulons.dat"
ZSCORE_OUT_CSV = "microglia_scenic_zscore_df_heatmap.csv"

# --- Load data ---------------------------------------------------------------

df_motifs_microglia = load_motifs(MOTIF_CSV)
regulons_microglia = derive_regulons(df_motifs_microglia)
print(f"Number of regulons: {len(regulons_microglia)}")

auc_mtx_microglia = pd.read_csv(AUC_CSV, index_col=0)

# --- Save pickled regulons ---------------------------------------------------

with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons_microglia, f)

# --- Subset to microglia -----------------------------------------------------

microglia = all[all.obs.Celltype_markers == "microglia"]
add_scenic_metadata(microglia, auc_mtx_microglia, regulons_microglia)

# --- UMAP from AUCell representation ----------------------------------------

sc.pp.neighbors(microglia, use_rep='X_aucell')
sc.tl.umap(microglia)

sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
sc.pl.umap(
    microglia,
    color=['Sex', 'Diagnosis2', 'Tissue'],
    title=['sex', 'Diagnosis', 'Tissue'],
    ncols=3,
    save='-dementia_scenic_microglia.pdf',
)

# --- Sort AUC matrix by mean activity ---------------------------------------

auc_mtx_microglia.loc["MEAN_ROW"] = auc_mtx_microglia.mean(axis=0)
auc_mtx_microglia = auc_mtx_microglia[
    auc_mtx_microglia.loc["MEAN_ROW"].sort_values(ascending=False).index
]

# --- Compute regulon Z-scores -----------------------------------------------

df_obs = microglia.obs

# Select only regulon-related numeric columns
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(
    filter(lambda s: s.startswith('Regulon('), signature_column_names)
)

# Sum regulon scores across all cells
df_regulon_sums = df_obs[signature_column_names].sum(axis=0).to_frame(name="sum_score")

# Compute Z-scores for regulons
df_regulon_sums["Z"] = (
    (df_regulon_sums["sum_score"] - df_regulon_sums["sum_score"].mean())
    / df_regulon_sums["sum_score"].std()
)

# Extract regulon names (removing "Regulon()" wrapper)
df_regulon_sums.index = df_regulon_sums.index.str[8:-1]
df_regulon_sums = df_regulon_sums.reset_index().rename(columns={"index": "regulon"})

df_results_microglia = df_regulon_sums

# --- Save results ------------------------------------------------------------

df_results_microglia.to_csv(ZSCORE_OUT_CSV)
print(f"Z-score table saved to {ZSCORE_OUT_CSV}")
