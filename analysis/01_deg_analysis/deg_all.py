"""
deg_all.py
==========
Wilcoxon DEG analysis: all condition groups vs REST.

Runs a single Wilcoxon rank-sum test across all groups in 'Cluster_or_dx',
saving combined and per-group CSVs and PDFs.

Inputs:
  - H5AD: AnnData object with condition annotations in 'Cluster_or_dx'
  - META: CSV with per-sample metadata for merging subgroup assignments

Outputs:
  - OUTDIR: directory containing:
    - Combined all-vs-REST DEG CSV and PDF
    - Per-group DEG CSVs and top-20 PDF plots

Usage:
  python deg_all.py
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# --- Configuration ---

BASE_DIR = "./"
H5AD = os.path.join(BASE_DIR, "data/h5ad/13AD_matched_label_final.h5ad")
META = "20_seeds_final_subgroups_with_meta.csv"

OUTDIR = "deg_results2"
KEY = "Cluster_or_dx"

# --- Load data ---

adata = sc.read_h5ad(H5AD)
df4 = pd.read_csv(META)

# Merge metadata onto adata.obs by Sample/Samples
obs_key = "Sample"
df_key = "Samples"

df4_ = (df4
        .drop_duplicates(subset=[df_key])
        .set_index(df_key))
adata.obs[obs_key] = adata.obs[obs_key].astype(str)
df4_.index = df4_.index.astype(str)

adata.obs = adata.obs.join(df4_, on=obs_key, how="left", rsuffix="_df4")

# --- Preprocessing ---

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

os.makedirs(OUTDIR, exist_ok=True)
adata.obs[KEY] = adata.obs[KEY].astype("category")
cats = list(adata.obs[KEY].cat.categories)

# --- Analysis: all groups vs REST ---

sc.tl.rank_genes_groups(
    adata,
    groupby=KEY,
    reference="rest",
    method="wilcoxon",
    use_raw=False,
    tie_correct=True,
)

# --- Save results ---

# Combined CSV
sc.get.rank_genes_groups_df(adata).to_csv(
    f"{OUTDIR}/DEG_{KEY}_ALL_vs_REST.csv", index=False
)

# Per-group CSVs
for g in cats:
    sc.get.rank_genes_groups_df(adata, group=g).to_csv(
        f"{OUTDIR}/DEG_{KEY}_{g}_vs_REST.csv", index=False
    )

# Combined plot (top-20 per group)
sc.pl.rank_genes_groups(adata, n_genes=20, show=False)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/DEG_{KEY}_ALL_vs_REST.pdf", bbox_inches="tight")
plt.close()

# Per-group plots (top-20)
for g in cats:
    sc.pl.rank_genes_groups(adata, groups=[g], n_genes=20, show=False)
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/DEG_{KEY}_{g}_top20.pdf", bbox_inches="tight")
    plt.close()

# --- Done ---

print("Saved CSVs (combined + per-group) and PDFs in:", OUTDIR)
