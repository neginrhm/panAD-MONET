"""
deg_ad1-4.py
============
Wilcoxon DEG analysis restricted to AD subgroups 1-4.

Runs pairwise comparisons (both directions) and each-vs-REST among
subgroups {1, 2, 3, 4} using all cells (not per cell type).

Inputs:
  - H5AD: AnnData object with condition annotations in 'Cluster_or_dx'
  - META: CSV with per-sample metadata for merging subgroup assignments

Outputs:
  - OUTDIR: directory containing:
    - Pairwise DEG CSVs and PDFs for subgroup pairs
    - Each-vs-REST DEG CSVs and PDFs

Usage:
  python deg_ad1-4.py
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# --- Configuration ---

BASE_DIR = "./"
H5AD = os.path.join(BASE_DIR, "data/h5ad/12AD_matched_label_final.h5ad")
META = "20_seeds_final_subgroups_with_meta.csv"

OUTDIR = "deg_results_1-4"
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
assert KEY in adata.obs.columns, f"'{KEY}' not found in adata.obs.columns"
adata.obs[KEY] = adata.obs[KEY].astype("category")
sc.settings.autoshow = False
use_raw = False

# --- Select subgroups 1-4 ---

want = {"1", "2", "3", "4"}
four = []
for c in adata.obs[KEY].cat.categories:
    s = str(c)
    if s.endswith(".0"):
        s = s[:-2]
    if s in want:
        four.append(c)

if len(four) < 2:
    raise ValueError(f"Need at least two of {{1,2,3,4}} in '{KEY}'. Found: {four}")

# --- Analysis: pairwise comparisons (both directions) ---

for i in range(len(four)):
    for j in range(i + 1, len(four)):
        g1, g2 = four[i], four[j]

        # g1 vs g2
        sc.tl.rank_genes_groups(
            adata, groupby=KEY, groups=[g1], reference=g2,
            method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
        )
        sc.get.rank_genes_groups_df(adata, group=g1).to_csv(
            f"{OUTDIR}/DEG_{KEY}_{g1}_vs_{g2}.csv", index=False
        )
        sc.pl.rank_genes_groups(adata, groups=[g1], n_genes=20, show=False)
        plt.savefig(f"{OUTDIR}/DEG_{KEY}_{g1}_vs_{g2}.pdf", bbox_inches="tight")
        plt.close("all")
        if "rank_genes_groups" in adata.uns:
            del adata.uns["rank_genes_groups"]

        # g2 vs g1
        sc.tl.rank_genes_groups(
            adata, groupby=KEY, groups=[g2], reference=g1,
            method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
        )
        sc.get.rank_genes_groups_df(adata, group=g2).to_csv(
            f"{OUTDIR}/DEG_{KEY}_{g2}_vs_{g1}.csv", index=False
        )
        sc.pl.rank_genes_groups(adata, groups=[g2], n_genes=20, show=False)
        plt.savefig(f"{OUTDIR}/DEG_{KEY}_{g2}_vs_{g1}.pdf", bbox_inches="tight")
        plt.close("all")
        if "rank_genes_groups" in adata.uns:
            del adata.uns["rank_genes_groups"]

# --- Analysis: each vs REST (within subgroups 1-4 only) ---

tmp = "tmp_1234"
adat_vals = adata.obs[KEY].values
adata.obs[tmp] = pd.Categorical(
    [v if v in set(four) else np.nan for v in adat_vals],
    categories=four
)

for g in four:
    sc.tl.rank_genes_groups(
        adata, groupby=tmp, groups=[g], reference="rest",
        method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
    )
    sc.get.rank_genes_groups_df(adata, group=g).to_csv(
        f"{OUTDIR}/DEG_1234_{g}_vs_REST.csv", index=False
    )
    sc.pl.rank_genes_groups(adata, groups=[g], n_genes=20, show=False)
    plt.savefig(f"{OUTDIR}/DEG_1234_{g}_vs_REST.pdf", bbox_inches="tight")
    plt.close("all")
    if "rank_genes_groups" in adata.uns:
        del adata.uns["rank_genes_groups"]

del adata.obs[tmp]

# --- Done ---

print("Done. Results in", OUTDIR)
