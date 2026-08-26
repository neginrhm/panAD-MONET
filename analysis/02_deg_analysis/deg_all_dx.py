"""
deg_all_dx.py
=============
Wilcoxon DEG analysis across all diagnosis groups (1-4, MCI, Control).

Runs pairwise comparisons (both directions) and each-vs-REST for the
selected condition groups using all cells (not per cell type).

Inputs:
  - H5AD: AnnData object with condition annotations in 'Cluster_or_dx'
  - META: CSV with per-sample metadata for merging subgroup assignments

Outputs:
  - OUTDIR: directory containing:
    - Pairwise DEG CSVs and PDFs for all group pairs
    - Each-vs-REST DEG CSVs and PDFs

Usage:
  python deg_all_dx.py
"""

import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# --- Configuration ---

BASE_DIR = "./"
H5AD = os.path.join(BASE_DIR, "data/h5ad/12AD_matched_label_final.h5ad")
META = "20_seeds_final_subgroups_with_meta.csv"

OUTDIR = "deg_results_1-4_MCI_Control"
KEY = "Cluster_or_dx"
CANDIDATES = {"1", "2", "3", "4", "MCI", "Control"}

# --- Helper functions ---

def normalize_label(x):
    """Normalize label strings (e.g., '1.0' -> '1')."""
    s = str(x).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


def safe_name(x):
    """Create filesystem-safe name from label."""
    return normalize_label(x).replace(os.sep, "_").replace(" ", "_")


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

# --- Select groups ---

label_map = {c: normalize_label(c) for c in adata.obs[KEY].cat.categories}
selected = [c for c, norm in label_map.items() if norm in CANDIDATES]

if len(selected) < 2:
    raise ValueError(f"Need at least two of {sorted(CANDIDATES)} in '{KEY}'. Found: {selected}")

# --- Analysis: pairwise comparisons (both directions) ---

for i in range(len(selected)):
    for j in range(i + 1, len(selected)):
        g1, g2 = selected[i], selected[j]
        n1, n2 = safe_name(g1), safe_name(g2)

        # g1 vs g2
        sc.tl.rank_genes_groups(
            adata, groupby=KEY, groups=[g1], reference=g2,
            method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
        )
        sc.get.rank_genes_groups_df(adata, group=g1).to_csv(
            f"{OUTDIR}/DEG_{KEY}_{n1}_vs_{n2}.csv", index=False
        )
        sc.pl.rank_genes_groups(adata, groups=[g1], n_genes=20, show=False)
        plt.savefig(f"{OUTDIR}/DEG_{KEY}_{n1}_vs_{n2}.pdf", bbox_inches="tight")
        plt.close("all")
        adata.uns.pop("rank_genes_groups", None)

        # g2 vs g1
        sc.tl.rank_genes_groups(
            adata, groupby=KEY, groups=[g2], reference=g1,
            method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
        )
        sc.get.rank_genes_groups_df(adata, group=g2).to_csv(
            f"{OUTDIR}/DEG_{KEY}_{n2}_vs_{n1}.csv", index=False
        )
        sc.pl.rank_genes_groups(adata, groups=[g2], n_genes=20, show=False)
        plt.savefig(f"{OUTDIR}/DEG_{KEY}_{n2}_vs_{n1}.pdf", bbox_inches="tight")
        plt.close("all")
        adata.uns.pop("rank_genes_groups", None)

# --- Analysis: each vs REST (within selected set) ---

tmp = "tmp_selected"
vals = adata.obs[KEY].values
selset = set(selected)
adata.obs[tmp] = pd.Categorical([v if v in selset else np.nan for v in vals],
                                categories=selected)

for g in selected:
    ng = safe_name(g)
    sc.tl.rank_genes_groups(
        adata, groupby=tmp, groups=[g], reference="rest",
        method="wilcoxon", use_raw=use_raw, n_genes=adata.shape[1], tie_correct=True
    )
    sc.get.rank_genes_groups_df(adata, group=g).to_csv(
        f"{OUTDIR}/DEG_selected_{ng}_vs_REST.csv", index=False
    )
    sc.pl.rank_genes_groups(adata, groups=[g], n_genes=20, show=False)
    plt.savefig(f"{OUTDIR}/DEG_selected_{ng}_vs_REST.pdf", bbox_inches="tight")
    plt.close("all")
    adata.uns.pop("rank_genes_groups", None)

del adata.obs[tmp]

# --- Done ---

print("Done. Results in", OUTDIR, "Groups:", [safe_name(g) for g in selected])
