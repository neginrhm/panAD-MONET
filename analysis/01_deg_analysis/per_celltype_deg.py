"""
per_celltype_deg.py
===================
Per-cell-type Wilcoxon DEG analysis across AD subgroups (1-4), MCI, and Control.

For each cell type, runs:
  - Pairwise Wilcoxon rank-sum tests (both directions) between all condition pairs
  - Each condition vs REST (within the selected set)

Inputs:
  - H5AD: AnnData object with cell-type and condition annotations
  - META: CSV with per-sample metadata (optional; set to None to skip)

Outputs:
  - Per-cell-type directories under OUTROOT, each containing:
    - DEG CSV files for all pairwise comparisons
    - DEG PDF plots (top 20 genes per comparison)

Usage:
  python per_celltype_deg.py
"""

import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# --- Configuration ---

BASE_DIR = "./"
H5AD = os.path.join(BASE_DIR, "data/h5ad/13AD_matched_label_final.h5ad")
META = "20_seeds_final_subgroups_with_meta.csv"   # optional; set to None to skip merge
OUTROOT = "deg_celltype_wilcoxon"

SAMPLE_KEY = "Sample"
COND_KEY = "Cluster_or_dx"
CELL_KEY = "Celltype_markers"

CANDIDATES = {"MCI", "Control", "1", "2", "3", "4"}

# --- Helper functions ---

def normalize_label(x):
    """Normalize label strings (e.g., '1.0' -> '1')."""
    s = str(x).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


def safe_name(x):
    """Create filesystem-safe name from label."""
    return normalize_label(x).replace("/", "_").replace(" ", "_")


def run_pair_and_rest(ad, groups, outdir, groupby):
    """Run pairwise both directions + each vs REST within 'groups'."""
    os.makedirs(outdir, exist_ok=True)
    use_raw = False

    # Pairwise both directions
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            n1, n2 = safe_name(g1), safe_name(g2)

            sc.tl.rank_genes_groups(
                ad, groupby=groupby, groups=[g1], reference=g2,
                method="wilcoxon", use_raw=use_raw, n_genes=ad.shape[1], tie_correct=True
            )
            sc.get.rank_genes_groups_df(ad, group=g1).to_csv(
                f"{outdir}/DEG_{groupby}_{n1}_vs_{n2}.csv", index=False
            )
            sc.pl.rank_genes_groups(ad, groups=[g1], n_genes=20, show=False)
            plt.savefig(f"{outdir}/DEG_{groupby}_{n1}_vs_{n2}.pdf", bbox_inches="tight")
            plt.close("all")
            ad.uns.pop("rank_genes_groups", None)

            sc.tl.rank_genes_groups(
                ad, groupby=groupby, groups=[g2], reference=g1,
                method="wilcoxon", use_raw=use_raw, n_genes=ad.shape[1], tie_correct=True
            )
            sc.get.rank_genes_groups_df(ad, group=g2).to_csv(
                f"{outdir}/DEG_{groupby}_{n2}_vs_{n1}.csv", index=False
            )
            sc.pl.rank_genes_groups(ad, groups=[g2], n_genes=20, show=False)
            plt.savefig(f"{outdir}/DEG_{groupby}_{n2}_vs_{n1}.pdf", bbox_inches="tight")
            plt.close("all")
            ad.uns.pop("rank_genes_groups", None)

    # Each vs REST (within these groups only)
    tmp = "_tmp_groups_"
    vals = ad.obs[groupby].values
    gset = set(groups)
    ad.obs[tmp] = pd.Categorical([v if v in gset else np.nan for v in vals], categories=groups)

    for g in groups:
        ng = safe_name(g)
        sc.tl.rank_genes_groups(
            ad, groupby=tmp, groups=[g], reference="rest",
            method="wilcoxon", use_raw=use_raw, n_genes=ad.shape[1], tie_correct=True
        )
        sc.get.rank_genes_groups_df(ad, group=g).to_csv(
            f"{outdir}/DEG_{safe_name('REST_WITHIN_' + groupby)}_{ng}_vs_REST.csv", index=False
        )
        sc.pl.rank_genes_groups(ad, groups=[g], n_genes=20, show=False)
        plt.savefig(f"{outdir}/DEG_{safe_name('REST_WITHIN_' + groupby)}_{ng}_vs_REST.pdf", bbox_inches="tight")
        plt.close("all")
        ad.uns.pop("rank_genes_groups", None)

    del ad.obs[tmp]


# --- Load data ---

adata = sc.read_h5ad(H5AD)

# (optional) merge per-sample metadata
if META:
    df = pd.read_csv(META)
    df_ = df.drop_duplicates("Samples").set_index("Samples")
    adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype(str)
    df_.index = df_.index.astype(str)
    adata.obs = adata.obs.join(df_, on=SAMPLE_KEY, how="left", rsuffix="_meta")

# --- Preprocessing ---

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

for k in [COND_KEY, CELL_KEY]:
    adata.obs[k] = adata.obs[k].astype("category")

os.makedirs(OUTROOT, exist_ok=True)

# --- Analysis: iterate cell types ---

for ct in adata.obs[CELL_KEY].cat.categories:
    ad = adata[adata.obs[CELL_KEY] == ct].copy()
    present = [c for c in ad.obs[COND_KEY].cat.categories if c in set(ad.obs[COND_KEY])]

    # Keep only target groups (with normalization)
    groups = []
    for c in present:
        if normalize_label(c) in CANDIDATES:
            groups.append(c)
    if len(groups) < 2:
        continue

    outdir = os.path.join(OUTROOT, safe_name(ct))
    print(f"[{ct}] groups: {list(map(safe_name, groups))} -> {outdir}")
    run_pair_and_rest(ad, groups, outdir, COND_KEY)

# --- Done ---

print("Done. Results in", OUTROOT)
