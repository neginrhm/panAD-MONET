#!/usr/bin/env python3
"""
07_clinical_stats.py
=====================
Description:
    Braak/CERAD/MMSE comparison across Control, MCI, and AD1-4, matching
    Fig. 7c / Suppl. Fig. 6b. Extracted and cleaned from the statistical cells
    of boxplots_pathology_traits_subtypes.ipynb (the rest of that notebook is
    exploratory data loading/merging and is not reproduced here).

    Note: 06_clinical_boxplots.py plots AD1-4 only, from a pre-merged CSV that
    does not include Control/MCI. This script instead operates on the final
    atlas h5ad's per-cell obs metadata (Control/MCI/AD1-4 all present) and
    produces both the boxplot input CSV and the omnibus/pairwise statistics
    tables behind Fig. 7c.

Method (Methods, "Clinical trait comparison across subtypes"):
    - CERAD harmonization: ROSMAP's raw scale (1=Definite AD ... 4=No AD) is
      inverted via CERAD_fixed = 4 - CERAD_raw to match SEA-AD's convention
      (0=Absent ... 3=Frequent). SEA-AD's CERAD is left unchanged.
    - Comparison performed at the donor level (one row per Sample/Study, not
      per cell) -- required for valid statistics.
    - Kruskal-Wallis omnibus test per trait (MMSE, CERAD_fixed, Braak) across
      all groups with >=2 donors, then pairwise two-sided Mann-Whitney U
      between every pair of groups, with Benjamini-Hochberg FDR correction
      applied within each trait.

Inputs:
    data/h5ad/13ad_3m_5k_atlas_final.h5ad -- final atlas (output of
    ../00_preprocessing_qc/), with obs columns: Sample, Study, CERAD, Braak,
    cts_mmse30_lv, Cluster_or_dx (values: "Control", "MCI", "1.0".."4.0")

Outputs:
    data/clinical_combined_ROSMAP_SEAAD_full.csv -- sample-level table (one row per Sample/Study) with
        CERAD_fixed, ready for boxplotting across Control/MCI/AD1-4
    figures/clinical_stats/stat_tests_omnibus_kruskalwallis.csv
    figures/clinical_stats/stat_tests_pairwise_mannwhitney.csv

Usage:
    python 07_clinical_stats.py
"""

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT    = Path(__file__).parent.parent.parent
H5AD    = ROOT / "data" / "h5ad" / "13ad_3m_5k_atlas_final.h5ad"
OUT_CSV = ROOT / "data" / "clinical_combined_ROSMAP_SEAAD_full.csv"
OUT_DIR = ROOT / "figures" / "clinical_stats"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GROUP_COL = "Cluster_or_dx"
GROUP_ORDER = ["Control", "MCI", "1.0", "2.0", "3.0", "4.0"]

# ---------------------------------------------------------------------------
# 1) Load obs metadata and standardize CERAD onto SEA-AD's 0=Absent...3=Frequent scale
# ---------------------------------------------------------------------------
print(f"Reading {H5AD}...")
adata = sc.read_h5ad(H5AD, backed="r")
obs = adata.obs.copy()

is_rosmap = obs["Study"] == "ROSMAP"
cerad_num = pd.to_numeric(obs["CERAD"], errors="coerce")

# ROSMAP raw CERAD (ceradsc): 1=Definite AD (most severe) ... 4=No AD (least severe)
# SEA-AD raw CERAD:            0=Absent (least severe) ... 3=Frequent (most severe)
cerad_fixed = cerad_num.copy()
cerad_fixed[is_rosmap] = 4 - cerad_num[is_rosmap]
obs["CERAD_fixed"] = cerad_fixed
obs["CERAD_label_fixed"] = obs["CERAD_fixed"].map(
    {0: "Absent", 1: "Sparse", 2: "Moderate", 3: "Frequent"}
)

# ---------------------------------------------------------------------------
# 2) Sample-level table (one row per Sample/Study) -- required for valid statistics
# ---------------------------------------------------------------------------
donor = obs.drop_duplicates(subset=["Sample", "Study"]).copy()
donor = donor[donor[GROUP_COL].notna()]
donor.to_csv(OUT_CSV, index=False)
print(f"Saved sample-level clinical table ({len(donor)} samples) -> {OUT_CSV}")

# ---------------------------------------------------------------------------
# 3) Kruskal-Wallis omnibus + pairwise Mann-Whitney U with BH-FDR
# ---------------------------------------------------------------------------
def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n)
    out[order] = fdr
    return out

test_cols = {
    "cts_mmse30_lv": "MMSE",
    "CERAD_fixed":   "CERAD (fixed scale, 0=Absent...3=Frequent)",
    "Braak":         "Braak",
}

omnibus_rows = []
pairwise_rows = []

for col, label in test_cols.items():
    sub = donor[[GROUP_COL, col]].copy()
    sub[col] = pd.to_numeric(sub[col], errors="coerce")
    sub = sub.dropna()

    groups_present = [g for g in GROUP_ORDER if g in sub[GROUP_COL].unique()]
    group_data = {g: sub.loc[sub[GROUP_COL] == g, col].values for g in groups_present}
    group_n = {g: len(v) for g, v in group_data.items()}
    valid_groups = [g for g in groups_present if group_n[g] >= 2]

    if len(valid_groups) >= 2:
        stat, p = stats.kruskal(*[group_data[g] for g in valid_groups])
    else:
        stat, p = np.nan, np.nan

    omnibus_rows.append({
        "trait": label, "column": col, "test": "Kruskal-Wallis",
        "groups": ",".join(valid_groups),
        "n_total": sum(group_n[g] for g in valid_groups),
        "statistic": stat, "p_value": p,
    })

    for g1, g2 in itertools.combinations(valid_groups, 2):
        u, p2 = stats.mannwhitneyu(group_data[g1], group_data[g2], alternative="two-sided")
        pairwise_rows.append({
            "trait": label, "column": col, "test": "Mann-Whitney U",
            "group1": g1, "n1": group_n[g1], "median1": np.median(group_data[g1]),
            "group2": g2, "n2": group_n[g2], "median2": np.median(group_data[g2]),
            "statistic": u, "p_value": p2,
        })

omnibus_df = pd.DataFrame(omnibus_rows)
pairwise_df = pd.DataFrame(pairwise_rows)

pairwise_df["p_fdr_bh"] = np.nan
for trait in pairwise_df["trait"].unique():
    mask = pairwise_df["trait"] == trait
    pairwise_df.loc[mask, "p_fdr_bh"] = benjamini_hochberg(pairwise_df.loc[mask, "p_value"].values)
pairwise_df["significant_fdr05"] = pairwise_df["p_fdr_bh"] < 0.05

omnibus_df.to_csv(OUT_DIR / "stat_tests_omnibus_kruskalwallis.csv", index=False)
pairwise_df.to_csv(OUT_DIR / "stat_tests_pairwise_mannwhitney.csv", index=False)

print(omnibus_df)
print()
print(pairwise_df.sort_values("p_fdr_bh").to_string())
print(f"\nSaved statistics tables -> {OUT_DIR}")
