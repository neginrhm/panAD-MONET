"""
CMap/CLUE Heatmap Generator (hierarchically clustered, all compounds)
======================================================================
Description:
    Generates a horizontal, hierarchically-clustered heatmap of CLUE/CMap
    connectivity scores (norm_cs) for the full top-N-per-signature union of
    compounds across the four AD subtypes (S1=AD1 ... S4=AD4), with vector
    (editable-text) PDF/SVG output.

    AD1 (S1) has no compounds meeting the significance/connectivity-score
    threshold applied in 04b_parse_clue_gct.R (see that script's
    norm_cs < -1.2 & fdr_q_nlog10 > 1.3 filter) -- its column is populated
    from the unfiltered top-ranked hits instead, per Methods ("Drug repurposing analysis"), and its
    input file is named accordingly (s1_nothreshold_...).

    NOTE ON FIG. 8C / SUPPL. TABLE 4: this script reproduces the full,
    unfiltered heatmap exactly as generated computationally. The published
    figure panel was then further curated by hand in Adobe Illustrator
    (removing a subset of compound labels for visual clarity/space) -- that
    curation step is editorial, not a computed filter, and is not scripted
    here. Do not treat any hardcoded "keep only these compounds" list as a
    computational result; none is applied by this script.

Inputs (relative to repo root):
    data/cmap/s1_nothreshold_cmap_top_compounds.csv
    data/cmap/s2_cmap_top_compounds.csv
    data/cmap/s3_cmap_top_compounds.csv
    data/cmap/s4_cmap_top_compounds.csv
    (produced by 04b_parse_clue_gct.R from the raw CLUE query .gct output --
    not included in this repo; see that script's docstring.)

Outputs:
    figures/cmap_heatmap_clean_colors.png / .pdf / .svg
    (open the .pdf/.svg in Illustrator/Inkscape to reproduce the manual
    compound-label curation applied for the published figure)

Usage:
    python 05_cmap_heatmap.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from matplotlib.colors import LinearSegmentedColormap

# Make exported text editable in Illustrator/Inkscape (TrueType, not outlines)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# -----------------------------
# CONFIG
# -----------------------------
ROOT     = Path(__file__).parent.parent
DATA_DIR = ROOT / "data" / "cmap"
OUTDIR   = ROOT / "figures"
OUTDIR.mkdir(exist_ok=True)

files = {
    "S1": str(DATA_DIR / "s1_nothreshold_cmap_top_compounds.csv"),
    "S2": str(DATA_DIR / "s2_cmap_top_compounds.csv"),
    "S3": str(DATA_DIR / "s3_cmap_top_compounds.csv"),
    "S4": str(DATA_DIR / "s4_cmap_top_compounds.csv"),
}

score_col = "norm_cs"
top_n_per_signature = 30
drug_col_override = None

# Interpretation: +1 = promoting disease-like state, -1 = treating/reversing
more_effective_direction = "more_negative"
ascending = True if more_effective_direction == "more_negative" else False

# Optional: require presence in >= K signatures
min_signatures = 1

output_png = str(OUTDIR / "cmap_heatmap_clean_colors.png")
output_pdf = str(OUTDIR / "cmap_heatmap_clean_colors.pdf")
output_svg = str(OUTDIR / "cmap_heatmap_clean_colors.svg")

# -----------------------------
# HELPERS
# -----------------------------
def guess_drug_col(df):
    candidates = [
        "pert_iname", "pert_name", "compound", "compound_name",
        "drug", "drug_name", "perturbagen", "cmap_name",
        "name", "pert_desc", "id"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError("Could not identify drug column.")

def load_one(signature, path):
    df = pd.read_csv(path)
    drug_col = drug_col_override or guess_drug_col(df)
    df = df[[drug_col, score_col]].copy()
    df.rename(columns={drug_col: "drug"}, inplace=True)
    df["signature"] = signature
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=["drug", score_col])
    return df

# -----------------------------
# LOAD DATA
# -----------------------------
dfs = []
for sig, fname in files.items():
    if not os.path.exists(fname):
        raise FileNotFoundError(f"Missing file: {fname}")
    dfs.append(load_one(sig, fname))

all_df = pd.concat(dfs, ignore_index=True)
grouped = all_df.groupby(["drug", "signature"], as_index=False)[score_col].mean()

# -----------------------------
# SELECT DRUGS (TOP-N UNION -- no further filtering)
# -----------------------------
top_sets = []
for sig in files.keys():
    top_sig = (
        grouped[grouped["signature"] == sig]
        .sort_values(score_col, ascending=ascending)
        .head(top_n_per_signature)
    )
    top_sets.append(set(top_sig["drug"]))

drug_union = sorted(set.union(*top_sets))
plot_df = grouped[grouped["drug"].isin(drug_union)].copy()
mat = plot_df.pivot_table(index="drug", columns="signature", values=score_col)

if min_signatures > 1:
    present = mat.notna().sum(axis=1)
    mat = mat.loc[present >= min_signatures]

# -----------------------------
# RESCALE TO [-1, 1]
# -----------------------------
max_abs = np.nanmax(np.abs(mat.values))
mat_scaled = (mat / max_abs).clip(-1, 1)

# -----------------------------
# HIERARCHICAL REORDERING (NO TREE)
# -----------------------------
mat_for_clust = mat_scaled.fillna(0.0)
dist = pdist(mat_for_clust.values, metric="euclidean")
Z = linkage(dist, method="average")
order = leaves_list(Z)
mat_scaled = mat_scaled.iloc[order]

# -----------------------------
# HORIZONTAL HEATMAP
# -----------------------------
mat_h = mat_scaled.T
heat_values = np.ma.masked_invalid(mat_h.values)

cmap = LinearSegmentedColormap.from_list(
    "blue_white_red",
    [
        (0.00, "#08306B"),
        (0.25, "#2171B5"),
        (0.50, "#FFFFFF"),
        (0.75, "#CB181D"),
        (1.00, "#67000D"),
    ],
    N=256
)
cmap.set_bad(color="#D1D5DB")

fig, ax = plt.subplots(figsize=(max(18, 0.25 * mat_h.shape[1]), 4))

# pcolormesh so every cell is a vector polygon (editable in Illustrator/Inkscape)
n_rows, n_cols = mat_h.shape
im = ax.pcolormesh(
    np.arange(n_cols + 1),
    np.arange(n_rows + 1),
    heat_values,
    cmap=cmap,
    vmin=-1,
    vmax=1,
    shading="flat",
    edgecolors="none",
    rasterized=False,
)
ax.invert_yaxis()
ax.set_aspect("auto")

ax.set_yticks(np.arange(n_rows) + 0.5)
ax.set_yticklabels(mat_h.index)
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels(mat_h.columns, rotation=90, fontsize=7)
ax.set_xlim(0, n_cols)
ax.set_ylim(n_rows, 0)

cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("Disease-state score")
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels(["Treating (-1)", "Neutral (0)", "Promoting (+1)"])
cbar.solids.set_rasterized(False)
cbar.solids.set_edgecolor("face")

ax.set_title(
    f"Clustered CLUE Heatmap\n"
    f"Top {top_n_per_signature} union/signature | Scale [-1,1]"
)

plt.tight_layout()
plt.savefig(output_png, dpi=300, bbox_inches="tight")
plt.savefig(output_pdf, bbox_inches="tight")
plt.savefig(output_svg, bbox_inches="tight")
plt.show()

print("Saved:", os.path.abspath(output_png))
print("Saved:", os.path.abspath(output_pdf))
print("Drugs plotted:", mat_h.shape[1])
