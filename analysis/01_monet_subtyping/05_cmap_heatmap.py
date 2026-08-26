"""
CMAP Heatmap Generator
======================
Description:
    Generates a horizontal heatmap of CMap connectivity scores (norm_cs)
    for the union of top-N compounds across multiple signatures (S1-S4).

Inputs:
    - CSV files per signature (s1/s2/s3/s4_cmap_top_compounds.csv) with
      columns: drug identifier, norm_cs, moa

Outputs:
    - cmap_heatmap_horizontal.png  (300 dpi)
    - cmap_heatmap_horizontal.pdf

Usage:
    python cmap_heatmap.py
    (Run from the directory containing the signature CSV files.)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---

from pathlib import Path
DATA_DIR = Path(__file__).parent.parent / "data" / "cmap"

files = {
    "S1": str(DATA_DIR / "s1_nothrehsold_cmap_top_compounds.csv"),
    "S2": str(DATA_DIR / "s2_cmap_top_compounds.csv"),
    "S3": str(DATA_DIR / "s3_cmap_top_compounds.csv"),
    "S4": str(DATA_DIR / "s4_cmap_top_compounds.csv"),
}

score_col = "norm_cs"
kind_col = "moa"
top_n_per_signature = 30
drug_col_override = None

# Reversal case: more negative = better
more_effective_direction = "more_negative"
ascending = True if more_effective_direction == "more_negative" else False

output_png = "cmap_heatmap_horizontal.png"
output_pdf = "cmap_heatmap_horizontal.pdf"

# --- Helpers ---

def guess_drug_col(df: pd.DataFrame) -> str:
    """Auto-detect the column containing drug/compound names."""
    candidates = [
        "pert_iname", "pert_name", "compound", "compound_name",
        "drug", "drug_name", "perturbagen", "cmap_name", "name", "pert_desc"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    if "id" in df.columns:
        return "id"
    raise ValueError(f"Could not identify drug column. Columns: {list(df.columns)}")


def load_one(signature: str, path: str) -> pd.DataFrame:
    """Load a single signature CSV file, returning drug/score/signature columns."""
    df = pd.read_csv(path)
    drug_col = drug_col_override or guess_drug_col(df)

    cols = [drug_col, score_col]
    if kind_col in df.columns:
        cols.append(kind_col)

    df = df[cols].copy()
    df.rename(columns={drug_col: "drug"}, inplace=True)
    df["signature"] = signature

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=["drug", score_col])

    return df

# --- Load data ---

dfs = []
for sig, fname in files.items():
    if not os.path.exists(fname):
        raise FileNotFoundError(f"Missing file: {fname}")
    dfs.append(load_one(sig, fname))

all_df = pd.concat(dfs, ignore_index=True)

# --- Aggregate duplicates (mean norm_cs per drug per signature) ---

group_cols = ["drug", "signature"]
if kind_col in all_df.columns:
    grouped = all_df.groupby(group_cols, as_index=False).agg(
        **{
            score_col: (score_col, "mean"),
            kind_col: (kind_col, "first"),
        }
    )
else:
    grouped = all_df.groupby(group_cols, as_index=False)[score_col].mean()

# --- Select drugs (top-N union) ---

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

# Pivot to matrix: rows=drugs, cols=signatures
mat = plot_df.pivot_table(
    index="drug",
    columns="signature",
    values=score_col,
    aggfunc="mean"
)

# Order drugs by average effect
mat["avg"] = mat.mean(axis=1, skipna=True)
mat = mat.sort_values("avg", ascending=ascending).drop(columns="avg")

# --- Horizontal heatmap ---

mat_h = mat.T  # signatures on Y, drugs on X
masked = np.ma.masked_invalid(mat_h.values)

# Symmetric color scale around zero
vmax = np.nanmax(np.abs(mat_h.values))
vmin = -vmax

plt.figure(figsize=(max(18, 0.22 * mat_h.shape[1]), 4))

im = plt.imshow(
    masked,
    aspect="auto",
    cmap="RdBu_r",
    vmin=vmin,
    vmax=vmax
)

plt.yticks(range(mat_h.shape[0]), mat_h.index)
plt.xticks(
    range(mat_h.shape[1]),
    mat_h.columns,
    rotation=90,
    fontsize=7
)

plt.colorbar(im, label=score_col)
plt.title(
    f"CLUE heatmap ({score_col})\n"
    f"Union of top {top_n_per_signature} drugs per signature"
)

plt.tight_layout()
plt.savefig(output_png, dpi=300, bbox_inches="tight")
plt.savefig(output_pdf, bbox_inches="tight")
plt.show()

print("Saved:", os.path.abspath(output_png))
print("Saved:", os.path.abspath(output_pdf))
print("Drugs plotted:", mat_h.shape[1])
