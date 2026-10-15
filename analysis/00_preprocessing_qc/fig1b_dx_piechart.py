"""
Fig. 1b: diagnosis distribution pie charts (nuclei and individuals).

Input (working directory):
  13ad_3m_5k_atlas_final_renamed_endo_meta_filled_no_bc.h5ad
Output:
  figures/dx_piechart.pdf (+ .png)
"""

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

H5AD = "13ad_3m_5k_atlas_final_renamed_endo_meta_filled_no_bc.h5ad"
OUT_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)


# ── helpers ────────────────────────────────────────────────────────────────────
def read_col(f, key):
    grp = f["obs"][key]
    if "categories" in grp:
        cats = [c.decode() if isinstance(c, bytes) else c for c in grp["categories"][:]]
        codes = grp["codes"][:]
        return pd.Categorical.from_codes(codes, cats)
    vals = grp[:]
    if vals.dtype.kind in ("S", "O"):
        vals = [v.decode() if isinstance(v, bytes) else v for v in vals]
    return vals


# ── load just the two columns we need ─────────────────────────────────────────
print("Reading h5ad …")
with h5py.File(H5AD, "r") as f:
    dx    = read_col(f, "dx")
    donor = read_col(f, "external_donor_name_label")

obs = pd.DataFrame({"dx": dx, "donor": donor})
obs = obs[obs["dx"] != "unknown"].copy()

# ── counts ────────────────────────────────────────────────────────────────────
cell_counts  = obs["dx"].value_counts()
indiv_counts = obs.drop_duplicates("donor")["dx"].value_counts()

# ── palette ───────────────────────────────────────────────────────────────────
PALETTE = {
    "Control":  "#4393c3",   # blue
    "MCI":      "#f4a582",   # salmon
    "Dementia": "#d6604d",   # brick red
}

ORDER = ["Control", "MCI", "Dementia"]

# ── layout ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
fig.patch.set_facecolor("white")

def draw_pie(ax, counts, title, unit):
    labels  = [dx for dx in ORDER if dx in counts]
    sizes   = [counts[dx] for dx in labels]
    colors  = [PALETTE[dx] for dx in labels]
    total   = sum(sizes)

    def autopct(pct):
        n = int(round(pct / 100 * total))
        return f"{pct:.1f}%\n({n:,})"

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        autopct=autopct,
        startangle=90,
        pctdistance=0.72,
        wedgeprops=dict(linewidth=1.5, edgecolor="white"),
        textprops=dict(fontsize=11),
    )

    for at in autotexts:
        at.set_fontsize(10)
        at.set_color("white")
        at.set_fontweight("bold")

    # centre label
    ax.text(0, 0, f"n = {total:,}\n{unit}", ha="center", va="center",
            fontsize=11, fontweight="bold", color="#333333")

    ax.set_title(title, fontsize=14, fontweight="bold", pad=14)

    # legend
    patches = [mpatches.Patch(color=PALETTE[dx], label=f"{dx}  ({counts[dx]:,})") for dx in labels]
    ax.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5, -0.12),
              frameon=False, fontsize=10, ncol=len(labels))


draw_pie(axes[0], cell_counts,  "Diagnosis — Cells",       "cells")
draw_pie(axes[1], indiv_counts, "Diagnosis — Individuals", "individuals")

fig.suptitle("AD Atlas — Diagnosis Distribution", fontsize=16, fontweight="bold", y=1.01)
fig.tight_layout()

for ext in ("pdf", "png"):
    path = os.path.join(OUT_DIR, f"dx_piechart.{ext}")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"Saved → {path}")

plt.close(fig)
print("Done.")
