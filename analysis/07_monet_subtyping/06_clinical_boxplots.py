#!/usr/bin/env python3
"""
Generate CERAD, BRAAK, and MMSE30 IV boxplots by AD subtype (AD1–AD4).

Inputs (relative to repo root):
    data/clinical_combined_ROSMAP_SEAAD.csv

Outputs:
    figures/clinical_boxplots.png / .pdf
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT   = Path(__file__).parent.parent
DATA   = ROOT / "data" / "clinical_combined_ROSMAP_SEAAD.csv"
OUTDIR = ROOT / "figures"
OUTDIR.mkdir(exist_ok=True)

df     = pd.read_csv(DATA)
COLORS = {"AD1": "#7fbfff", "AD2": "#90ee90", "AD3": "#ffb6b6", "AD4": "#cc3333"}
ORDER  = ["AD1", "AD2", "AD3", "AD4"]

panels = [
    ("CERAD", "CERAD",     [0,1,2,3,4],       ["Absent","Sparse","Moderate","Frequent",""]),
    ("Braak", "BRAAK",     [1,2,3,4,5,6],      ["1","2","3","4","5","6"]),
    ("MMSE",  "MMSE30 IV", [0,5,10,15,20,25,30],["0","5","10","15","20","25","30"]),
]

fig, axes = plt.subplots(1, 3, figsize=(13, 5))

for ax, (col, ylabel, yticks, ylabels) in zip(axes, panels):
    groups = [df.loc[df["subtype"] == s, col].dropna().values for s in ORDER]
    bp = ax.boxplot(groups, patch_artist=True, widths=0.5,
                    medianprops=dict(color="black", linewidth=1.5),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2),
                    flierprops=dict(marker="o", markersize=3, alpha=0.4))
    for patch, s in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[s])
        patch.set_alpha(0.85)

    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(ORDER, fontsize=11)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
    for i, (s, g) in enumerate(zip(ORDER, groups), 1):
        ax.text(i, ax.get_ylim()[0] - yrange * 0.09,
                f"n={len(g)}", ha="center", va="top", fontsize=8, color="#444")

fig.suptitle("Clinical Traits by AD Subtype", fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(OUTDIR / "clinical_boxplots.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTDIR / "clinical_boxplots.pdf", bbox_inches="tight")
print("Saved figures/clinical_boxplots.png/.pdf")
