#!/usr/bin/env python3
"""Combined supplementary figure: all 4 QC violin+box panels (total_counts,
n_genes_by_counts, pct_counts_MT, doublet_score) by Study, ALL 13 studies,
stacked into ONE figure sized exactly for Nature Neuroscience (7.2 x 8 in).
Same as make_qc_violin_supp_figure.py but with inner="box" (quartile box +
whiskers + median dot drawn inside each violin).

Font: Arial if installed, otherwise Liberation Sans (the metric-compatible
open-source substitute - identical glyph widths/kerning),
so swapping the font family to Arial in Illustrator afterwards causes no
reflow. Text is saved as true editable vector type (pdf.fonttype=42), not
outlines, so it is selectable/editable in Adobe Illustrator. All font sizes
are >=7pt per journal spec.
"""

import os
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

H5AD = "13ad_3m_5k_atlas_final_renamed_endo_meta_filled_no_bc_fig1_colors.h5ad"
OUTDIR = "qc_supp1_v5_scanpy_13studies"
os.makedirs(OUTDIR, exist_ok=True)

FONT_AVAILABLE = {f.name for f in fm.fontManager.ttflist}
FONT = "Arial" if "Arial" in FONT_AVAILABLE else "Liberation Sans" if "Liberation Sans" in FONT_AVAILABLE else "DejaVu Sans"
print(f"Using font: {FONT} (Arial installed: {'Arial' in FONT_AVAILABLE})")

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,   # TrueType-embedded -> real editable text in Illustrator, not outlines
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.6,
})
sc.set_figure_params(facecolor="white", frameon=True)
plt.rcParams["grid.linewidth"] = 0.35

STUDY_ORDER = [
    "AD_Grubman_2019", "AD_Lau_2020", "AD_Leng_2021", "AD_Mathys_2019",
    "AD_Otero-Garcia_2020", "AD_Swarup", "AD_Zhou_2020", "ADvas_Yang_2021",
    "ADDS", "Morabito_2020", "Olah_microglia_AD", "ROSMAP", "SEA-AD",
]
DISPLAY = {
    "AD_Grubman_2019": "Grubman et al., 2019",
    "AD_Lau_2020": "Lau et al., 2020",
    "AD_Leng_2021": "Leng et al., 2021",
    "AD_Mathys_2019": "Mathys et al., 2019",
    "AD_Otero-Garcia_2020": "Otero-Garcia et al., 2020",
    "AD_Swarup": "Swarup et al., 2021",
    "AD_Zhou_2020": "Zhou et al., 2020",
    "ADvas_Yang_2021": "Yang et al., 2021",
    "ADDS": "DSAD, 2024",
    "Morabito_2020": "Morabito et al., 2020",
    "Olah_microglia_AD": "Olah et al., 2020",
    "ROSMAP": "ROSMAP",
    "SEA-AD": "SEA-AD",
}
COLORS = {
    "ROSMAP": "#1F78B4", "SEA-AD": "#33A02C", "ADDS": "#E31A1C",
    "AD_Leng_2021": "#FF7F00", "AD_Zhou_2020": "#6A3D9A", "ADvas_Yang_2021": "#B15928",
    "AD_Mathys_2019": "#A6CEE3", "AD_Swarup": "#B2DF8A", "Morabito_2020": "#FB9A99",
    "AD_Otero-Garcia_2020": "#FDBF6F", "AD_Grubman_2019": "#CAB2D6",
    "Olah_microglia_AD": "#FFFF99", "AD_Lau_2020": "#B2E2E2",
}
METRICS = [
    ("A", "total_counts", "Total counts", (0, 80000)),
    ("B", "n_genes_by_counts", "Genes detected", (0, 10000)),
    ("C", "pct_counts_MT", "% mitochondrial counts", (0, 10)),
    ("D", "doublet_score", "Doublet score", (0, 0.10)),
]

print("Loading obs columns via h5py...")
f = h5py.File(H5AD, "r")
study_g = f["obs"]["Study"]
study_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in study_g["categories"][:]])
study_codes = study_g["codes"][:]
study = study_cats[study_codes]

qc_data = {"Study": study}
for _, col, _, _ in METRICS:
    qc_data[col] = f["obs"][col][:]
qc = pd.DataFrame(qc_data)
qc = qc[qc["Study"].isin(STUDY_ORDER)].copy()

label_order = [DISPLAY[c] for c in STUDY_ORDER]
palette = {DISPLAY[c]: COLORS[c] for c in STUDY_ORDER}
qc["Study"] = pd.Categorical(qc["Study"].map(DISPLAY), categories=label_order)
qc = qc.reset_index(drop=True)
print("n cells used:", len(qc))

X = sp.csr_matrix((len(qc), 1), dtype="float32")
qc_adata = ad.AnnData(X=X, obs=qc)

# 4 stacked panels; only the bottom panel carries the (rotated) study x-tick
# labels, so it gets a taller height ratio to leave room for them.
fig, axes = plt.subplots(
    4, 1, figsize=(7.2, 8.0), constrained_layout=True,
    gridspec_kw={"height_ratios": [1.0, 1.0, 1.0, 1.35]},
)

for ax, (letter, col, ylabel, ylim) in zip(axes, METRICS):
    # inner="box" draws a mini box-and-whisker (thin IQR bar, whisker line,
    # white median dot) inside each violin
    ax = sc.pl.violin(
        qc_adata, keys=col, groupby="Study", rotation=45, stripplot=False,
        cut=0, palette=palette, show=False, ax=ax, inner="box",
    )
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_xlabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    ax.tick_params(axis="both", width=0.6, length=2.5, labelsize=7)
    ax.text(-0.06, 1.05, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="right")

# only bottom panel keeps study-name tick labels; upper 3 panels get blank
# labels (categories still line up 1:1 across panels since Study order/count
# is identical) to cut redundant clutter in a stacked multi-panel figure
for ax in axes[:-1]:
    ax.set_xticklabels([])
axes[-1].set_xticklabels(label_order, rotation=45, ha="right", fontsize=7)

out_pdf = f"{OUTDIR}/qc_violin_box_supp_figure_13studies.pdf"
out_png = f"{OUTDIR}/qc_violin_box_supp_figure_13studies.png"
# no bbox_inches="tight" -- constrained_layout keeps the canvas at exactly
# the requested 7.2x8.0in; "tight" would crop to content and shrink the page
# below the journal-specified size (see project_qc_supp1_figure_regen memory)
fig.savefig(out_pdf)
fig.savefig(out_png, dpi=300)
plt.close(fig)
print("Saved:", out_pdf)
print("Saved:", out_png)
print("Done.")
