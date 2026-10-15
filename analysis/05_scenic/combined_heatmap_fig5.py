#!/usr/bin/env python3
"""
combined_heatmap_v1.py
======================
Combined TF x TG heatmap across all 5 cell types for Nature journal.

Rows   : top-6 TFs per cell type (30 rows total), grouped by cell type
         with a left-side color strip and cell-type label.
Columns: up to ~25 TGs: all AD-GWAS TGs found in any cell type's top-6-TF
         neighbourhood, then highest-scoring shared non-GWAS TGs.
         Ordered by waterfall (dominant cell type left → right).
Color  : NetRem weight_combined (plasma + PowerNorm); 0 where pair absent.

Output: scenic_v11_nature/celltype_outputs/combined_heatmap_all_celltypes.pdf/.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"]  = 42
matplotlib.rcParams["ps.fonttype"]   = 42
matplotlib.rcParams["font.family"]   = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Helvetica"]
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT     = Path(__file__).parent.parent.parent  # repo root
DATA_DIR = ROOT / "figures" / "scenic_v11_nature" / "celltype_outputs"  # written by regenerate_network_figures_fig5.py
OUT_DIR  = DATA_DIR          # same folder as other outputs

# ── Constants ─────────────────────────────────────────────────────────────────
AD_GWAS = {
    "APOE","SORL1","PICALM","CLU","BIN1","CR1","ABCA7","MEF2C","PTK2B",
    "NR4A1","NR4A2","RORA","RORB","BMP2K","RIN3","ATP8B4","ADAM10",
    "PLCG2","TREM2","GRN","MAPT","RUNX1","MAFB","OLIG1","OLIG2",
}

CELL_TYPES = ["AST","EXC","INH","MIC","ODC"]
CT_LABELS  = {
    "AST": "Astrocytes",
    "EXC": "Excitatory\nNeurons",
    "INH": "Inhibitory\nNeurons",
    "MIC": "Microglia",
    "ODC": "Oligodendrocytes",
}
CT_COLORS  = {
    "AST": "#3d9970",
    "EXC": "#e07b39",
    "INH": "#7b4fa6",
    "MIC": "#c0392b",
    "ODC": "#2980b9",
}

N_TF_PER_CT            = 6    # TFs shown per cell type (DTL + top-5 specific)
N_TG_COLS              = 28   # total TG columns
FORCE_TFS              = ["DTL"]   # always included if in that CT's top-12
HEATMAP_CMAP           = "plasma"
HEATMAP_GAMMA          = 0.65
GWAS_COLOR             = "#00e5ff"

# Match the exact v11 heatmap pipeline parameters
TOP_TF_TG_INTERSECTION = 300   # take top-300 pairs (CSV sorted by rank_score_netrem)
TOP_TF_PER_TG          = 3     # max 3 TFs per TG
TOP_TFS_HEATMAP        = 12    # top-12 TFs per CT (same as individual heatmaps)

# ── Load data — exact same pipeline as scenic_v11_nature ──────────────────────
def top_k_per_group(df, by_col, group_col, k):
    return (df.sort_values(by_col, ascending=False)
              .groupby(group_col, as_index=False, group_keys=False)
              .head(k))

def load_ct_heatmap_pool(ct):
    """Replicate the pair pool used by the individual v11 heatmaps exactly."""
    df = pd.read_csv(DATA_DIR / f"Shared_TF_TG_pairs_NetRem_x_pySCENIC_{ct}.csv")
    df["weight_combined"] = pd.to_numeric(df["weight_combined"], errors="coerce")
    df = df.dropna(subset=["weight_combined"])
    df = df.head(TOP_TF_TG_INTERSECTION)                      # top-300 by rank_score_netrem
    df = top_k_per_group(df, "weight_combined", "TG", TOP_TF_PER_TG)   # max 3 TFs/TG
    top12 = (df.groupby("TF")["weight_combined"].sum()
               .sort_values(ascending=False).head(TOP_TFS_HEATMAP).index.tolist())
    df = df[df["TF"].isin(top12)].copy()
    return df, top12

ct_data    = {}
ct_top12   = {}   # full top-12 per CT (used for TF specificity scoring)
for ct in CELL_TYPES:
    pool, top12 = load_ct_heatmap_pool(ct)
    ct_data[ct]  = pool
    ct_top12[ct] = top12

EPS = 1e-6

# ── Step 1: select TG columns from the FULL top-12 pool ──────────────────────
# Use ALL top-12 TF pairs to score TGs so no activity is lost by TF pruning.
# Score = gwas_bonus
#         + 4.0 × max_ct_weight      (absolute strength = interesting)
#         + 3.0 × max_ct_weight × specificity  (cell-type concentrated)
tg_ct_w_full = {}
for ct in CELL_TYPES:
    for _, r in ct_data[ct].iterrows():   # ct_data still has full top-12 here
        tg = r["TG"]
        if tg not in tg_ct_w_full:
            tg_ct_w_full[tg] = {c: 0.0 for c in CELL_TYPES}
        tg_ct_w_full[tg][ct] += r["weight_combined"]

tg_score = {}
for tg, cw in tg_ct_w_full.items():
    weights     = list(cw.values())
    total_w     = sum(weights) + EPS
    max_w       = max(weights)
    specificity = max_w / total_w
    gwas_bonus  = 1000.0 if tg.upper() in AD_GWAS else 0.0
    tg_score[tg] = gwas_bonus + 4.0 * max_w + 3.0 * max_w * specificity

# ── Step 2: select TFs per CT that best cover the selected TGs ───────────────
# For each CT, rank its top-12 TFs by total weight TO the selected TGs.
# This ensures rows are chosen to actually explain the visible columns.
# DTL is always forced in first.
ranked_tgs = sorted(tg_score, key=lambda t: -tg_score[t])
gwas_tgs   = [t for t in ranked_tgs if t.upper() in AD_GWAS]
other_tgs  = [t for t in ranked_tgs if t.upper() not in AD_GWAS]
selected_tgs_set = set((gwas_tgs + other_tgs)[:N_TG_COLS])

# Guarantee top-2 targets of each FORCE_TF per CT are included.
# DTL's targets are uniformly shared → low specificity score → otherwise excluded.
N_FORCE_TG_PER_CT = 2
for tf in FORCE_TFS:
    for ct in CELL_TYPES:
        sub = ct_data[ct][ct_data[ct]["TF"] == tf].sort_values(
            "weight_combined", ascending=False).head(N_FORCE_TG_PER_CT)
        for tg in sub["TG"]:
            selected_tgs_set.add(tg)

ct_top_tfs = {}
for ct in CELL_TYPES:
    # Weight each top-12 TF by its total connection to the selected TGs
    tf_tg_weight = {}
    for _, r in ct_data[ct].iterrows():
        if r["TG"] in selected_tgs_set:
            tf_tg_weight[r["TF"]] = tf_tg_weight.get(r["TF"], 0.0) + r["weight_combined"]

    forced    = [tf for tf in FORCE_TFS if tf in ct_top12[ct]]
    remaining = [tf for tf in ct_top12[ct] if tf not in forced]
    ranked_by_coverage = sorted(remaining, key=lambda t: -tf_tg_weight.get(t, 0.0))
    n_fill    = N_TF_PER_CT - len(forced)
    ct_top_tfs[ct] = forced + ranked_by_coverage[:n_fill]

# Filter ct_data to selected TFs now (after TG selection)
for ct in CELL_TYPES:
    ct_data[ct] = ct_data[ct][ct_data[ct]["TF"].isin(ct_top_tfs[ct])].copy()

# Final ordered TG list: global top-N + guaranteed FORCE_TF targets
selected_tgs = list(dict.fromkeys(
    (gwas_tgs + other_tgs)[:N_TG_COLS]
    + [t for t in ranked_tgs if t in selected_tgs_set and t not in (gwas_tgs + other_tgs)[:N_TG_COLS]]
))

# Waterfall ordering: GWAS TGs first, then non-GWAS grouped by dominant CT
def dominant_ct(tg):
    best_ct, best_w = CELL_TYPES[0], -1.0
    for ct in CELL_TYPES:
        w = tg_ct_w_full.get(tg, {}).get(ct, 0.0)
        if w > best_w:
            best_ct, best_w = ct, w
    return best_ct

gwas_cols  = [t for t in selected_tgs if t.upper() in AD_GWAS]
other_cols = [t for t in selected_tgs if t.upper() not in AD_GWAS]
other_cols_sorted = sorted(
    other_cols,
    key=lambda t: (CELL_TYPES.index(dominant_ct(t)),
                   -tg_ct_w_full.get(t, {}).get(dominant_ct(t), 0.0))
)
final_tgs = gwas_cols + other_cols_sorted
n_cols = len(final_tgs)

# ── Build combined matrix ─────────────────────────────────────────────────────
# Rows: flattened list of (ct, tf) in CT order; cols: final_tgs
row_labels = []   # TF names for y-axis
row_ct     = []   # which CT each row belongs to
for ct in CELL_TYPES:
    for tf in ct_top_tfs[ct]:
        row_labels.append(tf)
        row_ct.append(ct)

n_rows = len(row_labels)
M = np.zeros((n_rows, n_cols))

# Build lookup: (ct, tf, tg) -> weight
lookup = {}
for ct in CELL_TYPES:
    for _, r in ct_data[ct].iterrows():
        lookup[(ct, r["TF"], r["TG"])] = float(r["weight_combined"])

for i, (tf, ct) in enumerate(zip(row_labels, row_ct)):
    for j, tg in enumerate(final_tgs):
        M[i, j] = lookup.get((ct, tf, tg), 0.0)

# ── Colour scaling ────────────────────────────────────────────────────────────
flat = M[M > 0].flatten() if (M > 0).any() else np.array([0.0, 1.0])
vmin = float(np.quantile(flat, 0.05))
vmax = float(np.quantile(flat, 0.95))
if vmax <= vmin:
    vmin, vmax = flat.min(), flat.max()
norm = PowerNorm(gamma=HEATMAP_GAMMA, vmin=vmin, vmax=vmax)

# ── Figure layout ─────────────────────────────────────────────────────────────
# Main heatmap + left color strip (narrow axis) + colorbar
fig = plt.figure(figsize=(14, 9))

# Axes: [left, bottom, width, height] in figure fraction
strip_w = 0.018   # color strip
gap     = 0.005
left0   = 0.13    # space for row labels
hm_w    = 0.72
hm_h    = 0.70
bottom0 = 0.22    # space for rotated column labels + legend

ax_strip = fig.add_axes([left0, bottom0, strip_w, hm_h])
ax_hm    = fig.add_axes([left0 + strip_w + gap, bottom0, hm_w, hm_h])
ax_cb    = fig.add_axes([left0 + strip_w + gap + hm_w + 0.015, bottom0, 0.018, hm_h])

# ── Draw color strip ──────────────────────────────────────────────────────────
ct_boundaries = []   # (start_row, end_row, ct)
i = 0
for ct in CELL_TYPES:
    n = N_TF_PER_CT
    ct_boundaries.append((i, i + n, ct))
    i += n

for start, end, ct in ct_boundaries:
    ax_strip.add_patch(mpatches.Rectangle(
        (0, n_rows - end), 1, end - start,
        facecolor=CT_COLORS[ct], edgecolor="white", linewidth=0.5,
    ))
ax_strip.set_xlim(0, 1)
ax_strip.set_ylim(0, n_rows)
ax_strip.axis("off")

# Cell-type labels on strip (centered vertically on each block, left of strip)
for start, end, ct in ct_boundaries:
    mid_y = (n_rows - start + n_rows - end) / 2   # center in imshow coordinates
    fig.text(
        left0 - 0.005, bottom0 + hm_h * (mid_y / n_rows),
        CT_LABELS[ct], ha="right", va="center",
        fontsize=8, fontweight="bold", color=CT_COLORS[ct],
        transform=fig.transFigure,
    )

# ── Draw heatmap ──────────────────────────────────────────────────────────────
im = ax_hm.imshow(M, aspect="auto", cmap=HEATMAP_CMAP, norm=norm,
                  interpolation="nearest")

# Row (TF) labels
ax_hm.set_yticks(range(n_rows))
ax_hm.set_yticklabels(row_labels, fontsize=7.5, fontfamily="Arial")
# Bold + GWAS frame for TFs that are GWAS genes
gwas_bbox = dict(boxstyle="round,pad=0.12", facecolor="none",
                 edgecolor=GWAS_COLOR, linewidth=1.0)
for lbl, tf in zip(ax_hm.get_yticklabels(), row_labels):
    if tf.upper() in AD_GWAS:
        lbl.set_fontweight("bold")
        lbl.set_bbox(gwas_bbox)

# Column (TG) labels — 45°, colored by GWAS status
ax_hm.set_xticks(range(n_cols))
ax_hm.set_xticklabels(final_tgs, rotation=45, ha="right",
                       fontsize=7.5, fontfamily="Arial")
for lbl, tg in zip(ax_hm.get_xticklabels(), final_tgs):
    if tg.upper() in AD_GWAS:
        lbl.set_fontweight("bold")
        lbl.set_color(GWAS_COLOR)

ax_hm.tick_params(axis="both", which="both", length=0)

# Horizontal separator lines between cell-type blocks
for start, end, ct in ct_boundaries[:-1]:
    ax_hm.axhline(end - 0.5, color="white", linewidth=1.5, zorder=3)

# Thin border
for spine in ax_hm.spines.values():
    spine.set_visible(False)

# ── Colorbar ──────────────────────────────────────────────────────────────────
cb = fig.colorbar(im, cax=ax_cb)
cb.set_label("NetRem weight\n(PowerNorm γ=0.65)", fontsize=7.5)
cb.ax.tick_params(labelsize=7)

# ── Legend ────────────────────────────────────────────────────────────────────
ct_patches = [mpatches.Patch(facecolor=CT_COLORS[ct], label=CT_LABELS[ct].replace("\n"," "))
              for ct in CELL_TYPES]
gwas_handle = Line2D([0],[0], marker="s", color="none",
                     markerfacecolor="none", markeredgecolor=GWAS_COLOR,
                     markeredgewidth=1.2, markersize=8, label="AD GWAS gene")
fig.legend(handles=ct_patches + [gwas_handle],
           loc="lower center", ncol=6, frameon=True,
           fontsize=8, bbox_to_anchor=(0.5, 0.01))

fig.suptitle(
    "Combined TF → Target Gene regulatory weights across cell types\n"
    "(top-6 TFs per cell type · top-25 TGs by GWAS priority + cross-cell-type sharing)",
    fontsize=9, y=0.97,
)

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = OUT_DIR / "combined_heatmap_all_celltypes.pdf"
out_png = OUT_DIR / "combined_heatmap_all_celltypes.png"
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved:\n  {out_pdf}\n  {out_png}")
print(f"\nMatrix: {n_rows} TF rows × {n_cols} TG columns")
print(f"GWAS TGs shown: {[t for t in final_tgs if t.upper() in AD_GWAS]}")
print(f"Cell-type-specific non-GWAS TGs: {[t for t in final_tgs if t.upper() not in AD_GWAS][:15]}")
