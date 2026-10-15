#!/usr/bin/env python3
"""
fig1b_umap_substates.py

Fig. 1b: UMAP of the final atlas coloured by the 33 fine-grained cell
substates (EXC, INH, ODC, ASC, MIC, OPC, ENDO lineages). Labels are placed
per connected UMAP blob; colours are grouped by major cell type and shaded
within that type's hue family (with a slight global warm tint).

Input (working directory):
  13ad_3m_5k_atlas_final_renamed.h5ad   (read with h5py: obsm X_umap + substate labels)

Output:
  figures/endo_recomputed_umap/umap_subtypes_labeled_v10.{png,pdf}
"""

import colorsys
import warnings
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.rcParams["pdf.fonttype"]    = 42
matplotlib.rcParams["ps.fonttype"]     = 42
matplotlib.rcParams["font.family"]     = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb, to_hex, to_rgba
from scipy import ndimage
from collections import defaultdict, OrderedDict
from pathlib import Path
from datetime import datetime
warnings.filterwarnings("ignore")

BASE    = Path.cwd()  # run from the directory holding the atlas h5ad
H5AD    = BASE / "13ad_3m_5k_atlas_final_renamed.h5ad"
FIG_DIR = BASE / "figures" / "endo_recomputed_umap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OBS_COL        = "subcluster_annot_labeled_renamed"
MIN_CELLS_BLOB = 300
GRID_SIZE      = 150
POINT_SIZE     = 0.07
FONT_SIZE      = 6.5

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── 1. Load data ──────────────────────────────────────────────────────────────
log(f"Loading {H5AD.name}...")
with h5py.File(H5AD, "r") as f:
    umap     = f["obsm"]["X_umap"][:].astype(np.float32)
    cats_raw = [c.decode() if isinstance(c, bytes) else str(c)
                for c in f["obs"][OBS_COL]["categories"][:]]
    codes    = f["obs"][OBS_COL]["codes"][:]
    labels   = np.array(cats_raw)[codes]

log(f"  {len(labels):,} cells  |  {len(cats_raw)} subtypes")

counts = pd.Series(labels).value_counts()
pct    = (counts / len(labels) * 100)


# ── 2. Colour palette ─────────────────────────────────────────────────────────
prefix_base = {
    "EXC":           "#1a5fb4",   # royal blue    (~213°)
    "INH":           "#c2185b",   # hot pink      (~336°)
    "ODC":           "#2e7d32",   # forest green  (~125°)
    "OPC":           "#4527a0",   # deep indigo   (~271°)
    "ASC":           "#00838f",   # cyan-teal     (~185°)
    "MIC":           "#b86e00",   # golden amber  (~46°)
    "ENDO_dropped":  "#78909c",   # steel grey
    "ENDO":          "#b71c1c",   # dark crimson  (~0°)
}

# Per-group overrides (missing keys fall back to the function defaults)
prefix_direction = {
    "EXC": "up",        # darkest anchored at base blue, no teal bleed at dark end
    "INH": "symmetric", # centred at hot pink; [309°–3°] = pink family throughout
    "ODC": "down",      # darkest drifts toward olive (~71°), away from ASC teal
    # ASC, MIC, OPC, ENDO: "symmetric" (default)
}
# Per-group lightness floor: raise lo so dark shades are brighter/more visible
prefix_lo = {
    "EXC": 0.30,   # EXC_IT → bright royal blue, not near-black navy
    "ASC": 0.28,   # darkest ASC shade remains visibly teal, not dark green
}
# Per-group hue spread
prefix_spread = {
    "EXC": 0.09,   # [213°, 245°] — all solidly blue
    "ASC": 0.10,   # [167°, 203°] — teal-turquoise-cyan; spread < 0.15 prevents green bleed
}

def vary_hue_and_lightness(base_hex, n, lo=0.18, hi=0.82,
                            hue_spread=0.15, hue_direction="symmetric"):
    """
    Generate n shades dark (lo) → light (hi), also sweeping hue.
    hue_direction controls which end of the sweep is darkest:
      "symmetric": [base-spread/2 … base+spread/2]
      "up":        [base          … base+spread  ]   darkest = base hue
      "down":      [base-spread   … base         ]   lightest = base hue
    """
    if n == 0:
        return []
    r, g, b = to_rgb(base_hex)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if n == 1:
        return [to_hex(colorsys.hls_to_rgb(h, (lo + hi) / 2, min(s * 1.05, 0.97)))]
    ls = np.linspace(lo, hi, n)
    if hue_direction == "up":
        hs = np.linspace(h, h + hue_spread, n) % 1.0
    elif hue_direction == "down":
        hs = np.linspace(h - hue_spread, h, n) % 1.0
    else:  # symmetric
        hs = np.linspace(h - hue_spread / 2, h + hue_spread / 2, n) % 1.0
    return [to_hex(colorsys.hls_to_rgb(hv, lv, min(s * 1.05, 0.97)))
            for hv, lv in zip(hs, ls)]

WARMTH = 0.07   # 0 = no change, 1 = fully warm (red/no-blue)

def apply_warmth(hex_color, warmth=WARMTH):
    """Nudge R up and B down to add a subtle warm tint without changing lightness."""
    r, g, b = to_rgb(hex_color)
    r = min(1.0, r + warmth * (1.0 - r))
    b = max(0.0, b - warmth * b)
    return to_hex((r, g, b))

cats_sorted = sorted(set(labels))
groups = defaultdict(list)
for c in cats_sorted:
    matched = None
    for p in sorted(prefix_base, key=len, reverse=True):
        if c.startswith(p):
            matched = p
            break
    groups[matched if matched else "other"].append(c)

palette = {}
for prefix, members in groups.items():
    base      = prefix_base.get(prefix, "#bcbd22")
    direction = prefix_direction.get(prefix, "symmetric")
    spread    = prefix_spread.get(prefix, 0.15)
    lo        = prefix_lo.get(prefix, 0.18)
    shades    = vary_hue_and_lightness(base, len(members),
                                       lo=lo, hue_spread=spread, hue_direction=direction)
    # Apply warm tint to every group except the neutral grey
    if prefix != "ENDO_dropped":
        shades = [apply_warmth(s) for s in shades]
    members_by_count = sorted(members, key=lambda c: -counts.get(c, 0))
    for m, shade in zip(members_by_count, shades):
        palette[m] = shade


# ── 3. Scatter plot ───────────────────────────────────────────────────────────
log("Plotting scatter...")
color_arr = np.array([to_rgba(palette[c]) for c in labels])
order     = np.argsort([counts[c] for c in labels])[::-1]

fig, ax = plt.subplots(figsize=(15, 12))
ax.scatter(
    umap[order, 0], umap[order, 1],
    c=color_arr[order],
    s=POINT_SIZE, alpha=0.55, rasterized=True, linewidths=0,
)


# ── 4. Per-blob label placement ───────────────────────────────────────────────
log("Finding blobs and placing labels...")

x_min, x_max = float(umap[:, 0].min()), float(umap[:, 0].max())
y_min, y_max = float(umap[:, 1].min()), float(umap[:, 1].max())

def get_blob_centers(xy_sub, min_cells=MIN_CELLS_BLOB, grid=GRID_SIZE):
    h, xedges, yedges = np.histogram2d(
        xy_sub[:, 0], xy_sub[:, 1],
        bins=grid,
        range=[[x_min, x_max], [y_min, y_max]],
    )
    labeled, n_comp = ndimage.label(h > 0)
    blobs = []
    for cid in range(1, n_comp + 1):
        rows, cols = np.where(labeled == cid)
        w = h[rows, cols]
        n = int(w.sum())
        cx = float(np.average((xedges[rows] + xedges[rows + 1]) / 2.0, weights=w))
        cy = float(np.average((yedges[cols] + yedges[cols + 1]) / 2.0, weights=w))
        blobs.append((n, cx, cy))
    if not blobs:
        return [(float(np.median(xy_sub[:, 0])), float(np.median(xy_sub[:, 1])))]
    blobs.sort(key=lambda x: -x[0])
    centers = [(blobs[0][1], blobs[0][2])]
    for n, cx, cy in blobs[1:]:
        if n >= min_cells:
            centers.append((cx, cy))
    return centers

n_labels = 0
for subtype in cats_sorted:
    mask = labels == subtype
    if mask.sum() == 0:
        continue
    centers = get_blob_centers(umap[mask])
    for cx, cy in centers:
        ax.text(cx, cy, subtype,
                fontsize=FONT_SIZE, fontweight="bold",
                ha="center", va="center", color="black")
        n_labels += 1

log(f"  {n_labels} labels placed across {len(cats_sorted)} subtypes")


# ── 5. Grouped legend ─────────────────────────────────────────────────────────
GROUP_ORDER = ["EXC", "INH", "ODC", "ASC", "MIC", "OPC", "ENDO", "ENDO_dropped"]
GROUP_DISPLAY = {
    "EXC":          "Excitatory",
    "INH":          "Inhibitory",
    "ODC":          "Oligodendrocyte",
    "ASC":          "Astrocyte",
    "MIC":          "Microglia",
    "OPC":          "OPC",
    "ENDO":         "Endothelial",
    "ENDO_dropped": "Endothelial (excluded)",
}

def get_group(cat):
    for g in sorted(GROUP_ORDER, key=len, reverse=True):
        if cat.startswith(g):
            return g
    return "other"

grouped = OrderedDict()
for g in GROUP_ORDER:
    grouped[g] = []
for cat in cats_sorted:
    g = get_group(cat)
    if g in grouped:
        grouped[g].append(cat)
    else:
        grouped.setdefault("other", []).append(cat)

legend_handles = []
legend_labels  = []
for group, members in grouped.items():
    if not members:
        continue
    legend_handles.append(Patch(facecolor="none", edgecolor="none"))
    legend_labels.append(GROUP_DISPLAY.get(group, group))
    for cat in sorted(members, key=lambda c: -counts.get(c, 0)):
        legend_handles.append(Patch(facecolor=palette[cat], edgecolor="none"))
        legend_labels.append(f"  {cat}  ({pct[cat]:.1f}%, n={counts[cat]:,})")

leg = ax.legend(
    legend_handles, legend_labels,
    loc="center left", bbox_to_anchor=(1.02, 0.5),
    fontsize=7, frameon=False,
    title="Cell subtype", title_fontsize=9,
    handlelength=1.2, handleheight=1.0,
)
for text, handle in zip(leg.get_texts(), leg.legend_handles):
    if handle.get_facecolor()[3] == 0:
        text.set_fontweight("bold")
        text.set_fontsize(8)

# Palette log
log("Palette (dark → light, most → least abundant):")
for g in GROUP_ORDER:
    members = [c for c in cats_sorted if get_group(c) == g]
    dir_label = prefix_direction.get(g, "symmetric")
    log(f"  [{g} / {dir_label}]")
    for cat in sorted(members, key=lambda c: -counts.get(c, 0)):
        log(f"    {cat:40s}  {palette[cat]}  (n={counts.get(cat,0):,})")


# ── 6. Axes ───────────────────────────────────────────────────────────────────
ax.set_title(f"Cell subtypes  (n={len(labels):,})", fontsize=12, pad=6)
ax.set_xlabel("UMAP_1")
ax.set_ylabel("UMAP_2")
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
plt.tight_layout()


# ── 7. Save ───────────────────────────────────────────────────────────────────
for ext in ("png", "pdf"):
    out = FIG_DIR / f"umap_subtypes_labeled_v10.{ext}"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    log(f"  saved {out}")
plt.close(fig)
log("Done.")
