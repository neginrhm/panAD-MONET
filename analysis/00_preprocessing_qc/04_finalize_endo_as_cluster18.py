"""
finalize_endo_as_cluster18.py

The combined atlas + rescued ENDO h5ad already has:
  - Atlas cells: original X_umap, leiden 0–17, original Celltype_markers
  - Rescued ENDO cells (37,358): projected X_umap, leiden='unknown', Celltype_markers='ENDO'

This script:
  1. Assigns leiden='18' and Celltype_markers='Endothelial' to rescued ENDO cells
  2. Patches these two obs columns in-place in the h5ad using h5py (no full matrix load)
  3. Generates UMAP figures showing the original 17 clusters + ENDO as cluster 18

Output h5ad:  13ad_3m_5k_atlas_with_endo_final.h5ad
Output figs:  figures/atlas_with_endo/
"""

import os
import shutil
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from datetime import datetime

BASE  = os.getcwd()  # Set to project root (containing h5ad/, reports/) before running
SRC   = f"{BASE}/h5ad/13ad_3m_5k_atlas_plus_endo.h5ad"
OUT   = f"{BASE}/h5ad/13ad_3m_5k_atlas_with_endo_final.h5ad"
FIGS  = f"{BASE}/figures/atlas_with_endo"
os.makedirs(FIGS, exist_ok=True)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── Colour palette: same as the rest of the project, ENDO in a warm pink ───────
CELLTYPE_PALETTE = {
    "Excitatory":    "#1f77b4",
    "Inhibitory":    "#ff7f0e",
    "Oligodendrocyte": "#2ca02c",
    "Astrocyte":     "#d62728",
    "Microglia":     "#9467bd",
    "OPC":           "#8c564b",
    "Endothelial":   "#e377c2",
    "Unknown":       "#cccccc",
}

LEIDEN_ENDO_COLOR = "#e377c2"   # same pink for cluster 18


# ── Helpers ─────────────────────────────────────────────────────────────────────

def label_at_centroid(ax, xy, labels, fontsize=8):
    df = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "label": labels})
    centers = df.groupby("label", observed=True)[["x", "y"]].median()
    for name, row in centers.iterrows():
        t = ax.text(row["x"], row["y"], str(name),
                    fontsize=fontsize, fontweight="bold",
                    ha="center", va="center", color="black")
        t.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground="white")])


def umap_plot(xy, cats, title, fname, palette, point_size=0.08, label_fontsize=8):
    n = len(cats)
    counts = pd.Series(cats).value_counts()
    cat_names = counts.index.tolist()
    pct = (counts / n * 100).round(2)

    # fill any missing colours
    fb = plt.get_cmap("tab20")
    for i, c in enumerate(cat_names):
        if c not in palette:
            palette[c] = fb(i % 20)

    color_arr = np.array([to_rgba(palette[c]) for c in cats])
    order = np.argsort([counts[c] for c in cats])[::-1]

    fig, ax = plt.subplots(figsize=(15, 12))
    ax.scatter(xy[order, 0], xy[order, 1],
               c=color_arr[order], s=point_size,
               alpha=0.55, rasterized=True, linewidths=0)
    label_at_centroid(ax, xy, cats, fontsize=label_fontsize)

    handles = [Patch(facecolor=palette[c], edgecolor="none",
                     label=f"{c}  ({pct[c]:.2f}%, n={counts[c]:,})")
               for c in cat_names]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, frameon=False, title=title.split("\n")[0], title_fontsize=9)
    ax.set_title(f"{title}  (n={n:,})", fontsize=12)
    ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{FIGS}/{fname}.png", dpi=200, bbox_inches="tight")
    plt.savefig(f"{FIGS}/{fname}.pdf",            bbox_inches="tight")
    plt.close()
    log(f"  -> {FIGS}/{fname}.png/pdf")


def patch_obs_categorical(h5file, col, new_values):
    """
    Replace an obs categorical column in an open h5py.File (r+ mode).
    Writes the exact encoding anndata 0.2.0 expects:
      group attrs:      encoding-type='categorical', encoding-version='0.2.0', ordered=False
      categories attrs: encoding-type='string-array', encoding-version='0.2.0'
      codes attrs:      encoding-type='array',        encoding-version='0.2.0'
    """
    cats = sorted(set(new_values))
    cat_to_code = {c: i for i, c in enumerate(cats)}
    codes = np.array([cat_to_code[v] for v in new_values],
                     dtype=np.int8 if len(cats) <= 127 else np.int16)

    grp = h5file["obs"]
    if col in grp:
        del grp[col]
    sub = grp.create_group(col)
    sub.attrs["encoding-type"]    = "categorical"
    sub.attrs["encoding-version"] = "0.2.0"
    sub.attrs["ordered"]          = False

    dt = h5py.special_dtype(vlen=str)
    cats_ds = sub.create_dataset("categories", data=np.array(cats, dtype=object), dtype=dt)
    cats_ds.attrs["encoding-type"]    = "string-array"
    cats_ds.attrs["encoding-version"] = "0.2.0"

    codes_ds = sub.create_dataset("codes", data=codes)
    codes_ds.attrs["encoding-type"]    = "array"
    codes_ds.attrs["encoding-version"] = "0.2.0"


# ── Main ─────────────────────────────────────────────────────────────────────────

def main():
    # 1. Load obs + obsm (backed → X stays on disk)
    log(f"Reading obs/obsm from {SRC}")
    a = sc.read_h5ad(SRC, backed="r")
    obs = a.obs.copy()
    xy  = np.asarray(a.obsm["X_umap"])
    a.file.close()
    log(f"  n_obs={len(obs):,}  |  rescued_qc={obs['cell_origin'].eq('rescued_qc').sum():,}")

    # 2. Assign cluster 18 + Endothelial label to rescued ENDO cells
    mask = obs["cell_origin"] == "rescued_qc"

    # Convert to plain strings first to avoid categorical restriction issues
    leiden_str = obs["leiden"].astype(str)
    leiden_str[mask] = "18"
    obs["leiden"] = leiden_str

    ct_str = obs["Celltype_markers"].astype(str)
    ct_str[mask] = "Endothelial"
    # unify remaining ENDO variants in atlas cells too
    ct_str = ct_str.replace({"ENDO": "Endothelial", "ENDO_dropped_p2ry12_pos": "Endothelial"})
    obs["Celltype_markers"] = ct_str

    log("  Labels assigned: leiden='18', Celltype_markers='Endothelial'")

    # 3. Patch the h5ad (copy first, then update two obs columns)
    log(f"Copying {SRC} -> {OUT}")
    shutil.copy2(SRC, OUT)
    log("  Patching obs/leiden and obs/Celltype_markers ...")
    with h5py.File(OUT, "r+") as f:
        patch_obs_categorical(f, "leiden",           obs["leiden"].astype(str).values)
        patch_obs_categorical(f, "Celltype_markers", obs["Celltype_markers"].astype(str).values)
    log(f"  Saved: {OUT}")

    # 4. Verify the patch reads back correctly
    log("  Verifying patch ...")
    chk = sc.read_h5ad(OUT, backed="r")
    n_endo_leiden = (chk.obs["leiden"] == "18").sum()
    n_endo_ct     = (chk.obs["Celltype_markers"] == "Endothelial").sum()
    log(f"  leiden=='18': {n_endo_leiden:,}  |  Celltype_markers=='Endothelial': {n_endo_ct:,}")
    chk.file.close()

    # 5. UMAP: cell type (Celltype_markers)
    log("Plotting cell-type UMAP ...")
    ct = obs["Celltype_markers"].astype(str).values
    umap_plot(xy, ct, "Major cell type", "umap_celltype_with_endo",
              palette=CELLTYPE_PALETTE.copy(), point_size=0.08, label_fontsize=10)

    # 6. UMAP: leiden (0–17 + 18)
    log("Plotting leiden UMAP ...")
    leiden_vals = obs["leiden"].astype(str).values
    # build a palette: clusters 0-17 from tab20, cluster 18 in ENDO pink
    fb = plt.get_cmap("tab20")
    leiden_palette = {}
    for i in range(18):
        leiden_palette[str(i)] = fb(i % 20)
    leiden_palette["18"] = LEIDEN_ENDO_COLOR
    leiden_palette["unknown"] = "#cccccc"
    umap_plot(xy, leiden_vals, "Leiden cluster", "umap_leiden_with_endo",
              palette=leiden_palette, point_size=0.08, label_fontsize=7)

    # 7. UMAP: highlight only ENDO vs rest
    log("Plotting ENDO highlight UMAP ...")
    origin_label = np.where(mask.values, "Endothelial (rescued)", "Original atlas")
    umap_plot(xy, origin_label, "Cell origin", "umap_endo_highlight",
              palette={"Original atlas": "#dddddd", "Endothelial (rescued)": LEIDEN_ENDO_COLOR},
              point_size=0.05, label_fontsize=11)

    log("All done.")
    log(f"  Figures in: {FIGS}/")
    log(f"  Updated h5ad: {OUT}")


if __name__ == "__main__":
    main()
