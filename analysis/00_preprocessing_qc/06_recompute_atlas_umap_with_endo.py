"""
recompute_atlas_umap_with_endo.py

Recomputes the UMAP on all 3,183,466 cells (3.1M atlas + 37k ENDO) using
their X_scANVI (10D) latent coordinates.

Background
----------
The existing atlas UMAP was built on 3,146,108 cells with cluster 13
(Endothelial, n=319) as the only endothelial reference.  With only 319 cells,
no proper endothelial island formed, and adding 37k ENDO cells via kNN coord
transfer placed them on top of adjacent Excitatory clusters.

This script adds the 37k ENDO cells to the UMAP from scratch, giving them
first-class membership.  X_scANVI is already patched for ENDO rows in the
input h5ad (from project_endo_scanvi_umap.py attempt 12/13).

Parameters match the original atlas UMAP exactly:
  n_neighbors=15, metric="cosine", random_state=0, use_rep="X_scANVI"
  UMAP a=0.583030, b=1.334167  (min_dist ≈ 0.3, spread = 1.0)

Input:
  h5ad/13ad_3m_5k_atlas_with_endo_scanvi.h5ad   (3,183,466 cells, X_scANVI complete)

Output:
  h5ad/13ad_3m_5k_atlas_with_endo_umap_recomputed.h5ad   (new X_umap for all cells)
  figures/endo_recomputed_umap/umap_celltype.png/pdf
  figures/endo_recomputed_umap/umap_leiden.png/pdf
  figures/endo_recomputed_umap/umap_origin.png/pdf
  reports/endo_recomputed_umap_summary.md
"""

import os, sys, shutil
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import h5py
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba

BASE     = os.getcwd()  # Set to project root (containing h5ad/, reports/) before running
IN_H5AD  = f"{BASE}/h5ad/13ad_3m_5k_atlas_with_endo_scanvi.h5ad"
OUT_H5AD = f"{BASE}/h5ad/13ad_3m_5k_atlas_with_endo_umap_recomputed.h5ad"
FIGS     = f"{BASE}/figures/endo_recomputed_umap"
REPORT   = f"{BASE}/reports/endo_recomputed_umap_summary.md"

N_ATLAS    = 3_146_108
N_TOTAL    = 3_183_466
UMAP_A     = 0.583030019901822
UMAP_B     = 1.3341669931033755

os.makedirs(FIGS, exist_ok=True)

_log = []
def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    _log.append(line)


# ── 1. Load X_scANVI + obs metadata ──────────────────────────────────────────
log("Loading X_scANVI for all cells...")
with h5py.File(IN_H5AD, "r") as f:
    X_scanvi = f["obsm/X_scANVI"][:].astype(np.float32)   # (3183466, 10)
    assert not np.isnan(X_scanvi).any(), "NaN in X_scANVI — ENDO rows not patched"
    log(f"  shape: {X_scanvi.shape} | range [{X_scanvi.min():.3f}, {X_scanvi.max():.3f}]")

    leiden_cats  = f["obs/leiden/categories"][:].astype(str)
    leiden_codes = f["obs/leiden/codes"][:]
    leiden_vals  = leiden_cats[leiden_codes]

    ct_cats  = f["obs/Celltype_markers/categories"][:].astype(str)
    ct_codes = f["obs/Celltype_markers/codes"][:]
    ct_vals  = ct_cats[ct_codes]

    orig_cats  = f["obs/cell_origin/categories"][:].astype(str)
    orig_codes = f["obs/cell_origin/codes"][:]
    origin_vals = orig_cats[orig_codes]

log(f"  ENDO (leiden==18): {(leiden_vals=='18').sum():,} cells")
log(f"  NaN in ENDO X_scANVI: {np.isnan(X_scanvi[leiden_vals=='18']).any()}")


# ── 2. Build minimal AnnData for neighbors + UMAP ────────────────────────────
log("Building AnnData (obs only + X_scANVI)...")
obs_df = pd.DataFrame({
    "leiden":           leiden_vals,
    "Celltype_markers": ct_vals,
    "cell_origin":      origin_vals,
}, index=[str(i) for i in range(N_TOTAL)])

adata = ad.AnnData(obs=obs_df)
adata.obsm["X_scANVI"] = X_scanvi
del X_scanvi   # free memory; adata holds the reference

log(f"  AnnData: {adata.n_obs:,} cells")


# ── 3. Build neighbor graph (same params as original atlas) ───────────────────
log("Building neighbor graph (n_neighbors=15, metric=cosine, use_rep=X_scANVI)...")
log("  This will take 20-40 min on the full 3.18M cell set...")
sc.pp.neighbors(
    adata,
    use_rep="X_scANVI",
    n_neighbors=15,
    metric="cosine",
    random_state=0,
    n_pcs=None,
)
log("  Neighbor graph done.")


# ── 4. Run UMAP (same a/b as original atlas) ──────────────────────────────────
log("Running UMAP (a=0.583, b=1.334, same as original atlas)...")
log("  This will take 30-60 min on 3.18M cells...")
sc.tl.umap(adata, a=UMAP_A, b=UMAP_B, random_state=0)
log("  UMAP done.")

X_umap_new = adata.obsm["X_umap"].astype(np.float32)
log(f"  Full UMAP range: x=[{X_umap_new[:,0].min():.2f}, {X_umap_new[:,0].max():.2f}]"
    f"  y=[{X_umap_new[:,1].min():.2f}, {X_umap_new[:,1].max():.2f}]")

mask18 = (leiden_vals == "18")
log(f"  ENDO centroid: ({X_umap_new[mask18,0].mean():.2f}, {X_umap_new[mask18,1].mean():.2f})")
log(f"  ENDO x: [{X_umap_new[mask18,0].min():.2f}, {X_umap_new[mask18,0].max():.2f}]")
log(f"  ENDO y: [{X_umap_new[mask18,1].min():.2f}, {X_umap_new[mask18,1].max():.2f}]")


# ── 5. Copy input h5ad and overwrite X_umap ───────────────────────────────────
log(f"Copying input h5ad → {OUT_H5AD}")
shutil.copy2(IN_H5AD, OUT_H5AD)

log("Overwriting X_umap with recomputed embedding...")
with h5py.File(OUT_H5AD, "r+") as f:
    del f["obsm/X_umap"]
    f["obsm/X_umap"] = X_umap_new
log(f"Saved: {OUT_H5AD}")


# ── 6. Verify ─────────────────────────────────────────────────────────────────
log("Verifying...")
with h5py.File(OUT_H5AD, "r") as f:
    stored = f["obsm/X_umap"][:]
    nan_ok  = not np.isnan(stored).any()
    shape_ok = stored.shape == (N_TOTAL, 2)
    log(f"  Shape OK ({stored.shape}): {shape_ok}")
    log(f"  NaN-free: {nan_ok}")
    log(f"  ENDO (last {N_TOTAL - N_ATLAS:,}) centroid: "
        f"({stored[N_ATLAS:,0].mean():.2f}, {stored[N_ATLAS:,1].mean():.2f})")


# ── 7. Figures ────────────────────────────────────────────────────────────────
PALETTE = {
    "Excitatory":      "#1f77b4",
    "Inhibitory":      "#ff7f0e",
    "Oligodendrocyte": "#2ca02c",
    "Astrocyte":       "#d62728",
    "Microglia":       "#9467bd",
    "OPC":             "#8c564b",
    "Endothelial":     "#e377c2",
    "Unknown":         "#cccccc",
}

def umap_plot(xy, cats, title, fname, palette, ps=0.07, lfs=9):
    counts    = pd.Series(cats).value_counts()
    cat_names = counts.index.tolist()
    pct       = (counts / len(cats) * 100).round(2)
    fb = plt.get_cmap("tab20")
    for i, c in enumerate(cat_names):
        if c not in palette:
            palette[c] = fb(i % 20)
    color_arr = np.array([to_rgba(palette[c]) for c in cats])
    order     = np.argsort([counts[c] for c in cats])[::-1]
    fig, ax   = plt.subplots(figsize=(15, 12))
    ax.scatter(xy[order, 0], xy[order, 1],
               c=color_arr[order], s=ps, alpha=0.5, rasterized=True, linewidths=0)
    df = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "label": cats})
    for name, row in df.groupby("label", observed=True)[["x","y"]].median().iterrows():
        t = ax.text(row["x"], row["y"], str(name), fontsize=lfs,
                    fontweight="bold", ha="center", va="center", color="black")
        t.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground="white")])
    handles = [Patch(facecolor=palette[c], edgecolor="none",
                     label=f"{c}  ({pct[c]:.2f}%, n={counts[c]:,})")
               for c in cat_names]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, frameon=False, title=title, title_fontsize=9)
    ax.set_title(f"{title}  (n={len(cats):,})", fontsize=12)
    ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{FIGS}/{fname}.png", dpi=200, bbox_inches="tight")
    plt.savefig(f"{FIGS}/{fname}.pdf",            bbox_inches="tight")
    plt.close()
    log(f"  -> {fname}.png/pdf")

log("Plotting UMAPs...")
umap_plot(X_umap_new, ct_vals, "Cell type", "umap_celltype",
          palette=PALETTE.copy(), lfs=10)

fb = plt.get_cmap("tab20")
leiden_pal = {str(i): fb(i % 20) for i in range(18)}
leiden_pal["18"]      = "#e377c2"
leiden_pal["unknown"] = "#cccccc"
umap_plot(X_umap_new, leiden_vals, "Leiden cluster", "umap_leiden",
          palette=leiden_pal, lfs=7)

umap_plot(X_umap_new, origin_vals, "Cell origin", "umap_origin",
          palette={"atlas": "#dddddd", "rescued_qc": "#e377c2"},
          ps=0.05, lfs=11)


# ── 8. Summary report ─────────────────────────────────────────────────────────
endo_cent_x = X_umap_new[mask18, 0].mean()
endo_cent_y = X_umap_new[mask18, 1].mean()
report = f"""# ENDO Recomputed UMAP Summary

**Generated:** {datetime.now().isoformat(timespec='seconds')}

## Approach
Full UMAP recompute on all {N_TOTAL:,} cells (atlas + ENDO) from X_scANVI (10D).
ENDO cells participate as first-class members — not placed via kNN coord transfer.

## Parameters (match original atlas exactly)
- n_neighbors=15, metric=cosine, random_state=0
- UMAP a={UMAP_A}, b={UMAP_B}

## Input
`{IN_H5AD}`
- X_scANVI: 10D for all {N_TOTAL:,} cells (ENDO rows patched by attempt 12/13)

## Output
`{OUT_H5AD}`

## ENDO cells (leiden==18, n={mask18.sum():,})
- Centroid: ({endo_cent_x:.2f}, {endo_cent_y:.2f})
- x: [{X_umap_new[mask18,0].min():.2f}, {X_umap_new[mask18,0].max():.2f}]
- y: [{X_umap_new[mask18,1].min():.2f}, {X_umap_new[mask18,1].max():.2f}]

## Log
```
{chr(10).join(_log)}
```
"""
with open(REPORT, "w") as fh:
    fh.write(report)
log(f"Report written: {REPORT}")
log("All done.")
