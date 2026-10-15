"""Merge QC-rescued endothelial cells into the final 3.1M atlas.

Inputs:
  - Atlas: 13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad (3.1M × 5000)
  - Rescued: 13ad_endo_rescued_from_initial_qc.h5ad (93,750 × 25,251) — apply strict QC

Outputs:
  - 13ad_3m_5k_atlas_plus_endo.h5ad   (atlas + ~37k rescued)
  - figures/atlas_plus_endo_umap_*.png
  - reports/merge_endo_summary.md

The atlas is left untouched. New cells get a `cell_origin == 'rescued_qc'` flag.
"""

import os
import json
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.getcwd()  # Set to project root (containing h5ad/, reports/) before running
ATLAS = f"{BASE}/h5ad/13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad"
RESCUED = f"{BASE}/h5ad/13ad_endo_rescued_from_initial_qc.h5ad"
OUT_H5 = f"{BASE}/h5ad/13ad_3m_5k_atlas_plus_endo.h5ad"
OUT_REPORT = f"{BASE}/reports/merge_endo_summary.md"
FIG_DIR = f"{BASE}/figures"

# Strict QC for rescued cells (verified above to yield ~37k vascular cells)
QC_MIN_GENES = 800
QC_MAX_MT = 10.0
QC_MIN_ENDO_MARKERS = 5
QC_MIN_TOTAL_ENDO = 10

stats = {"started": datetime.now().isoformat()}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 1. Load rescued + apply strict QC
# ---------------------------------------------------------------------------
log(f"Loading rescued: {RESCUED}")
r = sc.read_h5ad(RESCUED)
log(f"  rescued shape: {r.shape}")

mask = (
    (r.obs["n_genes_by_counts"].values >= QC_MIN_GENES)
    & (r.obs["pct_counts_MT"].values < QC_MAX_MT)
    & (r.obs["n_endo_markers"].values >= QC_MIN_ENDO_MARKERS)
    & (r.obs["total_endo_count"].values >= QC_MIN_TOTAL_ENDO)
)
log(f"Strict QC: {int(mask.sum()):,} of {r.n_obs:,} pass")
stats["rescued_after_qc"] = int(mask.sum())
r = r[mask].copy()
log(f"  filtered shape: {r.shape}")
log(f"  per-study: {r.obs['Study'].value_counts().to_dict()}")
stats["rescued_per_study"] = r.obs["Study"].value_counts().to_dict()

# ---------------------------------------------------------------------------
# 2. Load atlas (full)
# ---------------------------------------------------------------------------
log(f"Loading atlas: {ATLAS}  (~5 min)")
a = sc.read_h5ad(ATLAS)
log(f"  atlas shape: {a.shape}")
log(f"  atlas layers: {list(a.layers.keys())}")
log(f"  atlas obsm: {list(a.obsm.keys())}")

# ---------------------------------------------------------------------------
# 3. Subset rescued to atlas's 5000 HVGs, in atlas gene order
# ---------------------------------------------------------------------------
log("Subsetting rescued to atlas's 5000 HVGs...")
common_genes = [g for g in a.var_names if g in r.var_names]
log(f"  common genes: {len(common_genes)} of atlas's {a.n_vars}")
stats["common_genes"] = len(common_genes)
missing_in_rescued = [g for g in a.var_names if g not in r.var_names]
log(f"  genes missing from rescued: {len(missing_in_rescued)}")

r = r[:, common_genes].copy()
# If some atlas genes are missing from rescued, pad with zero columns to match atlas gene order exactly
if len(common_genes) < a.n_vars:
    log("  padding rescued with zero columns for atlas-only genes...")
    pad_n = a.n_vars - r.n_vars
    pad = sparse.csr_matrix((r.n_obs, pad_n), dtype=np.float32)
    new_X = sparse.hstack([r.X, pad]).tocsr()
    new_var = pd.concat([r.var, a.var.loc[missing_in_rescued]])
    r = ad.AnnData(X=new_X, obs=r.obs.copy(), var=new_var)
    # Reorder to match atlas exactly
    r = r[:, a.var_names].copy()
log(f"  rescued shape after subset/reorder: {r.shape}")

# ---------------------------------------------------------------------------
# 4. Normalize rescued cells (matches atlas's `normalized` layer)
# ---------------------------------------------------------------------------
log("Normalizing rescued cells (sum=1e4 then log1p)...")
r.layers["raw_count"] = r.X.copy()
sc.pp.normalize_total(r, target_sum=1e4)
sc.pp.log1p(r)
r.layers["normalized"] = r.X.copy()

# ---------------------------------------------------------------------------
# 5. Project rescued cells onto atlas's UMAP via PCA + kNN (manual, avoids
#    scanpy ingest which doesn't handle sparse query data robustly)
# ---------------------------------------------------------------------------
from sklearn.decomposition import TruncatedSVD
import pynndescent

a.X = a.layers["normalized"]

log("Fitting TruncatedSVD (50 comps) on 200k-cell atlas sample...")
np.random.seed(42)
sample_n = min(200_000, a.n_obs)
sample_idx = np.random.choice(a.n_obs, size=sample_n, replace=False)
sample_idx.sort()
svd = TruncatedSVD(n_components=50, random_state=42, algorithm="randomized")
svd.fit(a.X[sample_idx])
log(f"  explained variance ratio (sum top 50): {svd.explained_variance_ratio_.sum():.3f}")

log("Transforming atlas to new shared PCA (chunked)...")
atlas_pca = np.zeros((a.n_obs, 50), dtype=np.float32)
CHUNK = 200_000
for start in range(0, a.n_obs, CHUNK):
    end = min(start + CHUNK, a.n_obs)
    atlas_pca[start:end] = svd.transform(a.X[start:end])
    if (start // CHUNK) % 4 == 0:
        log(f"  atlas transformed {end:,} / {a.n_obs:,}")

log("Transforming rescued to shared PCA...")
rescued_pca = svd.transform(r.X).astype(np.float32)

log("Building pynndescent index on atlas PCA (slow ~5-10 min)...")
nnd = pynndescent.NNDescent(atlas_pca, n_neighbors=15, random_state=42, n_jobs=-1)
nnd.prepare()
log("  index ready")

log("Querying nearest atlas neighbors for rescued cells...")
neighbors_idx, _ = nnd.query(rescued_pca, k=15)
log(f"  neighbors_idx shape: {neighbors_idx.shape}")

log("Median-projecting UMAP coords from neighbors...")
atlas_umap = a.obsm["X_umap"]
rescued_umap = np.zeros((r.n_obs, 2), dtype=np.float32)
for i in range(r.n_obs):
    rescued_umap[i] = np.median(atlas_umap[neighbors_idx[i]], axis=0)
log("  UMAP projection done")

# Stash in r.obsm
r.obsm["X_pca"] = rescued_pca
r.obsm["X_umap"] = rescued_umap

# Replace atlas X_pca with the new shared PCA for consistency in concatenated file
# (existing X_pca was from scVI / different fit — combined kNN downstream needs consistency)
a.obsm["X_pca_shared"] = atlas_pca
r.obsm["X_pca_shared"] = rescued_pca

# Checkpoint: save the projected rescued AnnData (small file) so if write fails we can resume
CHECKPOINT = OUT_H5.replace(".h5ad", "_rescued_projected_checkpoint.h5ad")
log(f"Checkpoint: writing projected rescued to {CHECKPOINT}")
r.write_h5ad(CHECKPOINT, compression="gzip")
log("  checkpoint written")

# ---------------------------------------------------------------------------
# 6. Build matching obs for rescued (same columns as atlas, filled appropriately)
# ---------------------------------------------------------------------------
log("Aligning obs columns to atlas...")
new_obs = pd.DataFrame(index=r.obs.index)
for col in a.obs.columns:
    if col in r.obs.columns:
        new_obs[col] = r.obs[col].values
    else:
        new_obs[col] = pd.NA

# Fill key annotation columns with the rescue label
ENDO_LABEL = "ENDO_rescued_qc"
for c in [
    "Celltype_markers",
    "subcluster_annot",
    "subcluster_annot_named",
    "subcluster_annot_labeled",
    "subcluster_annot_labeled_renamed",
    "subcluster_annot_labeled_renamed_neurons",
    "subcluster_annot_labeled_renamed_neurons2",
    "label_transfer_annot",
    "predictions_combined",
    "predictions_scanvi",
    "predictions_scanvi_grouped",
    "parent_cell_type",
    "parent_type",
]:
    if c in new_obs.columns:
        new_obs[c] = ENDO_LABEL

new_obs["major_type"] = "VASCULAR"

r_obs = new_obs
r.obs = r_obs

# Atlas + rescued provenance flag
a.obs["cell_origin"] = "atlas"
r.obs["cell_origin"] = "rescued_qc"

# ---------------------------------------------------------------------------
# 7. Concatenate
# ---------------------------------------------------------------------------
log("Concatenating atlas + rescued...")
combined = ad.concat(
    [a, r],
    axis=0,
    join="outer",
    merge="same",
    label=None,
    index_unique=None,
)
# concat doesn't always preserve obsm/var/layers — be explicit
combined.var = a.var.copy()
log(f"combined shape: {combined.shape}")
stats["combined_shape"] = list(combined.shape)
stats["atlas_cells"] = int(a.n_obs)
stats["rescued_cells_added"] = int(r.n_obs)

# Re-stitch obsm
log("Reattaching obsm (X_umap, X_pca, X_pca_shared)...")
umap = np.vstack([a.obsm["X_umap"], r.obsm["X_umap"]])
combined.obsm["X_umap"] = umap
# atlas X_pca is from scVI (50-d); rescued X_pca is from our SVD fit — different spaces,
# so we provide both: native atlas X_pca (rescued cells get NaN) and X_pca_shared (both)
atlas_xpca_padded = np.vstack(
    [a.obsm["X_pca"], np.full((r.n_obs, a.obsm["X_pca"].shape[1]), np.nan, dtype=np.float32)]
)
combined.obsm["X_pca"] = atlas_xpca_padded
combined.obsm["X_pca_shared"] = np.vstack([a.obsm["X_pca_shared"], r.obsm["X_pca_shared"]])

# Re-stitch layers
log("Reattaching layers (normalized, raw_count)...")
combined.layers["normalized"] = sparse.vstack([a.layers["normalized"], r.layers["normalized"]]).tocsr()
combined.layers["raw_count"] = sparse.vstack([a.layers["raw_count"], r.layers["raw_count"]]).tocsr()
combined.X = combined.layers["normalized"]

# ---------------------------------------------------------------------------
# 8. Clean obs columns for h5ad write (no pd.NA in object/categorical cols)
# ---------------------------------------------------------------------------
log("Cleaning obs columns (stringify objects, fill NA)...")
for col in combined.obs.columns:
    s = combined.obs[col]
    if pd.api.types.is_categorical_dtype(s):
        # Add 'unknown' as a category to absorb any NaN
        if s.isna().any():
            s = s.cat.add_categories(["unknown"]).fillna("unknown")
            combined.obs[col] = s
    elif s.dtype == object:
        # Stringify and replace NaN/<NA> with empty string
        combined.obs[col] = s.astype(str).replace({"nan": "", "<NA>": "", "None": ""})
    elif pd.api.types.is_bool_dtype(s):
        # Booleans with NaN end up as object — fill with False
        combined.obs[col] = s.fillna(False).astype(bool)
    elif pd.api.types.is_numeric_dtype(s):
        pass  # NaN ok in numeric
log("  obs cleaned")

# ---------------------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------------------
log(f"Writing {OUT_H5}  (~25 min)")
combined.write_h5ad(OUT_H5, compression="gzip")
log("h5ad written")

# ---------------------------------------------------------------------------
# 9. UMAP figure
# ---------------------------------------------------------------------------
log("Generating UMAP figure...")
os.makedirs(FIG_DIR, exist_ok=True)
xy = combined.obsm["X_umap"]
origin = combined.obs["cell_origin"].astype(str).values
fig, ax = plt.subplots(figsize=(11, 11))
ax.scatter(xy[origin == "atlas", 0], xy[origin == "atlas", 1],
           s=0.05, c="lightgray", alpha=0.4, rasterized=True)
ax.scatter(xy[origin == "rescued_qc", 0], xy[origin == "rescued_qc", 1],
           s=2.0, c="#d62728", label=f"Rescued ENDO (n={(origin=='rescued_qc').sum():,})", rasterized=True)
ax.set_title("Atlas (gray) + QC-rescued endothelial cells (red)")
ax.set_xlabel("UMAP_1")
ax.set_ylabel("UMAP_2")
ax.set_aspect("equal")
ax.legend(loc="best", markerscale=5)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/atlas_plus_endo_umap.png", dpi=200, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/atlas_plus_endo_umap.pdf", bbox_inches="tight")
plt.close()
log("  -> atlas_plus_endo_umap.png/pdf")

# ---------------------------------------------------------------------------
# 10. Markdown report
# ---------------------------------------------------------------------------
stats["finished"] = datetime.now().isoformat()
log(f"Writing report: {OUT_REPORT}")
lines = []
lines.append("# Atlas + ENDO Rescue (Merged) Summary")
lines.append(f"\n**Generated:** {stats['finished']}\n")
lines.append("## Inputs")
lines.append(f"- Atlas: `{ATLAS}` ({stats['atlas_cells']:,} cells)")
lines.append(f"- Rescued source: `{RESCUED}` ({stats.get('rescued_after_qc',0):,} after strict QC)")
lines.append("\n## Strict QC applied to rescued")
lines.append(f"- n_genes_by_counts ≥ {QC_MIN_GENES}")
lines.append(f"- pct_counts_MT < {QC_MAX_MT}")
lines.append(f"- n_endo_markers ≥ {QC_MIN_ENDO_MARKERS}")
lines.append(f"- total_endo_count ≥ {QC_MIN_TOTAL_ENDO}")
lines.append(f"- → **{stats['rescued_after_qc']:,} cells passed**")
lines.append("\n## Per-study rescued cells added")
lines.append("| Study | Cells |")
lines.append("|---|---:|")
for k, v in sorted(stats["rescued_per_study"].items(), key=lambda x: -x[1]):
    lines.append(f"| {k} | {v:,} |")
lines.append(f"\n## Combined dataset")
lines.append(f"- Final shape: {stats['combined_shape']}")
lines.append(f"- Output h5ad: `{OUT_H5}`")
lines.append(f"\n## How rescued cells are labeled in the merged file")
lines.append(f"- `cell_origin`: `'atlas'` for original 3.1M, `'rescued_qc'` for added cells")
lines.append(f"- `major_type`: `'VASCULAR'` for rescued")
lines.append(f"- `subcluster_annot_labeled_renamed` and related: `'{ENDO_LABEL}'` for rescued")
lines.append(f"- UMAP/PCA coords: projected via sc.tl.ingest using kNN to atlas neighbors")
lines.append("\n## Notes")
lines.append("- The original 3.1M atlas labels are unchanged. Adding cells, not relabeling.")
lines.append("- UMAP coords for rescued cells are projected (not native) — they sit at their nearest neighbors' positions in the atlas.")
lines.append("- For DEG / hdWGCNA, you can subset on `cell_origin == 'rescued_qc'` for ENDO-only analyses.")
lines.append("\n## Stats")
lines.append("```json")
lines.append(json.dumps(stats, indent=2, default=str))
lines.append("```")

with open(OUT_REPORT, "w") as f:
    f.write("\n".join(lines))
log("Report written")
log("DONE")
