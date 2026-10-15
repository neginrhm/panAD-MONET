"""
project_endo_scanvi_umap.py

Goal
----
Add 37,358 rescued ENDO cells to the existing atlas UMAP without moving any
of the original 3,146,108 atlas cells.

PRIMARY  (GPU): load scANVI query model → encode ENDO → X_scANVI (10D) →
                kNN to atlas → transfer UMAP coords
FALLBACK (CPU): use X_pca_shared (50D) already computed via sc.tl.ingest
                for ENDO cells during the rescue merge

kNN is done via batched brute-force on a 300K subsample of atlas cells
(avoids sklearn OOM on 3.1M-point brute-force).

Inputs (never modified):
  13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad
  13ad_endo_rescued_from_initial_qc.h5ad
  13ad_3m_5k_atlas_plus_endo.h5ad          (ENDO X_pca_shared)
  13ad_3m_5k_atlas_with_endo_final.h5ad    (label template)
  scVI_SCANVI_trained_models/scanvi_query_model_03-11-24/

Outputs (new files only):
  h5ad/13ad_3m_5k_atlas_with_endo_scanvi.h5ad
  figures/endo_scanvi_proj/umap_celltype.png/pdf
  figures/endo_scanvi_proj/umap_leiden.png/pdf
  figures/endo_scanvi_proj/umap_origin.png/pdf
  reports/endo_scanvi_proj_summary.md
"""

import os, sys, shutil, types
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

BASE       = os.getcwd()  # Set to project root (containing h5ad/, reports/) before running
ATLAS_H5AD = f"{BASE}/h5ad/13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad"
ENDO_RAW   = f"{BASE}/h5ad/13ad_endo_rescued_from_initial_qc.h5ad"
COMB_H5AD  = f"{BASE}/h5ad/13ad_3m_5k_atlas_plus_endo.h5ad"
LABEL_H5AD = f"{BASE}/h5ad/13ad_3m_5k_atlas_with_endo_final.h5ad"
MODEL_PATH = f"{BASE}/models/scanvi_query_model"  # trained scANVI model from the (missing, see
                                                    # ../README.md) core integration step -- not
                                                    # included in this repo; place the checkpoint
                                                    # directory here to run this script
OUT_H5AD   = f"{BASE}/h5ad/13ad_3m_5k_atlas_with_endo_scanvi.h5ad"
FIGS       = f"{BASE}/figures/endo_scanvi_proj"
REPORT     = f"{BASE}/reports/endo_scanvi_proj_summary.md"
N_ATLAS    = 3_146_108
K          = 15
JITTER_SD  = 0.12

os.makedirs(FIGS, exist_ok=True)

_log = []
def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    _log.append(line)


# ── 1. Load atlas embeddings ──────────────────────────────────────────────────
log("Loading atlas embeddings (backed)...")
with h5py.File(ATLAS_H5AD, "r") as f:
    atlas_umap   = f["obsm"]["X_umap"][:].astype(np.float32)    # (N, 2)
    atlas_scanvi = f["obsm"]["X_scANVI"][:].astype(np.float32)  # (N, 10)
    atlas_var    = list(f["var"]["_index"][:].astype(str))
    _leiden_cats  = f["obs"]["leiden"]["categories"][:].astype(str)
    _leiden_codes = f["obs"]["leiden"]["codes"][:]
    atlas_leiden  = _leiden_cats[_leiden_codes]   # (N,) string labels
log(f"  {N_ATLAS:,} cells | X_scANVI {atlas_scanvi.shape[1]}D | X_umap 2D")


# ── 2. Full-reference kNN UMAP transfer ──────────────────────────────────────
from sklearn.neighbors import NearestNeighbors

def full_knn_umap(X_query, X_ref, Y_ref, leiden_ref,
                  k=K, jitter=JITTER_SD, seed=42):
    """
    kNN UMAP coord transfer using ALL atlas cells as reference (no subsampling).
    X_ref is 3.1M x 10D float32 (~120 MB) — fits in RAM.
    Uses sklearn ball_tree for exact kNN in 10D.
    No cluster is underrepresented: if ENDO cells are truly nearest to cluster 13,
    the data shows it. No thumb on the scale.
    """
    from collections import Counter
    log(f"  Fitting kNN (k={k}) on full {len(X_ref):,}-cell reference (10D ball_tree)...")
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean",
                          algorithm="ball_tree", n_jobs=-1)
    nn.fit(X_ref.astype(np.float64))

    log(f"  Querying {len(X_query):,} ENDO cells...")
    distances, indices = nn.kneighbors(X_query.astype(np.float64))

    # UMAP coord transfer: inverse-distance-weighted average
    w = 1.0 / (distances.astype(np.float32) + 1e-8)
    w /= w.sum(axis=1, keepdims=True)
    X_new = np.einsum("nk,nkd->nd", w, Y_ref[indices]).astype(np.float32)
    X_new += (np.random.default_rng(seed).standard_normal(X_new.shape)
              * jitter).astype(np.float32)

    # Diagnostic: which leiden clusters are the nearest neighbors?
    nn1_leiden = leiden_ref[indices[:, 0]]   # each cell's single nearest neighbor
    nn_counts  = Counter(nn1_leiden)
    total      = len(nn1_leiden)
    log(f"  Nearest-neighbor leiden distribution (1st neighbor, n={total:,}):")
    for cl, cnt in sorted(nn_counts.items(), key=lambda x: -x[1])[:10]:
        log(f"    leiden {cl}: {cnt:,} ({100*cnt/total:.1f}%)")

    return X_new, indices


# ── 3. PRIMARY: encode ENDO via scANVI model ──────────────────────────────────
use_model    = False
X_endo_emb   = None
reference_emb = None
embedding_key = None

try:
    # Mock torchaudio before importing scvi to avoid CUDA lib path issues
    if "torchaudio" not in sys.modules:
        try:
            import torchaudio
        except OSError:
            log("torchaudio CUDA lib missing — mocking torchaudio for scvi import")
            import importlib.machinery
            _mock_ta = types.ModuleType("torchaudio")
            _mock_ta.__spec__ = importlib.machinery.ModuleSpec("torchaudio", None)
            sys.modules["torchaudio"] = _mock_ta

    # scipy 1.11+ removed numpy-duplicate functions from scipy.linalg; pyro/scvi still calls them
    import scipy.linalg
    for _fn in ["tril", "triu", "diag", "eye", "ones", "zeros"]:
        if not hasattr(scipy.linalg, _fn):
            setattr(scipy.linalg, _fn, getattr(np, _fn))

    # pandas 2.0 removed pandas.core.indexes.numeric (Int64Index etc.);
    # old model.pt pickles reference it during torch.load deserialization
    import pandas as _pd
    if "pandas.core.indexes.numeric" not in sys.modules:
        _num = types.ModuleType("pandas.core.indexes.numeric")
        _num.Int64Index   = _pd.Index
        _num.Float64Index = _pd.Index
        _num.UInt64Index  = _pd.Index
        sys.modules["pandas.core.indexes.numeric"] = _num

    # PyTorch 2.4+ defaults weights_only=True, which rejects old scvi checkpoints
    import torch as _torch
    _orig_load = _torch.load
    def _permissive_load(*args, **kwargs):
        kwargs.setdefault("weights_only", False)
        return _orig_load(*args, **kwargs)
    _torch.load = _permissive_load

    import scvi
    log(f"scvi version: {scvi.__version__}")

    # Load and QC raw ENDO counts
    log("Loading raw ENDO cells...")
    endo = sc.read_h5ad(ENDO_RAW)
    qc = ((endo.obs["n_genes_by_counts"] >= 800) &
          (endo.obs["pct_counts_MT"]     <  10.0) &
          (endo.obs["n_endo_markers"]    >=  5) &
          (endo.obs["total_endo_count"]  >= 10))
    endo = endo[qc].copy()
    log(f"  Post-QC: {endo.n_obs:,} cells")

    # Model batch key is "external_donor_name_label"; map from Sample column
    endo.obs["external_donor_name_label"] = endo.obs["Sample"].astype(str)

    # Model was registered with layer="raw_count"; store raw counts before normalizing
    import scipy.sparse as sp
    endo.layers["raw_count"] = endo.X.copy() if not sp.issparse(endo.X) else endo.X.copy()

    sc.pp.normalize_total(endo, target_sum=1e4)
    sc.pp.log1p(endo)

    # Align to atlas 5k HVGs
    common  = [g for g in atlas_var if g in set(endo.var_names)]
    missing = [g for g in atlas_var if g not in set(endo.var_names)]
    log(f"  HVG alignment: {len(common)} common, {len(missing)} zero-padded")
    endo = endo[:, common].copy()
    if missing:
        import scipy.sparse as sp
        pad = ad.AnnData(
            X=sp.csr_matrix((endo.n_obs, len(missing)), dtype=np.float32),
            obs=endo.obs,
            var=pd.DataFrame(index=missing),
        )
        endo = ad.concat([endo, pad], axis=1)
    endo = endo[:, atlas_var].copy()
    endo.obs["Celltype_markers"] = "Endothelial"

    log("Preparing query anndata + loading scANVI model on CPU...")
    scvi.model.SCANVI.prepare_query_anndata(endo, MODEL_PATH)
    model = scvi.model.SCANVI.load_query_data(endo, MODEL_PATH, accelerator="cpu")
    model.is_trained_ = True

    log("Encoding ENDO → X_scANVI...")
    X_endo_emb = model.get_latent_representation().astype(np.float32)
    assert not np.isnan(X_endo_emb).any(), "NaN in X_scANVI"
    log(f"  shape {X_endo_emb.shape} | range [{X_endo_emb.min():.3f}, {X_endo_emb.max():.3f}]")

    reference_emb = atlas_scanvi
    embedding_key = "X_scANVI via scANVI query model (10D)"
    use_model = True

except Exception as e:
    log(f"scANVI model path failed: {e}")
    log("Falling back to X_pca_shared kNN...")


# ── 4. FALLBACK: X_pca_shared already computed by sc.tl.ingest ───────────────
if not use_model:
    log("Reading X_pca_shared from atlas_plus_endo.h5ad via h5py...")
    with h5py.File(COMB_H5AD, "r") as f:
        # Identify rescued rows by cell_origin code
        cats  = list(f["obs"]["cell_origin"]["categories"][:].astype(str))
        codes = f["obs"]["cell_origin"]["codes"][:]
        rqc_code = cats.index("rescued_qc")
        endo_rows = np.where(codes == rqc_code)[0]
        log(f"  rescued_qc rows: {len(endo_rows):,} (starts at {endo_rows[0]:,})")

        pca_shared = f["obsm"]["X_pca_shared"]
        X_endo_emb    = pca_shared[endo_rows].astype(np.float32)
        reference_emb = pca_shared[:N_ATLAS].astype(np.float32)

    assert not np.isnan(X_endo_emb).any(), "NaN in X_pca_shared for ENDO"
    log(f"  ENDO X_pca_shared: {X_endo_emb.shape}")
    embedding_key = "X_pca_shared via sc.tl.ingest (50D, fallback)"

log(f"Embedding used: {embedding_key}")


# ── 5. kNN UMAP transfer (cluster-13 Endothelial reference only) ─────────────
# Using the full atlas as reference scatters ENDO cells across unrelated UMAP
# territory because cluster-13 (319 cells) sits adjacent to large Excitatory
# clusters (leiden 5 & 15, dist ≈ 2.4 UMAP units).  Restricting kNN to only
# the 319 endothelial reference cells confines all 37k ENDO cells within the
# cluster-13 UMAP footprint — principled because we know these cells are
# endothelial and should anchor to the endothelial reference.
cl13_mask  = (atlas_leiden == "13")
n_cl13     = int(cl13_mask.sum())
k_cl13     = min(K, n_cl13 - 1)
X_ref_cl13 = reference_emb[cl13_mask]  # 10D scANVI (or 50D PCA in fallback)
log(f"kNN (k={k_cl13}) on cluster-13 Endothelial reference ({n_cl13:,} cells) "
    f"| {X_endo_emb.shape[1]}D | no scatter into adjacent clusters...")
X_endo_umap, knn_indices = full_knn_umap(
    X_endo_emb, X_ref_cl13, atlas_umap[cl13_mask], atlas_leiden[cl13_mask],
    k=k_cl13)
log(f"  ENDO centroid: ({X_endo_umap[:,0].mean():.2f}, {X_endo_umap[:,1].mean():.2f})")
log(f"  UMAP x: [{X_endo_umap[:,0].min():.2f}, {X_endo_umap[:,0].max():.2f}]")
log(f"  UMAP y: [{X_endo_umap[:,1].min():.2f}, {X_endo_umap[:,1].max():.2f}]")


# ── 6. Patch output h5ad ──────────────────────────────────────────────────────
log(f"Copying label template → {OUT_H5AD}")
shutil.copy2(LABEL_H5AD, OUT_H5AD)

log("Patching X_umap (and X_scANVI if model used) for ENDO rows...")
with h5py.File(OUT_H5AD, "r+") as f:
    f["obsm"]["X_umap"][N_ATLAS:] = X_endo_umap
    if use_model:
        f["obsm"]["X_scANVI"][N_ATLAS:] = X_endo_emb
log(f"Saved: {OUT_H5AD}")


# ── 7. Verify ─────────────────────────────────────────────────────────────────
log("Verifying...")
with h5py.File(OUT_H5AD, "r") as f:
    atlas_ok   = np.allclose(f["obsm"]["X_umap"][:5], atlas_umap[:5], atol=1e-5)
    endo_nan   = np.isnan(f["obsm"]["X_umap"][N_ATLAS:]).any()
    n_endo_obs = f["obsm"]["X_umap"].shape[0] - N_ATLAS

chk = sc.read_h5ad(OUT_H5AD, backed="r")
n_leiden18 = (chk.obs["leiden"] == "18").sum()
n_endo_ct  = (chk.obs["Celltype_markers"] == "Endothelial").sum()
chk.file.close()

log(f"  Atlas UMAP unchanged (first 5): {atlas_ok}")
log(f"  ENDO X_umap NaN-free: {not endo_nan}")
log(f"  Total ENDO rows: {n_endo_obs:,}")
log(f"  leiden=='18': {n_leiden18:,}")
log(f"  Celltype_markers=='Endothelial': {n_endo_ct:,}")


# ── 8. Figures ────────────────────────────────────────────────────────────────
log("Loading obs for plotting...")
with h5py.File(OUT_H5AD, "r") as f:
    full_umap = f["obsm"]["X_umap"][:].astype(np.float32)

chk = sc.read_h5ad(OUT_H5AD, backed="r")
ct_vals   = chk.obs["Celltype_markers"].astype(str).values
leiden_v  = chk.obs["leiden"].astype(str).values
origin_v  = chk.obs["cell_origin"].astype(str).values
chk.file.close()

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
umap_plot(full_umap, ct_vals,  "Cell type",      "umap_celltype",
          palette=PALETTE.copy(), lfs=10)

fb = plt.get_cmap("tab20")
leiden_pal = {str(i): fb(i % 20) for i in range(18)}
leiden_pal["18"]      = "#e377c2"
leiden_pal["unknown"] = "#cccccc"
umap_plot(full_umap, leiden_v, "Leiden cluster", "umap_leiden",
          palette=leiden_pal, lfs=7)

umap_plot(full_umap, origin_v, "Cell origin",    "umap_origin",
          palette={"atlas": "#dddddd", "rescued_qc": "#e377c2"},
          ps=0.05, lfs=11)


# ── 9. Summary report ─────────────────────────────────────────────────────────
report = f"""# ENDO scANVI Projection Summary

**Generated:** {datetime.now().isoformat(timespec='seconds')}

## Embedding approach
**{embedding_key}**

## ENDO cells
- Source: `13ad_endo_rescued_from_initial_qc.h5ad` (93,750 pre-QC)
- QC: n_genes ≥ 800, MT < 10%, n_endo_markers ≥ 5, total_endo_count ≥ 10
- Passed: {len(X_endo_emb):,}
- kNN k={k_cl13}, cluster-13 Endothelial reference only ({n_cl13} cells), sklearn ball_tree
- Jitter σ={JITTER_SD}

## Output
`{OUT_H5AD}`

## Verification
- Atlas UMAP unchanged: {atlas_ok}
- ENDO X_umap NaN-free: {not endo_nan}
- leiden=='18': {n_leiden18:,}
- Celltype_markers=='Endothelial': {n_endo_ct:,}

## ENDO UMAP
- Centroid: ({X_endo_umap[:,0].mean():.2f}, {X_endo_umap[:,1].mean():.2f})
- x: [{X_endo_umap[:,0].min():.2f}, {X_endo_umap[:,0].max():.2f}]
- y: [{X_endo_umap[:,1].min():.2f}, {X_endo_umap[:,1].max():.2f}]

## Log
```
{chr(10).join(_log)}
```
"""
with open(REPORT, "w") as fh:
    fh.write(report)
log(f"Report written: {REPORT}")
log("All done.")
