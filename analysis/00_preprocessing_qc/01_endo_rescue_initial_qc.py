"""Marker-based ENDO/VASC rescue from the initial-QC h5ad.

The 276 GB initial-QC file holds 5.39M cells (vs 3.15M in the final dataset),
so 2.24M cells were dropped between initial QC and final. This script:

1. Streams the initial-QC h5ad in chunks (backed mode)
2. For each cell, counts how many endothelial markers have raw count > 0
   and the total counts across markers
3. Flags candidate ENDO cells using strict criteria (markers ≥ 4, total ≥ 5)
4. Cross-references against the final 3.1M dataset to find cells that:
     - are flagged as endo by markers, AND
     - are NOT in the final dataset (dropped during downstream QC)
5. Writes a sibling h5ad of just those rescued cells with raw counts.

Existing files are untouched.
"""

import os
import json
from datetime import datetime

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

BASE = os.getcwd()  # Set to project root (containing h5ad/, reports/) before running
INIT_QC = f"{BASE}/h5ad/13data_qcmetrics_initial_qc.h5ad"
FINAL = f"{BASE}/h5ad/13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad"
OUT_H5 = f"{BASE}/h5ad/13ad_endo_rescued_from_initial_qc.h5ad"
OUT_REPORT = f"{BASE}/reports/endo_rescue_initial_qc_summary.md"

# Markers — same as the prior rescue
ENDO_MARKERS = ["CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "ERG", "KDR", "ESAM"]
VASC_MARKERS_EXTRA = ["PDGFRB", "RGS5", "MYH11", "ACTA2", "DCN", "COL1A2"]  # mural / VLMC
ALL_MARKERS = ENDO_MARKERS + VASC_MARKERS_EXTRA

# Strict thresholds for un-annotated rescue
MIN_ENDO_MARKERS = 4         # ≥4 of 8 endothelial markers expressed (raw count > 0)
MIN_TOTAL_ENDO_COUNT = 5     # sum of raw counts across endo markers ≥ 5
MIN_GENES_RELAXED = 200      # very lenient cell-level QC for ENDO (vs typical 500-1000)

CHUNK = 250_000

stats = {"started": datetime.now().isoformat()}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 1. Build set of (Sample, org_barcode) keys already in the final dataset
# ---------------------------------------------------------------------------
log("Loading final dataset obs (backed)...")
fa = ad.read_h5ad(FINAL, backed="r")
fa_keys = (
    fa.obs["Sample"].astype(str) + "||" + fa.obs["org_barcode"].astype(str)
).values
final_key_set = set(fa_keys)
log(f"final keys: {len(final_key_set):,}")
del fa, fa_keys

# ---------------------------------------------------------------------------
# 2. Open initial-QC h5ad in backed mode and find endo marker var indices
# ---------------------------------------------------------------------------
log("Opening initial-QC h5ad (backed)...")
a = ad.read_h5ad(INIT_QC, backed="r")
log(f"initial-QC shape: {a.shape}")

present_endo = [g for g in ENDO_MARKERS if g in a.var_names]
present_extra = [g for g in VASC_MARKERS_EXTRA if g in a.var_names]
log(f"Endo markers present: {present_endo}")
log(f"Mural/VLMC markers present: {present_extra}")
endo_idx = np.array([a.var_names.get_loc(g) for g in present_endo])
extra_idx = np.array([a.var_names.get_loc(g) for g in present_extra])
all_idx = np.concatenate([endo_idx, extra_idx])

# ---------------------------------------------------------------------------
# 3. Stream in chunks; compute per-cell marker stats
# ---------------------------------------------------------------------------
n = a.n_obs
n_endo_expressed = np.zeros(n, dtype=np.int8)
total_endo_count = np.zeros(n, dtype=np.float32)
n_extra_expressed = np.zeros(n, dtype=np.int8)
total_extra_count = np.zeros(n, dtype=np.float32)

log("Streaming chunks for marker stats...")
processed = 0
for start in range(0, n, CHUNK):
    end = min(start + CHUNK, n)
    sub = a.X[start:end, :]  # backed slicing
    if sparse.issparse(sub):
        sub_endo = sub[:, endo_idx]
        sub_extra = sub[:, extra_idx]
        n_endo_expressed[start:end] = np.asarray((sub_endo > 0).sum(axis=1)).ravel()
        total_endo_count[start:end] = np.asarray(sub_endo.sum(axis=1)).ravel()
        n_extra_expressed[start:end] = np.asarray((sub_extra > 0).sum(axis=1)).ravel()
        total_extra_count[start:end] = np.asarray(sub_extra.sum(axis=1)).ravel()
    else:
        sub_endo = sub[:, endo_idx]
        sub_extra = sub[:, extra_idx]
        n_endo_expressed[start:end] = (sub_endo > 0).sum(axis=1)
        total_endo_count[start:end] = sub_endo.sum(axis=1)
        n_extra_expressed[start:end] = (sub_extra > 0).sum(axis=1)
        total_extra_count[start:end] = sub_extra.sum(axis=1)
    processed = end
    if processed % (CHUNK * 4) == 0 or processed == n:
        log(f"  processed {processed:,} / {n:,} ({100*processed/n:.1f}%)")

# ---------------------------------------------------------------------------
# 4. Build candidate flag using strict criteria + cell-level lenient QC
# ---------------------------------------------------------------------------
log("Building candidate flags...")
obs = a.obs.copy()
obs["n_endo_markers"] = n_endo_expressed
obs["total_endo_count"] = total_endo_count
obs["n_mural_markers"] = n_extra_expressed
obs["total_mural_count"] = total_extra_count

# Endo candidate: strict marker criterion + lenient cell-level QC
endo_candidate = (
    (n_endo_expressed >= MIN_ENDO_MARKERS)
    & (total_endo_count >= MIN_TOTAL_ENDO_COUNT)
    & (obs["n_genes_by_counts"].values >= MIN_GENES_RELAXED)
)

# Mural / vascular wall candidate: ≥3 of 6 mural markers AND ≥1 endo marker
mural_candidate = (
    (n_extra_expressed >= 3)
    & (n_endo_expressed >= 1)
    & (obs["n_genes_by_counts"].values >= MIN_GENES_RELAXED)
)

obs["endo_candidate"] = endo_candidate
obs["mural_candidate"] = mural_candidate
obs["any_vasc_candidate"] = endo_candidate | mural_candidate

stats["initial_qc_total_cells"] = int(n)
stats["endo_candidates"] = int(endo_candidate.sum())
stats["mural_candidates"] = int(mural_candidate.sum())
stats["any_vasc_candidates"] = int(obs["any_vasc_candidate"].sum())
log(f"Endo candidates (strict markers + n_genes>=200): {stats['endo_candidates']:,}")
log(f"Mural candidates: {stats['mural_candidates']:,}")
log(f"Any vascular candidates: {stats['any_vasc_candidates']:,}")

# ---------------------------------------------------------------------------
# 5. Identify which candidates are NOT in the final 3.1M dataset
# ---------------------------------------------------------------------------
log("Cross-referencing against final dataset keys...")
init_keys = (obs["Sample"].astype(str) + "||" + obs["org_barcode"].astype(str)).values
in_final = np.array([k in final_key_set for k in init_keys])
obs["in_final_dataset"] = in_final

new_endo = endo_candidate & (~in_final)
new_mural = mural_candidate & (~in_final)
new_any = obs["any_vasc_candidate"].values & (~in_final)

stats["candidates_already_in_final"] = int(obs["any_vasc_candidate"].values.sum() - new_any.sum())
stats["new_endo_rescued"] = int(new_endo.sum())
stats["new_mural_rescued"] = int(new_mural.sum())
stats["new_any_vasc_rescued"] = int(new_any.sum())
log(f"NEW endo rescued (not in final): {stats['new_endo_rescued']:,}")
log(f"NEW mural rescued: {stats['new_mural_rescued']:,}")
log(f"Per-study breakdown of NEW any-vasc rescue:")
per_study = obs.loc[new_any, "Study"].value_counts().to_dict()
stats["per_study_new_rescue"] = {str(k): int(v) for k, v in per_study.items()}
for k, v in per_study.items():
    log(f"  {k}: {v:,}")

# ---------------------------------------------------------------------------
# 6. Build & save the rescued AnnData (raw counts + obs metadata)
# ---------------------------------------------------------------------------
log(f"Building rescued AnnData with {new_any.sum():,} cells...")
rescue_idx = np.where(new_any)[0]
# Slice X chunk-wise to save memory
log("Slicing X for rescued cells (chunked)...")
parts = []
for start in range(0, len(rescue_idx), CHUNK):
    end = min(start + CHUNK, len(rescue_idx))
    parts.append(a.X[rescue_idx[start:end], :])
    if (start // CHUNK) % 4 == 0:
        log(f"  sliced {end:,} / {len(rescue_idx):,}")
X_new = sparse.vstack(parts).tocsr() if sparse.issparse(parts[0]) else np.vstack(parts)

# Build obs subset
obs_new = obs.iloc[rescue_idx].copy()
obs_new["rescue_class"] = np.where(
    obs_new["endo_candidate"].values, "ENDO_marker_rescued",
    np.where(obs_new["mural_candidate"].values, "MURAL_marker_rescued", "OTHER")
)
obs_new["rescue_source"] = "initial_qc_marker"

new_adata = ad.AnnData(X=X_new, obs=obs_new, var=a.var.copy())
new_adata.obs_names = obs_new.index.astype(str).values
log(f"new_adata: {new_adata.shape}")

log(f"Writing {OUT_H5}...")
new_adata.write_h5ad(OUT_H5, compression="gzip")
log("h5ad written")

# ---------------------------------------------------------------------------
# 7. Markdown report
# ---------------------------------------------------------------------------
stats["finished"] = datetime.now().isoformat()
log(f"Writing report: {OUT_REPORT}")
lines = []
lines.append("# ENDO/VASC Rescue from Initial-QC h5ad")
lines.append("")
lines.append(f"**Generated:** {stats['finished']}")
lines.append("")
lines.append("**Method:** Strict marker-based detection on 5.39M cells in the initial-QC dataset, ")
lines.append("identifying vascular cells that were dropped during downstream QC.")
lines.append("")
lines.append("## Inputs / outputs")
lines.append(f"- Input initial-QC h5ad: `{INIT_QC}`")
lines.append(f"- Final 3.1M h5ad (used as exclusion key): `{FINAL}`")
lines.append(f"- Output rescued h5ad: `{OUT_H5}`")
lines.append("")
lines.append("## Detection criteria")
lines.append(f"- Endothelial: ≥{MIN_ENDO_MARKERS}/8 endo markers expressed (raw count >0) AND total endo counts ≥ {MIN_TOTAL_ENDO_COUNT}")
lines.append(f"  - Markers used: {', '.join(present_endo)}")
lines.append(f"- Mural / VLMC: ≥3/6 mural markers AND ≥1 endo marker")
lines.append(f"  - Markers used: {', '.join(present_extra)}")
lines.append(f"- Cell QC (relaxed): n_genes_by_counts ≥ {MIN_GENES_RELAXED}")
lines.append("")
lines.append("## Results")
lines.append(f"- Initial-QC total cells: {stats['initial_qc_total_cells']:,}")
lines.append(f"- Endo candidates (any QC): **{stats['endo_candidates']:,}**")
lines.append(f"- Mural candidates: **{stats['mural_candidates']:,}**")
lines.append(f"- Any vascular candidates: **{stats['any_vasc_candidates']:,}**")
lines.append(f"- Already in final 3.1M (not rescued): {stats['candidates_already_in_final']:,}")
lines.append(f"- **NEW (dropped during QC, now rescued): {stats['new_any_vasc_rescued']:,}**")
lines.append(f"  - of those, ENDO: {stats['new_endo_rescued']:,}")
lines.append(f"  - of those, MURAL: {stats['new_mural_rescued']:,}")
lines.append("")
lines.append("## Per-study NEW rescue counts")
lines.append("| Study | NEW vascular cells rescued |")
lines.append("|---|---|")
for k, v in sorted(stats["per_study_new_rescue"].items(), key=lambda x: -x[1]):
    lines.append(f"| {k} | {v:,} |")
lines.append("")
lines.append("## How to use the rescued h5ad")
lines.append("- File contains raw counts only (no normalization, PCA, UMAP).")
lines.append("- For ENDO-specific analyses (DEGs, hdWGCNA), normalize and HVG-select on this file.")
lines.append("- To put rescued cells on the existing UMAP, project via scArches/scvi-tools `load_query_data` using the trained scVI/scANVI model.")
lines.append("- Existing files unchanged — original 3.1M atlas remains the canonical reference.")
lines.append("")
lines.append("## Stats dump")
lines.append("```json")
lines.append(json.dumps(stats, indent=2, default=str))
lines.append("```")

os.makedirs(os.path.dirname(OUT_REPORT), exist_ok=True)
with open(OUT_REPORT, "w") as f:
    f.write("\n".join(lines))
log("Report written")
log("DONE")
