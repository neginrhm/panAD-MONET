#!/usr/bin/env python3
"""
finalize_remove_cl17_rename_cl18.py

1. Removes 8 leiden==17 cells (likely artifact).
2. Renames leiden==18 (rescued ENDO) → leiden==17 to keep numbering sequential.
3. Writes new h5ad — original is never touched.

Input:  h5ad/13ad_3m_5k_atlas_with_endo_umap_recomputed.h5ad  (3,183,466 cells)
Output: h5ad/13ad_3m_5k_atlas_final.h5ad                      (3,183,458 cells)
"""

import scanpy as sc
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE   = Path.cwd()  # Set to project root (containing h5ad/, reports/) before running
IN_H5  = BASE / "h5ad" / "13ad_3m_5k_atlas_with_endo_umap_recomputed.h5ad"
OUT_H5 = BASE / "h5ad" / "13ad_3m_5k_atlas_final.h5ad"

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

log(f"Reading {IN_H5.name}...")
adata = sc.read_h5ad(IN_H5, backed="r")
log(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

n_cl17 = (adata.obs["leiden"] == "17").sum()
n_cl18 = (adata.obs["leiden"] == "18").sum()
log(f"  leiden==17: {n_cl17} cells (removing)")
log(f"  leiden==18: {n_cl18} cells (will become leiden==17)")

# Filter and load into memory
adata = adata[adata.obs["leiden"] != "17"].to_memory()
log(f"  {adata.n_obs:,} cells after removing cluster 17")

# Rename 18 → 17
adata.obs["leiden"] = (
    adata.obs["leiden"]
    .astype(str)
    .replace({"18": "17"})
    .astype("category")
)

counts = adata.obs["leiden"].value_counts()
log("Leiden cluster counts after rename:")
for cl in sorted(counts.index, key=lambda x: int(x) if x.isdigit() else 999):
    log(f"  leiden {cl}: {counts[cl]:,}")

log(f"Writing {OUT_H5.name}...")
adata.write_h5ad(OUT_H5, compression="gzip")
log(f"Done. {OUT_H5}")
