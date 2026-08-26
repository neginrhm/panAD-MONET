#!/usr/bin/env python3
"""
Build per-gene, per-sample 7-element cell-type expression vectors.

Cleaned up from Merging_regressed_pseudobulked_values_vectors.ipynb (Cell 47/49).

Input:
    Regressed_data_DEM_{CT}_thisexp~Study+Age+Sex.csv  (7 files, one per cell type)
    produced by 01_regress_covariates.r

Output:
    input_dnn_celltypes_pseudobulked_mapped_VectorsofFloats.csv
    Shape: genes x samples
    Each cell = list [AST, MIC, ENDO, ODC, OPC, INH, EXC] float values
    NaN in a position = that cell type had no data for that gene/sample
"""

import os
import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path

DATA_DIR = Path(".")
OUT_FILE = Path("input_dnn_celltypes_pseudobulked_mapped_VectorsofFloats.csv")

# Fixed cell type order — must match VAE preprocessing in 03_vae_ensemble_pipeline.py
CELLTYPE_MAP = {
    "AST":  "A",
    "MIC":  "M",
    "ENDO": "E",
    "ODC":  "D",
    "OPC":  "P",
    "INH":  "I",
    "EXC":  "X",
}
ORDERED_INITIALS = list(CELLTYPE_MAP.values())  # [A, M, E, D, P, I, X]


def extract_celltype_initial(filename):
    for ct, initial in CELLTYPE_MAP.items():
        if f"_{ct}_" in filename:
            return initial
    return None


# --- load all 7 regressed CSVs ---
csv_files = sorted(glob(str(DATA_DIR / "Regressed_data_DEM_*_thisexp~Study+Age+Sex.csv")))
if len(csv_files) != 7:
    print(f"WARNING: expected 7 CSVs, found {len(csv_files)}")

gene_sample_matrix = {}  # gene -> sample -> {initial: float}

for fpath in csv_files:
    ct_initial = extract_celltype_initial(os.path.basename(fpath))
    if ct_initial is None:
        print(f"Skipping unrecognized file: {fpath}")
        continue
    print(f"Loading {os.path.basename(fpath)} -> position [{ct_initial}]")
    df = pd.read_csv(fpath, index_col=0)
    for gene in df.index:
        for sample in df.columns:
            value = df.at[gene, sample]
            if pd.notna(value):
                gene_sample_matrix \
                    .setdefault(gene, {}) \
                    .setdefault(sample, {})[ct_initial] = float(value)

all_samples = sorted({s for gd in gene_sample_matrix.values() for s in gd})
print(f"\nGenes: {len(gene_sample_matrix)}  |  Samples: {len(all_samples)}")

# --- assemble output: each cell = [A, M, E, D, P, I, X] list ---
output_df = pd.DataFrame(index=gene_sample_matrix.keys(), columns=all_samples)

for gene, sample_dict in gene_sample_matrix.items():
    for sample in all_samples:
        cell_data = sample_dict.get(sample, {})
        vec = [cell_data.get(ct, np.nan) for ct in ORDERED_INITIALS]
        output_df.at[gene, sample] = vec

output_df.to_csv(OUT_FILE)
print(f"Saved: {OUT_FILE}  ({output_df.shape[0]} genes x {output_df.shape[1]} samples)")
