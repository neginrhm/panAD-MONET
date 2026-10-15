#!/usr/bin/env python3
"""Sample metadata harmonization and major cell-type annotation (Fig. 1).

Continues from 00_scvi_scanvi_allen_reference_mapping.py, in the order of the
original run:

1. Sample metadata: ADDS sample IDs get an "ADDS_" prefix, then the per-cohort
   metadata CSVs (one per study, column `Samples`) are concatenated and
   left-joined onto the cells by sample. A missing age (sample R6231318) is
   set to 90, and the "FeMaleale" typo in Sex is fixed.
2. Diagnosis2 = Diagnosis with NCI -> Control; cells whose diagnosis is
   "Not Applicable" are dropped (3,556,850 -> 3,316,964 cells).
3. Two independent major-class labels:
     Primary_Function  scANVI Allen subclass prediction mapped to its major class
                       ("unknown" -> Unclassified)
     Celltype_markers  Leiden cluster (resolution 0.3) annotated by canonical
                       markers (dot plot of MARKERS below)
   Unclassified cells take their marker-based label. Only cells where the two
   labels agree are kept (3,316,964 -> 3,158,097 cells).
   `label_transfer_annot` keeps the pre-filter label-transfer annotation.
4. One h5ad per major class.

The 33 substates and the final pre-ENDO atlas (3,146,108 cells) come from
subclustering these per-class objects; that step was interactive and is not
scripted here. The ENDO rescue (scripts 01-07) follows.

Usage:
    python 00b_metadata_and_celltype_annotation.py [--input H5AD] [--metadata-dir DIR]

Inputs:
    13datasets_after_label_transfer_leiden_doublet0.1.h5ad   output of script 00
    metadata/*.csv                                           per-cohort sample metadata;
                                                             columns incl. Samples, Diagnosis,
                                                             Sex, Age

Outputs (working directory):
    13ad_Dataset_leiden_scr0.1_wMeta.h5ad        all cells + metadata
    13AD_matched_cellmarker_labelt.h5ad          concordant cells (3,158,097)
    13ad_{Celltype}.h5ad                         one per major class
    figures/dotplot_markers_by_leiden.pdf, figures/umap_major_celltypes.pdf
"""

import argparse
import glob
import os

import pandas as pd
import scanpy as sc

# Allen subclass (scANVI prediction) -> major class
SUBCLASS_TO_MAJOR = {
    "IT": "Excitatory", "L4 IT": "Excitatory", "L5/6 IT Car3": "Excitatory",
    "L5/6 NP": "Excitatory", "L6 CT": "Excitatory",
    "LAMP5": "Inhibitory", "PVALB": "Inhibitory", "SST": "Inhibitory", "VIP": "Inhibitory",
    "Astrocyte": "Astrocyte", "Microglia": "Microglia", "OPC": "OPC",
    "Oligodendrocyte": "Oligodendrocyte",
    "unknown": "Unclassified",
}

# Leiden cluster (resolution 0.3) -> major class, from the marker dot plot
LEIDEN_TO_MAJOR = {
    "0": "Excitatory", "1": "Oligodendrocyte", "2": "Astrocyte", "3": "Inhibitory",
    "4": "Excitatory", "5": "Excitatory", "6": "Inhibitory", "7": "Microglia",
    "8": "Inhibitory", "9": "Excitatory", "10": "OPC", "11": "Oligodendrocyte",
    "12": "Inhibitory", "13": "ENDO", "14": "Inhibitory", "15": "Excitatory",
    "16": "Excitatory", "17": "Inhibitory",
}

MARKERS = {
    "Oligodendrocyte": ["MOBP", "MOG", "MBP", "PLP1"],
    "Astrocyte": ["GFAP", "SLC1A2", "AQP4", "S100B"],
    "Microglia": ["P2RY12", "CX3CR1", "CD4", "CD68", "PTPRC"],
    "Excitatory": ["RORB", "THEMIS", "SLC17A6", "SLC1A2", "TLE4", "CSRP2", "CUX1", "CUX2"],
    "Inhibitory": ["GAD2", "SST", "PVALB", "CALB2", "GAD1", "VIP", "ID2"],
    "OPC": ["PDGFRA", "CSPG4", "MBP"],
    "Endothelial cells": ["PDGFRB", "PECAM1", "VWF", "ITGA1", "FLT1", "CLDN5", "ABCB1", "ATP10A"],
}

CELL_TYPES = ["Astrocyte", "ENDO", "Excitatory", "Inhibitory", "Microglia", "OPC", "Oligodendrocyte"]


def add_sample_metadata(adata, metadata_dir):
    obs = adata.obs
    is_adds = obs["Study"] == "ADDS"
    obs["Sample"] = obs["Sample"].astype(str)
    obs.loc[is_adds, "Sample"] = "ADDS_" + obs.loc[is_adds, "Sample"]
    obs["Samples"] = obs["Sample"].astype(str)
    obs["barcode_index"] = obs.index

    files = sorted(glob.glob(os.path.join(metadata_dir, "*.csv")))
    if not files:
        raise SystemExit(f"no metadata CSVs in {metadata_dir}")
    meta = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    meta.loc[meta["Samples"] == "R6231318", "Age"] = "90"
    meta = meta.drop_duplicates(subset="Samples")
    meta["Samples"] = meta["Samples"].astype(str)

    # keep only metadata columns not already present on the cells
    meta = meta[["Samples"] + [c for c in meta.columns if c not in obs.columns]]
    merged = obs.merge(meta, on="Samples", how="left")
    merged.index = merged["barcode_index"].values
    merged.index.name = None
    merged["Sex"] = merged["Sex"].replace("FeMaleale", "Female")
    merged["Age"] = merged["Age"].astype(str)
    adata.obs = merged
    print(f"metadata from {len(files)} files; samples without metadata: "
          f"{adata.obs.loc[adata.obs['Diagnosis'].isna(), 'Samples'].nunique()}")
    return adata


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input", default="13datasets_after_label_transfer_leiden_doublet0.1.h5ad")
    p.add_argument("--metadata-dir", default="metadata")
    args = p.parse_args()
    sc.settings.figdir = "figures"

    adata = sc.read_h5ad(args.input)
    print(f"input: {adata.n_obs:,} cells")  # original: 3,556,850

    # 1-2. metadata and diagnosis filter
    adata = add_sample_metadata(adata, args.metadata_dir)
    adata.write_h5ad("13ad_Dataset_leiden_scr0.1_wMeta.h5ad")
    adata.obs["Diagnosis2"] = adata.obs["Diagnosis"].replace("NCI", "Control")
    adata = adata[adata.obs["Diagnosis2"] != "Not Applicable"].copy()
    print(f"after removing Diagnosis 'Not Applicable': {adata.n_obs:,} cells")  # original: 3,316,964

    # 3. two-way annotation
    adata.obs["Primary_Function"] = adata.obs["predictions_scanvi"].astype(str).map(SUBCLASS_TO_MAJOR)
    adata.obs["Celltype_markers"] = adata.obs["leiden"].astype(str).map(LEIDEN_TO_MAJOR)
    sc.tl.dendrogram(adata, "leiden", use_rep="X_scANVI")
    sc.pl.dotplot(adata, MARKERS, "leiden", dendrogram=True, standard_scale="var",
                  swap_axes=True, cmap="RdPu", save="markers_by_leiden.pdf")
    adata.obs["label_transfer_annot"] = adata.obs["Primary_Function"]

    pf = adata.obs["Primary_Function"].astype(object)
    unclassified = pf == "Unclassified"
    pf[unclassified] = adata.obs.loc[unclassified, "Celltype_markers"].astype(object)
    adata.obs["Primary_Function"] = pd.Categorical(pf)
    print(f"Unclassified -> marker label: {unclassified.sum():,} cells")

    agree = adata.obs["Primary_Function"].astype(str) == adata.obs["Celltype_markers"].astype(str)
    print(f"label transfer vs markers: {agree.sum():,} agree, {(~agree).sum():,} disagree")
    adata = adata[agree].copy()
    print(f"kept concordant cells: {adata.n_obs:,}")  # original: 3,158,097
    print(adata.obs["Primary_Function"].value_counts().to_string())
    sc.pl.umap(adata, color="Primary_Function", frameon=False, save="_major_celltypes.pdf")
    adata.write_h5ad("13AD_matched_cellmarker_labelt.h5ad")

    # 4. one object per major class
    for ct in CELL_TYPES:
        sub = adata[adata.obs["Primary_Function"] == ct].copy()
        sub.obs_names_make_unique()
        sub.write_h5ad(f"13ad_{ct}.h5ad")
        print(f"13ad_{ct}.h5ad: {sub.n_obs:,} cells")


if __name__ == "__main__":
    main()
