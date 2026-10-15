#!/usr/bin/env python3
"""pySCENIC per cell type (Fig. 5 prior GRN): loom -> grn -> ctx -> aucell.

For one major cell type, splits nuclei by diagnosis (Control / MCI /
Dementia), writes a loom per group, and runs the three pySCENIC CLI steps:
  1. pyscenic grn     GRNBoost2 co-expression TF -> target adjacencies
  2. pyscenic ctx     cisTarget motif pruning into regulons (--mask_dropouts)
  3. pyscenic aucell  per-cell regulon activity

The Dementia-group regulons are the "prior GRN" that NetREm refined for Fig. 5
(NetREm itself was run by a collaborator and is not part of this repository).

Settings follow the original runs: hg38 refseq-r80 500bp-up/100bp-down mc9nr
ranking database, motifs-v9-nr.hgnc annotations, human TF list
allTFs_hg38.txt. GRNBoost2 is stochastic and was not seeded, so a rerun
reproduces the regulons approximately, not edge-for-edge.

Usage:
    python 00a_pyscenic_per_celltype.py <CellType> [--groups dementia control mci]
                                        [--steps grn ctx aucell] [--workers 50]
    e.g. python 00a_pyscenic_per_celltype.py Oligodendrocyte

Inputs (working directory):
    13ad_{CellType}.h5ad      raw counts in .X; obs 'Diagnosis'
    allTFs_hg38.txt
    hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    motifs-v9-nr.hgnc-m0.001-o0.0.tbl

Outputs (working directory), per group g in {dementia, control, mci}:
    {CellType}_scenic_adata_{g}.loom
    {CellType}_adj_adata_{g}.csv           (grn)
    {CellType}_{dem|ctrl|mci}_reg.csv      (ctx)
    {CellType}_auc_matrix.csv              (aucell, dementia; original file name)
    {CellType}_{ctrl|mci}_auc_matrix.csv   (aucell, control / MCI)
"""

import argparse
import os
import subprocess

import loompy as lp
import numpy as np
import scanpy as sc

TF_FILE = "allTFs_hg38.txt"
RANKING_DB = "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_TBL = "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

DEMENTIA = ["ADDS", "AD", "AD, Vascular Dementia", "Dementia", "DSAD", "PD", "Vascular Dementia"]
CONTROL = ["No Dementia", "NCI", "Control"]

# group -> (dx label, short tag used in the original ctx/aucell file names)
GROUPS = {"dementia": ("Dementia", "dem"), "control": ("Control", "ctrl"), "mci": ("MCI", "mci")}


def assign_dx(adata):
    d = adata.obs["Diagnosis"]
    adata.obs["dx"] = np.where(d.isin(DEMENTIA), "Dementia",
                               np.where(d.isin(CONTROL), "Control",
                                        np.where(d == "MCI", "MCI", None)))
    return adata


def create_loom_file(adata, loom_path):
    row_attrs = {"Gene": np.array(adata.var_names)}
    col_attrs = {
        "CellID": np.array(adata.obs_names),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
    }
    lp.create(loom_path, adata.X.transpose(), row_attrs, col_attrs)


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("cell_type", help="e.g. Oligodendrocyte, Microglia, Astrocyte, Excitatory, Inhibitory")
    p.add_argument("--groups", nargs="+", default=["dementia", "control", "mci"],
                   choices=["dementia", "control", "mci"])
    p.add_argument("--steps", nargs="+", default=["grn", "ctx", "aucell"],
                   choices=["grn", "ctx", "aucell"])
    p.add_argument("--workers", type=int, default=50)
    args = p.parse_args()
    ct = args.cell_type

    adata = sc.read_h5ad(f"13ad_{ct}.h5ad")
    adata = adata[adata.obs["Diagnosis"] != "Brain Cancer"].copy()
    adata = assign_dx(adata)
    print(adata.obs["dx"].value_counts(dropna=False))

    for g in args.groups:
        loom = f"{ct}_scenic_adata_{g}.loom"
        adj = f"{ct}_adj_adata_{g}.csv"
        dx_label, short = GROUPS[g]
        reg = f"{ct}_{short}_reg.csv"
        auc = f"{ct}_auc_matrix.csv" if g == "dementia" else f"{ct}_{short}_auc_matrix.csv"

        if "grn" in args.steps:
            sub = adata[adata.obs["dx"] == dx_label].copy()
            print(f"{ct} {g}: {sub.n_obs:,} cells")
            create_loom_file(sub, loom)
            run(["pyscenic", "grn", loom, TF_FILE, "-o", adj, "--num_workers", str(args.workers)])
        if "ctx" in args.steps:
            run(["pyscenic", "ctx", adj, RANKING_DB, "--annotations_fname", MOTIF_TBL,
                 "--expression_mtx_fname", loom, "--output", reg, "--mask_dropouts",
                 "--num_workers", str(args.workers)])
        if "aucell" in args.steps:
            run(["pyscenic", "aucell", loom, reg, "--output", auc, "--num_workers", str(args.workers)])


if __name__ == "__main__":
    main()
