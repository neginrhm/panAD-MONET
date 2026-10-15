#!/usr/bin/env python3
"""pySCENIC GRN step for the sex-stratified analysis (Fig. 6c-d): loom -> grn.

Writes one loom per sex x diagnosis group and runs `pyscenic grn` (GRNBoost2)
to produce the TF -> target adjacency lists. Next steps:
  01_filter_adjacency_noncoding.py  -> adata_{group}_adj_scenic_filteredGenes.csv
  scenic_ctx_aucell_viz.py          -> ctx, aucell, plots

GRNBoost2 is stochastic and was not seeded, so a rerun reproduces the
networks approximately (in the original runs, two GRN runs on the same
group shared ~54-56% of edges).

Usage:
    python 00b_pyscenic_grn_sex_stratified.py [--groups dm_male dm_female ...] [--workers 40]

Inputs (working directory):
    13adata_{ct,mci,dem}_{male,female}.h5ad   raw counts in .X
    allTFs_hg38.txt

Outputs (working directory), per group (ct_male, ct_female, mci_male,
mci_female, dm_male, dm_female):
    adata_{group}_scenic.loom
    adata_{group}_adj.csv
"""

import argparse
import subprocess

import loompy as lp
import numpy as np
import scanpy as sc

TF_FILE = "allTFs_hg38.txt"
GROUP_H5AD = {
    "ct_male": "13adata_ct_male.h5ad",
    "ct_female": "13adata_ct_female.h5ad",
    "mci_male": "13adata_mci_male.h5ad",
    "mci_female": "13adata_mci_female.h5ad",
    "dm_male": "13adata_dem_male.h5ad",
    "dm_female": "13adata_dem_female.h5ad",
}


def create_loom_file(adata, loom_path):
    row_attrs = {"Gene": np.array(adata.var_names)}
    col_attrs = {
        "CellID": np.array(adata.obs_names),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
    }
    lp.create(loom_path, adata.X.transpose(), row_attrs, col_attrs)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--groups", nargs="+", default=list(GROUP_H5AD), choices=list(GROUP_H5AD))
    p.add_argument("--workers", type=int, default=40)
    args = p.parse_args()

    for g in args.groups:
        adata = sc.read_h5ad(GROUP_H5AD[g])
        loom, adj = f"adata_{g}_scenic.loom", f"adata_{g}_adj.csv"
        print(f"{g}: {adata.n_obs:,} cells -> {loom}", flush=True)
        create_loom_file(adata, loom)
        cmd = ["pyscenic", "grn", loom, TF_FILE, "-o", adj, "--num_workers", str(args.workers)]
        print("+", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        del adata


if __name__ == "__main__":
    main()
