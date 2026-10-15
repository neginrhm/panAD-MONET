#!/usr/bin/env python3
"""Build the up/down gene lists submitted to CLUE (Connectivity Map) per subtype.

For each subtype DEG table (`DEG_*_vs_REST.csv`, scanpy rank_genes_groups
output with columns names / logfoldchanges / pvals_adj), keeps genes with
adjusted p < 0.05 and writes the 150 most up- and 150 most down-regulated
genes (by log fold change) as plain-text lists. These lists were uploaded to
the CLUE L1000 query tool (clue.io); CLUE keeps only genes in its BING space,
so each submitted query contained ~70-110 genes per direction. The returned
.gct results are parsed by 04b_parse_clue_gct.R.

The input tables are scanpy Wilcoxon rank-sum DEGs of each subtype versus
all other samples (rank_genes_groups, one table per subtype).

Usage:
    python 04a_make_clue_gene_lists.py [--input-dir DIR] [--top-n 150] [--padj 0.05]

Inputs:  DEG_*_vs_REST.csv in --input-dir (default: current directory)
Outputs: DEG_*_vs_REST_UP_top150.txt, DEG_*_vs_REST_DOWN_top150.txt (same directory)
"""

import argparse
import glob
import os

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input-dir", default=".")
    p.add_argument("--top-n", type=int, default=150)
    p.add_argument("--padj", type=float, default=0.05)
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "DEG_*_vs_REST.csv")))
    if not files:
        raise SystemExit(f"no DEG_*_vs_REST.csv files in {args.input_dir}")

    for file in files:
        df = pd.read_csv(file)
        required = {"names", "logfoldchanges", "pvals_adj"}
        if not required.issubset(df.columns):
            raise ValueError(f"missing one of {required} in {file}")

        sig = df[df["pvals_adj"] < args.padj]
        up = sig.sort_values("logfoldchanges", ascending=False).head(args.top_n)["names"]
        down = sig.sort_values("logfoldchanges", ascending=True).head(args.top_n)["names"]

        base = os.path.splitext(file)[0]
        up.to_csv(f"{base}_UP_top{args.top_n}.txt", index=False, header=False)
        down.to_csv(f"{base}_DOWN_top{args.top_n}.txt", index=False, header=False)
        print(f"{os.path.basename(file)}: {len(sig)} significant -> {len(up)} up, {len(down)} down")


if __name__ == "__main__":
    main()
