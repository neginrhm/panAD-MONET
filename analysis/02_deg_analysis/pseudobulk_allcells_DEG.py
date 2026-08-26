"""
pseudobulk_allcells_DEG.py
==========================
Pseudobulk DEG analysis aggregating all cells per sample.

1) Creates pseudobulk count matrix by summing counts per sample (all cell types).
2) Runs DESeq2 via rpy2 (if available) for pairwise and each-vs-REST comparisons
   across {MCI, Control, 1, 2, 3, 4}.

Inputs:
  - H5AD: AnnData object with raw counts and sample/condition annotations
  - META: CSV with per-sample metadata (optional; set to None to skip)

Outputs:
  - OUTDIR: directory containing:
    - counts_allcells.tsv: pseudobulk count matrix (genes x samples)
    - coldata_allcells.tsv: sample metadata
    - DESeq2_results/: pairwise and vs-REST DE result CSVs (if rpy2 available)

Usage:
  python pseudobulk_allcells_DEG.py
"""

import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# --- Configuration ---

BASE_DIR = "./"
H5AD = os.path.join(BASE_DIR, "data/h5ad/12AD_matched_label_final.h5ad")
META = "20_seeds_final_subgroups_with_meta.csv"   # set to None to skip
SAMPLE_KEY = "Sample"
COND_KEY = "Cluster_or_dx"
OUTDIR = "pseudobulk_allcells_DEG"

CANDIDATES = {"MCI", "Control", "1", "2", "3", "4"}

# --- Helper functions ---

def normalize_label(x):
    """Normalize label strings (e.g., '1.0' -> '1')."""
    s = str(x).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


os.makedirs(OUTDIR, exist_ok=True)

# --- Load data ---

adata = sc.read_h5ad(H5AD)

if META:
    df = pd.read_csv(META)
    df_ = df.drop_duplicates("Samples").set_index("Samples")
    adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype(str)
    df_.index = df_.index.astype(str)
    adata.obs = adata.obs.join(df_, on=SAMPLE_KEY, how="left", rsuffix="_meta")

assert SAMPLE_KEY in adata.obs, f"{SAMPLE_KEY} not in adata.obs"
assert COND_KEY in adata.obs, f"{COND_KEY} not in adata.obs"

# --- Select counts matrix (prefer raw) ---

if "counts" in adata.layers:
    X = adata.layers["counts"]
elif adata.raw is not None:
    X = adata.raw.X
else:
    X = adata.X
    print("WARNING: using adata.X (likely normalized). Prefer raw counts for pseudobulk.")

samples = adata.obs[SAMPLE_KEY].astype(str).values
conds = adata.obs[COND_KEY].astype(str).values

# Map sample -> one condition (assumes fixed per sample)
sample_to_cond = pd.Series(conds, index=samples).groupby(level=0).first().to_dict()

# --- Build pseudobulk: one row per sample (all cells) ---

labels, inv = np.unique(samples, return_inverse=True)
n_cells, n_genes = adata.n_obs, adata.n_vars

if sp.issparse(X):
    X = X.tocsr()
    M = sp.csr_matrix((np.ones(n_cells), (inv, np.arange(n_cells))), shape=(len(labels), n_cells))
    PB = (M @ X).A   # samples x genes
else:
    PB = np.zeros((len(labels), n_genes), dtype=float)
    for i in range(len(labels)):
        PB[i] = X[inv == i].sum(axis=0)

# Sample metadata
meta = pd.DataFrame({
    "sample": labels,
    "condition": [normalize_label(sample_to_cond[s]) for s in labels],
})

# Keep only target groups actually present
meta = meta[meta["condition"].isin(CANDIDATES)].copy()
keep_idx = [i for i, s in enumerate(labels) if s in set(meta["sample"])]
PB = PB[keep_idx, :]
labels_kept = [labels[i] for i in keep_idx]
assert list(meta["sample"]) == labels_kept

# --- Save pseudobulk matrices ---

cts = pd.DataFrame(PB.T, index=adata.var_names, columns=meta["sample"].values)
cts_path = os.path.join(OUTDIR, "counts_allcells.tsv")
col_path = os.path.join(OUTDIR, "coldata_allcells.tsv")
cts.to_csv(cts_path, sep="\t")
meta.to_csv(col_path, sep="\t", index=False)
print(f"[pseudobulk] wrote:\n  {cts_path}\n  {col_path}")

# --- DESeq2 via rpy2 (if available) ---

def run_deseq2_with_rpy2(cts_df, col_df, out_dir):
    """Run DESeq2 pairwise and vs-REST comparisons using rpy2."""
    from rpy2 import robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    r = ro.r

    r('suppressPackageStartupMessages({library(DESeq2);})')
    r('suppressWarnings(quietly <- requireNamespace("apeglm", quietly=TRUE))')

    # Send data to R
    r.assign("cts_py", pandas2ri.py2rpy(cts_df))
    r.assign("col_py", pandas2ri.py2rpy(col_df))

    r('cts <- as.matrix(cts_py)')
    r('coldata <- as.data.frame(col_py)')
    r('rownames(coldata) <- coldata$sample')
    r('cts <- cts[, rownames(coldata), drop=FALSE]')
    r('coldata$condition <- droplevels(factor(coldata$condition))')

    # Design without intercept so each level has a coefficient
    r('dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ 0 + condition)')
    r('keep <- rowSums(counts(dds)) >= 10; dds <- dds[keep,]')
    r('dds <- DESeq(dds)')

    groups = list(col_df["condition"].unique())
    os.makedirs(out_dir, exist_ok=True)

    def write_res(r_var, filename):
        r('res_df <- as.data.frame(%s); res_df$gene <- rownames(%s)' % (r_var, r_var))
        df = pandas2ri.rpy2py(r('res_df'))
        df.to_csv(os.path.join(out_dir, filename), index=False)

    # Pairwise both directions
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            r.assign("g1", f"condition{g1}")
            r.assign("g2", f"condition{g2}")

            r('res <- results(dds, contrast = list(g1, g2))')
            r('if ("apeglm" %in% rownames(installed.packages())) { '
              '  nm <- g1; if (nm %in% resultsNames(dds)) { '
              '    suppressWarnings(res <- lfcShrink(dds, coef=nm, type="apeglm", res=res)) } }')
            write_res("res", f"DE_allcells_{g1[9:]}_vs_{g2[9:]}.csv")

            r('res2 <- results(dds, contrast = list(g2, g1))')
            r('if ("apeglm" %in% rownames(installed.packages())) { '
              '  nm <- g2; if (nm %in% resultsNames(dds)) { '
              '    suppressWarnings(res2 <- lfcShrink(dds, coef=nm, type="apeglm", res=res2)) } }')
            write_res("res2", f"DE_allcells_{g2[9:]}_vs_{g1[9:]}.csv")

    # Each group vs REST (mean of others)
    for g in groups:
        others = [x for x in groups if x != g]
        if not others:
            continue
        r.assign("contrast_list", ro.StrVector([f"condition{g}"] + [f"condition{o}" for o in others]))
        r.assign("weights", ro.FloatVector([1.0] + [-1.0 / len(others)] * len(others)))
        r('res_rest <- results(dds, contrast = list(contrast_list, weights))')
        r('if ("apeglm" %in% rownames(installed.packages())) { '
          '  nm <- paste0("condition", "%s"); if (nm %in% resultsNames(dds)) { '
          '    suppressWarnings(res_rest <- lfcShrink(dds, coef=nm, type="apeglm", res=res_rest)) } }' % g)
        write_res("res_rest", f"DE_allcells_{g}_vs_REST.csv")

    print(f"[DESeq2] results written to: {out_dir}")


# Try rpy2; if missing, print instructions for running in R
try:
    out_res = os.path.join(OUTDIR, "DESeq2_results")
    run_deseq2_with_rpy2(cts, meta, out_res)
except Exception as e:
    print(f"[info] rpy2/DESeq2 step not run ({e}).")
    print("Run in R instead:\n"
          "  library(DESeq2)\n"
          f"  cts <- as.matrix(read.table('{cts_path}', sep='\\t', header=TRUE, row.names=1, check.names=FALSE))\n"
          f"  col <- read.table('{col_path}', sep='\\t', header=TRUE, check.names=FALSE)\n"
          "  rownames(col) <- col$sample; cts <- toupper(cts*1) * 1 # noop to ensure numeric\n"
          "  col$condition <- factor(col$condition); dds <- DESeqDataSetFromMatrix(cts, col, design=~0+condition)\n"
          "  keep <- rowSums(counts(dds)) >= 10; dds <- dds[keep,]; dds <- DESeq(dds)\n"
          "  # now use results(dds, contrast=list('conditionMCI','conditionControl')) etc.\n")
