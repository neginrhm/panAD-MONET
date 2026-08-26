#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd

# non-interactive backend for headless servers/terminals
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def pick_col(cols, candidates):
    """
    Given a list-like of column names and a list of candidate names (case-insensitive),
    return the first matching column from `cols`. Raise KeyError if none found.
    """
    lc = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lc:
            return lc[cand.lower()]
    raise KeyError(f"Could not find any of {candidates} in columns: {list(cols)}")


def load_final_labels(final_csv, id_col=None, label_col=None):
    """
    Load final labels with flexible column detection.
    """
    df = pd.read_csv(final_csv)
    # if user didn't specify, auto-detect common names
    if id_col is None:
        id_col = pick_col(
            df.columns,
            ["SampleID", "sample_id", "sample", "id", "Sample", "ID"],
        )
    if label_col is None:
        label_col = pick_col(
            df.columns,
            ["FinalCluster", "final", "cluster", "label", "final_label", "Cluster"],
        )
    ser = df.set_index(id_col)[label_col]
    ser.index = ser.index.astype(str)
    return ser


def main(pairwise_csv, final_csv, out_dir, id_col=None, label_col=None):
    os.makedirs(out_dir, exist_ok=True)

    # ---------- Load ----------
    M = pd.read_csv(pairwise_csv, index_col=0)
    # normalize index to str
    M.index = M.index.astype(str)
    M.columns = M.columns.astype(str)

    final = load_final_labels(final_csv, id_col=id_col, label_col=label_col)

    # ---------- Align ----------
    common = M.index.intersection(M.columns).intersection(final.index)
    if len(common) == 0:
        raise ValueError(
            "No overlapping sample IDs between pairwise matrix and final labels.\n"
            f"pairwise rows n={len(M.index)}, cols n={len(M.columns)}; "
            f"final n={len(final.index)}"
        )
    if len(common) < len(final.index):
        print(f"[warn] Only {len(common)}/{len(final)} samples overlap. Proceeding with intersection.")

    M = M.loc[common, common]
    final = final.loc[common]

    # ---------- Sanity checks ----------
    A = M.values
    assert (M.index == M.columns).all(), "Row/column labels must match and be in the same order."
    if not np.allclose(A, A.T, equal_nan=True):
        raise ValueError("Co-assignment matrix is not symmetric.")
    if not np.allclose(np.diag(A), 1.0):
        raise ValueError("Diagonal should be 1.0 in co-assignment matrix.")
    if not (((A >= 0) & (A <= 1)).all()):
        raise ValueError("Co-assignment entries must be within [0,1].")

    # ---------- Reorder by FINAL labels ----------
    order = final.sort_values(kind="stable").index  # stable sort preserves ID order within clusters
    M_ord = M.loc[order, order]
    final_ord = final.loc[order]
    M_ord.to_csv(os.path.join(out_dir, "pairwise_coassignment_ordered_by_FINAL.csv"))

    # ---------- 1) Heatmap ----------
    plt.figure(figsize=(7.5, 6.5), dpi=200)
    im = plt.imshow(M_ord.values, vmin=0, vmax=1, cmap="viridis", interpolation="nearest")
    cbar = plt.colorbar(im)
    cbar.set_label("Co-assignment fraction", rotation=90)
    plt.title("Pairwise co-assignment (ordered by FINAL clusters)")
    plt.xticks([], [])
    plt.yticks([], [])
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "coassignment_heatmap.png"))
    plt.savefig(os.path.join(out_dir, "coassignment_heatmap.pdf"))
    plt.close()

    # ---------- 2) Within vs Between distributions ----------
    labels = final_ord.to_numpy()
    within, between = [], []
    n = len(M_ord)
    for i in range(n):
        for j in range(i + 1, n):
            (within if labels[i] == labels[j] else between).append(M_ord.iat[i, j])
    within = np.array(within, dtype=float)
    between = np.array(between, dtype=float)

    w_med = float(np.nanmedian(within)) if within.size else float("nan")
    b_med = float(np.nanmedian(between)) if between.size else float("nan")
    delta = w_med - b_med if np.isfinite(w_med) and np.isfinite(b_med) else float("nan")

    plt.figure(figsize=(4.8, 4.6), dpi=200)
    plt.boxplot([within, between], labels=["Within FINAL", "Between FINAL"], showfliers=False)
    plt.ylabel("Co-assignment fraction")
    plt.title(f"Within vs Between (Δ={delta:.2f}; med_w={w_med:.2f}, med_b={b_med:.2f})")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "within_vs_between_boxplot.png"))
    plt.savefig(os.path.join(out_dir, "within_vs_between_boxplot.pdf"))
    plt.close()

    # ---------- 3) Per-sample agreement ----------
    agree_vals = []
    for i, s in enumerate(M_ord.index):
        same = (final_ord == final_ord.loc[s]).to_numpy()
        same[i] = False
        val = np.nan if same.sum() == 0 else M_ord.values[i, same].mean()
        agree_vals.append(val)
    agree = pd.Series(agree_vals, index=M_ord.index, name="agreement_fraction")
    agree.to_csv(os.path.join(out_dir, "per_sample_agreement_fraction.csv"))

    thr = 0.70
    pct_ge = float((agree >= thr).mean() * 100.0)

    plt.figure(figsize=(5.6, 4.6), dpi=200)
    plt.hist(agree.dropna().values, bins=20, range=(0, 1))
    plt.xlabel("Per-sample agreement fraction")
    plt.ylabel("Count")
    plt.title(f"Agreement histogram (≥{thr:.2f}: {pct_ge:.1f}% of samples)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "agreement_hist.png"))
    plt.savefig(os.path.join(out_dir, "agreement_hist.pdf"))
    plt.close()

    # ---------- 4) Summary CSV ----------
    summary = pd.Series(
        {
            "n_samples": int(len(M_ord)),
            "median_within": w_med,
            "median_between": b_med,
            "delta_medians": delta,
            "pct_agreement_ge_0.70": pct_ge,
        },
        name="summary",
    )
    summary.to_csv(os.path.join(out_dir, "stability_summary.csv"))

    print("✅ Wrote figures & summaries to:", os.path.abspath(out_dir))
    print("    Files:")
    for f in [
        "coassignment_heatmap.png",
        "coassignment_heatmap.pdf",
        "within_vs_between_boxplot.png",
        "within_vs_between_boxplot.pdf",
        "agreement_hist.png",
        "agreement_hist.pdf",
        "pairwise_coassignment_ordered_by_FINAL.csv",
        "per_sample_agreement_fraction.csv",
        "stability_summary.csv",
    ]:
        print("   -", os.path.join(out_dir, f))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="QC plots for consensus stability from pairwise co-assignment matrix."
    )
    ap.add_argument("--pairwise", required=True, help="Path to pairwise_coassignment.csv")
    ap.add_argument("--final", required=True, help="Path to final_assignments.csv")
    ap.add_argument("--out", default="consensus_qc", help="Output directory")
    ap.add_argument("--id-col", default=None, help="ID column name in final_assignments.csv (optional)")
    ap.add_argument("--label-col", default=None, help="Label column name in final_assignments.csv (optional)")
    args = ap.parse_args()

    main(args.pairwise, args.final, args.out, id_col=args.id_col, label_col=args.label_col)
