#!/usr/bin/env python3
"""
scatter_with_gene_labels.py
===========================
Generate annotated scatter plots comparing DESeq2 vs MAST log2 fold-change
estimates for matching contrasts and cell types.  Points are colored by
significance category (both, DESeq2-only, MAST-only, non-sig) and gene
names are optionally labeled.

Inputs
------
- Directory of DESeq2 result CSVs  (--deseq_dir, default: deseq2/)
- Directory of MAST result CSVs    (--mast_dir,  default: mast/)

Outputs
-------
- Per-contrast/celltype scatter PNGs in --outdir (default: scatter_annotated/)

Usage
-----
    python scatter_with_gene_labels.py --deseq_dir deseq2 --mast_dir mast
    python scatter_with_gene_labels.py --contrast dementia_vs_control --ct exc
    python scatter_with_gene_labels.py --label_mode topN --top_n 30
"""

import os
import re
import glob
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Helpers: normalize inputs -----------------------------------------------

CELLTYPE_TOKENS = {"ast", "exc", "inh", "mic", "odc", "opc", "endo", "oligo"}
_fn_end_ct_re = re.compile(r"(.+)_([A-Za-z0-9]+)\.csv$", re.IGNORECASE)


def _first_present(cols, names):
    """Return the first column name from *names* that exists in *cols*."""
    for n in names:
        if n in cols:
            return n
    return None


def load_results_csv(fpath):
    """Load one results CSV and standardize columns: gene, log2FC, padj, pvalue."""
    df = pd.read_csv(fpath)

    # Promote first column to 'gene' if needed
    if 'gene' not in df.columns:
        c0 = df.columns[0]
        guess = {"Gene": "gene", "SYMBOL": "gene", "index": "gene"}
        df = df.rename(columns=guess)
        if 'gene' not in df.columns:
            if str(c0).lower().startswith("unnamed") or str(c0).lower() in {"", "x"}:
                df = df.rename(columns={c0: "gene"})
    if 'gene' not in df.columns:
        df = df.reset_index().rename(columns={"index": "gene"})

    ren = {}
    lfc = _first_present(df.columns,
                         ["log2FoldChange", "avg_log2FC", "avg_logFC", "log2FC", "logFC"])
    padj = _first_present(df.columns,
                          ["padj", "p_adj", "p_val_adj", "qval", "FDR"])
    pval = _first_present(df.columns,
                          ["pvalue", "p_val", "p"])
    if lfc:
        ren[lfc] = "log2FC"
    if padj:
        ren[padj] = "padj"
    if pval:
        ren[pval] = "pvalue"
    df = df.rename(columns=ren)

    for c in ["log2FC", "padj", "pvalue"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df["gene_key"] = df["gene"].astype(str).str.upper()
    df = df.drop_duplicates(subset="gene_key", keep="first")
    return df


def parse_name(fname, method_prefixes=("deseq2_", "filtered_", "mast_")):
    """Return (contrast, celltype) from filename patterns."""
    base = Path(fname).name
    for p in method_prefixes:
        if base.lower().startswith(p):
            base = base[len(p):]
            break
    m = _fn_end_ct_re.match(base)
    if m:
        left, last = m.group(1), m.group(2)
        if last.lower() in CELLTYPE_TOKENS:
            return (left.lower(), last.lower())
        parts = left.split("_")
        if parts[0].lower() in CELLTYPE_TOKENS:
            ct = parts[0].lower()
            contrast = "_".join(parts[1:] + [last])
            return (contrast.lower(), ct)
    parts = base.replace(".csv", "").split("_")
    if parts[0].lower() in CELLTYPE_TOKENS:
        return ("_".join(parts[1:]).lower(), parts[0].lower())
    return (base.replace(".csv", "").lower(), "unknown")


# --- Plotting core -----------------------------------------------------------

def make_scatter(de, ms, label, out_png, padj_thr=0.05, label_mode="both_sig",
                 top_n=50, min_abs_lfc_for_label=0.0, show_categories=True):
    """
    Create an annotated scatter of DESeq2 vs MAST log2FC values.

    Parameters
    ----------
    label_mode : str
        'all'       - label every point
        'both_sig'  - label genes significant in both methods (default)
        'topN'      - label top-N by combined significance
    """
    merged = (de[["gene", "gene_key", "log2FC", "padj"]]
              .rename(columns={"log2FC": "log2FC_DESeq2", "padj": "padj_DESeq2"})
              .merge(ms[["gene_key", "log2FC", "padj"]],
                     on="gene_key", how="inner")
              .rename(columns={"log2FC": "log2FC_MAST", "padj": "padj_MAST"}))
    if merged.empty:
        print(f"[{label}] No genes in common -- skipping.")
        return

    # Significance categories
    sig_de = (merged["padj_DESeq2"].notna() & (merged["padj_DESeq2"] <= padj_thr))
    sig_ms = (merged["padj_MAST"].notna() & (merged["padj_MAST"] <= padj_thr))
    both_sig = sig_de & sig_ms
    de_only = sig_de & ~sig_ms
    ms_only = sig_ms & ~sig_de

    x = merged["log2FC_DESeq2"]
    y = merged["log2FC_MAST"]

    # Symmetric axis limits
    lim = np.nanmax(np.abs(pd.concat([x, y])))
    if not np.isfinite(lim) or lim == 0:
        lim = 1.0
    lim = float(np.ceil(lim * 1.05))

    fig, ax = plt.subplots(figsize=(7.5, 6))

    # Quadrant shading
    ax.fill_between([-lim, 0], 0, lim, alpha=0.08)
    ax.fill_between([0, lim], 0, lim, alpha=0.08)
    ax.fill_between([-lim, 0], -lim, 0, alpha=0.08)
    ax.fill_between([0, lim], -lim, 0, alpha=0.08)

    # Reference lines
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.axhline(0, linestyle="--", linewidth=1)
    ax.plot([-lim, lim], [-lim, lim], linewidth=1)

    # Scatter points
    if show_categories:
        ax.scatter(x[~(sig_de | sig_ms)], y[~(sig_de | sig_ms)], s=10, label="non-sig")
        ax.scatter(x[de_only], y[de_only], s=14, label="DESeq2-only")
        ax.scatter(x[ms_only], y[ms_only], s=14, label="MAST-only")
        ax.scatter(x[both_sig], y[both_sig], s=18, label="both sig")
        ax.legend(frameon=False)
    else:
        ax.scatter(x, y, s=10)

    ax.set_xlabel("log2FC (DESeq2)")
    ax.set_ylabel("log2FC (MAST)")
    ax.set_title(label.replace("_", " "))

    # Decide which points to label
    to_label = pd.Series(False, index=merged.index)
    if label_mode == "all":
        to_label[:] = True
    elif label_mode == "both_sig":
        to_label = both_sig
        if min_abs_lfc_for_label > 0:
            to_label &= ((merged["log2FC_DESeq2"].abs().clip(lower=0) >= min_abs_lfc_for_label) |
                          (merged["log2FC_MAST"].abs().clip(lower=0) >= min_abs_lfc_for_label))
    elif label_mode == "topN":
        best_p = merged[["padj_DESeq2", "padj_MAST"]].min(axis=1)
        mag = merged[["log2FC_DESeq2", "log2FC_MAST"]].abs().max(axis=1)
        order = (best_p.rank(method="first") +
                 (mag.rank(ascending=False, method="first") / 10_000)).sort_values()
        to_label.loc[order.index[:top_n]] = True
    else:
        raise ValueError("label_mode must be one of: all, both_sig, topN")

    if min_abs_lfc_for_label > 0 and label_mode != "both_sig":
        to_label &= (merged[["log2FC_DESeq2", "log2FC_MAST"]].abs().max(axis=1)
                      >= min_abs_lfc_for_label)

    # Add text labels (use adjustText if available)
    use_adjust = False
    try:
        from adjustText import adjust_text
        use_adjust = True
    except Exception:
        pass

    texts = []
    for _, row in merged[to_label].iterrows():
        texts.append(ax.text(row["log2FC_DESeq2"], row["log2FC_MAST"],
                             str(row["gene"]), fontsize=7))

    if use_adjust and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", lw=0.5),
                    only_move={'points': 'y', 'text': 'y'})

    # Quadrant annotations
    ax.text(0.75 * lim, 0.85 * lim, "up in both", ha="center", va="center")
    ax.text(-0.75 * lim, -0.85 * lim, "down in both", ha="center", va="center")
    ax.text(-0.75 * lim, 0.85 * lim, "discordant", ha="center", va="center")
    ax.text(0.75 * lim, -0.85 * lim, "discordant", ha="center", va="center")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"Wrote {out_png} (labeled {to_label.sum()} genes of {len(merged)})")


def index_files(root, method):
    """Index CSV files in a directory by (contrast, celltype) tuple."""
    idx = {}
    for f in sorted(glob.glob(os.path.join(root, "*.csv"))):
        contrast, ct = parse_name(
            f, method_prefixes=(f"{method}_", f"{method.upper()}_",
                                f"{method.lower()}_", "filtered_", "mast_"))
        idx[(contrast, ct)] = f
    return idx


# --- CLI entry point ---------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Annotated DESeq2 vs MAST scatter with quadrants and gene labels.")
    ap.add_argument("--deseq_dir", default="deseq2")
    ap.add_argument("--mast_dir", default="mast")
    ap.add_argument("--contrast", default=None,
                    help="Contrast name (e.g. dementia_vs_control). Omit to do all.")
    ap.add_argument("--ct", default=None,
                    help="Cell type token (e.g. ast, exc). Omit to do all.")
    ap.add_argument("--padj", type=float, default=0.05,
                    help="Adjusted p-value threshold for significance.")
    ap.add_argument("--label_mode", choices=["all", "both_sig", "topN"],
                    default="both_sig")
    ap.add_argument("--top_n", type=int, default=50,
                    help="Top-N to label when --label_mode=topN")
    ap.add_argument("--min_abs_lfc_for_label", type=float, default=0.0,
                    help="Only label genes with max(|lfc|) >= this")
    ap.add_argument("--no_categories", action="store_true",
                    help="Do not color by significance categories")
    ap.add_argument("--outdir", default="scatter_annotated")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    de_idx = index_files(args.deseq_dir, "deseq2")
    ms_idx = index_files(args.mast_dir, "mast")

    # Select keys to render
    keys = sorted(set(de_idx.keys()) & set(ms_idx.keys()))
    if args.contrast:
        keys = [k for k in keys if k[0] == args.contrast.lower()]
    if args.ct:
        keys = [k for k in keys if k[1] == args.ct.lower()]

    if not keys:
        raise SystemExit(
            "No matching (contrast, celltype) pairs found. "
            "Check names or use --contrast/--ct.")

    for (contrast, ct) in keys:
        de = load_results_csv(de_idx[(contrast, ct)])
        ms = load_results_csv(ms_idx[(contrast, ct)])
        label = f"{ct} - {contrast}"
        out_png = os.path.join(args.outdir,
                               f"scatter_annotated_{ct}_{contrast}.png")
        make_scatter(de, ms, label, out_png,
                     padj_thr=args.padj,
                     label_mode=args.label_mode,
                     top_n=args.top_n,
                     min_abs_lfc_for_label=args.min_abs_lfc_for_label,
                     show_categories=not args.no_categories)


if __name__ == "__main__":
    main()
