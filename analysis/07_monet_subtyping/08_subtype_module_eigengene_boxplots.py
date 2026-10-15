#!/usr/bin/env python3
"""hdWGCNA module eigengenes across Control, MCI and AD1-AD4 (Fig. 8a).

For each cell type, merges the per-sample module eigengene table with the
sample-level group labels (Control / MCI / consensus subtypes 1-4) and draws
box + strip plots of every module by group. The published Fig. 8a panel shows
a curated subset of these modules (ASC-M5/M10, MIC-M5/M4, ODC-M1,
EXC-M1/M2, INH-M7/M10).

Usage:
    python 08_subtype_module_eigengene_boxplots.py [--outdir figures/subtype_module_eigengenes]

Inputs (working directory):
    {cell}_MEs_by_sample_manual.csv          cell in mic, ast, exc, inh, odc;
                                             columns: Sample + one column per module
    20_seeds_final_subgroups_with_meta.csv   columns: Samples, Cluster_or_dx
                                             (values Control / MCI / 1.0-4.0)

Outputs:
    {outdir}/final_clustering_20seeds_{cell}_MLP_eigengene_by_module_boxplot.pdf
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

CELL_TYPES = ["mic", "ast", "exc", "inh", "odc"]
EIG_FMT = "{cell}_MEs_by_sample_manual.csv"
SUBGROUPS_CSV = "20_seeds_final_subgroups_with_meta.csv"
GROUPS = ["Control", "MCI", "AD1", "AD2", "AD3", "AD4"]


def norm_id(s: pd.Series) -> pd.Series:
    """Harmonize sample IDs between the eigengene and metadata tables."""
    return (s.astype(str).str.strip()
            .str.replace(r"\.bam$", "", regex=True)
            .str.replace(r"_R\d$", "", regex=True)
            .str.replace(r"\s+", "", regex=True)
            .str.replace("-", "_", regex=False))


def to_group(series: pd.Series) -> pd.Series:
    """Map Control/MCI and 1.0-4.0 to Control/MCI/AD1-AD4."""
    s = series.astype(str).str.strip().str.replace(".0", "", regex=False)
    return s.replace({
        "control": "Control", "CONTROL": "Control",
        "mci": "MCI", "Mci": "MCI",
        "Not Applicable": np.nan, "NA": np.nan, "nan": np.nan,
        "1": "AD1", "2": "AD2", "3": "AD3", "4": "AD4",
    })


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--outdir", default="figures/subtype_module_eigengenes")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = pd.read_csv(SUBGROUPS_CSV)[["Samples", "Cluster_or_dx"]].copy()
    meta["Sample"] = norm_id(meta["Samples"])
    meta["Group"] = pd.Categorical(to_group(meta["Cluster_or_dx"]), categories=GROUPS, ordered=True)
    meta = meta.dropna(subset=["Group"]).drop_duplicates("Sample")
    print("Samples per group:\n" + meta["Group"].value_counts().to_string())

    palette = sns.color_palette("Paired", len(GROUPS))
    hue_kws = dict(hue="Group", hue_order=GROUPS, palette=palette)

    for cell in CELL_TYPES:
        eig = pd.read_csv(EIG_FMT.format(cell=cell)).drop(columns=["Unnamed: 0"], errors="ignore")
        eig["Sample"] = norm_id(eig["Sample"])
        merged = eig.merge(meta[["Sample", "Group"]], on="Sample", how="inner")
        if merged.empty:
            print(f"{cell}: no overlapping samples, skipped")
            continue

        modules = [c for c in merged.columns if c not in {"Sample", "Group"}]
        long_df = merged.melt(id_vars=["Sample", "Group"], value_vars=modules,
                              var_name="Module", value_name="Eigengene")
        long_df["Eigengene"] = pd.to_numeric(long_df["Eigengene"], errors="coerce")
        long_df = long_df.dropna(subset=["Eigengene"])

        plt.figure(figsize=(14, 6))
        ax = sns.boxplot(data=long_df, x="Module", y="Eigengene", **hue_kws)
        sns.stripplot(data=long_df, x="Module", y="Eigengene", dodge=True, alpha=0.5,
                      jitter=0.25, size=3, linewidth=0, legend=False, **hue_kws)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:len(GROUPS)], labels[:len(GROUPS)], title="Group",
                  bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.title(f"{cell.upper()} module eigengene by group")
        plt.xlabel("Module"); plt.ylabel("Eigengene expression")
        plt.xticks(rotation=45, ha="right"); plt.tight_layout()
        out = os.path.join(args.outdir, f"final_clustering_20seeds_{cell}_MLP_eigengene_by_module_boxplot.pdf")
        plt.savefig(out); plt.close()
        print(f"{cell}: {merged.shape[0]} samples, {len(modules)} modules -> {out}")


if __name__ == "__main__":
    main()
