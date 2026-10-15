"""
UpSet plots of subtype DEGs (Fig. 7d): sets = MCI, AD1-AD4 (each vs. Control),
one page per cell type, for all / up / down DEGs.

Filter: MAST p_val_adj < 0.05 and |avg_log2FC| > 1; intersections with fewer
than 5 genes are not shown. PDFs use pdf.fonttype=42 (editable text).

Usage (run from the MAST subtype DE output directory):
    python plot_upset_subtypes.py

Inputs (working directory):
    {MCI,AD1,AD2,AD3,AD4}_vs_Control/{celltype}_*.csv   MAST FindMarkers output
                                                        (from 01_deg_analysis/subtype_deg*.R)
Outputs:
    plots/upset_subtypes/upset_subtypes_{all,up,down}.pdf
"""
import os
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType — editable in Adobe AI
matplotlib.rcParams['ps.fonttype']  = 42
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

BASE = "."
OUT  = os.path.join(BASE, "plots", "upset_subtypes")
os.makedirs(OUT, exist_ok=True)

COMPARISONS = ["MCI_vs_Control", "AD1_vs_Control", "AD2_vs_Control", "AD3_vs_Control", "AD4_vs_Control"]
COMP_LABELS  = ["MCI", "AD1", "AD2", "AD3", "AD4"]

CT_ORDER  = ["ast", "exc", "inh", "mic", "odc", "opc"]   # endo excluded
CT_LABELS = {
    "ast": "Astrocytes",
    "exc": "Exc. Neurons",
    "inh": "Inh. Neurons",
    "mic": "Microglia",
    "odc": "Oligodendrocytes",
    "opc": "OPCs",
}

COMP_COLORS = {
    "MCI": "#7570b3",
    "AD1": "#1b9e77",
    "AD2": "#d95f02",
    "AD3": "#e7298a",
    "AD4": "#e41a1c",
}

PADJ = 0.05
LFC  = 1.0


def load_sets_for_celltype(ct_prefix, direction="all"):
    sets = {}
    for comp_dir_name, label in zip(COMPARISONS, COMP_LABELS):
        folder = os.path.join(BASE, comp_dir_name)
        hits = glob.glob(os.path.join(folder, f"{ct_prefix}_*.csv"))
        if not hits:
            continue
        df = pd.read_csv(hits[0], index_col=0)
        sig = df[(df["p_val_adj"] < PADJ) & (df["avg_log2FC"].abs() > LFC)]
        if direction == "up":
            sig = sig[sig["avg_log2FC"] > 0]
        elif direction == "down":
            sig = sig[sig["avg_log2FC"] < 0]
        sets[label] = set(sig.index.tolist())
    return sets


for dirkey, dirlabel in [("all", "All DEGs"), ("up", "Upregulated only"), ("down", "Downregulated only")]:
    pdf_path = os.path.join(OUT, f"upset_subtypes_{dirkey}.pdf")

    with PdfPages(pdf_path) as pdf:
        for ct_prefix in CT_ORDER:
            ct_name = CT_LABELS[ct_prefix]
            sets_clean = {k: v for k, v in load_sets_for_celltype(ct_prefix, direction=dirkey).items()
                          if len(v) > 0}

            if not sets_clean:
                print(f"  Skipping {ct_name} ({dirkey}) — no sets with data")
                continue

            data = from_contents(sets_clean)
            upset = UpSet(
                data,
                subset_size="count",
                min_subset_size=5,
                sort_by="cardinality",
                sort_categories_by=None,
                show_counts=True,
                totals_plot_elements=3,
            )

            fig = plt.figure(figsize=(14, 5))
            fig.patch.set_facecolor("#f9f9f9")
            d = upset.plot(fig=fig)

            for bar, label in zip(d["totals"].patches, list(sets_clean.keys())):
                bar.set_facecolor(COMP_COLORS.get(label, "#888888"))
                bar.set_alpha(0.9)

            title = (f"{ct_name} — {dirlabel}\n"
                     f"Sets = AD subtypes  |  padj < {PADJ}, |log2FC| > {LFC}")
            fig.suptitle(title, fontsize=12, fontweight="bold", y=1.02)

            legend_handles = [Patch(facecolor=COMP_COLORS[k], label=k, alpha=0.9)
                              for k in sets_clean.keys()]
            d["intersections"].legend(
                handles=legend_handles, title="Subtype (set size bar color)",
                fontsize=8, title_fontsize=8.5,
                loc="upper right", framealpha=0.8,
            )

            pdf.savefig(fig, bbox_inches="tight", facecolor=fig.get_facecolor())
            plt.close(fig)
            print(f"  Added page: {ct_name} ({dirkey})")

    print(f"Saved: {pdf_path}")
