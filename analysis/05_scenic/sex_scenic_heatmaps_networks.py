#!/usr/bin/env python3
"""
sex_scenic_heatmaps_networks.py
===============================
pySCENIC Male vs Female AD: Regulon heatmaps and TF-TG network plots.

Steps:
  1. Per-cell-type regulon Z-score heatmaps (M vs F side-by-side)
  2. Differential regulon heatmap (M - F)
  3. TF-TG network plots from pySCENIC adjacencies (per sex)
  4. Sex-specific TF-TG edge identification
  5. TF importance comparison bar chart

Inputs:
    - AUCell matrices:    <SCENIC_DIR>/dm_{male,female}_auc_matrix.csv
    - h5ad metadata:      <H5AD_DIR>/13adata_dem_{male,female}.h5ad
    - pySCENIC adjacencies: <SCENIC_DIR>/adata_dm_{male,female}_adj.csv

Outputs:
    - Regulon Z-score clustermaps (per sex and combined)
    - Differential regulon heatmap (Male - Female)
    - TF-TG network PDFs (overall and sex-specific)
    - TF importance sex-bias bar chart
    - CSV tables of sex-specific edges and TF importance

Usage:
    python sex_scenic_heatmaps_networks.py
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.lines import Line2D

# --- Configuration -----------------------------------------------------------

BASE_DIR = Path(".")
SCENIC_DIR = BASE_DIR / "scenic"
H5AD_DIR = BASE_DIR / "h5ad"

OUT = SCENIC_DIR / "sex_comparison"
OUT.mkdir(exist_ok=True)

# AUCell matrices
AUC_MALE = SCENIC_DIR / "dm_male_auc_matrix.csv"
AUC_FEMALE = SCENIC_DIR / "dm_female_auc_matrix.csv"

# h5ad for cell type annotations
H5AD_MALE = H5AD_DIR / "13adata_dem_male.h5ad"
H5AD_FEMALE = H5AD_DIR / "13adata_dem_female.h5ad"

# pySCENIC adjacencies
ADJ_MALE = SCENIC_DIR / "adata_dm_male_adj.csv"
ADJ_FEMALE = SCENIC_DIR / "adata_dm_female_adj.csv"

# Cell types to analyze
CELLTYPES = ["Excitatory", "Inhibitory", "Astrocyte", "Microglia", "Oligodendrocyte"]
CELLTYPE_COL = "Primary_Function"

# --- Load data ---------------------------------------------------------------

print("Loading AUCell matrices...")
auc_m = pd.read_csv(AUC_MALE, index_col=0)
auc_f = pd.read_csv(AUC_FEMALE, index_col=0)
print(f"  Male:   {auc_m.shape[0]} cells, {auc_m.shape[1]} regulons")
print(f"  Female: {auc_f.shape[0]} cells, {auc_f.shape[1]} regulons")

print("Loading h5ad metadata...")
adata_m = sc.read_h5ad(str(H5AD_MALE), backed="r")
adata_f = sc.read_h5ad(str(H5AD_FEMALE), backed="r")

# Extract cell type info
meta_m = adata_m.obs[[CELLTYPE_COL]].copy()
meta_f = adata_f.obs[[CELLTYPE_COL]].copy()

# Align AUC barcodes with metadata
common_m = auc_m.index.intersection(meta_m.index)
common_f = auc_f.index.intersection(meta_f.index)
print(f"  Male cells matched:   {len(common_m)}/{auc_m.shape[0]}")
print(f"  Female cells matched: {len(common_f)}/{auc_f.shape[0]}")

auc_m = auc_m.loc[common_m]
meta_m = meta_m.loc[common_m]
auc_f = auc_f.loc[common_f]
meta_f = meta_f.loc[common_f]

# Shared regulons
shared_regulons = sorted(set(auc_m.columns) & set(auc_f.columns))
print(f"  Shared regulons: {len(shared_regulons)}")

# --- Analysis: per-cell-type regulon Z-scores --------------------------------

print("\n=== Computing per-cell-type regulon Z-scores ===")


def compute_celltype_zscores(auc_df, meta_df, celltypes, regulons):
    """Compute mean AUC per cell type, then Z-score across cell types."""
    results = {}
    for ct in celltypes:
        mask = meta_df[CELLTYPE_COL] == ct
        if mask.sum() < 10:
            continue
        ct_auc = auc_df.loc[mask, regulons]
        results[ct] = ct_auc.mean(axis=0)

    mean_df = pd.DataFrame(results)  # regulons x celltypes
    # Z-score each regulon across cell types
    z_df = mean_df.apply(
        lambda row: (row - row.mean()) / row.std() if row.std() > 0 else row * 0,
        axis=1,
    )
    return mean_df, z_df


mean_m, z_m = compute_celltype_zscores(auc_m, meta_m, CELLTYPES, shared_regulons)
mean_f, z_f = compute_celltype_zscores(auc_f, meta_f, CELLTYPES, shared_regulons)

# Save Z-score tables
z_m.to_csv(OUT / "male_regulon_zscore_by_celltype.csv")
z_f.to_csv(OUT / "female_regulon_zscore_by_celltype.csv")

# --- Heatmaps ---------------------------------------------------------------

print("\n=== Generating heatmaps ===")

# Side-by-side Z-score clustermaps (M vs F)
for sex_label, z_df in [("Male", z_m), ("Female", z_f)]:
    z_sorted = z_df.loc[z_df.mean(axis=1).sort_values(ascending=False).index]

    g = sns.clustermap(
        z_sorted,
        cmap="vlag",
        figsize=(8, max(6, len(z_sorted) * 0.18)),
        yticklabels=True,
        col_cluster=True,
        row_cluster=True,
        vmin=-2, vmax=2,
        linewidths=0.2,
        cbar_kws={"label": "Z-score"},
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=5)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=9)
    g.fig.suptitle(f"Regulon Activity Z-scores — {sex_label} AD", y=1.02, fontsize=14)

    outpath = OUT / f"regulon_zscore_clustermap_{sex_label}.pdf"
    g.savefig(outpath, dpi=300, bbox_inches="tight")
    g.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath.name}")

# Differential heatmap (Male - Female)
diff_z = z_m - z_f
diff_z = diff_z.dropna(how="all")
diff_sorted = diff_z.loc[diff_z.abs().mean(axis=1).sort_values(ascending=False).index]

# Top 40 most differential
top_diff = diff_sorted.head(40)

fig, ax = plt.subplots(figsize=(8, max(6, len(top_diff) * 0.25)))
vmax = max(abs(top_diff.values.min()), abs(top_diff.values.max()), 1)
sns.heatmap(
    top_diff, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
    yticklabels=True, xticklabels=True, linewidths=0.3,
    cbar_kws={"label": "Z-score difference (Male - Female)"},
    ax=ax,
)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=9, rotation=45, ha="right")
ax.set_title("Differential Regulon Activity (Male − Female AD)\nTop 40 most variable", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "regulon_diff_heatmap_MvsF.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "regulon_diff_heatmap_MvsF.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: regulon_diff_heatmap_MvsF.pdf")

# Combined heatmap (M and F side by side per cell type)
combined_cols = {}
for ct in CELLTYPES:
    if ct in z_m.columns and ct in z_f.columns:
        combined_cols[f"{ct}\nMale"] = z_m[ct]
        combined_cols[f"{ct}\nFemale"] = z_f[ct]

combined_df = pd.DataFrame(combined_cols)
combined_sorted = combined_df.loc[combined_df.abs().mean(axis=1).sort_values(ascending=False).index]
top_combined = combined_sorted.head(50)

g2 = sns.clustermap(
    top_combined,
    cmap="vlag",
    figsize=(12, max(8, len(top_combined) * 0.2)),
    yticklabels=True,
    col_cluster=False,
    row_cluster=True,
    vmin=-2, vmax=2,
    linewidths=0.2,
    cbar_kws={"label": "Z-score"},
)
g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), fontsize=5)
g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), fontsize=8, rotation=45, ha="right")
g2.fig.suptitle("Regulon Activity: Male vs Female AD (top 50 regulons)", y=1.02, fontsize=14)

g2.savefig(OUT / "regulon_zscore_combined_MvsF.pdf", dpi=300, bbox_inches="tight")
g2.savefig(OUT / "regulon_zscore_combined_MvsF.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: regulon_zscore_combined_MvsF.pdf")

# --- TF-TG network plots (pySCENIC adjacencies) -----------------------------

print("\n=== Building TF-TG network plots ===")

NET_OUT = OUT / "networks"
NET_OUT.mkdir(exist_ok=True)

# Load adjacencies
adj_m = pd.read_csv(ADJ_MALE)
adj_f = pd.read_csv(ADJ_FEMALE)
print(f"  Male adjacencies:   {len(adj_m)}")
print(f"  Female adjacencies: {len(adj_f)}")

# Clean column names
for df in [adj_m, adj_f]:
    df.columns = [c.strip() for c in df.columns]


def build_tf_tg_network(adj_df, top_n=200, title="", outpath="network.pdf"):
    """Build and plot a TF-TG network from pySCENIC adjacencies using networkx."""

    top_edges = adj_df.nlargest(top_n, "importance").copy()
    if top_edges.empty:
        print(f"  No edges for {title}, skipping.")
        return

    tfs = set(top_edges["TF"].unique())
    tgs = set(top_edges["target"].unique()) - tfs

    G = nx.DiGraph()
    for node in sorted(tfs):
        G.add_node(node, node_type="TF")
    for node in sorted(tgs):
        G.add_node(node, node_type="TG")

    w_arr = top_edges["importance"].values
    w_min, w_max = w_arr.min(), w_arr.max()

    for _, row in top_edges.iterrows():
        w_norm = (row["importance"] - w_min) / (w_max - w_min) if w_max > w_min else 1.0
        G.add_edge(row["TF"], row["target"], weight=row["importance"], w_norm=w_norm)

    # Layout
    pos = nx.spring_layout(G, k=1.8 / np.sqrt(len(G.nodes)), iterations=80, seed=42)

    fig, ax = plt.subplots(figsize=(16, 12))

    # Draw TG nodes
    tg_nodes = [n for n in G.nodes if G.nodes[n]["node_type"] == "TG"]
    nx.draw_networkx_nodes(
        G, pos, nodelist=tg_nodes, node_color="#80B1D3",
        node_size=80, node_shape="d", edgecolors="#333",
        linewidths=0.5, ax=ax, alpha=0.8,
    )

    # Draw TF nodes
    tf_nodes = [n for n in G.nodes if G.nodes[n]["node_type"] == "TF"]
    nx.draw_networkx_nodes(
        G, pos, nodelist=tf_nodes, node_color="#ffe666",
        node_size=350, node_shape="o", edgecolors="#333",
        linewidths=1.5, ax=ax,
    )

    # Draw edges with varying alpha/width
    for u, v, data in G.edges(data=True):
        w = data["w_norm"]
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)],
            width=1 + w * 5,
            alpha=0.15 + 0.85 * w,
            edge_color=[(1.0, 0.35, 0.15)],
            arrows=True, arrowsize=8, ax=ax,
            connectionstyle="arc3,rad=0.05",
        )

    # Labels
    nx.draw_networkx_labels(
        G, pos, labels={n: n for n in tf_nodes},
        font_size=8, font_weight="bold", ax=ax,
    )
    nx.draw_networkx_labels(
        G, pos, labels={n: n for n in tg_nodes},
        font_size=4, font_color="#333", ax=ax,
    )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ffe666",
               markeredgecolor="k", markersize=12, label="TF"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#80B1D3",
               markeredgecolor="k", markersize=8, label="Target gene"),
        Line2D([0], [0], color=(1, 0.35, 0.15), lw=3,
               label="TF->TG (width/alpha ~ importance)"),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=9, frameon=True)

    outpath = Path(outpath)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath.name}")

    return G


# Build networks for each sex (top 200 overall edges)
for sex_label, adj_df in [("Male", adj_m), ("Female", adj_f)]:
    build_tf_tg_network(
        adj_df, top_n=200,
        title=f"Top TF->TG Network — {sex_label} AD (pySCENIC)",
        outpath=str(NET_OUT / f"tf_tg_network_{sex_label}_top200.pdf"),
    )

# --- Identify sex-specific TF-TG interactions -------------------------------

print("\n=== Identifying sex-specific TF-TG interactions ===")

# Top 500 per sex
top_m = adj_m.nlargest(500, "importance")[["TF", "target"]].assign(
    pair=lambda x: x["TF"] + "->" + x["target"]
)
top_f = adj_f.nlargest(500, "importance")[["TF", "target"]].assign(
    pair=lambda x: x["TF"] + "->" + x["target"]
)

male_only = set(top_m["pair"]) - set(top_f["pair"])
female_only = set(top_f["pair"]) - set(top_m["pair"])
shared = set(top_m["pair"]) & set(top_f["pair"])

print(f"  Top 500 edges — Male-only: {len(male_only)} | Female-only: {len(female_only)} | Shared: {len(shared)}")

# Save sex-specific edges
male_specific = adj_m[adj_m.apply(lambda r: f"{r['TF']}->{r['target']}" in male_only, axis=1)]
female_specific = adj_f[adj_f.apply(lambda r: f"{r['TF']}->{r['target']}" in female_only, axis=1)]

male_specific.to_csv(OUT / "male_specific_top_TF_TG.csv", index=False)
female_specific.to_csv(OUT / "female_specific_top_TF_TG.csv", index=False)

# Build sex-specific networks
if len(male_only) > 0:
    male_spec_adj = adj_m[adj_m.apply(lambda r: f"{r['TF']}->{r['target']}" in male_only, axis=1)]
    build_tf_tg_network(
        male_spec_adj, top_n=min(100, len(male_spec_adj)),
        title="Male-Specific TF->TG Network (AD)",
        outpath=str(NET_OUT / "tf_tg_network_Male_specific.pdf"),
    )

if len(female_only) > 0:
    female_spec_adj = adj_f[adj_f.apply(lambda r: f"{r['TF']}->{r['target']}" in female_only, axis=1)]
    build_tf_tg_network(
        female_spec_adj, top_n=min(100, len(female_spec_adj)),
        title="Female-Specific TF->TG Network (AD)",
        outpath=str(NET_OUT / "tf_tg_network_Female_specific.pdf"),
    )

# --- TF importance comparison bar chart --------------------------------------

print("\n=== TF importance comparison ===")

# Sum importance per TF
tf_imp_m = adj_m.groupby("TF")["importance"].sum().rename("Male")
tf_imp_f = adj_f.groupby("TF")["importance"].sum().rename("Female")

tf_imp = pd.concat([tf_imp_m, tf_imp_f], axis=1).fillna(0)
tf_imp["diff"] = tf_imp["Male"] - tf_imp["Female"]
tf_imp = tf_imp.sort_values("diff")

# Top 15 male-biased and female-biased TFs
top_male_biased = tf_imp.tail(15)
top_female_biased = tf_imp.head(15)
plot_df = pd.concat([top_female_biased, top_male_biased])

fig, ax = plt.subplots(figsize=(8, 8))
colors = ["#E8938A" if d < 0 else "#7BAFD4" for d in plot_df["diff"]]
ax.barh(plot_df.index, plot_df["diff"], color=colors, edgecolor="grey", linewidth=0.3)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("Sum Importance (Male) - Sum Importance (Female)", fontsize=10)
ax.set_title("Sex-Biased Transcription Factors in AD\n(pySCENIC GRN importance)", fontsize=12)
ax.tick_params(axis="y", labelsize=8)

legend_elems = [
    Line2D([0], [0], color="#7BAFD4", lw=8, label="Male-biased"),
    Line2D([0], [0], color="#E8938A", lw=8, label="Female-biased"),
]
ax.legend(handles=legend_elems, loc="lower right", fontsize=9)

fig.tight_layout()
fig.savefig(OUT / "TF_importance_sex_bias.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "TF_importance_sex_bias.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: TF_importance_sex_bias.pdf")

# --- Save results ------------------------------------------------------------

tf_imp.to_csv(OUT / "TF_importance_comparison.csv")

# --- Summary -----------------------------------------------------------------

print("\n" + "=" * 60)
print("Done! Outputs in:", OUT)
print("  - regulon_zscore_clustermap_Male/Female.pdf")
print("  - regulon_diff_heatmap_MvsF.pdf")
print("  - regulon_zscore_combined_MvsF.pdf")
print("  - networks/tf_tg_network_Male/Female_top200.pdf")
print("  - networks/tf_tg_network_Male/Female_specific.pdf")
print("  - TF_importance_sex_bias.pdf")
print("=" * 60)
