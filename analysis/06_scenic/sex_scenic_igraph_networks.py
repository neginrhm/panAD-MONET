"""
sex_scenic_igraph_networks.py
=============================
TF x TG heatmap + TF->TG network (igraph) for Male vs Female AD using
pySCENIC adjacency data (importance scores).

Generates per-sex heatmaps and igraph networks, plus a comparison network
showing sex-specific edges. Optionally highlights GWAS genes.

Inputs:
    - pySCENIC adjacencies: <SCENIC_DIR>/adata_dm_{male,female}_adj.csv
    - GWAS gene list (optional): <SCENIC_DIR>/gwas_genes_ad_ftd.xlsx

Outputs:
    - Per-sex TF x TG heatmaps (PDF + PNG)
    - Per-sex igraph network plots (PDF + PNG) with legends
    - Comparison network (male-specific / female-specific / shared edges)
    - CSV tables of sex-specific and shared top edges

Usage:
    python sex_scenic_igraph_networks.py
"""

import os
from pathlib import Path
from typing import Union, Tuple, Set, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import igraph as ig
from matplotlib.lines import Line2D
from matplotlib.colors import PowerNorm

# --- Configuration -----------------------------------------------------------

BASE_DIR = Path(".")
SCENIC_DIR = BASE_DIR / "scenic"
OUT_DIR = SCENIC_DIR / "sex_comparison" / "igraph_networks"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PATH_MALE_ADJ = SCENIC_DIR / "adata_dm_male_adj.csv"
PATH_FEMALE_ADJ = SCENIC_DIR / "adata_dm_female_adj.csv"

PATH_GWAS = SCENIC_DIR / "gwas_genes_ad_ftd.xlsx"
GWAS_GENE_COL = "GeneName"

# --- Tunable parameters ------------------------------------------------------

TOP_EDGES_PER_SEX = 300
TOP_TF_PER_TG = 3
TOP_TFS_HEATMAP = 15
MIN_WEIGHT_HEATMAP = 0.0

FIGSIZE_HEATMAP = (26, 9)

# --- Colour constants (matching Fig6 pySCENIC/NetRem GRN palette) -----------

COLOR_TG_NODE = "#ff4400"
COLOR_TF_NODE = "#ffe666"

EDGE_BASE_RGB = (1.00, 0.35, 0.05)

# Comparison network: sex-specific edge colours
EDGE_BASE_RGB_MALE = (0.05, 0.80, 0.95)
EDGE_BASE_RGB_FEMALE = (0.60, 0.35, 0.80)
EDGE_BASE_RGB_SHARED = (1.00, 0.35, 0.05)

COLOR_EDGE_MALE = "#0dccf2"
COLOR_EDGE_FEMALE = "#9959cc"
COLOR_EDGE_SHARED = "#ff5900"

# GWAS highlight
GWAS_HIGHLIGHT_COLOR = "#00e5ff"
GWAS_FRAME_WIDTH = 4

# --- Pop / contrast controls -------------------------------------------------

HEATMAP_CMAP = "viridis"
HEATMAP_QUANTILE_CLIP = (0.05, 0.95)
HEATMAP_GAMMA = 0.65

EDGE_ALPHA_MIN = 0.15
EDGE_ALPHA_MAX = 1.0
NODE_ALPHA_MIN = 0.25
NODE_ALPHA_MAX = 1.0
EDGE_ALPHA_GAMMA = 1.5
NODE_ALPHA_GAMMA = 0.6

TF_TG_WIDTH_MIN = 1.0
TF_TG_WIDTH_MAX = 14.0

# --- Helper functions --------------------------------------------------------


def hex_to_rgb(hex_color: str) -> Tuple[float, float, float]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def normalize_0_1(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    if x.dropna().empty:
        return pd.Series(np.zeros(len(x)), index=x.index)
    mn, mx = x.min(), x.max()
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - mn) / (mx - mn)


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def scaled_alpha(w01: float, amin: float, amax: float, gamma: float) -> float:
    w01 = clamp01(w01)
    return amin + (amax - amin) * (w01 ** gamma)


def scaled_width(w01: float, wmin: float, wmax: float, gamma: float) -> float:
    w01 = clamp01(w01)
    return wmin + (wmax - wmin) * (w01 ** gamma)


def top_k_per_group(df: pd.DataFrame, by_col: str, group_col: str, k: int) -> pd.DataFrame:
    if k is None:
        return df
    return (df.sort_values(by_col, ascending=False)
              .groupby(group_col, as_index=False, group_keys=False)
              .head(k))


def edges_to_matrix(edges: pd.DataFrame, row_key: str, col_key: str, weight_key: str,
                    row_order=None, col_order=None):
    rows = row_order if row_order is not None else sorted(edges[row_key].unique())
    cols = col_order if col_order is not None else sorted(edges[col_key].unique())
    r_index = {r: i for i, r in enumerate(rows)}
    c_index = {c: j for j, c in enumerate(cols)}
    M = np.zeros((len(rows), len(cols)))
    for _, r in edges.iterrows():
        i = r_index.get(r[row_key])
        j = c_index.get(r[col_key])
        if i is not None and j is not None:
            M[i, j] = max(M[i, j], float(r[weight_key]))
    return M, rows, cols


def load_gwas_genes(path: Path, gene_col: str) -> Set[str]:
    if not Path(path).exists():
        print(f"[GWAS] file not found: {path}")
        return set()
    df = pd.read_excel(path)
    if gene_col not in df.columns:
        print(f"[GWAS] column '{gene_col}' not found")
        return set()
    genes = set(df[gene_col].dropna().astype(str).str.strip().str.upper().tolist())
    print(f"[GWAS] loaded {len(genes)} genes")
    return genes


def is_gwas_gene(gene: str, gwas_genes_upper: Set[str]) -> bool:
    return str(gene).strip().upper() in gwas_genes_upper


# --- Loader ------------------------------------------------------------------


def load_pyscenic_adj(path: Path) -> pd.DataFrame:
    """Load and normalize a pySCENIC adjacency CSV."""
    df = pd.read_csv(path)
    df = df.rename(columns={"target": "TG"})
    df["TF"] = df["TF"].astype(str).str.strip()
    df["TG"] = df["TG"].astype(str).str.strip()
    df["importance"] = pd.to_numeric(df["importance"], errors="coerce")
    df = (df.groupby(["TF", "TG"], as_index=False)
            .agg(importance=("importance", "max")))
    df["weight"] = normalize_0_1(df["importance"])
    df = df.sort_values("importance", ascending=False).reset_index(drop=True)
    return df


# --- Heatmap -----------------------------------------------------------------


def draw_heatmap(tf_tg: pd.DataFrame, title: str, outpath: Path,
                 gwas_genes_upper: Set[str]) -> pd.DataFrame:
    """Draw a TF x TG importance heatmap with optional GWAS gene highlighting."""
    edges = tf_tg.dropna(subset=["weight"]).copy()
    if MIN_WEIGHT_HEATMAP > 0:
        edges = edges[edges["weight"] >= MIN_WEIGHT_HEATMAP]
    if edges.empty:
        print("[heatmap] no edges -> skipping.")
        return edges

    edges = top_k_per_group(edges, "weight", "TG", TOP_TF_PER_TG)
    if edges.empty:
        return edges

    tf_strength = edges.groupby("TF")["weight"].sum().sort_values(ascending=False)
    keep_tfs = tf_strength.head(TOP_TFS_HEATMAP).index.tolist()
    edges = edges[edges["TF"].isin(keep_tfs)]
    if edges.empty:
        return edges

    M, rows_tf, cols_tg = edges_to_matrix(
        edges, "TF", "TG", "weight", row_order=keep_tfs, col_order=None
    )

    # Waterfall ordering of TGs
    primary_tf_idx = M.argmax(axis=0)
    max_vals = M.max(axis=0)
    col_info = [(j, int(primary_tf_idx[j]), float(max_vals[j])) for j in range(M.shape[1])]
    col_info_sorted = sorted(col_info, key=lambda x: (x[1], -x[2]))
    col_order = [j for (j, _, _) in col_info_sorted]
    M = M[:, col_order]
    cols_tg = [cols_tg[j] for j in col_order]

    flat = M.flatten()
    flat = flat[np.isfinite(flat)]
    if flat.size > 0:
        qlo, qhi = HEATMAP_QUANTILE_CLIP
        vmin = float(np.quantile(flat, qlo))
        vmax = float(np.quantile(flat, qhi))
        if vmax <= vmin:
            vmin, vmax = float(np.min(flat)), float(np.max(flat))
    else:
        vmin, vmax = 0.0, 1.0

    norm = PowerNorm(gamma=HEATMAP_GAMMA, vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=FIGSIZE_HEATMAP)
    im = ax.imshow(M, aspect="auto", cmap=HEATMAP_CMAP, norm=norm)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(
        f"pySCENIC importance (norm 0-1; PowerNorm gamma={HEATMAP_GAMMA}; "
        f"clipped {HEATMAP_QUANTILE_CLIP})"
    )

    ax.set_yticks(range(len(rows_tf)))
    ax.set_yticklabels(rows_tf, fontsize=11)
    ax.set_xticks(range(len(cols_tg)))
    ax.set_xticklabels(cols_tg, rotation=90, fontsize=8)

    gwas_bbox = dict(
        boxstyle="round,pad=0.25", facecolor="none",
        edgecolor=GWAS_HIGHLIGHT_COLOR, linewidth=2.5,
    )

    for lbl, gene in zip(ax.get_yticklabels(), rows_tf):
        if is_gwas_gene(gene, gwas_genes_upper):
            lbl.set_fontweight("bold")
            lbl.set_bbox(gwas_bbox)

    for lbl, gene in zip(ax.get_xticklabels(), cols_tg):
        if is_gwas_gene(gene, gwas_genes_upper):
            lbl.set_fontweight("bold")
            lbl.set_bbox(gwas_bbox)

    ax.set_title(title, fontsize=14)
    fig.tight_layout()

    outpath = Path(outpath)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"[heatmap] saved to {outpath}")
    plt.close(fig)
    return edges.copy()


# --- igraph network ----------------------------------------------------------


def build_igraph(tf_tg_subset: pd.DataFrame,
                 gwas_genes_upper: Set[str],
                 edge_base_rgb: Tuple[float, float, float]) -> ig.Graph:
    """Build an igraph network from TF-TG edges with styled nodes and edges."""

    tg_max_strength: Dict[str, float] = {}
    for _, r in tf_tg_subset.iterrows():
        tg = r["TG"]
        w = float(r.get("weight", 0.0))
        if np.isfinite(w):
            tg_max_strength[tg] = max(tg_max_strength.get(tg, 0.0), w)

    tfs = set(tf_tg_subset["TF"])
    tgs = set(tf_tg_subset["TG"])
    vertices = sorted(tfs | tgs)

    g = ig.Graph()
    g.add_vertices(vertices)

    TF_SIZE = 100
    TG_SIZE = 35
    TF_LABEL_SIZE = 26
    TG_LABEL_SIZE = 22
    TG_LABEL_DIST = 1.3

    TG_BASE_RGB = hex_to_rgb(COLOR_TG_NODE)

    vtype, vshape, vcolor, vsize = [], [], [], []
    label_color, label_size, label_dist, label_angle = [], [], [], []
    frame_color, frame_width = [], []

    for name in vertices:
        hit = is_gwas_gene(name, gwas_genes_upper)

        if name in tfs:
            vtype.append("TF")
            vshape.append("circle")
            vsize.append(TF_SIZE)
            vcolor.append(COLOR_TF_NODE)
            label_color.append("#000000")
            label_size.append(TF_LABEL_SIZE)
            label_dist.append(0.0)
            label_angle.append(0.0)
        else:
            vtype.append("TG")
            vshape.append("diamond")
            vsize.append(TG_SIZE)

            strength = clamp01(tg_max_strength.get(name, 0.0))
            alpha_node = scaled_alpha(strength, NODE_ALPHA_MIN, NODE_ALPHA_MAX, NODE_ALPHA_GAMMA)
            vcolor.append((TG_BASE_RGB[0], TG_BASE_RGB[1], TG_BASE_RGB[2], alpha_node))
            label_color.append("#000000")
            label_size.append(TG_LABEL_SIZE)
            label_dist.append(TG_LABEL_DIST)
            label_angle.append(np.pi / 2)

        if hit:
            frame_color.append(GWAS_HIGHLIGHT_COLOR)
            frame_width.append(GWAS_FRAME_WIDTH)
        else:
            frame_color.append("#00000000")
            frame_width.append(0)

    g.vs["name"] = vertices
    g.vs["type"] = vtype
    g.vs["shape"] = vshape
    g.vs["color"] = vcolor
    g.vs["size"] = vsize
    g.vs["label"] = vertices
    g.vs["label_color"] = label_color
    g.vs["label_size"] = label_size
    g.vs["label_dist"] = label_dist
    g.vs["label_angle"] = label_angle
    g.vs["frame_color"] = frame_color
    g.vs["frame_width"] = frame_width

    idx = {name: i for i, name in enumerate(vertices)}

    edges_list, edge_colors, edge_widths = [], [], []
    for _, r in tf_tg_subset.iterrows():
        tf_name = r["TF"]
        tg_name = r["TG"]
        if tf_name not in idx or tg_name not in idx:
            continue

        edges_list.append((idx[tf_name], idx[tg_name]))
        w = clamp01(float(r.get("weight", 0.5)))
        alpha = scaled_alpha(w, EDGE_ALPHA_MIN, EDGE_ALPHA_MAX, EDGE_ALPHA_GAMMA)
        edge_colors.append((edge_base_rgb[0], edge_base_rgb[1], edge_base_rgb[2], alpha))
        edge_widths.append(scaled_width(w, TF_TG_WIDTH_MIN, TF_TG_WIDTH_MAX, EDGE_ALPHA_GAMMA))

    g.add_edges(edges_list)
    g.es["color"] = edge_colors
    g.es["width"] = edge_widths

    return g


def build_igraph_comparison(male_specific: pd.DataFrame,
                            female_specific: pd.DataFrame,
                            shared_edges: pd.DataFrame,
                            gwas_genes_upper: Set[str]) -> ig.Graph:
    """Build a comparison network with male-specific, female-specific, and shared edges."""

    all_edges = pd.concat([male_specific, female_specific, shared_edges], ignore_index=True)
    tg_max_strength: Dict[str, float] = {}
    for _, r in all_edges.iterrows():
        tg = r["TG"]
        w = float(r.get("weight", 0.0))
        if np.isfinite(w):
            tg_max_strength[tg] = max(tg_max_strength.get(tg, 0.0), w)

    tfs = set(all_edges["TF"])
    tgs = set(all_edges["TG"])
    vertices = sorted(tfs | tgs)

    g = ig.Graph()
    g.add_vertices(vertices)

    TF_SIZE = 100
    TG_SIZE = 35
    TF_LABEL_SIZE = 26
    TG_LABEL_SIZE = 22
    TG_LABEL_DIST = 1.3
    TG_BASE_RGB = hex_to_rgb(COLOR_TG_NODE)

    vtype, vshape, vcolor, vsize = [], [], [], []
    label_color, label_size, label_dist, label_angle = [], [], [], []
    frame_color, frame_width = [], []

    for name in vertices:
        hit = is_gwas_gene(name, gwas_genes_upper)
        if name in tfs:
            vtype.append("TF")
            vshape.append("circle")
            vsize.append(TF_SIZE)
            vcolor.append(COLOR_TF_NODE)
            label_color.append("#000000")
            label_size.append(TF_LABEL_SIZE)
            label_dist.append(0.0)
            label_angle.append(0.0)
        else:
            vtype.append("TG")
            vshape.append("diamond")
            vsize.append(TG_SIZE)
            strength = clamp01(tg_max_strength.get(name, 0.0))
            alpha_node = scaled_alpha(strength, NODE_ALPHA_MIN, NODE_ALPHA_MAX, NODE_ALPHA_GAMMA)
            vcolor.append((TG_BASE_RGB[0], TG_BASE_RGB[1], TG_BASE_RGB[2], alpha_node))
            label_color.append("#000000")
            label_size.append(TG_LABEL_SIZE)
            label_dist.append(TG_LABEL_DIST)
            label_angle.append(np.pi / 2)

        if hit:
            frame_color.append(GWAS_HIGHLIGHT_COLOR)
            frame_width.append(GWAS_FRAME_WIDTH)
        else:
            frame_color.append("#00000000")
            frame_width.append(0)

    g.vs["name"] = vertices
    g.vs["type"] = vtype
    g.vs["shape"] = vshape
    g.vs["color"] = vcolor
    g.vs["size"] = vsize
    g.vs["label"] = vertices
    g.vs["label_color"] = label_color
    g.vs["label_size"] = label_size
    g.vs["label_dist"] = label_dist
    g.vs["label_angle"] = label_angle
    g.vs["frame_color"] = frame_color
    g.vs["frame_width"] = frame_width

    idx = {name: i for i, name in enumerate(vertices)}

    def add_edges_from(df, base_rgb, sex_label):
        elist, ecol, ewid, esex = [], [], [], []
        for _, r in df.iterrows():
            tf, tg = r["TF"], r["TG"]
            if tf not in idx or tg not in idx:
                continue
            elist.append((idx[tf], idx[tg]))
            w = clamp01(float(r.get("weight", 0.5)))
            alpha = scaled_alpha(w, EDGE_ALPHA_MIN, EDGE_ALPHA_MAX, EDGE_ALPHA_GAMMA)
            ecol.append((base_rgb[0], base_rgb[1], base_rgb[2], alpha))
            ewid.append(scaled_width(w, TF_TG_WIDTH_MIN, TF_TG_WIDTH_MAX, EDGE_ALPHA_GAMMA))
            esex.append(sex_label)
        return elist, ecol, ewid, esex

    m_e, m_c, m_w, m_s = add_edges_from(male_specific, EDGE_BASE_RGB_MALE, "Male-specific")
    f_e, f_c, f_w, f_s = add_edges_from(female_specific, EDGE_BASE_RGB_FEMALE, "Female-specific")
    s_e, s_c, s_w, s_s = add_edges_from(shared_edges, EDGE_BASE_RGB_SHARED, "Shared")

    g.add_edges(m_e + f_e + s_e)
    g.es["color"] = m_c + f_c + s_c
    g.es["width"] = m_w + f_w + s_w
    g.es["sex"] = m_s + f_s + s_s

    return g


def plot_igraph_network(g: ig.Graph, outfile: Union[str, Path]) -> Path:
    """Plot an igraph network to PDF and PNG."""
    layout = g.layout_fruchterman_reingold()

    visual_style = {
        "layout": layout,
        "vertex_color": g.vs["color"],
        "vertex_size": g.vs["size"],
        "vertex_shape": g.vs["shape"],
        "vertex_label": g.vs["name"],
        "vertex_label_size": g.vs["label_size"],
        "vertex_label_color": g.vs["label_color"],
        "vertex_label_dist": g.vs["label_dist"],
        "vertex_label_angle": g.vs["label_angle"],
        "vertex_frame_color": g.vs["frame_color"],
        "vertex_frame_width": g.vs["frame_width"],
        "edge_color": g.es["color"],
        "edge_width": g.es["width"],
        "bbox": (2600, 1900),
        "margin": 100,
    }

    out = Path(outfile)
    ig.plot(g, target=str(out.with_suffix(".png")), **visual_style)
    ig.plot(g, target=str(out.with_suffix(".pdf")), **visual_style)
    print(f"[igraph] saved to {out.with_suffix('.png')} and .pdf")
    return out.with_suffix(".png")


def make_figure_with_legend(png_path: Path, title: str, legend_items: list,
                            outfile_base: Union[str, Path]) -> Path:
    """Overlay a matplotlib legend on an igraph-rendered PNG."""
    img = plt.imread(str(png_path))

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(title, fontsize=16)

    if legend_items:
        fig.subplots_adjust(bottom=0.10)
        fig.legend(
            handles=legend_items,
            loc="lower center",
            ncol=min(3, len(legend_items)),
            frameon=True, fontsize=10,
            bbox_to_anchor=(0.5, 0.02),
        )

    fig.tight_layout(rect=(0.02, 0.08, 0.98, 0.97))

    outfile_base = Path(outfile_base)
    fig.savefig(outfile_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outfile_base.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[legend] saved to {outfile_base.with_suffix('.png')} and .pdf")
    return outfile_base.with_suffix(".png")


# --- Main --------------------------------------------------------------------


def main():
    print("=" * 60)
    print("Sex-stratified pySCENIC TF-TG igraph networks — AD")
    print("=" * 60)

    gwas_genes = load_gwas_genes(PATH_GWAS, GWAS_GENE_COL)

    # Load adjacencies
    print("\nLoading male adjacency...")
    male_adj = load_pyscenic_adj(PATH_MALE_ADJ)
    print(f"  Total edges: {len(male_adj)}, TFs: {male_adj['TF'].nunique()}")

    print("Loading female adjacency...")
    female_adj = load_pyscenic_adj(PATH_FEMALE_ADJ)
    print(f"  Total edges: {len(female_adj)}, TFs: {female_adj['TF'].nunique()}")

    # Top edges for each sex
    male_top = male_adj.head(TOP_EDGES_PER_SEX).copy()
    female_top = female_adj.head(TOP_EDGES_PER_SEX).copy()

    # --- Male heatmap + network ---
    print("\n--- MALE ---")
    male_heatmap_edges = draw_heatmap(
        male_top,
        title=f"TF x TG heatmap — Male AD (top {TOP_EDGES_PER_SEX} pySCENIC edges)",
        outpath=OUT_DIR / "heatmap_tf_x_tg_Male.pdf",
        gwas_genes_upper=gwas_genes,
    )

    print("Building Male igraph network...")
    g_male = build_igraph(male_heatmap_edges, gwas_genes, EDGE_BASE_RGB)
    raw_male = plot_igraph_network(g_male, OUT_DIR / "igraph_network_Male_raw")

    legend_male = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_TF_NODE,
               markeredgecolor="k", markersize=12, label="TF (yellow circle)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_TG_NODE,
               markeredgecolor="k", markersize=8,
               label="TG (diamond; alpha ~ importance)"),
        Line2D([0], [0], color=EDGE_BASE_RGB, lw=3,
               label="TF->TG activator (width & alpha ~ importance)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor=GWAS_HIGHLIGHT_COLOR, markeredgewidth=3,
               markersize=12, label="GWAS hit (cyan outline)"),
    ]
    make_figure_with_legend(
        raw_male,
        f"TF-TG Network — Male AD (pySCENIC top {TOP_EDGES_PER_SEX})",
        legend_male,
        OUT_DIR / "igraph_network_Male_with_legend",
    )

    # --- Female heatmap + network ---
    print("\n--- FEMALE ---")
    female_heatmap_edges = draw_heatmap(
        female_top,
        title=f"TF x TG heatmap — Female AD (top {TOP_EDGES_PER_SEX} pySCENIC edges)",
        outpath=OUT_DIR / "heatmap_tf_x_tg_Female.pdf",
        gwas_genes_upper=gwas_genes,
    )

    print("Building Female igraph network...")
    g_female = build_igraph(female_heatmap_edges, gwas_genes, EDGE_BASE_RGB)
    raw_female = plot_igraph_network(g_female, OUT_DIR / "igraph_network_Female_raw")

    legend_female = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_TF_NODE,
               markeredgecolor="k", markersize=12, label="TF (yellow circle)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_TG_NODE,
               markeredgecolor="k", markersize=8,
               label="TG (diamond; alpha ~ importance)"),
        Line2D([0], [0], color=EDGE_BASE_RGB, lw=3,
               label="TF->TG activator (width & alpha ~ importance)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor=GWAS_HIGHLIGHT_COLOR, markeredgewidth=3,
               markersize=12, label="GWAS hit (cyan outline)"),
    ]
    make_figure_with_legend(
        raw_female,
        f"TF-TG Network — Female AD (pySCENIC top {TOP_EDGES_PER_SEX})",
        legend_female,
        OUT_DIR / "igraph_network_Female_with_legend",
    )

    # --- Comparison network ---
    print("\n--- COMPARISON: Male-specific / Female-specific / Shared ---")

    TOP_COMPARISON = 500
    male_comp = male_adj.head(TOP_COMPARISON).copy()
    female_comp = female_adj.head(TOP_COMPARISON).copy()

    male_pairs = set(zip(male_comp["TF"], male_comp["TG"]))
    female_pairs = set(zip(female_comp["TF"], female_comp["TG"]))

    shared_pairs = male_pairs & female_pairs
    male_only_pairs = male_pairs - female_pairs
    female_only_pairs = female_pairs - male_pairs

    print(f"  Top {TOP_COMPARISON} per sex:")
    print(f"    Shared:          {len(shared_pairs)}")
    print(f"    Male-specific:   {len(male_only_pairs)}")
    print(f"    Female-specific: {len(female_only_pairs)}")

    # Build subsets
    male_specific_df = male_comp[
        male_comp.apply(lambda r: (r["TF"], r["TG"]) in male_only_pairs, axis=1)
    ].head(150).copy()
    female_specific_df = female_comp[
        female_comp.apply(lambda r: (r["TF"], r["TG"]) in female_only_pairs, axis=1)
    ].head(150).copy()
    shared_df = male_comp[
        male_comp.apply(lambda r: (r["TF"], r["TG"]) in shared_pairs, axis=1)
    ].head(100).copy()

    g_comp = build_igraph_comparison(
        male_specific_df, female_specific_df, shared_df, gwas_genes
    )
    raw_comp = plot_igraph_network(g_comp, OUT_DIR / "igraph_network_comparison_raw")

    legend_comp = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_TF_NODE,
               markeredgecolor="k", markersize=12, label="TF (yellow circle)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_TG_NODE,
               markeredgecolor="k", markersize=8, label="TG (diamond)"),
        Line2D([0], [0], color=COLOR_EDGE_MALE, lw=3,
               label="Male-specific TF->TG"),
        Line2D([0], [0], color=COLOR_EDGE_FEMALE, lw=3,
               label="Female-specific TF->TG"),
        Line2D([0], [0], color=COLOR_EDGE_SHARED, lw=3,
               label="Shared TF->TG"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor=GWAS_HIGHLIGHT_COLOR, markeredgewidth=3,
               markersize=12, label="GWAS hit (cyan outline)"),
    ]
    make_figure_with_legend(
        raw_comp,
        f"TF-TG Network — Male vs Female AD (top {TOP_COMPARISON} edges)",
        legend_comp,
        OUT_DIR / "igraph_network_comparison_with_legend",
    )

    # --- Save results --------------------------------------------------------

    male_specific_df.to_csv(OUT_DIR / "male_specific_top_edges.csv", index=False)
    female_specific_df.to_csv(OUT_DIR / "female_specific_top_edges.csv", index=False)
    shared_df.to_csv(OUT_DIR / "shared_top_edges.csv", index=False)

    print("\n" + "=" * 60)
    print("Done. All outputs in:", OUT_DIR)
    print("=" * 60)


if __name__ == "__main__":
    main()
