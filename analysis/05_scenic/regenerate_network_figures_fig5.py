#!/usr/bin/env python
"""
regenerate_editable_pdfs_v10.py -- per-TF quotas (1 <= per TF <= 2)
============================================================================

Goal: keep the v3-style full coverage (top 300 pairs, all heatmap TGs
as columns, all TGs as diamond nodes in the network), but annotate FAR
FEWER genes than the previous v10 pass. Per-TF quotas tightened to
1-2 labels per TF: every TF gets at least 1 named target (and 2 when
data and shared-TG constraints permit), capped at 2. Result: a very
sparse label set where AD-GWAS TGs and the highest-priority targets per
TF dominate.

Quotas:
  * Phase 0: every AD-GWAS TG that lands in the heatmap is pinned --
    these bypass MAX so no AD-relevant target is ever dropped because a
    sibling TF was already saturated.
  * MIN_LABELS_PER_TF = 1  (Phase 1 ensures each TF reaches this floor;
    GWAS pinning often satisfies this for free).
  * MAX_LABELS_PER_TF = 2  (Phase 2 may add a 2nd non-GWAS label per TF
    if the score allows and no quota would be violated; TFs that already
    exceed MAX from Phase 0 GWAS pins simply receive nothing further).

Score (GWAS hit [AD/FTD], pySCENIC global degree, connectivity to top
heatmap TFs, in-degree, total weight) ranks candidates inside the
quota -- AD-GWAS hits dominate via the 1000x bonus, so they are
labeled whenever they appear under any heatmap TF.

v10 (current): per-TF quota 1-2 (this pass).
v10 (previous): per-TF quota 2-5.
v9 lineage: per-TF quotas, MIN/MAX-driven, supplemental sheets from v5.

Based on v8. v9 changes:
  * New per-TF quota selection in select_label_indices (Phase 1: ensure
    each TF hits MIN_LABELS_PER_TF; Phase 2: fill more by score under
    MAX_LABELS_PER_TF cap).
  * Soft global ceiling N_TGS_LABELED bumped to 250 -- per-TF cap is
    the real driver.
  * Supplemental sheets unchanged from v5.

Based on v6. v8 changes:
  * TOP_TF_TG_INTERSECTION 100 -> 300 (full v3 coverage, "show all").
  * N_TGS_LABELED 200 -> 50 (selective for legibility).
  * Wider canvas (28 in) and 10pt Arial labels with 1.5-column spacing
    so the sparse labels read cleanly.
  * Network unchanged in structure: uses tf_tg_used (heatmap-filtered),
    so all heatmap TGs and all top-12 TFs become network nodes; labels
    filtered to TFs + heatmap-labeled TGs.
  * Supplemental sheets unchanged from v5.

Based on v5. v6 changes:
  * Heatmap: dense labeling (N_TGS_LABELED 40 -> 200 effective "all"), 9pt
    Arial labels (>= 8pt Nature minimum), wider canvas (22 in) so labels
    don't go sparse. Tight 1-column spacing.
  * Network: keeps ALL nodes (TFs + every TG diamond) and ALL edges from
    the focused panel. Only the NAMES are filtered: TFs are always labeled;
    TGs are labeled only if their name appears in the heatmap. Unlabeled
    TG diamonds stay visible -- they just don't have text next to them.
  * Supplemental sheets unchanged from v5.

Based on v4. v5 changes:
  * Fewer genes visualized: TOP_TF_TG_INTERSECTION 300 -> 100, N_TGS_LABELED
    130 -> 40. Heatmap shows only the most-important TGs (GWAS hits + high
    pySCENIC degree + connectivity-to-top-TFs first).
  * Network plot drops the unlabeled TG diamonds entirely. Only the focused
    (labeled) TGs appear as diamond nodes; TFs are unchanged.
  * Larger label font (10pt) since fewer labels means more breathing room.
  * For each cell type, writes a Nature-style supplementary workbook
    (Supplementary_Table_<STEM>.xlsx) with sheets covering: shared TF->TG
    intersection, TF-TF coordination, TFs/TGs in heatmap with all selection-
    score components, the focused TG subset shown in the network, and the
    full NetREm + pySCENIC TF->TG tables for that cell type. Also mirrored
    to per-sheet CSVs.

v4 lineage: full heatmap + selective connector labels. Differences from v1:
  * Colormap is "plasma" (perceptually uniform, colorblind-friendly) --
    replaces viridis. Matches the hue family used in v2's custom anchor
    palette, which was a 5-point resampling of plasma.
  * Heatmap still shows ALL shared TGs as columns -- no columns dropped.
  * X-axis labels are drawn only for a SELECTED subset of TGs, placed in
    a row below the heatmap and connected to their column with an elbow
    (L-shaped) connector. Font size is kept >= 8 so labels stay legible.
  * TG selection score prioritizes, in order:
      1. GWAS hits           (always labeled)
      2. Global pySCENIC in-degree of the TG
      3. Connectivity to the top heatmap TFs
      4. In-degree inside the heatmap
      5. Total edge weight  (tiebreak)
  * Normalization is unchanged from v1 (PowerNorm with quantile clip),
    preserving the per-plot contrast "pop".
  * Networks are unchanged from v1.

Output directory: scenic/scenic_v3_connector_labels/celltype_outputs/
"""

import os
from pathlib import Path
from typing import Union, Tuple, Set, Dict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import igraph as ig
from matplotlib.lines import Line2D
from matplotlib.colors import PowerNorm

# =====================================================================
# EDITABLE TEXT IN PDF / ADOBE ILLUSTRATOR
# =====================================================================
matplotlib.rcParams['pdf.fonttype']  = 42
matplotlib.rcParams['ps.fonttype']   = 42
matplotlib.rcParams['svg.fonttype']  = 'none'
matplotlib.rcParams['font.family']   = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']

# =====================================================================
# PATHS
# =====================================================================
ROOT     = Path(__file__).parent.parent.parent  # repo root
BASE_DIR = ROOT / "data" / "scenic" / "netrem_tables"   # NetREm+pySCENIC TF-TG CSVs (external input, see README)
OUT_DIR  = ROOT / "figures" / "scenic_v11_nature" / "celltype_outputs"
SUPP_DIR = ROOT / "figures" / "scenic_v11_nature" / "supplemental_sheets"
SUPP_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

PATH_GWAS     = ROOT / "data" / "scenic" / "gwas_genes_ad_ftd.xlsx"
GWAS_GENE_COL = "GeneName"

CELL_TYPES = [
    ("AST", "Ast", "TF_TF_BRAIN_AST.csv",
     "NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_Ast_Nov_1_2025_bySaniya_Final.csv",
     "ast_adj_adata_dementia_filtered_for_Saniya.csv",
     None, 0),
    ("EXC", "ExN", "TF_TF_BRAIN_EXN.csv",
     "NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_ExN_Nov_1_2025_bySaniya_Final.csv",
     "exc_adj_adata_dementia_filtered_for_Saniya.csv",
     None, 0),
    ("INH", "InN", "TF_TF_BRAIN_INH.csv",
     "NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_InN_Nov_1_2025_bySaniya_Final.csv",
     "inh_adj_adata_dementia_filtered_for_Saniya.csv",
     None, 0),
    ("MIC", "Mic", "TF_TF_BRAIN_MIC.csv",
     "NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_Mic_Nov_1_2025_bySaniya_Final.csv",
     "mic_adj_adata_dementia2_filtered_for_Saniya.csv",
     "SPI1", 30),
    ("ODC", "Oli", "TF_TF_BRAIN_OLI.csv",
     "NetREm_TF_to_TG_regulatory_network_links_UCI_b1_a0_Dementia_pyscenicPriorGRN_Oli_Nov_1_2025_bySaniya_Final.csv",
     "odc_adj_adata_dementia_filtered_for_Saniya.csv",
     None, 0),
]

# =====================================================================
# GLOBAL CONFIG
# =====================================================================
TOP_TF_TG_INTERSECTION = 300        # full shared-pair coverage
TOP_TF_PER_TG          = 3
TOP_TFS_HEATMAP        = 12
MIN_WEIGHT_HEATMAP     = 0.0
TOP_TF_TG_PER_METHOD   = 5000
FIGSIZE_HEATMAP        = (10, 6)    # v11: Nature-friendly width, reduced height

# v11: cap displayed TG columns to top-N by selection score (GWAS first).
# All displayed columns get a direct 45-degree x-axis tick label -- no
# elbow connectors. Keeps the heatmap readable at journal column width.
HEATMAP_DISPLAY_COLS        = 15

N_TGS_LABELED               = 15    # matches HEATMAP_DISPLAY_COLS
MIN_LABELS_PER_TF           = 1
MAX_LABELS_PER_TF           = 2
N_TOP_TFS_FOR_CONNECTIVITY  = 3
LABEL_FONTSIZE              = 8     # v11: Nature minimum ≥7pt
LABEL_SPACE_ROWS            = 0.0   # no connector-label row; direct xtick labels used
# Bracket connectors: when >= BRACKET_MIN_SIZE selected columns sit on
# directly-consecutive heatmap indices, group them under a single horizontal
# bracket. At high label density the bracket bar runs across most of the
# heatmap bottom and reads as "gene names chained by a line", so we disable
# it by default (effectively-infinite threshold). Lower the threshold to
# re-enable for "hot block" callouts.
BRACKET_MIN_SIZE            = 9999  # disabled: use pure elbow connectors
BRACKET_MAX_GAP             = 1     # max index gap between consecutive members
BRACKET_LW                  = 0.9   # slightly heavier than elbow
CONNECTOR_COLOR             = "#555555"
CONNECTOR_LW                = 0.7
# Minimum horizontal spacing between adjacent labels, in data units (column
# widths). Labels are anchored at their actual column index and only pushed
# sideways if the neighbour is closer than this. Keep slightly above the
# rotated-label footprint so labels don't touch.
MIN_LABEL_SPACING           = 1.5   # v8: spread labels (fewer of them, more room)
# Selection-time spacing floor (data-unit columns). Selection picks highest-
# score TG first and only accepts subsequent candidates whose column is
# >= LABEL_TARGET_SPACING_MIN away from every already-picked column. 1.0
# means adjacent columns are allowed (brackets will group them); higher
# values force near-vertical connectors at the cost of coverage.
LABEL_TARGET_SPACING_MIN    = 1.0
# If the horizontal offset between a label and its column is smaller than
# this (data units), draw a direct vertical tick instead of an elbow.
DIRECT_ANCHOR_TOL           = 0.4

# =====================================================================
# COLOUR CONSTANTS
# =====================================================================
COLOR_TG_UP      = "#ff4400"
COLOR_TG_DOWN    = "#6d23b6"
COLOR_TG_NEUTRAL = "#48bfe3"
COLOR_TF_NODE    = "#ffe666"

COLOR_TF_TF_COOP  = "#2ecc71"
COLOR_TF_TF_ANTAG = "#6f3ce7"

EDGE_BASE_ACT_RGB = (1.00, 0.35, 0.05)
EDGE_BASE_REP_RGB = (0.60, 0.35, 0.80)
EDGE_BASE_NET_RGB = (0.05, 0.80, 0.95)

GWAS_HIGHLIGHT_COLOR = "#00e5ff"
GWAS_FRAME_WIDTH     = 4

# Heatmap contrast (v1 behavior -- PowerNorm with quantile clip)
HEATMAP_CMAP          = "plasma"          # v3: plasma, colorblind-friendly
HEATMAP_QUANTILE_CLIP = (0.05, 0.95)
HEATMAP_GAMMA         = 0.65

# Network contrast (unchanged from v1)
EDGE_ALPHA_MIN   = 0.10
EDGE_ALPHA_MAX   = 1.0
NODE_ALPHA_MIN   = 0.25
NODE_ALPHA_MAX   = 1.0
EDGE_ALPHA_GAMMA = 1.5
NODE_ALPHA_GAMMA = 0.6

TF_TG_WIDTH_MIN  = 1.0
TF_TG_WIDTH_MAX  = 14.0
TF_TF_WIDTH_MIN  = 5.0
TF_TF_WIDTH_MAX  = 14.0

# =====================================================================
# CELL TYPE FILTERING
# =====================================================================
CELLTYPE_KEYWORDS_MAP = {
    "AST": ["ast", "astro", "astrocyte"],
    "OL":  ["ol", "oli", "oligo", "oligodendrocyte"],
    "OLI": ["ol", "oli", "oligo", "oligodendrocyte"],
    "ODC": ["ol", "oli", "oligo", "oligodendrocyte", "odc"],
    "MIC": ["mic", "micro", "microglia", "microglial"],
    "EX":  ["ex", "excitatory"],
    "EXC": ["ex", "excitatory", "exc"],
    "INH": ["inh", "inhibitory"],
}
CELLTYPE_COLS = [
    "Final_Cell_type", "cell_type", "CellType", "celltype",
    "info_type", "cluster", "group", "celltype_label"
]

def _get_keywords(stem):
    return CELLTYPE_KEYWORDS_MAP.get(stem.upper(), [stem.lower()])

def filter_celltype(df, stem):
    keywords = _get_keywords(stem)
    pat = "|".join(keywords)
    for c in CELLTYPE_COLS:
        if c in df.columns:
            mask = df[c].astype(str).str.lower().str.contains(pat, na=False)
            return df.loc[mask].copy()
        for col in df.columns:
            if col.lower() == c.lower():
                mask = df[col].astype(str).str.lower().str.contains(pat, na=False)
                return df.loc[mask].copy()
    return df.copy()

# =====================================================================
# HELPERS
# =====================================================================
def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def normalize_0_1(x):
    x = pd.to_numeric(x, errors="coerce")
    if x.dropna().empty:
        return pd.Series(np.zeros(len(x)), index=x.index)
    mn, mx = x.min(), x.max()
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - mn) / (mx - mn)

def clamp01(x):
    return max(0.0, min(1.0, float(x)))

def scaled_alpha(w01, amin, amax, gamma):
    w01 = clamp01(w01)
    return amin + (amax - amin) * (w01 ** gamma)

def scaled_width(w01, wmin, wmax, gamma):
    w01 = clamp01(w01)
    return wmin + (wmax - wmin) * (w01 ** gamma)

def top_k_per_group(df, by_col, group_col, k):
    if k is None:
        return df
    return (df.sort_values(by_col, ascending=False)
              .groupby(group_col, as_index=False, group_keys=False)
              .head(k))

def edges_to_matrix(edges, row_key, col_key, weight_key, row_order=None, col_order=None):
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

def choose_score_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce")
            if not s.dropna().empty:
                return s
    return pd.Series(np.nan, index=df.index)

def load_gwas_genes(path, gene_col):
    if not Path(path).exists():
        print(f"[GWAS] file not found: {path}")
        return set()
    df = pd.read_excel(path)
    if gene_col not in df.columns:
        print(f"[GWAS] column '{gene_col}' not found")
        return set()
    genes = df[gene_col].dropna().astype(str).str.strip().str.upper().tolist()
    genes_set = set(genes)
    print(f"[GWAS] loaded {len(genes_set)} genes")
    return genes_set

def is_gwas(gene, gwas_set):
    return str(gene).strip().upper() in gwas_set

# =====================================================================
# LOADERS
# =====================================================================
def load_tf_tf_netrem(path, stem):
    df = pd.read_csv(path)
    df = filter_celltype(df, stem)

    req = ["TF1", "TF2", "overall_cell_type_coordination_relationship",
           "Percentile_cooperative", "Percentile_antagonistic"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"TF-TF file missing columns: {missing}")

    tf1 = df["TF1"].astype(str).str.strip()
    tf2 = df["TF2"].astype(str).str.strip()

    rel_txt = df["overall_cell_type_coordination_relationship"].astype(str).str.strip().str.lower()
    relation = np.where(rel_txt.str.contains("antagon", na=False), "antagonistic", "cooperative")

    pc = pd.to_numeric(df["Percentile_cooperative"], errors="coerce")
    pa = pd.to_numeric(df["Percentile_antagonistic"], errors="coerce")
    w_raw = np.where(relation == "antagonistic", pa.values, pc.values)
    w = pd.Series(w_raw, index=df.index)
    w = (w / 100.0).clip(lower=0.0, upper=1.0).fillna(0.0)

    out = pd.DataFrame({"TF1": tf1, "TF2": tf2, "relation": relation, "weight": w})
    return out.dropna(subset=["TF1", "TF2"])

def load_tf_tg_netrem(path):
    df = pd.read_csv(path)
    for col in ("coef", "TF", "TG"):
        if col not in df.columns:
            raise ValueError(f"NetRem TF-TG missing column '{col}'")

    df["TF"] = df["TF"].astype(str).str.strip()
    df["TG"] = df["TG"].astype(str).str.strip()
    coef = pd.to_numeric(df["coef"], errors="coerce")

    out = pd.DataFrame({
        "TF": df["TF"], "TG": df["TG"],
        "coef": coef, "rank_score_netrem": coef.abs(),
    })
    out["weight_netrem"] = normalize_0_1(out["rank_score_netrem"])
    sign = np.sign(out["coef"])
    out["role_netrem"] = np.where(sign > 0, "activator",
                          np.where(sign < 0, "repressor", "neutral"))

    out = (out.groupby(["TF", "TG"], as_index=False)
              .agg(coef=("coef", "mean"),
                   rank_score_netrem=("rank_score_netrem", "max"),
                   weight_netrem=("weight_netrem", "max"),
                   role_netrem=("role_netrem", "first")))
    return out

def load_tf_tg_pyscenic(path, stem):
    df = pd.read_csv(path)
    df = filter_celltype(df, stem)
    for col in ("TF", "target", "importance"):
        if col not in df.columns:
            raise ValueError(f"pySCENIC missing column '{col}'")

    df = df.rename(columns={"target": "TG"})
    df["TF"] = df["TF"].astype(str).str.strip()
    df["TG"] = df["TG"].astype(str).str.strip()

    score_raw = choose_score_column(df, [
        "Final_Cell_TF_TG_score_beta", "Final_Cell_TF_TG_overall_score_beta",
        "Final_Cell_TF_TG_score", "Final_Cell_TF_TG_score_abs", "importance",
    ])
    if score_raw.dropna().empty:
        score_raw = pd.to_numeric(df["importance"], errors="coerce")

    out = pd.DataFrame({"TF": df["TF"], "TG": df["TG"], "rank_score_pyscenic": score_raw})
    out["weight_pyscenic"] = normalize_0_1(out["rank_score_pyscenic"])
    out = (out.groupby(["TF", "TG"], as_index=False)
              .agg(rank_score_pyscenic=("rank_score_pyscenic", "max"),
                   weight_pyscenic=("weight_pyscenic", "max")))
    return out

# =====================================================================
# SHARED TF-TG
# =====================================================================
def build_shared_tf_tg_all(net_top, scen_top):
    merged = net_top.merge(
        scen_top[["TF", "TG", "rank_score_pyscenic", "weight_pyscenic"]],
        on=["TF", "TG"], how="inner")
    if merged.empty:
        print("[shared] WARNING: no shared TF-TG pairs.")
        merged["weight_combined"] = np.nan
        return merged
    merged["weight_combined"] = merged["weight_netrem"]
    merged = merged.sort_values("rank_score_netrem", ascending=False).reset_index(drop=True)
    return merged

def force_shared_pairs_for_tf(net_all, scen_all, tf_name):
    tf_name = str(tf_name).strip()
    net_tf = net_all[net_all["TF"].astype(str).str.strip() == tf_name].copy()
    scen_tf = scen_all[scen_all["TF"].astype(str).str.strip() == tf_name].copy()
    if net_tf.empty or scen_tf.empty:
        return pd.DataFrame()
    merged = net_tf.merge(
        scen_tf[["TF", "TG", "rank_score_pyscenic", "weight_pyscenic"]],
        on=["TF", "TG"], how="inner")
    if merged.empty:
        return pd.DataFrame()
    merged["weight_combined"] = merged["weight_netrem"]
    return merged.sort_values("rank_score_netrem", ascending=False).reset_index(drop=True)

def get_tg_status_from_role(tf_tg_df):
    status = {}
    if "role_netrem" not in tf_tg_df.columns:
        for g in tf_tg_df["TG"].unique():
            status[g] = "neutral"
        return status
    tmp = tf_tg_df.copy()
    if "weight_combined" in tmp.columns:
        tmp = tmp.sort_values("weight_combined", ascending=False)
    for tg, sub in tmp.groupby("TG"):
        role = str(sub.iloc[0]["role_netrem"]).lower()
        if "activator" in role:
            status[tg] = "up"
        elif "repressor" in role:
            status[tg] = "down"
        else:
            status[tg] = "neutral"
    return status

# =====================================================================
# LABEL POSITION ANCHORING
# =====================================================================
def _anchor_label_positions(targets, n_cols, min_spacing):
    """Return label x-positions that stay as close as possible to each
    target column index while respecting a minimum horizontal spacing
    and the [0, n_cols-1] extent of the heatmap.

    Two-pass algorithm (left-to-right then right-to-left) with a final
    left clamp pass -- standard adjusted-text placement.
    """
    xs = list(targets)
    n = len(xs)
    if n == 0:
        return []

    # Left-to-right: push each label right if it crowds its left neighbour
    for i in range(1, n):
        if xs[i] < xs[i - 1] + min_spacing:
            xs[i] = xs[i - 1] + min_spacing

    # Clamp right edge and cascade leftward if we overshot
    right = float(n_cols - 1)
    if xs[-1] > right:
        xs[-1] = right
        for i in range(n - 2, -1, -1):
            if xs[i] > xs[i + 1] - min_spacing:
                xs[i] = xs[i + 1] - min_spacing
            else:
                break

    # Clamp left edge and cascade rightward if we undershot
    if xs[0] < 0.0:
        xs[0] = 0.0
        for i in range(1, n):
            if xs[i] < xs[i - 1] + min_spacing:
                xs[i] = xs[i - 1] + min_spacing
            else:
                break

    return xs

# =====================================================================
# TG SELECTION FOR LABELING
# =====================================================================
def select_label_indices(M, cols_tg, gwas_set, scen_all_for_degree):
    """Rank TG columns for labeling. Returns sorted-by-column-index list.

    Priority (highest to lowest):
      1. GWAS hits            -- always labeled (large additive bonus)
      2. Global pySCENIC in-degree of the TG (normalized)
      3. Connectivity to the top-N heatmap TFs
      4. Heatmap in-degree
      5. Total edge weight    (tiebreak)
    """
    n_rows, n_cols = M.shape
    if n_cols == 0:
        return []

    # Heatmap in-degree
    in_degree = (M > 0).sum(axis=0).astype(float)

    # Total weight (tiebreaker)
    total_w = M.sum(axis=0)

    # GWAS hits
    is_gwas_col = np.array([is_gwas(g, gwas_set) for g in cols_tg], dtype=float)

    # Connectivity to top heatmap TFs (by row sum)
    tf_strength = M.sum(axis=1)
    n_top_tf = min(N_TOP_TFS_FOR_CONNECTIVITY, n_rows)
    if n_top_tf > 0:
        top_tf_idx = np.argsort(-tf_strength)[:n_top_tf]
        connected_to_top_tf = (M[top_tf_idx, :] > 0).any(axis=0).astype(float)
    else:
        connected_to_top_tf = np.zeros(n_cols)

    # Global pySCENIC degree of each TG (across full pySCENIC table)
    scenic_degree = np.zeros(n_cols)
    if scen_all_for_degree is not None and "TG" in scen_all_for_degree.columns:
        deg_map = scen_all_for_degree.groupby("TG").size().to_dict()
        for j, tg in enumerate(cols_tg):
            scenic_degree[j] = deg_map.get(tg, 0)
        if scenic_degree.max() > 0:
            scenic_degree = scenic_degree / scenic_degree.max()

    score = (
        1000.0 * is_gwas_col
        + 5.0  * scenic_degree
        + 3.0  * connected_to_top_tf
        + 1.0  * (in_degree / max(in_degree.max(), 1.0))
        + 0.001 * total_w
    )

    n_labels_target = min(N_TGS_LABELED, n_cols)
    if n_labels_target == 0:
        return []

    # v10: per-TF quota-based selection with GWAS pin.
    # Each TF row in the heatmap should have AT LEAST MIN_LABELS_PER_TF
    # of its connected TGs labeled, and AT MOST MAX_LABELS_PER_TF.
    # A TG can be "attached" to multiple TFs (M[i,j] > 0 for several i),
    # so we track per-TF counts and reject any TG whose addition would
    # push any of its connected TFs over the max.
    #
    # Phase 0: pin every AD-GWAS TG present in the heatmap (these are the
    #   "important TGs in human AD"; they bypass MAX_LABELS_PER_TF so no
    #   GWAS hit is dropped due to a sibling TF being saturated).
    # Phase 1: walk TFs in priority order (weakest first); for each TF,
    #   greedily pick its highest-scoring un-selected TGs until either it
    #   reaches MIN_LABELS_PER_TF or it runs out of admissible candidates.
    # Phase 2: walk all TGs in score order and add any that still fit
    #   under the per-TF max (this fills the heatmap up to the max
    #   without crossing it).
    tg_to_tfs = [list(np.where(M[:, j] > 0)[0]) for j in range(n_cols)]
    tf_counts = np.zeros(n_rows, dtype=int)
    selected_set = set()

    def _try_add(j_int):
        """Add TG j if no connected TF would exceed MAX_LABELS_PER_TF."""
        for tf_i in tg_to_tfs[j_int]:
            if tf_counts[tf_i] >= MAX_LABELS_PER_TF:
                return False
        selected_set.add(j_int)
        for tf_i in tg_to_tfs[j_int]:
            tf_counts[tf_i] += 1
        return True

    # Phase 0: always label every AD-GWAS hit (cap-bypassing).
    for j in range(n_cols):
        if is_gwas_col[j] > 0:
            selected_set.add(j)
            for tf_i in tg_to_tfs[j]:
                tf_counts[tf_i] += 1

    # Phase 1: ensure each TF has MIN_LABELS_PER_TF labels (where possible).
    # Process WEAKEST TFs first -- if we did strongest-first, strong TFs hit
    # MAX before we ever try the weak ones, blocking _try_add for shared TGs
    # and starving the weak TFs. Doing weakest-first means strong TFs are
    # still well below MAX when we reach them, so they happily accept
    # whatever shared TGs we picked for weak TFs.
    tf_priority = list(np.argsort(tf_strength))   # weakest -> strongest
    for tf_i in tf_priority:
        connected_idx = np.where(M[tf_i, :] > 0)[0]
        if len(connected_idx) == 0:
            continue
        ranked = sorted(connected_idx.tolist(), key=lambda j: -score[j])
        for j in ranked:
            if tf_counts[tf_i] >= MIN_LABELS_PER_TF:
                break
            if j in selected_set:
                continue
            _try_add(int(j))   # may skip if a sibling TF is already at MAX

    # Phase 2: fill remaining capacity by global score (still respecting MAX)
    order = np.argsort(-score)
    for j in order:
        j_int = int(j)
        if j_int in selected_set:
            continue
        if len(selected_set) >= n_labels_target:
            break
        _try_add(j_int)

    return sorted(selected_set)

# =====================================================================
# HEATMAP -- v3: full matrix + selective connector labels
# =====================================================================
def draw_heatmap(tf_tg, title, outpath, gwas_set,
                 force_tf_row=None, scen_all_for_degree=None):
    tg_status = get_tg_status_from_role(tf_tg)
    edges = tf_tg.dropna(subset=["weight_combined"]).copy()
    if MIN_WEIGHT_HEATMAP > 0:
        edges = edges[edges["weight_combined"] >= MIN_WEIGHT_HEATMAP]
    if edges.empty:
        print("[heatmap] no edges -> skipping.")
        return edges, set()

    edges = top_k_per_group(edges, "weight_combined", "TG", TOP_TF_PER_TG)
    if edges.empty:
        return edges, set()

    tf_strength = edges.groupby("TF")["weight_combined"].sum().sort_values(ascending=False)
    keep_tfs = tf_strength.head(TOP_TFS_HEATMAP).index.tolist()

    if force_tf_row:
        ftf = str(force_tf_row).strip()
        if ftf and (ftf in edges["TF"].unique()) and (ftf not in keep_tfs):
            keep_tfs.append(ftf)

    edges = edges[edges["TF"].isin(keep_tfs)]
    if edges.empty:
        return edges, set()

    M, rows_tf, cols_tg = edges_to_matrix(
        edges, "TF", "TG", "weight_combined", row_order=keep_tfs)

    # Waterfall ordering of TG columns by dominant TF row
    primary_tf_idx = M.argmax(axis=0)
    max_vals = M.max(axis=0)
    col_info = [(j, int(primary_tf_idx[j]), float(max_vals[j])) for j in range(M.shape[1])]
    col_info_sorted = sorted(col_info, key=lambda x: (x[1], -x[2]))
    col_order = [j for (j, _, _) in col_info_sorted]
    M = M[:, col_order]
    cols_tg = [cols_tg[j] for j in col_order]

    n_rows, n_cols = M.shape

    # v11: cap displayed columns to HEATMAP_DISPLAY_COLS best TGs (GWAS first,
    # then by composite selection score). After waterfall ordering, we score
    # every column and keep the top-N while preserving the waterfall sequence
    # so the block structure (TG groups under each dominant TF) is unchanged.
    if HEATMAP_DISPLAY_COLS > 0 and n_cols > HEATMAP_DISPLAY_COLS:
        # Compute column scores exactly as select_label_indices does
        in_degree   = (M > 0).sum(axis=0).astype(float)
        total_w     = M.sum(axis=0)
        is_gwas_col = np.array([is_gwas(g, gwas_set) for g in cols_tg], dtype=float)
        tf_strength_rows = M.sum(axis=1)
        n_top = min(N_TOP_TFS_FOR_CONNECTIVITY, n_rows)
        top_tf_idx = np.argsort(-tf_strength_rows)[:n_top] if n_top > 0 else []
        connected_to_top = (M[top_tf_idx, :] > 0).any(axis=0).astype(float) if n_top > 0 \
                           else np.zeros(n_cols)
        scenic_degree = np.zeros(n_cols)
        if scen_all_for_degree is not None and "TG" in scen_all_for_degree.columns:
            deg_map = scen_all_for_degree.groupby("TG").size().to_dict()
            for j, tg in enumerate(cols_tg):
                scenic_degree[j] = deg_map.get(tg, 0)
            if scenic_degree.max() > 0:
                scenic_degree /= scenic_degree.max()
        col_score = (1000.0 * is_gwas_col
                     + 5.0  * scenic_degree
                     + 3.0  * connected_to_top
                     + 1.0  * (in_degree / max(in_degree.max(), 1.0))
                     + 0.001 * total_w)
        # Select top-N by score, but preserve waterfall (column) order
        top_col_idx = sorted(np.argsort(-col_score)[:HEATMAP_DISPLAY_COLS])
        M       = M[:, top_col_idx]
        cols_tg = [cols_tg[j] for j in top_col_idx]
        n_cols  = len(cols_tg)

    n_rows, n_cols = M.shape   # recompute after possible cap

    # Contrast scaling (v1 behavior)
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

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=FIGSIZE_HEATMAP)
    im = ax.imshow(M, aspect="auto", cmap=HEATMAP_CMAP, norm=norm)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"NetRem weight  (PowerNorm g={HEATMAP_GAMMA})", fontsize=8)

    # TF (y-axis) labels -- bold + GWAS frame
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(rows_tf, fontsize=8, fontfamily="Arial")
    gwas_bbox = dict(boxstyle="round,pad=0.15", facecolor="none",
                     edgecolor=GWAS_HIGHLIGHT_COLOR, linewidth=1.2)
    for lbl, gene in zip(ax.get_yticklabels(), rows_tf):
        if is_gwas(gene, gwas_set):
            lbl.set_fontweight("bold")
            lbl.set_bbox(gwas_bbox)

    # v11: direct 45° x-axis tick labels -- no elbow connectors.
    # All HEATMAP_DISPLAY_COLS columns are labeled. GWAS TGs bold + colored.
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(cols_tg, rotation=45, ha="right",
                       fontsize=LABEL_FONTSIZE, fontfamily="Arial")
    statuses_present = set()
    for lbl, tg in zip(ax.get_xticklabels(), cols_tg):
        status = tg_status.get(tg, "neutral")
        statuses_present.add(status)
        if status == "up":
            lbl.set_color(COLOR_TG_UP)
        elif status == "down":
            lbl.set_color(COLOR_TG_DOWN)
        else:
            lbl.set_color(COLOR_TG_NEUTRAL)
        if is_gwas(tg, gwas_set):
            lbl.set_fontweight("bold")

    ax.tick_params(axis="both", which="both", length=0)
    for _side in ("top", "right"):
        ax.spines[_side].set_visible(False)

    ax.set_title(title, fontsize=10, fontfamily="Arial")

    legend_handles = []
    if "up" in statuses_present:
        legend_handles.append(Line2D([0], [0], marker="s", color="none",
            markerfacecolor=COLOR_TG_UP, markersize=8,
            label="Activated (coef > 0)"))
    if "down" in statuses_present:
        legend_handles.append(Line2D([0], [0], marker="s", color="none",
            markerfacecolor=COLOR_TG_DOWN, markersize=8,
            label="Repressed (coef < 0)"))
    if "neutral" in statuses_present:
        legend_handles.append(Line2D([0], [0], marker="s", color="none",
            markerfacecolor=COLOR_TG_NEUTRAL, markersize=8,
            label="Neutral / mixed"))

    if legend_handles:
        fig.legend(handles=legend_handles, loc="lower center",
                   ncol=min(len(legend_handles), 3),
                   frameon=True, fontsize=8,
                   title="TG x-axis color")
        fig.subplots_adjust(bottom=0.18)
    else:
        fig.tight_layout()

    outpath = Path(outpath)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    labeled_tg_set = set(cols_tg)   # all displayed columns count as "labeled"
    print(f"  [heatmap] saved: {outpath}  ({n_cols} TGs displayed)")
    return edges.copy(), labeled_tg_set

# =====================================================================
# IGRAPH NETWORK (unchanged from v1)
# =====================================================================
def build_igraph(tf_tg_subset, tf_tf_subset, gwas_set, tg_label_whitelist=None):
    tg_status_map = get_tg_status_from_role(tf_tg_subset)

    tg_max_strength = {}
    for _, r in tf_tg_subset.iterrows():
        tg = r["TG"]
        w = float(r.get("weight_combined", r.get("weight_netrem", 0.0)))
        if np.isfinite(w):
            tg_max_strength[tg] = max(tg_max_strength.get(tg, 0.0), w)

    tfs = set(tf_tg_subset["TF"]) | set(tf_tf_subset["TF1"]) | set(tf_tf_subset["TF2"])
    tgs = set(tf_tg_subset["TG"])
    vertices = sorted(tfs | tgs)

    g = ig.Graph()
    g.add_vertices(vertices)

    TF_SIZE = 100; TG_SIZE = 35
    TF_LABEL_SIZE = 52; TG_LABEL_SIZE = 38   # v11: larger for Nature panel scale
    TG_LABEL_DIST = 1.3

    TG_BASE_UP_RGB   = hex_to_rgb(COLOR_TG_UP)
    TG_BASE_DOWN_RGB = hex_to_rgb(COLOR_TG_DOWN)
    TG_BASE_NEUT_RGB = hex_to_rgb(COLOR_TG_NEUTRAL)

    vtype, vstatus, vshape, vcolor, vsize = [], [], [], [], []
    label_color, label_size, label_dist, label_angle = [], [], [], []
    frame_color, frame_width = [], []

    for name in vertices:
        is_g = is_gwas(name, gwas_set)
        if name in tfs:
            vtype.append("TF"); vstatus.append("TF")
            vshape.append("circle"); vsize.append(TF_SIZE)
            vcolor.append(COLOR_TF_NODE)
            label_color.append("#000000")
            label_size.append(TF_LABEL_SIZE)
            label_dist.append(0.0); label_angle.append(0.0)
        else:
            vtype.append("TG")
            status = tg_status_map.get(name, "neutral")
            vstatus.append(status)
            vshape.append("diamond"); vsize.append(TG_SIZE)

            strength = clamp01(tg_max_strength.get(name, 0.0))
            alpha_node = scaled_alpha(strength, NODE_ALPHA_MIN, NODE_ALPHA_MAX, NODE_ALPHA_GAMMA)

            if status == "up":
                base = TG_BASE_UP_RGB
            elif status == "down":
                base = TG_BASE_DOWN_RGB
            else:
                base = TG_BASE_NEUT_RGB

            vcolor.append((base[0], base[1], base[2], alpha_node))
            label_color.append("#000000")
            label_size.append(TG_LABEL_SIZE)
            label_dist.append(TG_LABEL_DIST)
            label_angle.append(np.pi / 2)

        if is_g:
            frame_color.append(GWAS_HIGHLIGHT_COLOR)
            frame_width.append(GWAS_FRAME_WIDTH)
        else:
            frame_color.append("#00000000")
            frame_width.append(0)

    if tg_label_whitelist is None:
        labels_per_vertex = list(vertices)
    else:
        labels_per_vertex = [
            name if (vt == "TF" or name in tg_label_whitelist) else ""
            for name, vt in zip(vertices, vtype)
        ]

    # v11: bold Arial for TFs, regular Arial for TGs (Nature font requirements)
    label_font = ["Arial Bold" if vt == "TF" else "Arial" for vt in vtype]

    g.vs["name"] = vertices; g.vs["type"] = vtype; g.vs["status"] = vstatus
    g.vs["shape"] = vshape; g.vs["color"] = vcolor; g.vs["size"] = vsize
    g.vs["label"] = labels_per_vertex; g.vs["label_color"] = label_color
    g.vs["label_size"] = label_size; g.vs["label_dist"] = label_dist
    g.vs["label_angle"] = label_angle; g.vs["label_font"] = label_font
    g.vs["frame_color"] = frame_color; g.vs["frame_width"] = frame_width

    idx = {name: i for i, name in enumerate(vertices)}

    tf_tf_edges, tf_tf_color, tf_tf_width = [], [], []
    for _, r in tf_tf_subset.iterrows():
        a, b = r["TF1"], r["TF2"]
        if a not in idx or b not in idx:
            continue
        tf_tf_edges.append((idx[a], idx[b]))
        w = clamp01(float(r["weight"]))
        relation = str(r.get("relation", "cooperative")).lower()
        col = COLOR_TF_TF_ANTAG if "antag" in relation else COLOR_TF_TF_COOP
        tf_tf_color.append(col)
        tf_tf_width.append(TF_TF_WIDTH_MIN + (TF_TF_WIDTH_MAX - TF_TF_WIDTH_MIN) * w)

    tf_tg_edges, tf_tg_color, tf_tg_width, tf_tg_roles = [], [], [], []
    for _, r in tf_tg_subset.iterrows():
        tf_name, tg_name = r["TF"], r["TG"]
        if tf_name not in idx or tg_name not in idx:
            continue
        tf_tg_edges.append((idx[tf_name], idx[tg_name]))

        w = clamp01(float(r.get("weight_combined", r.get("weight_netrem", 0.5))))
        alpha = scaled_alpha(w, EDGE_ALPHA_MIN, EDGE_ALPHA_MAX, EDGE_ALPHA_GAMMA)

        role = str(r.get("role_netrem", "neutral")).lower()
        if "activator" in role:
            base = EDGE_BASE_ACT_RGB; role_label = "activator"
        elif "repressor" in role:
            base = EDGE_BASE_REP_RGB; role_label = "repressor"
        else:
            base = EDGE_BASE_NET_RGB; role_label = "neutral"

        tf_tg_roles.append(role_label)
        tf_tg_color.append((base[0], base[1], base[2], alpha))
        tf_tg_width.append(scaled_width(w, TF_TG_WIDTH_MIN, TF_TG_WIDTH_MAX, EDGE_ALPHA_GAMMA))

    g.add_edges(tf_tf_edges + tf_tg_edges)
    n_tf_tf = len(tf_tf_edges)
    g.es["kind"]  = ["TF-TF"] * n_tf_tf + ["TF-TG"] * len(tf_tg_edges)
    g.es["color"] = tf_tf_color + tf_tg_color
    g.es["width"] = tf_tf_width + tf_tg_width
    g.es["role"]  = [None] * n_tf_tf + tf_tg_roles

    return g

def plot_igraph_network(g, outfile_png):
    layout = g.layout_fruchterman_reingold()
    visual_style = {
        "layout": layout,
        "vertex_color":       g.vs["color"],
        "vertex_size":        g.vs["size"],
        "vertex_shape":       g.vs["shape"],
        "vertex_label":       g.vs["label"],
        "vertex_label_size":  g.vs["label_size"],
        "vertex_label_color": g.vs["label_color"],
        "vertex_label_dist":  g.vs["label_dist"],
        "vertex_label_angle": g.vs["label_angle"],
        "vertex_label_font":  g.vs["label_font"],   # v11: Arial / Arial Bold
        "vertex_frame_color": g.vs["frame_color"],
        "vertex_frame_width": g.vs["frame_width"],
        "edge_color": g.es["color"],
        "edge_width": g.es["width"],
        "bbox": (2600, 1900),
        "margin": 100,
    }
    out = Path(outfile_png)
    ig.plot(g, target=str(out), **visual_style)
    print(f"  [igraph] saved raw PNG: {out}")

    out_pdf = out.with_suffix(".pdf")
    ig.plot(g, target=str(out_pdf), **visual_style)
    print(f"  [igraph] saved raw PDF: {out_pdf}")
    return out

def make_legend_figure(png_path, g, outfile_base, out_stem):
    img = plt.imread(str(png_path))
    fig, ax = plt.subplots(figsize=(5.5, 4.5))   # v11: Nature panel size
    ax.imshow(img); ax.axis("off")
    ax.set_title(f"TF-TF + TF->TG network - {out_stem}", fontsize=16)

    has_tf      = any(v["type"] == "TF" for v in g.vs)
    has_tg_up   = any(v["type"] == "TG" and v["status"] == "up" for v in g.vs)
    has_tg_down = any(v["type"] == "TG" and v["status"] == "down" for v in g.vs)
    has_tg_neut = any(v["type"] == "TG" and v["status"] == "neutral" for v in g.vs)
    has_act  = any(e["kind"] == "TF-TG" and e["role"] == "activator" for e in g.es)
    has_rep  = any(e["kind"] == "TF-TG" and e["role"] == "repressor" for e in g.es)
    has_net  = any(e["kind"] == "TF-TG" and e["role"] == "neutral" for e in g.es)
    has_coop = any(e["kind"] == "TF-TF" and e["color"] == COLOR_TF_TF_COOP for e in g.es)
    has_ant  = any(e["kind"] == "TF-TF" and e["color"] == COLOR_TF_TF_ANTAG for e in g.es)
    has_gwas_nodes = ("frame_color" in g.vs.attributes()) and any(
        str(fc).lower() == str(GWAS_HIGHLIGHT_COLOR).lower() for fc in g.vs["frame_color"])

    elems = []
    if has_tf:
        elems.append(Line2D([0],[0], marker="o", color="w",
            markerfacecolor=COLOR_TF_NODE, markeredgecolor="k",
            markersize=12, label="TF (yellow circle)"))
    if has_tg_up:
        elems.append(Line2D([0],[0], marker="D", color="w",
            markerfacecolor=COLOR_TG_UP, markeredgecolor="k",
            markersize=8, label="TG up (diamond; alpha = strength)"))
    if has_tg_down:
        elems.append(Line2D([0],[0], marker="D", color="w",
            markerfacecolor=COLOR_TG_DOWN, markeredgecolor="k",
            markersize=8, label="TG down (diamond; alpha = strength)"))
    if has_tg_neut:
        elems.append(Line2D([0],[0], marker="D", color="w",
            markerfacecolor=COLOR_TG_NEUTRAL, markeredgecolor="k",
            markersize=8, label="TG neutral (diamond; alpha = strength)"))
    if has_gwas_nodes:
        elems.append(Line2D([0],[0], marker="o", color="w",
            markerfacecolor="none", markeredgecolor=GWAS_HIGHLIGHT_COLOR,
            markeredgewidth=3, markersize=12, label="GWAS hit (cyan outline)"))
    if has_act:
        elems.append(Line2D([0],[0], color=EDGE_BASE_ACT_RGB, lw=3,
            label="TF->TG activator (width & alpha = strength)"))
    if has_rep:
        elems.append(Line2D([0],[0], color=EDGE_BASE_REP_RGB, lw=3,
            label="TF->TG repressor (width & alpha = strength)"))
    if has_net:
        elems.append(Line2D([0],[0], color=EDGE_BASE_NET_RGB, lw=3,
            label="TF->TG neutral (width & alpha = strength)"))
    if has_coop:
        elems.append(Line2D([0],[0], color=COLOR_TF_TF_COOP, lw=3,
            label="TF-TF cooperative"))
    if has_ant:
        elems.append(Line2D([0],[0], color=COLOR_TF_TF_ANTAG, lw=3,
            label="TF-TF antagonistic"))

    if elems:
        fig.legend(handles=elems, loc="lower center",
                   ncol=min(3, len(elems)), frameon=True, fontsize=10,
                   bbox_to_anchor=(0.5, 0.03))

    fig.tight_layout(rect=(0.02, 0.08, 0.98, 0.97))

    base = Path(outfile_base)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [legend] saved: {base.with_suffix('.pdf')}")

# =====================================================================
# SUPPLEMENTAL SHEETS  (Nature-style multi-sheet xlsx + tidy CSVs)
# =====================================================================
SUPP_README_ROWS = [
    ("Sheet", "Description"),
    ("README", "This sheet -- column dictionary and provenance for every other sheet."),
    ("01_TF_TG_shared",
     "TF->TG edges supported by BOTH NetREm and pySCENIC (intersection). "
     "One row per TF-TG pair. Includes coefficient sign (role), per-method "
     "rank scores, normalized weights, and GWAS flag for the TG."),
    ("02_TF_TF_coordination",
     "TF-TF cooperative or antagonistic relationships from NetREm "
     "(filtered to this cell type). One row per TF1-TF2 pair."),
    ("03_TFs_in_heatmap",
     "The 12 TFs visualized as rows in the heatmap and as yellow circles "
     "in the network. Includes total NetREm-weighted out-strength and "
     "number of distinct TGs each TF connects to in the focused panel."),
    ("04_TGs_in_heatmap",
     "Every TG plotted as a column in the heatmap. Columns include the "
     "selection score components (GWAS, global pySCENIC degree, "
     "connectivity to top heatmap TFs, in-degree, total weight) used to "
     "rank candidates for labeling."),
    ("05_TGs_visualized_in_network",
     "The focused subset of TGs that appear as labeled diamond nodes in "
     "the network plot AND as labeled columns in the heatmap. Same "
     "selection-score columns as sheet 04."),
    ("06_NetREm_TF_TG_full",
     "ALL NetREm TF->TG edges for this cell type. Useful for follow-up "
     "analyses that go beyond the visualized intersection."),
    ("07_pySCENIC_TF_TG_full",
     "ALL pySCENIC TF->TG edges for this cell type."),
    ("", ""),
    ("Definitions", ""),
    ("role_netrem",
     "Sign of NetREm coefficient: activator (>0), repressor (<0), neutral (0)."),
    ("rank_score_netrem",
     "Per-method ranking score for NetREm (|coef|)."),
    ("rank_score_pyscenic",
     "Per-method ranking score for pySCENIC (importance / final score)."),
    ("weight_netrem / weight_pyscenic / weight_combined",
     "Min-max normalized weights in [0,1] used for plotting."),
    ("is_gwas",
     "True if the gene is in the AD/FTD GWAS list at "
     "scenic/gwas_genes_ad_ftd.xlsx (column GeneName)."),
    ("selection_score",
     "Composite score used to pick which TG columns get visible labels in "
     "the heatmap (and therefore which TGs appear in the focused network). "
     "= 1000*is_gwas + 5*scenic_degree_norm + 3*connected_to_top_TF "
     "+ in_degree_norm + 0.001*total_weight."),
]

def _is_gwas_series(genes, gwas_set):
    return [is_gwas(g, gwas_set) for g in genes]

def _heatmap_tf_summary(tf_tg_used):
    if tf_tg_used.empty:
        return pd.DataFrame(columns=["TF", "n_TGs_in_heatmap",
                                     "sum_weight_combined",
                                     "mean_weight_combined"])
    g = (tf_tg_used.groupby("TF")
         .agg(n_TGs_in_heatmap=("TG", "nunique"),
              sum_weight_combined=("weight_combined", "sum"),
              mean_weight_combined=("weight_combined", "mean"))
         .reset_index()
         .sort_values("sum_weight_combined", ascending=False))
    return g

def _heatmap_tg_summary(tf_tg_used, tf_tg_scen_all, gwas_set, top_tfs):
    if tf_tg_used.empty:
        return pd.DataFrame()
    tg_status_map = get_tg_status_from_role(tf_tg_used)
    in_degree = tf_tg_used.groupby("TG")["TF"].nunique()
    total_w = tf_tg_used.groupby("TG")["weight_combined"].sum()
    in_degree_norm = in_degree / max(in_degree.max(), 1)

    if tf_tg_scen_all is not None and "TG" in tf_tg_scen_all.columns:
        scen_deg = tf_tg_scen_all.groupby("TG").size()
        if scen_deg.max() > 0:
            scen_deg_norm = scen_deg / scen_deg.max()
        else:
            scen_deg_norm = scen_deg.astype(float)
    else:
        scen_deg = pd.Series(dtype=int)
        scen_deg_norm = pd.Series(dtype=float)

    top_tf_set = set(top_tfs[:N_TOP_TFS_FOR_CONNECTIVITY])
    connected_top = (tf_tg_used[tf_tg_used["TF"].isin(top_tf_set)]
                     .groupby("TG").size() > 0)

    rows = []
    for tg in sorted(tf_tg_used["TG"].unique()):
        gwas_flag = is_gwas(tg, gwas_set)
        sd = float(scen_deg.get(tg, 0))
        sdn = float(scen_deg_norm.get(tg, 0.0)) if len(scen_deg_norm) else 0.0
        ctt = bool(connected_top.get(tg, False))
        idg = int(in_degree.get(tg, 0))
        idn = float(in_degree_norm.get(tg, 0.0))
        tw = float(total_w.get(tg, 0.0))
        sel = (1000.0 * (1.0 if gwas_flag else 0.0)
               + 5.0 * sdn
               + 3.0 * (1.0 if ctt else 0.0)
               + 1.0 * idn
               + 0.001 * tw)
        rows.append({
            "TG": tg,
            "is_gwas": gwas_flag,
            "regulatory_status_in_heatmap": tg_status_map.get(tg, "neutral"),
            "in_degree_in_heatmap": idg,
            "in_degree_norm": idn,
            "scenic_global_degree": int(sd),
            "scenic_global_degree_norm": sdn,
            "connected_to_top_heatmap_TF": ctt,
            "total_weight_in_heatmap": tw,
            "selection_score": sel,
        })
    out = (pd.DataFrame(rows)
           .sort_values(["is_gwas", "selection_score"], ascending=[False, False])
           .reset_index(drop=True))
    return out

def _write_supplemental_workbook(stem, supp_dir, gwas_set,
                                 tf_tg_shared, tf_tg_net_all,
                                 tf_tg_scen_all, tf_tf_all,
                                 tf_tg_used, heatmap_labeled_tgs):
    """Write a Nature-style multi-sheet XLSX with all the underlying data
    behind the focused heatmap+network. Also drops a tidy CSV per sheet
    for reviewers who prefer plain text."""
    xlsx_path = supp_dir / f"Supplementary_Table_{stem}.xlsx"
    csv_dir   = supp_dir / f"{stem}_csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    # ---- Build the dataframes ----
    shared = tf_tg_shared.copy()
    if not shared.empty:
        shared["is_gwas_TG"] = _is_gwas_series(shared["TG"], gwas_set)
        shared["is_gwas_TF"] = _is_gwas_series(shared["TF"], gwas_set)
        shared = shared.sort_values(
            ["TF", "rank_score_netrem"], ascending=[True, False]
        ).reset_index(drop=True)

    tftf = tf_tf_all.copy()
    if not tftf.empty:
        tftf["is_gwas_TF1"] = _is_gwas_series(tftf["TF1"], gwas_set)
        tftf["is_gwas_TF2"] = _is_gwas_series(tftf["TF2"], gwas_set)
        tftf = tftf.sort_values(
            ["relation", "weight"], ascending=[True, False]
        ).reset_index(drop=True)

    if not tf_tg_used.empty:
        tf_strength = (tf_tg_used.groupby("TF")["weight_combined"]
                       .sum().sort_values(ascending=False))
        top_tfs_ordered = tf_strength.index.tolist()
    else:
        top_tfs_ordered = []

    tf_summary = _heatmap_tf_summary(tf_tg_used)
    if not tf_summary.empty:
        tf_summary["is_gwas"] = _is_gwas_series(tf_summary["TF"], gwas_set)

    tg_summary = _heatmap_tg_summary(tf_tg_used, tf_tg_scen_all,
                                     gwas_set, top_tfs_ordered)

    tg_focused = tg_summary[tg_summary["TG"].isin(heatmap_labeled_tgs)].copy() \
        if not tg_summary.empty else tg_summary

    netrem_full = tf_tg_net_all.sort_values(
        "rank_score_netrem", ascending=False
    ).reset_index(drop=True).copy()
    netrem_full["is_gwas_TG"] = _is_gwas_series(netrem_full["TG"], gwas_set)

    pyscenic_full = tf_tg_scen_all.sort_values(
        "rank_score_pyscenic", ascending=False
    ).reset_index(drop=True).copy()
    pyscenic_full["is_gwas_TG"] = _is_gwas_series(pyscenic_full["TG"], gwas_set)

    readme_df = pd.DataFrame(SUPP_README_ROWS[1:], columns=SUPP_README_ROWS[0])

    sheet_payload = [
        ("README",                     readme_df),
        ("01_TF_TG_shared",            shared),
        ("02_TF_TF_coordination",      tftf),
        ("03_TFs_in_heatmap",          tf_summary),
        ("04_TGs_in_heatmap",          tg_summary),
        ("05_TGs_visualized_in_network", tg_focused),
        ("06_NetREm_TF_TG_full",       netrem_full),
        ("07_pySCENIC_TF_TG_full",     pyscenic_full),
    ]

    with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as xw:
        for sheet_name, df in sheet_payload:
            df.to_excel(xw, sheet_name=sheet_name, index=False)
            ws = xw.sheets[sheet_name]
            wb = xw.book
            header_fmt = wb.add_format({
                "bold": True, "bg_color": "#222222", "font_color": "white",
                "border": 1, "align": "left",
            })
            wrap_fmt = wb.add_format({"text_wrap": True, "valign": "top"})
            ws.freeze_panes(1, 0)
            for ci, col in enumerate(df.columns):
                ws.write(0, ci, str(col), header_fmt)
                series = df[col].astype(str)
                width = min(60, max(12, int(series.map(len).max() if len(series) else 12) + 2,
                                   len(str(col)) + 2))
                ws.set_column(ci, ci, width, wrap_fmt if sheet_name == "README" else None)
            ws.autofilter(0, 0, max(0, len(df)), max(0, len(df.columns) - 1))

        # Title sheet listing celltype + counts
        info_df = pd.DataFrame({
            "field": [
                "cell_type", "n_shared_TF_TG_pairs",
                "n_TFs_in_heatmap", "n_TGs_in_heatmap",
                "n_TGs_in_focused_network", "n_TF_TF_pairs_total",
                "n_NetREm_TF_TG_full", "n_pySCENIC_TF_TG_full",
                "n_GWAS_genes_loaded",
            ],
            "value": [
                stem, len(shared),
                len(tf_summary), len(tg_summary),
                len(tg_focused), len(tftf),
                len(netrem_full), len(pyscenic_full),
                len(gwas_set),
            ],
        })
        info_df.to_excel(xw, sheet_name="00_summary", index=False)

    # Mirror to CSV
    for sheet_name, df in sheet_payload:
        df.to_csv(csv_dir / f"{sheet_name}.csv", index=False)

    print(f"  [supp] wrote {xlsx_path.name}  "
          f"({len(sheet_payload)} sheets) and {len(sheet_payload)} CSVs")

# =====================================================================
# MAIN
# =====================================================================
def main():
    gwas_set = load_gwas_genes(PATH_GWAS, GWAS_GENE_COL)

    for (stem, subdir, tf_tf_file, netrem_file, pyscenic_file,
         force_tf, force_tf_count) in CELL_TYPES:

        data_dir = BASE_DIR / subdir
        print(f"\n{'='*60}")
        print(f"  Cell type: {stem}  (dir: {subdir})")
        print(f"  v3: full heatmap + selective connector labels")
        print(f"{'='*60}")

        path_tf_tf    = data_dir / tf_tf_file
        path_netrem   = data_dir / netrem_file
        path_pyscenic = data_dir / pyscenic_file

        missing = [str(p) for p in [path_tf_tf, path_netrem, path_pyscenic]
                   if not p.exists()]
        if missing:
            print(f"  SKIPPING {stem}: missing files:")
            for m in missing:
                print(f"    {m}")
            continue

        try:
            print(f"  Loading TF-TF...")
            tf_tf = load_tf_tf_netrem(path_tf_tf, stem)
            print(f"    rows: {len(tf_tf)}, relations: "
                  f"{dict(tf_tf['relation'].value_counts())}")

            print(f"  Loading NetRem TF-TG...")
            tf_tg_net_all = load_tf_tg_netrem(path_netrem)
            print(f"    total: {len(tf_tg_net_all)}")

            print(f"  Loading pySCENIC TF-TG...")
            tf_tg_scen_all = load_tf_tg_pyscenic(path_pyscenic, stem)
            print(f"    total: {len(tf_tg_scen_all)}")

            tf_tg_net_top = (tf_tg_net_all
                .sort_values("rank_score_netrem", ascending=False)
                .head(TOP_TF_TG_PER_METHOD).copy())
            tf_tg_scen_top = (tf_tg_scen_all
                .sort_values("rank_score_pyscenic", ascending=False)
                .head(TOP_TF_TG_PER_METHOD).copy())

            print(f"  Building shared TF-TG (NetRem x pySCENIC)...")
            tf_tg_shared = build_shared_tf_tg_all(tf_tg_net_top, tf_tg_scen_top)
            print(f"    shared pairs: {len(tf_tg_shared)}")

            shared_csv = OUT_DIR / f"Shared_TF_TG_pairs_NetRem_x_pySCENIC_{stem}.csv"
            tf_tg_shared.to_csv(shared_csv, index=False)

            tf_tg_plot = tf_tg_shared.head(TOP_TF_TG_INTERSECTION).copy()

            force_tf_row = None
            if force_tf and force_tf_count > 0:
                forced = force_shared_pairs_for_tf(
                    tf_tg_net_all, tf_tg_scen_all, force_tf)
                if not forced.empty:
                    forced_top = forced.head(force_tf_count).copy()
                    tf_tg_plot = pd.concat([tf_tg_plot, forced_top],
                                           ignore_index=True)
                    tf_tg_plot = tf_tg_plot.drop_duplicates(
                        subset=["TF", "TG"]).copy()
                    tf_tg_plot = tf_tg_plot.sort_values(
                        "rank_score_netrem", ascending=False
                    ).reset_index(drop=True)
                    print(f"    [force] Added {len(forced_top)} {force_tf} pairs")
                    force_tf_row = force_tf
                else:
                    print(f"    [force] No shared pairs for {force_tf}")

            print(f"    plotting TFs: {len(tf_tg_plot['TF'].unique())}, "
                  f"TGs: {len(tf_tg_plot['TG'].unique())}")

            print(f"  Drawing heatmap...")
            tf_tg_used, heatmap_labeled_tgs = draw_heatmap(
                tf_tg_plot,
                title=(f"TF x TG heatmap (shared NetRem x pySCENIC) - {stem}\n"
                       f"(plasma + PowerNorm g={HEATMAP_GAMMA}; "
                       f"top {N_TGS_LABELED} TGs labeled via elbow/bracket connectors)"),
                outpath=OUT_DIR / f"heatmap_tf_x_tg_shared_{stem}.pdf",
                gwas_set=gwas_set,
                force_tf_row=force_tf_row,
                scen_all_for_degree=tf_tg_scen_all,
            )

            if tf_tg_used.empty:
                print(f"  No heatmap edges for {stem}, skipping network.")
                _write_supplemental_workbook(
                    stem, SUPP_DIR, gwas_set,
                    tf_tg_shared, tf_tg_net_all, tf_tg_scen_all, tf_tf,
                    tf_tg_used, heatmap_labeled_tgs)
                continue

            heatmap_tfs = sorted(tf_tg_used["TF"].unique())
            tf_tf_used = tf_tf[
                tf_tf["TF1"].isin(heatmap_tfs) & tf_tf["TF2"].isin(heatmap_tfs)
            ].copy()
            n_tg_total = tf_tg_used["TG"].nunique()
            print(f"    heatmap TFs: {len(heatmap_tfs)}, "
                  f"TF-TF edges: {len(tf_tf_used)}, "
                  f"all TGs in network: {n_tg_total} "
                  f"({len(heatmap_labeled_tgs)} get name labels, "
                  f"{n_tg_total - len(heatmap_labeled_tgs)} stay as unlabeled diamonds)")

            print(f"  Building igraph...")
            # v6: keep ALL TG nodes/edges in the network. Only restrict which
            # TGs get a NAME drawn (TFs always labeled; TGs labeled only if
            # they were also labeled in the heatmap).
            g = build_igraph(tf_tg_used, tf_tf_used, gwas_set,
                             tg_label_whitelist=heatmap_labeled_tgs)

            print(f"  Plotting network...")
            raw_png = plot_igraph_network(
                g, outfile_png=str(OUT_DIR / f"igraph_network_{stem}_raw.png"))

            print(f"  Adding legend...")
            make_legend_figure(
                raw_png, g,
                outfile_base=str(OUT_DIR / f"igraph_network_{stem}_with_legend"),
                out_stem=stem)

            print(f"  Writing supplemental workbook...")
            _write_supplemental_workbook(
                stem, SUPP_DIR, gwas_set,
                tf_tg_shared, tf_tg_net_all, tf_tg_scen_all, tf_tf,
                tf_tg_used, heatmap_labeled_tgs)

        except Exception as e:
            print(f"  ERROR for {stem}: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*60}")
    print(f"All done! Output: {OUT_DIR}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
