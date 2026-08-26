"""
scenic_ctx_aucell_viz.py
========================
Full pySCENIC pipeline: ctx (regulon prediction), aucell (cellular enrichment),
and downstream visualizations (UMAP, heatmaps, Z-score clustermaps).

Starts from pre-computed filtered adjacency CSVs per sex-diagnosis group.

Inputs:
    - Filtered adjacency CSVs:  <SCENIC_DIR>/adata_{group}_adj_scenic_filteredGenes.csv
    - Loom files:               <SCENIC_DIR>/adata_{group}_scenic.loom
    - Ranking database:         <SCENIC_DIR>/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    - Motif annotation:         <SCENIC_DIR>/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
    - h5ad metadata files:      <H5AD_DIR>/12adata_{group}.h5ad

Outputs:
    - {group}_reg.csv           (ctx regulon CSVs)
    - {group}_auc_matrix.csv    (aucell output)
    - {group}_regulons.dat      (pickled regulon objects)
    - UMAP PDFs, heatmap PDFs, Z-score clustermap PDF

Usage:
    python scenic_ctx_aucell_viz.py
"""

import os
import pickle
import operator as op

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from cytoolz import compose

from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons

# --- Configuration -----------------------------------------------------------

BASE_DIR = "."
SCENIC_DIR = BASE_DIR + "/scenic"
ADJ_DIR = SCENIC_DIR
H5AD_DIR = BASE_DIR + "/h5ad"

F_DB = os.path.join(SCENIC_DIR, "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
F_MOTIF = os.path.join(SCENIC_DIR, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

H5AD_FILES = {
    "ct_male":    os.path.join(H5AD_DIR, "12adata_ct_male.h5ad"),
    "ct_female":  os.path.join(H5AD_DIR, "12adata_ct_female.h5ad"),
    "mci_male":   os.path.join(H5AD_DIR, "12adata_mci_male.h5ad"),
    "mci_female": os.path.join(H5AD_DIR, "12adata_mci_female.h5ad"),
    "dm_male":    os.path.join(H5AD_DIR, "12adata_dem_male.h5ad"),
    "dm_female":  os.path.join(H5AD_DIR, "12adata_dem_female.h5ad"),
}

GROUPS = ["ct_male", "ct_female", "mci_male", "mci_female", "dm_male", "dm_female"]
NUM_WORKERS = 50
TOP_N_REGULONS = 30

sc.settings.figdir = SCENIC_DIR
sc.settings.verbosity = 1

# --- Helper: load h5ad efficiently ------------------------------------------

def load_h5ad_light(path):
    """
    Load an h5ad in backed mode to avoid reading the full expression matrix.
    Extracts obs (metadata) and obsm (embeddings like UMAP) into a lightweight
    AnnData, then closes the file handle. Fast even for million-cell datasets.
    """
    print(f"  Loading (backed): {os.path.basename(path)}")
    adata_backed = sc.read_h5ad(path, backed="r")
    adata_light = ad.AnnData(
        obs=adata_backed.obs.copy(),
        obsm={k: np.array(v) for k, v in adata_backed.obsm.items()},
    )
    adata_backed.file.close()
    return adata_light

# --- Step 2+3: pyscenic ctx (regulon prediction) ----------------------------

print("=" * 60)
print("STEP 2+3: pyscenic ctx")
print("=" * 60)

for group in GROUPS:
    adj_csv = os.path.join(ADJ_DIR, f"adata_{group}_adj_scenic_filteredGenes.csv")
    loom = os.path.join(SCENIC_DIR, f"adata_{group}_scenic.loom")
    reg_csv = os.path.join(SCENIC_DIR, f"{group}_reg.csv")

    if os.path.exists(reg_csv):
        print(f"[SKIP] {group}_reg.csv already exists")
        continue

    if not os.path.exists(adj_csv):
        print(f"[SKIP] adjacency file not found: {adj_csv}")
        continue

    cmd = (
        f"pyscenic ctx {adj_csv} {F_DB} "
        f"--annotations_fname {F_MOTIF} "
        f"--expression_mtx_fname {loom} "
        f"--output {reg_csv} "
        f"--mask_dropouts "
        f"--num_workers {NUM_WORKERS}"
    )
    print(f"Running ctx for {group} ...")
    os.system(cmd)
    print(f"Done: {reg_csv}")

# --- Step 4: pyscenic aucell (cellular enrichment) --------------------------

print("\n" + "=" * 60)
print("STEP 4: pyscenic aucell")
print("=" * 60)

for group in GROUPS:
    loom = os.path.join(SCENIC_DIR, f"adata_{group}_scenic.loom")
    reg_csv = os.path.join(SCENIC_DIR, f"{group}_reg.csv")
    auc_csv = os.path.join(SCENIC_DIR, f"{group}_auc_matrix.csv")

    if os.path.exists(auc_csv):
        print(f"[SKIP] {group}_auc_matrix.csv already exists")
        continue

    if not os.path.exists(reg_csv):
        print(f"[SKIP] regulon file not found: {reg_csv}")
        continue

    cmd = (
        f"pyscenic aucell {loom} {reg_csv} "
        f"--output {auc_csv} "
        f"--num_workers {NUM_WORKERS}"
    )
    print(f"Running aucell for {group} ...")
    os.system(cmd)
    print(f"Done: {auc_csv}")

# --- Helper: derive regulons from motif CSV ---------------------------------

def derive_regulons(motifs):
    """Filter motifs and derive regulons with NES >= 3.0 and >= 10 target genes."""
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains("weight>75.0%")), motifs.Context), dtype=bool) &
        np.fromiter(map(contains("activating"), motifs.Context), dtype=bool)
    ]

    regulons = list(filter(
        lambda r: len(r) >= 10,
        df2regulons(motifs[
            (motifs["NES"] >= 3.0) &
            (
                (motifs["Annotation"] == "gene is directly annotated") |
                (
                    motifs["Annotation"].str.startswith("gene is orthologous to") &
                    motifs["Annotation"].str.endswith("which is directly annotated for motif")
                )
            )
        ])
    ))
    return [r.rename(r.transcription_factor) for r in regulons]

# --- Visualization functions -------------------------------------------------

def plot_umap_by_regulons(auc_mtx, group, top_n=TOP_N_REGULONS):
    """
    Color the existing h5ad UMAP by top N regulon AUC scores.
    If the h5ad is unavailable, falls back to computing a new UMAP from the
    AUC matrix. Saves one PDF per batch of 10 regulons.
    """
    h5ad_path = H5AD_FILES.get(group)

    if h5ad_path and os.path.exists(h5ad_path):
        adata_auc = load_h5ad_light(h5ad_path)
        if "X_umap" not in adata_auc.obsm:
            print(f"  WARNING: no X_umap in {h5ad_path}, computing from AUC matrix instead")
            adata_auc = None
    else:
        print(f"  h5ad not found for {group}, computing UMAP from AUC matrix ...")
        adata_auc = None

    if adata_auc is None:
        adata_auc = ad.AnnData(
            X=auc_mtx.values,
            obs=pd.DataFrame(index=auc_mtx.index),
            var=pd.DataFrame(index=auc_mtx.columns),
        )
        sc.pp.neighbors(adata_auc, use_rep="X")
        sc.tl.umap(adata_auc)
    else:
        # Align cell barcodes between AUC matrix and h5ad
        common_cells = adata_auc.obs_names.intersection(auc_mtx.index)
        adata_auc = adata_auc[common_cells].copy()
        auc_mtx = auc_mtx.loc[common_cells]

    # Add AUC scores to obs for coloring
    for col in auc_mtx.columns:
        adata_auc.obs[col] = auc_mtx[col].values

    # Top N regulons by mean activity
    top_regulons = auc_mtx.mean(axis=0).nlargest(top_n).index.tolist()

    for i in range(0, len(top_regulons), 10):
        batch = top_regulons[i:i + 10]
        sc.pl.umap(
            adata_auc,
            color=batch,
            ncols=5,
            show=False,
            save=f"_{group}_regulons_umap_batch{i // 10 + 1}.pdf",
        )
    print(f"  UMAP plots saved for {group}")


def plot_regulon_activity_heatmap(auc_mtx, group, top_n=50):
    """
    Heatmap: rows = top N regulons (TFs), columns = subsampled cells,
    values = AUC activity scores.
    """
    top_regulons = auc_mtx.mean(axis=0).nlargest(top_n).index.tolist()
    sub = auc_mtx[top_regulons]

    # Subsample cells if large
    if len(sub) > 5000:
        sub = sub.sample(n=5000, random_state=42)

    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(
        sub.T,
        cmap="YlOrRd",
        xticklabels=False,
        yticklabels=True,
        ax=ax,
    )
    ax.set_title(f"{group} — top {top_n} regulon AUC scores (cells subsampled)")
    ax.set_xlabel("Cells")
    ax.set_ylabel("Regulon (TF)")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    out_path = os.path.join(SCENIC_DIR, f"{group}_regulon_activity_heatmap.pdf")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Regulon activity heatmap saved: {out_path}")


def plot_regulon_size_heatmap(regulons, auc_mtx, group, top_n=30):
    """
    Heatmap summarising each TF regulon:
      - number of target genes
      - mean AUC activity score
    Rows = top N TFs by mean AUC.
    """
    records = []
    for r in regulons:
        tf = r.transcription_factor
        col = next((c for c in auc_mtx.columns if c.startswith(tf + "(")), None)
        records.append({
            "TF": tf,
            "n_targets": len(r.genes),
            "mean_AUC": auc_mtx[col].mean() if col else np.nan,
        })

    df = (
        pd.DataFrame(records)
        .dropna(subset=["mean_AUC"])
        .sort_values("mean_AUC", ascending=False)
        .head(top_n)
        .set_index("TF")
    )

    fig, ax = plt.subplots(figsize=(5, top_n * 0.35 + 1))
    sns.heatmap(
        df,
        cmap="YlOrRd",
        annot=True,
        fmt=".3g",
        linewidths=0.4,
        ax=ax,
    )
    ax.set_title(f"{group} — top {top_n} regulons")
    out_path = os.path.join(SCENIC_DIR, f"{group}_regulon_summary_heatmap.pdf")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Regulon summary heatmap saved: {out_path}")


def plot_zscore_clustermap(auc_dict):
    """
    Clustermap: rows = regulons, columns = groups,
    values = Z-scored mean AUC across cells.
    Shows top 50 most variable regulons across groups.
    """
    zscore_series = {}
    for group, auc_mtx in auc_dict.items():
        mean_auc = auc_mtx.mean(axis=0)
        z = (mean_auc - mean_auc.mean()) / mean_auc.std()
        zscore_series[group] = z

    merged = pd.DataFrame(zscore_series).fillna(0)

    # Keep top 50 most variable regulons
    top_var = merged.var(axis=1).nlargest(50).index
    merged = merged.loc[top_var]

    g = sns.clustermap(
        merged,
        cmap="vlag",
        center=0,
        figsize=(10, 18),
        yticklabels=True,
        col_cluster=True,
        row_cluster=True,
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
    g.ax_heatmap.set_title("Regulon Z-score across groups", pad=10)
    out_path = os.path.join(SCENIC_DIR, "all_groups_regulon_zscore_clustermap.pdf")
    g.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Z-score clustermap saved: {out_path}")

# --- Main: load results and generate all plots ------------------------------

print("\n" + "=" * 60)
print("VISUALIZATIONS")
print("=" * 60)

auc_dict = {}

for group in GROUPS:
    reg_csv = os.path.join(SCENIC_DIR, f"{group}_reg.csv")
    auc_csv = os.path.join(SCENIC_DIR, f"{group}_auc_matrix.csv")

    if not os.path.exists(reg_csv) or not os.path.exists(auc_csv):
        print(f"[SKIP] Missing output files for {group} — run ctx/aucell first")
        continue

    print(f"\n--- {group} ---")

    # Load AUC matrix
    auc_mtx = pd.read_csv(auc_csv, index_col=0)
    auc_dict[group] = auc_mtx
    print(f"  AUC matrix: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")

    # Derive regulons
    df_motifs = load_motifs(reg_csv)
    regulons = derive_regulons(df_motifs)
    print(f"  {len(regulons)} regulons derived")

    # Save pickled regulons
    dat_path = os.path.join(SCENIC_DIR, f"{group}_regulons.dat")
    with open(dat_path, "wb") as fh:
        pickle.dump(regulons, fh)

    # Generate plots
    plot_umap_by_regulons(auc_mtx, group)
    plot_regulon_activity_heatmap(auc_mtx, group)
    plot_regulon_size_heatmap(regulons, auc_mtx, group)

# --- Cross-group Z-score clustermap -----------------------------------------

if len(auc_dict) > 1:
    print("\n--- Cross-group Z-score clustermap ---")
    plot_zscore_clustermap(auc_dict)

print("\nAll done!")
