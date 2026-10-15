#!/usr/bin/env python3
"""scVI / scANVI reference mapping of the 13 cohorts onto the Allen brain reference (Fig. 1).

Maps the merged, per-cohort-QC'd snRNA-seq data from the 13 cohorts onto the
Allen Human Multiple Cortical Areas SMART-seq reference (49,417 nuclei,
3 donors, `subclass_label`) with scVI -> scANVI -> scArches query surgery, then
clusters the scANVI latent space with Leiden.

Two modes
---------
--mode tutorial (default)
    Follows the scvi-tools 1.3.3 tutorials "Integration and label transfer with
    Tabula Muris" and "Reference mapping with scvi-tools" (scArches):
    1. Reference restricted to genes shared with the query.
    2. SMART-seq gene-length normalization of the reference counts
       (counts / gene length * median gene length, rounded), as the Tabula
       Muris tutorial does for its Smart-seq2 data. 10x query counts are UMIs
       and are not length-normalized. Gene lengths come from a GTF (--gtf) or a
       two-column table (--gene-lengths); see --length-type.
    3. HVGs: flavor="seurat_v3" on the counts layer, batch-aware over the Allen
       donors (N_HVG genes).
    4. scVI on the counts with the scArches architecture (use_layer_norm="both",
       use_batch_norm="none", encode_covariates=True, dropout_rate=0.2,
       n_layers=2), default epochs.
    5. scANVI from the scVI model. Reference cells without a subclass are set to
       the unlabeled category "Unknown" (so they are not a class);
       max_epochs=20, n_samples_per_label=100.
    6. Query surgery: prepare_query_anndata / load_query_data, each query sample
       its own batch, max_epochs=100, weight_decay=0, check_val_every_n_epoch=10.
    7. Hard and soft predictions; reference self-prediction confusion matrix
       (the tutorial's evaluation).
    8. Doublet filter, then neighbors (cosine) on X_scANVI, UMAP, Leiden.
    scvi.settings.seed is fixed (--seed).

--mode original
    Reproduces the run behind the manuscript (scvi-tools 0.20.3, March 2024):
    no gene-length normalization; HVGs with the default flavor on
    log-normalized reference expression; reference cells without a subclass
    labeled "unknown" (a labeled class, because unlabeled_category="Unknown");
    scANVI 30 epochs; query 30 epochs; Leiden before the doublet filter; no
    seed. Use with --pretrained to load the saved models
    (scvi_model_allen_smartseq_03-11-24/, scanvi_model_allen_03-11-24/,
    scanvi_query_model_03-11-24/ in --model-dir) and reproduce the published
    latent space. The saved models load only in scvi-tools 0.20.3 (pandas < 2).

Environments: tutorial mode was tested with scvi-tools 1.3.3 / torch 2.5.1
(CUDA 12.1); original mode with --pretrained with scvi-tools 0.20.3 /
torch 1.13.1 / pandas 1.5.3 / anndata 0.9.1 / scanpy 1.9.3.

Original-run reference values:
    reference after HVG selection   49,417 x 5,000
    query                           4,807,280 cells, 913 samples
    query-model batches             916 (3 Allen donors + 913 query samples)
    Leiden, resolution 0.3          18 clusters (0-17); cluster 17 = 9,412 cells
    after doublet filter            3,556,850 cells
Cluster 17 is removed in 07_finalize_remove_cl17_rename_cl18.py, which leaves
17 Leiden clusters (0-16). A tutorial-mode run gives a new embedding, so its
cluster numbers and the cluster -> cell type map in script 00b must be
re-checked against the marker dot plot.

Usage:
    python 00_scvi_scanvi_allen_reference_mapping.py --gtf genes.gtf
    python 00_scvi_scanvi_allen_reference_mapping.py --mode original --pretrained
    python 00_scvi_scanvi_allen_reference_mapping.py --gtf genes.gtf --smoke-test   # quick CPU check

Inputs:
    allen_countmat.h5ad                   Allen Human Multiple Cortical Areas SMART-seq
                                          (49,417 x 50,281); raw read counts (exon + intron)
                                          in layers['raw_count']; obs: 'subclass_label',
                                          'external_donor_name_label'
    13ad_filtered_samples200&15000.h5ad   merged 13-cohort data after per-cohort QC
                                          (4,807,280 x 25,251); raw UMI counts in .X;
                                          obs: 'Sample', 'Study', 'doublet', 'doublet_score'
    genes.gtf                             (tutorial mode) GTF with gene_name attributes, e.g.
                                          the Cell Ranger GRCh38-2020-A genes.gtf

Outputs (--outdir / --model-dir):
    scVI, scANVI and query models (when training)
    13datasets_after_label_transfer_leiden_doublet0.1.h5ad
    reference_self_prediction_confusion.csv / .pdf      (tutorial mode)
    umap_leiden_study_predictions.png
"""

import argparse
import os
import re
from collections import defaultdict
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# This script trains in a single process. Under a multi-task SLURM allocation
# without srun, Lightning would set up distributed training and fail; a job name
# of "bash" makes it skip SLURM cluster detection.
if int(os.environ.get("SLURM_NTASKS", "1")) > 1:
    os.environ["SLURM_JOB_NAME"] = "bash"

import scvi  # noqa: E402
import torch  # noqa: E402

BATCH_KEY = "external_donor_name_label"  # Allen donor field; query samples are written into it
REF_LABEL = "subclass_label"
SCANVI_LABELS = "labels_scanvi"
UNLABELED = "Unknown"
N_HVG = 5000
LEIDEN_RES = 0.3
DOUBLET_MAX = 0.1

SETTINGS = {
    #            counts layer   gene-length norm  HVG flavor   missing label  scANVI ep  per-label  query ep
    "tutorial": dict(layer="counts", length_norm=True, hvg_flavor="seurat_v3", missing_label=UNLABELED,
                     scanvi_epochs=20, n_samples_per_label=100, query_epochs=100, doublets_first=True),
    "original": dict(layer="raw_count", length_norm=False, hvg_flavor="seurat", missing_label="unknown",
                     scanvi_epochs=30, n_samples_per_label=None, query_epochs=30, doublets_first=False),
}

REF_H5AD = "allen_countmat.h5ad"
QUERY_H5AD = "13ad_filtered_samples200&15000.h5ad"
OUT_H5AD = "13datasets_after_label_transfer_leiden_doublet0.1.h5ad"
MODEL_NAMES = {
    "tutorial": ("scvi_model_allen_tutorial", "scanvi_model_allen_tutorial", "scanvi_query_model_tutorial"),
    "original": ("scvi_model_allen_smartseq_03-11-24", "scanvi_model_allen_03-11-24", "scanvi_query_model_03-11-24"),
}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["tutorial", "original"], default="tutorial")
    p.add_argument("--ref", default=REF_H5AD, help=f"Allen SMART-seq reference h5ad (default: {REF_H5AD})")
    p.add_argument("--query", default=QUERY_H5AD, help=f"merged 13-cohort query h5ad (default: {QUERY_H5AD})")
    p.add_argument("--gtf", help="GTF with gene_name attributes, for gene lengths (tutorial mode)")
    p.add_argument("--gene-lengths", help="alternative to --gtf: TSV with columns gene, length")
    p.add_argument("--length-type", choices=["exonic", "genomic"], default="exonic",
                   help="exonic = union of exons per gene (as the tutorial's transcript lengths); "
                        "genomic = gene span, closer to exon+intron read counts (default: exonic)")
    p.add_argument("--n-hvg", type=int, default=N_HVG)
    p.add_argument("--query-epochs", type=int, help="override query epochs (tutorial 100, original 30)")
    p.add_argument("--seed", type=int, default=0, help="scvi.settings.seed (tutorial mode)")
    p.add_argument("--outdir", default=".")
    p.add_argument("--model-dir", default=".")
    p.add_argument("--pretrained", action="store_true", help="load saved models from --model-dir instead of training")
    p.add_argument("--final-atlas", help="optional final atlas h5ad for a gene/batch consistency check")
    p.add_argument("--smoke-test", action="store_true",
                   help="subsample (5k reference / 20k query cells), 1 epoch per model; checks the pipeline runs")
    args = p.parse_args()
    if args.pretrained and int(scvi.__version__.split(".")[0]) >= 1:
        p.error("the saved models were trained with scvi-tools 0.20.3 and do not load in scvi-tools >= 1.0 "
                "(scANVI label layout changed); run --pretrained in a scvi-tools 0.20.3 environment")
    if args.pretrained and args.mode != "original":
        p.error("--pretrained loads the original models; use it with --mode original")
    if args.mode == "tutorial" and not (args.gtf or args.gene_lengths):
        p.error("tutorial mode needs --gtf or --gene-lengths for SMART-seq gene-length normalization")
    return args


def gene_lengths_from_gtf(gtf, length_type):
    """Per gene_name length: union of exon intervals (exonic) or gene span (genomic)."""
    name_re = re.compile(r'gene_name "([^"]+)"')
    exons, spans = defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] not in ("exon", "gene"):
                continue
            m = name_re.search(f[8])
            if not m:
                continue
            g, start, end = m.group(1), int(f[3]), int(f[4])
            if f[2] == "gene":
                spans[g] = max(spans.get(g, 0), end - start + 1)
            else:
                exons[g].append((start, end))
    if length_type == "genomic":
        return pd.Series(spans, dtype=float)
    lengths = {}
    for g, iv in exons.items():
        iv.sort()
        total, cs, ce = 0, *iv[0]
        for s, e in iv[1:]:
            if s > ce:
                total += ce - cs + 1
                cs, ce = s, e
            else:
                ce = max(ce, e)
        lengths[g] = total + ce - cs + 1
    return pd.Series(lengths, dtype=float)


def length_normalize(ref, layer, lengths):
    """Tabula Muris tutorial: counts / gene length * median length, rounded (SMART-seq only)."""
    lengths = lengths.reindex(ref.var_names).dropna()
    lengths = lengths[lengths > 0]
    dropped = ref.n_vars - lengths.size
    ref = ref[:, lengths.index].copy()
    X = sp.csr_matrix(ref.layers[layer], dtype=np.float32)
    scale = (np.median(lengths.values) / lengths.values).astype(np.float32)
    X = X.multiply(scale[None, :]).tocsr()
    X.data = np.rint(X.data)
    X.eliminate_zeros()
    ref.layers[layer] = X
    log(f"gene-length normalized reference: {lengths.size:,} genes ({dropped:,} without a length dropped), "
        f"median length {np.median(lengths.values):,.0f} bp")
    return ref


def load_query(backed, genes, rows=None, chunk=250_000):
    """Read the backed query in row chunks, keeping only `genes` (as float32)."""
    cols = backed.var_names.get_indexer(genes)
    if (cols < 0).any():
        raise ValueError(f"{(cols < 0).sum()} reference genes missing from the query")
    rows = np.arange(backed.n_obs) if rows is None else rows
    parts = []
    for i in range(0, rows.size, chunk):
        r = rows[i:i + chunk]
        X = backed.X[r[0]:r[-1] + 1] if np.all(np.diff(r) == 1) else backed.X[r]
        parts.append(sp.csr_matrix(X[:, cols], dtype=np.float32))
        log(f"  query rows {i + r.size:,} / {rows.size:,}")
    X = sp.vstack(parts, format="csr")
    obs = backed.obs.iloc[rows].copy()
    q = ad.AnnData(X=X, obs=obs, var=backed.var.iloc[cols].copy())
    backed.file.close()
    log(f"query loaded: {q.shape}, {X.nnz:,} nonzeros")
    return q


def train_or_load(cls, path, adata, pretrained, build, train):
    if pretrained:
        return cls.load(path, adata=adata)
    model = build()
    train(model)
    model.save(path, overwrite=True)
    return model


def main():
    args = parse_args()
    cfg = dict(SETTINGS[args.mode])
    if args.query_epochs:
        cfg["query_epochs"] = args.query_epochs
    epochs = (lambda n: 1) if args.smoke_test else (lambda n: n)
    layer = cfg["layer"]
    scvi_name, scanvi_name, query_name = (os.path.join(args.model_dir, n) for n in MODEL_NAMES[args.mode])
    if args.smoke_test and not args.pretrained:
        scvi_name, scanvi_name, query_name = (n + "_smoketest" for n in (scvi_name, scanvi_name, query_name))
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.model_dir, exist_ok=True)
    sc.settings.figdir = args.outdir
    if args.mode == "tutorial":
        scvi.settings.seed = args.seed
    log(f"mode={args.mode} | scvi-tools {scvi.__version__} | torch {torch.__version__} | "
        f"CUDA {torch.cuda.is_available()} | {'pretrained' if args.pretrained else 'training'}")

    # --- 1. Load; reference restricted to shared genes ---------------------------------
    # The query (4.8M cells, ~16.5e9 nonzeros as float64) is opened backed and only
    # read after HVG selection, restricted to the reference genes (see load_query).
    ref = sc.read_h5ad(args.ref)
    query_backed = ad.read_h5ad(args.query, backed="r")
    rng = np.random.default_rng(0)
    query_rows = None
    if args.smoke_test:
        ref = ref[np.sort(rng.choice(ref.n_obs, min(5000, ref.n_obs), replace=False))].copy()
        query_rows = np.sort(rng.choice(query_backed.n_obs, min(20000, query_backed.n_obs), replace=False))
    log(f"reference {ref.shape}, query {query_backed.shape}")

    common = ref.var_names.intersection(query_backed.var_names)
    ref = ref[:, common].copy()
    log(f"shared genes: {len(common):,}")
    if layer != "raw_count":
        ref.layers[layer] = ref.layers.pop("raw_count")

    # --- 2. SMART-seq gene-length normalization (tutorial) ---------------------------------
    if cfg["length_norm"]:
        if args.gtf:
            lengths = gene_lengths_from_gtf(args.gtf, args.length_type)
        else:
            t = pd.read_csv(args.gene_lengths, sep="\t")
            lengths = pd.Series(t.iloc[:, 1].values, index=t.iloc[:, 0].astype(str), dtype=float)
        lengths = lengths[~lengths.index.duplicated()]
        ref = length_normalize(ref, layer, lengths)

    # --- 3. HVGs on the reference -----------------------------------------------------------
    if args.pretrained:
        ref_before_hvg, full_ref_genes = ref.copy(), ref.var_names
    if cfg["hvg_flavor"] == "seurat_v3":
        sc.pp.highly_variable_genes(ref, n_top_genes=args.n_hvg, flavor="seurat_v3",
                                    layer=layer, batch_key=BATCH_KEY, subset=True)
    else:
        sc.pp.normalize_total(ref, target_sum=1e4)
        sc.pp.log1p(ref)
        sc.pp.highly_variable_genes(ref, n_top_genes=args.n_hvg, batch_key=BATCH_KEY, subset=True)
    log(f"reference after HVG selection: {ref.shape}")  # original: (49417, 5000)
    if args.pretrained:
        # use exactly the genes the saved models were trained on
        saved = pd.Index(torch.load(os.path.join(scvi_name, "model.pt"), map_location="cpu")["var_names"])
        log(f"recomputed HVGs identical to saved model genes: {ref.var_names.equals(saved)}")
        missing = saved.difference(full_ref_genes)
        if len(missing):
            raise ValueError(f"{len(missing)} saved-model genes missing from reference/query")
        ref = ref_before_hvg[:, saved].copy()

    # --- 4. scVI on the reference ------------------------------------------------------------
    scvi.model.SCVI.setup_anndata(ref, batch_key=BATCH_KEY, layer=layer)
    scvi_model = train_or_load(
        scvi.model.SCVI, scvi_name, ref, args.pretrained,
        build=lambda: scvi.model.SCVI(ref, use_layer_norm="both", use_batch_norm="none",
                                      encode_covariates=True, dropout_rate=0.2, n_layers=2),
        train=lambda m: m.train(max_epochs=1) if args.smoke_test else m.train(),
    )

    # --- 5. scANVI on the reference ----------------------------------------------------------
    labels = ref.obs[REF_LABEL].astype(object).where(ref.obs[REF_LABEL].notna(), cfg["missing_label"])
    ref.obs[SCANVI_LABELS] = pd.Categorical(labels.astype(str))
    log(f"reference cells without a subclass -> '{cfg['missing_label']}': {(labels == cfg['missing_label']).sum():,}")
    scanvi_train = {"max_epochs": epochs(cfg["scanvi_epochs"])}
    if cfg["n_samples_per_label"]:
        scanvi_train["n_samples_per_label"] = cfg["n_samples_per_label"]
    scanvi_model = train_or_load(
        scvi.model.SCANVI, scanvi_name, ref, args.pretrained,
        build=lambda: scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category=UNLABELED,
                                                        labels_key=SCANVI_LABELS),
        train=lambda m: m.train(**scanvi_train),
    )

    if args.mode == "tutorial":
        # tutorial evaluation: reference predictions vs. known labels
        known = ref.obs[SCANVI_LABELS] != UNLABELED
        pred = pd.Series(scanvi_model.predict(ref), index=ref.obs_names)
        conf = pd.crosstab(ref.obs.loc[known, SCANVI_LABELS], pred[known], normalize="index")
        conf.to_csv(os.path.join(args.outdir, "reference_self_prediction_confusion.csv"))
        acc = (pred[known] == ref.obs.loc[known, SCANVI_LABELS].astype(str)).mean()
        log(f"reference self-prediction accuracy: {acc:.3f}")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.pcolor(conf.values, cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(np.arange(conf.shape[1]) + 0.5, conf.columns, rotation=90)
        ax.set_yticks(np.arange(conf.shape[0]) + 0.5, conf.index)
        ax.set_xlabel("Predicted"); ax.set_ylabel("Allen subclass")
        fig.colorbar(im, ax=ax); fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, "reference_self_prediction_confusion.pdf")); plt.close(fig)

    # --- 6. Query surgery (batch = sample) -----------------------------------------------------
    query = load_query(query_backed, ref.var_names, query_rows)
    query.layers[layer] = query.X
    query.obs[BATCH_KEY] = query.obs["Sample"].astype(str).values
    scvi.model.SCANVI.prepare_query_anndata(query, scanvi_model)
    if args.pretrained:
        # load_query_data adds this (all unlabeled) when training; a saved query model needs it present
        query.obs[SCANVI_LABELS] = pd.Categorical([UNLABELED] * query.n_obs)
    log(f"query after gene alignment: {query.shape}")  # original: (4807280, 5000)
    query_model = train_or_load(
        scvi.model.SCANVI, query_name, query, args.pretrained,
        build=lambda: scvi.model.SCANVI.load_query_data(query, scanvi_model),
        train=lambda m: m.train(max_epochs=epochs(cfg["query_epochs"]),
                                plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10),
    )

    # --- 7. Latent space and predictions --------------------------------------------------------
    query.obsm["X_scANVI"] = query_model.get_latent_representation()
    soft = query_model.predict(soft=True)
    query.obs["predictions_scanvi"] = pd.Categorical(soft.idxmax(axis=1).values)
    query.obs["predictions_scanvi_prob"] = soft.max(axis=1).values
    log("predicted labels:\n" + query.obs["predictions_scanvi"].value_counts().to_string())
    log(f"median max prediction probability: {np.median(query.obs['predictions_scanvi_prob']):.3f}")

    # --- 8. Doublets, neighbors, UMAP, Leiden ----------------------------------------------------
    def drop_doublets(a):
        a = a[a.obs["doublet_score"] <= DOUBLET_MAX].copy()
        log(f"after doublet filter (<= {DOUBLET_MAX}): {a.n_obs:,} cells")  # original: 3,556,850
        return a

    if cfg["doublets_first"]:
        query = drop_doublets(query)
    sc.pp.neighbors(query, use_rep="X_scANVI", metric="cosine")
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=LEIDEN_RES)
    counts = query.obs["leiden"].value_counts().sort_index(key=lambda s: s.astype(int))
    log(f"{counts.size} Leiden clusters at resolution {LEIDEN_RES}:\n{counts.to_string()}")  # original: 18 (0-17)
    if not cfg["doublets_first"]:
        query = drop_doublets(query)
    query.obs["doublet"] = query.obs["doublet"].astype(str)

    out_h5 = os.path.join(args.outdir, OUT_H5AD.replace(".h5ad", "_smoketest.h5ad") if args.smoke_test else OUT_H5AD)
    query.write_h5ad(out_h5)
    log(f"wrote {out_h5}")
    sc.pl.umap(query, color=["leiden", "Study", "predictions_scanvi"], frameon=False, ncols=1,
               legend_fontoutline=1, legend_fontsize=9, show=False, save="_leiden_study_predictions.png")

    # --- Optional consistency check against the final atlas ----------------------------------------
    if args.final_atlas:
        atlas = sc.read_h5ad(args.final_atlas, backed="r")
        log(f"genes identical & same order as final atlas: {list(atlas.var_names) == list(query.var_names)}")
        model_batches = set(query_model.adata_manager.get_state_registry("batch")["categorical_mapping"])
        atlas_batches = set(atlas.obs[BATCH_KEY].astype(str)) - {"unknown"}  # "unknown" = rescued ENDO cells
        log(f"final-atlas batches found in query model: {len(atlas_batches & model_batches)} / {len(atlas_batches)}")
        atlas.file.close()


if __name__ == "__main__":
    main()
