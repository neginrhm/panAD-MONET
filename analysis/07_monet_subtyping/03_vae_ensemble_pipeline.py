#!/usr/bin/env python3
"""
VAE Ensemble Pipeline
=====================
1. Train masked VAE across 20 seeds
2. GMM clustering (k=4) per seed; silhouette scores for k=2-10 recorded per
   seed (k=2-3 scored highest; k=4 was used, see manuscript Methods)
3. Consensus clustering -> final labels
4. Sankey diagrams: per-seed -> final mapping

Requirements:
    pip install torch numpy pandas scikit-learn scipy matplotlib \
                plotly "kaleido==0.2.1" umap-learn tqdm

Usage:
    Provide sample_vectors_padded (dict[sample_id -> padded vector])
    and sample_vectors_mask (dict[sample_id -> mask vector]) before
    calling run_pipeline(), or modify load_data() to load from disk.
"""

import os
import re
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm
import plotly.graph_objects as go

# =====================================================================
# CONFIG
# =====================================================================
OUTPUT_DIR = "vae_ensemble_output"
K_FINAL = 4
K_RANGE_SILHOUETTE = range(2, 11)
SEEDS = list(range(20))
LATENT_DIM = 128
BATCH_SIZE = 32
LR = 1e-4
MAX_EPOCHS = 200
DELTA = 1e-4

# Sankey aesthetics
GRID_ROWS = 4
GRID_COLS = 5
SINGLE_W = 1000
SINGLE_H = 600
GRID_W = GRID_COLS * 550
GRID_H = GRID_ROWS * 430
FONT_FAMILY = "Helvetica, Arial, sans-serif"
FONT_SIZE = 10
NODE_THICK = 12
NODE_PAD = 10
LINK_ALPHA = 0.55
SEED_NODE_ALPHA = 0.30

# Okabe-Ito colorblind-safe palette (final clusters 1..4)
PALETTE = {
    1: "#0072B2",   # blue
    2: "#E69F00",   # orange
    3: "#D55E00",   # vermilion
    4: "#009E73",   # bluish green
}


# =====================================================================
# MODEL
# =====================================================================
class VAE(nn.Module):
    def __init__(self, input_dim, latent_dim=128):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 8192), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(8192, 4096), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(4096, 2048), nn.ReLU(),
            nn.Linear(2048, 1024), nn.ReLU(),
            nn.Linear(1024, 512), nn.ReLU(),
            nn.Linear(512, 256), nn.ReLU(),
        )
        self.fc_mu = nn.Linear(256, latent_dim)
        self.fc_logvar = nn.Linear(256, latent_dim)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256), nn.ReLU(),
            nn.Linear(256, 512), nn.ReLU(),
            nn.Linear(512, 1024), nn.ReLU(),
            nn.Linear(1024, 2048), nn.ReLU(),
            nn.Linear(2048, 4096), nn.ReLU(),
            nn.Linear(4096, 8192), nn.ReLU(),
            nn.Linear(8192, input_dim),
        )

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar


class GeneDataset(Dataset):
    def __init__(self, data, mask):
        self.data = torch.tensor(data, dtype=torch.float32)
        self.mask = torch.tensor(mask, dtype=torch.float32)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.mask[idx]


def masked_vae_loss(recon_x, x, mask, mu, logvar):
    diff = (recon_x - x) * mask
    mse = torch.sum(diff ** 2) / torch.sum(mask)
    kld = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return mse + kld


# =====================================================================
# STEP 1: TRAIN VAE + GMM OVER 20 SEEDS
# =====================================================================
def train_all_seeds(X, mask_array, sample_ids):
    """Train VAE for each seed, cluster with GMM, return per-seed labels."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    vae_dir = os.path.join(OUTPUT_DIR, "vae_runs")
    os.makedirs(vae_dir, exist_ok=True)

    # columns = per-seed cluster labels, rows = samples
    all_labels = pd.DataFrame({"SampleID": sample_ids})
    silhouette_rows = []  # per-seed silhouette scores across candidate k

    for seed in SEEDS:
        print(f"\n=== SEED {seed} ===")
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

        dataset = GeneDataset(X, mask_array)
        train_size = int(0.8 * len(dataset))
        val_size = len(dataset) - train_size
        train_ds, val_ds = random_split(dataset, [train_size, val_size])
        train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)

        model = VAE(input_dim=X.shape[1], latent_dim=LATENT_DIM).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=LR)

        best_loss = float("inf")
        best_state = None
        train_losses = []

        for epoch in range(MAX_EPOCHS):
            model.train()
            epoch_loss = 0.0
            for batch_data, batch_mask in train_loader:
                batch_data, batch_mask = batch_data.to(device), batch_mask.to(device)
                optimizer.zero_grad()
                recon, mu, logvar = model(batch_data)
                loss = masked_vae_loss(recon, batch_data, batch_mask, mu, logvar)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                epoch_loss += loss.item() * batch_data.size(0)

            avg_loss = epoch_loss / len(train_loader.dataset)
            train_losses.append(avg_loss)
            if (epoch + 1) % 20 == 0:
                print(f"  Epoch {epoch+1}/{MAX_EPOCHS} | Loss: {avg_loss:.6f}")

            if best_loss - avg_loss > DELTA:
                best_loss = avg_loss
                best_state = {k: v.detach().cpu().clone()
                              for k, v in model.state_dict().items()}

        # save checkpoint
        if best_state is not None:
            torch.save(best_state, os.path.join(vae_dir, f"vae_seed_{seed}.pth"))

        # save loss curve
        plt.figure()
        plt.plot(train_losses)
        plt.xlabel("Epoch"); plt.ylabel("Training Loss")
        plt.title(f"Seed {seed}"); plt.tight_layout()
        plt.savefig(os.path.join(vae_dir, f"loss_seed_{seed}.png"))
        plt.close()

        # --- embed with best model ---
        if best_state is not None:
            model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            data_t = torch.tensor(X, dtype=torch.float32, device=device)
            mu_all, logvar_all = model.encode(data_t)
            z = model.reparameterize(mu_all, logvar_all).cpu().numpy()

        # --- silhouette across candidate k (reported; not used to pick k) ---
        for k in K_RANGE_SILHOUETTE:
            labels_k = GaussianMixture(n_components=k, random_state=42).fit_predict(z)
            silhouette_rows.append({"seed": seed, "k": k,
                                    "silhouette": silhouette_score(z, labels_k)})

        # --- GMM k=4 ---
        gmm = GaussianMixture(n_components=K_FINAL, random_state=42)
        labels = gmm.fit_predict(z)
        all_labels[f"Seed_{seed}"] = labels

    all_labels.to_csv(os.path.join(OUTPUT_DIR, "all_seed_labels.csv"), index=False)

    sil = pd.DataFrame(silhouette_rows)
    sil.to_csv(os.path.join(OUTPUT_DIR, "silhouette_by_seed_k.csv"), index=False)
    summary = sil.groupby("k")["silhouette"].agg(["mean", "std"])
    print("Silhouette by k (mean over seeds):\n" + summary.to_string())
    best_k = sil.loc[sil.groupby("seed")["silhouette"].idxmax(), "k"].value_counts().sort_index()
    print("Best k per seed (count of seeds):\n" + best_k.to_string())
    plt.figure(figsize=(5, 3.5))
    plt.errorbar(summary.index, summary["mean"], yerr=summary["std"], marker="o", capsize=3)
    plt.axvline(K_FINAL, color="grey", ls="--", lw=0.8)
    plt.xlabel("Number of GMM components (k)"); plt.ylabel("Silhouette score")
    plt.title(f"Silhouette across {len(SEEDS)} seeds"); plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "silhouette_by_k.pdf")); plt.close()
    return all_labels


# =====================================================================
# STEP 2: CONSENSUS CLUSTERING -> FINAL LABELS
# =====================================================================
def consensus_clustering(all_labels):
    """Build co-assignment matrix across seeds, cluster into K_FINAL."""
    sample_ids = all_labels["SampleID"].astype(str).tolist()
    seed_cols = [c for c in all_labels.columns if c.startswith("Seed_")]
    L = all_labels[seed_cols].to_numpy(dtype=float)
    N, S = L.shape

    # pairwise co-assignment matrix
    M_num = np.zeros((N, N), dtype=float)
    M_den = np.zeros((N, N), dtype=float)
    for s in range(S):
        lab = L[:, s]
        valid = ~np.isnan(lab)
        if valid.sum() < 2:
            continue
        both = np.outer(valid, valid)
        eq = np.equal.outer(lab, lab) & both
        M_num += eq.astype(float)
        M_den += both.astype(float)

    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.divide(M_num, M_den, out=np.zeros_like(M_num), where=M_den > 0)
    np.fill_diagonal(M, 1.0)
    D = 1.0 - M

    # agglomerative clustering on distance matrix
    try:
        agg = AgglomerativeClustering(
            n_clusters=K_FINAL, linkage="average", metric="precomputed")
    except TypeError:
        agg = AgglomerativeClustering(
            n_clusters=K_FINAL, linkage="average", affinity="precomputed")
    final0 = agg.fit_predict(D)          # 0-indexed
    final_labels = (final0 + 1).astype(int)  # 1-indexed

    # per-sample agreement fraction
    agree = np.zeros(N, dtype=float)
    for si in range(S):
        seed_lab = L[:, si]
        valid = ~np.isnan(seed_lab)
        if valid.sum() == 0:
            continue
        seed_vals = np.unique(seed_lab[valid])
        conf = np.zeros((len(seed_vals), K_FINAL), dtype=int)
        for r, v in enumerate(seed_vals):
            mask = valid & (seed_lab == v)
            conf[r, :] = np.bincount(final0[mask], minlength=K_FINAL)
        ri, ci = linear_sum_assignment(-conf)
        mapping = {seed_vals[r]: c for r, c in zip(ri, ci)}
        mapped = -np.ones(N, dtype=int)
        mapped[valid] = [mapping[x] for x in seed_lab[valid]]
        agree += (mapped == final0).astype(float)
    agree /= max(1, S)

    # save outputs
    pd.DataFrame(M, index=sample_ids, columns=sample_ids).to_csv(
        os.path.join(OUTPUT_DIR, "pairwise_coassignment.csv"))

    final_df = pd.DataFrame({
        "SampleID": sample_ids,
        "FinalCluster": final_labels,
        "agreement_fraction": np.round(agree, 4),
    })
    final_df.to_csv(os.path.join(OUTPUT_DIR, "final_assignments.csv"), index=False)

    print(f"\nConsensus clustering done. Cluster sizes:")
    for k in range(1, K_FINAL + 1):
        print(f"  Cluster {k}: {(final_labels == k).sum()}")

    return final_df, final0, M


# =====================================================================
# STEP 3: SANKEY DIAGRAMS
# =====================================================================
def _hex_to_rgba(h, a):
    h = h.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{a})"


def _evenly_spaced(n):
    if n <= 0:
        return []
    step = 1.0 / (n + 1)
    return [step * (i + 1) for i in range(n)]


def build_sankeys(all_labels, final0):
    """Sankey diagrams mapping each seed's clusters to final consensus."""
    sankey_dir = os.path.join(OUTPUT_DIR, "sankey")
    os.makedirs(sankey_dir, exist_ok=True)

    sample_ids = all_labels["SampleID"].astype(str).tolist()
    seed_cols = [c for c in all_labels.columns if c.startswith("Seed_")]
    L = all_labels[seed_cols].to_numpy(dtype=float)
    N, S = L.shape

    final_order = list(range(K_FINAL))
    final_y = _evenly_spaced(K_FINAL)

    # --- per-seed figures ---
    for si, col in enumerate(seed_cols):
        a = L[:, si]
        valid = ~np.isnan(a)
        if valid.sum() == 0:
            continue
        seed_short = f"S{si}"
        avals = np.unique(a[valid]).astype(int)

        # sort seed clusters by majority final destination
        cluster_meta = []
        for av in avals:
            mask = valid & (a == av)
            size = int(mask.sum())
            dest_counts = np.bincount(final0[mask], minlength=K_FINAL)
            maj_dest = int(dest_counts.argmax())
            cluster_meta.append((av, maj_dest, size))
        cluster_meta.sort(key=lambda t: (t[1], -t[2]))

        node_labels, node_colors, node_x, node_y = [], [], [], []
        node_idx = {}

        def _add(label, color, x, y):
            i = len(node_labels)
            node_labels.append(label); node_colors.append(color)
            node_x.append(x); node_y.append(y)
            node_idx[label] = i
            return i

        # seed nodes (left)
        sy = _evenly_spaced(len(cluster_meta))
        seed_node = {}
        for pos, (av, maj, size) in enumerate(cluster_meta):
            lbl = f"{seed_short}:C{av} (n={size})"
            seed_node[av] = _add(lbl, _hex_to_rgba(PALETTE[maj+1], SEED_NODE_ALPHA), 0.0, sy[pos])

        # final nodes (right)
        final_node = {}
        for j, f in enumerate(final_order):
            size_f = int(np.sum(final0 == f))
            final_node[f] = _add(f"Final:{f+1} (n={size_f})", PALETTE[f+1], 1.0, final_y[j])

        sources, targets, values, link_colors = [], [], [], []
        for av, _, _ in cluster_meta:
            for f in final_order:
                cnt = int(np.sum(valid & (a == av) & (final0 == f)))
                if cnt > 0:
                    sources.append(seed_node[av])
                    targets.append(final_node[f])
                    values.append(cnt)
                    link_colors.append(_hex_to_rgba(PALETTE[f+1], LINK_ALPHA))

        fig = go.Figure(go.Sankey(
            arrangement="fixed",
            node=dict(label=node_labels, color=node_colors,
                      x=node_x, y=node_y,
                      thickness=NODE_THICK, pad=NODE_PAD,
                      line=dict(color="rgba(0,0,0,0.35)", width=0.5)),
            link=dict(source=sources, target=targets, value=values, color=link_colors),
        ))
        fig.update_layout(
            title=dict(text=f"{seed_short} -> Final (K={K_FINAL})", x=0.02, y=0.98),
            font=dict(family=FONT_FAMILY, size=FONT_SIZE),
            margin=dict(l=10, r=10, t=40, b=10),
            paper_bgcolor="white", width=SINGLE_W, height=SINGLE_H,
        )
        base = os.path.join(sankey_dir, f"{col}_to_final")
        fig.write_image(base + ".png", scale=2)
        fig.write_image(base + ".pdf", scale=2)

    # --- grid figure (all seeds) ---
    fig_grid = go.Figure()
    for si, col in enumerate(seed_cols[:GRID_ROWS * GRID_COLS]):
        a = L[:, si]; valid = ~np.isnan(a)
        if valid.sum() == 0:
            continue
        seed_short = f"S{si}"
        avals = np.unique(a[valid]).astype(int)

        cluster_meta = []
        for av in avals:
            mask = valid & (a == av)
            size = int(mask.sum())
            dest_counts = np.bincount(final0[mask], minlength=K_FINAL)
            cluster_meta.append((av, int(dest_counts.argmax()), size))
        cluster_meta.sort(key=lambda t: (t[1], -t[2]))

        node_labels, node_colors, node_x, node_y, node_idx = [], [], [], [], {}
        def _add_g(label, color, x, y):
            i = len(node_labels)
            node_labels.append(label); node_colors.append(color)
            node_x.append(x); node_y.append(y)
            node_idx[label] = i
            return i

        sy = _evenly_spaced(len(cluster_meta))
        for pos, (av, maj, sz) in enumerate(cluster_meta):
            _add_g(f"{seed_short}:C{av}", _hex_to_rgba(PALETTE[maj+1], SEED_NODE_ALPHA), 0.0, sy[pos])
        for j, f in enumerate(final_order):
            _add_g(f"Final:{f+1}", PALETTE[f+1], 1.0, final_y[j])

        seed_map = {av: node_idx[f"{seed_short}:C{av}"] for av,_,_ in cluster_meta}
        final_map = {f: node_idx[f"Final:{f+1}"] for f in final_order}

        sources, targets, values, link_colors = [], [], [], []
        for av, _, _ in cluster_meta:
            for f in final_order:
                cnt = int(np.sum(valid & (a == av) & (final0 == f)))
                if cnt > 0:
                    sources.append(seed_map[av])
                    targets.append(final_map[f])
                    values.append(cnt)
                    link_colors.append(_hex_to_rgba(PALETTE[f+1], LINK_ALPHA))

        r_idx = si // GRID_COLS; c_idx = si % GRID_COLS
        x0 = c_idx / GRID_COLS; x1 = (c_idx + 1) / GRID_COLS
        y_top = 1 - (r_idx / GRID_ROWS); y_bot = 1 - ((r_idx + 1) / GRID_ROWS)
        domain = dict(x=[x0 + 0.02/GRID_COLS, x1 - 0.02/GRID_COLS],
                      y=[y_bot + 0.03/GRID_ROWS, y_top - 0.03/GRID_ROWS])

        fig_grid.add_trace(go.Sankey(
            domain=domain, arrangement="fixed",
            node=dict(label=node_labels, color=node_colors, x=node_x, y=node_y,
                      thickness=NODE_THICK, pad=NODE_PAD,
                      line=dict(color="rgba(0,0,0,0.35)", width=0.4)),
            link=dict(source=sources, target=targets, value=values, color=link_colors),
            name=seed_short,
        ))

    fig_grid.update_layout(
        title=dict(text=f"All seeds -> Final (K={K_FINAL})", x=0.02, y=0.98),
        font=dict(family=FONT_FAMILY, size=FONT_SIZE),
        margin=dict(l=10, r=10, t=50, b=10),
        paper_bgcolor="white", width=GRID_W, height=GRID_H,
    )
    base = os.path.join(sankey_dir, "all_seeds_grid")
    fig_grid.write_image(base + ".png", scale=2)
    fig_grid.write_image(base + ".pdf", scale=2)

    print(f"Sankey figures saved to: {sankey_dir}")


# =====================================================================
# DATA LOADING  (from mlp_masked-Copy2_ensemble.ipynb, Cell 101 + 125)
# =====================================================================
def load_data(csv_path):
    """
    Parse VectorsofFloats.csv into padded sample vectors and binary masks.

    Each cell in the CSV is a 7-element string list [AST,MIC,ENDO,ODC,OPC,INH,EXC].
    - Genes with >1 NaN across cell types are skipped.
    - Genes with exactly 1 NaN are mean-imputed (mean across all samples for that
      gene/position).
    - Each sample vector is zero-padded to 35,000 elements; mask=0 for padded positions.
    """
    from collections import defaultdict

    print(f"Loading: {csv_path}")
    output_df = pd.read_csv(csv_path, index_col=0)
    samples = output_df.columns
    genes   = output_df.index

    # --- Cell 101: build position_map ---
    sample_vectors = {s: [] for s in samples}
    position_map   = []

    for gene in genes:
        for sample in samples:
            raw = output_df.at[gene, sample]
            if pd.isna(raw):
                continue
            if isinstance(raw, str):
                parts  = [v.strip() for v in raw.strip("[]").split(",")]
                values = [float(v) if v.lower() != "nan" else np.nan for v in parts]
            else:
                values = [float(v) if pd.notna(v) else np.nan for v in raw]

            if len(values) != 7:
                continue
            nan_count = sum(np.isnan(v) for v in values)
            if nan_count <= 1:
                for i, val in enumerate(values):
                    if not np.isnan(val):
                        sample_vectors[sample].append(val)
                        position_map.append((gene, sample, i, val))

    sample_vectors = {k: v for k, v in sample_vectors.items() if len(v) > 0}

    # --- Cell 125: mean imputation + padding ---
    celltype_position_values = defaultdict(list)
    for gene, sample, pos, val in position_map:
        celltype_position_values[(gene, pos)].append(val)
    celltype_position_means = {
        k: np.mean(v) for k, v in celltype_position_values.items()
    }

    sample_vectors_padded = {}
    sample_vectors_mask   = {}

    for sample in sample_vectors:
        vector = []
        mask   = []

        for gene in genes:
            raw = output_df.at[gene, sample]
            if pd.isna(raw):
                continue
            if isinstance(raw, str):
                parts  = [v.strip() for v in raw.strip("[]").split(",")]
                values = [float(v) if v.lower() != "nan" else np.nan for v in parts]
            else:
                values = [float(v) if pd.notna(v) else np.nan for v in raw]

            if len(values) != 7:
                continue
            nan_count = int(np.isnan(values).sum())
            if nan_count == 0:
                vector.extend(values)
                mask.extend([1] * 7)
            elif nan_count == 1:
                nan_idx            = int(np.where(np.isnan(values))[0][0])
                values[nan_idx]    = celltype_position_means.get((gene, nan_idx), 0.0)
                vector.extend(values)
                mask.extend([1 if not np.isnan(v) else 0 for v in values])
            # genes with >1 NaN are skipped

        # zero-pad to 35,000
        if len(vector) < 35000:
            diff   = 35000 - len(vector)
            vector = vector + [0] * diff
            mask   = mask   + [0] * diff

        sample_vectors_padded[sample] = vector
        sample_vectors_mask[sample]   = np.array(mask, dtype=int)

    print(f"Loaded {len(sample_vectors_padded)} samples, "
          f"vector length {len(next(iter(sample_vectors_padded.values())))}")
    return sample_vectors_padded, sample_vectors_mask


# =====================================================================
# MAIN
# =====================================================================
def run_pipeline(sample_vectors_padded, sample_vectors_mask):
    """
    Run the full pipeline.

    Parameters
    ----------
    sample_vectors_padded : dict[str, np.ndarray]
        sample_id -> padded feature vector
    sample_vectors_mask : dict[str, np.ndarray]
        sample_id -> binary mask (1 = valid, 0 = padded)
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    sample_ids = list(sample_vectors_padded.keys())

    X = np.array([sample_vectors_padded[s] for s in sample_ids], dtype=np.float32)
    mask_array = np.array([sample_vectors_mask[s] for s in sample_ids], dtype=np.float32)

    # Step 1: VAE + GMM
    all_labels = train_all_seeds(X, mask_array, sample_ids)

    # Step 2: consensus -> final clusters
    final_df, final0, M = consensus_clustering(all_labels)

    # Step 3: Sankey diagrams
    build_sankeys(all_labels, final0)

    print(f"\nPipeline complete. All outputs in: {OUTPUT_DIR}")
    return final_df


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="input_dnn_celltypes_pseudobulked_mapped_VectorsofFloats.csv",
        help="Path to VectorsofFloats CSV (output of 02_build_feature_vectors.py)"
    )
    args = parser.parse_args()

    sample_vectors_padded, sample_vectors_mask = load_data(args.input)
    run_pipeline(sample_vectors_padded, sample_vectors_mask)
