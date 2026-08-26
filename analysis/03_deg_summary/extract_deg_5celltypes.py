"""
extract_deg_5celltypes.py
=========================
Extract and filter differentially expressed genes (DEGs) across five major
cell types (EXC, ODC, INH, AST, MIC), overlap them with hdWGCNA hub genes,
and produce summary heatmaps (single-condition and multi-condition).

Inputs
------
- DEG marker CSVs per cell type (MAST output), e.g.
      ./data/control_vs_dementia/mast_{celltype}_control_vs_dementia_markers.csv
- hdWGCNA hub-gene tables per cell type, e.g.
      ./data/celltype_hub_genes/12ad_{celltype}_modules_hub_df_table.csv

Outputs
-------
- deg_heatmap_across_celltypes.png       Single-condition heatmap (control vs dementia)
- comprehensive_deg_heatmap.png          Multi-condition heatmap (3 contrasts x 5 cell types)
- scaled_comprehensive_deg_clustermap_fc0.5.png   Z-score clustermap

Usage
-----
    python extract_deg_5celltypes.py
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

# --- Configuration -----------------------------------------------------------

BASE_DIR = "./data"

cell_types = ['exc', 'odc', 'inh', 'ast', 'mic']

deg_path_template = f"{BASE_DIR}/control_vs_dementia/mast_{{celltype}}_control_vs_dementia_markers.csv"
module_path_template = f"{BASE_DIR}/celltype_hub_genes/12ad_{{celltype}}_modules_hub_df_table.csv"

# --- Load data ---------------------------------------------------------------

deg_data = {}
module_data = {}

for celltype in cell_types:
    # Load DEG data
    deg_file = deg_path_template.format(celltype=celltype)
    try:
        deg_df = pd.read_csv(deg_file)
        deg_df.rename(columns={deg_df.columns[0]: 'gene_name'}, inplace=True)
        deg_data[celltype] = deg_df
        print(f"Successfully loaded and renamed DEG data for {celltype.upper()}")
    except FileNotFoundError:
        print(f"DEG file for {celltype.upper()} not found at {deg_file}")

    # Load module hub gene data
    module_file = module_path_template.format(celltype=celltype)
    try:
        module_data[celltype] = pd.read_csv(module_file)
        print(f"Successfully loaded module data for {celltype.upper()}")
    except FileNotFoundError:
        print(f"Module file for {celltype.upper()} not found at {module_file}")

# Verify loaded data
for celltype in cell_types:
    if celltype in deg_data:
        print(f"DEG data for {celltype.upper()} cell type:")
        print(deg_data[celltype].head(), "\n")

# --- Overlap: all DEGs vs hub genes ------------------------------------------

overlap_genes = defaultdict(set)

for celltype in deg_data:
    deg_genes = set(deg_data[celltype]['gene_name'])
    module_genes = set(module_data[celltype]['gene_name'])
    overlap_genes[celltype] = deg_genes.intersection(module_genes)

for celltype, genes in overlap_genes.items():
    print(f"Number of overlapping genes in {celltype.upper()}: {len(genes)}")

# --- Overlap: filtered DEGs vs hub genes -------------------------------------

p_val_threshold = 0.05
logfc_threshold = 0.25

overlap_genes = defaultdict(set)

for celltype in deg_data:
    filtered_deg = deg_data[celltype][
        (deg_data[celltype]['p_val_adj'] < p_val_threshold) &
        ((deg_data[celltype]['avg_log2FC'] > logfc_threshold) |
         (deg_data[celltype]['avg_log2FC'] < -logfc_threshold))
    ]
    deg_genes = set(filtered_deg['gene_name'])
    module_genes = set(module_data[celltype]['gene_name'])
    overlap_genes[celltype] = deg_genes.intersection(module_genes)

for celltype, genes in overlap_genes.items():
    print(f"Number of filtered overlapping genes in {celltype.upper()}: {len(genes)}")

# --- Heatmap: single condition (control vs dementia) -------------------------

deg_filtered_data = {}
all_degs = set()

for celltype in deg_data:
    filtered_deg = deg_data[celltype][
        (deg_data[celltype]['p_val_adj'] < p_val_threshold) &
        ((deg_data[celltype]['avg_log2FC'] > logfc_threshold) |
         (deg_data[celltype]['avg_log2FC'] < -logfc_threshold))
    ]
    deg_filtered_data[celltype] = filtered_deg[['gene_name', 'avg_log2FC']]
    all_degs.update(filtered_deg['gene_name'])

heatmap_data = pd.DataFrame(index=list(all_degs))

for celltype, filtered_df in deg_filtered_data.items():
    logfc_data = filtered_df.set_index('gene_name')['avg_log2FC']
    heatmap_data = heatmap_data.join(logfc_data.rename(celltype.upper()), how='left')

heatmap_data.fillna(0, inplace=True)

plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, cmap='coolwarm', annot=False,
            cbar_kws={'label': 'Log2 Fold Change'})
plt.title('Heatmap of Filtered DEGs (Control vs Dementia) Across Cell Types')
plt.xlabel('Cell Types')
plt.ylabel('DEG Genes')
plt.savefig('deg_heatmap_across_celltypes.png')
plt.show()

# --- Heatmap: multi-condition (3 contrasts x 5 cell types) -------------------

logfc_threshold_multi = 0.5
conditions = ['control_vs_dementia', 'control_vs_mci', 'mci_vs_dementia']
directories = {cond: f"{BASE_DIR}/{cond}/" for cond in conditions}

heatmap_data = pd.DataFrame()

for condition, dir_path in directories.items():
    for celltype in cell_types:
        file_path = f"{dir_path}mast_{celltype}_{condition}_markers.csv"
        try:
            data = pd.read_csv(file_path)
            data.rename(columns={data.columns[0]: 'gene_name'}, inplace=True)

            filtered_data = data[
                (data['p_val_adj'] < p_val_threshold) &
                (abs(data['avg_log2FC']) > logfc_threshold_multi)
            ]

            col_name = f"{condition}_{celltype}"
            filtered_data.set_index('gene_name', inplace=True)
            heatmap_data[col_name] = filtered_data['avg_log2FC']

        except FileNotFoundError:
            print(f"File not found: {file_path}")
            continue

heatmap_data.fillna(0, inplace=True)

scaled_data = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Multi-condition heatmap
plt.figure(figsize=(18, 10))
sns.heatmap(heatmap_data, cmap='coolwarm', annot=False,
            cbar_kws={'label': 'Log2 Fold Change'})
plt.title('Comprehensive Heatmap of DEGs Across Conditions and Cell Types')
plt.xlabel('Conditions and Cell Types')
plt.ylabel('DEG Genes')
plt.savefig('comprehensive_deg_heatmap.png')
plt.show()

# --- Clustermap: z-score scaled ---------------------------------------------

g = sns.clustermap(
    scaled_data, cmap='coolwarm', annot=False,
    cbar_kws={'label': 'Z-Score (Scaled Log2 Fold Change'},
    figsize=(18, 10), metric='euclidean', method='average'
)
plt.title('Comprehensive Clustermap of scaled DEGs Across Conditions and Cell Types fc0.5')
g.savefig('scaled_comprehensive_deg_clustermap_fc0.5.png')
plt.show()
