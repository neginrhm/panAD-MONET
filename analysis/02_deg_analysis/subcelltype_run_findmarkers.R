#!/usr/bin/env Rscript
# ==============================================================================
# subcelltype_run_findmarkers.R
# ==============================================================================
# Description:
#   SLURM array-compatible MAST DEG analysis for sub-cell-type h5ad files.
#   Takes a file path and cell type name as arguments. Uses pre-normalized data
#   in the 'data' slot for FindMarkers.
#
# Inputs:
#   - file_path: path to a .h5ad file for one sub-cell-type
#   - celltype:  name used as output file prefix
#
# Outputs:
#   - MAST DEG CSVs:
#     {celltype}_subcelltype_dementia_vs_control_markers.csv
#     {celltype}_subcelltype_mci_vs_control_markers.csv
#     {celltype}_subcelltype_dementia_vs_mci_markers.csv
#
# Usage:
#   Rscript subcelltype_run_findmarkers.R <file_path> <celltype>
# ==============================================================================

# --- Parse arguments ---

args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
celltype <- args[2]

# --- Load libraries ---

library(Seurat)
library(anndata)
library(reticulate)
library(dplyr)

set.seed(12345)

# --- Load data ---

adata <- read_h5ad(file_path)

# Convert AnnData to Seurat object
seurat_obj <- CreateSeuratObject(
  counts = t(adata$X),
  meta.data = as.data.frame(adata$obs)
)

# Use pre-normalized data in 'X' as 'data' slot
seurat_obj <- SetAssayData(
  object = seurat_obj,
  assay = "RNA",
  layer = "data",
  new.data = GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
)

# --- Run MAST differential expression ---

Idents(seurat_obj) <- seurat_obj$dx

dementia_vs_control <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
mci_vs_control      <- FindMarkers(seurat_obj, ident.1 = "MCI", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
dementia_vs_mci     <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "MCI", test.use = "MAST", min.pct = 0, logfc.threshold = 0)

# --- Save results ---

write.csv(dementia_vs_control, paste0(celltype, "_subcelltype_dementia_vs_control_markers.csv"))
write.csv(mci_vs_control, paste0(celltype, "_subcelltype_mci_vs_control_markers.csv"))
write.csv(dementia_vs_mci, paste0(celltype, "_subcelltype_dementia_vs_mci_markers.csv"))
