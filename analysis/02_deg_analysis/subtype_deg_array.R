#!/usr/bin/env Rscript
# ==============================================================================
# subtype_deg_array.R
# ==============================================================================
# Description:
#   SLURM array-compatible version of subtype DEG analysis. Processes a single
#   .h5ad file passed as a command-line argument and runs MAST differential
#   expression: Dementia vs Control, MCI vs Control, Dementia vs MCI.
#
# Inputs:
#   - A single .h5ad file (command-line argument) with 'dx' column in obs
#
# Outputs:
#   - MAST DEG CSVs in the working directory:
#     {basename}_mast_dementia_vs_control.csv
#     {basename}_mast_mci_vs_control.csv
#     {basename}_mast_dementia_vs_mci.csv
#
# Usage:
#   Rscript subtype_deg_array.R <path_to_h5ad_file>
# ==============================================================================

# --- Parse arguments ---

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

cat("Processing:", file, "\n")

# --- Load libraries ---

library(Seurat)
library(anndata)
library(SeuratDisk)
library(reticulate)
library(dplyr)

set.seed(12345)

# --- Load data ---

adata <- read_h5ad(file)

# Create Seurat object (transpose needed)
seurat_obj <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
seurat_obj <- NormalizeData(seurat_obj)

# Check if 'dx' column exists
if (!"dx" %in% colnames(seurat_obj@meta.data)) {
  stop("No 'dx' column found in metadata.")
}

# Set identity class
Idents(seurat_obj) <- seurat_obj$dx
levels_present <- levels(Idents(seurat_obj))

# Output file name prefix
out_prefix <- tools::file_path_sans_ext(basename(file))

# --- Run MAST comparisons (only if groups exist) ---

if (all(c("Dementia", "Control") %in% levels_present)) {
  cat("Running Dementia vs Control\n")
  dementia_vs_control <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  write.csv(dementia_vs_control, file = paste0(out_prefix, "_mast_dementia_vs_control.csv"))
}

if (all(c("MCI", "Control") %in% levels_present)) {
  cat("Running MCI vs Control\n")
  mci_vs_control <- FindMarkers(seurat_obj, ident.1 = "MCI", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  write.csv(mci_vs_control, file = paste0(out_prefix, "_mast_mci_vs_control.csv"))
}

if (all(c("Dementia", "MCI") %in% levels_present)) {
  cat("Running Dementia vs MCI\n")
  dementia_vs_mci <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "MCI", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  write.csv(dementia_vs_mci, file = paste0(out_prefix, "_mast_dementia_vs_mci.csv"))
}

cat("Finished processing", file, "\n")
