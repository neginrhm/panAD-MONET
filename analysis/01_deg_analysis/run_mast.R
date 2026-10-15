#!/usr/bin/env Rscript
# ==============================================================================
# run_mast.R
# ==============================================================================
# Description:
#   Runs MAST differential expression analysis on a single cell-type h5ad file.
#   Performs three comparisons: AD vs Control, MCI vs Control, AD vs MCI.
#
# Inputs:
#   - A single .h5ad file containing cells of one cell type with a 'dx' column
#     in obs (e.g., unique_ons_all_Celltype_Astrocyte_385040cells.h5ad)
#
# Outputs:
#   - Three CSV files per cell type with MAST DEG results:
#     {celltype}_mast_AD_vs_Control.csv
#     {celltype}_mast_MCI_vs_Control.csv
#     {celltype}_mast_AD_vs_MCI.csv
#
# Usage:
#   Rscript run_mast.R unique_ons_all_Celltype_Astrocyte_385040cells.h5ad
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
infile <- args[[1]]

# --- Load libraries ---

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(reticulate)
  library(Matrix)
})

# Thread configuration
nthreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
set.seed(12345)

# Reticulate Python configuration (optional)
if (nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)
}

# --- Load data ---

message(sprintf("Reading: %s", infile))
data <- read_h5ad(infile)

# --- Create Seurat object and normalize ---

data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
seurat_obj <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(seurat_obj) <- seurat_obj$dx

# --- Run MAST differential expression ---
# min.pct=0 matches Methods ("Differential expression analysis") and every other MAST call in the study
# (subtype/sex/hdWGCNA-embedded DE); an earlier draft of this script used 0.1.

ad_vs_control  <- FindMarkers(seurat_obj, ident.1 = "AD",  ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
mci_vs_control <- FindMarkers(seurat_obj, ident.1 = "MCI", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
ad_vs_mci      <- FindMarkers(seurat_obj, ident.1 = "AD",  ident.2 = "MCI",     test.use = "MAST", min.pct = 0, logfc.threshold = 0)

# --- Save results ---

# Derive cell type from filename pattern: unique_ons_all_Celltype_{celltype}_{N}cells.h5ad
celltype <- sub("^unique_ons_all_Celltype_([^_]+)_.*$", "\\1", basename(infile))

write.csv(ad_vs_control,  sprintf("%s_mast_AD_vs_Control.csv",  celltype), quote = TRUE)
write.csv(mci_vs_control, sprintf("%s_mast_MCI_vs_Control.csv", celltype), quote = TRUE)
write.csv(ad_vs_mci,      sprintf("%s_mast_AD_vs_MCI.csv",      celltype), quote = TRUE)

message("Done.")
