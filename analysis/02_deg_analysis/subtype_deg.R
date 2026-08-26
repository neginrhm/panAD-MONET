# ==============================================================================
# subtype_deg.R
# ==============================================================================
# Description:
#   Loops over all .h5ad files in a directory and runs MAST differential
#   expression for each cell subtype: Dementia vs Control, MCI vs Control,
#   Dementia vs MCI. Skips files without a 'dx' column.
#
# Inputs:
#   - Directory of .h5ad files, each containing a cell subtype with 'dx' in obs
#
# Outputs:
#   - Per-subtype MAST DEG CSVs in output_dir:
#     {basename}_mast_dementia_vs_control.csv
#     {basename}_mast_mci_vs_control.csv
#     {basename}_mast_dementia_vs_mci.csv
#
# Usage:
#   Rscript subtype_deg.R
# ==============================================================================

# --- Load libraries ---

library(Seurat)
library(SeuratDisk)
library(anndata)
library(reticulate)
library(dplyr)

# --- Configuration ---

theme_set(ggplot2::theme_bw())
set.seed(12345)

h5ad_dir <- "./"
output_dir <- "differential_expression_results"
dir.create(output_dir, showWarnings = FALSE)

# --- Process each h5ad file ---

h5ad_files <- list.files(h5ad_dir, pattern = "\\.h5ad$", full.names = TRUE)

for (file in h5ad_files) {
  cat("Processing", file, "\n")

  # Load .h5ad and convert to Seurat
  ad <- read_h5ad(file)
  seurat_obj <- CreateSeuratObject(counts = t(ad$X), meta.data = ad$obs)
  seurat_obj <- NormalizeData(seurat_obj)

  # Set identities
  if (!"dx" %in% colnames(seurat_obj@meta.data)) {
    cat("Skipping", file, "- 'dx' column not found\n")
    next
  }

  Idents(seurat_obj) <- seurat_obj$dx
  levels_present <- levels(Idents(seurat_obj))
  base <- tools::file_path_sans_ext(basename(file))

  # --- Run MAST comparisons (only if both groups exist) ---

  if (all(c("Dementia", "Control") %in% levels_present)) {
    res <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
    write.csv(res, file = file.path(output_dir, paste0(base, "_mast_dementia_vs_control.csv")))
  }

  if (all(c("MCI", "Control") %in% levels_present)) {
    res <- FindMarkers(seurat_obj, ident.1 = "MCI", ident.2 = "Control", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
    write.csv(res, file = file.path(output_dir, paste0(base, "_mast_mci_vs_control.csv")))
  }

  if (all(c("Dementia", "MCI") %in% levels_present)) {
    res <- FindMarkers(seurat_obj, ident.1 = "Dementia", ident.2 = "MCI", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
    write.csv(res, file = file.path(output_dir, paste0(base, "_mast_dementia_vs_mci.csv")))
  }
}
