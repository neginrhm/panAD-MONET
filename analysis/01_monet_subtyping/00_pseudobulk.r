###############################################################################
# Pseudobulk Matrix Generation — all 7 major cell types
#
# Logic: identical to all_cells_pseudobulking.r (scripts/ and fig7/panAD/utils/)
# Extended here to include OPC and ENDO, which were missing from the original.
#
# Input:  per-cell-type Seurat RDS files (RNA assay, raw counts layer)
# Output: pseudobulk_{CellType}.csv  — genes x samples, raw summed counts
#         (these are the files read by 01_regress_covariates.r / regree_out.r)
#
# NOTE: RDS paths below are HPC paths inferred from the existing 5-cell-type
#       script. Verify OPC and ENDO paths before running.
###############################################################################

library(Seurat)
library(dplyr)

create_pseudobulk <- function(rds_file, output_file) {
  obj <- readRDS(rds_file)

  counts_matrix  <- GetAssayData(obj, assay = "RNA", layer = "counts")
  sample_metadata <- obj@meta.data

  unique_samples <- unique(sample_metadata$Sample)
  pseudobulk_matrix <- matrix(
    0,
    nrow = nrow(counts_matrix),
    ncol = length(unique_samples)
  )
  rownames(pseudobulk_matrix) <- rownames(counts_matrix)
  colnames(pseudobulk_matrix) <- unique_samples

  for (sample in unique_samples) {
    cell_idx <- which(sample_metadata$Sample == sample)
    pseudobulk_matrix[, sample] <- rowSums(
      counts_matrix[, cell_idx, drop = FALSE]
    )
  }

  write.csv(pseudobulk_matrix, output_file)
  cat("Saved:", output_file, "\n")
}

# --- 5 original cell types (paths confirmed in all_cells_pseudobulking.r) ---

create_pseudobulk(
  "./Microglia/12ad_microglia_states.rds",
  "pseudobulk_Microglia.csv"
)
create_pseudobulk(
  "./Astrocyte/12ad_astrocye_hdWGCNA_object.rds",
  "pseudobulk_Astrocyte.csv"
)
create_pseudobulk(
  "./ODC/12ad_odc_hdWGCNA_object.rds",
  "pseudobulk_Oligodendrocyte.csv"
)
create_pseudobulk(
  "./Excitatory/12ad_excitatory_hdWGCNA_object.rds",
  "pseudobulk_Excitatory.csv"
)
create_pseudobulk(
  "./Inhibitory/12ad_INH_hdWGCNA_object.rds",
  "pseudobulk_Inhibitory.csv"
)

# --- 2 additional cell types (paths inferred — verify before running) ---

create_pseudobulk(
  "./OPC/12ad_opc_hdWGCNA_object.rds",
  "pseudobulk_OPC.csv"
)
create_pseudobulk(
  "./ENDO/12ad_endo_hdWGCNA_object.rds",
  "pseudobulk_ENDO.csv"
)
