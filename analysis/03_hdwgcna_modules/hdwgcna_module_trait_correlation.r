# ============================================================
# hdwgcna_module_trait_correlation.r
# ============================================================
# Compute module-trait correlations for Inhibitory neurons using
# hdWGCNA's built-in ModuleTraitCorrelation function.  Traits
# include clinical pathology scores (CERAD, Braak, MMSE) and
# technical covariates.
#
# Inputs
# ------
# - hdWGCNA Seurat object for Inhibitory neurons (RDS)
# - Donor metadata CSV (rosmap_pathology_cluster_ctrl_mci_ad.csv)
#
# Outputs
# -------
# - hwcgna_inhh_module_trait_corr2.pdf (heatmap of correlations)
#
# Usage
# -----
#   Rscript hdwgcna_module_trait_correlation.r
# ============================================================

# --- Load libraries ----------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(WGCNA)
library(hdWGCNA)

theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 8)
set.seed(12345)

# --- Configuration -----------------------------------------------------------

BASE_DIR <- "./data"

# --- Load data ---------------------------------------------------------------

seurat_obj <- readRDS(file.path(BASE_DIR, "CellTypes/Inhibitory/13ad_INH_hdWGCNA_object.rds"))

meta <- read.csv(file.path(BASE_DIR, "csv/rosmap_pathology_cluster_ctrl_mci_ad.csv"))

# --- Prepare metadata --------------------------------------------------------

# Convert metadata types
seurat_obj$Sex <- as.factor(seurat_obj$Sex)
seurat_obj$Study <- as.factor(seurat_obj$Study)
seurat_obj$Age <- as.numeric(seurat_obj$Age)
seurat_obj$CERAD <- as.numeric(seurat_obj$CERAD)
seurat_obj$Braak <- as.numeric(seurat_obj$Braak)

# Match sample IDs between Seurat object and donor metadata
meta$sample_id <- as.character(meta$sample_id)
seurat_obj$Sample <- as.character(seurat_obj$Sample)

# Transfer traits from donor metadata to Seurat object
traits_to_add <- c("CERAD", "Braak", "cts_mmse30_first_ad_dx", "cts_mmse30_lv")

for (trait in traits_to_add) {
  trait_vector <- setNames(meta[[trait]], meta$sample_id)
  values <- trait_vector[seurat_obj$Sample]

  if (length(values) != ncol(seurat_obj)) {
    stop(paste("Length mismatch for trait:", trait))
  }

  names(values) <- NULL
  seurat_obj[[trait]] <- values
}

# --- Analysis: module-trait correlation --------------------------------------

cur_traits <- c('cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'CERAD', 'Braak',
                'Sex', 'Age', 'doublet_score', 'nCount_RNA', 'nFeature_RNA',
                'total_counts_MT')

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by = 'label_transfer_annot'
)

mt_cor <- GetModuleTraitCorrelation(seurat_obj)

# --- Save results ------------------------------------------------------------

pdf("hwcgna_inhh_module_trait_corr2.pdf")
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = '#FFA07A',
  mid_color = '#F5F5F5',
  low_color = '#20B2AA',
  plot_max = 0.2,
  combine = TRUE
)
dev.off()
