# ============================================================
# signed_cor_hdwgcna_modules_celltypes.r
# ============================================================
# Donor-level module-trait Spearman correlations using
# pseudobulked (median) module eigengenes across seven cell types.
# Results are BH-FDR corrected and plotted as faceted bar charts
# per trait.
#
# Inputs
# ------
# - hdWGCNA Seurat objects per cell type (RDS files)
# - Donor metadata CSV (rosmap_pathology_cluster_ctrl_mci_ad.csv)
#
# Outputs
# -------
# - correlation_summary_DONORLEVEL.csv
# - DONOR_trait_{trait}_module_trait_spearman.pdf (per trait)
# - DONOR_ALL_TRAITS_panel_module_trait_spearman.pdf (combined)
#
# Usage
# -----
#   Rscript signed_cor_hdwgcna_modules_celltypes.r
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(cowplot)
  library(tidytext)
  library(rlang)
  library(glue)
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(scales)
})

theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 8)
set.seed(12345)

# --- Configuration -----------------------------------------------------------

BASE_DIR <- "./data"

celltype_paths <- c(
  Microglia   = file.path(BASE_DIR, "CellTypes/Microglia/hdWGCNA/13ad_microglia_hdWGCNA_object.rds"),
  Inhibitory  = file.path(BASE_DIR, "CellTypes/Inhibitory/13ad_INH_hdWGCNA_object.rds"),
  Excitatory  = file.path(BASE_DIR, "CellTypes/Excitatory/13ad_excitatory_post_dme_enrichr.rds"),
  ODC         = file.path(BASE_DIR, "CellTypes/ODC/13ad_odc_hdWGCNA_object.rds"),
  Astrocyte   = file.path(BASE_DIR, "CellTypes/Astrocyte/astrocyte_moduleUmap_hdwgcna.rds"),
  OPC         = file.path(BASE_DIR, "CellTypes/opc/13ad_opcc_hdWGCNA_object.rds"),
  Endothelial = file.path(BASE_DIR, "CellTypes/Endo/13ad_endo_hdWGCNA_object.rds")
)

meta <- read.csv(file.path(BASE_DIR, "csv/rosmap_pathology_cluster_ctrl_mci_ad.csv"),
                 stringsAsFactors = FALSE)

traits_wanted <- c("Braak", "CERAD", "delta_mmse", "Sex_numeric")

trait_colors <- c(
  "delta_mmse"  = "#8DA0CB",
  "Braak"       = "#4C72B0",
  "CERAD"       = "#66C2A5",
  "Sex_numeric" = "#C27BA0"
)

alpha_adj <- 0.05

# --- Helper functions --------------------------------------------------------

find_best_id_columns <- function(seurat_meta, meta_df) {
  seurat_candidates <- c("Sample", "sample", "SampleID", "sample_id", "orig.ident",
                         "donor", "Donor", "Subject", "subject", "individual",
                         "Individual", "projid", "proj_id")
  meta_candidates   <- c("sample_id", "Sample", "sample", "SampleID", "orig.ident",
                         "donor", "Donor", "Subject", "subject", "individual",
                         "Individual", "projid", "proj_id")
  seurat_cands <- intersect(seurat_candidates, colnames(seurat_meta))
  meta_cands   <- intersect(meta_candidates,   colnames(meta_df))

  best <- NULL; best_n <- 0
  for (sc in seurat_cands) {
    sv <- unique(as.character(seurat_meta[[sc]]))
    for (mc in meta_cands) {
      mv <- unique(as.character(meta_df[[mc]]))
      n <- length(intersect(sv, mv))
      if (n > best_n) best <- list(seurat_col = sc, meta_col = mc, overlap = n)
    }
  }
  if (is.null(best)) best <- list(seurat_col = NA_character_, meta_col = NA_character_, overlap = 0)
  best
}

harmonize_ids <- function(x) {
  x <- as.character(x)
  x <- sub("^[^_]+_", "", x)
  x <- sub("-[0-9]+$", "", x)
  x
}

# --- Donor-level correlations for one cell type ------------------------------

compute_cor_df_donor <- function(ct_label, rds_path, meta_df,
                                 traits_keep = traits_wanted,
                                 agg_fun = c("mean", "median")) {

  agg_fun <- match.arg(agg_fun)
  seurat_obj <- readRDS(rds_path)

  # Ensure Sex_numeric exists
  if (!"Sex_numeric" %in% colnames(seurat_obj@meta.data) && "Sex" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$Sex_numeric <- ifelse(seurat_obj@meta.data$Sex == "Male", 1, 0)
  }

  # Ensure delta_mmse exists in meta
  if (!"delta_mmse" %in% colnames(meta_df)) {
    if (all(c("cts_mmse30_lv", "cts_mmse30_first_ad_dx") %in% colnames(meta_df))) {
      meta_df$delta_mmse <- meta_df$cts_mmse30_lv - meta_df$cts_mmse30_first_ad_dx
    } else {
      meta_df$delta_mmse <- NA_real_
    }
  }

  # Detect join keys
  id_match <- find_best_id_columns(seurat_obj@meta.data, meta_df)
  if (is.na(id_match$seurat_col) || is.na(id_match$meta_col) || id_match$overlap == 0) {
    message(glue("[{ct_label}] No sample-ID overlap."))
    return(NULL)
  }
  seurat_key <- id_match$seurat_col
  meta_key   <- id_match$meta_col

  # Per-cell donor key
  cell_meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(DonorID = as.character(.data[[seurat_key]])) %>%
    select(cell, DonorID)

  # Donor traits (one row per donor)
  traits_keep <- traits_keep[traits_keep %in% colnames(meta_df)]
  if (!length(traits_keep)) return(NULL)

  meta_keyed <- meta_df %>%
    mutate(DonorID = as.character(.data[[meta_key]])) %>%
    distinct(DonorID, .keep_all = TRUE) %>%
    select(DonorID, all_of(traits_keep))

  cell_meta_full <- inner_join(cell_meta, meta_keyed, by = "DonorID")
  if (nrow(cell_meta_full) == 0) return(NULL)

  # Per-cell module eigengenes
  me <- GetMEs(seurat_obj, harmonized = FALSE)
  if (is.null(me) || nrow(me) == 0) {
    message(glue("[{ct_label}] No module eigengenes found."))
    return(NULL)
  }
  me <- as.data.frame(me)
  me$cell <- rownames(me)

  # Keep variable ME columns (drop grey)
  me_cols <- setdiff(colnames(me), c("cell", "grey"))
  if (!length(me_cols)) return(NULL)

  keep <- vapply(me_cols, function(cl) {
    v <- me[[cl]]
    any(!is.na(v)) && sd(v, na.rm = TRUE) > 1e-12
  }, logical(1))
  me_keep <- me_cols[keep]
  if (!length(me_keep)) return(NULL)

  me <- me[, c("cell", me_keep), drop = FALSE]

  # Harmonize cell IDs if needed
  if (length(intersect(me$cell, cell_meta_full$cell)) == 0) {
    me$cell <- harmonize_ids(me$cell)
    cell_meta_full$cell <- harmonize_ids(cell_meta_full$cell)
  }

  merged <- inner_join(me, cell_meta_full, by = "cell")
  if (nrow(merged) == 0) return(NULL)

  # Pseudobulk: aggregate within donor
  summarizer <- if (agg_fun == "mean") ~ mean(.x, na.rm = TRUE) else ~ median(.x, na.rm = TRUE)

  donor_me <- merged %>%
    group_by(DonorID) %>%
    summarise(across(all_of(me_keep), summarizer), .groups = "drop")

  donor_traits <- merged %>%
    select(DonorID, all_of(traits_keep)) %>%
    distinct(DonorID, .keep_all = TRUE)

  donor_df <- inner_join(donor_me, donor_traits, by = "DonorID")
  if (nrow(donor_df) < 3) return(NULL)

  # Spearman correlations across donors
  res <- lapply(me_keep, function(mod) {
    lapply(traits_keep, function(tr) {
      x <- donor_df[[mod]]
      y <- donor_df[[tr]]
      ok <- stats::complete.cases(x, y)
      if (sum(ok) < 3) return(NULL)
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      data.frame(
        celltype = ct_label,
        module   = mod,
        trait    = tr,
        rho      = unname(ct$estimate),
        pval     = ct$p.value,
        n_used   = sum(ok),
        agg      = agg_fun,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  res
}

# --- Run over all cell types -------------------------------------------------

res_list <- lapply(names(celltype_paths), function(ct) {
  compute_cor_df_donor(ct, celltype_paths[[ct]], meta, agg_fun = "median")
})
res_list <- Filter(Negate(is.null), res_list)
if (!length(res_list)) stop("No donor-level correlations computed.")

all_cor <- bind_rows(res_list)

# Multiple-testing correction
all_cor <- all_cor %>%
  mutate(
    pval_adj = p.adjust(pval, method = "BH"),
    R2       = as.numeric(rho)^2
  )

# --- Save results ------------------------------------------------------------

write.csv(all_cor, "correlation_summary_DONORLEVEL.csv", row.names = FALSE)
message("Wrote: correlation_summary_DONORLEVEL.csv")

# --- Plotting: bar charts per trait ------------------------------------------

traits_to_plot <- intersect(traits_wanted, unique(all_cor$trait))
trait_colors <- trait_colors[names(trait_colors) %in% traits_to_plot]

# Only plot significant results
plot_df <- all_cor %>% filter(!is.na(pval_adj) & pval_adj < alpha_adj)

if (nrow(plot_df) == 0) {
  warning(glue("No donor-level results pass FDR < {alpha_adj}. Plotting top 10 by |rho| per trait x celltype instead."))
  plot_df <- all_cor %>%
    group_by(trait, celltype) %>%
    slice_max(order_by = abs(rho), n = 10, with_ties = FALSE) %>%
    ungroup()
}

present_pairs <- plot_df %>% distinct(celltype, module)

label_cutoff <- 0.12
effect_lines_df <- data.frame(
  y  = c(-0.5, -0.3, 0, 0.3, 0.5),
  lt = c("dashed", "dashed", "solid", "dashed", "dashed"),
  lw = c(0.4, 0.4, 0.7, 0.4, 0.4),
  col = c("grey70", "grey70", "grey35", "grey70", "grey70")
)

make_trait_plot <- function(tr) {
  df_t <- plot_df %>%
    filter(trait == tr) %>%
    semi_join(present_pairs, by = c("celltype", "module")) %>%
    group_by(celltype, module) %>%
    summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(rho)) %>%
    group_by(celltype) %>%
    mutate(
      module_in_ct = tidytext::reorder_within(module, abs(rho), celltype, .desc = TRUE),
      hjust_lbl    = ifelse(rho >= 0, -0.15, 1.15)
    ) %>%
    ungroup() %>%
    mutate(module_in_ct = forcats::fct_drop(module_in_ct))

  ggplot(df_t, aes(x = module_in_ct, y = rho)) +
    geom_col(width = 0.7, fill = trait_colors[tr] %||% "grey70") +
    geom_segment(
      data = df_t,
      aes(x = module_in_ct, xend = module_in_ct, y = 0, yend = rho),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey20",
      linewidth = 0.6,
      lineend = "round",
      alpha = 0.9
    ) +
    geom_hline(
      data = effect_lines_df,
      aes(yintercept = y, linetype = lt, linewidth = lw, colour = col),
      inherit.aes = FALSE
    ) +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    scale_colour_identity() +
    geom_text(
      data = dplyr::filter(df_t, abs(rho) >= label_cutoff),
      aes(label = scales::number(rho, accuracy = 0.01), hjust = hjust_lbl),
      size = 3
    ) +
    facet_grid(celltype ~ ., scales = "free_y", space = "free_y") +
    tidytext::scale_x_reordered() +
    scale_x_discrete(drop = TRUE) +
    scale_y_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, by = 0.25),
                       minor_breaks = seq(-1, 1, by = 0.05)) +
    labs(
      title = paste0("DONOR-level Spearman correlations: ", tr),
      subtitle = paste0("Module eigengenes pseudobulked per donor within each cell type; BH-FDR<", alpha_adj, " for plotted set."),
      x = "Module",
      y = "Spearman \u03C1 (across donors)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      strip.text.y = element_text(face = "bold"),
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 9),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_flip(clip = "off")
}

plots_by_trait <- lapply(traits_to_plot, make_trait_plot)
names(plots_by_trait) <- traits_to_plot

# --- Save individual trait plots ---------------------------------------------

for (tr in traits_to_plot) {
  ggsave(
    filename = paste0("DONOR_trait_", tr, "_module_trait_spearman.pdf"),
    plot = plots_by_trait[[tr]],
    width = 10, height = 14, dpi = 300, bg = "white"
  )
  message("Saved: ", paste0("DONOR_trait_", tr, "_module_trait_spearman.pdf"))
}

# --- Save combined panel -----------------------------------------------------

if (length(plots_by_trait) > 0) {
  combined <- cowplot::plot_grid(plotlist = plots_by_trait[traits_to_plot], ncol = 2, align = "hv")
  ggsave(
    filename = "DONOR_ALL_TRAITS_panel_module_trait_spearman.pdf",
    plot = combined, width = 18, height = 14, dpi = 300, bg = "white"
  )
  message("Saved: DONOR_ALL_TRAITS_panel_module_trait_spearman.pdf")
}
