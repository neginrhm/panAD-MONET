# ============================================================
# signed_cor_hdwgcna_modules_LAYOUT3.r
# ============================================================
# Donor-level (statistically correct) module-trait correlations
# using pseudobulked module eigengenes (mean + median) with
# Spearman correlation across donors and BH-FDR correction.
#
# Generates grouped publication-ready figures (2x2 trait panels):
#   - Neurons:  Inhibitory + Excitatory
#   - Glia:     Microglia + ODC + Astrocyte
#   - VascOPC:  Endothelial + OPC
# Each figure shows all four traits (Braak, CERAD, delta_mmse,
# Sex_numeric) with shared x-axis scale.
#
# Inputs
# ------
# - hdWGCNA Seurat objects per cell type (RDS files)
# - Donor metadata CSV (rosmap_pathology_cluster_ctrl_mci_ad.csv)
#
# Outputs
# -------
# - correlation_DONORLEVEL_all_mean+median.csv
# - correlation_DONORLEVEL_sigONLY_mean+median.csv
# - donor_pseudobulk_MEs_{celltype}_{agg}.csv  (per cell type)
# - ALL4TRAITS_{group}_DONORLEVEL_{agg}_SIGONLY.pdf (per group x agg)
#
# Usage
# -----
#   Rscript signed_cor_hdwgcna_modules_LAYOUT3.r
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(glue)
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(scales)
  library(tibble)
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
aggs_to_run   <- c("mean", "median")

alpha_adj  <- 0.05
min_donors <- 10
PAD_TO_N   <- 12

trait_colors <- c(
  "delta_mmse"  = "#8DA0CB",
  "Braak"       = "#4C72B0",
  "CERAD"       = "#66C2A5",
  "Sex_numeric" = "#C27BA0"
)

celltype_groups <- list(
  Neurons = c("Inhibitory", "Excitatory"),
  Glia    = c("Microglia", "ODC", "Astrocyte"),
  VascOPC = c("Endothelial", "OPC")
)

# --- Ensure Sex_numeric exists in donor metadata -----------------------------

if (!"Sex_numeric" %in% colnames(meta)) {
  sex_cols <- intersect(c("Sex", "sex", "gender", "Gender"), colnames(meta))
  if (length(sex_cols) > 0) {
    sc <- sex_cols[1]
    meta$Sex_numeric <- dplyr::case_when(
      tolower(as.character(meta[[sc]])) %in% c("male", "m")   ~ 1,
      tolower(as.character(meta[[sc]])) %in% c("female", "f") ~ 0,
      TRUE ~ NA_real_
    )
  } else {
    meta$Sex_numeric <- NA_real_
  }
}

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
      if (n > best_n) {
        best_n <- n
        best <- list(seurat_col = sc, meta_col = mc, overlap = n)
      }
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

tight_lim <- function(max_abs) {
  if (!is.finite(max_abs) || is.na(max_abs)) return(0.20)
  lim <- max_abs * 1.08
  lim <- max(lim, 0.06)
  lim <- min(lim, 0.40)
  ceiling(lim * 100) / 100
}

# --- Pad and center modules for uniform bar thickness ------------------------

pad_center_one_celltype <- function(df_ct, pad_to = 12) {
  df_ct <- df_ct %>%
    mutate(module = as.character(module), abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho), desc(rho), module)

  k <- nrow(df_ct)
  if (k == 0) return(NULL)

  if (k > pad_to) {
    df_ct <- df_ct %>% slice_head(n = pad_to)
    k <- pad_to
  }

  n_pad   <- pad_to - k
  n_above <- floor(n_pad / 2)
  n_below <- n_pad - n_above

  pads_above <- tibble(
    module  = paste0("__PAD__A__", seq_len(n_above)),
    rho     = 0,
    trait   = df_ct$trait[1],
    .is_pad = TRUE
  )
  pads_below <- tibble(
    module  = paste0("__PAD__B__", seq_len(n_below)),
    rho     = 0,
    trait   = df_ct$trait[1],
    .is_pad = TRUE
  )

  df_ct <- df_ct %>% mutate(.is_pad = FALSE)

  d_all <- bind_rows(pads_above, df_ct, pads_below)
  level_order <- c(pads_above$module, df_ct$module, pads_below$module)

  d_all <- d_all %>%
    mutate(
      module_plot = factor(module, levels = level_order),
      label_rho   = ifelse(.is_pad, "", scales::number(rho, accuracy = 0.01)),
      hjust_rho   = ifelse(rho >= 0, -0.08, 1.08)
    )

  real_breaks <- d_all %>% filter(!.is_pad) %>% pull(module_plot)
  list(data = d_all, breaks = real_breaks)
}

# --- Right-side celltype strip label -----------------------------------------

strip_right <- function(label) {
  ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
              fill = "white", color = "black", linewidth = 0.6) +
    geom_text(aes(x = 0.5, y = 0.5, label = label),
              angle = 90, fontface = "bold", size = 4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# --- Core: donor-level correlations ------------------------------------------

compute_cor_df_donor <- function(ct_label, rds_path, meta_df,
                                 traits_keep = traits_wanted,
                                 agg_fun = c("mean", "median")) {

  agg_fun <- match.arg(agg_fun)
  seurat_obj <- readRDS(rds_path)

  # Create delta_mmse in donor meta if needed
  if (!"delta_mmse" %in% colnames(meta_df)) {
    if (all(c("cts_mmse30_lv", "cts_mmse30_first_ad_dx") %in% colnames(meta_df))) {
      meta_df$delta_mmse <- meta_df$cts_mmse30_lv - meta_df$cts_mmse30_first_ad_dx
    } else {
      meta_df$delta_mmse <- NA_real_
    }
  }

  id_match <- find_best_id_columns(seurat_obj@meta.data, meta_df)
  if (is.na(id_match$seurat_col) || is.na(id_match$meta_col) || id_match$overlap == 0) {
    message(glue("[{ct_label}] No sample-ID overlap."))
    return(NULL)
  }
  seurat_key <- id_match$seurat_col
  meta_key   <- id_match$meta_col

  cell_meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(DonorID = as.character(.data[[seurat_key]])) %>%
    select(cell, DonorID)

  traits_keep <- traits_keep[traits_keep %in% colnames(meta_df)]
  if (!length(traits_keep)) return(NULL)

  donor_traits <- meta_df %>%
    mutate(DonorID = as.character(.data[[meta_key]])) %>%
    distinct(DonorID, .keep_all = TRUE) %>%
    select(DonorID, all_of(traits_keep))

  cell_meta_full <- inner_join(cell_meta, donor_traits, by = "DonorID")
  if (nrow(cell_meta_full) == 0) return(NULL)

  me <- GetMEs(seurat_obj, harmonized = FALSE)
  if (is.null(me) || nrow(me) == 0) return(NULL)

  me <- as.data.frame(me)
  me$cell <- rownames(me)

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

  summarizer <- if (agg_fun == "mean") ~ mean(.x, na.rm = TRUE) else ~ median(.x, na.rm = TRUE)

  donor_me <- merged %>%
    group_by(DonorID) %>%
    summarise(across(all_of(me_keep), summarizer), .groups = "drop")

  donor_df <- inner_join(donor_me, donor_traits, by = "DonorID")
  if (nrow(donor_df) < min_donors) return(NULL)

  write.csv(donor_df,
            paste0("donor_pseudobulk_MEs_", ct_label, "_", agg_fun, ".csv"),
            row.names = FALSE)

  res <- lapply(me_keep, function(mod) {
    lapply(traits_keep, function(tr) {
      x <- donor_df[[mod]]
      y <- donor_df[[tr]]
      ok <- stats::complete.cases(x, y)
      n  <- sum(ok)
      if (n < min_donors) return(NULL)
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      data.frame(
        celltype = ct_label,
        module   = mod,
        trait    = tr,
        agg      = agg_fun,
        rho      = unname(ct$estimate),
        pval     = ct$p.value,
        n_donors = n,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  res
}

# --- Run correlations --------------------------------------------------------

all_cor_list <- list()
for (agg in aggs_to_run) {
  res_list <- lapply(names(celltype_paths), function(ct) {
    compute_cor_df_donor(ct, celltype_paths[[ct]], meta, agg_fun = agg)
  })
  res_list <- Filter(Negate(is.null), res_list)
  if (!length(res_list)) stop("No donor-level correlations computed for agg = ", agg)
  all_cor_list[[agg]] <- bind_rows(res_list)
}

all_cor <- bind_rows(all_cor_list) %>%
  group_by(celltype, trait, agg) %>%
  mutate(pval_adj = p.adjust(pval, method = "BH")) %>%
  ungroup()

# --- Save results ------------------------------------------------------------

write.csv(all_cor, "correlation_DONORLEVEL_all_mean+median.csv", row.names = FALSE)

sig_only <- all_cor %>% filter(!is.na(pval_adj) & pval_adj < alpha_adj)
write.csv(sig_only, "correlation_DONORLEVEL_sigONLY_mean+median.csv", row.names = FALSE)

message("Wrote: correlation_DONORLEVEL_all_mean+median.csv")
message("Wrote: correlation_DONORLEVEL_sigONLY_mean+median.csv")

# --- Plotting: grouped figures with 2x2 trait panels ------------------------

make_trait_panel <- function(df_sig, tr, agg_choice, ct_order_full, lim, brks) {

  df_tr <- df_sig %>%
    filter(trait == tr, agg == agg_choice, celltype %in% ct_order_full) %>%
    mutate(celltype = as.character(celltype), module = as.character(module))

  if (nrow(df_tr) == 0) return(NULL)

  ct_keep <- intersect(ct_order_full, unique(df_tr$celltype))
  df_tr <- df_tr %>% filter(celltype %in% ct_keep)

  panels <- lapply(ct_keep, function(ct) {
    df_ct <- df_tr %>% filter(celltype == ct)
    pad_info <- pad_center_one_celltype(df_ct, pad_to = PAD_TO_N)
    if (is.null(pad_info)) return(NULL)

    d_all <- pad_info$data
    real_breaks <- pad_info$breaks

    p_main <- ggplot(d_all, aes(x = module_plot, y = rho)) +
      geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4) +
      geom_col(aes(fill = trait, alpha = !.is_pad), width = 0.60, color = NA) +
      geom_text(
        aes(label = label_rho, hjust = hjust_rho, alpha = !.is_pad),
        size = 3.0, color = "grey10", show.legend = FALSE
      ) +
      scale_x_discrete(drop = FALSE, breaks = real_breaks) +
      scale_y_continuous(limits = c(-lim, lim), breaks = brks) +
      scale_fill_manual(values = trait_colors, name = "Trait") +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") +
      coord_flip(clip = "off") +
      theme_classic(base_size = 11) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 6)
      )

    cowplot::plot_grid(
      p_main,
      strip_right(ct),
      nrow = 1,
      rel_widths = c(1, 0.075)
    )
  })

  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) return(NULL)

  # x tick labels only on bottom celltype in this trait-panel
  if (length(panels) > 1) {
    for (i in seq_len(length(panels) - 1)) {
      panels[[i]] <- panels[[i]] +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
  }

  stacked <- cowplot::plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "l")

  cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(tr, x = 0, hjust = 0, fontface = "bold", size = 14),
    stacked,
    ncol = 1,
    rel_heights = c(0.06, 0.94)
  )
}

make_group_alltraits_figure <- function(group_name, ct_order_full, agg_choice) {

  df_g <- sig_only %>%
    filter(agg == agg_choice, celltype %in% ct_order_full, trait %in% traits_wanted) %>%
    mutate(celltype = as.character(celltype))

  if (nrow(df_g) == 0) return(NULL)

  # Shared x-limits across all four traits in this group
  lim  <- tight_lim(max(abs(df_g$rho), na.rm = TRUE))
  brks <- pretty(c(-lim, lim), n = 6)

  trait_order <- c("Braak", "CERAD", "delta_mmse", "Sex_numeric")

  panels <- lapply(trait_order, function(tr) {
    make_trait_panel(df_g, tr, agg_choice, ct_order_full, lim = lim, brks = brks)
  })

  # Empty placeholder for traits with no significant results
  panels <- lapply(seq_along(panels), function(i) {
    if (!is.null(panels[[i]])) return(panels[[i]])
    tr <- trait_order[i]
    ggdraw() +
      draw_label(paste0(tr, "\n(no FDR<", alpha_adj, ")"), size = 12, fontface = "bold", x = 0.5, y = 0.6) +
      draw_label("No significant modules", size = 10, x = 0.5, y = 0.45)
  })

  grid_traits <- cowplot::plot_grid(plotlist = panels, ncol = 2, align = "hv", axis = "tblr")

  # Legend
  leg <- ggplot(data.frame(trait = factor(trait_order, levels = trait_order),
                           x = 1:4, y = 1:4),
                aes(x = x, y = y, fill = trait)) +
    geom_col() +
    scale_fill_manual(values = trait_colors, name = "Trait") +
    theme_void(base_size = 11) +
    theme(legend.position = "right")

  legend <- cowplot::get_legend(leg)

  header <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("DONOR-level Spearman correlations (FDR<", alpha_adj, "): ", group_name, "  [", agg_choice, "]"),
      x = 0, hjust = 0, fontface = "bold", size = 16
    ) +
    cowplot::draw_label(
      "Only significant modules shown (celltype-specific). Bars padded to fixed slots for uniform thickness. Dotted line = 0.",
      x = 0, y = 0.55, hjust = 0, size = 11
    )

  body <- cowplot::plot_grid(grid_traits, legend, ncol = 2, rel_widths = c(1, 0.12))

  out <- cowplot::plot_grid(header, body, ncol = 1, rel_heights = c(0.10, 0.90))

  # Global axis labels
  out <- ggdraw(out) +
    draw_label("Spearman \u03c1 across donors", x = 0.5, y = 0.015, hjust = 0.5, size = 12) +
    draw_label("Module (significant only)", x = 0.012, y = 0.5, angle = 90, hjust = 0.5, size = 12)

  out
}

# --- Save group figures ------------------------------------------------------

for (agg_choice in aggs_to_run) {
  for (group_name in names(celltype_groups)) {
    ct_order_full <- celltype_groups[[group_name]]
    fig <- make_group_alltraits_figure(group_name, ct_order_full, agg_choice)
    if (is.null(fig)) next
    fn <- paste0("ALL4TRAITS_", group_name, "_DONORLEVEL_", agg_choice, "_SIGONLY.pdf")
    ggsave(fn, fig, width = 18, height = 12, dpi = 300, bg = "white")
    message("Saved: ", fn)
  }
}

message("DONE. Generated 3 group figures per agg (Neurons/Glia/VascOPC), each containing 4 traits on the same figure with shared x-scale.")
