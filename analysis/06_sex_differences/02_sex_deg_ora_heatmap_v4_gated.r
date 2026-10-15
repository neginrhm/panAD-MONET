###############################################################################
# Sex-DEG x hdWGCNA module ORA heatmaps -- FULL modules, v4 (publication, gated)
#
# v4 differences vs v3:
#   1) STAT ANNOTATIONS GATED. Stars / OR / % labels render only for cells with
#         FDR < 0.05  AND  overlap q >= 4  AND  module overlap >= 2%
#      Failing cells keep their color but show no text. n.t. cells unchanged.
#   2) Combined Supplementary Table (sex_ORA_supplementary_table_v4.csv) with
#      every (cell_type, direction, module) tuple's full statistics, plus a
#      passes_heatmap_gate flag. Per-cell-type CSVs are extended likewise.
#   3) Color scale, palette, and layout match v3.
#
# Output dir: deg_module_enrichment_v3/
#   sex_<ct>_ORA_heatmap_capped_fullmodules_v4.pdf       (per cell type)
#   sex_combined_ORA_heatmap_capped_fullmodules_v4.pdf   (unified panel)
#   sex_<ct>_ORA_results_fullmodules_v4.csv              (per cell type)
#   sex_ORA_supplementary_table_v4.csv
#   sex_ORA_scale_summary_v4.txt
###############################################################################

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(grid)

BASE_DIR <- getwd()  # Set to project root before running
source(file.path(BASE_DIR, "analysis", "utils", "ORA.R"))

MIN_DEG            <- 10
OR_LABEL_THRESHOLD <- 2
GATE_FDR           <- 0.05
GATE_Q             <- 4
GATE_PCT           <- 2
deg_dir  <- file.path(BASE_DIR, "sex_DEG_analysis/pipeline_v2/ad_only_degs")  # written by 01_sex_deg_pipeline.r
out_dir  <- file.path(BASE_DIR, "sex_DEG_analysis/module_enrichment_v4_gated")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ct_root <- file.path(BASE_DIR, "CellTypes")  # full (non-hub-only) hdWGCNA module tables per cell type
celltype_info <- list(
  ast = list(label = "Astrocytes",
             deg_file = file.path(deg_dir, "ast_AD_only_sex_DEGs.csv"),
             mod_file = file.path(ct_root, "Astrocyte/13ad_astrocye_modules_table.csv")),
  exc = list(label = "Excitatory neurons",
             deg_file = file.path(deg_dir, "exc_AD_only_sex_DEGs.csv"),
             mod_file = file.path(ct_root, "Excitatory/13ad_excitatory_modules_table.csv")),
  inh = list(label = "Inhibitory neurons",
             deg_file = file.path(deg_dir, "inh_AD_only_sex_DEGs.csv"),
             mod_file = file.path(ct_root, "Inhibitory/13ad_INH_modules_table.csv")),
  mic = list(label = "Microglia",
             deg_file = file.path(deg_dir, "mic_AD_only_sex_DEGs.csv"),
             mod_file = file.path(ct_root, "Microglia/hdWGCNA/13ad_microglia_modules_table.csv")),
  odc = list(label = "Oligodendrocytes",
             deg_file = file.path(deg_dir, "odc_AD_only_sex_DEGs.csv"),
             mod_file = file.path(ct_root, "ODC/13ad_odc_modules_table.csv"))
)
celltype_order <- names(celltype_info)

module_rename <- list(
  ast = c(brown = "ASC-M1", turquoise = "ASC-M2", pink = "ASC-M3",
          blue = "ASC-M4", yellow = "ASC-M5", magenta = "ASC-M6",
          purple = "ASC-M7", black = "ASC-M8", greenyellow = "ASC-M9",
          green = "ASC-M10", red = "ASC-M11"),
  mic = c(brown = "MIC-M1", turquoise = "MIC-M2", blue = "MIC-M4",
          yellow = "MIC-M5", green = "MIC-M10"),
  odc = c(brown = "ODC-M1", turquoise = "ODC-M2", blue = "ODC-M4"),
  exc = c(green = "EXC-M1", turquoise = "EXC-M2", yellow = "EXC-M3",
          blue = "EXC-M4", red = "EXC-M5", brown = "EXC-M6"),
  inh = c(brown = "INH-M1", turquoise = "INH-M2", pink = "INH-M3",
          blue = "INH-M4", yellow = "INH-M5", magenta = "INH-M6",
          purple = "INH-M7", black = "INH-M8", greenyellow = "INH-M9",
          green = "INH-M10", red = "INH-M11", tan = "INH-M12")
)
relabel_modules <- function(ct_key, colors) {
  mapping <- module_rename[[ct_key]]
  if (is.null(mapping)) return(as.character(colors))
  new_names <- unname(mapping[as.character(colors)])
  miss <- is.na(new_names); new_names[miss] <- as.character(colors[miss])
  new_names
}

get_gene_col <- function(df) {
  if ("gene" %in% colnames(df))      df$gene
  else if ("Gene" %in% colnames(df)) df$Gene
  else if ("X" %in% colnames(df))    df$X
  else                               rownames(df)
}

# ---------------------------------------------------------------------------
# Pass 1: compute ORAs + apply MIN_DEG filter
# ---------------------------------------------------------------------------
results_per_ct <- list()
all_finite_ORs <- c()
drop_log <- character()
supp_rows <- list()

for (ct in celltype_order) {
  info <- celltype_info[[ct]]
  if (!file.exists(info$deg_file) || !file.exists(info$mod_file)) {
    cat(sprintf("[skip] %s: missing input\n", ct)); next
  }
  module_data <- read.csv(info$mod_file, stringsAsFactors = FALSE)
  background   <- as.character(module_data$gene_name)
  datKME       <- subset(module_data, module != "grey",
                         select = c("gene_name", "module"))
  module_lists <- split(as.character(datKME$gene_name), datKME$module)
  module_sizes <- sapply(module_lists, length)

  degs <- read.csv(info$deg_file, stringsAsFactors = FALSE)
  degs$gene <- get_gene_col(degs)
  male_genes   <- as.character(degs$gene[degs$avg_log2FC > 0])
  female_genes <- as.character(degs$gene[degs$avg_log2FC < 0])

  male_in_bg   <- length(intersect(male_genes,   background))
  female_in_bg <- length(intersect(female_genes, background))

  keep_dirs <- c()
  if (male_in_bg   >= MIN_DEG) keep_dirs <- c(keep_dirs, "Up_in_Male")
  else if (male_in_bg > 0)
    drop_log <- c(drop_log, sprintf("%s Up_in_Male dropped (k=%d < %d)",
                                    ct, male_in_bg, MIN_DEG))
  if (female_in_bg >= MIN_DEG) keep_dirs <- c(keep_dirs, "Up_in_Female")
  else if (female_in_bg > 0)
    drop_log <- c(drop_log, sprintf("%s Up_in_Female dropped (k=%d < %d)",
                                    ct, female_in_bg, MIN_DEG))
  if (length(keep_dirs) == 0) {
    cat(sprintf("  [skip] %s: no directions survive MIN_DEG\n", ct)); next
  }

  deg_lists <- list()
  if ("Up_in_Male"   %in% keep_dirs) deg_lists$Up_in_Male   <- male_genes
  if ("Up_in_Female" %in% keep_dirs) deg_lists$Up_in_Female <- female_genes
  n_in_bg_vec <- c(Up_in_Male = male_in_bg, Up_in_Female = female_in_bg)

  ORA.OR <- matrix(NA_real_, length(deg_lists), length(module_lists),
                   dimnames = list(names(deg_lists), names(module_lists)))
  ORA.P <- ORA.OR; pctMat <- ORA.OR; qMat <- ORA.OR
  for (i in seq_along(module_lists)) {
    for (j in seq_along(deg_lists)) {
      r <- ORA(deg_lists[[j]], module_lists[[i]], background, background)
      ORA.OR[j, i] <- as.numeric(r["OR"])
      ORA.P [j, i] <- as.numeric(r["Fisher p"])
      q <- as.numeric(r["Overlap"])
      ms <- length(module_lists[[i]])
      qMat[j, i]   <- q
      pctMat[j, i] <- if (!is.na(q) && ms > 0) round((q/ms)*100, 1) else NA
    }
  }
  FDRmat <- matrix(p.adjust(ORA.P, method = "fdr"),
                   nrow = nrow(ORA.P), dimnames = dimnames(ORA.P))
  all_finite_ORs <- c(all_finite_ORs, ORA.OR[is.finite(ORA.OR)])

  module_labels <- relabel_modules(ct, colnames(ORA.OR))
  long_df <- data.frame(
    cell_type            = info$label,
    direction            = rep(rownames(ORA.OR), times = ncol(ORA.OR)),
    n_DEGs_in_background = rep(n_in_bg_vec[rownames(ORA.OR)], times = ncol(ORA.OR)),
    module_color         = rep(colnames(ORA.OR), each = nrow(ORA.OR)),
    module_label         = rep(module_labels, each = nrow(ORA.OR)),
    module_size          = rep(module_sizes[colnames(ORA.OR)], each = nrow(ORA.OR)),
    overlap_q            = as.vector(qMat),
    odds_ratio           = as.vector(ORA.OR),
    pct_module_overlap   = as.vector(pctMat),
    fisher_p             = as.vector(ORA.P),
    fdr                  = as.vector(FDRmat),
    stringsAsFactors     = FALSE
  )
  long_df$passes_heatmap_gate <- with(long_df,
    !is.na(fdr) & fdr < GATE_FDR &
    !is.na(overlap_q) & overlap_q >= GATE_Q &
    !is.na(pct_module_overlap) & pct_module_overlap >= GATE_PCT
  )

  write.csv(long_df, file.path(out_dir,
            sprintf("sex_%s_ORA_results_fullmodules_v4.csv", ct)),
            row.names = FALSE)
  supp_rows[[ct]] <- long_df

  rownames(ORA.OR) <- paste0(ct, "_", rownames(ORA.OR))
  rownames(FDRmat) <- rownames(pctMat) <- rownames(qMat) <- rownames(ORA.OR)
  results_per_ct[[ct]] <- list(or_raw = ORA.OR, fdr = FDRmat,
                               pct = pctMat, q = qMat)
}

if (length(all_finite_ORs) == 0) stop("No finite ORs to calibrate scale.")

# ---------------------------------------------------------------------------
# Data-driven color scale (v3: p90 top / p95 cap, median mid)
# ---------------------------------------------------------------------------
p50 <- quantile(all_finite_ORs, 0.50, na.rm = TRUE, names = FALSE)
p90 <- quantile(all_finite_ORs, 0.90, na.rm = TRUE, names = FALSE)
p95 <- quantile(all_finite_ORs, 0.95, na.rm = TRUE, names = FALSE)
or_max <- max(all_finite_ORs, na.rm = TRUE)

OR_TOP <- max(2,         ceiling(p90 * 2) / 2)
OR_CAP <- max(OR_TOP + 0.5, ceiling(p95 * 2) / 2)
OR_MID <- min(OR_TOP - 0.5, max(1.2, round(p50, 1)))

col_fun <- colorRamp2(
  c(0,         1,         OR_MID,    OR_TOP,    OR_CAP),
  c("#2E0854", "#6E1B8C", "#E8419A", "#FFD93D", "#B9F8D3")
)

summary_lines <- c(
  "Data-driven OR scale for sex-DEG heatmaps (v4, full modules, gated)",
  sprintf("Finite OR cells: n = %d across retained rows", length(all_finite_ORs)),
  sprintf("  median (p50) = %.3f, p90 = %.3f, p95 = %.3f, max = %.3f",
          p50, p90, p95, or_max),
  sprintf("Anchors: 0, 1, %.1f (mid~p50), %.1f (top=p90), %.1f (cap=p95/Inf)",
          OR_MID, OR_TOP, OR_CAP),
  sprintf("MIN_DEG filter (row admission): %d DEGs in background", MIN_DEG),
  sprintf("Heatmap text gate: FDR < %.2f AND q >= %d AND module overlap >= %.0f%%",
          GATE_FDR, GATE_Q, GATE_PCT),
  "Cells failing the gate retain their OR-derived color but display no text.",
  "Palette: #2E0854 -> #6E1B8C -> #E8419A -> #FFD93D -> #B9F8D3",
  "",
  "Rows dropped for low power:",
  if (length(drop_log)) paste("  ", drop_log, collapse = "\n") else "  (none)"
)
writeLines(summary_lines, file.path(out_dir, "sex_ORA_scale_summary_v4.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")

legend_param <- list(
  title  = "Odds Ratio",
  at     = c(0, 1, OR_MID, OR_TOP, OR_CAP),
  labels = c("0", "1",
             sprintf("%.1f", OR_MID),
             sprintf("%.1f (p90)", OR_TOP),
             sprintf("max\n(Inf, >= %.1f)", OR_CAP))
)

# ---------------------------------------------------------------------------
# Build text annotations + cell drawing (gated)
# ---------------------------------------------------------------------------
build_txt <- function(fdr_mat, pct_mat, or_raw, disp, q_mat) {
  txt <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat),
                dimnames = dimnames(fdr_mat))
  gate <- !is.na(fdr_mat) & fdr_mat < GATE_FDR &
          !is.na(q_mat)   & q_mat   >= GATE_Q  &
          !is.na(pct_mat) & pct_mat >= GATE_PCT
  for (i in seq_len(nrow(fdr_mat))) for (j in seq_len(ncol(fdr_mat))) {
    if (!isTRUE(gate[i, j])) next
    f <- fdr_mat[i, j]
    stars <- if (is.na(f))       ""
             else if (f < 0.005) "***"
             else if (f < 0.01)  "**"
             else if (f < 0.05)  "*"
             else                ""
    or_val <- disp[i, j]; raw <- or_raw[i, j]
    num_str <- if (is.na(or_val)) {
      if (!is.na(pct_mat[i, j])) paste0(pct_mat[i, j], "%") else ""
    } else if (or_val >= OR_LABEL_THRESHOLD) {
      if (!is.finite(raw)) paste0("OR>=", OR_CAP)
      else paste0("OR=", formatC(raw, format = "f", digits = 1))
    } else {
      if (!is.na(pct_mat[i, j])) paste0(pct_mat[i, j], "%") else ""
    }
    txt[i, j] <- if (nzchar(num_str)) paste0(stars, "\n", num_str) else stars
  }
  txt
}

make_cell_fun <- function(txtMat, valMat) {
  function(j, i, x, y, width, height, fill) {
    v <- valMat[i, j]
    if (is.na(v)) {
      grid.text("n.t.", x, y,
                gp = gpar(fontsize = 6, fontface = "italic", col = "grey60"))
      return(invisible())
    }
    label <- txtMat[i, j]
    if (!nzchar(label)) return(invisible())  # gate failed -- color only
    parts <- unlist(strsplit(label, "\n"))
    stars <- if (length(parts) >= 1) parts[1] else ""
    num   <- if (length(parts) >= 2) parts[2] else ""
    text_col <- if (v >= OR_MID) "black" else "white"
    if (nzchar(stars))
      grid.text(stars, x, y + unit(1.5, "mm"),
                gp = gpar(fontsize = 20, fontface = "bold", col = text_col))
    if (nzchar(num))
      grid.text(num, x, y - unit(1.5, "mm"),
                gp = gpar(fontsize = 6, col = text_col))
  }
}

direction_colors <- c(Up_in_Male = "steelblue", Up_in_Female = "indianred")

gate_footer <- sprintf(
  "Annotated: FDR<%.2f & q>=%d & overlap>=%.0f%%",
  GATE_FDR, GATE_Q, GATE_PCT
)

# ---------------------------------------------------------------------------
# Per-cell-type heatmaps
# ---------------------------------------------------------------------------
all_results <- list()
for (ct in names(results_per_ct)) {
  R <- results_per_ct[[ct]]
  disp <- R$or_raw
  disp[is.infinite(disp)] <- OR_CAP
  disp[is.nan(disp)]      <- NA
  txt <- build_txt(R$fdr, R$pct, R$or_raw, disp, R$q)
  all_results[[ct]] <- list(or = disp, fdr = R$fdr, pct = R$pct,
                            or_raw = R$or_raw, txt = txt, q = R$q)

  meta_ct <- data.frame(
    row = rownames(disp),
    direction = ifelse(grepl("Up_in_Male", rownames(disp)),
                       "Up_in_Male", "Up_in_Female"),
    stringsAsFactors = FALSE
  )
  row_anno <- rowAnnotation(
    direction = meta_ct$direction,
    col = list(direction = direction_colors),
    show_annotation_name = TRUE
  )

  pdf_file <- file.path(out_dir,
                        sprintf("sex_%s_ORA_heatmap_capped_fullmodules_v4.pdf", ct))
  pdf(pdf_file, width = 10, height = 4)
  draw(Heatmap(
    matrix = disp, name = "Odds Ratio", col = col_fun, na_col = "grey92",
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = make_cell_fun(txt, disp),
    row_names_side = "left", column_names_rot = 45,
    column_labels = relabel_modules(ct, colnames(disp)),
    row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10),
    left_annotation = row_anno,
    column_title = celltype_info[[ct]]$label,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = legend_param
  ), heatmap_legend_side = "right")
  grid.text(gate_footer,
            x = unit(1, "npc") - unit(2, "mm"), y = unit(2, "mm"),
            just = c("right", "bottom"),
            gp = gpar(fontsize = 7, col = "grey30", fontface = "italic"))
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(pdf_file)))
}

# ---------------------------------------------------------------------------
# Combined heatmap across all cell types
# ---------------------------------------------------------------------------
if (length(all_results) > 0) {
  all_colors <- unique(unlist(lapply(all_results, function(x) colnames(x$or))))
  dispMat_all <- NULL; txtMat_all <- NULL
  for (ct in names(all_results)) {
    dm <- all_results[[ct]]$or; tm <- all_results[[ct]]$txt
    dm_full <- matrix(NA_real_, nrow = nrow(dm), ncol = length(all_colors),
                      dimnames = list(rownames(dm), all_colors))
    tm_full <- matrix("", nrow = nrow(tm), ncol = length(all_colors),
                      dimnames = list(rownames(tm), all_colors))
    shared <- intersect(colnames(dm), all_colors)
    dm_full[, shared] <- dm[, shared]
    tm_full[, shared] <- tm[, shared]
    dispMat_all <- rbind(dispMat_all, dm_full)
    txtMat_all  <- rbind(txtMat_all,  tm_full)
  }

  metadata_ora <- data.frame(
    row = rownames(dispMat_all),
    celltype = gsub("_Up_in_Male$|_Up_in_Female$", "", rownames(dispMat_all)),
    direction = ifelse(grepl("Up_in_Male", rownames(dispMat_all)),
                       "Up_in_Male", "Up_in_Female"),
    stringsAsFactors = FALSE
  )
  ct_levels <- unique(metadata_ora$celltype)
  n_ct <- max(3, length(ct_levels))
  ct_colors <- setNames(brewer.pal(n_ct, "Set3")[seq_along(ct_levels)], ct_levels)

  row_anno_all <- rowAnnotation(
    celltype = metadata_ora$celltype,
    direction = metadata_ora$direction,
    col = list(celltype = ct_colors, direction = direction_colors),
    show_annotation_name = TRUE
  )

  cts_present <- names(all_results)
  compact_col_labels <- vapply(all_colors, function(color) {
    parts <- character(0)
    for (ct in cts_present) {
      mapping <- module_rename[[ct]]
      if (is.null(mapping)) next
      nm <- mapping[color]
      if (!is.na(nm)) parts <- c(parts, unname(nm))
    }
    if (length(parts) == 0) return(color)
    paste(parts, collapse = "\n")
  }, character(1))

  pdf_file <- file.path(out_dir, "sex_combined_ORA_heatmap_capped_fullmodules_v4.pdf")
  pdf(pdf_file, width = 11, height = 7)
  draw(Heatmap(
    matrix = dispMat_all, name = "Odds Ratio",
    col = col_fun, na_col = "grey92",
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = make_cell_fun(txtMat_all, dispMat_all),
    row_names_side = "left", column_names_side = "bottom",
    column_names_rot = 0, column_names_centered = TRUE,
    column_labels = compact_col_labels,
    row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 7),
    left_annotation = row_anno_all,
    column_title = "AD-only sex DEGs x hdWGCNA modules",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = legend_param
  ), heatmap_legend_side = "right")
  grid.text(gate_footer,
            x = unit(1, "npc") - unit(2, "mm"), y = unit(2, "mm"),
            just = c("right", "bottom"),
            gp = gpar(fontsize = 7, col = "grey30", fontface = "italic"))
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(pdf_file)))
}

# ---------------------------------------------------------------------------
# Combined Supplementary Table
# ---------------------------------------------------------------------------
if (length(supp_rows) > 0) {
  supp <- do.call(rbind, supp_rows)
  rownames(supp) <- NULL
  supp <- supp[order(supp$cell_type, supp$direction, supp$module_label), ]
  write.csv(supp,
            file.path(out_dir, "sex_ORA_supplementary_table_v4.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved supplementary table: %d rows\n", nrow(supp)))
}

cat(sprintf("\nDone. Output: %s\n", out_dir))
