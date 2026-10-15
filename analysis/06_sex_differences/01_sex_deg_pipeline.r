###############################################################################
# Sex-Differential Analysis Pipeline v2: AD vs Control
#
# Steps:
#   1. Remove sex-chromosome (X/Y) genes from all DEG CSV files & save cleaned
#   2. Filter AD DEGs p_adj < 0.05, Venn diagram vs Control p_adj < 0.05,
#      subset to AD-only DEGs
#   3. ME comparison by sex in Dementia samples (Wilcoxon + boxplots)
#   4. ORA enrichment: AD-only sex DEGs against hdWGCNA modules (heatmaps)
#   5. Compact bar chart of ME sex differences across all cell types
###############################################################################

library(dplyr)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(12345)

# ==============================================================================
# Output directories
# ==============================================================================
BASE_DIR <- getwd()  # Set to project root before running
source(file.path(BASE_DIR, "analysis", "utils", "ORA.R"))
out_base   <- file.path(BASE_DIR, "sex_DEG_analysis/pipeline_v2")
clean_dir  <- file.path(out_base, "cleaned_DEGs")
venn_dir   <- file.path(out_base, "venn_diagrams")
deg_dir    <- file.path(out_base, "ad_only_degs")
me_dir     <- file.path(out_base, "ME_by_sex")
ora_dir    <- file.path(out_base, "module_enrichment")
plot_dir   <- file.path(out_base, "plots")

for (d in c(out_base, clean_dir, venn_dir, deg_dir, me_dir, ora_dir, plot_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ==============================================================================
# STEP 0: Build X/Y chromosome gene list
# ==============================================================================
cat("Building X/Y gene filter...\n")
xy_annot <- AnnotationDbi::select(org.Hs.eg.db, keys = c("X", "Y"),
                                   keytype = "CHR", columns = c("SYMBOL", "CHR"))
xy_genes <- unique(xy_annot$SYMBOL)
cat(sprintf("  %d X/Y genes will be excluded\n", length(xy_genes)))
write.csv(data.frame(gene = xy_genes, stringsAsFactors = FALSE),
          file.path(out_base, "xy_genes_excluded.csv"), row.names = FALSE)

# ==============================================================================
# Cell type configuration
# ==============================================================================
mast_dir <- file.path(BASE_DIR, "sex_DEG_analysis/MAST_results")

celltype_info <- list(
  ast = list(
    label     = "Astrocyte",
    dem_csv   = file.path(mast_dir, "Dementia/Dementia_ast_male_vs_female_normalized.csv"),
    ctrl_csv  = file.path(mast_dir, "Control/control_ast_male_vs_female_normalized.csv"),
    me_file   = file.path(BASE_DIR, "CellTypes/Astrocyte/ast_MEs_by_sample_manual.csv"),
    hub_file  = file.path(BASE_DIR, "CellTypes/Astrocyte/13ad_astrocye_modules_hub_df_table.csv"),
    mod_col   = "module"
  ),
  exc = list(
    label     = "Excitatory",
    dem_csv   = file.path(mast_dir, "Dementia/Dementia_exc_male_vs_female_normalized.csv"),
    ctrl_csv  = file.path(mast_dir, "Control/control_exc_male_vs_female_normalized.csv"),
    me_file   = file.path(BASE_DIR, "CellTypes/Excitatory/exc_MEs_by_sample.csv"),
    hub_file  = file.path(BASE_DIR, "CellTypes/Excitatory/13ad_excitatory_modules_hub_df_table.csv"),
    mod_col   = "module"
  ),
  inh = list(
    label     = "Inhibitory",
    dem_csv   = file.path(mast_dir, "Dementia/Dementia_inh_male_vs_female_normalized.csv"),
    ctrl_csv  = file.path(mast_dir, "Control/control_inh_male_vs_female_normalized.csv"),
    me_file   = file.path(BASE_DIR, "CellTypes/Inhibitory/inh_MEs_by_sample_manual.csv"),
    hub_file  = file.path(BASE_DIR, "CellTypes/Inhibitory/13ad_INH_modules_hub_df_table.csv"),
    mod_col   = "module"
  ),
  mic = list(
    label     = "Microglia",
    dem_csv   = file.path(mast_dir, "Dementia/Dementia_mic_male_vs_female_normalized.csv"),
    ctrl_csv  = file.path(mast_dir, "Control/control_mic_male_vs_female_normalized.csv"),
    me_file   = file.path(BASE_DIR, "CellTypes/Microglia/hdWGCNA/mic_MEs_by_sample_manual.csv"),
    hub_file  = file.path(BASE_DIR, "CellTypes/Microglia/hdWGCNA/13ad_microglia_modules_hub_df_table.csv"),
    mod_col   = "module"
  ),
  odc = list(
    label     = "Oligodendrocyte",
    dem_csv   = file.path(mast_dir, "Dementia/Dementia_odc_male_vs_female_normalized.csv"),
    ctrl_csv  = file.path(mast_dir, "Control/control_odc_male_vs_female_normalized.csv"),
    me_file   = file.path(BASE_DIR, "CellTypes/ODC/odc_MEs_by_sample_manual.csv"),
    hub_file  = file.path(BASE_DIR, "CellTypes/ODC/13ad_odc_modules_hub_df_table.csv"),
    mod_col   = "module"
  )
)

# ==============================================================================
# Load metadata: one row per sample with Sex and dx
# ==============================================================================
cat("Loading metadata...\n")
meta_full <- read.csv(file.path(BASE_DIR, "ipynb/not_ipynb/13ad_metadata.csv"),
                      stringsAsFactors = FALSE)
sample_meta <- meta_full %>%
  dplyr::select(Sample, Sex, dx) %>%
  distinct(Sample, .keep_all = TRUE)
colnames(sample_meta)[3] <- "Diagnosis"
cat(sprintf("  %d unique samples\n", nrow(sample_meta)))

# Load module GO annotations for bar chart labels
go_annot <- tryCatch({
  read.csv(file.path(BASE_DIR,
    "sex_DEG_analysis/first_try/results/module_enrichment/module_GO_annotations_top1.csv"),
    stringsAsFactors = FALSE)
}, error = function(e) NULL)

# ==============================================================================
# Helper: extract gene column from DEG data frame (handles X vs Gene column)
# ==============================================================================
get_gene_col <- function(df) {
  if ("Gene" %in% colnames(df)) {
    return(df$Gene)
  } else if ("X" %in% colnames(df)) {
    return(df$X)
  } else {
    # Fall back to row names or first column
    return(rownames(df))
  }
}

# ==============================================================================
# Venn diagram helper
# ==============================================================================
draw_venn2 <- function(set_a, set_b, label_a, label_b, title, file_path) {
  n_a    <- length(set_a)
  n_b    <- length(set_b)
  n_both <- length(intersect(set_a, set_b))
  n_a_only <- n_a - n_both
  n_b_only <- n_b - n_both

  pdf(file_path, width = 7, height = 6)
  grid.newpage()
  grid.circle(x = 0.38, y = 0.5, r = 0.22,
              gp = gpar(fill = adjustcolor("#E74C3C", 0.4), col = "#E74C3C", lwd = 2))
  grid.circle(x = 0.62, y = 0.5, r = 0.22,
              gp = gpar(fill = adjustcolor("#3498DB", 0.4), col = "#3498DB", lwd = 2))
  grid.text(as.character(n_a_only), x = 0.28, y = 0.5,
            gp = gpar(fontsize = 18, fontface = "bold"))
  grid.text(as.character(n_both), x = 0.50, y = 0.5,
            gp = gpar(fontsize = 18, fontface = "bold"))
  grid.text(as.character(n_b_only), x = 0.72, y = 0.5,
            gp = gpar(fontsize = 18, fontface = "bold"))
  grid.text(paste0(label_a, "\n(n=", n_a, ")"), x = 0.25, y = 0.82,
            gp = gpar(fontsize = 11, col = "#E74C3C"))
  grid.text(paste0(label_b, "\n(n=", n_b, ")"), x = 0.75, y = 0.82,
            gp = gpar(fontsize = 11, col = "#3498DB"))
  grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
}

###############################################################################
# STEP 1: Remove sex-chromosome genes from all DEG files & save cleaned CSVs
###############################################################################
cat("\n=== STEP 1: Remove X/Y genes from all DEG CSV files ===\n")

for (ct in names(celltype_info)) {
  info <- celltype_info[[ct]]
  for (cond in c("dem", "ctrl")) {
    csv_path <- if (cond == "dem") info$dem_csv else info$ctrl_csv
    cond_label <- if (cond == "dem") "Dementia" else "Control"

    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    df$gene <- get_gene_col(df)
    n_before <- nrow(df)
    df_clean <- df %>% filter(!gene %in% xy_genes)
    n_removed <- n_before - nrow(df_clean)

    out_name <- paste0(ct, "_", cond_label, "_male_vs_female_no_sex_genes.csv")
    write.csv(df_clean, file.path(clean_dir, out_name), row.names = FALSE)
    cat(sprintf("  %s %s: %d -> %d genes (%d X/Y removed)\n",
                info$label, cond_label, n_before, nrow(df_clean), n_removed))
  }
}

###############################################################################
# STEP 2: Filter + Venn + AD-only DEGs
###############################################################################
cat("\n=== STEP 2: p_adj<0.05 filtering, Venn diagrams, AD-only DEGs ===\n")

ad_only_degs   <- list()
venn_summary   <- data.frame()

for (ct in names(celltype_info)) {
  info <- celltype_info[[ct]]
  cat(sprintf("\n  %s:\n", info$label))

  # Load cleaned DEGs
  dem_clean  <- read.csv(file.path(clean_dir,
                  paste0(ct, "_Dementia_male_vs_female_no_sex_genes.csv")),
                  stringsAsFactors = FALSE)
  ctrl_clean <- read.csv(file.path(clean_dir,
                  paste0(ct, "_Control_male_vs_female_no_sex_genes.csv")),
                  stringsAsFactors = FALSE)

  # Filter p_adj < 0.05
  dem_sig  <- dem_clean  %>% filter(p_val_adj < 0.05)
  ctrl_sig <- ctrl_clean %>% filter(p_val_adj < 0.05)

  dem_genes  <- dem_sig$gene
  ctrl_genes <- ctrl_sig$gene

  # Venn categories
  ad_only   <- setdiff(dem_genes, ctrl_genes)
  shared    <- intersect(dem_genes, ctrl_genes)
  ctrl_only <- setdiff(ctrl_genes, dem_genes)

  cat(sprintf("    AD sig: %d | Ctrl sig: %d\n", length(dem_genes), length(ctrl_genes)))
  cat(sprintf("    AD-only: %d | Shared: %d | Ctrl-only: %d\n",
              length(ad_only), length(shared), length(ctrl_only)))

  # Store AD-only DEGs with full info
  ad_only_df <- dem_sig %>% filter(gene %in% ad_only)
  ad_only_degs[[ct]] <- ad_only_df

  write.csv(ad_only_df, file.path(deg_dir, paste0(ct, "_AD_only_sex_DEGs.csv")),
            row.names = FALSE)
  write.csv(dem_sig %>% filter(gene %in% shared),
            file.path(deg_dir, paste0(ct, "_shared_sex_DEGs.csv")),
            row.names = FALSE)
  write.csv(ctrl_sig %>% filter(gene %in% ctrl_only),
            file.path(deg_dir, paste0(ct, "_control_only_sex_DEGs.csv")),
            row.names = FALSE)

  venn_summary <- rbind(venn_summary, data.frame(
    celltype = info$label, AD_sig = length(dem_genes), Ctrl_sig = length(ctrl_genes),
    AD_only = length(ad_only), Shared = length(shared), Control_only = length(ctrl_only)
  ))

  # Venn diagram
  draw_venn2(dem_genes, ctrl_genes,
             "M vs F\n(AD, p<0.05)", "M vs F\n(Control, p<0.05)",
             paste0(info$label, ": Sex-Differential Genes (X/Y removed)"),
             file.path(venn_dir, paste0(ct, "_sex_DEG_venn.pdf")))
}

write.csv(venn_summary, file.path(out_base, "venn_summary.csv"), row.names = FALSE)
cat("\nVenn summary:\n")
print(venn_summary)

###############################################################################
# STEP 3: ME comparison by sex (Dementia samples)
###############################################################################
cat("\n=== STEP 3: Module Eigengene comparison by Sex (Dementia) ===\n")

me_test_all <- data.frame()

for (ct in names(celltype_info)) {
  info <- celltype_info[[ct]]
  cat(sprintf("\n  %s MEs...\n", info$label))

  if (!file.exists(info$me_file)) {
    cat("    ME file not found, skipping.\n")
    next
  }

  me_df <- read.csv(info$me_file, stringsAsFactors = FALSE)

  # Merge with metadata
  me_df$Sample <- as.character(me_df$Sample)
  sample_meta$Sample <- as.character(sample_meta$Sample)
  me_merged <- merge(me_df, sample_meta, by = "Sample", all.x = TRUE)

  # Filter to Dementia
  me_dem <- me_merged %>% filter(Diagnosis == "Dementia")
  n_m <- sum(me_dem$Sex == "Male", na.rm = TRUE)
  n_f <- sum(me_dem$Sex == "Female", na.rm = TRUE)
  cat(sprintf("    %d Dementia samples (M: %d, F: %d)\n", nrow(me_dem), n_m, n_f))

  if (n_m < 3 || n_f < 3) {
    cat("    Too few samples, skipping.\n")
    next
  }

  # Module columns (exclude metadata cols and grey)
  module_cols <- setdiff(colnames(me_dem),
                         c("X", "Sample", "Sex", "Diagnosis", "Age", "grey"))

  for (mod in module_cols) {
    male_vals   <- me_dem %>% filter(Sex == "Male")   %>% pull(!!sym(mod))
    female_vals <- me_dem %>% filter(Sex == "Female") %>% pull(!!sym(mod))

    male_vals   <- male_vals[!is.na(male_vals)]
    female_vals <- female_vals[!is.na(female_vals)]

    if (length(male_vals) < 3 || length(female_vals) < 3) next

    wt <- wilcox.test(male_vals, female_vals)
    diff_median <- median(male_vals) - median(female_vals)

    me_test_all <- rbind(me_test_all, data.frame(
      celltype     = info$label,
      ct_key       = ct,
      module       = mod,
      median_male  = median(male_vals),
      median_female = median(female_vals),
      diff_median  = diff_median,
      p_value      = wt$p.value,
      stringsAsFactors = FALSE
    ))
  }

  # Boxplot per cell type
  me_long <- me_dem %>%
    dplyr::select(Sample, Sex, all_of(module_cols)) %>%
    pivot_longer(cols = all_of(module_cols), names_to = "Module", values_to = "ME")

  p <- ggplot(me_long, aes(x = Module, y = ME, fill = Sex)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("Male" = "steelblue", "Female" = "indianred")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = paste0(info$label, ": Module Eigengenes by Sex (Dementia)"),
         x = "Module", y = "Module Eigengene") +
    stat_compare_means(aes(group = Sex), method = "wilcox.test",
                       label = "p.signif", size = 3, hide.ns = TRUE)

  ggsave(file.path(me_dir, paste0(ct, "_ME_by_sex_boxplot.pdf")),
         p, width = 12, height = 6)
  cat("    Boxplot saved.\n")
}

# FDR correction
if (nrow(me_test_all) > 0) {
  me_test_all$fdr <- p.adjust(me_test_all$p_value, method = "fdr")
  me_test_all <- me_test_all %>% arrange(p_value)
  write.csv(me_test_all, file.path(me_dir, "ME_sex_wilcoxon_results.csv"),
            row.names = FALSE)

  cat("\nSignificant ME differences (FDR < 0.05):\n")
  sig_me <- me_test_all %>% filter(fdr < 0.05)
  if (nrow(sig_me) > 0) print(sig_me) else cat("  None at FDR < 0.05\n")

  cat("\nNominally significant (p < 0.05):\n")
  nom_me <- me_test_all %>% filter(p_value < 0.05)
  if (nrow(nom_me) > 0) print(nom_me) else cat("  None at p < 0.05\n")
}

###############################################################################
# STEP 4: ORA Enrichment of AD-only Sex DEGs in hdWGCNA Modules
###############################################################################
cat("\n=== STEP 4: ORA Enrichment (AD-only DEGs vs hdWGCNA Modules) ===\n")

all_ora_results <- list()

for (ct in names(celltype_info)) {
  info <- celltype_info[[ct]]
  cat(sprintf("\n  %s...\n", info$label))

  if (!file.exists(info$hub_file)) {
    cat("    Hub file not found, skipping.\n")
    next
  }

  # Module data
  module_data <- read.csv(info$hub_file, stringsAsFactors = FALSE)
  datKME <- subset(module_data, module_data[[info$mod_col]] != "grey",
                   select = c("gene_name", info$mod_col))
  colnames(datKME) <- c("gene_name", "module")
  background    <- as.character(module_data$gene_name)
  module_lists  <- split(as.character(datKME$gene_name), datKME$module)

  # AD-only DEGs
  degs <- ad_only_degs[[ct]]
  if (is.null(degs) || nrow(degs) == 0) {
    cat("    No AD-only DEGs, skipping.\n")
    next
  }

  deg_lists <- list(
    Up_in_Male   = as.character(degs$gene[degs$avg_log2FC > 0]),
    Up_in_Female = as.character(degs$gene[degs$avg_log2FC < 0])
  )

  # ORA matrices
  ORA.OR     <- matrix(NA, nrow = length(deg_lists), ncol = length(module_lists))
  ORA.P      <- matrix(NA, nrow = length(deg_lists), ncol = length(module_lists))
  pctOverlap <- matrix(NA, nrow = length(deg_lists), ncol = length(module_lists))
  colnames(ORA.OR) <- colnames(ORA.P) <- colnames(pctOverlap) <- names(module_lists)
  rownames(ORA.OR) <- rownames(ORA.P) <- rownames(pctOverlap) <- names(deg_lists)

  for (i in seq_along(module_lists)) {
    for (j in seq_along(deg_lists)) {
      result <- ORA(deg_lists[[j]], module_lists[[i]], background, background)
      ORA.OR[j, i]     <- as.numeric(result["OR"])
      ORA.P[j, i]      <- as.numeric(result["Fisher p"])
      q <- as.numeric(result["Overlap"])
      k <- as.numeric(result["Input List"])
      pctOverlap[j, i] <- ifelse(!is.na(q) && k > 0, round((q / k) * 100, 1), NA)
    }
  }

  # FDR
  FDRmat <- matrix(p.adjust(ORA.P, method = "fdr"), nrow = nrow(ORA.P))
  rownames(FDRmat) <- rownames(ORA.P)
  colnames(FDRmat) <- colnames(ORA.P)

  # Text matrix
  txtMat <- matrix("", nrow = nrow(FDRmat), ncol = ncol(FDRmat))
  rownames(txtMat) <- rownames(FDRmat)
  colnames(txtMat) <- colnames(FDRmat)
  for (i in 1:nrow(FDRmat)) {
    for (j in 1:ncol(FDRmat)) {
      stars <- ifelse(FDRmat[i, j] < 0.005, "***",
               ifelse(FDRmat[i, j] < 0.01, "**",
               ifelse(FDRmat[i, j] < 0.05, "*", "")))
      txtMat[i, j] <- ifelse(!is.na(pctOverlap[i, j]),
                             paste0(stars, "\n", pctOverlap[i, j], "%"),
                             stars)
    }
  }

  dispMat <- ORA.OR
  dispMat[is.infinite(dispMat) | is.nan(dispMat)] <- NA
  rownames(dispMat) <- paste0(ct, "_", rownames(dispMat))
  rownames(txtMat)  <- paste0(ct, "_", rownames(txtMat))

  # Save CSV
  ora_csv <- data.frame(
    celltype    = info$label,
    direction   = rep(names(deg_lists), each = ncol(ORA.OR)),
    module      = rep(colnames(ORA.OR), times = nrow(ORA.OR)),
    odds_ratio  = as.vector(t(ORA.OR)),
    p_value     = as.vector(t(ORA.P)),
    fdr         = as.vector(t(FDRmat)),
    pct_overlap = as.vector(t(pctOverlap))
  )
  write.csv(ora_csv, file.path(ora_dir, paste0(ct, "_ORA_results.csv")),
            row.names = FALSE)

  all_ora_results[[ct]] <- list(dispMat = dispMat, txtMat = txtMat)

  # Per-celltype heatmap
  col_fun <- colorRamp2(c(0, 2.5, 5), viridis(3))
  direction_colors <- c(Up_in_Male = "steelblue", Up_in_Female = "indianred")

  meta_ct <- data.frame(
    row = rownames(dispMat),
    direction = ifelse(grepl("Up_in_Male", rownames(dispMat)),
                       "Up_in_Male", "Up_in_Female")
  )
  row_anno_ct <- rowAnnotation(
    direction = meta_ct$direction,
    col = list(direction = direction_colors),
    show_annotation_name = TRUE
  )

  pdf(file.path(ora_dir, paste0(ct, "_ORA_heatmap.pdf")), width = 10, height = 4)
  draw(Heatmap(
    matrix = dispMat, name = "Odds Ratio", col = col_fun, na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      label <- txtMat[i, j]
      parts <- unlist(strsplit(label, "\n"))
      stars <- parts[1]
      perc  <- ifelse(length(parts) > 1, parts[2], "")
      if (nzchar(stars)) grid.text(stars, x, y + unit(1.5, "mm"),
                                   gp = gpar(fontsize = 20, fontface = "bold"))
      if (nzchar(perc))  grid.text(perc, x, y - unit(1.5, "mm"),
                                   gp = gpar(fontsize = 6))
    },
    row_names_side = "left", column_names_rot = 45,
    row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10),
    left_annotation = row_anno_ct,
    column_title = paste0(info$label, ": AD-Only Sex DEGs in hdWGCNA Modules")
  ), heatmap_legend_side = "right")
  dev.off()

  cat(sprintf("    ORA heatmap saved (%d modules tested)\n", length(module_lists)))
}

# --- Combined ORA heatmap across all cell types ---
if (length(all_ora_results) > 0) {
  all_modules <- unique(unlist(lapply(all_ora_results, function(x) colnames(x$dispMat))))
  dispMat_all <- NULL
  txtMat_all  <- NULL

  for (ct in names(all_ora_results)) {
    dm <- all_ora_results[[ct]]$dispMat
    tm <- all_ora_results[[ct]]$txtMat
    dm_full <- matrix(NA, nrow = nrow(dm), ncol = length(all_modules))
    tm_full <- matrix("", nrow = nrow(tm), ncol = length(all_modules))
    colnames(dm_full) <- colnames(tm_full) <- all_modules
    rownames(dm_full) <- rownames(dm)
    rownames(tm_full) <- rownames(tm)
    shared_cols <- intersect(colnames(dm), all_modules)
    dm_full[, shared_cols] <- dm[, shared_cols]
    tm_full[, shared_cols] <- tm[, shared_cols]
    dispMat_all <- rbind(dispMat_all, dm_full)
    txtMat_all  <- rbind(txtMat_all, tm_full)
  }

  metadata_ora <- data.frame(
    row       = rownames(dispMat_all),
    celltype  = gsub("_Up_in_Male$|_Up_in_Female$", "", rownames(dispMat_all)),
    direction = ifelse(grepl("Up_in_Male", rownames(dispMat_all)),
                       "Up_in_Male", "Up_in_Female"),
    stringsAsFactors = FALSE
  )

  col_fun <- colorRamp2(c(0, 2.5, 5), viridis(3))
  ct_levels <- unique(metadata_ora$celltype)
  n_ct <- max(3, length(ct_levels))
  ct_colors <- setNames(brewer.pal(n_ct, "Set3")[seq_along(ct_levels)], ct_levels)
  direction_colors <- c(Up_in_Male = "steelblue", Up_in_Female = "indianred")

  row_anno_all <- rowAnnotation(
    celltype  = metadata_ora$celltype,
    direction = metadata_ora$direction,
    col = list(celltype = ct_colors, direction = direction_colors),
    show_annotation_name = TRUE
  )

  pdf(file.path(ora_dir, "combined_ORA_heatmap.pdf"), width = 14, height = 8)
  draw(Heatmap(
    matrix = dispMat_all, name = "Odds Ratio", col = col_fun, na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      label <- txtMat_all[i, j]
      parts <- unlist(strsplit(label, "\n"))
      stars <- parts[1]
      perc  <- ifelse(length(parts) > 1, parts[2], "")
      if (nzchar(stars)) grid.text(stars, x, y + unit(1.5, "mm"),
                                   gp = gpar(fontsize = 20, fontface = "bold"))
      if (nzchar(perc))  grid.text(perc, x, y - unit(1.5, "mm"),
                                   gp = gpar(fontsize = 6))
    },
    row_names_side = "left", column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10),
    left_annotation = row_anno_all,
    column_title = "AD-Only Sex DEGs Enrichment in hdWGCNA Modules (All Cell Types)"
  ), heatmap_legend_side = "right")
  dev.off()
  cat("\nCombined ORA heatmap saved.\n")
}

###############################################################################
# STEP 5: Compact bar chart — ME sex differences across all cell types
#         (inspired by Trait-Module Spearman correlation plot)
###############################################################################
cat("\n=== STEP 5: Compact bar chart of ME sex differences ===\n")

if (nrow(me_test_all) > 0) {
  # Use p < 0.05 (nominal) for the bar chart; note FDR in label
  sig_for_plot <- me_test_all %>% filter(p_value < 0.05)

  if (nrow(sig_for_plot) == 0) {
    cat("  No nominally significant ME differences to plot.\n")
    # Fall back to top results
    sig_for_plot <- me_test_all %>% arrange(p_value) %>% head(20)
    plot_title <- "ME Sex Differences by Cell Type (Top 20, Dementia)"
  } else {
    plot_title <- "ME Sex Differences (p<0.05, Dementia)"
  }

  # Add module annotation labels if available
  sig_for_plot$go_label <- ""
  if (!is.null(go_annot)) {
    for (i in 1:nrow(sig_for_plot)) {
      match_row <- go_annot %>%
        filter(ct_key == sig_for_plot$ct_key[i] &
               module == sig_for_plot$module[i])
      if (nrow(match_row) > 0) {
        sig_for_plot$go_label[i] <- match_row$term_short[1]
      }
    }
  }

  # Create y-axis labels: CellType-Module\nGO term
  sig_for_plot <- sig_for_plot %>%
    mutate(
      module_label = ifelse(
        nzchar(go_label),
        paste0(toupper(ct_key), "-", module, "\n", go_label),
        paste0(toupper(ct_key), "-", module)
      ),
      stars = case_when(
        fdr < 0.005 ~ "***",
        fdr < 0.01  ~ "**",
        fdr < 0.05  ~ "*",
        TRUE        ~ ""
      )
    )

  # Assign colors by cell type
  celltype_palette <- c(
    "Astrocyte"       = "#8DD3C7",
    "Excitatory"      = "#80B1D3",
    "Inhibitory"      = "#B3DE69",
    "Microglia"       = "#FDB462",
    "Oligodendrocyte" = "#FB8072"
  )

  # Order by cell type, then by diff_median
  sig_for_plot <- sig_for_plot %>%
    arrange(celltype, diff_median)
  sig_for_plot$module_label <- factor(sig_for_plot$module_label,
                                      levels = sig_for_plot$module_label)

  # Build cell type grouping annotation for right side
  ct_groups <- sig_for_plot %>%
    group_by(celltype) %>%
    summarise(ymin = min(as.numeric(module_label)) - 0.4,
              ymax = max(as.numeric(module_label)) + 0.4,
              ymid = mean(as.numeric(module_label)),
              .groups = "drop")

  p_bar <- ggplot(sig_for_plot,
                  aes(x = diff_median, y = module_label, fill = celltype)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(label = round(diff_median, 2)),
              hjust = ifelse(sig_for_plot$diff_median >= 0, -0.1, 1.1),
              size = 2.8) +
    scale_fill_manual(values = celltype_palette) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 7, lineheight = 0.8),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = plot_title,
      subtitle = "Median(Male) - Median(Female) ME values",
      x     = "Difference in Median ME (Male - Female)",
      y     = "",
      fill  = "Cell Type"
    )

  # Add cell type group separators (dashed lines between groups)
  if (nrow(ct_groups) > 1) {
    boundaries <- c()
    ct_order <- unique(sig_for_plot$celltype)
    for (i in 1:(length(ct_order) - 1)) {
      grp1 <- sig_for_plot %>% filter(celltype == ct_order[i])
      grp2 <- sig_for_plot %>% filter(celltype == ct_order[i + 1])
      boundary <- (max(as.numeric(grp1$module_label)) +
                   min(as.numeric(grp2$module_label))) / 2
      boundaries <- c(boundaries, boundary)
    }
    for (b in boundaries) {
      p_bar <- p_bar + geom_hline(yintercept = b, linetype = "dashed",
                                   color = "grey50", linewidth = 0.3)
    }
  }

  ggsave(file.path(plot_dir, "ME_sex_differences_barchart.pdf"),
         p_bar, width = 10, height = max(6, nrow(sig_for_plot) * 0.4))
  ggsave(file.path(plot_dir, "ME_sex_differences_barchart.png"),
         p_bar, width = 10, height = max(6, nrow(sig_for_plot) * 0.4), dpi = 300)
  cat("  Bar chart saved.\n")

  # Also make a version with ALL modules (not just significant)
  me_test_all_plot <- me_test_all %>%
    mutate(
      go_label = "",
      sig_label = case_when(
        fdr < 0.005 ~ "***",
        fdr < 0.01  ~ "**",
        fdr < 0.05  ~ "*",
        p_value < 0.05 ~ "+",
        TRUE        ~ ""
      )
    )

  if (!is.null(go_annot)) {
    for (i in 1:nrow(me_test_all_plot)) {
      match_row <- go_annot %>%
        filter(ct_key == me_test_all_plot$ct_key[i] &
               module == me_test_all_plot$module[i])
      if (nrow(match_row) > 0) {
        me_test_all_plot$go_label[i] <- match_row$term_short[1]
      }
    }
  }

  me_test_all_plot <- me_test_all_plot %>%
    mutate(
      module_label = ifelse(
        nzchar(go_label),
        paste0(toupper(ct_key), "-", module, "\n", go_label),
        paste0(toupper(ct_key), "-", module)
      )
    ) %>%
    arrange(celltype, diff_median)
  me_test_all_plot$module_label <- factor(me_test_all_plot$module_label,
                                          levels = me_test_all_plot$module_label)

  # Color: significant modules get cell type color, NS get grey
  me_test_all_plot <- me_test_all_plot %>%
    mutate(bar_color = ifelse(p_value < 0.05, celltype, "NS"))

  bar_palette <- c(celltype_palette, "NS" = "grey80")

  p_bar_all <- ggplot(me_test_all_plot,
                      aes(x = diff_median, y = module_label, fill = bar_color)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(label = sig_label),
              hjust = ifelse(me_test_all_plot$diff_median >= 0, -0.3, 1.3),
              size = 3, fontface = "bold") +
    scale_fill_manual(values = bar_palette,
                      breaks = names(celltype_palette),
                      labels = names(celltype_palette)) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 6, lineheight = 0.8),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "All Module Eigengene Sex Differences (Dementia)",
      subtitle = "Median(Male) - Median(Female) | + p<0.05, * FDR<0.05, ** FDR<0.01, *** FDR<0.005",
      x     = "Difference in Median ME (Male - Female)",
      y     = "",
      fill  = "Cell Type"
    )

  # Add cell type separators
  ct_order_all <- unique(me_test_all_plot$celltype)
  if (length(ct_order_all) > 1) {
    for (i in 1:(length(ct_order_all) - 1)) {
      grp1 <- me_test_all_plot %>% filter(celltype == ct_order_all[i])
      grp2 <- me_test_all_plot %>% filter(celltype == ct_order_all[i + 1])
      boundary <- (max(as.numeric(grp1$module_label)) +
                   min(as.numeric(grp2$module_label))) / 2
      p_bar_all <- p_bar_all + geom_hline(yintercept = boundary, linetype = "dashed",
                                           color = "grey50", linewidth = 0.3)
    }
  }

  ggsave(file.path(plot_dir, "ME_sex_differences_all_modules_barchart.pdf"),
         p_bar_all, width = 11, height = max(8, nrow(me_test_all_plot) * 0.35))
  ggsave(file.path(plot_dir, "ME_sex_differences_all_modules_barchart.png"),
         p_bar_all, width = 11, height = max(8, nrow(me_test_all_plot) * 0.35), dpi = 300)
  cat("  Full bar chart (all modules) saved.\n")
}

###############################################################################
# Summary
###############################################################################
cat("\n")
cat("=============================================\n")
cat("Pipeline v2 complete! Outputs:\n")
cat(sprintf("  Cleaned DEGs:         %s\n", clean_dir))
cat(sprintf("  Venn diagrams:        %s\n", venn_dir))
cat(sprintf("  AD-only DEGs:         %s\n", deg_dir))
cat(sprintf("  ME by sex:            %s\n", me_dir))
cat(sprintf("  Module enrichment:    %s\n", ora_dir))
cat(sprintf("  Compact bar charts:   %s\n", plot_dir))
cat(sprintf("  Venn summary:         %s\n", file.path(out_base, "venn_summary.csv")))
cat("=============================================\n")
