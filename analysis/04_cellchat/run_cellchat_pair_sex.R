#!/usr/bin/env Rscript
# ============================================================================
# CellChat Differential Analysis: Male vs. Female
#
# Description:
#   Generic, parameterized male-vs-female CellChat comparison, run once per
#   diagnosis group (Control, MCI, Dementia) against the corresponding pair
#   of pre-computed per-sex CellChat RDS objects (see the per-condition/sex
#   scripts in this directory). Merges the two objects, rebuilding
#   net/netP inference on either side if missing, then generates the full
#   comparison suite: compareInteractions, differential interaction circle
#   and heatmap plots, per-dataset circle plots, signaling-role scatter and
#   heatmaps (incoming/outgoing), rankNet info-flow bar charts, bubble
#   plots, functional/structural similarity embeddings, and per-pathway
#   circle/heatmap/chord plots for every shared signaling pathway.
#
#   object.list is built as list(female = ..., male = ...), so throughout
#   CellChat's standard diverging red/blue outputs (netVisual_diffInteraction,
#   netVisual_heatmap), RED = higher in male (dataset 2) and BLUE = higher in
#   female (dataset 1). rankNet uses an explicit categorical palette instead
#   (female = "#CC79A7" pink, male = "#0072B2" blue) -- that one encodes which
#   sex a bar belongs to, not direction of change.
#
#   This is the analysis underlying the manuscript's sex-stratified CellChat
#   supplementary figures (differential interaction heatmaps and rankNet
#   plots, Control/MCI/AD).
#
# Inputs:
#   --female_rds, --male_rds  paths to the per-sex CellChat RDS objects for
#                              one diagnosis group (from cellchat_ct_fe.R /
#                              cellchat_ct_male.R, cellchat_mci_fe.R /
#                              cellchat_mci_male.R, or the dementia/female
#                              equivalent)
#
# Outputs (written to --outdir, suffixed "_female_vs_male"):
#   int_compare, diff_interactions_circle, diff_interactions_heatmap,
#   circle_plot_per_dataset, cell_population_sending_receiving,
#   differential_signaling_scatter, ranknets, signalingRole_heatmap_out/in,
#   netVisual_bubble, netVisual_embeddingPairwise/ZoomIn + rankSimilarity
#   (functional and/or structural), per-pathway plots in
#   cellchat_pathway_plots/, and cellchat_female_vs_male_final.rds
#
# Usage:
#   Rscript run_cellchat_pair_sex.R \
#     --female_rds cellchat_by_ct_fe.rds --male_rds cellchat_by_ct_male.rds \
#     --outdir cellchat_diff_results_CT_female_vs_male
#
#   Run once per diagnosis group (CT / MCI / Dementia) pointing at that
#   group's female/male RDS pair.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(ComplexHeatmap)
  library(grid)        # unit()
  library(patchwork)   # wrap_plots
  library(ggplot2)
})

# ---------------- CLI ----------------
opt_list <- list(
  make_option("--female_rds", type="character", default="cellchat_by_mci_fe.rds",
              help="path to female CellChat RDS"),
  make_option("--male_rds", type="character", default="cellchat_by_mci_male.rds",
              help="path to male CellChat RDS"),
  make_option("--pos_dataset", type="character", default="male",
              help="DEG mapping positive dataset: 'male' or 'female'"),
  make_option("--outdir", type="character", default=".",
              help="output directory"),
  make_option("--min_cells", type="integer", default=10,
              help="min.cells for filterCommunication() when rebuilding"),
  make_option("--rebuild_if_missing", type="logical", default=TRUE,
              help="if TRUE, rebuild net/netP when missing so ALL plots can run"),
  make_option("--cache_rebuilt", type="logical", default=TRUE,
              help="if TRUE, save & reuse rebuilt per-dataset objects in outdir/cache"),
  make_option("--do_parallel", type="logical", default=FALSE,
              help="CellChat internal parallelization flag for netClustering; computeCommunProb uses future if you configured it externally"),
  make_option("--similarity_types", type="character", default="both",
              help="functional|structural|both")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(opt$pos_dataset %in% c("male", "female"))
stopifnot(opt$similarity_types %in% c("functional", "structural", "both"))

outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pathway_dir <- file.path(outdir, "cellchat_pathway_plots")
dir.create(pathway_dir, showWarnings = FALSE, recursive = TRUE)

cache_dir <- file.path(outdir, "cache_rebuilt_objects")
if (isTRUE(opt$cache_rebuilt)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

comparison_label <- "female_vs_male"

# ---------------- Logging helpers ----------------
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) message("[", ts(), "] ", paste0(..., collapse = ""))

# ---------------- PDF helper ----------------
safe_write_pdf <- function(fname, expr, width = NULL, height = NULL) {
  fname_full <- file.path(outdir, fname)
  if (!is.null(width) && !is.null(height)) pdf(fname_full, width = width, height = height) else pdf(fname_full)
  ok <- TRUE
  tryCatch({ force(expr) },
           error = function(e) { logmsg("SKIP (", basename(fname_full), "): ", conditionMessage(e)); ok <<- FALSE })
  try(dev.off(), silent = TRUE)
  invisible(ok)
}

out_pdf <- function(stem, ext = "pdf") sprintf("%s_%s.%s", stem, comparison_label, ext)

# ---------------- Validity helpers ----------------
is_valid_mat <- function(x, min_rows = 2, min_cols = 2) {
  if (is.null(x)) return(FALSE)
  m <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  d <- dim(m)
  if (is.null(d) || length(d) != 2) return(FALSE)
  isTRUE(d[1] >= min_rows && d[2] >= min_cols)
}

has_valid_net <- function(obj) {
  is_valid_mat(obj@net$count, 2, 2) || is_valid_mat(obj@net$weight, 2, 2)
}

has_2plus_groups <- function(obj) {
  lvl <- levels(obj@idents)
  !is.null(lvl) && length(lvl) >= 2
}

get_pathways <- function(obj) {
  if (!is.null(obj@netP$pathways)) return(obj@netP$pathways)
  if (!is.null(obj@netP$prob) && length(dimnames(obj@netP$prob)) >= 3) return(dimnames(obj@netP$prob)[[3]])
  character(0)
}

all_cluster_idxs <- function(obj) {
  lvl <- levels(obj@idents)
  if (!is.null(lvl) && length(lvl) > 0) return(seq_along(lvl))
  integer(0)
}

# ---------------- Rebuild pipeline (single dataset) ----------------
# This mirrors the *single-dataset* inference steps needed before comparison analysis.
rebuild_cellchat_if_needed <- function(obj, nm, min_cells = 10, rebuild_if_missing = TRUE) {

  # Tutorial recommendation: update objects made with older versions
  logmsg("updateCellChat(", nm, ")")
  obj <- updateCellChat(obj)

  if (!has_2plus_groups(obj)) {
    stop("Dataset '", nm, "' has <2 identity groups. ALL comparison plots (centrality/similarity) require >=2 groups.")
  }

  # If we already have valid net + netP, we can skip rebuild
  netp_ok <- (!is.null(obj@netP$prob) || !is.null(obj@netP$pathways))
  if (has_valid_net(obj) && netp_ok) {
    logmsg(nm, ": net and netP look OK; skipping rebuild.")
    return(obj)
  }

  if (!isTRUE(rebuild_if_missing)) {
    stop(nm, ": missing/invalid net/netP but --rebuild_if_missing=FALSE. Cannot make ALL plots.")
  }

  logmsg(nm, ": net/netP missing or invalid → rebuilding full inference pipeline.")
  logmsg(nm, ": steps: identifyOEGenes → identifyOEInteractions → computeCommunProb → filterCommunication → computeCommunProbPathway → aggregateNet")

  # These steps can be slow. We log before/after each to avoid “silent hangs”.
  t0 <- Sys.time()

  obj <- identifyOverExpressedGenes(obj)
  logmsg(nm, ": identifyOverExpressedGenes done (", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 2), " min)")

  obj <- identifyOverExpressedInteractions(obj)
  logmsg(nm, ": identifyOverExpressedInteractions done")

  obj <- computeCommunProb(obj)
  logmsg(nm, ": computeCommunProb done")

  obj <- filterCommunication(obj, min.cells = min_cells)
  logmsg(nm, ": filterCommunication done (min.cells=", min_cells, ")")

  obj <- computeCommunProbPathway(obj)
  logmsg(nm, ": computeCommunProbPathway done")

  obj <- aggregateNet(obj)
  logmsg(nm, ": aggregateNet done")

  # Validate rebuild
  if (!has_valid_net(obj)) {
    stop(nm, ": after rebuild, net$count/net$weight still not a valid >=2x2 matrix. Input object likely lacks usable underlying data.")
  }
  netp_ok2 <- (!is.null(obj@netP$prob) || !is.null(obj@netP$pathways))
  if (!netp_ok2) {
    stop(nm, ": after rebuild, netP still missing. computeNetSimilarityPairwise will not work.")
  }

  logmsg(nm, ": rebuild complete (total ", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 2), " min)")
  obj
}

compute_centrality_required <- function(obj, nm) {
  logmsg("netAnalysis_computeCentrality(", nm, ")")
  obj <- netAnalysis_computeCentrality(obj)
  obj
}

# ---------------- Run ----------------
options(future.globals.maxSize = 300000 * 1024^2)

logmsg(">>> Running: female vs male")
logmsg("female_rds: ", normalizePath(opt$female_rds))
logmsg("male_rds:   ", normalizePath(opt$male_rds))
logmsg("outdir:     ", normalizePath(outdir))

# Cache filenames (so reruns aren’t painful)
female_cache <- file.path(cache_dir, "female_rebuilt.rds")
male_cache   <- file.path(cache_dir, "male_rebuilt.rds")

# Load (and optionally reuse cached rebuilt objects)
if (isTRUE(opt$cache_rebuilt) && file.exists(female_cache)) {
  logmsg("Loading cached rebuilt female: ", female_cache)
  cellchat.female <- readRDS(female_cache)
} else {
  cellchat.female <- readRDS(opt$female_rds)
  cellchat.female <- rebuild_cellchat_if_needed(cellchat.female, "female",
                                                min_cells = opt$min_cells,
                                                rebuild_if_missing = opt$rebuild_if_missing)
  cellchat.female <- compute_centrality_required(cellchat.female, "female")
  if (isTRUE(opt$cache_rebuilt)) {
    saveRDS(cellchat.female, female_cache)
    logmsg("Saved rebuilt female cache: ", female_cache)
  }
}

if (isTRUE(opt$cache_rebuilt) && file.exists(male_cache)) {
  logmsg("Loading cached rebuilt male: ", male_cache)
  cellchat.male <- readRDS(male_cache)
} else {
  cellchat.male <- readRDS(opt$male_rds)
  cellchat.male <- rebuild_cellchat_if_needed(cellchat.male, "male",
                                              min_cells = opt$min_cells,
                                              rebuild_if_missing = opt$rebuild_if_missing)
  cellchat.male <- compute_centrality_required(cellchat.male, "male")
  if (isTRUE(opt$cache_rebuilt)) {
    saveRDS(cellchat.male, male_cache)
    logmsg("Saved rebuilt male cache: ", male_cache)
  }
}

# Build & merge (tutorial-style)
object.list <- list(female = cellchat.female, male = cellchat.male)
logmsg("Merging CellChat objects...")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, file.path(outdir, "checkpoint_merged_cellchat.rds"))
logmsg("Saved checkpoint: checkpoint_merged_cellchat.rds")

# Colors
PALETTE <- c(female = "#CC79A7", male = "#0072B2")
color.use.ranknet <- setNames(unname(PALETTE[c("female","male")]), c("female","male"))

# ---------------- Plots (your original set) ----------------

# compareInteractions
safe_write_pdf(
  out_pdf("int_compare"),
  {
    gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2),
                               color.use = unname(color.use.ranknet))
    gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2),
                               color.use = unname(color.use.ranknet), measure = "weight")
    print(gg1 + gg2)
  }
)

# diff interactions
safe_write_pdf(
  out_pdf("diff_interactions_circle"),
  {
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
  },
  width = 14, height = 12
)

# heatmap (merged object)
safe_write_pdf(
  out_pdf("diff_interactions_heatmap"),
  {
    g1 <- netVisual_heatmap(cellchat)
    g2 <- netVisual_heatmap(cellchat, measure = "weight")
    print(g1 + g2)
  },
  width = 11, height = 5.5
)

# per-dataset circles
safe_write_pdf(
  out_pdf("circle_plot_per_dataset"),
  {
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in seq_along(object.list)) {
      netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                       edge.weight.max = weight.max[2], edge.width.max = 12,
                       title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
  },
  width = 14, height = 10
)

# signaling role scatter (per dataset) — requires centrality already computed separately (we did)
safe_write_pdf(
  out_pdf("cell_population_sending_receiving"),
  {
    num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
    weight.MinMax <- range(num.link)
    gg <- lapply(seq_along(object.list), function(i)
      netAnalysis_signalingRole_scatter(object.list[[i]],
                                       title = names(object.list)[i],
                                       weight.MinMax = weight.MinMax))
    print(wrap_plots(plots = gg))
  },
  width = 11, height = 6
)

# signaling changes scatter (centrality-dependent)
safe_write_pdf(
  out_pdf("differential_signaling_scatter"),
  {
    plots <- list()
    idents_try <- c("Microglia","Astrocyte")
    for (idn in idents_try) {
      plots[[length(plots)+1]] <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = idn,
                                                                      signaling.exclude = c("OPC"))
    }
    print(wrap_plots(plots = plots))
  }
)

# rankNet
safe_write_pdf(
  out_pdf("ranknets"),
  {
    gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight",
                   stacked = TRUE,  do.stat = TRUE, color.use = color.use.ranknet)
    gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight",
                   stacked = FALSE, do.stat = TRUE, color.use = color.use.ranknet)
    print(gg1 + gg2)
  }
)

# Signaling role heatmaps (union of available)
pw1 <- get_pathways(object.list$female)
pw2 <- get_pathways(object.list$male)
pathway.union <- union(pw1, pw2)

safe_write_pdf(
  out_pdf("signalingRole_heatmap_out"),
  {
    ht1 <- netAnalysis_signalingRole_heatmap(object.list$female, pattern = "outgoing",
                                             signaling = pathway.union, title = "female", width = 8, height = 10)
    ht2 <- netAnalysis_signalingRole_heatmap(object.list$male,   pattern = "outgoing",
                                             signaling = pathway.union, title = "male", width = 8, height = 10)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  },
  width = 10, height = 12
)

safe_write_pdf(
  out_pdf("signalingRole_heatmap_in"),
  {
    ht1 <- netAnalysis_signalingRole_heatmap(object.list$female, pattern = "incoming",
                                             signaling = pathway.union, title = "female", width = 8, height = 10,
                                             color.heatmap = "GnBu")
    ht2 <- netAnalysis_signalingRole_heatmap(object.list$male,   pattern = "incoming",
                                             signaling = pathway.union, title = "male", width = 8, height = 10,
                                             color.heatmap = "GnBu")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  },
  width = 10, height = 10
)

# Bubble
safe_write_pdf(
  out_pdf("netVisual_bubble"),
  {
    g <- netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
    g <- g + theme(axis.text.x = element_text(size = 8))
    print(g)
  },
  width = 30, height = 10
)

# ---------------- Similarity analysis (tutorial) ----------------
# We run functional + structural (or whichever you requested).
run_similarity <- function(type) {
  logmsg("Running computeNetSimilarityPairwise(type='", type, "')")
  cellchat2 <- computeNetSimilarityPairwise(cellchat, type = type)
  logmsg("Running netEmbedding(type='", type, "')")
  cellchat2 <- netEmbedding(cellchat2, type = type)
  logmsg("Running netClustering(type='", type, "')")
  cellchat2 <- netClustering(cellchat2, type = type, do.parallel = opt$do_parallel)

  safe_write_pdf(out_pdf(paste0("netVisual_embeddingPairwise_", type)),
                 netVisual_embeddingPairwise(cellchat2, type = type, label.size = 3.5))
  safe_write_pdf(out_pdf(paste0("netVisual_embeddingZoomIn_", type)),
                 netVisual_embeddingZoomIn(cellchat2, type = type))
  safe_write_pdf(out_pdf(paste0("rankSimilarity_", type)),
                 rankSimilarity(cellchat2, type = type))

  saveRDS(cellchat2, file.path(outdir, paste0("checkpoint_cellchat_similarity_", type, ".rds")))
  logmsg("Saved checkpoint: checkpoint_cellchat_similarity_", type, ".rds")
  invisible(cellchat2)
}

if (opt$similarity_types %in% c("functional","both")) run_similarity("functional")
if (opt$similarity_types %in% c("structural","both")) run_similarity("structural")

# ---------------- Dynamic pathway loop (intersection) ----------------
pathway.list <- intersect(get_pathways(object.list$female), get_pathways(object.list$male))
logmsg("Dynamic pathway intersection size: ", length(pathway.list))

writeLines(sort(get_pathways(object.list$female)), file.path(pathway_dir, "pathways_female.txt"))
writeLines(sort(get_pathways(object.list$male)),   file.path(pathway_dir, "pathways_male.txt"))
writeLines(sort(pathway.list),                     file.path(pathway_dir, "pathways_intersection.txt"))

for (pathway in pathway.list) {

  safe_write_pdf(
    file.path("cellchat_pathway_plots", paste0("circle_shared_signaling_", pathway, ".pdf")),
    {
      weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
      par(mfrow = c(1, 2), xpd = TRUE)
      for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle",
                            edge.weight.max = weight.max[1], edge.width.max = 10,
                            signaling.name = paste(pathway, names(object.list)[i]))
      }
    },
    width = 10, height = 6
  )

  safe_write_pdf(
    file.path("cellchat_pathway_plots", paste0("heatmap_signaling_", pathway, ".pdf")),
    {
      ht <- vector("list", length(object.list))
      for (i in seq_along(object.list)) {
        ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathway,
                                     color.heatmap = "Reds",
                                     title.name = paste(pathway, "signaling", names(object.list)[i]))
      }
      ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    },
    width = 10, height = 6
  )

  safe_write_pdf(
    file.path("cellchat_pathway_plots", paste0("chord_signaling_", pathway, ".pdf")),
    {
      par(mfrow = c(1, 2), xpd = TRUE)
      for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "chord",
                            signaling.name = paste(pathway, names(object.list)[i]))
      }
    },
    width = 10, height = 6
  )
}

# ---------------- Save final merged object ----------------
saveRDS(cellchat, file.path(outdir, "cellchat_female_vs_male_final.rds"))
logmsg("Done: female vs male → ", outdir)
