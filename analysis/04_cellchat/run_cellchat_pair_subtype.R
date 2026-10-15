#!/usr/bin/env Rscript
# ==============================================================================
# run_cellchat_pair_subtype.R
# ==============================================================================
# Generic, parameterized pairwise CellChat differential-comparison script for
# one AD subtype vs. Control (or any two pre-built CellChat .rds objects).
# Run once per subtype to produce the AD{n}-vs-Control comparison suite behind
# Fig. 8b / Suppl. Fig. 6a: compareInteractions, differential interaction
# circle/heatmap/chord plots, signaling-role scatter/heatmaps, rankNet
# info-flow bar charts, bubble plots, functional/structural similarity
# embeddings, and a full per-pathway (including ANGPTL, RA, TGF-beta)
# circle/heatmap/chord sweep.
#
# Usage (run once per subtype):
#   Rscript run_cellchat_pair_subtype.R \
#     --a_label ad1 --a_rds cellchat_by_ad1.rds \
#     --b_label ct  --b_rds cellchat_by_control.rds \
#     --outdir ad1_vs_ct
#   ... repeat for ad2, ad3, ad4 ...
#
# See cellchat_subtype_panel_fig8b.R for the combined Control+AD1-4 rankNet
# panel built from the 5 CellChat .rds objects these comparisons share.
# ==============================================================================

suppressPackageStartupMessages({
  library(hdf5r)
  library(CellChat)
  library(patchwork)
  library(anndata)
  library(ComplexHeatmap)
  library(optparse)
  library(grid)   # for unit()
})

# -------- CLI --------
opt_list <- list(
  make_option("--a_label", type="character", help="dataset A label (e.g., ad1)"),
  make_option("--a_rds",   type="character", help="dataset A .rds file"),
  make_option("--b_label", type="character", help="dataset B label (e.g., ct)"),
  make_option("--b_rds",   type="character", help="dataset B .rds file"),
  make_option("--outdir",  type="character", default=NULL, help="output directory (default: <a_label>_vs_<b_label>)")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$a_label), !is.null(opt$a_rds), !is.null(opt$b_label), !is.null(opt$b_rds))

a_lab <- opt$a_label
b_lab <- opt$b_label
a_rds <- opt$a_rds
b_rds <- opt$b_rds
outdir <- if (is.null(opt$outdir)) sprintf("%s_vs_%s", a_lab, b_lab) else opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------- Load data --------
cellchat.ct <- readRDS(b_rds)
cellchat.ad <- readRDS(a_rds)

object.list <- list()
object.list[[b_lab]] <- cellchat.ct
object.list[[a_lab]] <- cellchat.ad

# -------- Ensure idents & dimnames exist so scPalette/netVisual_circle work --------
for (i in seq_along(object.list)) {
  obj <- object.list[[i]]

  if (is.null(obj@idents) || length(obj@idents) == 0) {
    rn <- rownames(obj@net$count)
    if (is.null(rn)) stop("CellChat object has no idents and net$count has no rownames.")
    obj@idents <- factor(rn, levels = rn)
  } else if (!is.factor(obj@idents)) {
    obj@idents <- factor(obj@idents)
  }

  lv <- levels(obj@idents)
  if (is.null(rownames(obj@net$count)) || is.null(colnames(obj@net$count))) {
    rownames(obj@net$count)  <- lv
    colnames(obj@net$count)  <- lv
  }
  if (!is.null(obj@net$weight)) {
    if (is.null(rownames(obj@net$weight)) || is.null(colnames(obj@net$weight))) {
      rownames(obj@net$weight) <- lv
      colnames(obj@net$weight) <- lv
    }
  }

  object.list[[i]] <- obj
}

all_groups <- unique(unlist(lapply(object.list, function(x) levels(x@idents))))
pal <- tryCatch(CellChat::scPalette(length(all_groups)), error = function(e) grDevices::rainbow(length(all_groups)))
names(pal) <- all_groups

for (i in seq_along(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

color.use <- setNames(c("#20B2AA", "#B39DDB"), c(b_lab, a_lab))

f <- function(name) file.path(outdir, name)

# -------- Overview plots --------
pdf(f(sprintf("int_%s_%s.pdf", a_lab, b_lab)))
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), color.use = c("#4C9DAF", "#A593E0"))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), color.use = c("#4C9DAF", "#A593E0"), measure = "weight")
gg1 + gg2
dev.off()

pdf(f(sprintf("diff_interactions_%s_%s.pdf", a_lab, b_lab)), width = 14, height = 12)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

pdf(f(sprintf("diff_interactions_%s_%s_heatmap.pdf", a_lab, b_lab)), width = 11, height = 5.5)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()

pdf(f(sprintf("diff_interactions_%s_%s_circle_plot.pdf", a_lab, b_lab)), width = 14, height = 10)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in seq_along(object.list)) {
  grp_cols <- pal[levels(object.list[[i]]@idents)]
  netVisual_circle(object.list[[i]]@net$count,
    weight.scale = TRUE, label.edge = FALSE,
    edge.weight.max = weight.max[2], edge.width.max = 12,
    color.use = grp_cols,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}
dev.off()

# -------- Rank / similarity --------
pdf(f(sprintf("%s_%s_ranknets.pdf", a_lab, b_lab)))
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE,  do.stat = TRUE, color.use = color.use)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE, color.use = color.use)
gg1 + gg2
dev.off()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
options(future.globals.maxSize = 300000 * 1024^2)
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)

pdf(f(sprintf("%s_%s_netVisual_embeddingPairwise.pdf", a_lab, b_lab)))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

# -------- Role heatmaps --------
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

pdf(f(sprintf("%s_%s_netAnalysis_signalingRole_heatmap.pdf", a_lab, b_lab)), width = 10, height = 12)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],   pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i],   width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# -------- Per-pathway sweep (covers ANGPTL, RA, TGF-beta among all shared pathways) --------
availA <- tryCatch(dimnames(object.list[[1]]@netP$prob)[[3]], error = function(e) character())
availB <- tryCatch(dimnames(object.list[[2]]@netP$prob)[[3]], error = function(e) character())
pathway.list <- sort(unique(c(availA, availB)))

loop.dir <- file.path(outdir, "cellchat_pathway_plots")
dir.create(loop.dir, showWarnings = FALSE, recursive = TRUE)

for (pathway in pathway.list) {
  message("Processing pathway: ", pathway)

  pathway.exists <- sapply(object.list, function(obj)
    pathway %in% tryCatch(dimnames(obj@netP$prob)[[3]], error=function(e) character()))
  if (!all(pathway.exists)) {
    warning("Skipping pathway '", pathway, "' because it's missing in one or more datasets.")
    next
  }

  pdf(file.path(loop.dir, paste0("circle_shared_signaling_", pathway, ".pdf")), width = 10, height = 6)
  weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in seq_along(object.list)) {
    grp_cols <- pal[levels(object.list[[i]]@idents)]
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle",
                        edge.weight.max = weight.max[1], edge.width.max = 10,
                        color.use = grp_cols,
                        signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()

  pdf(file.path(loop.dir, paste0("heatmap_signaling_", pathway, ".pdf")), width = 10, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  ht <- list()
  for (i in seq_along(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathway, color.heatmap = "Reds",
                                 title.name = paste(pathway, "signaling", names(object.list)[i]))
  }
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()

  pdf(file.path(loop.dir, paste0("chord_signaling_", pathway, ".pdf")), width = 10, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in seq_along(object.list)) {
    grp_cols <- pal[levels(object.list[[i]]@idents)]
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "chord",
                        color.use = grp_cols,
                        signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
}

save(object.list, file = f(sprintf("cellchat_differential_object.list_%s_%s.RData", a_lab, b_lab)))
save(cellchat,    file = f(sprintf("cellchat_differential_merged_%s_%s.RData", a_lab, b_lab)))
