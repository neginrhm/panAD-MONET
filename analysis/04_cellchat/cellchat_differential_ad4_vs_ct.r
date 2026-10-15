# ============================================================================
# CellChat Differential Analysis: AD Subtype 4 vs Control
# ============================================================================
# Description:
#   Compares cell-cell communication between AD subtype 4 and control using
#   CellChat. Builds CellChat objects from h5ad files (control, AD4), runs
#   the full CellChat pipeline (normalization, overexpressed genes/interactions,
#   communication probability), merges objects, and generates differential
#   interaction plots, signaling role analysis, functional similarity,
#   ranked networks, up/down-regulated LR pairs, FGF pathway visualizations,
#   and batch per-pathway plots.
#
# Inputs:
#   - unique_ons_control_adata_final.h5ad (control h5ad)
#   - unique_ons_4_adata_final.h5ad (AD subtype 4 h5ad)
#   - Pre-computed CellChat RDS objects (dementia, MCI, control)
#
# Outputs:
#   - Multiple PDF plots comparing AD4 vs control signaling
#   - RData files: cellchat_differential_object.list_AD4_ct.RData,
#                  cellchat_differential_merged_AD4_ct.RData
#   - Per-pathway PDFs in cellchat_pathway_plots/ directory
#
# Usage:
#   Rscript cellchat_differential_sexes.r
# ============================================================================

library(hdf5r)
library(CellChat)
library(patchwork)
library(anndata)
library(ComplexHeatmap)
library(grid)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 300000 * 1024^2)

BASE_DIR <- getwd()  # Set to project root before running

# --- Load and process control h5ad ---

ptm <- Sys.time()

ct <- read_h5ad("unique_ons_control_adata_final.h5ad")
unique(ct$dx_or_cluster)
ct <- ct[ct$obs$predictions_scanvi != "unknown", ]

counts <- t(as.matrix(ct$X))
library.size <- Matrix::colSums(counts)
data.input <- as(log1p(Matrix::t(Matrix::t(counts) / library.size) * 10000), "dgCMatrix")
meta <- ct$obs
meta$labels <- meta[["dx_or_cluster"]]

cellchat_ct <- createCellChat(object = data.input, meta = meta, group.by = "predictions_scanvi")
CellChatDB <- CellChatDB.human
cellchat_ct@DB <- CellChatDB
cellchat_ct <- subsetData(cellchat_ct)

future::plan("multisession", workers = 4)

ptm <- Sys.time()
cellchat_ct <- computeCommunProb(cellchat_ct, type = "triMean")
cellchat_ct <- filterCommunication(cellchat_ct, min.cells = 5)
cellchat_ct <- computeCommunProbPathway(cellchat_ct)

cellchat_ct <- identifyOverExpressedGenes(cellchat_ct)
cellchat_ct <- identifyOverExpressedInteractions(cellchat_ct)

save(cellchat_ct, file = "cellchat_control_object.RData")

# --- Load and process AD4 h5ad ---

AD4 <- read_h5ad("unique_ons_4_adata_final.h5ad")
AD4_filtered <- AD4[AD4$obs$predictions_scanvi != "unknown", ]
table(AD4_filtered$obs$predictions_scanvi)

counts2 <- t(as.matrix(AD4_filtered$X))
library.size2 <- Matrix::colSums(counts2)
data.input2 <- as(log1p(Matrix::t(Matrix::t(counts2) / library.size2) * 10000), "dgCMatrix")
meta2 <- AD4_filtered$obs

cellchat_AD4 <- createCellChat(object = data.input2, meta = meta2, group.by = "predictions_scanvi")
cellchat_AD4@DB <- CellChatDB
cellchat_AD4 <- subsetData(cellchat_AD4)

cellchat_AD4 <- computeCommunProb(cellchat_AD4, type = "triMean")
cellchat_AD4 <- filterCommunication(cellchat_AD4, min.cells = 5)
cellchat_AD4 <- computeCommunProbPathway(cellchat_AD4)

options(future.globals.maxSize = 300000 * 1024^2)

cellchat_AD4 <- identifyOverExpressedGenes(cellchat_AD4)
cellchat_AD4 <- identifyOverExpressedInteractions(cellchat_AD4)

# --- Merge and time ---

object.list <- list(dem = cellchat_AD4, ct = cellchat_ct)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
execution.time <- Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# --- Load pre-computed RDS and re-merge for downstream analysis ---

cellchat.AD4 <- readRDS(file.path(BASE_DIR, "CellChat/cellchat_dementia_Redo_scanvi_predictions/cellchat_dementia_zellekonverter.rds"))
cellchat.mci <- readRDS(file.path(BASE_DIR, "cellchat_mci_zellekonverter.rds"))
cellchat.ct <- readRDS(file.path(BASE_DIR, "CellChat/cellchat_control_redo/cellchat_control_zellekonverter.rds"))

object.list <- list(ct = cellchat.ct, AD4 = cellchat.AD4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

execution.time <- Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# --- Compare total interactions ---

pdf("int_AD4_ct_reddo.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), color.use = c("#D1495B", "#A593E0"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), color.use = c("#D1495B", "#A593E0"), measure = "weight")
gg1 + gg2
dev.off()

# --- Differential interaction network ---

pdf("diff_interactions_AD4_ct_reddo.pdf", width = 14, height = 12)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

# --- Differential interaction heatmap ---

pdf("diff_interactions_AD4_ct_heatmap_reddo.pdf", width = 11, height = 5.5)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()

# --- Circle plot: interaction counts ---

pdf("diff_interactions_AD4_ct_circle_plot_reddo.pdf", width = 14, height = 10)
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = F,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# --- Circle plot: coarse cell types ---

group.cellType <- c(rep("Microglia", 4), rep("Excitatory", 4), rep("Inhibitory", 4))
group.cellType <- factor(group.cellType, levels = c("Microglia", "Excitatory", "Inhibitory"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf("diff_interactions_AD4_ct_circle_plot_coarse_celltypes_reddo.pdf", width = 14, height = 10)
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"),
                           attribute = c("idents", "count", "count.merged"))
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge = T,
                   edge.weight.max = weight.max[3], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# --- Differential number of interactions (merged) ---

pdf("diff_interactions_AD4_ct_circle_plot_differential_num_reddo.pdf", width = 14, height = 10)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()

# --- Signaling role scatter ---

for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}

pdf("cell_population_sending_receiving_AD4_ct_reddo.pdf", width = 11, height = 6)
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
               title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()

# --- Differential signaling scatter for specific cell types ---

pdf("AD4_ct_differential_signaling_scatter_reddo.pdf")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Microglia", signaling.exclude = "OPC")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astrocyte", signaling.exclude = c("OPC"))
patchwork::wrap_plots(plots = list(gg1, gg2))
dev.off()

# --- Functional similarity embedding ---

ptm <- Sys.time()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")

options(future.globals.maxSize = 300000 * 1024^2)
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)

pdf("AD4_ct_netVisual_embeddingPairwise.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf("AD4_ct_zoomin_netVisual_embeddingPairwise.pdf")
netVisual_embeddingZoomIn(cellchat, type = "functional")
dev.off()

pdf("AD4_ct_rankSimilarityc.pdf")
rankSimilarity(cellchat, type = "functional")
dev.off()

# --- Ranked signaling networks ---

pdf("AD4_ct_ranknets.pdf")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight",
               sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight",
               sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

# --- Signaling role heatmaps ---

pdf("AD4_ct_netAnalysis_signalingRole_heatmap.pdf", width = 10, height = 12)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing",
         signaling = pathway.union, title = names(object.list)[i], width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing",
         signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf("AD4_ct_netAnalysis_signalingRole_heatmap2.pdf", width = 10, height = 10)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming",
         signaling = pathway.union, title = names(object.list)[i],
         width = 8, height = 10, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming",
         signaling = pathway.union, title = names(object.list)[i + 1],
         width = 8, height = 10, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# --- Bubble plots ---

pdf("AD4_ct_netVisual_bubblee.pdf", width = 10, height = 10)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),
                 comparison = c(1, 2), angle.x = 45)
dev.off()

pdf("AD4_ct_netVisual_bubblee2.pdf", width = 10, height = 10)
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),
         comparison = c(1, 2), max.dataset = 2,
         title.name = "Increased signaling in ct", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),
         comparison = c(1, 2), max.dataset = 1,
         title.name = "Decreased signaling in ct", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()

# --- Differential expression and LR pair extraction ---

pos.dataset <- "AD4"
features.name <- paste0(pos.dataset, ".merged")

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
              pos.dataset = pos.dataset, features.name = features.name,
              only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,
              thresh.p = 0.05, group.DE.combined = FALSE)

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
              pos.dataset = pos.dataset, features.name = features.name,
              only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05)

net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "AD4",
            ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "ct",
              ligand.logFC = -0.1, receptor.logFC = -0.1)

# --- Save objects ---

save(object.list, file = "cellchat_differential_object.list_AD4_ct.RData")
save(cellchat, file = "cellchat_differential_merged_AD4_ct.RData")

# --- Up/down-regulated LR pair bubble plots ---

pdf("up_down_regulated_in_AD4_up_heatmap.pdf", width = 40, height = 15)
pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, font.size = 15, pairLR.use = pairLR.use.up,
         comparison = c(1, 2), angle.x = 45, remove.isolate = T,
         title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down <- net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
         comparison = c(1, 2), angle.x = 90, remove.isolate = T,
         title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1
dev.off()

pdf("up_down_regulated_in_AD4_down_heatmap.pdf", width = 40, height = 15)
pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, font.size = 15, pairLR.use = pairLR.use.up,
         comparison = c(1, 2), angle.x = 45, remove.isolate = T,
         title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down <- net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
         comparison = c(1, 2), angle.x = 45, remove.isolate = T,
         title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg2
dev.off()

# --- Up/down-regulated chord diagrams ---

pdf("up_down_regulated_in_AD4_chord.pdf", width = 14, height = 14)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11),
  slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5,
  title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11),
  slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5,
  title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

# --- Enrichment word clouds ---

pdf("up_down_regulated_in_AD4_words_down_AD4.pdf")
computeEnrichmentScore(net.down, species = 'human')
dev.off()

pdf("up_down_regulated_in_AD4_words_net_up_AD4.pdf")
computeEnrichmentScore(net.up, species = 'human')
dev.off()

# --- FGF pathway: circle, heatmap, chord ---

pdf("AD4_ct_fgf_circle_plot.pdf", width = 12, height = 12)
pathways.show <- c("FGF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle",
    edge.weight.max = weight.max[1], edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf("AD4_ct_fgf_heatmap.pdf", width = 12, height = 12)
pathways.show <- c("FGF")
par(mfrow = c(1, 2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
    color.heatmap = "Reds",
    title.name = paste(pathways.show, "signaling ", names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

pdf("AD4_ct_fgf_chord.pdf", width = 12, height = 12)
pathways.show <- c("FGF")
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",
    signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# --- FGF gene expression distribution ---

pdf("gene_exp_dist_FGF.pdf")
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("ct", "AD4"))
plotGeneExpression(cellchat, signaling = "FGF", split.by = "datasets", colors.ggplot = T)
dev.off()

# --- FGF hierarchy signaling circle ---

pdf("gene_exp_hierarchy_signaling_FGF.pdf", width = 10, height = 6)
pathways.show <- c("FGF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle",
    edge.weight.max = weight.max[1], edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# --- FGF heatmap (alternate dimensions) ---

pdf("AD4_ct_heatmap_signaling_FGF.pdf", width = 10, height = 6)
pathways.show <- c("FGF")
par(mfrow = c(1, 2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
    color.heatmap = "Reds",
    title.name = paste(pathways.show, "signaling ", names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# --- FGF chord (alternate dimensions) ---

pdf("AD4_ct_chord_signaling_FGF.pdf", width = 10, height = 6)
pathways.show <- c("FGF")
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",
    signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# --- Batch pathway visualization loop ---

pathway.list <- c(
  "GRN", "GALECTIN", "MHC-II", "HGF", "PECAM1", "FLRT", "EGF", "PERIOSTIN",
  "UNC5", "RA", "RELN", "FN1", "CD45", "SPP1", "TAC", "Netrin", "CNTN",
  "TENASCIN", "TGFb", "BMP", "CD39", "Cholesterol", "PSAP", "SLITRK", "GAS",
  "DHEA", "PROS", "SEMA4", "THBS", "MAG", "CDH", "LAMININ", "COLLAGEN", "FGF",
  "PDGF", "MPZ", "ANGPTL", "IGF", "ANGPT", "SLIT", "NGL", "JAM", "CLDN",
  "MHC-I", "NRG", "PTPRM", "PTN", "COMPLEMENT", "GAP", "CADM"
)

if (is.null(names(object.list))) {
  names(object.list) <- paste0("dataset", seq_along(object.list))
}

output.dir <- "cellchat_pathway_plots"
dir.create(output.dir, showWarnings = FALSE)

for (pathway in pathway.list) {
  message("Processing pathway: ", pathway)

  pathway.exists <- sapply(object.list, function(obj) pathway %in% dimnames(obj@netP$prob)[[3]])
  if (!all(pathway.exists)) {
    warning("Skipping pathway '", pathway, "' because it's missing in one or more datasets.")
    next
  }

  # Shared weight circle plot
  pdf(file.path(output.dir, paste0("circle_shared_signaling_", pathway, ".pdf")), width = 10, height = 6)
  weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle",
      edge.weight.max = weight.max[1], edge.width.max = 10,
      signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()

  # Individual weight circle plot
  pdf(file.path(output.dir, paste0("circle_individual_signaling_", pathway, ".pdf")), width = 10, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle",
      edge.weight.max = weight.max[i], edge.width.max = 10,
      signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()

  # Heatmap
  pdf(file.path(output.dir, paste0("heatmap_signaling_", pathway, ".pdf")), width = 10, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathway,
      color.heatmap = "Reds",
      title.name = paste(pathway, "signaling", names(object.list)[i]))
  }
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()

  # Chord diagram
  pdf(file.path(output.dir, paste0("chord_signaling_", pathway, ".pdf")), width = 10, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "chord",
      signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
}
