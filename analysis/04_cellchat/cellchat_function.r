# ============================================================================
# CellChat Analysis Function
# ============================================================================
# Description:
#   Reusable function to run the full CellChat pipeline on an h5ad file.
#   Loads an anndata object, filters unknown cells, normalizes counts,
#   creates a CellChat object, identifies overexpressed genes/interactions,
#   computes communication probabilities, aggregates the network, generates
#   standard visualizations, and saves the CellChat RDS object.
#
# Inputs:
#   - h5ad_file: path to an .h5ad file with predictions_scanvi annotations
#   - output_prefix: (unused, kept for API compatibility)
#
# Outputs:
#   - Per-dataset output directory containing:
#       <dataset>_cellchat_log.txt
#       <dataset>_cellchat.rds
#       PDF plots: num_int.pdf, int_weights.pdf, CXCL_*.pdf,
#                  network_centrality.pdf
#
# Usage:
#   source("cellchat_function.r")
#   cellchat <- run_cellchat("my_data.h5ad", "prefix")
# ============================================================================

library(CellChat)
library(patchwork)
library(anndata)
library(Matrix)
library(reticulate)

# --- Main CellChat pipeline function ---

run_cellchat <- function(h5ad_file, output_prefix) {
  options(stringsAsFactors = FALSE)
  options(future.globals.maxSize = 900000 * 1024^2)

  # Create output directory based on the h5ad filename
  dataset_name <- sub("\\.h5ad$", "", basename(h5ad_file))
  output_dir <- file.path(getwd(), paste0(dataset_name, "_cellchat_outputs"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- Start logging ---

  log_file <- file.path(output_dir, paste0(dataset_name, "_cellchat_log.txt"))
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "message")

  message("====================================================")
  message("Starting CellChat analysis for: ", h5ad_file)
  ptm <- Sys.time()

  # --- Load and filter data ---

  message("Loading the .h5ad file...")
  data <- read_h5ad(h5ad_file)

  message("Filtering cells labeled as 'unknown'...")
  data_filtered <- data[data$obs$predictions_scanvi != "unknown", ]

  message("Filtered cells summary:")
  print(table(data_filtered$obs$predictions_scanvi))

  # --- Normalize count data ---

  message("Accessing and normalizing the count data matrix...")
  counts <- t(as.matrix(data_filtered$X))
  library.size <- Matrix::colSums(counts)
  data_input <- as(log1p(Matrix::t(Matrix::t(counts) / library.size) * 10000), "dgCMatrix")

  # --- Access metadata ---

  message("Accessing metadata...")
  meta <- data_filtered$obs
  meta$labels <- meta[["clusters"]]

  # --- Create CellChat object ---

  message("Creating the CellChat object...")
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "predictions_scanvi")

  message("Setting the signaling database...")
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)

  # --- Run pipeline ---

  message("Running parallel processes...")
  future::plan("multisession", workers = 4)

  message("Identifying overexpressed genes...")
  cellchat <- identifyOverExpressedGenes(cellchat)

  message("Identifying overexpressed interactions...")
  cellchat <- identifyOverExpressedInteractions(cellchat)

  message("Computing communication probabilities...")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  execution_time <- Sys.time() - ptm
  message("Execution time (seconds): ", as.numeric(execution_time, units = "secs"))

  # --- Helper for saving PDFs ---

  save_pdf <- function(filename) {
    pdf(file.path(output_dir, filename))
    on.exit(dev.off())
  }

  # --- Generate visualizations ---

  message("Generating visualizations...")

  save_pdf("num_int.pdf")
  netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")

  save_pdf("int_weights.pdf")
  netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction weights/strength")

  # CXCL pathway visualizations
  pathways.show <- c("CXCL")
  vertex.receiver <- seq(1, 4)

  save_pdf("CXCL_hierarchy.pdf")
  netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)

  save_pdf("CXCL_circle.pdf")
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

  save_pdf("CXCL_chord.pdf")
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

  save_pdf("CXCL_heatmap.pdf")
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

  # Network centrality
  message("Computing and visualizing network centrality...")
  save_pdf("network_centrality.pdf")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                    width = 8, height = 2.5, font.size = 10)

  # --- Save CellChat object ---

  message("Saving the CellChat object...")
  saveRDS(cellchat, file = file.path(output_dir, paste0(dataset_name, "_cellchat.rds")))

  message("CellChat analysis completed for: ", h5ad_file)

  # --- Close logging ---

  sink(type = "message")
  close(log_conn)

  return(cellchat)
}
