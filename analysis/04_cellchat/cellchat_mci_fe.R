# ============================================================================
# CellChat analysis: MCI Female samples
#
# Description: Run CellChat cell-cell communication inference on MCI
#              female snRNA-seq data grouped by scANVI-predicted cell types.
#
# Inputs:  h5ad/13adata_mci_female.h5ad (AnnData with predictions_scanvi labels)
# Outputs: cellchat_by_mci_fe.rds       (CellChat object)
#          celltype_global_circlesmci_fe.pdf (circle plot of interactions)
#
# Usage: Rscript cellchat_mci_fe.R
# ============================================================================

library(CellChat)
library(reticulate)
library(Matrix)
library(future)

# --- Configuration ---
BASE_DIR <- getwd()
options(future.globals.maxSize = 300000 * 1024^2)  # 100 GB

# --- Load data ---
ad <- reticulate::import("anndata", convert = FALSE)
obj <- ad$read_h5ad(file.path(BASE_DIR, "h5ad", "13adata_mci_female.h5ad"))

obj <- obj[obj$obs$predictions_scanvi != "unknown", ]

X   <- reticulate::py_to_r(obj$X)
mat <- if (inherits(X, c("dgCMatrix", "dgTMatrix"))) Matrix::t(X) else t(as.matrix(X))
var <- reticulate::py_to_r(obj$var)
obs <- reticulate::py_to_r(obj$obs)
rownames(mat) <- rownames(var)
colnames(mat) <- rownames(obs)

# Make a clean factor for cell types
stopifnot("predictions_scanvi" %in% colnames(obs))
obs$predictions_scanvi <- factor(as.character(obs$predictions_scanvi))

# Drop NA/empty labels if any
keep <- !is.na(obs$predictions_scanvi) & nzchar(as.character(obs$predictions_scanvi))
mat  <- mat[, keep]
obs  <- obs[keep, , drop = FALSE]

# --- Create CellChat object ---
cellchat <- createCellChat(object = mat, meta = obs, group.by = "predictions_scanvi")
cellchat@DB <- CellChatDB.human

# --- Run CellChat pipeline ---
future::plan("multisession", workers = 8)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Verify cell type labels
levels(cellchat@idents)
table(cellchat@idents)

# --- Visualize ---
groupSize <- as.numeric(table(cellchat@idents))

pdf("celltype_global_circlesmci_fe.pdf")
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

# --- Save results ---
saveRDS(cellchat, "cellchat_by_mci_fe.rds")
