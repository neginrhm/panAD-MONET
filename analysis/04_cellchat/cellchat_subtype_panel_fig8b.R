#!/usr/bin/env Rscript
# ==============================================================================
# cellchat_subtype_panel_fig8b.R
# ==============================================================================
# Combines the 5 per-group CellChat .rds objects (Control, AD1, AD2, AD3, AD4
# -- each built the same way as cellchat_function.r in ../04_cellchat, grouped
# by cell type) into 4 Control-vs-AD{n} merged objects and produces the
# unified rankNet info-flow panel with a consistent Control+AD1-4 color
# palette, matching Fig. 8b / Suppl. Fig. 6a (e.g. ANGPTL/RA/TGF-beta
# information-flow bars across Control and all four subtypes).
#
# Extracted and cleaned from the tail of an exploratory per-subtype build
# script; the object-construction boilerplate (identical to
# cellchat_function.r) and one-off sanity-check plots were removed.
#
# Inputs:
#   cellchat_by_control.rds, cellchat_by_ad1.rds, ..., cellchat_by_ad4.rds
#   (one CellChat object per group, built via cellchat_function.r)
#
# Outputs:
#   all_ranknets_subtypes.pdf -- 4-panel rankNet comparison (Control vs each
#   subtype), sharing one Control+AD1-4 palette across panels
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
})

BASE_DIR <- getwd()  # Set to project root (containing the 5 .rds objects) before running

ct  <- readRDS(file.path(BASE_DIR, "cellchat_by_control.rds"))
ad1 <- readRDS(file.path(BASE_DIR, "cellchat_by_ad1.rds"))
ad2 <- readRDS(file.path(BASE_DIR, "cellchat_by_ad2.rds"))
ad3 <- readRDS(file.path(BASE_DIR, "cellchat_by_ad3.rds"))
ad4 <- readRDS(file.path(BASE_DIR, "cellchat_by_ad4.rds"))

object.list1 <- list(ct = ct, ad1 = ad1)
cellchat1 <- mergeCellChat(object.list1, add.names = names(object.list1))

object.list2 <- list(ct = ct, ad2 = ad2)
cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2))

object.list3 <- list(ct = ct, ad3 = ad3)
cellchat3 <- mergeCellChat(object.list3, add.names = names(object.list3))

object.list4 <- list(ct = ct, ad4 = ad4)
cellchat4 <- mergeCellChat(object.list4, add.names = names(object.list4))

# Consistent Control+AD1-4 palette across all 4 panels
PALETTE <- c(
  ct  = "#20B2AA",  # control (teal)
  ad1 = "#6C5CE7",  # indigo
  ad2 = "#E69F00",  # orange
  ad3 = "#D55E00",  # vermilion
  ad4 = "#CC79A7"   # pink
)

pdf(file.path(BASE_DIR, "all_ranknets_subtypes.pdf"))
gg1 <- rankNet(cellchat1, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, color.use = PALETTE)
gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, color.use = PALETTE)
gg3 <- rankNet(cellchat3, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, color.use = PALETTE)
gg4 <- rankNet(cellchat4, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, color.use = PALETTE)
gg1 + gg2 + gg3 + gg4
dev.off()
