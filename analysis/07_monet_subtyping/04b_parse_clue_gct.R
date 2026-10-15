#!/usr/bin/env Rscript
# ==============================================================================
# 04b_parse_clue_gct.R
# ==============================================================================
# Description:
#   Parses a raw CLUE/CMap query result (.gct) for one AD-subtype signature
#   (S1=AD1 ... S4=AD4, each queried against CLUE separately from the
#   subtype-vs-Control pseudobulk MAST DEG signature) and writes the filtered
#   top-compound table consumed by 05_cmap_heatmap.py.
#
#   Filter: norm_cs < -1.2 & fdr_q_nlog10 > 1.3 (i.e. FDR q < 0.05), matching
#   Methods, "Drug repurposing analysis". AD1 (S1) had no compounds pass this filter -- its run uses
#   THRESHOLD <- FALSE, which skips the filter and writes the unfiltered,
#   rank-ordered hits instead ("shown for comparison only" per Methods, "Drug repurposing analysis");
#   its output file is named accordingly (s1_nothreshold_...).
#
# Inputs:
#   A CLUE sig_queryl1k_tool result directory containing arfs/TAG/query_result.gct
#   (one such directory per signature S1-S4; not included in this repo -- CLUE
#   query results must be regenerated from the subtype DEG signatures via the
#   CLUE website/API, see https://clue.io)
#
# Outputs:
#   data/cmap/s{N}_cmap_top_compounds.csv (or s1_nothreshold_... for S1)
#
# Usage:
#   Run once per signature, editing SIGNATURE/GCT_DIR/THRESHOLD below (or
#   pass as Rscript args: SIGNATURE GCT_DIR THRESHOLD), e.g.:
#     Rscript 04b_parse_clue_gct.R S4 ./s4_my_analysis.sig_queryl1k_tool.<hash> TRUE
# ==============================================================================

library(cmapR)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
SIGNATURE <- if (length(args) >= 1) args[[1]] else "S4"
GCT_DIR   <- if (length(args) >= 2) args[[2]] else "./s4_my_analysis.sig_queryl1k_tool.<hash>"
THRESHOLD <- if (length(args) >= 3) as.logical(args[[3]]) else TRUE   # FALSE for S1/AD1

OUT_DIR <- file.path(dirname(dirname(getwd())), "data", "cmap")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sig_lower <- tolower(SIGNATURE)
out_file <- if (THRESHOLD) {
  file.path(OUT_DIR, sprintf("%s_cmap_top_compounds.csv", sig_lower))
} else {
  file.path(OUT_DIR, sprintf("%s_nothreshold_cmap_top_compounds.csv", sig_lower))
}

## ----- 1) Read the .gct using cmapR -----
gct_path <- file.path(GCT_DIR, "arfs", "TAG", "query_result.gct")
ds <- parse_gctx(gct_path)

## ----- 2) Build a tidy data frame (row descriptors + scores) -----
df <- cbind(ds@rdesc, as.data.frame(ds@mat, check.names = FALSE)) %>% as_tibble()
df <- as.data.frame(df)
df[df == "-666"] <- NA

df2 <- subset(df, qc_pass == 1)
df2 <- subset(df2, pert_type == "trt_cp")

if (THRESHOLD) {
  ## fdr_q_nlog10 > 1.3 == q < 0.05
  df3 <- subset(df2, norm_cs < -1.2 & fdr_q_nlog10 > 1.3)
} else {
  ## No compounds passed threshold for this signature (e.g. AD1) --
  ## write the unfiltered, rank-ordered hits "for comparison only" (Methods, "Drug repurposing analysis")
  df3 <- df2[order(df2$norm_cs, decreasing = FALSE), ]
}

hits_sets <- subset(df3, moa != "<NA>")
hits_sets <- hits_sets[order(hits_sets$norm_cs, decreasing = FALSE), ]

write.csv(hits_sets, file = out_file)
message(sprintf("Wrote %d compounds to %s", nrow(hits_sets), out_file))
