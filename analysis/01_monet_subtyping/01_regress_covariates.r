###############################################################################
# Covariate Regression — all 7 cell types
#
# Generalised from regree_out.r (scripts/regree_out.r).
# Runs bootstrap linear regression per gene (R=1000) to remove Study, Age,
# and Sex effects. Coefficients are summarised with median across bootstraps.
#
# Input:
#   pseudobulk_{CellType}.csv          (output of 00_pseudobulk.r)
#   dementia_with_metadata_final.csv   (columns: Sample, dx, Age, Sex, Study)
#
# Output:
#   Regressed_data_DEM_{CT}_thisexp~Study+Age+Sex.csv  (one per cell type)
#   (these filenames are parsed by 02_build_feature_vectors.py)
#
# Runtime: ~30–60 min per cell type on HPC (5000 genes x 1000 bootstraps)
###############################################################################

library(boot)
library(dplyr)

options(stringsAsFactors = FALSE)

numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices, ]
  return(coef(lm(formula, data = d)))
}

# --- metadata ---
targets       <- read.csv("../data/dementia_with_metadata_final.csv")
targets$Sample <- targets$Unnamed..0
targets_unique <- targets[!duplicated(targets$Sample), ]

# --- cell type → input pseudobulk CSV mapping ---
cell_types <- list(
  list(input = "pseudobulk_Astrocyte.csv",      abbrev = "AST"),
  list(input = "pseudobulk_Microglia.csv",       abbrev = "MIC"),
  list(input = "pseudobulk_ENDO.csv",            abbrev = "ENDO"),
  list(input = "pseudobulk_Oligodendrocyte.csv", abbrev = "ODC"),
  list(input = "pseudobulk_OPC.csv",             abbrev = "OPC"),
  list(input = "pseudobulk_Inhibitory.csv",      abbrev = "INH"),
  list(input = "pseudobulk_Excitatory.csv",      abbrev = "EXC")
)

for (ct in cell_types) {
  cat("\n=== Processing", ct$abbrev, "===\n")

  datExpr <- read.csv(ct$input, row.names = 1)
  colnames(datExpr) <- gsub("\\.", "-", colnames(datExpr))

  common_samples <- intersect(colnames(datExpr), targets_unique$Sample)
  datExpr        <- datExpr[, common_samples]
  targets_bulk   <- targets_unique[match(common_samples, targets_unique$Sample), ]
  stopifnot(all(colnames(datExpr) == targets_bulk$Sample))

  datExpr  <- datExpr[, sapply(datExpr, is.numeric)]
  normExpr <- log2(datExpr + 1)

  Study   <- as.numeric(factor(targets_bulk$Study))
  Age     <- as.numeric(targets_bulk$Age)
  Sex     <- as.numeric(factor(targets_bulk$Sex))
  regvars <- data.frame(Study, Age, Sex)

  normExpr.reg <- matrix(NA, nrow = nrow(normExpr), ncol = ncol(normExpr))
  rownames(normExpr.reg) <- rownames(normExpr)
  colnames(normExpr.reg) <- colnames(normExpr)

  set.seed(8675309)
  for (i in seq_len(nrow(normExpr))) {
    if (i %% 500 == 0) cat("  gene", i, "of", nrow(normExpr), "\n")
    thisexp    <- as.numeric(normExpr[i, ])
    bs.results <- boot(
      data      = data.frame(thisexp, regvars),
      statistic = bs,
      R         = numboot,
      formula   = thisexp ~ Study + Age + Sex
    )
    # median across bootstrap resamples (robust coefficient estimate)
    bs.stats <- sapply(seq_len(ncol(bs.results$t)), function(n)
      median(na.omit(bs.results$t[, n]))
    )
    normExpr.reg[i, ] <- thisexp -
      bs.stats[2] * regvars$Study -
      bs.stats[3] * regvars$Age -
      bs.stats[4] * regvars$Sex
  }

  out_file <- paste0("Regressed_data_DEM_", ct$abbrev, "_thisexp~Study+Age+Sex.csv")
  write.csv(normExpr.reg, out_file)
  cat("Saved:", out_file, "\n")
}
