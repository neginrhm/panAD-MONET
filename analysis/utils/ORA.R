###############################################################################
# Over-Representation Analysis (ORA) using Fisher's Exact Test
#
# Description:
#   Provides two functions for over-representation analysis:
#   - OR(): Runs Fisher's exact test on a 2x2 contingency table and returns
#           odds ratio, p-value, and 95% confidence interval.
#   - ORA(): Computes overlap between a test gene set and a reference gene set,
#            then calls OR() to assess enrichment significance.
#
# Author: Vivek Swarup, PhD, UC Irvine (vswarup@uci.edu)
#
# Inputs (ORA function):
#   - testpath:       character vector of genes in the test set (e.g., DEGs)
#   - refpath:        character vector of genes in the reference set (e.g., a module)
#   - testbackground: character vector of all genes in the test background
#   - refbackground:  character vector of all genes in the reference background
#
# Outputs:
#   Named character vector with: OR, Fisher p, -95%CI, +95%CI,
#   Overlap, Reference List, Input List, Background, % List Overlap
#
# Usage:
#   source("ORA.R")
#   result <- ORA(my_degs, module_genes, all_test_genes, all_ref_genes)
###############################################################################

OR <- function(q, k, m, t) {
  # q: intersection of test list and reference list (white balls drawn)
  # m: all genes in reference list (number of draws)
  # k: all genes in test list (white balls)
  # t: total genes assessed (black + white balls)

  fisher.out <- fisher.test(
    matrix(c(q, k - q, m - q, t - m - k + q), 2, 2),
    conf.int = TRUE
  )

  output <- c(fisher.out$estimate, fisher.out$p.value,
              fisher.out$conf.int[1], fisher.out$conf.int[2])
  names(output) <- c("OR", "Fisher p", "-95%CI", "+95%CI")
  return(output)
}

ORA <- function(testpath, refpath, testbackground, refbackground) {
  q <- length(intersect(testpath, refpath))
  k <- length(intersect(refpath, testbackground))
  m <- length(intersect(testpath, refbackground))
  t <- length(intersect(testbackground, refbackground))

  empvals <- OR(q, k, m, t)

  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals, q, k, m, t, 100 * signif(q / k, 3)))
  names(empvals) <- c(tmpnames, "Overlap", "Reference List", "Input List",
                      "Background", "% List Overlap")
  return(empvals)
}
