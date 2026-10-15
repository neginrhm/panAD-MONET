#!/usr/bin/env Rscript
# ============================================================================
# 00_hdwgcna_construct_network.r
# ============================================================================
# Builds the per-cell-type hdWGCNA co-expression network (Figs. 2-3,
# Suppl. Tables 2-3): Seurat object from the cell-type h5ad -> Harmony ->
# metacells -> soft-power selection -> network/modules -> module eigengenes,
# kME, hub genes -> GO enrichment (enrichR).
#
# Consolidated from the per-cell-type working scripts used for the
# manuscript (microglia, oligodendrocytes, astrocytes), which share one
# template and differ only in the network/enrichment settings in
# CT_SETTINGS below. The excitatory/inhibitory/OPC/endothelial networks were
# built with the same template, but their original settings were not
# preserved; for those cell types this script uses the hdWGCNA defaults
# (automatic soft power, deepSplit=4, minModuleSize=50, mergeCutHeight=0.2).
#
# Steps (as in the original runs):
#   1. Read 13ad_{CellType}.h5ad (counts in X), drop "Brain Cancer" donors,
#      map Diagnosis -> dx (Control / MCI / Dementia)
#   2. LogNormalize (scale factor 1e4), ScaleData, PCA, Harmony on Study,
#      neighbors/UMAP on the Harmony embedding (the working scripts also ran
#      an earlier Harmony on Sample, which the Study run overwrites)
#   3. SetupForWGCNA (genes expressed in >= 5% of cells)
#   4. MetacellsByGroups (Celltype_markers x Sample, k = 25, max_shared = 10,
#      KNN on Harmony), NormalizeMetacells, SetDatExpr
#   5. TestSoftPowers (signed network), ConstructNetwork
#   6. ModuleEigengenes (harmonized by Sample), ModuleConnectivity (kME),
#      GetHubGenes (top 10 per module)
#   7. RunEnrichr on module genes (GO BP/CC/MF)
#
# Usage:
#   Rscript 00_hdwgcna_construct_network.r <CellType>
#   CellType: Microglia | Oligodendrocyte | Astrocyte | Excitatory |
#             Inhibitory | OPC | ENDO
#
# Inputs (working directory):
#   13ad_{CellType}.h5ad
#
# Outputs (working directory, prefix 13ad_{tag}_):
#   power_table.csv, hMEs.csv, MEs.csv, modules_table.csv,
#   modules_hub_df_table.csv, hdwgcna_enrichr.csv, dendrogram.pdf,
#   soft_powers.pdf, hdWGCNA_object.rds
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(harmony)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(WGCNA)
  library(hdWGCNA)
  library(enrichR)
})

theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 8)
set.seed(12345)

# --- Per-cell-type settings ---------------------------------------------------
# tag: file prefix / wgcna_name; group_name: value in Celltype_markers.
# soft_power NULL = hdWGCNA automatic selection from TestSoftPowers.
# network_args NULL = hdWGCNA ConstructNetwork defaults.
CT_SETTINGS <- list(
  Microglia = list(tag = "microglia", wgcna_name = "mic", group_name = "Microglia",
                   soft_power = NULL, network_args = NULL,
                   enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                                   "GO_Molecular_Function_2023"),
                   enrichr_max_genes = Inf),
  Oligodendrocyte = list(tag = "odc", wgcna_name = "odc", group_name = "odc",
                         soft_power = 14,
                         network_args = list(deepSplit = 4, minModuleSize = 25, mergeCutHeight = 0.25),
                         enrichr_dbs = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                         "GO_Molecular_Function_2021"),
                         enrichr_max_genes = 100),
  Astrocyte = list(tag = "astrocye", wgcna_name = "ast", group_name = "Astrocyte",
                   soft_power = NULL,
                   network_args = list(deepSplit = 4, minModuleSize = 25, mergeCutHeight = 0.25),
                   enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                                   "GO_Molecular_Function_2023"),
                   enrichr_max_genes = 100),
  Excitatory = list(tag = "excitatory", wgcna_name = "exc", group_name = "Excitatory",
                    soft_power = NULL, network_args = NULL,
                    enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                                    "GO_Molecular_Function_2023"),
                    enrichr_max_genes = 100),
  Inhibitory = list(tag = "INH", wgcna_name = "inh", group_name = "Inhibitory",
                    soft_power = NULL, network_args = NULL,
                    enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                                    "GO_Molecular_Function_2023"),
                    enrichr_max_genes = 100),
  OPC = list(tag = "opc", wgcna_name = "opc", group_name = "OPC",
             soft_power = NULL, network_args = NULL,
             enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                             "GO_Molecular_Function_2023"),
             enrichr_max_genes = 100),
  ENDO = list(tag = "endo", wgcna_name = "endo", group_name = "ENDO",
              soft_power = NULL, network_args = NULL,
              enrichr_dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                              "GO_Molecular_Function_2023"),
              enrichr_max_genes = 100)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !(args[1] %in% names(CT_SETTINGS))) {
  stop("Usage: Rscript 00_hdwgcna_construct_network.r <",
       paste(names(CT_SETTINGS), collapse = "|"), ">")
}
cell_type <- args[1]
cfg <- CT_SETTINGS[[cell_type]]
out <- function(x) paste0("13ad_", cfg$tag, "_", x)

# --- 1. Load and harmonize diagnosis --------------------------------------------
data <- read_h5ad(paste0("13ad_", cell_type, ".h5ad"))
seurat_obj <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
rm(data); gc()

seurat_obj <- subset(seurat_obj, subset = Diagnosis != "Brain Cancer")
seurat_obj$dx <- ifelse(
  seurat_obj$Diagnosis %in% c("ADDS", "AD", "AD, Vascular Dementia", "Dementia",
                              "DSAD", "PD", "Vascular Dementia"), "Dementia",
  ifelse(seurat_obj$Diagnosis %in% c("No Dementia", "NCI", "Control"), "Control",
         ifelse(seurat_obj$Diagnosis == "MCI", "MCI", NA)))
print(table(seurat_obj$dx, useNA = "ifany"))

# --- 2. Normalize, PCA, Harmony, UMAP --------------------------------------------
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)
seurat_obj <- RunPCA(seurat_obj, features = all_genes)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Study", plot_convergence = TRUE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony")

# --- 3-4. hdWGCNA setup and metacells ------------------------------------------------
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = "fraction", fraction = 0.05,
                            wgcna_name = cfg$wgcna_name)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Celltype_markers", "Sample"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "Celltype_markers"
)
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(seurat_obj, group_name = cfg$group_name,
                         group.by = "Celltype_markers", assay = "RNA", slot = "data")

# --- 5. Soft power and network -----------------------------------------------------
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
pdf(out("soft_powers.pdf"))
print(wrap_plots(PlotSoftPowers(seurat_obj), ncol = 2))
dev.off()
write_csv(GetPowerTable(seurat_obj), out("power_table.csv"))

net_args <- c(list(seurat_obj, tom_name = cfg$wgcna_name), cfg$network_args)
if (!is.null(cfg$soft_power)) net_args$soft_power <- cfg$soft_power
seurat_obj <- do.call(ConstructNetwork, net_args)

pdf(out("dendrogram.pdf"))
PlotDendrogram(seurat_obj, main = paste(cell_type, "hdWGCNA Dendrogram"))
dev.off()

# --- 6. Eigengenes, connectivity, hub genes ------------------------------------------
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "Sample")
write.csv(GetMEs(seurat_obj), out("hMEs.csv"))
write.csv(GetMEs(seurat_obj, harmonized = FALSE), out("MEs.csv"))

seurat_obj <- ModuleConnectivity(seurat_obj, group.by = "Celltype_markers",
                                 group_name = cfg$group_name)
write.csv(GetModules(seurat_obj), out("modules_table.csv"))
write.csv(GetHubGenes(seurat_obj, n_hubs = 10), out("modules_hub_df_table.csv"))

# --- 7. GO enrichment of module genes -------------------------------------------------
seurat_obj <- RunEnrichr(seurat_obj, dbs = cfg$enrichr_dbs, max_genes = cfg$enrichr_max_genes)
write_csv(GetEnrichrTable(seurat_obj), out("hdwgcna_enrichr.csv"))

saveRDS(seurat_obj, file = out("hdWGCNA_object.rds"))
message("Done: ", out("hdWGCNA_object.rds"))
