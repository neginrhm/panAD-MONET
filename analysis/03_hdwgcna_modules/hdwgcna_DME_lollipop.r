# ============================================================
# hdwgcna_DME_lollipop.r
# ============================================================
# Find differentially-expressed modules (DMEs) between diagnosis
# groups (Dementia vs Control, Dementia vs MCI, MCI vs Control)
# for Microglia, then produce lollipop plots of DME log2FC
# across all major cell types.
#
# Inputs
# ------
# - hdWGCNA Seurat object for Microglia (RDS)
# - Pre-computed DME CSVs per cell type (e.g. exc_DME_dementia_vs_control.csv)
#
# Outputs
# -------
# - mic_DME_dementia_vs_control.csv
# - mic_DME_dementia_vs_mci.csv
# - mic_DME_mci_vs_control.csv
# - Lollipop plot PDFs
#
# Usage
# -----
#   Rscript hdwgcna_DME_lollipop.r
# ============================================================

# --- Load libraries ----------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(forcats)

theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 8)
set.seed(12345)

# --- Configuration -----------------------------------------------------------

BASE_DIR <- "./data"

# --- Compute Microglia DMEs --------------------------------------------------

r <- readRDS(file.path(BASE_DIR, "CellTypes/Microglia/hdWGCNA/13ad_microglia_hdWGCNA_object.rds"))

group1 <- r@meta.data %>% subset(dx == "Control") %>% rownames
group2 <- r@meta.data %>% subset(dx == "Dementia") %>% rownames
group3 <- r@meta.data %>% subset(dx == "MCI") %>% rownames

# Dementia vs Control
DMEs <- FindDMEs(
  r,
  barcodes1 = group2,
  barcodes2 = group1,
  test.use = 'wilcox',
  wgcna_name = 'mic'
)
write.csv(DMEs, "mic_DME_dementia_vs_control.csv")
rm(DMEs)

# Dementia vs MCI
DMEs <- FindDMEs(
  r,
  barcodes1 = group2,
  barcodes2 = group3,
  test.use = 'wilcox',
  wgcna_name = 'mic'
)
write.csv(DMEs, "mic_DME_dementia_vs_mci.csv")
rm(DMEs)

# MCI vs Control
DMEs <- FindDMEs(
  r,
  barcodes1 = group3,
  barcodes2 = group1,
  test.use = 'wilcox',
  wgcna_name = 'mic'
)
write.csv(DMEs, "mic_DME_mci_vs_control.csv")
rm(r, DMEs)

# --- Lollipop plot: 5 cell types (basic) ------------------------------------

file_paths_5ct <- c(
  "exc_DME_dementia_vs_control.csv",
  "INH_DME_dementia_vs_control.csv",
  "odc_DME_dementia_vs_control.csv",
  "ast_DME_dementia_vs_control.csv",
  "mic_DME_dementia_vs_control.csv"
)

cell_types_5ct <- c("Excitatory", "Inhibitory", "Oligodendrocytes", "Astrocytes", "Microglia")

dme_data <- lapply(seq_along(file_paths_5ct), function(i) {
  data <- read.csv(file_paths_5ct[i])
  data$CellType <- cell_types_5ct[i]
  return(data)
}) %>%
  bind_rows()

dme_data$module <- factor(dme_data$module, levels = unique(dme_data$module))

pdf("lollipop Plot of avg_log2FC Across Modules and Cell Types 3.pdf")
ggplot(dme_data, aes(x = avg_log2FC, y = module, color = CellType)) +
  geom_segment(aes(x = 0, xend = avg_log2FC, y = module, yend = module), size = 0.5) +
  geom_point(size = 3) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 1) +
  labs(
    title = "Lollipop Plot of avg_log2FC Across Modules and Cell Types",
    x = "avg_log2FC",
    y = "Modules"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 8, margin = margin(r = 12)),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 11),
    panel.spacing.y = unit(1.5, "lines"),
    legend.position = "none"
  )
dev.off()

# --- Lollipop plot: 7 cell types (sorted, with gradient) --------------------

file_paths_7ct <- c(
  "exc_DME_dementia_vs_control.csv",
  "inh_DME_dementia_vs_control.csv",
  "odc_DME_dementia_vs_control.csv",
  "ast_DME_dementia_vs_control.csv",
  "mic_DME_dementia_vs_control.csv",
  "opc_DME_dementia_vs_control.csv",
  "endo_DME_dementia_vs_control.csv"
)

cell_types_7ct <- c("Excitatory", "Inhibitory", "Oligodendrocytes",
                    "Astrocytes", "Microglia", "OPC", "Endo")
cell_type_initials <- c("E", "I", "O", "A", "M", "OP", "EN")

dme_data <- lapply(seq_along(file_paths_7ct), function(i) {
  data <- read.csv(file_paths_7ct[i])
  data$CellType <- cell_types_7ct[i]
  data$Module_Updated <- paste0(data$module, "-", cell_type_initials[i])
  return(data)
}) %>%
  bind_rows()

# Reorder modules within each CellType by avg_log2FC
dme_data <- dme_data %>%
  group_by(CellType) %>%
  arrange(avg_log2FC, .by_group = TRUE) %>%
  mutate(Module_Updated = factor(Module_Updated,
                                 levels = Module_Updated[order(avg_log2FC)])) %>%
  ungroup()

# --- Save sorted lollipop plot -----------------------------------------------

pdf("Final_Sorted_Lollipop_Plot_with_Updated_Module_Names_2.pdf")
ggplot(dme_data, aes(x = avg_log2FC, y = Module_Updated,
                     color = avg_log2FC, size = abs(avg_log2FC))) +
  geom_segment(aes(x = 0, xend = avg_log2FC,
                   y = Module_Updated, yend = Module_Updated), size = 1.5) +
  geom_point() +
  geom_point(shape = 21, color = "grey30", stroke = 0.4, fill = NA) +
  expand_limits(x = c(min(dme_data$avg_log2FC) - 0.5,
                       max(dme_data$avg_log2FC) + 0.5)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(15, 15, 15, 15)) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 1) +
  labs(
    title = "Lollipop Plot of avg_log2FC Across Modules and Cell Types with Unique Module Names",
    x = "avg_Log2FC",
    y = "Modules"
  ) +
  scale_color_gradient(low = "#51c8c1", high = "#ff6538") +
  scale_size(range = c(2, 7)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 2, face = "bold"),
    axis.text.y = element_text(size = 8, margin = margin(r = 14)),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 9),
    panel.spacing.y = unit(0.2, "lines"),
    legend.position = "right"
  )
dev.off()
