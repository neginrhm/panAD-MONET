
library(Seurat)
library(hdf5r)
library(SeuratData)
library(SeuratDisk)
library(anndata)
library(harmony)
library(reticulate)
library(dplyr)
library(ggplot2)

library(tidyverse)
library(cowplot)
library(patchwork)

library(WGCNA)
library(hdWGCNA)
library(UCell)

library(enrichR)


theme_set(theme_cowplot())

enableWGCNAThreads(nThreads = 8)
set.seed(12345)


# Update this path to your own conda environment before running.
reticulate::use_condaenv(Sys.getenv("PYR_CONDAENV", "pyr-env"), required = TRUE)


data <- read_h5ad("12ad_Microglia.h5ad")

data <- read_h5ad("12ad_3m_5x_sparse.h5ad")

'''
> data
AnnData object with n_obs × n_vars = 164638 × 5000
    obs: 'Sample', 'Study', 'org_barcode', 'n_genes_by_counts', 'total_counts', 'total_counts_MT', 'pct_counts_MT', 'n_genes', 'n_counts', 'barcode2', 'doublet', 'doublet_score', 'external_dono_donor_name_label', '_scvi_batch', 'labels_scanvi', '_scvi_labels', 'predictions_scanvi', 'UMAP_1', 'UMAP_2', 'leiden', 'Doublet', 'Samples', 'barcode_index', 'Age', 'Sex', 'Diagnosis', 'Tsu'Dissue', 'Diagnosis2', 'Primary_Function', 'Celltype_markers', 'label_transfer_annot'
    var: 'gene_ids-0', 'gene_ids-1', 'gene_names-1', 'MT', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'
    uns: 'Celltype_markers_colors', 'Diagnosis2_colors', 'Primary_Function_colors', 'Sex_colors', 'Study_colors', '_scvi_manager_uuid', '_scvi_uuid', 'dendrogram_leiden', 'doublet_colors','lab 'labels_scanvi_colors', 'leiden', 'leiden_colors', 'leiden_sizes', 'neighbors', 'paga', 'predictions_scanvi_colors', 'umap'
    obsm: 'X_pca', 'X_scANVI', 'X_umap'
    layers: 'raw_count'
    obsp: 'connectivities', 'distances'
'''

data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)


# Assuming 'data' is your Seurat object and it has a metadata column named 'Diagnosis'
library(Seurat)
# Subset the Seurat object to exclude rows where 'Diagnosis' is "Brain Cancer"
data <- subset(data, subset = Diagnosis != "Brain Cancer")
# Optionally, print a summary to confirm the rows have been removed
table(data$Diagnosis)

## adjust the dx
data@meta.data$dx <- ifelse(data@meta.data$Diagnosis %in% c("ADDS", "AD", "AD, Vascular Dementia", "Dementia", "DSAD", "PD", "Vascular Dementia"), "Dementia",
                            ifelse(data@meta.data$Diagnosis %in% c("No Dementia", "NCI","Control"), "Control",
                                   ifelse(data@meta.data$Diagnosis == "MCI", "MCI", NA)))  # Ensure MCI and other cases are handled

# Verify by looking at the updated unique values
unique(data@meta.data$dx)


umap_data <- FetchData(data, c("UMAP_1", "UMAP_2", "Celltype_markers","dx","Study"))

png("12ad_microglia_scatter_study.png")

# Now plot using this extracted data
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Study)) +
  geom_point(alpha = 0.5, size = 1.5) +
  labs(title = "Scatter Plot of UMAP Coordinates", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_viridis_d()  # Using a color-blind friendly palette
dev.off()


seurat_obj<- data



###################################################
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Extracting cluster and study data from the Seurat object
data_frame <- data.frame(
  Cluster = Idents(seurat_obj),  # Extracts the cluster identity
  Study = seurat_obj$Study  # Assumes 'Study' is a metadata column in your Seurat object
)

# Create a summary table with counts and then calculate percentages
cluster_study_counts <- data_frame %>%
  group_by(Cluster, Study) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Cluster) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)  # Calculate percentage contribution

# Define a custom color palette with enough colors for each study
# You can add more colors to this list if you have more than 10 studies
colors <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
            "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
            "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
            "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
            "#ffffff", "#000000")

# Create the stacked percentage bar plot with custom colors
p <- ggplot(cluster_study_counts, aes(x = as.factor(Cluster), y = Percentage, fill = Study)) +
  geom_bar(stat = "identity", position = "fill") +  # Uses 'fill' to stack to 100%
  scale_fill_manual(values = colors) +  # Apply custom colors
  labs(title = "Percentage Distribution of Cells by Study and Cluster",
       x = "Cluster",
       y = "Percentage (%)",
       fill = "Study") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improves x-axis label readability

# Save the plot to a PDF file
# ggsave("microglia_Cluster_Study_Percentage_Distribution_Custom_Colors_harmony3.pdf", plot = p, width = 12, height = 8)
ggsave("Astrocyte_Cluster_Study_Percentage_Distribution_Custom_Colors_after_harmony3_final.pdf", plot = p, width = 12, height = 8)

# Print the plot to the console for immediate viewing
print(p)
#############################################################


seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = all.genes)

# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# seurat_obj <- RunHarmony(seurat_obj, "Sample")
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = 'Study',  plot_convergence = TRUE)
print("harmony is done")
#Warning: The following 1 features requested have zero variance; running reduction without them: OR7A5


seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction='harmony')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6, reduction='harmony')
# seurat_obj <- RunHarmony(seurat_obj, group.by.vars = 'Study',  plot_convergence = TRUE)
print("clustering is done")


seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = 'harmony')
print("RunUMAP is done")


png("12ad_microglia_normalized_scaled_umap.png")
DimPlot(seurat_obj, reduction = "umap")
dev.off()


png("12ad_microglia_normalized_scaled_umap_study_harmony3.png")
DimPlot(seurat_obj, reduction = "harmony", group.by = "Study")
dev.off()

saveRDS(seurat_obj, '12ad_microglia_post_umap2.rds')




##### performing hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "mic" # the name of the hdWGCNA experiment
)



seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Celltype_markers", "Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Celltype_markers' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = 'Microglia', # the name of the group of interest in the group.by column
  group.by='Celltype_markers', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

pdf("soft_powers_microgliac2.pdf")
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)

dev.off()


power_table <- GetPowerTable(seurat_obj)
head(power_table)
write_csv(power_table, "12ad_mic_power_table2.csv")

seurat_obj <- ConstructNetwork(seurat_obj)

# seurat_obj <- ConstructNetwork(seurat_obj, tom_name="MIC",soft_power=14)
# seurat_obj <- ConstructNetwork(seurat_obj, tom_name="MIC",soft_power=14, deepSplit =4, minModuleSize =25,mergeCutHeight = 0.25, overwrite_tom = TRUE)


pdf("12ad_microglia_wgcna_dendogram_softpower14_ds4_mm25_merge25_test.pdf")
PlotDendrogram(seurat_obj, main='MIC hdWGCNA Dendrogram')
dev.off()

TOM <- GetTOM(seurat_obj)


# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="Sample"
)



# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
write.csv(hMEs, "12ad_microglia_hMEs.csv")
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.csv(MEs, "12ad_microglia_MEs.csv")



# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Celltype_markers'
)


UpdateSeuratObject(object)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "MIC"
)



seurat_obj <- ModuleConnectivity(seurat_obj, group.by = 'Celltype_markers')
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')


pdf("12ad_microglia_plotKMEs_test2.pdf")
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=3)

p
dev.off()



# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])
write.csv(modules, "12ad_microglia_modules_table.csv")

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)
write.csv(hub_df, "12ad_microglia_modules_hub_df_table.csv")


saveRDS(seurat_obj, file='12ad_microglia_hdWGCNA_object.rds')


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)


png("12ad_microglia_ModuleFeaturePlot_hMEs.png")
# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
dev.off()



png("12ad_microglia_ModuleFeaturePlot_signatueScore.png")
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
dev.off()



png("12ad_microglia_ModuleCorrelogram.png")
# plot module correlagram
ModuleCorrelogram(seurat_obj)
dev.off()






# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)



# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'Celltype_markers')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

png("12ad_microglia_plot_dotplot_module_scores.png")

# plot output
p
dev.off()




png("12ad_microglia_modulenetworkplot.png")
ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks'
)
dev.off()



# Load necessary library
library(ggplot2)

# Assuming seurat_obj_mic$dx and seurat_obj_mic$blue are already loaded in your environment

# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$blue,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_blue_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)


#############################

# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$turquoise,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_turquoise_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)

#############################


# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$brown,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_brown_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)

#############################

# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$yellow,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_yellow_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)

#############################

# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$green,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_green_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)

#############################

# Create a data frame
data <- data.frame(
  Expression = seurat_obj_mic$grey,
  Diagnosis = seurat_obj_mic$dx
)

# Create the boxplot with explicit theme adjustments for a white background
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),  # Set the panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Set the plot background to white
  )

# Print the plot
print(plot)

# Save the plot
ggsave("mic_grey_module_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)

#############################


# Load packages
library(ggplot2)
library(dplyr)
library(ggsignif)

data <- data.frame(
  Expression = seurat_obj_mic$grey,
  Diagnosis = seurat_obj_mic$dx
)

# Assuming 'data' is your DataFrame
# Create the boxplot
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white", colour = "white"))

# Adding statistical annotations
# Adjust the comparisons based on your pairwise test results
plot <- plot + geom_signif(
  comparisons = list(c("Control", "Dementia"), c("Control", "MCI")),
  map_signif_level = TRUE,  # Automatically map significance levels to annotations
  test = "wilcox.test",     # Specify the test used for comparisons
  step_increase = 0.1       # Adjust the vertical spacing of the annotations
)

# Print the plot
print(plot)

# Save the plot
ggsave("annotated_expression_boxplot.png", plot, width = 10, height = 8, units = "in", dpi = 300)
#############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)


data <- data.frame(
  Expression = seurat_obj_mic$yellow,
  Diagnosis = seurat_obj_mic$dx
)
# Assuming 'data' is your DataFrame
# Subset data for each test
data_control_dementia <- subset(data, Diagnosis %in% c("Control", "Dementia"))
data_control_mci <- subset(data, Diagnosis %in% c("Control", "MCI"))

# Perform the statistical tests manually
test_control_dementia <- wilcox.test(Expression ~ Diagnosis, data = data_control_dementia)
test_control_mci <- wilcox.test(Expression ~ Diagnosis, data = data_control_mci)

# Function to format p-values for display
format_p_value <- function(p_value) {
  # Format p-values based on their magnitude
  if (p_value < 1e-10) {  # Use a stricter threshold for the smallest p-values
    "p < 1e-10"
  } else if (p_value < 0.001) {
    sprintf("p = %.3e", p_value)  # Scientific notation for small p-values
  } else {
    sprintf("p = %.3f", p_value)  # Regular formatting for larger p-values
  }
}

# Annotations for the plot based on test results
annotations <- c(
  format_p_value(test_control_dementia$p.value),
  format_p_value(test_control_mci$p.value)
)

# Create the boxplot
plot <- ggplot(data, aes(x = Diagnosis, y = Expression, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Expression Levels by Diagnostic Category", x = "Diagnosis", y = "Expression Level") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white", colour = "white"))

# Adding statistical annotations with formatted p-values
plot <- plot + geom_signif(
  comparisons = list(c("Control", "Dementia"), c("Control", "MCI")),
  annotations = annotations,
  tip_length = 0.05,      # Adjust the tip length of the signif bars
  step_increase = 0.1     # Adjust the vertical spacing of the annotations
)

# Print the plot
print(plot)

# Save the plot
ggsave("microglia_yellow_module_wilcox.png", plot, width = 10, height = 8, units = "in", dpi = 300)





#### run MAST
#In all these cases, a positive fold change indicates higher expression in the first group (ident.1) compared to the second group (ident.2),
#while a negative fold change indicates higher expression in the second group (ident.2) compared to the first group (ident.1).
#


# Set the identities to the diagnosis group
Idents(mic) <- mic$dx

# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "Control" group compared to the "Dementia" group.

# Control vs Dementia
control_vs_dementia <- FindMarkers(mic, ident.1 = "Dementia", ident.2 = "Control", test.use ="MAST", min.pct = 0,logfc.threshold = 0)


# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "Control" group compared to the "MCI" group.
# Control vs MCI
control_vs_mci <- FindMarkers(mic, ident.1 = "Control", ident.2 = "MCI", test.use ="MAST", min.pct = 0,logfc.threshold = 0)


# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "MCI" group compared to the "Dementia" group.
# MCI vs Dementia
mci_vs_dementia <- FindMarkers(mic, ident.1 = "MCI", ident.2 = "Dementia", test.use ="MAST", min.pct = 0,logfc.threshold = 0)


# Print the top results for each comparison
print(head(control_vs_dementia))
print(head(control_vs_mci))
print(head(mci_vs_dementia))

# Optionally, save the results to CSV files
write.csv(control_vs_dementia, file = "mast_control_vs_dementia_markers.csv")
write.csv(control_vs_mci, file = "mast_control_vs_mci_markers.csv")
write.csv(mci_vs_dementia, file = "mast_mci_vs_dementia_markers.csv")



----------------------------------------------------------------------------
# Filter results for p_val less than 0.05 and sort based on avg_log2FC
control_vs_dementia2 <- control_vs_dementia %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
control_vs_mci2 <- control_vs_mci %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
mci_vs_dementia2 <- mci_vs_dementia %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))

# Print the top results for each comparison
print(head(control_vs_dementia2))
print(head(control_vs_mci2))
print(head(mci_vs_dementia2))

# Optionally, save the results to CSV files
write.csv(control_vs_dementia2, file = "control_vs_dementia_markers_p0.05.csv")
write.csv(control_vs_mci2, file = "control_vs_mci_markers_p0.05.csv")
write.csv(mci_vs_dementia2, file = "mci_vs_dementia_markers_p0.05.csv")



----------------------------------------------------------------------------


# Define the patterns to remove
patterns_to_remove <- c("^AC0", "^MT-", "^AC1", "^AC026316")

# Combine the patterns into a single regular expression
combined_pattern <- paste(patterns_to_remove, collapse = "|")

# Remove rows where the gene names match the patterns
filtered_df <- control_vs_dementia2[!grepl(combined_pattern, rownames(control_vs_dementia)), ]

# Display the filtered dataframe
head(filtered_df)


# colnames(control_vs_dementia)[0] <- "gene"
control_vs_dementia2$gene <- rownames(control_vs_dementia2)
pdf("12ad_miccroglia_control_vs_dementia_volcanoplot2.pdf")
EnhancedVolcano(control_vs_dementia, lab = control_vs_dementia$gene,x = 'avg_log2FC',y = 'p_val_adj', title="control vs dementia")
dev.off()

------------------------------------------------------------------------------

control_vs_mci2$gene <- rownames(control_vs_mci2)
pdf("12ad_miccroglia_control_vs_mci2_volcanoplot_fc0.5.pdf")
EnhancedVolcano(control_vs_mci2, lab = control_vs_mci2$gene,x = 'avg_log2FC',y = 'p_val_adj', title="control vs MCI", subtitle = "", FCcutoff = 0.5, labSize = 4.0)
dev.off()

------------------------------------------------

mci_vs_dementia2$gene <- rownames(mci_vs_dementia2)
pdf("12ad_microglia_mci_vs_dementia2_volcanoplot_fc0.5.pdf")
EnhancedVolcano(mci_vs_dementia2, lab = mci_vs_dementia2$gene,x = 'avg_log2FC',y = 'p_val_adj', title="MCI vs Dementia", subtitle = "", FCcutoff = 0.5, labSize = 4.0)
dev.off()

------------------------------------------------
# performing enrichR

#downreg in control, up reg in dementia
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
enriched <- enrichr(control_vs_dementia2[control_vs_dementia2$avg_log2FC < 0, "gene"], dbs)

enriched[[3]] <- enriched[[3]][order(enriched[[3]]$P.value), ]
# print(head(enriched[[2]], 3))
# write.csv(enriched[[1]],"niched1.csv")
# Extract the first 10 terms and their combined scores

pdf("down_ENRICH_control_vs_dementia2.pdf", height =10, width=20)  # Adjust width and height as needed

if (websiteLive) {
    plotEnrich(enriched[[3]], showTerms = 20, numChar = 70, y = "Count", orderBy = "Combined.Score", title="Control vs Dementia - Downreg")
}
# Save the plot
dev.off()  # Closes the PNG



####### up reg in dementia

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

if (websiteLive) {
    enriched <- enrichr(control_vs_dementia2[control_vs_dementia2$avg_log2FC < 0, "gene"], dbs)
}

# if (websiteLive) enriched[["GO_Biological_Process_2021"]]

pdf("control_vs_dementia2_enrich_down_cellularprocess.pdf", height = 8, width = 11)  # Further increase width and height

# Adjust graphical parameters to increase font sizes
#par(mar = c(30, 50, 12, 15) + 0.1, cex.axis = 3, cex.lab = 2.5, cex.main = 3)  # Further increase text size

if (websiteLive) {
    plotEnrich(enriched[[3]], showTerms = 20, numChar = 70, y = "Count", orderBy = "Combined.Score", title="Control vs Dementia - Downreg")
}
dev.off()


'''
if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018",
           "GO_Biological_Process_2018")
  enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)
  # Plot top 20 terms from "GO_Biological_Process_2018" and ordered by P-value
  if (enrichRLive) {
    plotEnrich(enriched[[3]], showTerms = 20, numChar = 50, y = "Count",
               orderBy = "P.value")
  }
}

'''

-----------------------------------------------------------------------
#DME analysis comparing two groups

r <- seurat_obj
#group1 <- r@meta.data %>% subset(cell_type == 'INH' & dx == "Control") %>% rownames
group1 <- r@meta.data %>% subset(dx == "Control") %>% rownames
group2 <- r@meta.data %>% subset(dx == "Dementia") %>% rownames
group3 <- r@meta.data %>% subset(dx == "MCI") %>% rownames

head(group1)

DMEs <- FindDMEs(
  r,
  barcodes1 = group2,
  barcodes2 = group3,
  test.use='wilcox',
  wgcna_name='exc'
)

head(DMEs)

pdf("exc2_wgcna_Dementia_vs_control.pdf")
PlotDMEsVolcano(
  r,
  DMEs,
  wgcna_name = 'exc'
  #xlim_range = c(-2, 2)
)
dev.off()




pdf("exc2_PlotDMEsLollipop_exc2_wgcna_Dementia_vs_control.pdf")
PlotDMEsLollipop(
  r,
  DMEs,
  wgcna_name='exc',
  pvalue = "p_val_adj"
)
dev.off()


##hdwgcna enrichment

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = Inf # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

write_csv(enrich_df, "odcc_hdwgcna_enrichr.csv")

# make GO term plots:
EnrichrBarPlot(
seurat_obj,
outdir = "odc_enrichr_plots2", # name of output directory
n_terms = 12, # number of enriched terms to show (sometimes more show if there are ties!!!)
plot_size = c(5,7), # width, height of the output .pdfs
logscale=TRUE # do you want to show the enrichment as a log scale?
)


pdf("odc_dotplot_hdwgcna_goterm_modules.pdf")
# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=4, # number of terms for each module
  term_size=8, #
  p_adj = FALSE
) + scale_color_stepsn(colors=rev(viridis::magma(256)))
dev.off()



seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

umap_df <- GetModuleUMAP(seurat_obj)


pdf("ast_hdwgcna_ModuleUMAPPlot_umap.pdf")
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE

)
dev.off()


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Assuming r is your dataframe containing the data
# Create the contingency table
contingency_table <- table(r$Study, r$seurat_clusters)

# Convert the table to a data frame for easier plotting
df <- as.data.frame(contingency_table)
colnames(df) <- c("Study", "Cluster", "Count")

# Plot the stacked bar chart and save it as a PDF
pdf("Cells_Per_Study_Per_Cluster_jun12.pdf", height = 10, width = 15)

ggplot(df, aes(x = Cluster, y = Count, fill = Study)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells from Each Study in Each Cluster",
       x = "Cluster",
       y = "Number of Cells",
       fill = "Study") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()  # Close the PDF device


library(ggplot2)
library(dplyr)

# Assuming r is your dataframe containing the data
# Create the contingency table
contingency_table <- table(r$Study, r$seurat_clusters)

# Convert the table to a data frame for easier plotting
df <- as.data.frame(contingency_table)
colnames(df) <- c("Study", "Cluster", "Count")

# Calculate the percentage of cells from each study within each cluster
df <- df %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Define a custom color palette with at least 13 colors
custom_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31"
)

# Plot the stacked bar chart showing percentages and save it as a PDF
pdf("Cells_Per_Study_Per_Cluster_Percentage3.pdf", height = 10, width = 15)

ggplot(df, aes(x = Cluster, y = Percentage, fill = Study)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Percentage of Cells from Each Study in Each Cluster",
       x = "Cluster",
       y = "Percentage of Cells",
       fill = "Study") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors)  # Use custom colors

dev.off()  # Close the PDF device
  # Close the PDF device




library(ggplot2)
library(dplyr)

# Assuming r is your dataframe containing the data
# Create the contingency table
contingency_table <- table(r$Study, r$seurat_clusters)

# Convert the table to a data frame for easier plotting
df <- as.data.frame(contingency_table)
colnames(df) <- c("Study", "Cluster", "Count")

# Define a custom color palette with at least 13 colors
custom_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31"
)

# Plot the stacked bar chart showing actual cell counts and save it as a PDF
pdf("Cells_Per_Study_Per_Cluster_Count.pdf", height = 10, width = 15)

ggplot(df, aes(x = Cluster, y = Count, fill = Study)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Cells from Each Study in Each Cluster",
       x = "Cluster",
       y = "Number of Cells",
       fill = "Study") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors)  # Use custom colors

dev.off()  # Close the PDF device

----------------------------------------------------
## to assign a cluster identity, i am finding which genes are deg in one cluster vs rest

##A positive avg_log2FC (+ log2 fold change) indicates that the gene is more highly expressed in the cells
#of the cluster being compared (the "ident.1" group) relative to the rest of the cells (the "ident.2" group).
#A negative avg_log2FC (- log2 fold change) indicates that the gene is less expressed in the cells of the cluster
#being compared relative to the rest of the cells.
dir.create("differential_genes", showWarnings = FALSE)

# Loop through each cluster and save filtered and sorted results
for (cluster in unique(r$seurat_clusters)) {
  # Find markers for the cluster compared to all other cells
  markers <- FindMarkers(r, ident.1 = cluster, min.pct = 0.0, logfc.threshold = 0.5)
  # Filter markers by p-value
  filtered_markers <- markers %>%
    filter(p_val < 0.05) %>%
    arrange(desc(avg_log2FC))
  # Save the filtered and sorted markers to a CSV file
  write.csv(filtered_markers, file = paste0("differential_genes/cluster_", cluster, "_markers.csv"))
}


------------------------------------------------------

Idents(mic) <- mic$Sex
# subset_data <- subset(mic, dx %in% c("Control", "Dementia", "MCI"))

control1 <- subset(mic, dx %in% c("Control"))
control1 <- NormalizeData(control1, normalization.method = "LogNormalize", scale.factor = 10000)

Dementia1 <- subset(mic, dx %in% c("Dementia"))
Dementia1 <- NormalizeData(Dementia1, normalization.method = "LogNormalize", scale.factor = 10000)

MCI1 <- subset(mic, dx %in% c("MCI"))
MCI1 <- NormalizeData(MCI1, normalization.method = "LogNormalize", scale.factor = 10000)

# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "Control" group compared to the "Dementia" group.

#for inhibiory cells, normalizing was done

# Control vs Dementia
control_mic_male_vs_female1 <- FindMarkers(control1, ident.1 = "Male", ident.2 = "Female", test.use ="MAST", min.pct = 0,logfc.threshold = 0)
write.csv(control_mic_male_vs_female1, file = "control_mic_male_vs_female_normalized.csv")
head(order(control_mic_male_vs_female1$avg_log2FC),)

sorted_data <- control_microglia_male_vs_female[order(control_microglia_male_vs_female$avg_log2FC, decreasing = TRUE), ]


# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "Control" group compared to the "MCI" group.
# Control vs MCI
Dementia_mic_male_vs_female1 <- FindMarkers(Dementia1, ident.1 = "Male", ident.2 = "Female", test.use ="MAST", min.pct = 0,logfc.threshold = 0)
write.csv(Dementia_mic_male_vs_female1, file = "Dementia_ast_male_vs_female_normalized.csv")


# Positive Fold Change: A positive fold change for a gene means that this gene is more highly expressed in the "MCI" group compared to the "Dementia" group.
# MCI vs Dementia
MCI_mic_male_vs_female1 <- FindMarkers(MCI1, ident.1 = "Male", ident.2 = "Female", test.use ="MAST", min.pct = 0,logfc.threshold = 0)
write.csv(MCI_mic_male_vs_female1, file = "MCI_ast_male_vs_female_normalized.csv")


# Print the top results for each comparison
print(head(control_vs_dementia))
print(head(control_vs_mci))
print(head(mci_vs_dementia))

# Optionally, save the results to CSV files
write.csv(Dementia_microglia_male_vs_female, file = "Dementia_microglia_male_vs_female.csv")
write.csv(MCI_microglia_male_vs_female, file = "MCI_microglia_male_vs_female.csv")


----------------------------------------------------------------------------

mic_control_male_female1 <- read.csv("control_mic_male_vs_female_normalized.csv")
mic_dementia_male_female1 <- read.csv("Dementia_MIC_male_vs_female_normalized.csv")
mic_mci_male_female1 <- read.csv("MCI_mic_male_vs_female_normalized.csv")

----------------------------------------------------------------------------

control_ast_male_vs_female1$gene <- rownames(control_ast_male_vs_female1)
Dementia_ast_male_vs_female1$gene <- rownames(Dementia_ast_male_vs_female1)
MCI_mic_male_vs_female1$gene <- rownames(MCI_mic_male_vs_female1)

# Filter results for p_val less than 0.05 and sort based on avg_log2FC
mic_control_male_female2 <- control_mic_male_vs_female1 %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
mic_control_male_female2$gene <- mic_control_male_female2$X


mic_dementia_male_female2 <- Dementia_mic_male_vs_female1 %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
mic_dementia_male_female2$gene <- mic_dementia_male_female2$X


mic_mci_male_female2 <- MCI_mic_male_vs_female1 %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
mic_mci_male_female2$gene <- mic_mci_male_female2$X

# Print the top results for each comparison
print(head(mic_control_male_female2))
print(head(mic_dementia_male_female2))
print(head(mic_mci_male_female2))

# Optionally, save the results to CSV files
write.csv(mic_control_male_female2, file = "ast_control_male_female2_padj0.05.csv")
write.csv(mic_dementia_male_female2, file = "ast_dementia_male_female2_padj0.05.csv")
write.csv(mic_mci_male_female2, file = "ast_mci_male_female2_padj0.05.csv")


library(enrichR)


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive) {
    enriched <- enrichr(mic_mci_male_female2[ast_mci_male_female2$avg_log2FC < 0, "gene"], dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("down_ENRICH_mic_mci_male_female_normalized_upinfemale.pdf", width=10, height=10)
if (websiteLive) {
    plotEnrich(enriched[[3]], showTerms = 10, numChar = 60, y = "Ratio", orderBy = "Combined.Score", title = "GO_Biological_Process_2023")
}
dev.off()


library("EnhancedVolcano")

# mic_control_male_female2$gene <- rownames(mic_control_male_female2)

pdf("mic_control_male_female_volcanoplot_fc0.5_normalized.pdf")
EnhancedVolcano(mic_control_male_female2, lab = ast_control_male_female2$gene,x = 'avg_log2FC',y = 'p_val_adj', title="male  vs female", subtitle = "", FCcutoff = 0.5)
dev.off()


# control_vs_mci2$gene <- rownames(control_vs_mci2)
pdf("mic_dementia_male_female2_fc0.5_normalized.pdf")
EnhancedVolcano(mic_dementia_male_female2, lab = ast_dementia_male_female2$gene,x = 'avg_log2FC',y = 'p_val_adj', title="male vs female", subtitle = "", FCcutoff = 0.5)
dev.off()


# mci_vs_dementia2$gene <- rownames(mci_vs_dementia2)
pdf("mic_mci_male_female_volcanoplot_fc0.5_normalized.pdf")
EnhancedVolcano(mic_mci_male_female2, lab = ast_mci_male_female2$gene,x = 'avg_log2FC',y = 'p_val_adj', title="male vs female", subtitle = "", FCcutoff = 0.5, labSize = 4.0)
dev.off()
