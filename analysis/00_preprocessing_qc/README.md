# Preprocessing, Integration & Annotation (Fig. 1)

Atlas-construction pipeline behind Fig. 1: 13 published snRNA-seq cohorts, harmonized to Control/MCI/AD, integrated and annotated into 7 major cell classes / 33 substates, with a post-hoc rescue of a missed endothelial population as its own cluster.

**Per-cohort processing is not included** (CellBender ambient-RNA removal, Scanpy QC at >=200 genes/nucleus, >=250 UMI and <=15% MT, Scrublet doublet scores). This module starts from the 13 cohorts after that processing.

## Integration, label transfer and clustering (script 00)

`00_scvi_scanvi_allen_reference_mapping.py` maps the merged 13-cohort data onto the **Allen Human Multiple Cortical Areas SMART-seq** reference (49,417 nuclei, 3 donors, 19 subclasses) with scVI/scANVI (scArches-style reference mapping), then clusters the scANVI latent space with Leiden:

- Reference: 5,000 HVGs selected on genes shared with the query; scVI (`n_latent=10`, `n_layers=2`, `dropout_rate=0.2`, layer norm, `encode_covariates=True`, 162 epochs) -> scANVI on `subclass_label` (30 epochs)
- Query: 4,807,280 nuclei from 913 samples, mapped with `load_query_data` (30 epochs, `weight_decay=0`); **each sample is its own batch** (`external_donor_name_label := Sample`; 916 batches including the 3 Allen donors)
- Neighbors (cosine) on `X_scANVI` -> UMAP -> **Leiden, resolution 0.3** -> 18 clusters (0-17). Cluster 17 is a small residual cluster (9,412 cells, 8 by the final atlas) that is removed in script 07. That leaves **17 Leiden clusters (0-16)**, which are annotated into 7 major classes (EXC, INH, ODC, OPC, ASC, MIC, ENDO).
- Doublet filter: `doublet_score <= 0.1` -> 3,556,850 nuclei

The script has two modes:

- `--mode original` reproduces the run behind the manuscript (the settings above; scvi-tools 0.20.3). With `--pretrained` it loads the saved models (`scvi_model_allen_smartseq_03-11-24`, `scanvi_model_allen_03-11-24`, `scanvi_query_model_03-11-24`) and reproduces the published latent space. The query model's 5,000 genes and all 873 sample batches match the final atlas. Retraining in this mode gives a similar but not bit-identical embedding, because no random seed was fixed in the original run.
- `--mode tutorial` (default) follows the scvi-tools 1.3.3 tutorials ([Tabula Muris label transfer](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/tabula_muris.html), [scArches reference mapping](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/scarches_scvi_tools.html)). It differs from the original run in these ways: the SMART-seq reference counts are gene-length normalized (`--gtf`); HVGs use `flavor="seurat_v3"` on counts; reference cells without a subclass are left unlabeled (`Unknown`) instead of forming an `unknown` class; scANVI uses 20 epochs with `n_samples_per_label=100`; the query is trained for 100 epochs; doublets are removed before clustering; the seed is fixed; and the script writes soft prediction probabilities and a reference self-prediction confusion matrix. A run in this mode gives a new embedding, so the Leiden cluster -> cell type map in script 00b must be re-checked against the marker dot plot.

## Metadata and major cell-type annotation (script 00b)

`00b_metadata_and_celltype_annotation.py` continues from script 00:

- Joins the per-cohort sample metadata (ADDS samples prefixed `ADDS_`), sets `Diagnosis2` (NCI -> Control) and drops cells with diagnosis "Not Applicable" (-> 3,316,964 nuclei)
- Labels every cell twice: `Primary_Function` from the scANVI Allen subclass prediction, and `Celltype_markers` from the Leiden cluster annotated by canonical markers. Cells the label transfer left as "unknown" take the marker label.
- Keeps only cells where the two labels agree (-> 3,158,097 nuclei) and writes one h5ad per major class

**Between script 00b and script 01 (not scripted):** subclustering of the per-class objects into the 33 substates, done interactively (pre-ENDO atlas: 3,146,108 nuclei). The output is `13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad`, the input to script 01.

## ENDO rescue (scripts 01-07)

An internal step in the same workflow: a small batch of endothelial cells missed by the initial annotation is folded in as its own cluster (labeled 18 during the rescue, renamed to 17 in script 07) so the atlas ends up with a complete ENDO population, keeping the merge-through-final-atlas sequence consistent end to end. Run order, with the input and output files of each step:

| # | Script | Input | Output |
|---|--------|-------|--------|
| 0 | `00_scvi_scanvi_allen_reference_mapping.py` | Allen SMART-seq reference + merged 13-cohort h5ad | `13datasets_after_label_transfer_leiden_doublet0.1.h5ad` |
| 0b | `00b_metadata_and_celltype_annotation.py` | output of (0) + per-cohort metadata CSVs | `13AD_matched_cellmarker_labelt.h5ad`, `13ad_{Celltype}.h5ad` |
| 1 | `01_endo_rescue_initial_qc.py` | initial-QC metrics h5ad + the finished pre-ENDO atlas | `13ad_endo_rescued_from_initial_qc.h5ad` |
| 2 | `02_merge_endo_rescue_into_atlas.py` | atlas + rescued-ENDO h5ad from (1) | `13ad_3m_5k_atlas_plus_endo.h5ad` |
| 3 | `03_rename_endo_rescued_to_endo.py` | atlas+endo h5ad from (2), in place | relabels `ENDO_rescued_qc` -> `ENDO` in all annotation columns |
| 4 | `04_finalize_endo_as_cluster18.py` | atlas+endo h5ad from (3) | `13ad_3m_5k_atlas_with_endo_final.h5ad` |
| 5 | `05_project_endo_scanvi_umap.py` | output of (4) | `13ad_3m_5k_atlas_with_endo_scanvi.h5ad` |
| 6 | `06_recompute_atlas_umap_with_endo.py` | output of (5) | `13ad_3m_5k_atlas_with_endo_umap_recomputed.h5ad` |
| 7 | `07_finalize_remove_cl17_rename_cl18.py` | output of (6) | `13ad_3m_5k_atlas_final.h5ad` (final atlas used by all downstream analyses) |

## Figure scripts (Fig. 1b, Suppl. Fig. 1)

| Script | Output |
|--------|--------|
| `fig1b_umap_substates.py` | Fig. 1b UMAP of the final atlas coloured by the 33 substates |
| `fig1b_pie_charts_cells_individuals.py` | Fig. 1b pie charts of major cell types and donors (nuclei vs. individuals) |
| `fig1b_dx_piechart.py` | Fig. 1b diagnosis pie charts (nuclei and individuals) |
| `supp_fig1_qc_violins.py` | Suppl. Fig. 1 QC violin + box panels (total counts, genes, % MT, doublet score) by study |

Each reads the final atlas h5ad (file name in the script) from the working directory.

## Usage

`00_scvi_scanvi_allen_reference_mapping.py` takes its inputs as arguments (defaults are the original file names, resolved in the working directory; see `--help`). Training needs a GPU. Tutorial mode also needs a GTF with `gene_name` attributes, e.g. the Cell Ranger GRCh38-2020-A `genes.gtf`. Add `--smoke-test` (5k reference / 20k query cells, 1 epoch per model) to check that the environment and inputs work before a full run.

Environments: tutorial mode needs scvi-tools 1.3.3 (smoke-tested with torch 2.5.1 and CUDA 12.1 on an NVIDIA A30). `--mode original --pretrained` needs scvi-tools 0.20.3 with pandas < 2 (smoke-tested with torch 1.13.1, anndata 0.9.1, scanpy 1.9.3). The saved scANVI models do not load in scvi-tools >= 1.0, and the script stops with an explanation if you try.

Memory: the query is streamed from disk and only the 5,000 model genes are kept (float32), which keeps the query at roughly 25-30 GB in memory (estimated; loading the whole query h5ad at once needs > 250 GB). We request 300 GB for a full run. Under SLURM, request a single task (`--ntasks=1`). The script also guards against Lightning's multi-task SLURM detection.

Scripts 01-07 read/write relative to the current working directory (`BASE = os.getcwd()` / `Path.cwd()` at the top of each file) -- run from the repository root with `data/h5ad/` populated, in the numbered order above.

## Inputs

- `allen_countmat.h5ad` -- Allen Human Multiple Cortical Areas SMART-seq reference (raw counts), input to script 00
- `13ad_filtered_samples200&15000.h5ad` -- the 13 cohorts merged after per-cohort CellBender/QC/Scrublet, input to script 00
- `metadata/*.csv` -- per-cohort sample metadata (Samples, Diagnosis, Sex, Age, ...), input to script 00b
- `data/h5ad/13ad_3m_5k_subcluster_annot_labeled_renamed_neurons.h5ad` -- annotated pre-ENDO atlas (script 00b output after the interactive subclustering above), input to script 01
- `data/h5ad/13data_qcmetrics_initial_qc.h5ad` -- initial QC metrics used to identify rescuable ENDO cells

## Outputs

- `data/h5ad/13ad_3m_5k_atlas_final.h5ad` -- the final annotated atlas (7 classes, 33 substates; Leiden clusters 0-16 plus ENDO as cluster 17) consumed by every other analysis module in this repository
