#!/usr/bin/env Rscript

# @theme: T5 - Joint EC preprocessing, integration, and legacy score/UMAP exports
# @task_id: T5-02
# @description: Step-by-step script for reproducing the key cross-species brain EC preprocessing, QC, integration, endothelial refinement, cell-cycle scoring, and subtype-labeling workflow used in the manuscript.
# @inputs:
# - `path/to/brain_input_matrices/*.h5`: example filtered 10x gene-expression matrices for the late-stage mouse and guinea pig brain endothelial comparison.
# @outputs:
# - `path/to/brain_ec_public_repro/`: example release directory for checkpoint objects, QC plots, summary tables, and labeled UMAP exports.
# - `path/to/brain_ec_public_repro/example_sample_manifest.tsv`: example sample manifest used by the walkthrough.
# @key_params:
# - `expected_doublet_rate = 0.025`: expected per-sample doublet rate passed to `scDblFinder`.
# - `whole_brain_integration_features = 2000`: shared features used for the initial cross-species integration.
# - `endothelial_integration_features = 2000`: shared features used for endothelial-only reintegration.
# - `endothelial_doublet_score_max = 0.20`: upper bound applied during endothelial refinement.
# @dependencies:
# - `Seurat` (v4.4.0): preprocessing, integration, dimensional reduction, clustering, UMAP, and cell-cycle scoring.
# - `SeuratObject` (v5.0.1): Seurat object infrastructure used by the retained environment snapshot.
# - `scDblFinder` (v1.17.2): cluster-aware doublet scoring.
# - `BiocParallel` (v1.32.5): multicore execution for doublet scoring.
# - `gprofiler2` (v0.2.3): retained in the release environment used for manuscript support.
# - `ggplot2` (v3.5.0): visualization.
# - `patchwork` (v1.2.0): plot composition.
# - `viridis` (v0.6.5): color scaling for manuscript-supporting plots.
# @examples:
# Rscript path/to/reproduce_cross_species_brain_ec_pipeline.R
#
# Run this script block by block in the Seurat 4 environment. Each section saves
# the intermediate objects, tables, or plots needed for the next stage of the
# late-stage brain endothelial comparison. Replace the example paths and example
# file names below with the locations used in your own release package.
#
# Recorded environment snapshot from the seurat4 conda environment on 2026-03-29:
# - R: 4.2.0
# - platform: x86_64-conda-linux-gnu (64-bit)
# - running: CentOS Linux 7 (Core)
# - Seurat: 4.4.0
# - SeuratObject: 5.0.1
# - scDblFinder: 1.17.2
# - BiocParallel: 1.32.5
# - gprofiler2: 0.2.3
# - ggplot2: 3.5.0
# - patchwork: 1.2.0
# - viridis: 0.6.5

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(BiocParallel)
  library(gprofiler2)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

###############################################################################
# STEP 1. Define the input files, stage labels, and analysis parameters.
#
# The stage labels are used throughout the manuscript-facing outputs. Generic
# raw input IDs are retained as provenance placeholders.
###############################################################################

input_dir <- "path/to/brain_input_matrices"
output_dir <- "path/to/brain_ec_public_repro"
workers <- 4

sample_sheet <- data.frame(
  sample_key = c(
    "sample_01",
    "sample_02",
    "sample_03",
    "sample_04"
  ),
  sample_stage = c(
    "Mouse E18",
    "Mouse P0",
    "Guinea pig GD35 replicate 1",
    "Guinea pig GD35 replicate 2"
  ),
  legacy_input_id = c("raw_id_01", "raw_id_02", "raw_id_03", "raw_id_04"),
  species_code = c("mm", "mm", "gp", "gp"),
  species = c("Mouse", "Mouse", "Guinea pig", "Guinea pig"),
  input_file = c(
    file.path(input_dir, "mouse_stage_01_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "mouse_stage_02_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "guinea_pig_stage_01_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "guinea_pig_stage_02_filtered_feature_bc_matrix.h5")
  ),
  stringsAsFactors = FALSE
)

# The second mouse dataset is labeled as the P0 sample throughout this release
# script; the example input filename above is only a placeholder for the file
# you plan to distribute with the manuscript support package.

mouse_min_features <- 1000
mouse_min_counts <- 1000
mouse_max_percent_mt <- 12.5

gp_min_features <- 1000
gp_min_counts <- 1000
gp_min_percent_ribo <- 5

expected_doublet_rate <- 0.025
sample_pcs <- 30
sample_resolution <- 0.4

whole_brain_gene_min_cells <- 3
whole_brain_pcs <- 30
whole_brain_resolution <- 0.4
whole_brain_integration_features <- 2000
whole_brain_endothelial_clusters <- c("0", "1", "2", "3", "13")

endothelial_pcs <- 20
endothelial_resolution <- 0.4
endothelial_integration_features <- 2000
endothelial_doublet_score_max <- 0.20

bec_min_features <- 2000
bec_min_counts <- 3000
bec_min_percent_ribo <- 5
bec_pcs <- 30
bec_neighbor_dims <- 1:15
bec_resolution <- 0.4

cell_markers <- c(
  "Pecam1", "Cdh5", "Tek", "Cspg4", "Kcnj8", "Rbfox3", "Ptprc",
  "Aif1", "Mog", "Gfap", "Adgre1", "Neurod1", "Itga8", "Tagln"
)
fate_markers <- c(
  "Nrp2", "Vwf", "Aplnr", "Nr2f2", "Apln", "Meox1",
  "Dll4", "Igfbp3", "Gja4", "Vegfc", "Gja5", "Cxcr4"
)

fine_label_map <- c(
  "0" = "Cap3",
  "1" = "Vein",
  "2" = "Pre-artery",
  "3" = "Cap1",
  "4" = "Cycling",
  "5" = "Cap2",
  "6" = "Cycling",
  "7" = "Artery"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "qc_plots"), recursive = TRUE, showWarnings = FALSE)

missing_inputs <- sample_sheet$input_file[!file.exists(sample_sheet$input_file)]
if (length(missing_inputs) > 0) {
  stop(
    sprintf("Missing input files:\n- %s", paste(missing_inputs, collapse = "\n- ")),
    call. = FALSE
  )
}

write.table(
  sample_sheet,
  file = file.path(output_dir, "example_sample_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 2. Build one Seurat object per sample in its native feature space.
#
# Each sample is normalized and clustered first so QC thresholds and doublet
# scores are evaluated on the same sample-specific structure used in the
# original brain preprocessing workflow.
###############################################################################

sample_objects <- list()
sample_qc_summary <- vector("list", nrow(sample_sheet))
bp <- MulticoreParam(workers, RNGseed = 1234)

for (i in seq_len(nrow(sample_sheet))) {
  sample_key <- sample_sheet$sample_key[[i]]
  sample_stage <- sample_sheet$sample_stage[[i]]
  legacy_input_id <- sample_sheet$legacy_input_id[[i]]
  species_code <- sample_sheet$species_code[[i]]
  species_name <- sample_sheet$species[[i]]
  input_file <- sample_sheet$input_file[[i]]

  cat(sprintf("Reading %s (example input ID: %s) from %s\n", sample_stage, legacy_input_id, input_file))
  counts <- Read10X_h5(input_file)
  if (inherits(counts, "list")) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  input_cells <- ncol(counts)

  seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3)
  seurat_obj$sample_key <- sample_key
  seurat_obj$sample_stage <- sample_stage
  seurat_obj$legacy_input_id <- legacy_input_id
  seurat_obj$samples <- sample_stage
  seurat_obj$species <- species_name
  seurat_obj$species_code <- species_code

  # Mitochondrial and ribosomal fractions are computed for all samples, even
  # though the retained mouse and guinea pig filters use different metrics.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(Rp[ls]|RP[LS]|Rpl|Rps)")

  # Native-space normalization and clustering define the sample structure used
  # for preliminary QC inspection and cluster-aware doublet scoring.
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = sample_pcs, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:sample_pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = sample_resolution, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:sample_pcs, verbose = FALSE)

  # Sample-level QC thresholds remove low-complexity mouse libraries and low-
  # ribosomal-content guinea pig libraries exactly as encoded in the retained
  # preprocessing script.
  if (species_code == "mm") {
    keep_cells <- with(
      seurat_obj@meta.data,
      nFeature_RNA > mouse_min_features &
        nCount_RNA > mouse_min_counts &
        percent.mt < mouse_max_percent_mt
    )
  } else {
    keep_cells <- with(
      seurat_obj@meta.data,
      nFeature_RNA > gp_min_features &
        nCount_RNA > gp_min_counts &
        percent.ribo > gp_min_percent_ribo
    )
  }

  seurat_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[keep_cells])

  # After filtering, the neighborhood graph, cluster labels, and UMAP are
  # recomputed on the retained cells so the saved sample objects reflect the
  # post-QC cell population carried into integration.
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:sample_pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = sample_resolution, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:sample_pcs, verbose = FALSE)

  sce <- scDblFinder(
    GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
    clusters = Idents(seurat_obj),
    BPPARAM = bp,
    dbr = expected_doublet_rate
  )
  seurat_obj$scDblFinder.score <- sce$scDblFinder.score
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class

  sample_qc_plot <- (
    VlnPlot(
      seurat_obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
      ncol = 4,
      pt.size = 0
    ) + NoLegend()
  ) / (
    DimPlot(seurat_obj, reduction = "umap", group.by = "samples") + ggtitle(sample_stage)
  )
  ggsave(
    filename = file.path(output_dir, "qc_plots", paste0(sample_key, "_qc_panel.pdf")),
    plot = sample_qc_plot,
    width = 11,
    height = 7,
    units = "in"
  )

  saveRDS(seurat_obj, file.path(output_dir, paste0(sample_key, "_example_post_qc.rds")))

  sample_qc_summary[[i]] <- data.frame(
    sample_key = sample_key,
    sample_stage = sample_stage,
    legacy_input_id = legacy_input_id,
    species = species_name,
    input_cells = input_cells,
    post_qc_cells = ncol(seurat_obj),
    predicted_singlets_post_qc = sum(seurat_obj$scDblFinder.class == "singlet", na.rm = TRUE),
    predicted_doublets_post_qc = sum(seurat_obj$scDblFinder.class == "doublet", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  sample_objects[[sample_key]] <- seurat_obj
}

sample_qc_summary <- do.call(rbind, sample_qc_summary)
write.table(
  sample_qc_summary,
  file = file.path(output_dir, "example_sample_qc_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(sample_objects, file.path(output_dir, "example_post_qc_sample_objects.rds"))

###############################################################################
# STEP 3. Merge mouse and guinea pig samples separately, then restrict both
# species to the shared feature space used for cross-species whole-brain
# integration.
#
# The retained brain preprocessing script defines the matched feature set by the
# intersection of genes present in the mouse and guinea pig merged objects.
###############################################################################

mouse_sample_keys <- sample_sheet$sample_key[sample_sheet$species_code == "mm"]
gp_sample_keys <- sample_sheet$sample_key[sample_sheet$species_code == "gp"]

MM <- merge(
  x = sample_objects[[mouse_sample_keys[[1]]]],
  y = sample_objects[[mouse_sample_keys[[2]]]],
  add.cell.ids = mouse_sample_keys
)
GP <- merge(
  x = sample_objects[[gp_sample_keys[[1]]]],
  y = sample_objects[[gp_sample_keys[[2]]]],
  add.cell.ids = gp_sample_keys
)

features_MM <- rownames(MM)
features_GP <- rownames(GP)
shared_features <- intersect(features_MM, features_GP)

MM <- subset(MM, features = shared_features)
GP <- subset(GP, features = shared_features)

brain_by_sample <- c(
  SplitObject(MM, split.by = "sample_key"),
  SplitObject(GP, split.by = "sample_key")
)

# FilterGenes is retained here because the original script removes genes detected
# in too few cells before cross-sample integration.
FilterGenes <- function(seurat_obj, min_cells = 3) {
  gene_counts <- rowSums(seurat_obj@assays$RNA@data > 0)
  genes_to_keep <- names(gene_counts[gene_counts >= min_cells])
  subset(seurat_obj, features = genes_to_keep)
}

for (sample_id in names(brain_by_sample)) {
  brain_by_sample[[sample_id]] <- FilterGenes(brain_by_sample[[sample_id]], min_cells = whole_brain_gene_min_cells)
  brain_by_sample[[sample_id]] <- FindVariableFeatures(brain_by_sample[[sample_id]], verbose = FALSE)
}

whole_brain_features <- SelectIntegrationFeatures(
  object.list = brain_by_sample,
  nfeatures = whole_brain_integration_features
)

for (sample_id in names(brain_by_sample)) {
  brain_by_sample[[sample_id]] <- ScaleData(
    brain_by_sample[[sample_id]],
    features = whole_brain_features,
    verbose = FALSE
  )
  brain_by_sample[[sample_id]] <- RunPCA(
    brain_by_sample[[sample_id]],
    features = whole_brain_features,
    npcs = whole_brain_pcs,
    verbose = FALSE
  )
}

whole_brain_anchors <- FindIntegrationAnchors(
  object.list = brain_by_sample,
  anchor.features = whole_brain_features,
  dims = 1:whole_brain_pcs
)

whole_brain <- IntegrateData(anchorset = whole_brain_anchors, dims = 1:whole_brain_pcs)
DefaultAssay(whole_brain) <- "integrated"
whole_brain <- ScaleData(whole_brain, verbose = FALSE)
whole_brain <- RunPCA(whole_brain, npcs = whole_brain_pcs, verbose = FALSE)
whole_brain <- FindNeighbors(whole_brain, reduction = "pca", dims = 1:whole_brain_pcs, verbose = FALSE)
whole_brain <- FindClusters(whole_brain, resolution = whole_brain_resolution, verbose = FALSE)
whole_brain <- RunUMAP(whole_brain, reduction = "pca", dims = 1:whole_brain_pcs, verbose = FALSE)

saveRDS(whole_brain, file.path(output_dir, "example_whole_brain_integrated.rds"))

whole_brain_plot <- (
  DimPlot(whole_brain, reduction = "umap", group.by = "samples") +
    ggtitle("Whole-brain integrated atlas by sample")
) / (
  DimPlot(whole_brain, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Whole-brain integrated atlas by cluster")
)
ggsave(
  filename = file.path(output_dir, "qc_plots", "example_whole_brain_integrated_umap.pdf"),
  plot = whole_brain_plot,
  width = 8,
  height = 10,
  units = "in"
)

###############################################################################
# STEP 4. Isolate the endothelial compartment from the whole-brain atlas and
# reintegrate the retained endothelial cells across samples.
#
# The first integrated atlas is used only to locate the endothelial clusters.
# Those clusters are then re-integrated on their own to sharpen endothelial
# state structure before subtype annotation.
###############################################################################

Idents(whole_brain) <- "seurat_clusters"
ec <- subset(whole_brain, idents = whole_brain_endothelial_clusters)
DefaultAssay(ec) <- "RNA"
ec <- FilterGenes(ec, min_cells = whole_brain_gene_min_cells)

DefaultAssay(ec) <- "integrated"
ec <- ScaleData(ec, verbose = FALSE)
ec <- RunPCA(ec, npcs = endothelial_pcs, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = "pca", dims = 1:endothelial_pcs, verbose = FALSE)
ec <- FindClusters(ec, resolution = endothelial_resolution, verbose = FALSE)
ec <- RunUMAP(ec, reduction = "pca", dims = 1:endothelial_pcs, verbose = FALSE)

ggsave(
  filename = file.path(output_dir, "qc_plots", "example_endothelial_subset_before_reintegration_umap.pdf"),
  plot = DimPlot(ec, reduction = "umap", group.by = "samples") +
    ggtitle("Endothelial clusters before endothelial-only reintegration"),
  width = 8,
  height = 6,
  units = "in"
)

ec_by_sample <- SplitObject(ec, split.by = "sample_key")
for (sample_id in names(ec_by_sample)) {
  DefaultAssay(ec_by_sample[[sample_id]]) <- "RNA"
  ec_by_sample[[sample_id]] <- FindVariableFeatures(ec_by_sample[[sample_id]], verbose = FALSE)
}

endothelial_features <- SelectIntegrationFeatures(
  object.list = ec_by_sample,
  nfeatures = endothelial_integration_features
)

for (sample_id in names(ec_by_sample)) {
  ec_by_sample[[sample_id]] <- ScaleData(
    ec_by_sample[[sample_id]],
    features = endothelial_features,
    verbose = FALSE
  )
  ec_by_sample[[sample_id]] <- RunPCA(
    ec_by_sample[[sample_id]],
    features = endothelial_features,
    npcs = endothelial_pcs,
    verbose = FALSE
  )
}

endothelial_anchors <- FindIntegrationAnchors(
  object.list = ec_by_sample,
  anchor.features = endothelial_features,
  dims = 1:endothelial_pcs
)

ec <- IntegrateData(anchorset = endothelial_anchors, dims = 1:endothelial_pcs)
DefaultAssay(ec) <- "integrated"
ec <- ScaleData(ec, verbose = FALSE)
ec <- RunPCA(ec, npcs = endothelial_pcs, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = "pca", dims = 1:endothelial_pcs, verbose = FALSE)
ec <- FindClusters(ec, resolution = endothelial_resolution, verbose = FALSE)
ec <- RunUMAP(ec, reduction = "pca", dims = 1:endothelial_pcs, verbose = FALSE)
ec <- subset(ec, subset = scDblFinder.score < endothelial_doublet_score_max)

saveRDS(ec, file.path(output_dir, "example_endothelial_integrated_filtered.rds"))

###############################################################################
# STEP 5. Score endothelial cell-cycle phase and define the higher-stringency
# brain EC subset used for final subtype labeling.
#
# Cell-cycle scoring is performed on the endothelial object after reintegration.
# The higher-stringency BEC subset keeps the cells that passed the additional
# library-size and ribosomal-content thresholds used in the original workflow.
###############################################################################

DefaultAssay(ec) <- "RNA"
m.s.genes <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
m.g2m.genes <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

ec <- CellCycleScoring(
  ec,
  s.features = unique(m.s.genes),
  g2m.features = unique(m.g2m.genes),
  set.ident = TRUE
)

bec <- subset(
  ec,
  subset =
    nFeature_RNA > bec_min_features &
      nCount_RNA > bec_min_counts &
      percent.ribo > bec_min_percent_ribo
)

DefaultAssay(bec) <- "integrated"
bec <- ScaleData(bec, verbose = FALSE)
bec <- RunPCA(bec, npcs = bec_pcs, verbose = FALSE)
bec <- FindNeighbors(bec, reduction = "pca", dims = bec_neighbor_dims, verbose = FALSE)
bec <- FindClusters(bec, resolution = bec_resolution, verbose = FALSE)
bec <- RunUMAP(bec, reduction = "pca", dims = bec_neighbor_dims, verbose = FALSE)

saveRDS(bec, file.path(output_dir, "example_brain_ec_subset_before_labeling.rds"))

bec_plot <- (
  DimPlot(bec, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Brain EC subset by cluster")
) / (
  DimPlot(bec, reduction = "umap", group.by = "samples") +
    ggtitle("Brain EC subset by sample")
)
ggsave(
  filename = file.path(output_dir, "qc_plots", "example_brain_ec_subset_umap.pdf"),
  plot = bec_plot,
  width = 8,
  height = 10,
  units = "in"
)

###############################################################################
# STEP 6. Assign the final endothelial subtype labels.
#
# Cluster 8 is excluded before label transfer, following the retained script.
# The remaining clusters are labeled as artery, pre-artery, capillary, cycling,
# or venous endothelial states based on the original marker-guided annotation.
###############################################################################

Idents(bec) <- "seurat_clusters"
if ("8" %in% levels(Idents(bec))) {
  bec <- subset(bec, idents = "8", invert = TRUE)
}

present_clusters <- levels(Idents(bec))
if (!all(present_clusters %in% names(fine_label_map))) {
  stop(
    sprintf(
      "The final BEC clustering produced unmapped clusters: %s",
      paste(setdiff(present_clusters, names(fine_label_map)), collapse = ", ")
    ),
    call. = FALSE
  )
}

label_arguments <- as.list(unname(fine_label_map[present_clusters]))
names(label_arguments) <- present_clusters
bec <- do.call(RenameIdents, c(list(object = bec), label_arguments))
bec$cell.type.l1 <- Idents(bec)

saveRDS(bec, file.path(output_dir, "example_brain_ec_final_labeled.rds"))

label_counts <- as.data.frame(table(bec$samples, bec$cell.type.l1))
colnames(label_counts) <- c("sample_stage", "cell_type", "n_cells")
write.table(
  label_counts,
  file = file.path(output_dir, "example_brain_ec_label_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = file.path(output_dir, "qc_plots", "example_brain_ec_final_labels_umap.pdf"),
  plot = (
    DimPlot(bec, reduction = "umap", group.by = "cell.type.l1", label = TRUE) +
      ggtitle("Final brain EC subtype labels")
  ) / (
    FeaturePlot(
      bec,
      features = fate_markers,
      order = TRUE,
      cols = viridis(256, alpha = 0.8, begin = 0.1, option = "C")
    ) + plot_annotation(title = "Marker support for final labels")
  ),
  width = 11,
  height = 12,
  units = "in"
)

cat("Finished cross-species brain EC workflow template.\n")
cat(sprintf("Configured example output directory: %s\n", output_dir))
