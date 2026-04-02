#!/usr/bin/env Rscript

# @theme: T5 - Joint EC preprocessing, integration, and legacy score/UMAP exports
# @task_id: T5-01
# @description: Step-by-step script for reproducing the key cross-species heart EC preprocessing, ortholog remapping, integration, endothelial refinement, and label-assignment workflow used in the manuscript.
# @inputs:
# - `path/to/heart_input_matrices/*.h5`: example filtered 10x gene-expression matrices for three mouse and three guinea pig heart samples.
# - `path/to/ortholog_table.tsv`: example guinea pig-to-mouse ortholog table used for one-to-one remapping before integration.
# @outputs:
# - `path/to/heart_ec_public_repro/`: example release directory for checkpoint objects, QC plots, summary tables, and labeled UMAP exports.
# - `path/to/heart_ec_public_repro/example_sample_manifest.tsv`: example sample manifest used by the walkthrough.
# @key_params:
# - `expected_doublet_rate = 0.025`: expected per-sample doublet rate passed to `scDblFinder`.
# - `whole_heart_integration_features = 2000`: shared features used for whole-heart integration.
# - `endo_integration_features = 1000`: shared features used for endothelial-only reintegration.
# - `endo_neighbor_dims = 1:7`: PCA dimensions used for endothelial reclustering.
# @dependencies:
# - `Seurat` (v4.4.0): preprocessing, integration, dimensional reduction, clustering, UMAP, and cell-cycle scoring.
# - `SeuratObject` (v5.0.1): Seurat object infrastructure used by the retained environment snapshot.
# - `scDblFinder` (v1.17.2): cluster-aware doublet scoring.
# - `BiocParallel` (v1.32.5): multicore execution for doublet scoring.
# - `gprofiler2` (v0.2.3): ortholog-table support utilities retained in the release environment.
# - `ggplot2` (v3.5.0): visualization.
# - `patchwork` (v1.2.0): plot composition.
# @examples:
# Rscript path/to/reproduce_cross_species_heart_ec_pipeline.R
#
# Run this script block by block in the Seurat 4 environment. Each major section
# creates a checkpoint object, table, or figure that is used in the next stage
# of the heart endothelial-cell analysis. Replace the example paths and example
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

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(BiocParallel)
  library(gprofiler2)
  library(ggplot2)
  library(patchwork)
})

###############################################################################
# STEP 1. Define the input files, sample metadata, and analysis parameters.
#
# The sample sheet links stage labels to example 10x matrices. Downstream plots
# and metadata use the stage labels, while generic raw input IDs remain
# available as provenance placeholders.
###############################################################################

input_dir <- "path/to/heart_input_matrices"
ortholog_table_path <- "path/to/ortholog_table.tsv"
output_dir <- "path/to/heart_ec_public_repro"
workers <- 4

sample_sheet <- data.frame(
  sample_key = c(
    "sample_01",
    "sample_02",
    "sample_03",
    "sample_04",
    "sample_05",
    "sample_06"
  ),
  sample_stage = c(
    "Mouse E13.5",
    "Mouse E15.5",
    "Mouse E17.5",
    "Guinea pig GD25",
    "Guinea pig GD32",
    "Guinea pig GD35"
  ),
  legacy_input_id = c("raw_id_01", "raw_id_02", "raw_id_03", "raw_id_04", "raw_id_05", "raw_id_06"),
  species_code = c("mm", "mm", "mm", "gp", "gp", "gp"),
  species = c("Mouse", "Mouse", "Mouse", "Guinea pig", "Guinea pig", "Guinea pig"),
  input_file = c(
    file.path(input_dir, "mouse_stage_01_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "mouse_stage_02_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "mouse_stage_03_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "guinea_pig_stage_01_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "guinea_pig_stage_02_filtered_feature_bc_matrix.h5"),
    file.path(input_dir, "guinea_pig_stage_03_filtered_feature_bc_matrix.h5")
  ),
  stringsAsFactors = FALSE
)

# Sample-level cell filters applied before whole-heart integration.
mouse_min_features <- 750
mouse_min_counts <- 1000
mouse_max_percent_mt <- 10

# Guinea pig filtering in this workflow uses ribosomal-content burden together
# with library size and detected-gene thresholds.
gp_min_features <- 750
gp_min_counts <- 1000
gp_max_percent_ribo <- 20

# Shared settings for per-sample preprocessing, whole-heart integration, and
# endothelial-only reanalysis.
expected_doublet_rate <- 0.025
whole_heart_pcs <- 30
whole_heart_resolution <- 0.4
whole_heart_integration_features <- 2000
whole_heart_endothelial_clusters <- c("6")

endo_gene_min_cells <- 10
endo_pcs <- 20
endo_neighbor_dims <- 1:7
endo_resolution <- 0.4
endo_umap_neighbors <- 50
endo_integration_features <- 1000

# Cluster-to-label maps applied after endothelial reclustering.
fine_label_map <- c(
  "0" = "Cap4",
  "1" = "Cap1",
  "2" = "Arteriole",
  "3" = "Cap3",
  "4" = "Cap2",
  "5" = "Artery",
  "6" = "Vein",
  "7" = "Cap4"
)

coarse_label_map <- c(
  "0" = "Capillary",
  "1" = "Capillary",
  "2" = "Arteriole",
  "3" = "Capillary",
  "4" = "Capillary",
  "5" = "Artery",
  "6" = "Vein",
  "7" = "Capillary"
)

fate_markers <- c("Nrp2", "Vwf", "Aplnr", "Nr2f2", "Apln", "Bmx", "Cxcr4", "Gja4", "Gja5")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "qc_plots"), recursive = TRUE, showWarnings = FALSE)

missing_inputs <- sample_sheet$input_file[!file.exists(sample_sheet$input_file)]
if (length(missing_inputs) > 0) {
  stop(
    sprintf("Missing input files:\n- %s", paste(missing_inputs, collapse = "\n- ")),
    call. = FALSE
  )
}

if (!file.exists(ortholog_table_path)) {
  stop(sprintf("Ortholog table not found: %s", ortholog_table_path), call. = FALSE)
}

write.table(
  sample_sheet,
  file = file.path(output_dir, "example_sample_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 2. Import the guinea pig-to-mouse ortholog table used later to remap
# guinea pig features and define the matched feature space for integration.
#
# One-to-one ortholog pairs are carried forward as a lookup table during sample-
# level preprocessing and are applied immediately before cross-species integration.
###############################################################################

ortholog_table <- read.delim(ortholog_table_path, check.names = FALSE)
ortholog_table <- ortholog_table[
  ortholog_table[["Mouse homology type"]] == "ortholog_one2one",
  c("Gene name", "Mouse gene name"),
  drop = FALSE
]
ortholog_table <- ortholog_table[
  ortholog_table[["Gene name"]] != "" & ortholog_table[["Mouse gene name"]] != "",
  ,
  drop = FALSE
]
ortholog_table <- ortholog_table[!duplicated(ortholog_table[["Gene name"]]), , drop = FALSE]
ortholog_table <- ortholog_table[!duplicated(ortholog_table[["Mouse gene name"]]), , drop = FALSE]
ortholog_table <- ortholog_table[order(ortholog_table[["Gene name"]]), , drop = FALSE]

cat(sprintf("Loaded %d one-to-one ortholog pairs from %s\n\n", nrow(ortholog_table), ortholog_table_path))

###############################################################################
# STEP 3. Create one Seurat object per sample in its native feature space,
# compute sample-level QC metrics, score doublets, and save the filtered objects.
#
# Sample-level normalization, dimensional reduction, and doublet scoring are
# performed before cross-species feature restriction so each sample is processed
# from its original filtered count matrix.
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
  seurat_obj$sample <- sample_stage
  seurat_obj$sample_stage <- sample_stage
  seurat_obj$sample_key <- sample_key
  seurat_obj$legacy_input_id <- legacy_input_id
  seurat_obj$species <- species_name
  seurat_obj$species_code <- species_code

  # Compute mitochondrial and ribosomal fractions for cell-level filtering.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(MT-|mt-)")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(RPL|RPS|Rp[ls])")

  # Normalize and cluster each sample once so scDblFinder can use cluster-aware
  # doublet scoring and the sample can be inspected in a preliminary UMAP.
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = whole_heart_pcs, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = whole_heart_resolution, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)

  sce <- scDblFinder(
    GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
    clusters = Idents(seurat_obj),
    BPPARAM = bp,
    dbr = expected_doublet_rate
  )
  seurat_obj$scDblFinder.score <- sce$scDblFinder.score
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class

  # Apply species-specific cell filters to remove low-complexity libraries and
  # cells with high mitochondrial or ribosomal burden.
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
        percent.ribo < gp_max_percent_ribo
    )
  }

  seurat_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[keep_cells])

  # Recompute normalized expression, variable features, PCA, neighbors, clusters,
  # and UMAP after cell filtering so the saved object reflects the retained cells.
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = whole_heart_pcs, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = whole_heart_resolution, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)

  # Export a compact QC panel summarizing library-complexity metrics and sample
  # structure after filtering.
  sample_qc_plot <- (
    VlnPlot(
      seurat_obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
      ncol = 4,
      pt.size = 0
    ) + NoLegend()
  ) / (
    DimPlot(seurat_obj, reduction = "umap", group.by = "sample") + ggtitle(sample_stage)
  )
  ggsave(
    filename = file.path(output_dir, "qc_plots", paste0(sample_key, "_qc_panel.pdf")),
    plot = sample_qc_plot,
    width = 11,
    height = 7,
    units = "in"
  )

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
# STEP 4. Remap guinea pig genes to mouse ortholog names and restrict all
# samples to the shared one-to-one ortholog feature set used for integration.
#
# The filtered sample objects saved above are converted here into a common
# mouse-named expression space after sample-level QC and doublet annotation.
###############################################################################

sample_objects_ortholog <- list()
shared_features <- unique(ortholog_table[["Mouse gene name"]])

for (sample_id in names(sample_objects)) {
  sample_obj <- sample_objects[[sample_id]]
  counts_mat <- GetAssayData(sample_obj, assay = "RNA", slot = "counts")

  if (sample_obj$species_code[[1]] == "gp") {
    gp_keep <- ortholog_table[ortholog_table[["Gene name"]] %in% rownames(counts_mat), , drop = FALSE]
    gp_keep <- gp_keep[order(gp_keep[["Gene name"]]), , drop = FALSE]
    counts_mat <- counts_mat[gp_keep[["Gene name"]], , drop = FALSE]
    rownames(counts_mat) <- gp_keep[["Mouse gene name"]]
  } else {
    counts_mat <- counts_mat[rownames(counts_mat) %in% ortholog_table[["Mouse gene name"]], , drop = FALSE]
  }

  ortholog_obj <- CreateSeuratObject(
    counts = counts_mat,
    meta.data = sample_obj@meta.data,
    min.cells = 0,
    min.features = 0
  )
  sample_objects_ortholog[[sample_id]] <- ortholog_obj
  shared_features <- intersect(shared_features, rownames(ortholog_obj))
}

cat(sprintf("Retaining %d shared one-to-one ortholog features across all samples\n\n", length(shared_features)))

for (sample_id in names(sample_objects_ortholog)) {
  sample_objects_ortholog[[sample_id]] <- subset(sample_objects_ortholog[[sample_id]], features = shared_features)
}

sample_objects <- sample_objects_ortholog
saveRDS(sample_objects, file.path(output_dir, "example_ortholog_filtered_sample_objects.rds"))

###############################################################################
# STEP 5. Integrate the six heart datasets, cluster the shared whole-heart atlas,
# and identify the endothelial compartment for focused reanalysis.
###############################################################################

whole_heart_list <- sample_objects
for (sample_id in names(whole_heart_list)) {
  DefaultAssay(whole_heart_list[[sample_id]]) <- "RNA"
  whole_heart_list[[sample_id]] <- FindVariableFeatures(whole_heart_list[[sample_id]], verbose = FALSE)
}

whole_heart_features <- SelectIntegrationFeatures(
  object.list = whole_heart_list,
  nfeatures = whole_heart_integration_features
)

for (sample_id in names(whole_heart_list)) {
  whole_heart_list[[sample_id]] <- ScaleData(
    whole_heart_list[[sample_id]],
    features = whole_heart_features,
    verbose = FALSE
  )
  whole_heart_list[[sample_id]] <- RunPCA(
    whole_heart_list[[sample_id]],
    features = whole_heart_features,
    npcs = whole_heart_pcs,
    verbose = FALSE
  )
}

whole_heart_anchors <- FindIntegrationAnchors(
  object.list = whole_heart_list,
  anchor.features = whole_heart_features,
  dims = 1:whole_heart_pcs
)

whole_heart <- IntegrateData(anchorset = whole_heart_anchors, dims = 1:whole_heart_pcs)
DefaultAssay(whole_heart) <- "integrated"
whole_heart <- ScaleData(whole_heart, verbose = FALSE)
whole_heart <- RunPCA(whole_heart, npcs = whole_heart_pcs, verbose = FALSE)
whole_heart <- FindNeighbors(whole_heart, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)
whole_heart <- FindClusters(whole_heart, resolution = whole_heart_resolution, verbose = FALSE)
whole_heart <- RunUMAP(whole_heart, reduction = "pca", dims = 1:whole_heart_pcs, verbose = FALSE)

saveRDS(whole_heart, file.path(output_dir, "example_whole_heart_integrated.rds"))

whole_heart_plot <- (
  DimPlot(whole_heart, reduction = "umap", group.by = "sample", label = FALSE) +
    ggtitle("Whole-heart integrated atlas by sample")
) / (
  DimPlot(whole_heart, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Whole-heart integrated atlas by cluster")
)
ggsave(
  filename = file.path(output_dir, "qc_plots", "example_whole_heart_integrated_umap.pdf"),
  plot = whole_heart_plot,
  width = 8,
  height = 10,
  units = "in"
)

###############################################################################
# STEP 6. Subset endothelial cells from the whole-heart atlas and remove genes
# detected in too few endothelial cells to support the second integration round.
###############################################################################

Idents(whole_heart) <- "seurat_clusters"
endothelial <- subset(whole_heart, idents = whole_heart_endothelial_clusters)
DefaultAssay(endothelial) <- "RNA"

endothelial_counts <- GetAssayData(endothelial, assay = "RNA", slot = "counts")
genes_to_keep <- rownames(endothelial_counts)[Matrix::rowSums(endothelial_counts > 0) >= endo_gene_min_cells]
endothelial <- subset(endothelial, features = genes_to_keep)

saveRDS(endothelial, file.path(output_dir, "example_endothelial_subset_before_reintegration.rds"))

###############################################################################
# STEP 7. Reintegrate endothelial cells across samples, remove high-scoring
# doublets, and generate the refined endothelial embedding used for annotation.
###############################################################################

endothelial_list <- SplitObject(endothelial, split.by = "sample_key")
for (sample_key in names(endothelial_list)) {
  DefaultAssay(endothelial_list[[sample_key]]) <- "RNA"
  endothelial_list[[sample_key]] <- FindVariableFeatures(endothelial_list[[sample_key]], verbose = FALSE)
}

endothelial_features <- SelectIntegrationFeatures(
  object.list = endothelial_list,
  nfeatures = endo_integration_features
)

for (sample_key in names(endothelial_list)) {
  endothelial_list[[sample_key]] <- ScaleData(
    endothelial_list[[sample_key]],
    features = endothelial_features,
    verbose = FALSE
  )
  endothelial_list[[sample_key]] <- RunPCA(
    endothelial_list[[sample_key]],
    features = endothelial_features,
    npcs = endo_pcs,
    verbose = FALSE
  )
}

endothelial_anchors <- FindIntegrationAnchors(
  object.list = endothelial_list,
  anchor.features = endothelial_features,
  dims = 1:endo_pcs
)

endothelial_integrated <- IntegrateData(anchorset = endothelial_anchors, dims = 1:endo_pcs)
endothelial_integrated <- subset(endothelial_integrated, subset = scDblFinder.score < 0.2)

DefaultAssay(endothelial_integrated) <- "integrated"
endothelial_integrated <- ScaleData(endothelial_integrated, verbose = FALSE)
endothelial_integrated <- RunPCA(endothelial_integrated, npcs = endo_pcs, verbose = FALSE)
endothelial_integrated <- FindNeighbors(
  endothelial_integrated,
  reduction = "pca",
  dims = endo_neighbor_dims,
  verbose = FALSE
)
endothelial_integrated <- FindClusters(endothelial_integrated, resolution = endo_resolution, verbose = FALSE)
endothelial_integrated <- RunUMAP(
  endothelial_integrated,
  reduction = "pca",
  dims = endo_neighbor_dims,
  n.neighbors = endo_umap_neighbors,
  verbose = FALSE
)

saveRDS(endothelial_integrated, file.path(output_dir, "example_endothelial_integrated_before_cell_cycle.rds"))

###############################################################################
# STEP 8. Score cell-cycle state in the refined endothelial atlas.
#
# The canonical Seurat S-phase and G2/M-phase gene sets are converted from human
# to mouse ortholog names so cell-cycle state can be scored in the mouse-named
# expression matrix.
###############################################################################

DefaultAssay(endothelial_integrated) <- "RNA"

s_genes <- unique(gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name)

g2m_genes <- unique(gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name)

endothelial_integrated <- CellCycleScoring(
  endothelial_integrated,
  s.features = intersect(s_genes, rownames(endothelial_integrated)),
  g2m.features = intersect(g2m_genes, rownames(endothelial_integrated)),
  set.ident = FALSE
)

# Keep only scDblFinder singlets in the object that will receive final subtype labels.
endothelial_integrated <- subset(endothelial_integrated, subset = scDblFinder.class == "singlet")
saveRDS(endothelial_integrated, file.path(output_dir, "example_endothelial_integrated_prelabel.rds"))

###############################################################################
# STEP 9. Assign fine and coarse endothelial subtype labels from the final
# endothelial clusters and save the labeled object and metadata table.
###############################################################################

endothelial_integrated$cell.type.l1 <- unname(fine_label_map[as.character(endothelial_integrated$seurat_clusters)])
endothelial_integrated$cell.type.l2 <- unname(coarse_label_map[as.character(endothelial_integrated$seurat_clusters)])

if (any(is.na(endothelial_integrated$cell.type.l1))) {
  warning("Some final clusters were not assigned a fine label. Check `fine_label_map`.")
}
if (any(is.na(endothelial_integrated$cell.type.l2))) {
  warning("Some final clusters were not assigned a coarse label. Check `coarse_label_map`.")
}

endothelial_integrated$cell.type.l1 <- factor(
  endothelial_integrated$cell.type.l1,
  levels = c("Artery", "Arteriole", "Cap4", "Cap3", "Cap2", "Cap1", "Vein")
)
endothelial_integrated$cell.type.l2 <- factor(
  endothelial_integrated$cell.type.l2,
  levels = c("Artery", "Arteriole", "Capillary", "Vein")
)

saveRDS(endothelial_integrated, file.path(output_dir, "example_endothelial_integrated_labeled.rds"))

endothelial_metadata <- endothelial_integrated@meta.data
endothelial_metadata$cell <- rownames(endothelial_metadata)
write.table(
  endothelial_metadata,
  file = file.path(output_dir, "example_endothelial_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 10. Export summary plots for the final endothelial atlas.
###############################################################################

final_endothelial_plot <- (
  DimPlot(endothelial_integrated, reduction = "umap", group.by = "sample", label = FALSE) +
    ggtitle("Endothelial atlas by sample")
) / (
  DimPlot(endothelial_integrated, reduction = "umap", group.by = "cell.type.l1", label = TRUE, repel = TRUE) +
    ggtitle("Fine endothelial labels")
) / (
  DimPlot(endothelial_integrated, reduction = "umap", group.by = "cell.type.l2", label = TRUE, repel = TRUE) +
    ggtitle("Coarse endothelial labels")
)

ggsave(
  filename = file.path(output_dir, "qc_plots", "example_endothelial_final_umap_panels.pdf"),
  plot = final_endothelial_plot,
  width = 8,
  height = 14,
  units = "in"
)

feature_plot <- FeaturePlot(
  endothelial_integrated,
  features = fate_markers,
  ncol = 3,
  order = TRUE,
  raster = FALSE
)

ggsave(
  filename = file.path(output_dir, "qc_plots", "example_endothelial_fate_marker_featureplots.pdf"),
  plot = feature_plot,
  width = 10,
  height = 10,
  units = "in"
)

cat("Interactive heart EC reproduction script template completed.\n")
cat(sprintf("Configured example output directory: %s\n", output_dir))
