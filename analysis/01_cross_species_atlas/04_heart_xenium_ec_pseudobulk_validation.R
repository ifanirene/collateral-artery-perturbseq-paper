#!/usr/bin/env Rscript

# @theme: T2 - Xenium spatial endothelial pseudobulk
# @task_id: T2-01
# @description: Minimal step-by-step script for reproducing the maintained Xenium endothelial pseudobulk workflow used for manuscript support, including region discovery, Seurat v5 integration, endothelial masking, and section-cluster pseudobulk DE.
# @inputs:
# - `path/to/GP_arc_reanalysis/spatial_transcriptome/Xenium_Runs/output-XETG00277__*`: four vendor Xenium region-output folders used in the maintained reanalysis.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/cluster_to_celltype.csv`: optional manual cluster-to-celltype map for endothelial masking.
# @outputs:
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/xenium_section_manifest.tsv`: region-level manifest with run metadata.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/01_merged_qc.rds`: merged Xenium object after QC and centroid-based section-cluster assignment.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/02_integrated_sketch.rds`: Seurat v5 sketch/CCA-integrated Xenium object.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/03_ec_subset.rds`: endothelial-only Xenium object used for pseudobulk quantification.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/outputs/de/ec_pseudobulk_counts.tsv`: raw endothelial pseudobulk matrix.
# - `path/to/GP_arc_reanalysis/Docs/xenium_public_repro/outputs/de/ec_deseq2_results.tsv`: DESeq2 results for Group2 versus Group1.
# @key_params:
# - `qc_min_nCount = 100`, `qc_max_nCount = 2000`, `qc_min_nFeature = 50`: Xenium cell-level QC thresholds copied from the maintained branch.
# - `target_clusters_per_region = 4`: number of centroid-derived spatial replicate units per Xenium region output in this minimal reproduction script.
# - `sketch_ncells_per_region = 50000`: number of cells per region layer retained for Seurat v5 sketch integration.
# - `n_pcs = 30`, `cluster_resolution = 0.4`: dimensional-reduction and clustering settings for the integrated sketch.
# - `ec_fallback_min_markers_detected = 2`: endothelial marker threshold applied when no manual cluster labels are supplied.
# @dependencies:
# - `Seurat` (v5): Xenium loading, preprocessing, sketch integration, dimensional reduction, clustering, and visualization.
# - `SeuratObject`: Seurat object infrastructure and layer handling.
# - `DESeq2`: pseudobulk differential expression.
# - `BiocParallel`: parallel DESeq2 execution.
# - `Matrix`: sparse pseudobulk aggregation.
# - `future`: Seurat parallelization.
# - `jsonlite`: parsing Xenium `experiment.xenium` metadata.
# - `data.table`: importing `cells.csv.gz` metadata.
# - `ggplot2`: PCA and sample-distance diagnostics.
# @examples:
# Rscript path/to/reproduce_xenium_ec_pseudobulk_pipeline.R
#
# Run this script block by block in the Seurat v5 environment used for the
# Xenium branch. Replace the example paths below with the locations used in your
# own manuscript-support package. This script is intentionally compact and
# mirrors the maintained all-sections analysis at a reader-facing level rather
# than reproducing every internal checkpoint file in `xenium-ec-pseudobulk/`.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DESeq2)
  library(BiocParallel)
  library(Matrix)
  library(future)
  library(jsonlite)
  library(data.table)
  library(ggplot2)
})

cat("Live session information for this reproduction script:\n")
print(sessionInfo())
cat("\n")

###############################################################################
# STEP 1. Define the Xenium region outputs, output directory, and analysis
# parameters.
#
# The maintained branch processed four vendor Xenium region outputs and then
# defined spatial replicate units from cell-centroid structure within each
# region. This minimal script mirrors that branch but keeps all inputs and
# outputs in one reader-facing location.
###############################################################################

project_root <- "path/to/GP_arc_reanalysis"
output_dir <- file.path(project_root, "Docs", "xenium_public_repro")

region_dirs <- c(
  file.path(project_root, "spatial_transcriptome", "Xenium_Runs", "output-XETG00277__0038141__Region_1__20241113__194049"),
  file.path(project_root, "spatial_transcriptome", "Xenium_Runs", "output-XETG00277__0038141__Region_2__20241113__194049"),
  file.path(project_root, "spatial_transcriptome", "Xenium_Runs", "output-XETG00277__0038159__Region_1__20241113__194049"),
  file.path(project_root, "spatial_transcriptome", "Xenium_Runs", "output-XETG00277__0038159__Region_2__20241113__194049")
)

sample_id_to_group <- c(
  "0038141" = "Group1",
  "0038159" = "Group2"
)

qc_min_nCount <- 100
qc_max_nCount <- 2000
qc_min_nFeature <- 50

# This minimal script mirrors the maintained all-sections run by producing four
# spatial replicate units per region output from cell-centroid coordinates.
target_clusters_per_region <- 4

sketch_ncells_per_region <- 50000
n_pcs <- 30
cluster_resolution <- 0.4

ec_celltype_patterns <- c("endothelial", "^ec$")
ec_fallback_markers <- c("Flt1", "Kdr", "Tek", "Egfl7", "Klf2", "Klf4", "Dach1", "Notch1")
ec_fallback_min_markers_detected <- 2
min_ec_cells_warn <- 500

workers_seurat <- 4
workers_deseq <- 4

manifest_tsv <- file.path(output_dir, "xenium_section_manifest.tsv")
cluster_map_csv <- file.path(output_dir, "cluster_to_celltype.csv")

merged_rds <- file.path(output_dir, "01_merged_qc.rds")
integrated_rds <- file.path(output_dir, "02_integrated_sketch.rds")
ec_rds <- file.path(output_dir, "03_ec_subset.rds")

qc_dir <- file.path(output_dir, "outputs", "qc")
umap_dir <- file.path(output_dir, "outputs", "umap")
marker_dir <- file.path(output_dir, "outputs", "markers")
de_dir <- file.path(output_dir, "outputs", "de")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

missing_region_dirs <- region_dirs[!file.exists(region_dirs)]
if (length(missing_region_dirs) > 0) {
  stop(
    sprintf("Missing Xenium region directories:\n- %s", paste(missing_region_dirs, collapse = "\n- ")),
    call. = FALSE
  )
}

###############################################################################
# STEP 2. Discover one Xenium output per region folder and build a region-level
# manifest carrying the identifiers used downstream.
#
# The maintained branch looked for `experiment.xenium` plus `cell_feature_matrix`
# to define what counts as one vendor Xenium output. The same logic is used
# here, but the resulting manifest is written directly to the release directory.
###############################################################################

find_xenium_output_dirs <- function(region_dir) {
  experiment_paths <- list.files(region_dir, pattern = "experiment\\.xenium$", recursive = TRUE, full.names = TRUE)
  if (length(experiment_paths) == 0) {
    stop(sprintf("No experiment.xenium file found under %s", region_dir), call. = FALSE)
  }

  output_dirs <- unique(dirname(experiment_paths))
  output_dirs <- output_dirs[file.exists(file.path(output_dirs, "cell_feature_matrix"))]
  if (length(output_dirs) == 0) {
    stop(sprintf("No Xenium output directory with cell_feature_matrix found under %s", region_dir), call. = FALSE)
  }
  normalizePath(sort(output_dirs), winslash = "/", mustWork = TRUE)
}

maybe_use_outs_dir <- function(xenium_dir) {
  outs_dir <- file.path(xenium_dir, "outs")
  if (dir.exists(outs_dir)) outs_dir else xenium_dir
}

section_manifest <- list()

for (region_dir in region_dirs) {
  xenium_output_dirs <- find_xenium_output_dirs(region_dir)

  for (i in seq_along(xenium_output_dirs)) {
    xenium_dir <- xenium_output_dirs[[i]]
    experiment_path <- file.path(xenium_dir, "experiment.xenium")
    experiment_meta <- jsonlite::fromJSON(experiment_path)

    sample_id <- as.character(experiment_meta$slide_id)
    region_name <- as.character(experiment_meta$region_name)
    section_id <- sprintf("%s_%s_sec%01d", sample_id, region_name, i)

    section_manifest[[length(section_manifest) + 1]] <- data.frame(
      section_id = section_id,
      group = unname(sample_id_to_group[[sample_id]]),
      sample_id = sample_id,
      region = region_name,
      xenium_out_dir = normalizePath(maybe_use_outs_dir(xenium_dir), winslash = "/", mustWork = TRUE),
      analysis_sw_version = as.character(experiment_meta$analysis_sw_version),
      analysis_uuid = as.character(experiment_meta$analysis_uuid),
      experiment_uuid = as.character(experiment_meta$experiment_uuid),
      roi_uuid = as.character(experiment_meta$roi_uuid),
      panel_num_targets_custom = as.integer(experiment_meta$panel_num_targets_custom),
      stringsAsFactors = FALSE
    )
  }
}

section_manifest <- do.call(rbind, section_manifest)
write.table(section_manifest, file = manifest_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote Xenium section manifest: %s\n\n", manifest_tsv))

###############################################################################
# STEP 3. Load the Xenium region outputs into Seurat v5, append cell-level
# metrics from `cells.csv.gz`, and apply the maintained QC thresholds.
#
# The maintained branch disabled molecule-coordinate loading, retained the
# Xenium gene-expression assay, and appended per-cell centroid and control-count
# metadata before filtering to 100-2,000 Xenium transcripts and at least 50
# detected genes.
###############################################################################

safe_future_plan <- function(workers) {
  available <- tryCatch(parallelly::availableCores(), error = function(e) parallel::detectCores())
  workers_use <- min(workers, available)
  options(future.globals.maxSize = 30 * 1024^3)
  future::plan(future::multisession, workers = workers_use)
}

safe_future_plan(workers_seurat)

read_cells_metadata <- function(cells_csv_gz) {
  if (!file.exists(cells_csv_gz)) {
    stop(sprintf("Missing cells.csv.gz file: %s", cells_csv_gz), call. = FALSE)
  }

  data.table::fread(
    cmd = paste("zcat", shQuote(cells_csv_gz)),
    select = c(
      "cell_id",
      "x_centroid",
      "y_centroid",
      "transcript_counts",
      "control_probe_counts",
      "genomic_control_counts",
      "control_codeword_counts",
      "unassigned_codeword_counts",
      "deprecated_codeword_counts",
      "total_counts",
      "cell_area",
      "nucleus_area",
      "nucleus_count",
      "segmentation_method"
    ),
    showProgress = FALSE
  )
}

xenium_objects <- list()

for (i in seq_len(nrow(section_manifest))) {
  section_row <- section_manifest[i, , drop = FALSE]
  cat(sprintf("Loading %s from %s\n", section_row$section_id, section_row$xenium_out_dir))

  xenium_obj <- LoadXenium(
    data.dir = section_row$xenium_out_dir,
    assay = "Xenium",
    molecule.coordinates = FALSE
  )

  cells_csv_gz <- file.path(dirname(section_row$xenium_out_dir), "cells.csv.gz")
  if (!file.exists(cells_csv_gz)) {
    cells_csv_gz <- file.path(section_row$xenium_out_dir, "cells.csv.gz")
  }
  cell_meta <- as.data.frame(read_cells_metadata(cells_csv_gz), stringsAsFactors = FALSE)
  rownames(cell_meta) <- cell_meta$cell_id
  cell_meta$cell_id <- NULL

  denom <- pmax(cell_meta$total_counts, 1)
  cell_meta$control_probe_frac <- cell_meta$control_probe_counts / denom
  cell_meta$unassigned_codeword_frac <- cell_meta$unassigned_codeword_counts / denom

  xenium_obj$section_id <- section_row$section_id
  xenium_obj$group <- section_row$group
  xenium_obj$sample_id <- section_row$sample_id
  xenium_obj$region <- section_row$region
  xenium_obj$analysis_sw_version <- section_row$analysis_sw_version

  xenium_obj <- AddMetaData(xenium_obj, metadata = cell_meta)
  xenium_objects[[section_row$section_id]] <- xenium_obj
}

merged_xenium <- merge(
  x = xenium_objects[[1]],
  y = xenium_objects[-1],
  add.cell.ids = names(xenium_objects)
)

DefaultAssay(merged_xenium) <- "Xenium"
merged_xenium <- subset(
  merged_xenium,
  subset =
    nCount_Xenium >= qc_min_nCount &
    nCount_Xenium <= qc_max_nCount &
    nFeature_Xenium >= qc_min_nFeature
)

cat(sprintf("Cells retained after Xenium QC: %d\n\n", ncol(merged_xenium)))

###############################################################################
# STEP 4. Define spatial replicate units (`section_cluster`) within each region
# output using x-y cell centroids.
#
# The maintained branch stores downstream pseudobulk replicates as spatial
# section clusters rather than as whole-region outputs. This minimal script uses
# centroid-based k-means with four clusters per region output so the resulting
# replicate count matches the maintained all-sections run.
###############################################################################

meta <- merged_xenium@meta.data
meta$section_cluster_in_region <- NA_integer_
meta$section_cluster <- NA_character_

for (section_id in sort(unique(meta$section_id))) {
  section_idx <- which(meta$section_id == section_id)
  section_meta <- meta[section_idx, , drop = FALSE]
  coords <- as.matrix(section_meta[, c("x_centroid", "y_centroid")])

  km <- stats::kmeans(coords, centers = target_clusters_per_region, iter.max = 50)
  cluster_sizes <- sort(table(km$cluster), decreasing = TRUE)
  cluster_map <- setNames(seq_along(cluster_sizes), names(cluster_sizes))
  section_cluster_in_region <- as.integer(cluster_map[as.character(km$cluster)])

  meta$section_cluster_in_region[section_idx] <- section_cluster_in_region
  meta$section_cluster[section_idx] <- sprintf(
    "%s_sc%02d",
    section_id,
    section_cluster_in_region
  )
}

merged_xenium@meta.data <- meta
saveRDS(merged_xenium, file = merged_rds)

section_cluster_summary <- as.data.frame(table(merged_xenium$section_cluster), stringsAsFactors = FALSE)
colnames(section_cluster_summary) <- c("section_cluster", "n_cells")
write.table(
  section_cluster_summary,
  file = file.path(qc_dir, "section_cluster_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 5. Integrate the Xenium region outputs in Seurat v5 using sketch-based
# CCA integration and generate a blank cluster-to-celltype template.
#
# The maintained branch split the Xenium assay by region output, sketch-sampled
# each layer, integrated the sketch with CCA, then projected the integrated
# embedding and cluster labels back to all cells.
###############################################################################

gene_features <- rownames(merged_xenium[["Xenium"]])

merged_xenium[["Xenium"]] <- JoinLayers(merged_xenium[["Xenium"]], layers = "counts")
merged_xenium[["Xenium"]] <- split(merged_xenium[["Xenium"]], f = merged_xenium$section_id, layers = "counts")

merged_xenium <- NormalizeData(merged_xenium, assay = "Xenium", normalization.method = "LogNormalize", verbose = FALSE)
merged_xenium <- FindVariableFeatures(
  merged_xenium,
  assay = "Xenium",
  features = gene_features,
  nfeatures = min(200, length(gene_features))
)

sketch_method <- if (length(gene_features) < 50) "Uniform" else "LeverageScore"

merged_xenium <- SketchData(
  object = merged_xenium,
  assay = "Xenium",
  ncells = sketch_ncells_per_region,
  sketched.assay = "sketch",
  method = sketch_method,
  features = gene_features,
  over.write = TRUE,
  seed = 123,
  verbose = FALSE
)

DefaultAssay(merged_xenium) <- "sketch"
merged_xenium <- ScaleData(merged_xenium, assay = "sketch", features = gene_features, verbose = FALSE)
merged_xenium <- RunPCA(merged_xenium, assay = "sketch", features = gene_features, npcs = n_pcs, verbose = FALSE)

sketch_layers <- Layers(merged_xenium[["sketch"]], search = "data")
layer_cell_counts <- vapply(sketch_layers, function(layer_name) length(Cells(merged_xenium[["sketch"]], layer = layer_name)), integer(1))
k_weight_use <- min(100, max(10, min(layer_cell_counts) - 1))

merged_xenium <- IntegrateLayers(
  object = merged_xenium,
  method = CCAIntegration,
  assay = "sketch",
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:n_pcs,
  layers = sketch_layers,
  k.weight = k_weight_use,
  verbose = FALSE
)

merged_xenium <- RunUMAP(merged_xenium, reduction = "integrated.cca", dims = 1:n_pcs, return.model = TRUE, verbose = FALSE)
merged_xenium <- FindNeighbors(merged_xenium, reduction = "integrated.cca", dims = 1:n_pcs, verbose = FALSE)
merged_xenium <- FindClusters(merged_xenium, resolution = cluster_resolution, verbose = FALSE)

merged_xenium <- ProjectData(
  object = merged_xenium,
  assay = "Xenium",
  sketched.assay = "sketch",
  sketched.reduction = "integrated.cca",
  full.reduction = "integrated.cca",
  dims = 1:n_pcs,
  refdata = list(seurat_clusters = "seurat_clusters"),
  umap.model = "umap",
  recompute.neighbors = FALSE,
  recompute.weights = FALSE,
  verbose = FALSE
)

cluster_template <- data.frame(
  cluster = sort(unique(as.character(merged_xenium$seurat_clusters))),
  celltype = "",
  stringsAsFactors = FALSE
)
write.csv(cluster_template, file = cluster_map_csv, row.names = FALSE, quote = TRUE)

saveRDS(merged_xenium, file = integrated_rds)

###############################################################################
# STEP 6. Apply manual cell-type labels when provided, otherwise use the
# maintained branch's endothelial marker fallback mask.
#
# The saved all-sections run used the marker-based fallback mask because the
# cluster-to-celltype table remained empty. This script preserves the same
# default while still allowing manual labels to be supplied.
###############################################################################

apply_cluster_mapping <- function(obj, mapping_df) {
  map <- setNames(as.character(mapping_df$celltype), as.character(mapping_df$cluster))
  celltype <- unname(map[as.character(obj$seurat_clusters)])
  celltype[is.na(celltype) | celltype == ""] <- "Unassigned"
  obj$celltype <- celltype
  obj
}

get_counts_matrix <- function(obj, assay = "Xenium") {
  tryCatch({
    GetAssayData(obj, assay = assay, layer = "counts")
  }, error = function(e) {
    obj[[assay]] <- JoinLayers(obj[[assay]], layers = "counts")
    GetAssayData(obj, assay = assay, layer = "counts")
  })
}

select_endothelial_cells <- function(obj) {
  if (file.exists(cluster_map_csv)) {
    cluster_map <- read.csv(cluster_map_csv, stringsAsFactors = FALSE)
    if (all(c("cluster", "celltype") %in% colnames(cluster_map))) {
      obj <- apply_cluster_mapping(obj, cluster_map)

      keep <- rep(FALSE, ncol(obj))
      for (pattern in ec_celltype_patterns) {
        keep <- keep | grepl(pattern, obj$celltype, ignore.case = TRUE)
      }
      if (sum(keep) > 0) {
        return(list(object = obj, cells = colnames(obj)[keep], method = "manual_cluster_mapping"))
      }
    }
  }

  counts <- get_counts_matrix(obj, assay = "Xenium")
  markers_present <- intersect(ec_fallback_markers, rownames(counts))
  if (length(markers_present) == 0) {
    stop("No endothelial fallback markers were found in the Xenium assay.", call. = FALSE)
  }

  marker_detected <- Matrix::colSums(counts[markers_present, , drop = FALSE] > 0)
  keep <- marker_detected >= ec_fallback_min_markers_detected

  obj$celltype <- "Unassigned"
  obj$celltype[keep] <- "Endothelial"

  list(object = obj, cells = colnames(obj)[keep], method = "marker_fallback")
}

merged_xenium <- readRDS(integrated_rds)
endo_selection <- select_endothelial_cells(merged_xenium)
merged_xenium <- endo_selection$object

cat(sprintf("Endothelial-cell selection method: %s\n", endo_selection$method))
cat(sprintf("Endothelial cells retained: %d\n\n", length(endo_selection$cells)))

ec_xenium <- subset(merged_xenium, cells = endo_selection$cells)
saveRDS(ec_xenium, file = ec_rds)

###############################################################################
# STEP 7. Aggregate raw Xenium counts within spatial section-cluster replicates
# and perform DESeq2 differential testing for Group2 versus Group1.
#
# The maintained branch summed raw counts across endothelial cells within each
# `section_cluster`, exported control features separately for QC, and tested the
# group effect with DESeq2 using the design `~ region + group`.
###############################################################################

make_pseudobulk <- function(counts, replicate_ids) {
  replicate_factor <- factor(replicate_ids)
  design_matrix <- Matrix::sparse.model.matrix(~ 0 + replicate_factor)
  colnames(design_matrix) <- levels(replicate_factor)
  counts %*% design_matrix
}

counts_ec <- get_counts_matrix(ec_xenium, assay = "Xenium")
pb_gene <- as.matrix(make_pseudobulk(counts_ec, ec_xenium$section_cluster))

pb_gene_tsv <- file.path(de_dir, "ec_pseudobulk_counts.tsv")
write.table(
  data.frame(gene = rownames(pb_gene), pb_gene, check.names = FALSE),
  file = pb_gene_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

control_assays <- intersect(SeuratObject::Assays(ec_xenium), c("ControlProbe", "ControlCodeword", "BlankCodeword"))
if (length(control_assays) > 0) {
  control_pseudobulk <- list()
  for (control_assay in control_assays) {
    counts_control <- get_counts_matrix(ec_xenium, assay = control_assay)
    pb_control <- as.matrix(make_pseudobulk(counts_control, ec_xenium$section_cluster))
    rownames(pb_control) <- paste0(control_assay, "::", rownames(pb_control))
    control_pseudobulk[[control_assay]] <- pb_control
  }
  control_pseudobulk <- do.call(rbind, control_pseudobulk)
  write.table(
    data.frame(gene = rownames(control_pseudobulk), control_pseudobulk, check.names = FALSE),
    file = file.path(de_dir, "ec_pseudobulk_counts_controls.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

ec_counts_by_section_cluster <- as.data.frame(table(ec_xenium$section_cluster), stringsAsFactors = FALSE)
colnames(ec_counts_by_section_cluster) <- c("section_cluster", "n_ec_cells")
write.table(
  ec_counts_by_section_cluster,
  file = file.path(de_dir, "ec_cell_counts_by_section_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ec_counts_by_section <- as.data.frame(table(ec_xenium$section_id), stringsAsFactors = FALSE)
colnames(ec_counts_by_section) <- c("section_id", "n_ec_cells")
write.table(
  ec_counts_by_section,
  file = file.path(de_dir, "ec_cell_counts_by_section_id.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

low_ec_replicates <- ec_counts_by_section_cluster$section_cluster[ec_counts_by_section_cluster$n_ec_cells < min_ec_cells_warn]
if (length(low_ec_replicates) > 0) {
  warning(
    sprintf(
      "The following section clusters contain fewer than %d endothelial cells: %s",
      min_ec_cells_warn,
      paste(low_ec_replicates, collapse = ", ")
    ),
    call. = FALSE
  )
}

coldata <- unique(ec_xenium@meta.data[, c("section_cluster", "section_id", "group", "region", "sample_id")])
coldata <- coldata[match(colnames(pb_gene), coldata$section_cluster), , drop = FALSE]
rownames(coldata) <- coldata$section_cluster
coldata$group <- factor(coldata$group, levels = c("Group1", "Group2"))
coldata$region <- factor(coldata$region)
coldata$sample_id <- factor(coldata$sample_id)

bp <- MulticoreParam(min(workers_deseq, tryCatch(parallelly::availableCores(), error = function(e) parallel::detectCores())))
dds <- DESeqDataSetFromMatrix(
  countData = round(pb_gene),
  colData = coldata,
  design = ~ region + group
)
dds$group <- relevel(dds$group, ref = "Group1")
dds <- DESeq(dds, parallel = TRUE, BPPARAM = bp)

deseq_results_obj <- results(dds, contrast = c("group", "Group2", "Group1"))
deseq_results <- as.data.frame(deseq_results_obj)
deseq_results$gene <- rownames(deseq_results)
deseq_results <- deseq_results[, c("gene", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
deseq_results <- deseq_results[order(deseq_results$padj, na.last = TRUE), ]

results_tsv <- file.path(de_dir, "ec_deseq2_results.tsv")
write.table(
  deseq_results,
  file = results_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = c("group", "region"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

pdf(file.path(de_dir, "ec_pseudobulk_pca.pdf"), width = 7, height = 5, useDingbats = FALSE)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = region)) +
    geom_point(size = 3) +
    theme_classic() +
    xlab(sprintf("PC1: %s%% variance", percent_var[1])) +
    ylab(sprintf("PC2: %s%% variance", percent_var[2])) +
    ggtitle("Xenium endothelial pseudobulk PCA")
)
dev.off()

sample_distance <- as.matrix(dist(t(assay(vsd))))
sample_distance_df <- as.data.frame(as.table(sample_distance), stringsAsFactors = FALSE)
colnames(sample_distance_df) <- c("replicate_i", "replicate_j", "distance")

pdf(file.path(de_dir, "ec_sample_distance_heatmap.pdf"), width = 7, height = 6, useDingbats = FALSE)
print(
  ggplot(sample_distance_df, aes(x = replicate_i, y = replicate_j, fill = distance)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Xenium endothelial pseudobulk sample distances", x = "", y = "")
)
dev.off()

pdf(file.path(de_dir, "ec_ma_plot.pdf"), width = 7, height = 5, useDingbats = FALSE)
plotMA(deseq_results_obj, ylim = c(-3, 3))
dev.off()

volcano_df <- deseq_results
volcano_df$neglog10padj <- -log10(volcano_df$padj)
volcano_df$neglog10padj[is.infinite(volcano_df$neglog10padj)] <- NA_real_

pdf(file.path(de_dir, "ec_volcano.pdf"), width = 7, height = 6, useDingbats = FALSE)
print(
  ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_classic() +
    labs(title = "Xenium endothelial pseudobulk: Group2 versus Group1", x = "log2 fold change", y = "-log10 adjusted P")
)
dev.off()

cat("Minimal Xenium reproduction script completed.\n")
cat(sprintf("Integrated object: %s\n", integrated_rds))
cat(sprintf("Endothelial subset: %s\n", ec_rds))
cat(sprintf("Pseudobulk matrix: %s\n", pb_gene_tsv))
cat(sprintf("DESeq2 results: %s\n", results_tsv))
