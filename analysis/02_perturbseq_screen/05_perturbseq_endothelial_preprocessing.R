#!/usr/bin/env Rscript

# @description
# Publication-supporting Perturb-seq preprocessing workflow for the neonatal
# brain endothelial screen.
# This script reconstructs the manuscript-facing preprocessing path from the
# legacy Final Pool analysis: it combines resequenced RNA counts with the
# original guide-capture counts, applies RNA-level QC, performs manual
# endothelial-focused cluster pruning, removes doublets, applies MOI filtering,
# and exports matrices for downstream SCEPTRE and cNMF analyses.
#
# Key features:
# - Loads deeper gene-expression matrices together with original raw gRNA counts
# - Reproduces the manual cluster-pruning steps used to focus the analysis on
#   endothelial and closely related screen populations
# - Applies scDblFinder singlet filtering and MOI-based guide-complexity
#   filtering before downstream export
# - Writes manuscript-facing intermediate objects, QC plots, session info, and
#   SCEPTRE/cNMF-ready outputs to a dedicated Docs subdirectory
#
# @dependencies
# - Seurat / SeuratDisk: single-cell preprocessing, annotation, and export
# - Matrix / dplyr: sparse-matrix handling and metadata tables
# - scDblFinder / SingleCellExperiment / BiocParallel: doublet detection
# - gprofiler2: ortholog mapping for cell-cycle scoring
# - ggplot2 / patchwork / viridis: QC visualizations
#
# @examples
# VALIDATE_ONLY=1 conda run -p /oak/stanford/groups/kredhors/irenefan/seurat4 \
#   Rscript Docs/reproduce_perturb_seq_preprocessing_pipeline.R
# conda run -p /oak/stanford/groups/kredhors/irenefan/seurat4 \
#   Rscript Docs/reproduce_perturb_seq_preprocessing_pipeline.R

options(stringsAsFactors = FALSE)
options(Seurat.parallel = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(gprofiler2)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(getwd())
  }
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "perturb_seq_public_repro")
plot_dir <- file.path(output_dir, "qc_plots")
validate_only <- identical(Sys.getenv("VALIDATE_ONLY"), "1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# STEP 0. Record the environment snapshot and the live session info.
#
# The snapshot below matches the documented seurat4 environment used for the
# manuscript-facing Docs scripts in this project. The live session info is also
# written so any rerun can be traced to the exact runtime environment.
###############################################################################

cat("============================================================\n")
cat("Neonatal brain Perturb-seq preprocessing reproduction script\n")
cat("Recorded seurat4 environment snapshot (2026-03-29)\n")
cat("============================================================\n")
cat("R.version: R version 4.2.0 (2022-04-22)\n")
cat("platform: x86_64-conda-linux-gnu (64-bit)\n")
cat("running: CentOS Linux 7 (Core)\n")
cat("Seurat: 4.4.0\n")
cat("SeuratObject: 5.0.1\n")
cat("scDblFinder: 1.17.2\n")
cat("BiocParallel: 1.32.5\n")
cat("gprofiler2: 0.2.3\n")
cat("ggplot2: 3.5.0\n")
cat("patchwork: 1.2.0\n")
cat("------------------------------------------------------------\n")
cat("Live sessionInfo() from the current run\n")
cat("------------------------------------------------------------\n")
print(sessionInfo())
cat("============================================================\n\n")

###############################################################################
# STEP 1. Define the input locations and fixed preprocessing parameters.
#
# The script pairs each resequenced `*_gex` directory with the corresponding
# original Cell Ranger output so that the final object uses deeper RNA counts
# while retaining the original CRISPR guide-capture matrix.
###############################################################################

base_dir <- "/home/irenefan/oak/perturb-seq/Final_pool/CellRanger_output"
gRNA_feature_name <- "CRISPR Guide Capture"

qc_min_features <- 2000
qc_max_features <- 7500
qc_min_counts <- 3000
qc_max_counts <- 30000
qc_max_percent_mt <- 10

initial_pcs <- 30
initial_resolution <- 0.4
recluster_pcs <- 10
recluster_resolution <- 0.4

workers <- 4L
doublet_rate <- 0.10
moi_umi_threshold <- 5
moi_max <- 15
sceptre_variable_genes <- 15000
min_cells_per_export_gene <- 100
h5_export <- TRUE

clusters_to_drop_round1 <- c("7", "8", "9", "10", "11", "13", "14", "15", "17", "18")
clusters_to_drop_round2 <- c("8")
clusters_to_drop_round3 <- c("11", "13")
clusters_to_keep_final <- c("0", "1", "2", "3", "4", "5", "6", "7", "10", "12", "16", "20")

cluster_map_broad <- c(
  `0` = "capillary",
  `1` = "capillary",
  `2` = "pre-artery",
  `3` = "capillary",
  `4` = "capillary",
  `5` = "vein",
  `6` = "artery",
  `7` = "cycling",
  `8` = "pericyte",
  `9` = "erythrocyte",
  `10` = "capillary",
  `11` = "pericyte",
  `12` = "cycling",
  `13` = "microglia",
  `14` = "microglia",
  `15` = "neuron",
  `16` = "artery",
  `17` = "astrocyte",
  `18` = "astrocyte",
  `19` = "neuron",
  `20` = "capillary"
)

cluster_map_fine <- c(
  `0` = "cap_1",
  `1` = "cap_2",
  `2` = "pre-artery",
  `3` = "cap_3",
  `4` = "cap_4",
  `5` = "vein",
  `6` = "artery_1",
  `7` = "cycling",
  `8` = "pericyte_1",
  `9` = "erythrocyte",
  `10` = "cap_5",
  `11` = "pericyte_2",
  `12` = "cycling",
  `13` = "microglia_1",
  `14` = "microglia_2",
  `15` = "neuron_dev",
  `16` = "artery_2",
  `17` = "astrocyte_1",
  `18` = "astrocyte_2",
  `19` = "neuron_mature",
  `20` = "cap_5"
)

gex_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
gex_dirs <- gex_dirs[grepl("_gex$", basename(gex_dirs))]
gex_dirs <- sort(gex_dirs)

if (length(gex_dirs) == 0) {
  stop("No *_gex directories were found in the Cell Ranger output directory.", call. = FALSE)
}

sample_sheet <- data.frame(
  sample_key = gsub("_gex$", "", basename(gex_dirs)),
  gex_dir = gex_dirs,
  original_run_dir = file.path(base_dir, gsub("_gex$", "", basename(gex_dirs))),
  stringsAsFactors = FALSE
)

missing_filtered <- sample_sheet$gex_dir[
  !file.exists(file.path(sample_sheet$gex_dir, "filtered_feature_bc_matrix.h5"))
]
missing_raw <- sample_sheet$original_run_dir[
  !file.exists(file.path(sample_sheet$original_run_dir, "outs", "raw_feature_bc_matrix.h5"))
]

if (length(missing_filtered) > 0 || length(missing_raw) > 0) {
  stop(
    paste(
      c(
        if (length(missing_filtered) > 0) {
          sprintf("Missing filtered_feature_bc_matrix.h5 for:\n- %s", paste(missing_filtered, collapse = "\n- "))
        },
        if (length(missing_raw) > 0) {
          sprintf("Missing raw_feature_bc_matrix.h5 for:\n- %s", paste(missing_raw, collapse = "\n- "))
        }
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

write.table(
  sample_sheet,
  file = file.path(output_dir, "sample_sheet_used.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("Discovered %d resequenced batches.\n", nrow(sample_sheet)))
print(sample_sheet)

if (validate_only) {
  cat("VALIDATE_ONLY=1, so the script stops after configuration and environment checks.\n")
  quit(save = "no", status = 0)
}

###############################################################################
# STEP 2. Load one Seurat object per batch using deeper RNA counts and the
# original raw guide-capture counts.
#
# Each batch keeps the resequenced filtered RNA matrix together with the raw
# `CRISPR Guide Capture` counts from the matching original run. Barcodes are
# intersected so both assays refer to the same cells before merging.
###############################################################################

load_batch <- function(gex_dir, original_run_dir, sample_key, guide_feature_name) {
  cat(sprintf("[LOAD] %s\n", sample_key))

  gex_data <- Read10X_h5(file.path(gex_dir, "filtered_feature_bc_matrix.h5"))
  if (inherits(gex_data, "list")) {
    if ("Gene Expression" %in% names(gex_data)) {
      rna_counts <- gex_data[["Gene Expression"]]
    } else {
      rna_counts <- gex_data[[1]]
    }
  } else {
    rna_counts <- gex_data
  }

  raw_data <- Read10X_h5(file.path(original_run_dir, "outs", "raw_feature_bc_matrix.h5"))
  if (!guide_feature_name %in% names(raw_data)) {
    stop(
      sprintf(
        "Guide feature set '%s' was not found in %s. Available slots: %s",
        guide_feature_name,
        original_run_dir,
        paste(names(raw_data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  gRNA_counts <- raw_data[[guide_feature_name]]

  common_barcodes <- intersect(colnames(rna_counts), colnames(gRNA_counts))
  if (length(common_barcodes) == 0) {
    stop(sprintf("No shared barcodes were found for sample %s.", sample_key), call. = FALSE)
  }

  rna_counts <- rna_counts[, common_barcodes, drop = FALSE]
  gRNA_counts <- gRNA_counts[, common_barcodes, drop = FALSE]

  seurat_obj <- CreateSeuratObject(
    counts = rna_counts,
    min.cells = 3,
    names.field = 2,
    names.delim = "-"
  )
  seurat_obj[["gRNA"]] <- CreateAssayObject(counts = gRNA_counts)
  seurat_obj$batch <- sample_key
  seurat_obj$sample_key <- sample_key
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  cat(sprintf(
    "[LOAD] %s: retained %d shared barcodes, %d RNA features, %d guide features\n",
    sample_key,
    ncol(seurat_obj),
    nrow(seurat_obj[["RNA"]]),
    nrow(seurat_obj[["gRNA"]])
  ))

  seurat_obj
}

batch_objects <- vector("list", nrow(sample_sheet))
for (i in seq_len(nrow(sample_sheet))) {
  batch_objects[[i]] <- load_batch(
    gex_dir = sample_sheet$gex_dir[[i]],
    original_run_dir = sample_sheet$original_run_dir[[i]],
    sample_key = sample_sheet$sample_key[[i]],
    guide_feature_name = gRNA_feature_name
  )
}

combined <- merge(
  batch_objects[[1]],
  y = batch_objects[-1],
  add.cell.id = sample_sheet$sample_key
)

saveRDS(combined, file = file.path(output_dir, "perturb_seq_merged_unfiltered.rds"))
cat(sprintf("[MERGE] Combined object contains %d cells before QC.\n", ncol(combined)))

###############################################################################
# STEP 3. Apply RNA-level QC and build the initial clustering used for manual
# population pruning.
#
# These thresholds match the legacy preprocessing workflow that defined the
# starting pool for endothelial-focused reanalysis.
###############################################################################

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
combined[["percent.ribo"]] <- PercentageFeatureSet(combined, pattern = "^Rp[ls]")

qc_violin <- VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  log = TRUE,
  pt.size = 0,
  group.by = "batch"
) + NoLegend()
ggsave(
  filename = file.path(plot_dir, "initial_qc_violin_by_batch.pdf"),
  plot = qc_violin,
  width = 12,
  height = 5
)

combined <- subset(
  combined,
  subset =
    nFeature_RNA > qc_min_features &
    nFeature_RNA < qc_max_features &
    nCount_RNA > qc_min_counts &
    nCount_RNA < qc_max_counts &
    percent.mt < qc_max_percent_mt
)

cat(sprintf("[QC] %d cells remain after RNA-level QC filtering.\n", ncol(combined)))

combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = initial_pcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:initial_pcs, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:initial_pcs, verbose = FALSE)
combined <- FindClusters(combined, resolution = initial_resolution, verbose = FALSE)

initial_umap <- DimPlot(combined, reduction = "umap", label = TRUE) + NoLegend()
ggsave(
  filename = file.path(plot_dir, "initial_umap_after_qc.pdf"),
  plot = initial_umap,
  width = 8,
  height = 6
)

saveRDS(combined, file = file.path(output_dir, "perturb_seq_after_initial_qc_and_clustering.rds"))

###############################################################################
# STEP 4. Reproduce the manual endothelial-focused cluster pruning from the
# legacy screen preprocessing.
#
# The cluster IDs below come directly from the original exploratory script. They
# should be treated as an explicit frozen preprocessing choice for the released
# manuscript-supporting workflow rather than a generic recipe for new data.
###############################################################################

combined <- subset(combined, idents = clusters_to_drop_round1, invert = TRUE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:recluster_pcs, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:recluster_pcs, verbose = FALSE)
combined <- FindClusters(combined, resolution = recluster_resolution, verbose = FALSE)

combined <- subset(combined, idents = clusters_to_drop_round2, invert = TRUE)
combined <- subset(combined, idents = clusters_to_drop_round3, invert = TRUE)

m_s_genes <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
m_g2m_genes <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

combined <- CellCycleScoring(
  combined,
  s.features = unique(m_s_genes),
  g2m.features = unique(m_g2m_genes),
  set.ident = TRUE
)

Idents(combined) <- "seurat_clusters"
combined <- RenameIdents(combined, cluster_map_broad)
combined$cell.type.l1 <- as.character(Idents(combined))

Idents(combined) <- "seurat_clusters"
combined <- RenameIdents(combined, cluster_map_fine)
combined$cell.type.l2 <- as.character(Idents(combined))

Idents(combined) <- "seurat_clusters"
combined <- subset(combined, idents = clusters_to_keep_final)

focused_umap <- DimPlot(combined, reduction = "umap", group.by = "cell.type.l2", label = TRUE)
ggsave(
  filename = file.path(plot_dir, "endothelial_focused_umap_before_doublet_filtering.pdf"),
  plot = focused_umap,
  width = 10,
  height = 7
)

saveRDS(combined, file = file.path(output_dir, "perturb_seq_endothelial_focused_pre_doublet.rds"))
cat(sprintf("[FOCUS] %d cells remain after manual endothelial-focused pruning.\n", ncol(combined)))

###############################################################################
# STEP 5. Remove doublets with scDblFinder in batch-aware mode.
#
# Doublet scoring is run after the endothelial-focused pruning step, matching
# the sequence used in the legacy preprocessing workflow.
###############################################################################

sce <- as.SingleCellExperiment(combined)
sce <- scDblFinder(
  sce,
  samples = "batch",
  BPPARAM = MulticoreParam(workers),
  clusters = FALSE,
  dbr = doublet_rate
)

combined$scDblFinder.class <- colData(sce)$scDblFinder.class
combined$scDblFinder.score <- colData(sce)$scDblFinder.score

doublet_plot <- DimPlot(
  combined,
  group.by = c("scDblFinder.class", "cell.type.l2"),
  label = TRUE
) + NoLegend()
ggsave(
  filename = file.path(plot_dir, "doublet_classification_overview.pdf"),
  plot = doublet_plot,
  width = 12,
  height = 6
)

combined <- subset(combined, subset = scDblFinder.class == "singlet")
saveRDS(combined, file = file.path(output_dir, "perturb_seq_endothelial_singlets.rds"))
cat(sprintf("[DOUBLET] %d singlet cells remain after scDblFinder filtering.\n", ncol(combined)))

###############################################################################
# STEP 6. Estimate MOI, apply the final guide-complexity filter, and export the
# matrices used for downstream SCEPTRE and cNMF analyses.
#
# Guide assignment itself is not performed here. Instead, this script exports
# the filtered RNA and gRNA matrices plus covariates that are consumed by later
# SCEPTRE analyses.
###############################################################################

gRNA_counts <- GetAssayData(combined, assay = "gRNA", slot = "counts")
assigned_guides <- gRNA_counts > moi_umi_threshold
combined$MOI <- Matrix::colSums(assigned_guides)

moi_density <- ggplot(data.frame(MOI = combined$MOI), aes(x = MOI)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_minimal() +
  labs(x = "MOI", y = "Density", title = "Guide multiplicity before MOI filtering")
ggsave(
  filename = file.path(plot_dir, "moi_density_before_filtering.pdf"),
  plot = moi_density,
  width = 7,
  height = 5
)

combined <- subset(combined, subset = MOI < moi_max & nCount_gRNA > 0)
cat(sprintf("[MOI] %d cells remain after retaining nCount_gRNA > 0 and MOI < %d.\n", ncol(combined), moi_max))

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = sceptre_variable_genes,
  verbose = FALSE
)
variable_genes <- VariableFeatures(combined)

DefaultAssay(combined) <- "RNA"
response_matrix <- GetAssayData(combined, assay = "RNA", slot = "counts")
response_matrix <- response_matrix[variable_genes, , drop = FALSE]
genes_to_keep <- Matrix::rowSums(response_matrix > 0) > min_cells_per_export_gene
response_matrix <- response_matrix[genes_to_keep, , drop = FALSE]

gRNA_data <- GetAssayData(combined, assay = "gRNA", slot = "counts")

meta_data <- combined@meta.data %>%
  transmute(
    response_n_umis = nCount_RNA,
    response_n_nonzero = nFeature_RNA,
    grna_n_umis = nCount_gRNA,
    grna_n_nonzero = nFeature_gRNA,
    p_mito = percent.mt,
    batch = batch,
    cell.type.l1 = cell.type.l1,
    cell.type.l2 = cell.type.l2,
    doublet = scDblFinder.class,
    MOI = MOI
  )

save(
  response_matrix,
  gRNA_data,
  meta_data,
  file = file.path(output_dir, "neonatal_brain_perturbseq_sceptre_input.RData")
)

saveRDS(combined, file = file.path(output_dir, "neonatal_brain_perturbseq_preprocessed.rds"))

if (h5_export) {
  DefaultAssay(combined) <- "RNA"
  expression_counts <- GetAssayData(combined, assay = "RNA", slot = "counts")
  export_genes <- names(Matrix::rowSums(expression_counts > 0)[Matrix::rowSums(expression_counts > 0) > min_cells_per_export_gene])
  combined_export <- subset(combined, features = export_genes)
  combined_export@assays$RNA@scale.data <- matrix(nrow = 0, ncol = 0)
  combined_export$cell.type.l1 <- as.character(combined_export$cell.type.l1)
  combined_export$cell.type.l2 <- as.character(combined_export$cell.type.l2)
  combined_export$batch <- as.character(combined_export$batch)
  combined_export@assays$RNA@data <- combined_export@assays$RNA@counts

  h5seurat_path <- file.path(output_dir, "neonatal_brain_perturbseq_preprocessed.h5Seurat")
  SaveH5Seurat(combined_export, filename = h5seurat_path, overwrite = TRUE)
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
}

cell_summary <- data.frame(
  stage = c(
    "after_initial_merge",
    "after_rna_qc",
    "after_endothelial_pruning",
    "after_singlet_filtering",
    "after_moi_and_guide_filtering"
  ),
  n_cells = c(
    readRDS(file.path(output_dir, "perturb_seq_merged_unfiltered.rds")) |> ncol(),
    readRDS(file.path(output_dir, "perturb_seq_after_initial_qc_and_clustering.rds")) |> ncol(),
    readRDS(file.path(output_dir, "perturb_seq_endothelial_focused_pre_doublet.rds")) |> ncol(),
    readRDS(file.path(output_dir, "perturb_seq_endothelial_singlets.rds")) |> ncol(),
    ncol(combined)
  )
)
write.table(
  cell_summary,
  file = file.path(output_dir, "cell_count_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[DONE] Wrote manuscript-facing preprocessing outputs to:\n")
cat(sprintf("  %s\n", output_dir))
cat("[DONE] Key exports:\n")
cat("  - neonatal_brain_perturbseq_sceptre_input.RData\n")
cat("  - neonatal_brain_perturbseq_preprocessed.rds\n")
if (h5_export) {
  cat("  - neonatal_brain_perturbseq_preprocessed.h5Seurat\n")
  cat("  - neonatal_brain_perturbseq_preprocessed.h5ad\n")
}
