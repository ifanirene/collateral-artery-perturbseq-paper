#!/usr/bin/env Rscript
# """
# @description
# This script derives endothelial-state UCell scores from a filtered
# Perturb-seq Seurat object with endothelial cell-state labels. It regroups
# arterial subclusters, discovers de novo markers for each endothelial state,
# selects signature genes, computes per-cell UCell scores, and exports score
# tables for downstream perturbation testing.
#
# Key features:
# - Regroups `artery_1` and `artery_2` into a single `artery` state before
#   marker discovery.
# - Selects state-specific marker signatures with the same ranking and fallback
#   rules used in the maintained analysis.
# - Exports both a per-cell score table and a signatures-by-cells sparse matrix
#   for downstream SCEPTRE testing.
#
# @dependencies
# - Seurat: single-cell container, normalization, and marker discovery.
# - UCell: rank-based gene-set scoring.
# - dplyr, Matrix: data wrangling and sparse matrix conversion.
#
# @examples
# - Source this script in RStudio and update the user-editable input block.
# - Run section by section to inspect marker tables before computing UCell scores.
# """

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(UCell)
})

## -------------------------------------------------------------------------
## User-editable inputs and output names
## -------------------------------------------------------------------------
# Relative path to the filtered Perturb-seq Seurat object with EC state labels.
input_seurat_rds <- "data/filtered_perturbseq_seurat_object_with_ec_state_labels.rds"

# Metadata column that stores the EC state labels used for marker discovery.
state_label_column <- "cell.type.l2"

# Name of the regrouped label column that will be added to the Seurat object.
regrouped_state_column <- "cell.type.l2.regrouped"

endothelial_state_labels <- c("artery", "pre-artery", "cap_2", "cap_3", "cap_4", "vein")
priority_state_labels <- c("artery", "pre-artery", "vein")
top_n_default <- 30L
top_n_priority <- 50L
minimum_signature_size <- 10L
random_seed <- 1L

output_dir <- "output/endothelial_state_scoring"
markers_all_csv <- file.path(output_dir, "endothelial_state_markers_all.csv")
signature_markers_csv <- file.path(output_dir, "endothelial_state_signature_markers.csv")
signature_summary_csv <- file.path(output_dir, "endothelial_state_signature_summary.csv")
signature_gene_sets_rds <- file.path(output_dir, "endothelial_state_signature_gene_sets.rds")
scored_seurat_rds <- file.path(output_dir, "endothelial_state_scored_seurat.rds")
ucell_scores_matrix_rds <- file.path(output_dir, "endothelial_state_ucell_scores_matrix.rds")
ucell_scores_csv <- file.path(output_dir, "endothelial_state_ucell_scores.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_seurat_rds)) {
  stop(
    "Input Seurat object not found: ", input_seurat_rds,
    "\nUpdate `input_seurat_rds` so it points to the filtered Perturb-seq Seurat object with EC state labels.",
    call. = FALSE
  )
}

cat("Session info at script start:\n")
print(sessionInfo())

## -------------------------------------------------------------------------
## Load the Seurat object and define the regrouped endothelial states
## This block creates the identity labels used for endothelial-state scoring.
## -------------------------------------------------------------------------
set.seed(random_seed)
message("Reading Seurat object: ", input_seurat_rds)
seu <- readRDS(input_seurat_rds)

if (!inherits(seu, "Seurat")) {
  stop("Input RDS is not a Seurat object. Class: ", paste(class(seu), collapse = ", "), call. = FALSE)
}

metadata_columns <- colnames(seu@meta.data)
if (!(state_label_column %in% metadata_columns)) {
  stop(
    "Missing required metadata column '", state_label_column,
    "'. Available columns: ", paste(metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

original_state_labels <- as.character(seu[[state_label_column]][, 1])
regrouped_state_labels <- original_state_labels
regrouped_state_labels[original_state_labels %in% c("artery_1", "artery_2")] <- "artery"

cells_to_keep <- !is.na(regrouped_state_labels) & regrouped_state_labels != ""
if (!all(cells_to_keep)) {
  message("Dropping ", sum(!cells_to_keep), " cells with missing state labels after regrouping.")
  seu <- subset(seu, cells = colnames(seu)[cells_to_keep])
  regrouped_state_labels <- regrouped_state_labels[cells_to_keep]
}

if ("RNA" %in% names(seu@assays)) {
  DefaultAssay(seu) <- "RNA"
} else {
  DefaultAssay(seu) <- names(seu@assays)[[1]]
}

message("Using assay: ", DefaultAssay(seu))
message("Normalizing assay data before marker discovery.")
seu <- NormalizeData(seu, verbose = FALSE)

seu[[regrouped_state_column]] <- droplevels(factor(regrouped_state_labels))
Idents(seu) <- seu[[regrouped_state_column]][, 1]

available_state_labels <- levels(Idents(seu))
retained_state_labels <- intersect(endothelial_state_labels, available_state_labels)
missing_state_labels <- setdiff(endothelial_state_labels, retained_state_labels)

if (length(retained_state_labels) == 0L) {
  stop(
    "None of the requested endothelial states were found after regrouping. Requested: ",
    paste(endothelial_state_labels, collapse = ", "),
    ". Available labels: ", paste(available_state_labels, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_state_labels) > 0L) {
  message("Skipping missing endothelial states: ", paste(missing_state_labels, collapse = ", "))
}

message("Regrouped state counts:")
print(sort(table(as.character(Idents(seu))), decreasing = TRUE))
message("Endothelial states retained for scoring: ", paste(retained_state_labels, collapse = ", "))

## -------------------------------------------------------------------------
## Find markers across regrouped endothelial states
## This block reproduces the ROC-based marker discovery used for scoring.
## -------------------------------------------------------------------------
message("Running Seurat::FindAllMarkers(test.use = 'roc') on regrouped endothelial states.")
markers <- FindAllMarkers(
  object = seu,
  only.pos = TRUE,
  test.use = "roc",
  min.pct = 0.1,
  logfc.threshold = 0,
  random.seed = random_seed,
  verbose = TRUE
)

if (nrow(markers) == 0L) {
  stop("FindAllMarkers returned no markers.", call. = FALSE)
}

required_marker_columns <- c("cluster", "gene", "avg_log2FC", "pct.1")
missing_marker_columns <- setdiff(required_marker_columns, colnames(markers))
if (length(missing_marker_columns) > 0L) {
  stop(
    "FindAllMarkers output missing required columns: ",
    paste(missing_marker_columns, collapse = ", "),
    call. = FALSE
  )
}

if (!("myAUC" %in% colnames(markers))) {
  if ("AUC" %in% colnames(markers)) {
    markers$myAUC <- markers$AUC
  } else {
    stop("FindAllMarkers output missing myAUC/AUC required for ranking.", call. = FALSE)
  }
}

if (!("pct.2" %in% colnames(markers))) {
  markers$pct.2 <- NA_real_
}

markers <- markers %>%
  mutate(
    cluster = as.character(cluster),
    gene = as.character(gene),
    avg_log2FC = as.numeric(avg_log2FC),
    myAUC = as.numeric(myAUC),
    pct.1 = as.numeric(pct.1),
    pct.2 = as.numeric(pct.2)
  ) %>%
  filter(!is.na(cluster), !is.na(gene), gene != "")

markers_clean <- markers %>%
  filter(!grepl("^(mt-|rpl|rps)", gene, ignore.case = TRUE))

markers_all <- markers_clean %>%
  filter(cluster %in% retained_state_labels) %>%
  arrange(cluster, desc(myAUC), desc(avg_log2FC), desc(pct.1), gene)

write.csv(markers_all, markers_all_csv, row.names = FALSE, quote = FALSE)
message("Wrote full endothelial-state marker table: ", markers_all_csv)

## -------------------------------------------------------------------------
## Select signature genes for each endothelial state
## This block applies the maintained ranking, cutoff, and fallback rules.
## -------------------------------------------------------------------------
selected_marker_tables <- list()
signature_summary_tables <- list()

for (state_label in retained_state_labels) {
  top_n_target <- if (state_label %in% priority_state_labels) top_n_priority else top_n_default
  log2fc_cutoff <- if (state_label %in% c("cap_2", "cap_3", "cap_4")) 1.0 else 1.5

  state_markers_ranked <- markers_clean %>%
    filter(cluster == state_label) %>%
    arrange(desc(myAUC), desc(avg_log2FC), desc(pct.1), gene) %>%
    distinct(gene, .keep_all = TRUE)

  state_markers_thresholded <- state_markers_ranked %>%
    filter(avg_log2FC > log2fc_cutoff) %>%
    arrange(desc(myAUC), desc(avg_log2FC), desc(pct.1), gene) %>%
    distinct(gene, .keep_all = TRUE)

  selected_markers <- state_markers_thresholded %>% head(top_n_target)
  used_minimum_size_fallback <- FALSE
  used_top_n_fallback <- FALSE

  if (nrow(selected_markers) < minimum_signature_size) {
    n_needed <- minimum_signature_size - nrow(selected_markers)
    fill_markers <- state_markers_ranked %>%
      filter(!(gene %in% selected_markers$gene)) %>%
      head(n_needed)
    if (nrow(fill_markers) > 0L) {
      used_minimum_size_fallback <- TRUE
      selected_markers <- bind_rows(selected_markers, fill_markers)
    }
  }

  if (nrow(selected_markers) < top_n_target) {
    n_needed <- top_n_target - nrow(selected_markers)
    fill_markers <- state_markers_ranked %>%
      filter(!(gene %in% selected_markers$gene)) %>%
      head(n_needed)
    if (nrow(fill_markers) > 0L) {
      used_top_n_fallback <- TRUE
      selected_markers <- bind_rows(selected_markers, fill_markers)
    }
  }

  selected_markers <- selected_markers %>%
    distinct(gene, .keep_all = TRUE) %>%
    head(top_n_target)

  if (nrow(selected_markers) == 0L) {
    stop("No markers available for endothelial state: ", state_label, call. = FALSE)
  }

  selected_markers <- selected_markers %>%
    mutate(
      state_label = state_label,
      rank_selected = row_number(),
      top_n_target = top_n_target,
      threshold_cutoff = log2fc_cutoff,
      priority_label = state_label %in% priority_state_labels,
      fallback_minimum_signature_size = used_minimum_size_fallback,
      fallback_top_n = used_top_n_fallback,
      low_marker_flag = nrow(selected_markers) < minimum_signature_size
    )

  selected_marker_tables[[state_label]] <- selected_markers
  signature_summary_tables[[state_label]] <- data.frame(
    state_label = state_label,
    top_n_target = top_n_target,
    priority_label = state_label %in% priority_state_labels,
    threshold_cutoff = log2fc_cutoff,
    n_markers_available = nrow(state_markers_ranked),
    n_markers_thresholded = nrow(state_markers_thresholded),
    n_markers_selected = nrow(selected_markers),
    fallback_minimum_signature_size = used_minimum_size_fallback,
    fallback_top_n = used_top_n_fallback,
    stringsAsFactors = FALSE
  )
}

signature_markers <- bind_rows(selected_marker_tables) %>%
  select(
    state_label,
    gene,
    rank_selected,
    top_n_target,
    priority_label,
    myAUC,
    avg_log2FC,
    pct.1,
    pct.2,
    threshold_cutoff,
    fallback_minimum_signature_size,
    fallback_top_n,
    low_marker_flag
  )

signature_summary <- bind_rows(signature_summary_tables)
signature_gene_sets <- lapply(selected_marker_tables, function(marker_table) unique(marker_table$gene))

write.csv(signature_markers, signature_markers_csv, row.names = FALSE, quote = FALSE)
write.csv(signature_summary, signature_summary_csv, row.names = FALSE, quote = FALSE)
saveRDS(signature_gene_sets, signature_gene_sets_rds)

message("Wrote signature marker table: ", signature_markers_csv)
message("Wrote signature summary table: ", signature_summary_csv)
message("Wrote signature gene sets: ", signature_gene_sets_rds)

## -------------------------------------------------------------------------
## Compute per-cell UCell scores for each endothelial state
## This block exports both a per-cell score table and a sparse response matrix.
## -------------------------------------------------------------------------
message("Computing UCell scores for ", length(signature_gene_sets), " endothelial states.")
seu <- AddModuleScore_UCell(seu, features = signature_gene_sets)

available_score_columns <- colnames(seu@meta.data)
resolved_score_columns <- vapply(
  names(signature_gene_sets),
  function(state_label) {
    candidate_columns <- c(
      paste0(state_label, "_UCell"),
      paste0(make.names(state_label), "_UCell"),
      paste0(gsub("-", ".", state_label, fixed = TRUE), "_UCell")
    )
    matching_columns <- candidate_columns[candidate_columns %in% available_score_columns]
    if (length(matching_columns) == 0L) {
      return(NA_character_)
    }
    matching_columns[[1]]
  },
  FUN.VALUE = character(1)
)

missing_score_columns <- names(signature_gene_sets)[is.na(resolved_score_columns)]
if (length(missing_score_columns) > 0L) {
  stop(
    "Missing expected UCell score columns for endothelial states: ",
    paste(missing_score_columns, collapse = ", "),
    call. = FALSE
  )
}

expected_score_columns <- paste0(names(signature_gene_sets), "_UCell")
ucell_score_table <- seu@meta.data[, resolved_score_columns, drop = FALSE]
colnames(ucell_score_table) <- expected_score_columns

ucell_score_matrix <- Matrix(t(as.matrix(ucell_score_table)), sparse = TRUE)
rownames(ucell_score_matrix) <- expected_score_columns
colnames(ucell_score_matrix) <- rownames(ucell_score_table)

saveRDS(seu, scored_seurat_rds)
saveRDS(ucell_score_matrix, ucell_scores_matrix_rds)

ucell_score_export <- data.frame(
  cell_barcode = rownames(ucell_score_table),
  ucell_score_table,
  check.names = FALSE
)
write.csv(ucell_score_export, ucell_scores_csv, row.names = FALSE, quote = FALSE)

message("Wrote scored Seurat object: ", scored_seurat_rds)
message("Wrote signatures-by-cells UCell score matrix: ", ucell_scores_matrix_rds)
message("Wrote per-cell UCell score table: ", ucell_scores_csv)
