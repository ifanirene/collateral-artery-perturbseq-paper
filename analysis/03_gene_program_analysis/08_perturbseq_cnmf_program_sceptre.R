#!/usr/bin/env Rscript
# @description
# This script performs program-level SCEPTRE testing for cNMF-derived
# Perturb-seq endothelial gene programs.
# It is responsible for converting per-cell cNMF usage values into
# count-like program responses, aligning those responses with guide
# assignments and cell-level covariates, and running high-MOI SCEPTRE
# analyses across all endothelial cells or within cell-type subsets.
#
# Key features:
# - Uses descriptive relative input and output filenames for release
# - Converts row-normalized cNMF usage values into per-program count proxies
# - Reuses the gene-level SCEPTRE input bundle, guide annotations, and QC covariates
# - Supports all-cell and per-cell-type program-response testing
# - Drops unused factor levels after subsetting to avoid redundant formula terms
# - Records session information for reproducibility
#
# @dependencies
# - sceptre: CRISPR perturbation testing framework
# - Matrix: sparse matrix operations
# - dplyr: metadata manipulation
# - readr: guide annotation import
#
# @examples
# Rscript docs/fp_seq2_cnmf_program_level_sceptre_submission_clean_2026-04-01.R

suppressPackageStartupMessages({
  library(sceptre)
  library(Matrix)
  library(dplyr)
  library(readr)
})

message("Session information at script start:")
message(paste(capture.output(sessionInfo()), collapse = "\n"))

## ------------------
## Release inputs and outputs
## ------------------
input_dir <- "input"
output_root <- file.path("output", "perturbseq_program_level_sceptre")

gene_level_input_bundle_path <- file.path(
  input_dir,
  "perturbseq_gene_level_sceptre_input_bundle.RData"
)
guide_annotation_path <- file.path(
  input_dir,
  "perturbseq_guide_target_annotation.csv"
)
program_usage_matrix_path <- file.path(
  input_dir,
  "perturbseq_cnmf_program_usage_k100_dt02_consensus.tsv"
)

analysis_run_label <- "perturbseq_cnmf_program_sceptre_k100_dt02"
expected_program_count <- 100L
local_density_threshold <- 0.2
thread_count <- 12L
run_per_cell_type <- FALSE
cell_type_column <- "cell.type"
cell_types_to_run <- NULL
min_cells_per_cell_type <- 200L

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

## ------------------
## Utility helpers
## ------------------
row_normalize <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  sweep(mat, 1, row_sums, FUN = "/")
}

drop_unused_factor_columns <- function(df) {
  df[] <- lapply(df, function(column) {
    if (is.factor(column)) {
      droplevels(column)
    } else {
      column
    }
  })
  df
}

run_sceptre_program_analysis <- function(
  program_response_matrix,
  guide_count_matrix,
  guide_annotation_table,
  extra_covariate_table,
  model_covariate_table,
  threads,
  output_directory
) {
  message("[SCEPTRE] Importing program-level response matrix")
  sceptre_object <- import_data(
    response_matrix = program_response_matrix,
    grna_matrix = guide_count_matrix,
    grna_target_data_frame = guide_annotation_table,
    extra_covariates = extra_covariate_table,
    moi = "high"
  )

  sceptre_object@covariate_data_frame <- model_covariate_table

  positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
  discovery_pairs <- construct_trans_pairs(
    sceptre_object,
    positive_control_pairs,
    pairs_to_exclude = "none"
  )

  model_formula <- ~ log(response_n_nonzero) + log(response_n_umis) +
    log(grna_n_nonzero) + log(grna_n_umis) + batch

  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs,
      formula_object = model_formula,
      side = "both"
    ) |>
    assign_grnas(
      method = "thresholding",
      threshold = 10,
      parallel = TRUE,
      n_processors = threads
    ) |>
    run_qc() |>
    run_calibration_check(
      parallel = TRUE,
      n_processors = threads,
      output_amount = 2
    ) |>
    run_discovery_analysis(
      parallel = TRUE,
      n_processors = threads,
      output_amount = 2
    )

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  write_outputs_to_directory(sceptre_object, directory = output_directory)
  invisible(sceptre_object)
}

## ------------------
## Load release inputs
## ------------------
if (!file.exists(gene_level_input_bundle_path)) {
  stop("Missing gene-level SCEPTRE input bundle: ", gene_level_input_bundle_path)
}
if (!file.exists(guide_annotation_path)) {
  stop("Missing guide annotation table: ", guide_annotation_path)
}
if (!file.exists(program_usage_matrix_path)) {
  stop("Missing cNMF program usage matrix: ", program_usage_matrix_path)
}

load(gene_level_input_bundle_path)
expected_objects <- c("response_matrix", "gRNA_data", "meta_data")
missing_objects <- expected_objects[!vapply(expected_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Expected objects missing from gene-level input bundle: ",
    paste(missing_objects, collapse = ", ")
  )
}

guide_annotation_table <- read_csv(
  guide_annotation_path,
  show_col_types = FALSE
)
colnames(guide_annotation_table) <- c("grna_id", "grna_target")
guide_annotation_table$grna_id <- gsub("_", "-", guide_annotation_table$grna_id)
guide_annotation_table <- distinct(guide_annotation_table)

message("Reading program usage matrix from ", program_usage_matrix_path)
program_usage_table <- read.table(
  program_usage_matrix_path,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (ncol(program_usage_table) != expected_program_count) {
  message(
    "Program-usage matrix contains ",
    ncol(program_usage_table),
    " programs (expected ",
    expected_program_count,
    ")"
  )
}

## ------------------
## Align cells across matrices
## ------------------
raw_expression_matrix <- response_matrix
shared_cell_barcodes <- Reduce(
  intersect,
  list(
    rownames(program_usage_table),
    colnames(raw_expression_matrix),
    rownames(meta_data),
    colnames(gRNA_data)
  )
)
if (length(shared_cell_barcodes) == 0) {
  stop("No shared cell barcodes across program usage, expression, metadata, and guide counts")
}

program_usage_table <- program_usage_table[shared_cell_barcodes, , drop = FALSE]
raw_expression_matrix <- raw_expression_matrix[, shared_cell_barcodes, drop = FALSE]
meta_data <- meta_data[shared_cell_barcodes, , drop = FALSE]
gRNA_data <- gRNA_data[, shared_cell_barcodes, drop = FALSE]

message(
  "Aligned inputs: ",
  nrow(program_usage_table), " cells, ",
  ncol(program_usage_table), " programs, ",
  nrow(guide_annotation_table), " guide annotations"
)

## ------------------
## Build count-like program responses
## ------------------
gene_expression_umi_totals <- Matrix::colSums(raw_expression_matrix)
gene_expression_umi_totals <- gene_expression_umi_totals[
  match(rownames(program_usage_table), names(gene_expression_umi_totals))
]
stopifnot(!anyNA(gene_expression_umi_totals))

normalized_program_usage <- row_normalize(as.matrix(program_usage_table))
program_count_proxy_matrix <- round(normalized_program_usage * gene_expression_umi_totals)
program_response_matrix <- as(t(program_count_proxy_matrix), "dgCMatrix")

## ------------------
## Build covariate tables
## ------------------
if (!cell_type_column %in% colnames(meta_data)) {
  stop("Requested cell-type column not found in metadata: ", cell_type_column)
}

meta_data <- meta_data %>%
  mutate(
    batch = factor(batch),
    !!cell_type_column := factor(.data[[cell_type_column]])
  )

extra_covariate_table <- meta_data %>%
  transmute(
    p_mito = p_mito,
    batch = batch
  )

model_covariate_table <- meta_data

## ------------------
## Run all-cell or per-cell-type analyses
## ------------------
result_root <- file.path(output_root, analysis_run_label)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

if (!run_per_cell_type) {
  message(
    "[GLOBAL] Running program-level SCEPTRE on ",
    ncol(program_response_matrix),
    " endothelial cells"
  )
  run_sceptre_program_analysis(
    program_response_matrix = program_response_matrix,
    guide_count_matrix = gRNA_data,
    guide_annotation_table = guide_annotation_table,
    extra_covariate_table = extra_covariate_table,
    model_covariate_table = model_covariate_table,
    threads = thread_count,
    output_directory = file.path(result_root, "all_endothelial_cells")
  )
} else {
  cell_type_levels <- levels(meta_data[[cell_type_column]])
  if (!is.null(cell_types_to_run)) {
    cell_type_levels <- intersect(cell_type_levels, cell_types_to_run)
  }

  for (cell_type_label in cell_type_levels) {
    keep_cells <- meta_data[[cell_type_column]] == cell_type_label
    if (sum(keep_cells) < min_cells_per_cell_type) {
      message(
        "[", cell_type_label, "] skipped (",
        sum(keep_cells),
        " cells < ",
        min_cells_per_cell_type,
        ")"
      )
      next
    }

    message(
      "[", cell_type_label, "] running program-level SCEPTRE on ",
      sum(keep_cells),
      " cells"
    )

    subset_extra_covariates <- drop_unused_factor_columns(
      extra_covariate_table[keep_cells, , drop = FALSE]
    )
    subset_model_covariates <- drop_unused_factor_columns(
      model_covariate_table[keep_cells, , drop = FALSE]
    )
    if (cell_type_column %in% colnames(subset_model_covariates)) {
      subset_model_covariates[[cell_type_column]] <- NULL
    }

    run_sceptre_program_analysis(
      program_response_matrix = program_response_matrix[, keep_cells, drop = FALSE],
      guide_count_matrix = gRNA_data[, keep_cells, drop = FALSE],
      guide_annotation_table = guide_annotation_table,
      extra_covariate_table = subset_extra_covariates,
      model_covariate_table = subset_model_covariates,
      threads = thread_count,
      output_directory = file.path(
        result_root,
        paste0("cell_type_", make.names(cell_type_label))
      )
    )
  }
}

writeLines(
  capture.output(sessionInfo()),
  file.path(result_root, "session_info_program_level_sceptre.txt")
)

message(
  "Completed program-level SCEPTRE testing for ",
  nrow(program_response_matrix),
  " programs at density threshold ",
  local_density_threshold,
  ". Results written to ",
  result_root
)
