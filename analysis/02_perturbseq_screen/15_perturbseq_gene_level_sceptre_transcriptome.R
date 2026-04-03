#!/usr/bin/env Rscript
# @description
# This script reproduces the gene-level SCEPTRE analysis for the
# Final Pool MOI<15 Perturb-seq transcriptome screen.
# It is responsible for loading the maintained transcriptome and guide matrices,
# reconstructing the global and per-endothelial-subtype high-MOI models, and
# writing clearly labeled SCEPTRE outputs together with analysis summaries.
#
# Key features:
# - Uses the maintained MOI15 singlet transcriptome inputs shared across final SCEPTRE runs
# - Preserves the final guide-assignment threshold, covariates, calibration size, and subtype split
# - Writes descriptive run summaries alongside the standard SCEPTRE outputs
# - Supports a lightweight dry-run mode so the script can be validated without rerunning the full screen
#
# @dependencies
# - sceptre: CRISPR screen testing framework
# - Matrix: sparse count matrices
# - dplyr: metadata handling
# - parallelly: cgroup-aware core detection for safe parallel defaults (optional)
#
# @examples
# Rscript src/gene_level_sceptre_transcriptome_moi15.R \
#   --gene-guide-input path/to/gene_and_guide_counts.RData \
#   --technical-covariates path/to/technical_covariates.RData \
#   --guide-target-annotation path/to/guide_target_annotation.csv \
#   --output-root path/to/output_directory \
#   --dry-run
# Rscript src/gene_level_sceptre_transcriptome_moi15.R \
#   --gene-guide-input path/to/gene_and_guide_counts.RData \
#   --technical-covariates path/to/technical_covariates.RData \
#   --guide-target-annotation path/to/guide_target_annotation.csv \
#   --output-root path/to/output_directory

suppressPackageStartupMessages({
  library(sceptre)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
run_per_celltype <- !("--all-cells-only" %in% args)

# Final analysis settings recovered from the maintained transcriptome SCEPTRE runs.
assign_threshold <- 10L
n_calibration_pairs <- 200000L
min_cells_per_type <- 200L
threads <- if (requireNamespace("parallelly", quietly = TRUE)) {
  max(1L, min(16L, as.integer(parallelly::availableCores())))
} else {
  8L
}

get_option_value <- function(flag, args) {
  flag_index <- match(flag, args)
  if (is.na(flag_index)) {
    return(NULL)
  }
  if (flag_index == length(args)) {
    stop("Missing value for ", flag)
  }
  args[[flag_index + 1]]
}

gene_and_guide_count_input_rdata <- get_option_value("--gene-guide-input", args)
technical_covariate_input_rdata <- get_option_value("--technical-covariates", args)
guide_target_annotation_csv <- get_option_value("--guide-target-annotation", args)
transcriptome_sceptre_output_dir <- get_option_value("--output-root", args)

required_cli_options <- c(
  "--gene-guide-input" = gene_and_guide_count_input_rdata,
  "--technical-covariates" = technical_covariate_input_rdata,
  "--guide-target-annotation" = guide_target_annotation_csv,
  "--output-root" = transcriptome_sceptre_output_dir
)
missing_cli_options <- names(required_cli_options)[vapply(required_cli_options, is.null, logical(1))]
if (length(missing_cli_options) > 0) {
  stop(
    "Missing required argument(s): ",
    paste(missing_cli_options, collapse = ", "),
    "\nProvide explicit input and output locations on the command line."
  )
}

drop_unused_factor_levels <- function(df) {
  for (nm in names(df)) {
    if (is.factor(df[[nm]])) {
      df[[nm]] <- droplevels(df[[nm]])
    }
  }
  df
}

build_formula_object <- function(covars) {
  formula_terms <- c(
    "log(response_n_nonzero)",
    "log(response_n_umis)",
    "log(grna_n_nonzero)",
    "log(grna_n_umis)"
  )

  if ("batch" %in% colnames(covars)) {
    batch_n_levels <- if (is.factor(covars$batch)) {
      nlevels(droplevels(covars$batch))
    } else {
      dplyr::n_distinct(covars$batch)
    }
    if (batch_n_levels > 1) {
      formula_terms <- c(formula_terms, "batch")
    }
  }

  as.formula(paste("~", paste(formula_terms, collapse = " + ")))
}

write_design_summary <- function(out_dir, subset_label, payload, formula_object, positive_control_pairs = NA_integer_, discovery_pairs = NA_integer_) {
  summary_lines <- c(
    paste0("analysis_label: ", subset_label),
    paste0("response_matrix: gene-level transcriptome counts"),
    paste0("guide_matrix: gRNA UMI counts"),
    paste0("cells_entering_analysis: ", ncol(payload$response_mat)),
    paste0("response_genes_entering_analysis: ", nrow(payload$response_mat)),
    paste0("guides_entering_analysis: ", nrow(payload$gRNA_mat)),
    paste0("guide_assignment_threshold_umis: ", assign_threshold),
    paste0("negative_control_pairs_for_calibration: ", n_calibration_pairs),
    paste0("formula_object: ", paste(deparse(formula_object), collapse = "")),
    paste0("positive_control_pairs_before_pairwise_qc: ", positive_control_pairs),
    paste0("discovery_pairs_before_pairwise_qc: ", discovery_pairs),
    paste0("dry_run: ", dry_run)
  )
  writeLines(summary_lines, file.path(out_dir, "analysis_inputs_and_model_summary.txt"))
}

run_sceptre_subset <- function(subset_label, payload, guide_map, out_dir) {
  formula_object <- build_formula_object(payload$covars)

  message("[", subset_label, "] cells entering analysis: ", ncol(payload$response_mat))
  message("[", subset_label, "] response genes entering analysis: ", nrow(payload$response_mat))
  message("[", subset_label, "] model formula: ", paste(deparse(formula_object), collapse = ""))

  if (dry_run) {
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: import the maintained transcriptome and guide matrices.
  # Input: gene-by-cell transcriptome counts, guide-by-cell counts, and the curated guide-to-target table.
  # Why: SCEPTRE requires the response matrix, perturbation matrix, and guide annotation to define the testable design.
  # Output: a high-MOI sceptre object ready to receive the precomputed technical covariates.
  message("[", subset_label, "] import_data()")
  sce <- import_data(
    response_matrix = payload$response_mat,
    grna_matrix = payload$gRNA_mat,
    grna_target_data_frame = guide_map,
    extra_covariates = payload$extra_covars,
    moi = "high"
  )

  # Step 2: attach the precomputed technical covariates from preprocessing.
  # Input: the aligned covariate table containing RNA and guide library metrics plus batch.
  # Why: the final analysis regresses out library size and detection-rate terms already computed upstream.
  # Output: a sceptre object carrying the exact technical covariate table used for the final model.
  sce@covariate_data_frame <- payload$covars

  # Step 3: define the positive-control and discovery pair sets.
  # Input: the imported high-MOI sceptre object.
  # Why: cis guide-target pairs provide the on-target positive controls, while trans pairs define transcriptome-wide discovery testing.
  # Output: two pair tables that fully specify the hypothesis set for calibration and discovery.
  message("[", subset_label, "] construct_positive_control_pairs()")
  positive_control_pairs <- construct_positive_control_pairs(sce)
  message("[", subset_label, "] positive-control pairs before pairwise QC: ", nrow(positive_control_pairs))

  message("[", subset_label, "] construct_trans_pairs()")
  discovery_pairs <- construct_trans_pairs(
    sceptre_object = sce,
    positive_control_pairs = positive_control_pairs,
    pairs_to_exclude = "none"
  )
  message("[", subset_label, "] discovery pairs before pairwise QC: ", nrow(discovery_pairs))

  write_design_summary(
    out_dir = out_dir,
    subset_label = subset_label,
    payload = payload,
    formula_object = formula_object,
    positive_control_pairs = nrow(positive_control_pairs),
    discovery_pairs = nrow(discovery_pairs)
  )

  # Step 4: run guide assignment, QC, calibration, and discovery in the same order as the maintained final runs.
  # Input: the imported object plus the pair definitions and model formula.
  # Why: this recreates the final gene-level screen analysis used for transcriptome-wide perturbation testing.
  # Output: the standard sceptre result files, including calibration and discovery RDS outputs.
  message("[", subset_label, "] set_analysis_parameters()")
  sce <- sce |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs,
      formula_object = formula_object,
      side = "both"
    )

  message("[", subset_label, "] assign_grnas() with threshold = ", assign_threshold, " UMIs")
  sce <- sce |>
    assign_grnas(
      method = "thresholding",
      threshold = assign_threshold,
      parallel = TRUE,
      n_processors = threads
    ) |>
    run_qc()

  message("[", subset_label, "] run_calibration_check() with ", n_calibration_pairs, " negative-control pairs")
  sce <- sce |>
    run_calibration_check(
      parallel = TRUE,
      n_processors = threads,
      n_calibration_pairs = n_calibration_pairs,
      output_amount = 2
    )

  message("[", subset_label, "] run_discovery_analysis()")
  sce <- sce |>
    run_discovery_analysis(
      parallel = TRUE,
      n_processors = threads,
      output_amount = 2
    )

  message("[", subset_label, "] write_outputs_to_directory() -> ", out_dir)
  write_outputs_to_directory(sce, directory = out_dir)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

  invisible(sce)
}

message("Dry run: ", dry_run)
message("Run per cell type: ", run_per_celltype)
message("Threads: ", threads)
message("Gene-level transcriptome SCEPTRE session info:")
writeLines(capture.output(sessionInfo()))

if (!file.exists(gene_and_guide_count_input_rdata)) {
  stop("Gene-level count input RData not found at the provided path.")
}
if (!file.exists(technical_covariate_input_rdata)) {
  stop("Technical covariate RData not found at the provided path.")
}
if (!file.exists(guide_target_annotation_csv)) {
  stop("Guide annotation CSV not found at the provided path.")
}

# Step 0: load the maintained MOI15 transcriptome inputs.
# Input: the saved gene response matrix, guide matrix, metadata table, and covariate table.
# Why: these objects are the stable starting point shared across the final transcriptome SCEPTRE runs.
# Output: in-memory matrices and metadata ready for barcode alignment and subtype splitting.
message("Loading the maintained MOI15 transcriptome inputs")
load(gene_and_guide_count_input_rdata)
load(technical_covariate_input_rdata)

required_objects <- c("response_matrix", "gRNA_data", "meta_data", "covariate_data")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("Missing expected object(s): ", paste(missing_objects, collapse = ", "))
}

message("Loading the curated guide-to-target annotation")
guide_target_annotation <- read.csv(guide_target_annotation_csv, stringsAsFactors = FALSE)
colnames(guide_target_annotation) <- c("grna_id", "grna_target")
guide_target_annotation$grna_id <- gsub("_", "-", guide_target_annotation$grna_id)
guide_target_annotation <- distinct(guide_target_annotation)

# Rename the loaded objects so the script reads in terms of the scientific inputs rather than legacy object names.
gene_by_cell_count_matrix <- response_matrix
guide_by_cell_umi_matrix <- gRNA_data
cell_metadata <- meta_data
technical_covariate_table <- covariate_data

# Step 1: align barcodes across the transcriptome, guide, metadata, and covariate objects.
# Input: saved response, guide, metadata, and covariate tables from preprocessing.
# Why: the same cell set must be used consistently for transcriptome responses, perturbation labels, and model covariates.
# Output: aligned matrices and data frames with a shared barcode order.
message("Aligning shared barcodes across response, guide, metadata, and covariate tables")
common_barcodes <- Reduce(
  intersect,
  list(
    colnames(gene_by_cell_count_matrix),
    colnames(guide_by_cell_umi_matrix),
    rownames(cell_metadata),
    rownames(technical_covariate_table)
  )
)
if (length(common_barcodes) == 0) {
  stop("No shared barcodes across the count, guide, metadata, and covariate objects")
}

gene_by_cell_count_matrix <- gene_by_cell_count_matrix[, common_barcodes, drop = FALSE]
guide_by_cell_umi_matrix <- guide_by_cell_umi_matrix[, common_barcodes, drop = FALSE]
cell_metadata <- cell_metadata[common_barcodes, , drop = FALSE]
technical_covariate_table <- technical_covariate_table[common_barcodes, , drop = FALSE]

# Step 2: prepare the model covariates and subtype labels.
# Input: the aligned metadata and precomputed technical covariates.
# Why: the final transcriptome analysis uses precomputed RNA/guide library covariates plus batch, and cell.type.l2 defines the endothelial subtype strata.
# Output: aligned covariate tables for all cells and for later subtype-specific subsetting.
cell_metadata <- cell_metadata |>
  mutate(
    batch = factor(batch),
    cell.type.l2 = factor(cell.type.l2)
  )

technical_covariate_table <- technical_covariate_table |>
  mutate(batch = factor(batch))

if ("cell.type" %in% colnames(technical_covariate_table)) {
  technical_covariate_table$cell.type <- NULL
}
technical_covariate_table <- drop_unused_factor_levels(technical_covariate_table)

model_extra_covariates <- technical_covariate_table |>
  transmute(
    p_mito = p_mito,
    batch = batch
  )

message(
  "Prepared aligned matrices: ",
  nrow(gene_by_cell_count_matrix), " response genes x ",
  ncol(gene_by_cell_count_matrix), " cells; ",
  nrow(guide_by_cell_umi_matrix), " guides"
)

if (!dry_run) {
  dir.create(transcriptome_sceptre_output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(transcriptome_sceptre_output_dir, "sessionInfo_script_start.txt"))
}

all_cells_payload <- list(
  response_mat = gene_by_cell_count_matrix,
  gRNA_mat = guide_by_cell_umi_matrix,
  covars = technical_covariate_table,
  extra_covars = model_extra_covariates
)

run_sceptre_subset(
  subset_label = "AllCells_EndothelialScreen",
  payload = all_cells_payload,
  guide_map = guide_target_annotation,
  out_dir = file.path(transcriptome_sceptre_output_dir, "AllCells_EndothelialScreen")
)

if (run_per_celltype) {
  cell_type_counts <- sort(table(cell_metadata$cell.type.l2), decreasing = TRUE)
  message("Endothelial subtype counts before per-cell-type SCEPTRE runs:")
  print(cell_type_counts)

  eligible_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells_per_type])
  message(
    "Eligible endothelial subtypes (>= ", min_cells_per_type, " cells): ",
    paste(eligible_cell_types, collapse = ", ")
  )

  for (ct in eligible_cell_types) {
    idx <- cell_metadata$cell.type.l2 == ct

    # Step 3: subset to one endothelial subtype and drop unused factor levels.
    # Input: aligned transcriptome, guide, and covariate objects plus the cell.type.l2 label vector.
    # Why: the study reports both the global screen and subtype-specific transcriptome tests.
    # Output: one subtype-specific payload ready for SCEPTRE, with redundant factor levels removed.
    payload <- list(
      response_mat = gene_by_cell_count_matrix[, idx, drop = FALSE],
      gRNA_mat = guide_by_cell_umi_matrix[, idx, drop = FALSE],
      covars = drop_unused_factor_levels(technical_covariate_table[idx, , drop = FALSE]),
      extra_covars = drop_unused_factor_levels(model_extra_covariates[idx, , drop = FALSE])
    )

    subtype_label <- paste0("CellTypeL2_", ct)
    subtype_dir <- file.path(transcriptome_sceptre_output_dir, paste0("CellTypeL2_", make.names(ct)))

    run_sceptre_subset(
      subset_label = subtype_label,
      payload = payload,
      guide_map = guide_target_annotation,
      out_dir = subtype_dir
    )
  }
}

message("Gene-level transcriptome SCEPTRE analysis complete.")
