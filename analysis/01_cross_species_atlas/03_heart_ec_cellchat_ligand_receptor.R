#!/usr/bin/env Rscript

# @theme: T3 - CellChat ligand-receptor communication
# @task_id: T3-01
# @description: Step-by-step script for reproducing the CellChat comparison of mouse versus guinea pig heart signaling into endothelial cells, including global pathway summaries and EC-targeted differential interaction ranking.
# @inputs:
# - `path/to/heart_integrated_mm_gp_seurat.rds`: example integrated heart Seurat object with mouse and guinea pig cells plus `species`, `sample` or `samples`, and `cell.type.l2` metadata.
# @outputs:
# - `path/to/heart_ec_cellchat_release/rds/`: example release directory for species-specific and merged CellChat objects.
# - `path/to/heart_ec_cellchat_release/tables/`: example release directory for EC-receiver pathway summaries, EC-targeted interaction tables, and receptor-ranking tables.
# - `path/to/heart_ec_cellchat_release/plots/`: example release directory for overview plots.
# @key_params:
# - `target_celltype = "EC"`: receiver population used throughout the comparison.
# - `min_cells_per_group = 3`: relaxed CellChat cell-count filter used to retain low-frequency developmental cell groups.
# - `positive_dataset = "gp"`: dataset treated as the positive comparison for CellChat differential mapping.
# - `gp_ligand_logfc_min = 0.01` and `gp_receptor_logfc_min = 0.25`: thresholds used to label gp-enriched EC-targeted interactions.
# - `mm_ligand_logfc_max = -0.01` and `mm_receptor_logfc_max = -0.50`: thresholds used to label mm-enriched EC-targeted interactions.
# @dependencies:
# - `CellChat` (v2.1.2): ligand-receptor communication inference and differential communication mapping.
# - `Seurat` (v4.4.0): integrated heart atlas input and metadata handling.
# - `dplyr` (v1.1.4): table manipulation.
# - `tidyr` (v1.3.1): reshaping of summary tables.
# - `ggplot2` (v3.5.0): overview plots.
# @examples:
# Rscript path/to/reproduce_heart_ec_cellchat_comparison.R
#
# Run this script block by block in the Seurat 4 / CellChat environment. The
# workflow analyzes mouse-versus-guinea-pig heart signaling into endothelial
# cells and writes release-ready CellChat objects, tables, and overview plots.
#
# Recorded environment snapshot from the seurat4 conda environment on 2026-03-30:
# - R: 4.2.0
# - platform: x86_64-conda-linux-gnu (64-bit)
# - running: CentOS Linux 7 (Core)
# - CellChat: 2.1.2
# - Seurat: 4.4.0
# - dplyr: 1.1.4
# - tidyr: 1.3.1
# - ggplot2: 3.5.0

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

###############################################################################
# STEP 1. Define the input object, output locations, and analysis settings.
#
# The release script starts from the integrated heart Seurat object. Replace the
# example paths below with the locations bundled in your own release package.
###############################################################################

input_seurat_rds <- "Processed_data/heart.integrated.050524.rds"
output_dir <- "path/to/heart_ec_cellchat_release"

target_celltype <- "EC"
positive_dataset <- "gp"
negative_dataset <- "mm"

pathway_plot_top_n <- 20
receptor_plot_top_n <- 20
min_cells_per_group <- 3

gp_ligand_logfc_min <- 0.01
gp_receptor_logfc_min <- 0.25
mm_ligand_logfc_max <- -0.01
mm_receptor_logfc_max <- -0.50

canonical_pathways <- c("VEGF", "NOTCH", "CXCL", "WNT", "ncWNT", "TGFb")

rds_dir <- file.path(output_dir, "rds")
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_seurat_rds)) {
  stop(sprintf("Integrated heart Seurat object not found: %s", input_seurat_rds), call. = FALSE)
}

cat("CellChat release workflow for heart EC incoming signaling\n")
cat(sprintf("Input Seurat object: %s\n", input_seurat_rds))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Target receiver cell type: %s\n", target_celltype))
cat(sprintf("Differential comparison: %s vs %s\n\n", positive_dataset, negative_dataset))

cat("Live session info at script start:\n")
print(sessionInfo())

###############################################################################
# STEP 2. Load the integrated heart atlas, keep the mouse-versus-guinea-pig
# comparison, and confirm the metadata columns that define species, sample, and
# annotated cardiac cell classes.
###############################################################################

choose_metadata_col <- function(meta_df, candidates, label) {
  hits <- candidates[candidates %in% colnames(meta_df)]
  if (length(hits) == 0) {
    stop(sprintf("No `%s` metadata column found. Tried: %s", label, paste(candidates, collapse = ", ")), call. = FALSE)
  }
  hits[[1]]
}

heart <- readRDS(input_seurat_rds)
DefaultAssay(heart) <- "RNA"

heart_meta <- heart@meta.data
heart_meta <- heart_meta[, !grepl("^(pANN|DF)", colnames(heart_meta)), drop = FALSE]
heart@meta.data <- heart_meta

species_col <- choose_metadata_col(heart@meta.data, c("species", "datasets"), "species")
sample_col <- choose_metadata_col(heart@meta.data, c("samples", "sample"), "sample")
celltype_col <- choose_metadata_col(heart@meta.data, c("cell.type.l2", "celltype", "cell_type"), "cell type")

heart <- subset(
  heart,
  cells = rownames(heart@meta.data)[as.character(heart@meta.data[[species_col]]) %in% c(positive_dataset, negative_dataset)]
)

heart@meta.data[[species_col]] <- factor(as.character(heart@meta.data[[species_col]]), levels = c(negative_dataset, positive_dataset))
heart@meta.data[[sample_col]] <- factor(as.character(heart@meta.data[[sample_col]]))
heart@meta.data[[celltype_col]] <- factor(as.character(heart@meta.data[[celltype_col]]))
Idents(heart) <- species_col

cat(sprintf("Retained %d cells after limiting the analysis to %s and %s.\n", ncol(heart), negative_dataset, positive_dataset))
cat(sprintf("Using metadata columns: species=%s | sample=%s | celltype=%s\n\n", species_col, sample_col, celltype_col))

overview_summary <- bind_rows(
  tibble(section = "overview", variable = "cells_total", level = "all", count = ncol(heart), pct = 100),
  as.data.frame(table(heart@meta.data[[species_col]]), stringsAsFactors = FALSE) %>%
    transmute(
      section = "composition",
      variable = "species",
      level = as.character(Var1),
      count = as.integer(Freq),
      pct = round(100 * count / sum(count), 3)
    ),
  as.data.frame(table(heart@meta.data[[sample_col]]), stringsAsFactors = FALSE) %>%
    transmute(
      section = "composition",
      variable = "sample",
      level = as.character(Var1),
      count = as.integer(Freq),
      pct = round(100 * count / sum(count), 3)
    ),
  as.data.frame(table(heart@meta.data[[celltype_col]]), stringsAsFactors = FALSE) %>%
    transmute(
      section = "composition",
      variable = "celltype",
      level = as.character(Var1),
      count = as.integer(Freq),
      pct = round(100 * count / sum(count), 3)
    )
)

write.table(
  overview_summary,
  file = file.path(table_dir, "cellchat_input_summary_mm_gp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 3. Build one relaxed CellChat object per species with the annotated heart
# cell classes as sender/receiver groups.
#
# The relaxed `min.cells` filter is kept because several developmental sender
# groups are small but biologically important.
###############################################################################

build_species_cellchat <- function(seurat_obj, species_id, group_col, min_cells) {
  species_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj@meta.data[[species_col]]) == species_id]
  species_obj <- subset(seurat_obj, cells = species_cells)
  species_obj@meta.data[[group_col]] <- droplevels(factor(as.character(species_obj@meta.data[[group_col]])))

  cellchat_obj <- createCellChat(object = species_obj, group.by = group_col, assay = "RNA")
  cellchat_obj@DB <- CellChatDB.mouse

  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- projectData(cellchat_obj, PPI.mouse)
  cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj, thresh = 1.0)
  cellchat_obj <- aggregateNet(cellchat_obj, thresh = 1.0)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")

  cellchat_obj
}

all_celltypes <- levels(heart@meta.data[[celltype_col]])
species_ids <- c(negative_dataset, positive_dataset)

cellchat_by_species <- list()
for (species_id in species_ids) {
  cat(sprintf("Running CellChat for %s\n", species_id))
  cellchat_by_species[[species_id]] <- build_species_cellchat(
    seurat_obj = heart,
    species_id = species_id,
    group_col = celltype_col,
    min_cells = min_cells_per_group
  )
  cellchat_by_species[[species_id]] <- liftCellChat(cellchat_by_species[[species_id]], all_celltypes)
  saveRDS(
    cellchat_by_species[[species_id]],
    file = file.path(rds_dir, sprintf("cellchat_%s_relaxed.rds", species_id))
  )
}

cellchat_merged <- mergeCellChat(cellchat_by_species, add.names = names(cellchat_by_species))
saveRDS(cellchat_merged, file = file.path(rds_dir, "cellchat_merged_relaxed.rds"))

###############################################################################
# STEP 4. Summarize all incoming signaling to endothelial cells across pathways.
#
# This step asks whether whole-pathway input into ECs is globally different between
# guinea pig versus mouse. We therefore aggregate all source probabilities
# directed toward EC rather than focusing only on a few predefined
# ligand-receptor pairs.
###############################################################################

extract_pathway_source_probabilities <- function(cellchat_obj, species_id, target_cell) {
  pathway_array <- cellchat_obj@netP$prob
  target_levels <- dimnames(pathway_array)[[2]]
  if (!(target_cell %in% target_levels)) {
    stop(sprintf("Target cell type `%s` is absent from the `%s` CellChat object.", target_cell, species_id), call. = FALSE)
  }

  target_index <- match(target_cell, target_levels)
  sources <- dimnames(pathway_array)[[1]]
  pathways <- dimnames(pathway_array)[[3]]

  probability_table <- expand.grid(
    source = sources,
    pathway_name = pathways,
    stringsAsFactors = FALSE
  )
  probability_table$prob <- as.numeric(pathway_array[, target_index, pathways, drop = FALSE])
  probability_table$species <- species_id
  probability_table$target <- target_cell
  probability_table$pathway_class <- ifelse(probability_table$pathway_name %in% canonical_pathways, "canonical_arterial_pathway", "other")
  probability_table
}

incoming_pathway_long <- bind_rows(lapply(names(cellchat_by_species), function(species_id) {
  extract_pathway_source_probabilities(cellchat_by_species[[species_id]], species_id, target_celltype)
}))

incoming_pathway_summary <- incoming_pathway_long %>%
  group_by(species, target, pathway_name, pathway_class) %>%
  summarise(
    total_incoming_prob = sum(prob, na.rm = TRUE),
    mean_source_prob = mean(prob, na.rm = TRUE),
    max_source_prob = max(prob, na.rm = TRUE),
    top_source = source[which.max(prob)],
    n_active_senders = sum(prob > 0, na.rm = TRUE),
    .groups = "drop"
  )

incoming_pathway_diff <- incoming_pathway_summary %>%
  select(species, pathway_name, pathway_class, total_incoming_prob, mean_source_prob, max_source_prob, top_source, n_active_senders) %>%
  pivot_wider(
    names_from = species,
    values_from = c(total_incoming_prob, mean_source_prob, max_source_prob, top_source, n_active_senders),
    names_sep = "_"
  ) %>%
  mutate(
    gp_minus_mm_total_incoming = total_incoming_prob_gp - total_incoming_prob_mm,
    gp_over_mm_ratio = ifelse(total_incoming_prob_mm == 0, NA_real_, total_incoming_prob_gp / total_incoming_prob_mm)
  ) %>%
  arrange(desc(pmax(total_incoming_prob_gp, total_incoming_prob_mm, na.rm = TRUE)))

write.table(
  incoming_pathway_summary,
  file = file.path(table_dir, "ec_incoming_pathway_summary_by_species.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  incoming_pathway_diff,
  file = file.path(table_dir, "ec_incoming_pathway_totals_mm_vs_gp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pathway_levels_for_plot <- incoming_pathway_diff %>%
  slice_head(n = min(pathway_plot_top_n, n())) %>%
  pull(pathway_name)

top_pathways_for_plot <- incoming_pathway_diff %>%
  filter(pathway_name %in% pathway_levels_for_plot) %>%
  select(pathway_name, total_incoming_prob_mm, total_incoming_prob_gp) %>%
  pivot_longer(
    cols = c(total_incoming_prob_mm, total_incoming_prob_gp),
    names_to = "species",
    values_to = "total_incoming_prob"
  ) %>%
  mutate(
    species = factor(gsub("^total_incoming_prob_", "", species), levels = c("mm", "gp")),
    pathway_name = factor(pathway_name, levels = rev(pathway_levels_for_plot))
  )

top_pathway_plot <- ggplot(top_pathways_for_plot, aes(x = total_incoming_prob, y = pathway_name, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c(mm = "#4E79A7", gp = "#F28E2B")) +
  labs(
    title = "Incoming CellChat pathway strength to heart ECs",
    subtitle = "All inferred pathways ranked by total sender-to-EC probability",
    x = "Total incoming CellChat probability",
    y = "Pathway",
    fill = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(
  filename = file.path(plot_dir, "ec_incoming_pathway_totals_mm_gp.pdf"),
  plot = top_pathway_plot,
  width = 9,
  height = 7
)

###############################################################################
# STEP 5. Map cross-dataset differential expression onto the merged CellChat
# network and extract EC-targeted differential interactions.
#
# This step keeps the output transparent: we save the full EC-targeted
# differential table and then apply explicit thresholds to define gp-enriched
# and mm-enriched interaction subsets.
###############################################################################

differential_group_col <- choose_metadata_col(cellchat_merged@meta, c("datasets", "species"), "CellChat dataset group")

cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged,
  group.dataset = differential_group_col,
  pos.dataset = positive_dataset,
  features.name = paste0(positive_dataset, ".merged"),
  only.pos = FALSE,
  thresh.pc = 0,
  thresh.fc = 0,
  thresh.p = 1,
  group.DE.combined = FALSE,
  do.fast = FALSE
)

incoming_diff_network <- netMappingDEG(
  cellchat_merged,
  features.name = paste0(positive_dataset, ".merged"),
  variable.all = TRUE
)

incoming_diff_to_ec <- incoming_diff_network %>%
  filter(target == target_celltype) %>%
  mutate(
    pathway_class = ifelse(pathway_name %in% canonical_pathways, "canonical_arterial_pathway", "other")
  ) %>%
  arrange(datasets, pathway_name, desc(prob))

gp_enriched_to_ec <- incoming_diff_to_ec %>%
  filter(
    datasets == positive_dataset,
    ligand.logFC >= gp_ligand_logfc_min,
    receptor.logFC >= gp_receptor_logfc_min
  ) %>%
  arrange(desc(prob))

mm_enriched_to_ec <- incoming_diff_to_ec %>%
  filter(
    datasets == negative_dataset,
    ligand.logFC <= mm_ligand_logfc_max,
    receptor.logFC <= mm_receptor_logfc_max
  ) %>%
  arrange(desc(prob))

write.table(
  incoming_diff_to_ec,
  file = file.path(table_dir, "ec_incoming_interactions_differential_full.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gp_enriched_to_ec,
  file = file.path(table_dir, "ec_incoming_interactions_gp_higher.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mm_enriched_to_ec,
  file = file.path(table_dir, "ec_incoming_interactions_mm_higher.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# STEP 6. Collapse the EC-targeted differential interactions into receptor-level
# ranking tables for candidate nomination review.
#
# The ranking is intentionally transparent: each receptor is represented by the
# strongest EC-targeted interaction that passed the retained differential filter
# for that species-direction subset.
###############################################################################

rank_receptors <- function(interaction_df, dataset_id) {
  if (nrow(interaction_df) == 0) {
    return(tibble(
      dataset = character(),
      receptor = character(),
      pathway_name = character(),
      top_interaction_name = character(),
      top_interaction_name_2 = character(),
      top_ligand = character(),
      max_prob_to_ec = numeric(),
      receptor_logFC = numeric(),
      receptor_pvalue = numeric(),
      ligand_logFC = numeric(),
      ligand_pvalue = numeric(),
      n_supporting_interactions = integer()
    ))
  }

  interaction_df %>%
    group_by(receptor) %>%
    arrange(desc(prob), pathway_name, interaction_name_2, .by_group = TRUE) %>%
    summarise(
      dataset = dataset_id,
      pathway_name = pathway_name[[1]],
      top_interaction_name = interaction_name[[1]],
      top_interaction_name_2 = interaction_name_2[[1]],
      top_ligand = ligand[[1]],
      max_prob_to_ec = prob[[1]],
      receptor_logFC = receptor.logFC[[1]],
      receptor_pvalue = receptor.pvalues[[1]],
      ligand_logFC = ligand.logFC[[1]],
      ligand_pvalue = ligand.pvalues[[1]],
      n_supporting_interactions = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(max_prob_to_ec), receptor)
}

receptor_rank_gp <- rank_receptors(gp_enriched_to_ec, positive_dataset)
receptor_rank_mm <- rank_receptors(mm_enriched_to_ec, negative_dataset)
receptor_rank_all <- bind_rows(receptor_rank_mm, receptor_rank_gp)

write.table(
  receptor_rank_all,
  file = file.path(table_dir, "ec_incoming_receptor_candidate_ranking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_receptors_for_plot <- receptor_rank_all %>%
  group_by(dataset) %>%
  slice_head(n = min(receptor_plot_top_n, n())) %>%
  ungroup() %>%
  mutate(
    receptor_label = paste(receptor, pathway_name, sep = " | "),
    receptor_label = factor(receptor_label, levels = rev(unique(receptor_label)))
  )

if (nrow(top_receptors_for_plot) > 0) {
  receptor_plot <- ggplot(
    top_receptors_for_plot,
    aes(x = max_prob_to_ec, y = receptor_label, color = dataset, size = -log10(pmax(receptor_pvalue, 1e-300)))
  ) +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = c(mm = "#4E79A7", gp = "#F28E2B")) +
    labs(
      title = "Top EC-targeted differential receptors from the heart CellChat comparison",
      subtitle = "Each receptor is represented by its strongest retained incoming interaction to ECs",
      x = "Maximum CellChat probability into ECs",
      y = "Receptor | pathway",
      color = "Dataset",
      size = "-log10(receptor P)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    filename = file.path(plot_dir, "ec_incoming_receptor_candidates_mm_gp.pdf"),
    plot = receptor_plot,
    width = 9,
    height = 7
  )
}

###############################################################################
# STEP 7. Save a short manifest describing the main release objects produced by
# this script.
###############################################################################

release_manifest <- tibble(
  file = c(
    file.path(rds_dir, "cellchat_mm_relaxed.rds"),
    file.path(rds_dir, "cellchat_gp_relaxed.rds"),
    file.path(rds_dir, "cellchat_merged_relaxed.rds"),
    file.path(table_dir, "ec_incoming_pathway_totals_mm_vs_gp.tsv"),
    file.path(table_dir, "ec_incoming_interactions_differential_full.tsv"),
    file.path(table_dir, "ec_incoming_receptor_candidate_ranking.tsv"),
    file.path(plot_dir, "ec_incoming_pathway_totals_mm_gp.pdf")
  ),
  description = c(
    "Mouse-only relaxed CellChat object for the heart atlas",
    "Guinea pig-only relaxed CellChat object for the heart atlas",
    "Merged mouse-versus-guinea-pig relaxed CellChat object",
    "Pathway-level EC incoming signaling totals and gp-versus-mm differences",
    "Full EC-targeted differential interaction table from netMappingDEG",
    "Collapsed receptor-level ranking table for candidate review",
    "Overview plot of pathway-level EC incoming signaling totals"
  )
)

write.table(
  release_manifest,
  file = file.path(output_dir, "release_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCellChat release workflow complete.\n")
cat(sprintf("Objects and tables were written under: %s\n", output_dir))
