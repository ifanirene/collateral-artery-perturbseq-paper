#!/usr/bin/env python3
"""
@theme: T2 - Cross-species program and module differentials
@task_id: T2-02
@description
Step-by-step script for reproducing the cross-species Hallmark pathway
comparison between heart and brain endothelial cells used in Figure 2h.
It is responsible for restoring count-like matrices from processed AnnData
inputs, renormalizing counts, scoring Hallmark gene sets with AUCell,
testing guinea pig-versus-mouse pathway differences within each organ,
merging the organ-level pathway tables, and exporting the stringent
cross-organ comparison figure and summary tables.

@inputs
- `path/to/processed_heart_endothelial_cells.h5ad`: Processed heart endothelial
  AnnData object with `obs["species"]` and raw counts stored in
  `layers["counts"]` or `raw.X`.
- `path/to/processed_brain_endothelial_cells.h5ad`: Processed brain endothelial
  AnnData object with `obs["species"]` and raw counts stored in
  `layers["counts"]` or `raw.X`.
- `path/to/hallmark_gene_sets_mouse_symbols.csv`: Hallmark gene-set table
  with `source`, `target`, and optional `weight` columns.

@outputs
- `path/to/cross_species_endothelial_hallmark_pathway_release/heart_endothelial_pathway_differences_guinea_pig_vs_mouse.csv`
- `path/to/cross_species_endothelial_hallmark_pathway_release/brain_endothelial_pathway_differences_guinea_pig_vs_mouse.csv`
- `path/to/cross_species_endothelial_hallmark_pathway_release/shared_hallmark_pathway_differences_heart_and_brain.csv`
- `path/to/cross_species_endothelial_hallmark_pathway_release/cross_species_endothelial_hallmark_pathway_summary.tsv`
- `path/to/cross_species_endothelial_hallmark_pathway_release/cross_species_endothelial_hallmark_pathway_comparison.pdf`

@key_params
- `species_column = "species"`: Metadata column used for the guinea pig-versus-
  mouse comparison.
- `target_sum = 1e4`: Count normalization target used before AUCell scoring.
- `minimum_genes_per_pathway = 5`: Minimum overlap required for Hallmark scoring.
- `adjusted_p_value_threshold = 0.01`: Stringent threshold used for figure
  highlighting.
- `absolute_log2_fold_change_threshold = 0.25`: Effect-size threshold required
  in both organs for the stringent comparison.

@dependencies
- `Python` (v3.12.9): Script runtime.
- `anndata` (v0.11.4): Reads and rebuilds AnnData objects from count-like layers.
- `scanpy` (v1.11.3): Count normalization, log transform, and Wilcoxon testing.
- `decoupler` (v2.1.1): AUCell Hallmark pathway scoring.
- `pandas` (v2.2.3): Table handling and export.
- `numpy` (v2.2.4): Numeric utilities.
- `scipy` (v1.15.2): Pearson correlation.
- `matplotlib` (v3.10.1): Figure export.
- `seaborn` (v0.13.2): Regression overlay.
- `adjustText` (v1.3.0): Label collision handling.

@examples
python path/to/reproduce_cross_species_heart_brain_endothelial_hallmark_pathway_comparison.py
"""

# Run this script top to bottom in the `perturb2` environment. Each major
# section creates a table or figure that is used in the final cross-organ
# Hallmark pathway comparison. Replace the example paths and example file names
# below with the locations used in your own release package.
#
# Recorded environment snapshot from the perturb2 conda environment on 2026-04-03:
# - Python: 3.12.9
# - anndata: 0.11.4
# - scanpy: 1.11.3
# - decoupler: 2.1.1
# - pandas: 2.2.3
# - numpy: 2.2.4
# - scipy: 1.15.2
# - matplotlib: 3.10.1
# - seaborn: 0.13.2
# - adjustText: 1.3.0

from pathlib import Path
import sys

import anndata as ad
import decoupler as dc
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text
from scipy.stats import pearsonr

# Keep PDF text editable and prefer Arial when the system provides it.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


###############################################################################
# STEP 1. Define the input files, output files, and analysis parameters.
#
# The two processed AnnData objects represent the retained heart and brain
# endothelial comparisons. Both inputs are expected to carry species metadata
# and a recoverable count-like matrix so the pathway scoring can start from raw
# counts rather than a previously normalized expression matrix.
###############################################################################

heart_endothelial_input_path = Path("path/to/processed_heart_endothelial_cells.h5ad")
brain_endothelial_input_path = Path("path/to/processed_brain_endothelial_cells.h5ad")
hallmark_gene_set_table_path = Path("path/to/hallmark_gene_sets_mouse_symbols.csv")

output_directory = Path("path/to/cross_species_endothelial_hallmark_pathway_release")
heart_pathway_output_path = (
    output_directory / "heart_endothelial_pathway_differences_guinea_pig_vs_mouse.csv"
)
brain_pathway_output_path = (
    output_directory / "brain_endothelial_pathway_differences_guinea_pig_vs_mouse.csv"
)
shared_pathway_output_path = (
    output_directory / "shared_hallmark_pathway_differences_heart_and_brain.csv"
)
summary_output_path = output_directory / "cross_species_endothelial_hallmark_pathway_summary.tsv"
figure_output_path = output_directory / "cross_species_endothelial_hallmark_pathway_comparison.pdf"

species_column = "species"
guinea_pig_species_label = "gp"
mouse_species_label = "mm"

target_sum = 1e4
minimum_genes_per_pathway = 5
adjusted_p_value_threshold = 0.01
absolute_log2_fold_change_threshold = 0.25

pathways_excluded_from_stringent_highlight = [
    "Apical Surface",
    "Wnt Beta Catenin Signaling",
    "Epithelial Mesenchymal Transition",
]
pathways_forced_into_figure_labels = ["Hedgehog Signaling"]

dataset_manifest = pd.DataFrame(
    {
        "organ_name": ["heart", "brain"],
        "organ_display_name": ["Heart endothelial cells", "Brain endothelial cells"],
        "processed_input_path": [heart_endothelial_input_path, brain_endothelial_input_path],
        "organ_level_output_path": [heart_pathway_output_path, brain_pathway_output_path],
    }
)

output_directory.mkdir(parents=True, exist_ok=True)

missing_input_files = [
    str(path)
    for path in [heart_endothelial_input_path, brain_endothelial_input_path, hallmark_gene_set_table_path]
    if not path.exists()
]
if missing_input_files:
    raise FileNotFoundError(
        "The following input files are missing. Replace the example paths in the "
        "configuration block before running the script:\n- "
        + "\n- ".join(missing_input_files)
    )

print("Running the cross-species heart and brain endothelial Hallmark pathway comparison")
print(f"Python version: {sys.version.split()[0]}")
print(f"scanpy version: {sc.__version__}")
print(f"decoupler version: {dc.__version__}")
print(f"Heart endothelial input: {heart_endothelial_input_path}")
print(f"Brain endothelial input: {brain_endothelial_input_path}")
print(f"Hallmark gene-set table: {hallmark_gene_set_table_path}")
print(f"Output directory: {output_directory}")


###############################################################################
# STEP 2. Load the Hallmark gene-set table and confirm that the required
# columns are present.
#
# The retained workflow uses a long-format gene-set table with one row per
# pathway-gene pair. The `source` column stores the Hallmark pathway name and
# the `target` column stores the gene symbol used for scoring.
###############################################################################

hallmark_gene_set_table = pd.read_csv(hallmark_gene_set_table_path)
required_gene_set_columns = {"source", "target"}
missing_gene_set_columns = required_gene_set_columns.difference(hallmark_gene_set_table.columns)
if missing_gene_set_columns:
    raise ValueError(
        "The Hallmark gene-set table is missing required columns: "
        + ", ".join(sorted(missing_gene_set_columns))
    )

print(
    "Loaded "
    f"{hallmark_gene_set_table['source'].nunique()} Hallmark pathways from "
    f"{hallmark_gene_set_table_path.name}"
)


###############################################################################
# STEP 3. Define two small utility functions used for both organs.
#
# These functions keep the heart and brain branches synchronized: both recover
# counts from the same supported locations and both extract AUCell scores with
# the same fallback logic in case the storage key differs across decoupler
# versions.
###############################################################################


def select_count_like_matrix(processed_endothelial_data: ad.AnnData) -> tuple[object, pd.DataFrame, str]:
    """Return the retained count-like matrix, its gene table, and the source label."""
    if "counts" in processed_endothelial_data.layers:
        return (
            processed_endothelial_data.layers["counts"].copy(),
            processed_endothelial_data.var.copy(),
            "layers['counts']",
        )

    if processed_endothelial_data.raw is not None:
        raw_gene_table = processed_endothelial_data.raw.var.copy()
        sampled_raw_names = raw_gene_table.index.astype(str)[: min(10, len(raw_gene_table.index))]
        if (
            processed_endothelial_data.raw.shape[1] == processed_endothelial_data.n_vars
            and all(gene_name.isdigit() for gene_name in sampled_raw_names)
        ):
            raw_gene_table = processed_endothelial_data.var.copy()
        return processed_endothelial_data.raw.X.copy(), raw_gene_table, "raw.X"

    raise ValueError(
        "No count-like matrix was found. The processed AnnData object must store "
        "raw counts in `layers[\"counts\"]` or `raw.X`."
    )


def extract_aucell_score_table(
    scored_endothelial_data: ad.AnnData,
    hallmark_gene_set_table: pd.DataFrame,
) -> pd.DataFrame:
    """Return AUCell pathway scores as a cell-by-pathway DataFrame."""
    if "score_aucell" in scored_endothelial_data.obsm:
        aucell_scores = scored_endothelial_data.obsm["score_aucell"]
    else:
        aucell_score_keys = [
            key for key in scored_endothelial_data.obsm.keys() if "aucell" in key.lower()
        ]
        if not aucell_score_keys:
            raise ValueError(
                "AUCell scores were not found in `obsm`. Confirm that decoupler "
                "completed successfully before pathway extraction."
            )
        aucell_scores = scored_endothelial_data.obsm[aucell_score_keys[0]]

    if isinstance(aucell_scores, pd.DataFrame):
        return aucell_scores.copy()

    hallmark_pathway_names = pd.Index(
        hallmark_gene_set_table["source"].drop_duplicates().tolist(),
        name="pathway_name",
    )
    return pd.DataFrame(
        aucell_scores,
        index=scored_endothelial_data.obs_names.copy(),
        columns=hallmark_pathway_names,
    )


###############################################################################
# STEP 4. Process heart and brain endothelial cells with the same retained
# single-cell pathway workflow.
#
# For each organ, the script restores a count-like matrix, renormalizes counts,
# scores Hallmark pathways with AUCell, and then uses Scanpy's Wilcoxon test to
# compare guinea pig versus mouse pathway activity.
###############################################################################

organ_level_results = {}
organ_level_metadata = {}

for dataset_row in dataset_manifest.itertuples(index=False):
    organ_name = dataset_row.organ_name
    organ_display_name = dataset_row.organ_display_name
    processed_input_path = Path(dataset_row.processed_input_path)
    organ_level_output_path = Path(dataset_row.organ_level_output_path)

    print(f"\nProcessing {organ_display_name} from {processed_input_path.name}")
    processed_endothelial_data = sc.read_h5ad(processed_input_path)

    if species_column not in processed_endothelial_data.obs.columns:
        raise KeyError(
            f"{processed_input_path.name} does not contain the required "
            f"`obs[{species_column!r}]` metadata column."
        )

    count_like_matrix, count_like_gene_table, count_like_matrix_source = select_count_like_matrix(
        processed_endothelial_data
    )

    working_endothelial_data = ad.AnnData(
        X=count_like_matrix,
        obs=processed_endothelial_data.obs.copy(),
        var=count_like_gene_table.copy(),
    )
    working_endothelial_data.obs_names = processed_endothelial_data.obs_names.copy()
    working_endothelial_data.var_names = count_like_gene_table.index.astype(str)

    sc.pp.normalize_total(working_endothelial_data, target_sum=target_sum)
    sc.pp.log1p(working_endothelial_data)
    dc.mt.aucell(
        working_endothelial_data,
        hallmark_gene_set_table,
        tmin=minimum_genes_per_pathway,
        verbose=True,
    )

    pathway_activity_table = extract_aucell_score_table(
        working_endothelial_data,
        hallmark_gene_set_table,
    )
    pathway_activity_table.columns = (
        pd.Index(pathway_activity_table.columns)
        .astype(str)
        .str.replace("HALLMARK_", "", regex=False)
        .str.replace("_", " ", regex=False)
        .str.title()
    )

    pathway_activity_data = ad.AnnData(
        X=pathway_activity_table.to_numpy(),
        obs=working_endothelial_data.obs.loc[pathway_activity_table.index].copy(),
        var=pd.DataFrame(index=pathway_activity_table.columns),
    )

    sc.tl.rank_genes_groups(
        pathway_activity_data,
        species_column,
        method="wilcoxon",
        groups=[guinea_pig_species_label],
        reference=mouse_species_label,
        key_added="guinea_pig_vs_mouse_pathway_difference",
    )

    organ_pathway_results = sc.get.rank_genes_groups_df(
        pathway_activity_data,
        group=guinea_pig_species_label,
        key="guinea_pig_vs_mouse_pathway_difference",
    )
    organ_pathway_results = organ_pathway_results.loc[
        :, ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    ].rename(
        columns={
            "names": "pathway_name",
            "scores": "wilcoxon_score",
            "logfoldchanges": "log2_fold_change_guinea_pig_vs_mouse",
            "pvals": "p_value",
            "pvals_adj": "adjusted_p_value",
        }
    )

    organ_pathway_results.to_csv(organ_level_output_path, index=False)
    organ_level_results[organ_name] = organ_pathway_results
    organ_level_metadata[organ_name] = {
        "processed_input_path": str(processed_input_path),
        "count_like_matrix_source": count_like_matrix_source,
        "species_cell_counts": (
            pathway_activity_data.obs[species_column].astype(str).value_counts().to_dict()
        ),
        "cell_count": int(pathway_activity_data.n_obs),
        "pathway_count": int(pathway_activity_data.n_vars),
    }

    print(f"Saved {organ_display_name.lower()} pathway table to {organ_level_output_path.name}")


###############################################################################
# STEP 5. Merge the heart and brain pathway tables and define the stringent
# highlighted pathway set used in the final figure.
#
# A pathway is highlighted only if it passes the adjusted P-value and effect-
# size thresholds in both organs, after which the retained manual exclusions are
# applied so the shared figure emphasizes the pathways kept in the paper.
###############################################################################

shared_pathway_results = organ_level_results["heart"].merge(
    organ_level_results["brain"],
    on="pathway_name",
    suffixes=("_heart", "_brain"),
)

heart_passes_stringent_threshold = (
    shared_pathway_results["adjusted_p_value_heart"] < adjusted_p_value_threshold
) & (
    shared_pathway_results["log2_fold_change_guinea_pig_vs_mouse_heart"].abs()
    >= absolute_log2_fold_change_threshold
)
brain_passes_stringent_threshold = (
    shared_pathway_results["adjusted_p_value_brain"] < adjusted_p_value_threshold
) & (
    shared_pathway_results["log2_fold_change_guinea_pig_vs_mouse_brain"].abs()
    >= absolute_log2_fold_change_threshold
)

passes_both_organs = heart_passes_stringent_threshold & brain_passes_stringent_threshold
highlight_pathway_mask = passes_both_organs & (
    ~shared_pathway_results["pathway_name"].isin(pathways_excluded_from_stringent_highlight)
)

highlighted_pathway_results = shared_pathway_results.loc[highlight_pathway_mask].copy()
threshold_only_excluded_pathway_results = shared_pathway_results.loc[
    passes_both_organs
    & shared_pathway_results["pathway_name"].isin(pathways_excluded_from_stringent_highlight)
].copy()

shared_pathway_results.to_csv(shared_pathway_output_path, index=False)
print(f"\nSaved merged heart-and-brain pathway table to {shared_pathway_output_path.name}")


###############################################################################
# STEP 6. Create the cross-organ Hallmark pathway comparison figure.
#
# The figure places heart and brain pathway log2 fold changes on shared axes,
# shows the highlighted concordant pathways in red, and overlays the Pearson
# correlation calculated from the final highlighted pathway set.
###############################################################################

if len(highlighted_pathway_results) >= 2:
    pearson_correlation, pearson_p_value = pearsonr(
        highlighted_pathway_results["log2_fold_change_guinea_pig_vs_mouse_heart"],
        highlighted_pathway_results["log2_fold_change_guinea_pig_vs_mouse_brain"],
    )
    correlation_input = "highlighted_pathways"
else:
    pearson_correlation, pearson_p_value = pearsonr(
        shared_pathway_results["log2_fold_change_guinea_pig_vs_mouse_heart"],
        shared_pathway_results["log2_fold_change_guinea_pig_vs_mouse_brain"],
    )
    correlation_input = "all_pathways"

plt.figure(figsize=(10, 10))
plt.style.use("seaborn-v0_8-paper")

plt.scatter(
    shared_pathway_results.loc[
        ~highlight_pathway_mask, "log2_fold_change_guinea_pig_vs_mouse_heart"
    ],
    shared_pathway_results.loc[
        ~highlight_pathway_mask, "log2_fold_change_guinea_pig_vs_mouse_brain"
    ],
    color="lightgrey",
    alpha=0.5,
    s=40,
    label="Not highlighted",
)

plt.scatter(
    highlighted_pathway_results["log2_fold_change_guinea_pig_vs_mouse_heart"],
    highlighted_pathway_results["log2_fold_change_guinea_pig_vs_mouse_brain"],
    color="#d62728",
    alpha=0.8,
    s=60,
    label=(
        "Highlighted in both organs "
        f"(adjusted P value < {adjusted_p_value_threshold}, "
        f"absolute log2 fold change >= {absolute_log2_fold_change_threshold})"
    ),
)

forced_label_pathway_results = shared_pathway_results.loc[
    shared_pathway_results["pathway_name"].isin(pathways_forced_into_figure_labels)
].copy()
forced_label_pathway_results = forced_label_pathway_results.loc[
    ~forced_label_pathway_results["pathway_name"].isin(highlighted_pathway_results["pathway_name"])
]

if not forced_label_pathway_results.empty:
    plt.scatter(
        forced_label_pathway_results["log2_fold_change_guinea_pig_vs_mouse_heart"],
        forced_label_pathway_results["log2_fold_change_guinea_pig_vs_mouse_brain"],
        color="#d62728",
        alpha=0.8,
        s=60,
        zorder=4,
    )

label_table = pd.concat(
    [highlighted_pathway_results, forced_label_pathway_results],
    ignore_index=True,
)
label_text_objects = [
    plt.text(
        row["log2_fold_change_guinea_pig_vs_mouse_heart"],
        row["log2_fold_change_guinea_pig_vs_mouse_brain"],
        row["pathway_name"],
        fontsize=9,
        fontweight="bold",
    )
    for _, row in label_table.iterrows()
]
if label_text_objects:
    adjust_text(
        label_text_objects,
        arrowprops={"arrowstyle": "->", "color": "grey", "alpha": 0.5},
    )

if len(highlighted_pathway_results) >= 2:
    sns.regplot(
        data=highlighted_pathway_results,
        x="log2_fold_change_guinea_pig_vs_mouse_heart",
        y="log2_fold_change_guinea_pig_vs_mouse_brain",
        scatter=False,
        color="#d62728",
        line_kws={"linestyle": "--", "alpha": 0.5},
    )

plt.axhline(0, color="black", linestyle="--", alpha=0.3)
plt.axvline(0, color="black", linestyle="--", alpha=0.3)
plt.xlabel("Heart endothelial cells log2 fold change (guinea pig versus mouse)", fontsize=12)
plt.ylabel("Brain endothelial cells log2 fold change (guinea pig versus mouse)", fontsize=12)
plt.title(
    "Cross-species Hallmark pathway comparison across heart and brain endothelial cells",
    fontsize=14,
)
plt.text(
    0.05,
    0.95,
    f"Pearson r = {pearson_correlation:.3f}\n"
    f"P value = {pearson_p_value:.2e}\n"
    f"Highlighted pathways = {len(highlighted_pathway_results)}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.7},
)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(figure_output_path, bbox_inches="tight")
plt.close()

print(f"Saved cross-organ Hallmark pathway figure to {figure_output_path.name}")


###############################################################################
# STEP 7. Write a compact run summary for the release package.
#
# The summary records the organ-level input files, where the raw counts were
# restored from, the number of cells per species, how many pathways entered the
# merged comparison, and the exact thresholds used for the final highlighted set.
###############################################################################

summary_rows = [
    {"metric": "heart_input_file", "value": str(heart_endothelial_input_path)},
    {"metric": "brain_input_file", "value": str(brain_endothelial_input_path)},
    {"metric": "hallmark_gene_set_table", "value": str(hallmark_gene_set_table_path)},
    {
        "metric": "heart_count_like_matrix_source",
        "value": organ_level_metadata["heart"]["count_like_matrix_source"],
    },
    {
        "metric": "brain_count_like_matrix_source",
        "value": organ_level_metadata["brain"]["count_like_matrix_source"],
    },
    {
        "metric": "heart_species_cell_counts",
        "value": str(organ_level_metadata["heart"]["species_cell_counts"]),
    },
    {
        "metric": "brain_species_cell_counts",
        "value": str(organ_level_metadata["brain"]["species_cell_counts"]),
    },
    {
        "metric": "shared_hallmark_pathway_count",
        "value": str(len(shared_pathway_results)),
    },
    {
        "metric": "pathways_passing_thresholds_in_both_organs_before_manual_exclusions",
        "value": str(int(passes_both_organs.sum())),
    },
    {
        "metric": "highlighted_pathway_count_after_manual_exclusions",
        "value": str(len(highlighted_pathway_results)),
    },
    {
        "metric": "excluded_pathways_still_passing_thresholds",
        "value": ",".join(threshold_only_excluded_pathway_results["pathway_name"].tolist()),
    },
    {
        "metric": "pathways_forced_into_figure_labels",
        "value": ",".join(pathways_forced_into_figure_labels),
    },
    {
        "metric": "adjusted_p_value_threshold",
        "value": str(adjusted_p_value_threshold),
    },
    {
        "metric": "absolute_log2_fold_change_threshold",
        "value": str(absolute_log2_fold_change_threshold),
    },
    {"metric": "pearson_correlation", "value": f"{pearson_correlation:.6f}"},
    {"metric": "pearson_p_value", "value": f"{pearson_p_value:.6e}"},
    {"metric": "correlation_input", "value": correlation_input},
    {"metric": "figure_output_file", "value": str(figure_output_path)},
]

pd.DataFrame(summary_rows).to_csv(summary_output_path, sep="\t", index=False)
print(f"Saved run summary to {summary_output_path.name}")
print("\nDone. The release directory now contains organ-level pathway tables,")
print("the merged heart-and-brain pathway table, a compact run summary,")
print("and the final cross-organ Hallmark pathway comparison figure.")
