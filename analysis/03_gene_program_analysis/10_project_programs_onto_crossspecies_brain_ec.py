#!/usr/bin/env python3
"""Project Perturb-seq cNMF programs onto cross-species brain endothelial cells using STARcat.

Overview
--------
This script projects the 100 gene expression programs (GEPs) discovered by
consensus NMF (cNMF) in the mouse neonatal brain Perturb-seq endothelial
dataset onto an independently collected, integrated cross-species brain
endothelial cell (EC) dataset comprising mouse (E18, P0) and guinea pig (GD35)
brain ECs.

The projection uses STARcat v1.0.9 (Kotliar and Curtis, immunogenomics/starCAT)
with the same Perturb-seq cNMF reference used for the hypoxia/normoxia analysis.
Because the query includes guinea pig cells, whose gene symbols were mapped to
mouse orthologs during integration, only the subset of reference genes present
in the shared ortholog space is used for fitting (6,704 of 10,000 reference
genes).  Per-cell program usage scores are estimated by NNLS with spectra held
fixed and row-normalized to sum to 1 per cell.  Species annotation is retained
in the output to enable species-stratified comparison of program activity.

Prerequisites
-------------
1. The same Perturb-seq cNMF reference used for the hypoxia analysis (K=100,
   exported via `cnmf consensus --build-ref`).  See the configuration block
   below and project_cNMF_programs_onto_hypoxia_normoxia_STARcat.py for details.

2. The integrated cross-species brain EC AnnData object (cells × genes, raw
   UMI counts in .X or layers['counts']) with:
     obs['species']      — 'mm' (mouse) or 'gp' (guinea pig)
     obs['batch']        — sample/batch identifier
     obs['cell.type.l1'] — broad cell-type annotation

   Guinea pig gene symbols must already have been mapped to mouse ortholog
   symbols during the integration step (performed upstream; not repeated here).

Outputs
-------
- TSV: per-cell STARcat usage matrix (cells × 100 programs), rows = cell
  barcodes, columns = GEP1 … GEP100.
- AnnData (.h5ad): the cross-species brain EC dataset augmented with:
    obsm['X_starcat_usage_norm'] — row-normalized usage (sums to 1 per cell)
    uns['starcat_program_names'] — ordered list of program names (GEP1..GEP100)
  Species, batch, and cell-type metadata are preserved in .obs for downstream
  species-stratified program comparisons.

Dependencies
------------
Python >= 3.10, anndata, scanpy, numpy, pandas, starcatpy (pip install starcatpy)

Reference
---------
Kotliar D, et al. (2019) Identifying gene expression programs of cell-type
identity and cellular activity with single-cell RNA-seq. eLife 8:e43803.
https://doi.org/10.7554/eLife.43803

STARcat package: https://github.com/immunogenomics/starCAT
"""

from pathlib import Path

import numpy as np
import scanpy as sc
from starcat import starCAT

# ---------------------------------------------------------------------------
# Configuration — edit these paths before running
# ---------------------------------------------------------------------------

# Directory and run name for the completed Perturb-seq cNMF run.
# Must match the --output-dir and --name arguments used during cNMF.
# The run must have been finalized with `cnmf consensus --build-ref`.
CNMF_OUTPUT_DIR = Path("perturbseq_endothelial_cNMF")
CNMF_RUN_NAME   = "perturbseq_endothelial_k100"

# Number of programs (K) and local-density threshold used for consensus.
K = 100
DENSITY_THRESHOLD = "0_2"  # encodes 0.2 as used in the spectra filename

# Integrated cross-species brain EC AnnData (raw counts in .X or layers['counts']).
# Cells from mouse (obs['species'] == 'mm') and guinea pig (obs['species'] == 'gp').
# Guinea pig gene symbols must be pre-mapped to mouse orthologs.
QUERY_H5AD = Path("cross_species_brain_EC_mouse_guinea_pig_integrated.h5ad")

# Cache directory for STARcat intermediate files (created automatically).
STARCAT_CACHE_DIR = Path("starcat_projection_cache")

# Output files.
OUTPUT_USAGE_TSV = Path("cross_species_brain_EC_mouse_guinea_pig_STARcat_usage.tsv")
OUTPUT_H5AD      = Path("cross_species_brain_EC_mouse_guinea_pig_with_STARcat_usage.h5ad")

# ---------------------------------------------------------------------------
# Derive the STARcat reference file path from the cNMF run parameters.
# ---------------------------------------------------------------------------
STARCAT_SPECTRA = (
    CNMF_OUTPUT_DIR
    / CNMF_RUN_NAME
    / f"{CNMF_RUN_NAME}.starcat_spectra.k_{K}.dt_{DENSITY_THRESHOLD}.txt"
)

if not STARCAT_SPECTRA.exists():
    raise FileNotFoundError(
        f"STARcat reference not found: {STARCAT_SPECTRA}\n"
        "Re-run cnmf consensus with the --build-ref flag to generate it."
    )

# ---------------------------------------------------------------------------
# Load the cross-species brain EC dataset.
#
# Raw counts are required in .X.  If the object stores raw counts in a layer
# (e.g., layers['counts']), copy them to .X before proceeding.
# ---------------------------------------------------------------------------
print(f"Loading cross-species brain EC dataset from {QUERY_H5AD} …")
adata = sc.read_h5ad(QUERY_H5AD)

# Restore raw counts if the object was log-normalized before saving.
if adata.raw is not None:
    adata.X = adata.raw.X.copy()
elif "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()

print(f"  {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Species composition:\n{adata.obs['species'].value_counts().to_string()}")
print(f"  Cell-type composition:\n{adata.obs['cell.type.l1'].value_counts().to_string()}")

# ---------------------------------------------------------------------------
# Initialize STARcat with the Perturb-seq cNMF reference spectra.
#
# The reference is a 100-program × 10,000-gene matrix of TP10K-scaled gene
# loadings from the Perturb-seq screen (mouse gene space).
#
# Because this query includes guinea pig cells whose gene symbols were mapped
# to mouse orthologs during upstream integration, only the subset of reference
# genes present in the shared ortholog space enters the projection.  STARcat
# prints the number of overlapping genes; expect ~6,700 of 10,000 here
# (compared with ~9,750 for a pure mouse query).
#
# Projection steps performed by STARcat:
#   1. Intersect reference genes with query genes.
#   2. Column-standardize query counts (divide by gene-wise SD, no centering).
#   3. NNLS with spectra held fixed (update_H=False); scikit-learn
#      non_negative_factorization, beta_loss='frobenius', solver='mu',
#      tol=1e-4, max_iter=1000.
#   4. Row-normalize usage matrix so each cell sums to 1.
# ---------------------------------------------------------------------------
print(f"\nInitializing STARcat with reference: {STARCAT_SPECTRA}")
STARCAT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
tcat = starCAT(
    reference=str(STARCAT_SPECTRA),
    cachedir=str(STARCAT_CACHE_DIR),
)
print(f"Reference: {tcat.ref.shape[0]} programs × {tcat.ref.shape[1]} genes")

# ---------------------------------------------------------------------------
# Project cross-species brain ECs onto the Perturb-seq program space.
#
# fit_transform returns:
#   usage_norm — DataFrame (cells × programs), row-normalized to sum to 1.
#   _          — pre-built scores (None for custom references).
# ---------------------------------------------------------------------------
print("\nRunning STARcat projection …")
usage_norm, _ = tcat.fit_transform(adata)

print(f"  Output: {usage_norm.shape[0]} cells × {usage_norm.shape[1]} programs")
print(f"  Row sums (first 5, should be ≈1.0): "
      f"{usage_norm.sum(axis=1).head().values.round(4)}")

# Confirm gene overlap (printed by STARcat; expected ~6,700 for cross-species).
# Lower overlap than a pure mouse query is expected because guinea pig genes
# are restricted to the shared ortholog set.

# ---------------------------------------------------------------------------
# Save the usage matrix as a TSV for archival and downstream analysis.
# Rows = cell barcodes | Columns = GEP1 … GEP100
# ---------------------------------------------------------------------------
usage_norm.to_csv(OUTPUT_USAGE_TSV, sep="\t")
print(f"\nSaved per-cell usage TSV → {OUTPUT_USAGE_TSV}")

# ---------------------------------------------------------------------------
# Embed normalized usage scores and program names in the AnnData object.
#
# obsm['X_starcat_usage_norm'] — normalized usage (sums to 1 per cell);
#     used for species-stratified program comparisons and visualization.
# uns['starcat_program_names'] — program name list in column order.
#
# obs metadata (species, batch, cell.type.l1) is preserved for downstream
# stratified analyses comparing mouse vs guinea pig program activity.
# ---------------------------------------------------------------------------
usage_aligned = usage_norm.loc[adata.obs_names]

adata.obsm["X_starcat_usage_norm"] = usage_aligned.to_numpy()
adata.uns["starcat_program_names"]  = usage_aligned.columns.tolist()

adata.write_h5ad(OUTPUT_H5AD)
print(f"Saved augmented AnnData → {OUTPUT_H5AD}")
print("\nDone.  The output AnnData contains per-cell Perturb-seq program usage")
print("scores for both mouse and guinea pig brain ECs, enabling cross-species")
print("comparison of program activity alongside species and cell-type labels.")
