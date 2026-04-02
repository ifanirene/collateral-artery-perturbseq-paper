#!/usr/bin/env python3
"""Project Perturb-seq cNMF programs onto hypoxia/normoxia brain endothelial cells using STARcat.

Overview
--------
This script projects the 100 gene expression programs (GEPs) discovered by
consensus NMF (cNMF) in the Perturb-seq endothelial dataset onto a separate
scRNA-seq dataset of neonatal mouse brain endothelial cells exposed to either
normoxia (21% O₂) or hypoxia (11% O₂) for 7 days (P1–P8).

The projection uses STARcat v1.0.9 (Kotliar and Curtis, immunogenomics/starCAT),
which holds the cNMF consensus spectra fixed and estimates per-cell program usage
scores by non-negative least squares (NNLS).  Normalized usage scores (summing
to 1 per cell) are stored alongside batch and cell-type annotations for
downstream differential-usage analysis.

Prerequisites
-------------
1. A completed cNMF run on the Perturb-seq endothelial data (K=100, 200
   replicates, 10,000 overdispersed genes, local-density threshold 0.2) with
   the consensus spectra exported via:
       cnmf consensus --build-ref --output-dir <CNMF_OUTPUT_DIR> \
           --name <CNMF_RUN_NAME> --components 100 --local-density-threshold 0.2
   This produces the STARcat reference file at:
       <CNMF_OUTPUT_DIR>/<CNMF_RUN_NAME>/<CNMF_RUN_NAME>.starcat_spectra.k_100.dt_0_2.txt

2. The integrated hypoxia/normoxia AnnData object (cells × genes, raw UMI
   counts in .X) with obs['batch'] labeling each cell as 'Hypoxia' or
   'Normoxia'.

Outputs
-------
- TSV: per-cell STARcat usage matrix (cells × 100 programs), rows = cell
  barcodes, columns = GEP1 … GEP100.
- AnnData (.h5ad): the hypoxia/normoxia dataset augmented with:
    obsm['X_starcat_usage_norm'] — row-normalized usage (sums to 1 per cell)
    uns['starcat_program_names'] — ordered list of program names (GEP1..GEP100)

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

# Directory containing the completed cNMF run on the Perturb-seq dataset.
# The run must have been finalized with `cnmf consensus --build-ref`.
CNMF_OUTPUT_DIR = Path("perturbseq_endothelial_cNMF")

# Name used when calling cnmf prepare/factorize/consensus for this run.
# Must match the --name argument used during the cNMF run.
CNMF_RUN_NAME = "perturbseq_endothelial_k100"

# Number of programs (K) and local-density threshold used for consensus.
K = 100
DENSITY_THRESHOLD = "0_2"  # encodes 0.2 as used in the spectra filename

# Integrated hypoxia/normoxia brain EC AnnData (raw counts in .X).
# Expected: obs['batch'] = 'Hypoxia' or 'Normoxia', obs['cell.type.l1'].
QUERY_H5AD = Path("neonatal_brain_EC_hypoxia_normoxia_integrated.h5ad")

# Cache directory for STARcat intermediate files (created automatically).
STARCAT_CACHE_DIR = Path("starcat_projection_cache")

# Output files.
OUTPUT_USAGE_TSV = Path("neonatal_brain_EC_hypoxia_normoxia_STARcat_usage.tsv")
OUTPUT_H5AD      = Path("neonatal_brain_EC_hypoxia_normoxia_with_STARcat_usage.h5ad")

# ---------------------------------------------------------------------------
# Derive the STARcat reference file path from the cNMF run parameters.
# This file is written by `cnmf consensus --build-ref`.
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
# Load the hypoxia/normoxia brain EC dataset.
# Raw counts are required in .X; STARcat column-standardizes the query
# internally before NNLS fitting.
# ---------------------------------------------------------------------------
print(f"Loading hypoxia/normoxia brain EC dataset from {QUERY_H5AD} …")
adata = sc.read_h5ad(QUERY_H5AD)

# Ensure .X holds raw counts (not log-normalized values).
# If raw counts were stored in adata.raw, restore them here.
if adata.raw is not None:
    adata.X = adata.raw.X.copy()

print(f"  {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Batch composition:\n{adata.obs['batch'].value_counts().to_string()}")

# ---------------------------------------------------------------------------
# Initialize STARcat with the Perturb-seq cNMF reference spectra.
#
# The reference is a 100-program × 10,000-gene matrix of TP10K-scaled gene
# loadings representing the consensus GEPs from the Perturb-seq screen.
# STARcat will:
#   1. Find the intersection of reference genes with the query dataset.
#   2. Column-standardize query counts (divide by gene-wise SD, no centering)
#      to place query values on the same scale as the reference spectra.
#   3. Solve for per-cell usage by NNLS with spectra held fixed
#      (update_H=False), using scikit-learn non_negative_factorization with
#      beta_loss='frobenius', solver='mu', tol=1e-4, max_iter=1000.
#   4. Row-normalize the usage matrix so each cell's scores sum to 1.
# ---------------------------------------------------------------------------
print(f"\nInitializing STARcat with reference: {STARCAT_SPECTRA}")
STARCAT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
tcat = starCAT(
    reference=str(STARCAT_SPECTRA),
    cachedir=str(STARCAT_CACHE_DIR),
)
print(f"Reference: {tcat.ref.shape[0]} programs × {tcat.ref.shape[1]} genes")

# ---------------------------------------------------------------------------
# Project the hypoxia/normoxia cells onto the Perturb-seq program space.
#
# fit_transform returns:
#   usage_norm — DataFrame (cells × programs), row-normalized to sum to 1
#   _          — pre-built scores (None for custom references)
#
# STARcat prints the number of overlapping genes during fitting.
# ---------------------------------------------------------------------------
print("\nRunning STARcat projection …")
usage_norm, _ = tcat.fit_transform(adata)

print(f"  Output: {usage_norm.shape[0]} cells × {usage_norm.shape[1]} programs")
print(f"  Row sums (first 5, should be ≈1.0): "
      f"{usage_norm.sum(axis=1).head().values.round(4)}")

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
#     used for differential analysis and visualization.
# uns['starcat_program_names'] — program name list in column order; needed
#     by downstream scripts that reconstruct the usage DataFrame.
#
# Cell order is aligned explicitly to guard against any index mismatch.
# ---------------------------------------------------------------------------
usage_aligned = usage_norm.loc[adata.obs_names]

adata.obsm["X_starcat_usage_norm"] = usage_aligned.to_numpy()
adata.uns["starcat_program_names"]  = usage_aligned.columns.tolist()

adata.write_h5ad(OUTPUT_H5AD)
print(f"Saved augmented AnnData → {OUTPUT_H5AD}")
print("\nDone.  Next step: run hypoxia_vs_normoxia_program_diff.py to compare")
print("program usage between Hypoxia and Normoxia conditions.")
