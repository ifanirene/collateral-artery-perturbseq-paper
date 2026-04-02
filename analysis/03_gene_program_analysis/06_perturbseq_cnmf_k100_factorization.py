#!/usr/bin/env python
"""
@description
This script runs the retained Final Pool endothelial cNMF program-discovery workflow
without the earlier exploratory notebook steps.
It is responsible for restoring raw counts from the source AnnData object,
preparing the cNMF inputs, running the retained final no-batch-correction
`k = 100` factorization workflow,
and building the consensus outputs used downstream in the paper.

Key features:
- Restores raw UMI counts from an AnnData `raw` layer into a cNMF-ready `.h5ad`
- Runs only the retained production settings (`k = 100`, `n_iter = 200`,
  `num_highvar_genes = 10000`, `max_NMF_iter = 2000`, `dt = 0.2`)
- Supports per-worker factorization for parallel execution plus a sequential fallback
- Prints descriptive next-step commands and expected output filenames for debugging
- Makes the retained run identity explicit in the file name, run label, logs, and defaults

@dependencies
- scanpy: Reads and writes AnnData files and restores raw counts
- cnmf: Performs consensus NMF preparation, factorization, combination, and consensus
- argparse, pathlib: Provides a lightweight command-line interface

@examples
# Prepare the raw-count input file and cNMF job ledger
python src/final_pool_endothelial_cnmf_retained_k100_run_submission_clean_2026-04-01.py prepare \
  --source-h5ad data/final_pool_endothelial_cells.h5ad

# Run worker 3 of 12 in parallel
python src/final_pool_endothelial_cnmf_retained_k100_run_submission_clean_2026-04-01.py factorize-worker \
  --worker-index 3 --total-workers 12 --skip-completed-runs

# Finish the retained production run
python src/final_pool_endothelial_cnmf_retained_k100_run_submission_clean_2026-04-01.py combine
python src/final_pool_endothelial_cnmf_retained_k100_run_submission_clean_2026-04-01.py consensus
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

if TYPE_CHECKING:
    from cnmf import cNMF

DEFAULT_SOURCE_H5AD = Path("data/final_pool_endothelial_cells.h5ad")
DEFAULT_COUNTS_H5AD = Path("data/final_pool_endothelial_cells_raw_counts_for_cnmf.h5ad")
DEFAULT_OUTPUT_DIR = Path("output/cnmf_final_pool_endothelial_retained_k100_release")
DEFAULT_RUN_NAME = "final_pool_endothelial_no_batch_correction_k100"
DEFAULT_COMPONENTS = 100
DEFAULT_N_ITER = 200
DEFAULT_SEED = 123
DEFAULT_NUM_HIGHVAR_GENES = 10000
DEFAULT_MAX_NMF_ITER = 2000
DEFAULT_DENSITY_THRESHOLD = 0.2
DEFAULT_TOTAL_WORKERS = 12


def log(message: str) -> None:
    """Print a prefixed status message."""
    print(f"[final_pool_endothelial_cnmf] {message}", flush=True)


def density_threshold_tag(density_threshold: float) -> str:
    """Match cNMF's output naming convention for density thresholds."""
    return str(density_threshold).replace(".", "_")


def expected_outputs(
    output_dir: Path,
    run_name: str,
    components: int,
    density_threshold: float,
) -> list[Path]:
    """Return the main consensus outputs expected from the retained run."""
    density_tag = density_threshold_tag(density_threshold)
    run_dir = output_dir / run_name
    prefix = f"{run_name}"
    return [
        run_dir / f"{prefix}.clustering.k_{components}.dt_{density_tag}.png",
        run_dir / f"{prefix}.gene_spectra_score.k_{components}.dt_{density_tag}.txt",
        run_dir / f"{prefix}.gene_spectra_tpm.k_{components}.dt_{density_tag}.txt",
        run_dir / f"{prefix}.spectra.k_{components}.dt_{density_tag}.consensus.txt",
        run_dir / f"{prefix}.starcat_spectra.k_{components}.dt_{density_tag}.txt",
        run_dir / f"{prefix}.usages.k_{components}.dt_{density_tag}.consensus.txt",
    ]


def get_scanpy_module() -> Any:
    """Import scanpy only when a stage actually needs it."""
    import scanpy as sc

    return sc


def get_cnmf_class() -> Any:
    """Import cNMF only when a stage actually needs it."""
    from cnmf import cNMF

    return cNMF


def build_cnmf(output_dir: Path, run_name: str) -> Any:
    """Initialize the cNMF object for the retained Final Pool endothelial run."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cNMF = get_cnmf_class()
    return cNMF(output_dir=str(output_dir), name=run_name)


def restore_raw_counts_to_h5ad(
    source_h5ad: Path,
    counts_h5ad: Path,
    force: bool,
) -> Path:
    """Write a cNMF-ready `.h5ad` with raw UMI counts in `.X`."""
    if counts_h5ad.exists() and not force:
        log(f"Using existing raw-count AnnData: {counts_h5ad}")
        return counts_h5ad

    if not source_h5ad.exists():
        raise FileNotFoundError(
            f"Source AnnData not found: {source_h5ad}. "
            "Pass --source-h5ad to the endothelial AnnData file with a populated raw layer."
        )

    log(f"Reading source AnnData: {source_h5ad}")
    sc = get_scanpy_module()
    adata = sc.read_h5ad(source_h5ad)
    if adata.raw is None:
        raise ValueError(
            "The source AnnData object does not contain a raw layer. "
            "Provide a pre-restored counts file with --counts-h5ad or use a source file with adata.raw."
        )

    raw_adata = adata.raw.to_adata()
    raw_adata.obs = adata.obs.copy()
    if "_index" in raw_adata.var.columns and "features" not in raw_adata.var.columns:
        raw_adata.var = raw_adata.var.rename(columns={"_index": "features"})

    counts_h5ad.parent.mkdir(parents=True, exist_ok=True)
    log(f"Writing raw-count AnnData for cNMF: {counts_h5ad}")
    raw_adata.write_h5ad(counts_h5ad)
    return counts_h5ad


def print_expected_outputs(
    output_dir: Path,
    run_name: str,
    components: int,
    density_threshold: float,
) -> None:
    """Print the main output files expected from the consensus step."""
    log("Expected final consensus outputs:")
    for path in expected_outputs(output_dir, run_name, components, density_threshold):
        print(f"  - {path}")


def print_parallel_commands(
    script_path: Path,
    total_workers: int,
    common_args: Iterable[str],
) -> None:
    """Print example commands for running factorization workers in parallel."""
    display_script_path = Path("src") / script_path.name
    quoted_args = " ".join(common_args)
    log("Suggested retained-run commands for the no-batch-correction k=100 workflow:")
    print(f"  python {display_script_path.as_posix()} prepare {quoted_args}")
    print(
        "  parallel "
        f"python {display_script_path.as_posix()} factorize-worker {quoted_args} "
        f"--total-workers {total_workers} --worker-index {{}} --skip-completed-runs ::: "
        + " ".join(str(worker_i) for worker_i in range(total_workers))
    )
    print(f"  python {display_script_path.as_posix()} combine {quoted_args}")
    print(f"  python {display_script_path.as_posix()} consensus {quoted_args}")


def run_prepare(args: argparse.Namespace) -> None:
    """Restore raw counts when needed and build the cNMF job ledger."""
    counts_h5ad = restore_raw_counts_to_h5ad(
        source_h5ad=args.source_h5ad,
        counts_h5ad=args.counts_h5ad,
        force=args.force_rewrite_counts,
    )

    cnmf_obj = build_cnmf(args.output_dir, args.run_name)
    log(
        "Preparing retained final cNMF run "
        f"(run_name={args.run_name}, k={args.components}, n_iter={args.n_iter}, "
        f"seed={args.seed}, num_highvar_genes={args.num_highvar_genes}, "
        f"max_nmf_iter={args.max_nmf_iter}, density_threshold={args.density_threshold})"
    )
    cnmf_obj.prepare(
        counts_fn=str(counts_h5ad),
        components=[args.components],
        n_iter=args.n_iter,
        seed=args.seed,
        num_highvar_genes=args.num_highvar_genes,
        beta_loss="frobenius",
        init="random",
        max_NMF_iter=args.max_nmf_iter,
    )

    print_parallel_commands(
        script_path=Path(__file__),
        total_workers=args.total_workers,
        common_args=common_cli_args(args),
    )


def run_factorize_worker(args: argparse.Namespace) -> None:
    """Run the assigned cNMF worker shard."""
    cnmf_obj = build_cnmf(args.output_dir, args.run_name)
    log(
        f"Starting factorization worker {args.worker_index} of {args.total_workers} "
        f"(skip_completed_runs={args.skip_completed_runs})"
    )
    cnmf_obj.factorize(
        worker_i=args.worker_index,
        total_workers=args.total_workers,
        skip_completed_runs=args.skip_completed_runs,
    )


def run_combine(args: argparse.Namespace) -> None:
    """Merge per-worker spectra for the retained component count."""
    cnmf_obj = build_cnmf(args.output_dir, args.run_name)
    log(f"Combining retained factorization runs for k={args.components}")
    cnmf_obj.combine(components=[args.components])


def run_consensus(args: argparse.Namespace) -> None:
    """Build the retained consensus solution and write the final program outputs."""
    cnmf_obj = build_cnmf(args.output_dir, args.run_name)
    log(
        f"Running retained consensus step for k={args.components} "
        f"with density_threshold={args.density_threshold}"
    )
    cnmf_obj.consensus(
        k=args.components,
        density_threshold=args.density_threshold,
        show_clustering=True,
        close_clustergram_fig=True,
    )
    print_expected_outputs(
        output_dir=args.output_dir,
        run_name=args.run_name,
        components=args.components,
        density_threshold=args.density_threshold,
    )


def run_all_sequential(args: argparse.Namespace) -> None:
    """Run the full retained workflow serially as a simple one-command fallback."""
    run_prepare(args)
    for worker_index in range(args.total_workers):
        worker_args = argparse.Namespace(**vars(args))
        worker_args.worker_index = worker_index
        run_factorize_worker(worker_args)
    run_combine(args)
    run_consensus(args)


def common_cli_args(args: argparse.Namespace) -> list[str]:
    """Build the shared argument list used in printed example commands."""
    shared_args = [
        f"--counts-h5ad {args.counts_h5ad.as_posix()}",
        f"--output-dir {args.output_dir.as_posix()}",
        f"--run-name {args.run_name}",
        f"--components {args.components}",
        f"--n-iter {args.n_iter}",
        f"--seed {args.seed}",
        f"--num-highvar-genes {args.num_highvar_genes}",
        f"--max-nmf-iter {args.max_nmf_iter}",
        f"--density-threshold {args.density_threshold}",
    ]
    if hasattr(args, "source_h5ad"):
        shared_args.append(f"--source-h5ad {args.source_h5ad.as_posix()}")
    if getattr(args, "force_rewrite_counts", False):
        shared_args.append("--force-rewrite-counts")
    return shared_args


def add_shared_runtime_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared across retained-run stages."""
    parser.add_argument("--counts-h5ad", type=Path, default=DEFAULT_COUNTS_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-name", default=DEFAULT_RUN_NAME)
    parser.add_argument("--components", type=int, default=DEFAULT_COMPONENTS)
    parser.add_argument("--n-iter", type=int, default=DEFAULT_N_ITER)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--num-highvar-genes", type=int, default=DEFAULT_NUM_HIGHVAR_GENES)
    parser.add_argument("--max-nmf-iter", type=int, default=DEFAULT_MAX_NMF_ITER)
    parser.add_argument("--density-threshold", type=float, default=DEFAULT_DENSITY_THRESHOLD)
    parser.add_argument("--total-workers", type=int, default=DEFAULT_TOTAL_WORKERS)


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Standalone retained Final Pool endothelial cNMF script for the "
            "final no-batch-correction 100-program production run."
        )
    )
    subparsers = parser.add_subparsers(dest="stage", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Restore raw counts if needed and prepare the retained cNMF run.",
    )
    add_shared_runtime_arguments(prepare_parser)
    prepare_parser.add_argument("--source-h5ad", type=Path, default=DEFAULT_SOURCE_H5AD)
    prepare_parser.add_argument(
        "--force-rewrite-counts",
        action="store_true",
        help="Rewrite the derived raw-count `.h5ad` even if it already exists.",
    )

    factorize_parser = subparsers.add_parser(
        "factorize-worker",
        help="Run one retained cNMF factorization worker shard.",
    )
    add_shared_runtime_arguments(factorize_parser)
    factorize_parser.add_argument("--worker-index", type=int, required=True)
    factorize_parser.add_argument(
        "--skip-completed-runs",
        action="store_true",
        help="Skip worker jobs whose spectra files already exist.",
    )

    combine_parser = subparsers.add_parser(
        "combine",
        help="Combine retained factorization outputs for the final k=100 run.",
    )
    add_shared_runtime_arguments(combine_parser)

    consensus_parser = subparsers.add_parser(
        "consensus",
        help="Build the retained consensus solution and write final outputs.",
    )
    add_shared_runtime_arguments(consensus_parser)

    all_parser = subparsers.add_parser(
        "all-sequential",
        help="Run the full retained workflow serially as a simple fallback.",
    )
    add_shared_runtime_arguments(all_parser)
    all_parser.add_argument("--source-h5ad", type=Path, default=DEFAULT_SOURCE_H5AD)
    all_parser.add_argument(
        "--force-rewrite-counts",
        action="store_true",
        help="Rewrite the derived raw-count `.h5ad` even if it already exists.",
    )
    all_parser.add_argument(
        "--skip-completed-runs",
        action="store_true",
        help="Skip worker jobs whose spectra files already exist during serial factorization.",
    )

    commands_parser = subparsers.add_parser(
        "print-parallel-commands",
        help="Print example commands for the retained production run.",
    )
    add_shared_runtime_arguments(commands_parser)

    return parser


def main() -> None:
    """Parse arguments and dispatch to the requested retained-run stage."""
    parser = build_parser()
    args = parser.parse_args()

    if args.stage == "prepare":
        run_prepare(args)
    elif args.stage == "factorize-worker":
        run_factorize_worker(args)
    elif args.stage == "combine":
        run_combine(args)
    elif args.stage == "consensus":
        run_consensus(args)
    elif args.stage == "all-sequential":
        run_all_sequential(args)
    elif args.stage == "print-parallel-commands":
        print_parallel_commands(
            script_path=Path(__file__),
            total_workers=args.total_workers,
            common_args=common_cli_args(args),
        )
    else:
        parser.error(f"Unsupported stage: {args.stage}")


if __name__ == "__main__":
    main()
