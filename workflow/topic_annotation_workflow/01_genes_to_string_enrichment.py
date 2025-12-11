"""
/**
 * @description 
 * Unified CLI to create program-wise top-N gene lists from a cNMF loading CSV and
 * run STRING functional enrichment per program with optional figure generation.
 * 
 * It merges functionality from `annotate_cnmf_programs_string.py` and
 * `run_string_enrichment.py` into subcommands:
 * - extract:  read loading CSV → save JSON {program_id: [genes]} and overview CSV
 * - enrich:   read JSON → call STRING API → write full and filtered CSVs, optionally figures
 * - all:      extract → enrich (convenience)
 * 
 * Key features:
 * - Robust HTTP with retries and pacing
 * - Full unfiltered CSV and Process/KEGG filtered CSV (<500 background genes)
 * - Optional Process/KEGG figure generation per program
 * 
 * @dependencies
 * - pandas, numpy, requests, matplotlib
 * 
 * @examples
 * - Extract top 100 genes per program:
 *   python tools/topic_annotation_workflow/genes_to_string_enrichment.py extract \
 *     --input data/Hypoxia/p2_HnN_loading_gene_k30.csv \
 *     --n-top 100 \
 *     --json-out output/topic_annotations/hnn_k30_top100_genes.json \
 *     --csv-out output/topic_annotations/hnn_k30_top100_genes_overview.csv
 * 
 * - Run enrichment and figures:
 *   python tools/topic_annotation_workflow/genes_to_string_enrichment.py enrich \
 *     --genes-json output/topic_annotations/hnn_k30_top100_genes.json \
 *     --species 10090 \
 *     --out-csv-full output/topic_annotations/hnn_k30_string_full.csv \
 *     --out-csv-filtered output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
 *     --figures-dir output/topic_annotations/enrichment_figures
 * 
 * - End-to-end:
 *   python tools/topic_annotation_workflow/genes_to_string_enrichment.py all \
 *     --input data/Hypoxia/p2_HnN_loading_gene_k30.csv \
 *     --n-top 100 \
 *     --json-out output/topic_annotations/hnn_k30_top100_genes.json \
 *     --csv-out output/topic_annotations/hnn_k30_top100_genes_overview.csv \
 *     --species 10090 \
 *     --out-csv-full output/topic_annotations/hnn_k30_string_full.csv \
 *     --out-csv-filtered output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv
 */
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


# --------------------------- Extract top genes (CSV) --------------------------

def ensure_parent_dir(path_str: str) -> None:
    path = Path(path_str)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)


def extract_top_genes_by_program(df: pd.DataFrame, n_top: int) -> Dict[str, List[str]]:
    required_cols = {"Name", "Score", "RowID"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    top_map: Dict[str, List[str]] = {}
    for program_id, sub in df.groupby("RowID", sort=True):
        sub_sorted = sub.sort_values("Score", ascending=False).head(n_top)
        genes = [str(g) for g in sub_sorted["Name"].dropna().tolist()]
        seen = set()
        unique_genes: List[str] = []
        for g in genes:
            if g not in seen:
                seen.add(g)
                unique_genes.append(g)
        top_map[str(program_id)] = unique_genes
    return top_map


def build_overview_long_table(df: pd.DataFrame, top_map: Dict[str, List[str]]) -> pd.DataFrame:
    records = []
    sub_indexed_cache: Dict[int, pd.DataFrame] = {}
    for program_id_str, genes in top_map.items():
        program_id = int(program_id_str)
        if program_id not in sub_indexed_cache:
            sub_indexed_cache[program_id] = df[df["RowID"] == program_id].set_index("Name")
        sub_idx = sub_indexed_cache[program_id]
        for rank, gene in enumerate(genes, start=1):
            score = float(sub_idx.loc[gene, "Score"]) if gene in sub_idx.index else np.nan
            records.append({"program_id": program_id, "rank": rank, "gene": gene, "score": score})
    out_df = pd.DataFrame.from_records(records)
    if not out_df.empty:
        out_df.sort_values(["program_id", "rank"], inplace=True)
    return out_df


def cmd_extract(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input not found: {input_path}")
        return 2
    df = pd.read_csv(input_path)
    logger.info(f"Loaded input: {input_path} with shape {df.shape} and columns {list(df.columns)}")

    top_map = extract_top_genes_by_program(df=df, n_top=args.n_top)
    logger.info(f"Extracted gene lists for {len(top_map)} programs")

    ensure_parent_dir(args.json_out)
    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(top_map, f, indent=2)
    logger.info(f"Wrote JSON: {args.json_out}")

    overview_df = build_overview_long_table(df, top_map)
    ensure_parent_dir(args.csv_out)
    overview_df.to_csv(args.csv_out, index=False)
    logger.info(f"Wrote overview CSV: {args.csv_out} (rows={overview_df.shape[0]})")
    return 0


# ----------------------------- STRING enrichment -----------------------------

STRING_ENRICH_ENDPOINT = "https://string-db.org/api/json/enrichment"


def call_string_enrichment(genes: List[str], species: int, retries: int = 3, sleep_between: float = 0.6) -> List[Dict[str, Any]]:
    identifiers_value = "\r".join(genes)
    params = {"identifiers": identifiers_value, "species": species, "caller_identity": "topic_analysis_string_enrichment"}

    attempt = 0
    while attempt <= retries:
        try:
            response = requests.get(STRING_ENRICH_ENDPOINT, params=params, timeout=60)
            if response.status_code == 200:
                try:
                    data = response.json()
                except Exception as json_err:
                    logger.error(f"Failed to parse JSON (n={len(genes)}): {json_err}")
                    data = []
                return data if isinstance(data, list) else []
            else:
                logger.warning(f"STRING returned status {response.status_code}: {response.text[:200]}")
        except requests.RequestException as e:
            logger.warning(f"HTTP error on STRING request (attempt {attempt+1}/{retries+1}): {e}")

        attempt += 1
        time.sleep(min(2.0, sleep_between * (attempt + 1)))

    return []


def build_full_csv(program_to_results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for pid, terms in program_to_results.items():
        for t in terms:
            rows.append(
                {
                    "program_id": int(pid),
                    "category": str(t.get("category", "NA")),
                    "term": str(t.get("term", t.get("description", "NA"))),
                    "term_id": str(t.get("term_id", "NA")),
                    "description": str(t.get("description", t.get("term", "NA"))),
                    "fdr": float(t.get("fdr", float("nan"))),
                    "p_value": float(t.get("p_value", float("nan"))),
                    "number_of_genes": int(t.get("number_of_genes", 0)),
                    "number_of_genes_in_background": int(t.get("number_of_genes_in_background", 0)),
                    "ncbiTaxonId": int(t.get("ncbiTaxonId", 0)),
                    "inputGenes": "|".join(t.get("inputGenes", [])) if t.get("inputGenes") else "",
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(["program_id", "fdr", "p_value"], inplace=True)
    return df


def filter_process_kegg(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    category_mask = df["category"].str.contains("Process|KEGG", case=False, na=False)
    background_mask = df["number_of_genes_in_background"] < 500
    filtered_df = df[category_mask & background_mask].copy()
    if not filtered_df.empty:
        filtered_df.sort_values(["program_id", "fdr", "p_value"], inplace=True)
    return filtered_df


def cmd_enrich(args: argparse.Namespace) -> int:
    genes_path = Path(args.genes_json)
    if not genes_path.exists():
        logger.error(f"Genes JSON not found: {genes_path}")
        return 2

    program_to_genes: Dict[str, List[str]] = json.loads(genes_path.read_text(encoding="utf-8"))
    logger.info(f"Loaded gene lists for {len(program_to_genes)} programs from {genes_path}")

    program_to_results: Dict[str, List[Dict[str, Any]]] = {}
    total = len(program_to_genes)
    for idx, program_id in enumerate(sorted(program_to_genes.keys(), key=lambda x: int(x)), start=1):
        genes = [g for g in program_to_genes[program_id] if isinstance(g, str) and g.strip()]
        logger.info(f"[{idx}/{total}] STRING enrichment for program {program_id} with {len(genes)} genes ...")
        results = call_string_enrichment(genes=genes, species=args.species, retries=args.retries, sleep_between=args.sleep)
        logger.info(f"Program {program_id}: retrieved {len(results)} enriched terms")
        program_to_results[program_id] = results
        time.sleep(args.sleep)

    if args.out_csv_full:
        out_csv_full_path = Path(args.out_csv_full)
        df_full = build_full_csv(program_to_results)
        out_csv_full_path.parent.mkdir(parents=True, exist_ok=True)
        df_full.to_csv(out_csv_full_path, index=False)
        logger.info(f"Wrote full unfiltered CSV with {len(df_full)} rows to {out_csv_full_path}")

    if args.out_csv_filtered:
        out_csv_filtered_path = Path(args.out_csv_filtered)
        df_full = build_full_csv(program_to_results)
        df_filtered = filter_process_kegg(df_full)
        out_csv_filtered_path.parent.mkdir(parents=True, exist_ok=True)
        df_filtered.to_csv(out_csv_filtered_path, index=False)
        logger.info(f"Wrote filtered CSV (Process/KEGG) with {len(df_filtered)} rows to {out_csv_filtered_path}")

    # Optional simple figures (Process and KEGG top 10)
    if args.figures_dir:
        figures_dir = Path(args.figures_dir)
        figures_dir.mkdir(parents=True, exist_ok=True)

        def create_enrichment_plot(enrichment_data: List[Dict[str, Any]], category: str, output_path: Path, top_n: int = 10) -> bool:
            try:
                filtered_data = [term for term in enrichment_data if term.get("category", "").lower() == category.lower()]
                if not filtered_data:
                    return False
                sorted_data = sorted(filtered_data, key=lambda x: float(x.get("fdr", float("inf"))))[:top_n]
                if not sorted_data:
                    return False
                terms = []
                fdrs = []
                gene_counts = []
                for term in sorted_data:
                    term_name = str(term.get("term", term.get("description", "Unknown")))
                    if len(term_name) > 50:
                        term_name = term_name[:47] + "..."
                    terms.append(term_name)
                    fdrs.append(float(term.get("fdr", 1.0)))
                    gene_counts.append(int(term.get("number_of_genes", 0)))
                fig, ax = plt.subplots(figsize=(10, max(6, len(terms) * 0.4)))
                x_values = [-np.log10(max(fdr, 1e-10)) for fdr in fdrs]
                y_positions = range(len(terms))
                sizes = [max(20, count * 5) for count in gene_counts]
                scatter = ax.scatter(x_values, y_positions, s=sizes, alpha=0.7, c=x_values, cmap="viridis")
                ax.set_yticks(list(y_positions))
                ax.set_yticklabels(terms)
                ax.set_xlabel("-log10(FDR)")
                ax.set_title(f"{category} Enrichment (Top {len(terms)} Terms)")
                ax.grid(True, alpha=0.3)
                cbar = plt.colorbar(scatter)
                cbar.set_label("-log10(FDR)")
                plt.tight_layout()
                plt.savefig(output_path, dpi=150, bbox_inches="tight")
                plt.close()
                return True
            except Exception:
                return False

        for program_id, terms in program_to_results.items():
            ok_p = create_enrichment_plot(terms, "Process", figures_dir / f"program_{program_id}_process_enrichment.png")
            ok_k = create_enrichment_plot(terms, "KEGG", figures_dir / f"program_{program_id}_kegg_enrichment.png")
            logger.info(f"Program {program_id}: figures - Process={'✓' if ok_p else '✗'}, KEGG={'✓' if ok_k else '✗'}")

    return 0


# -------------------------------- Entry points -------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract program gene lists and run STRING enrichment.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # extract
    p_extract = subparsers.add_parser("extract", help="Extract top-N genes per program from loading CSV")
    p_extract.add_argument("--input", required=True, help="Path to cNMF loading CSV (columns: Name, Score, RowID)")
    p_extract.add_argument("--n-top", type=int, default=100, help="Number of top genes per program to extract")
    p_extract.add_argument("--json-out", required=True, help="Output JSON {program_id: [genes...]}")
    p_extract.add_argument("--csv-out", required=True, help="Output overview CSV")
    p_extract.set_defaults(func=cmd_extract)

    # enrich
    p_enrich = subparsers.add_parser("enrich", help="Run STRING enrichment for program gene lists from JSON")
    p_enrich.add_argument("--genes-json", required=True, help="Path to JSON mapping {program_id: [genes...]}")
    p_enrich.add_argument("--species", type=int, default=10090, help="NCBI/STRING species id (default: 10090 mouse)")
    p_enrich.add_argument("--out-csv-full", help="Full unfiltered CSV output path")
    p_enrich.add_argument("--out-csv-filtered", help="Filtered CSV (Process/KEGG only, background<500)")
    p_enrich.add_argument("--figures-dir", help="Directory to save enrichment figures")
    p_enrich.add_argument("--sleep", type=float, default=0.6, help="Sleep seconds between API calls")
    p_enrich.add_argument("--retries", type=int, default=3, help="Retries per program on HTTP failures")
    p_enrich.set_defaults(func=cmd_enrich)

    # all
    p_all = subparsers.add_parser("all", help="Run extract then enrich")
    # extract args
    p_all.add_argument("--input", required=True)
    p_all.add_argument("--n-top", type=int, default=100)
    p_all.add_argument("--json-out", required=True)
    p_all.add_argument("--csv-out", required=True)
    # enrich args
    p_all.add_argument("--species", type=int, default=10090)
    p_all.add_argument("--out-csv-full")
    p_all.add_argument("--out-csv-filtered")
    p_all.add_argument("--figures-dir")
    p_all.add_argument("--sleep", type=float, default=0.6)
    p_all.add_argument("--retries", type=int, default=3)

    def run_all(args: argparse.Namespace) -> int:
        rc = cmd_extract(args)
        if rc != 0:
            return rc
        enrich_args = argparse.Namespace(
            genes_json=args.json_out,
            species=args.species,
            out_csv_full=args.out_csv_full,
            out_csv_filtered=args.out_csv_filtered,
            figures_dir=args.figures_dir,
            sleep=args.sleep,
            retries=args.retries,
            func=cmd_enrich,
        )
        return cmd_enrich(enrich_args)

    p_all.set_defaults(func=run_all)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    rc = args.func(args)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()


