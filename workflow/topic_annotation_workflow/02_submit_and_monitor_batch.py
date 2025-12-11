"""
/**
 * @description 
 * Unified CLI to prepare Anthropic batch requests for topic annotations, submit them,
 * monitor job status, and download the results when complete.
 * 
 * It merges functionality from `01_submit_batch_annotation.py` and
 * `check_and_download_batch.py` into one consolidated interface with subcommands:
 * - prepare: build the batch request JSON only
 * - submit: build and submit the batch; prints Batch ID
 * - monitor: poll a Batch ID; optionally download results
 * - submit-monitor: build → submit → monitor → download (end-to-end)
 * 
 * Key features:
 * - Enrichment context (STRING KEGG/Process) can be optionally included in prompts
 * - Robust polling with informative logging and optional result download
 * - Consistent error handling and clear CLI ergonomics
 * 
 * @dependencies
 * - pandas: CSV loading and grouping
 * - requests: HTTP calls to Anthropic API
 * - python-dotenv: load API key from .env
 * 
 * @examples
 * - Prepare only (test first 2 topics):
 *   python tools/topic_annotation_workflow/submit_and_monitor_batch.py \
 *     prepare \
 *     --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
 *     --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
 *     --num-topics 2 \
 *     --output-file output/topic_annotations/hypoxia_k30/anthropic_batch_request_with_enrichment_test.json
 * 
 * - Submit now (all topics):
 *   python tools/topic_annotation_workflow/submit_and_monitor_batch.py \
 *     submit \
 *     --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
 *     --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv
 * 
 * - Monitor and download:
 *   python tools/topic_annotation_workflow/submit_and_monitor_batch.py \
 *     monitor \
 *     --batch-id <BATCH_ID> \
 *     --output-file output/topic_annotations/hypoxia_k30/batch_results.jsonl
 * 
 * - End-to-end (submit then wait and download):
 *   python tools/topic_annotation_workflow/submit_and_monitor_batch.py \
 *     submit-monitor \
 *     --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
 *     --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
 *     --download-to output/topic_annotations/hypoxia_k30/batch_results.jsonl \
 *     --poll-interval 30
 */
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from typing import Dict, List, Optional

import pandas as pd
import requests
from dotenv import load_dotenv


# Load environment variables from .env file (if present)
load_dotenv()

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


# ----------------------------- Anthropic helpers -----------------------------

ANTHROPIC_API_URL = "https://api.anthropic.com/v1/messages/batches"
ANTHROPIC_API_VERSION = "2023-06-01"
MODEL = "claude-3-7-sonnet-20250219"


def check_batch_status(batch_id: str, api_key: str) -> Optional[dict]:
    """Checks the status of a given batch job."""
    url = f"{ANTHROPIC_API_URL}/{batch_id}"
    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_API_VERSION,
        "content-type": "application/json",
    }
    try:
        response = requests.get(url, headers=headers, timeout=60)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Error checking batch status: {e}")
        if getattr(e, "response", None) is not None:
            logger.error(f"Response Body: {e.response.text}")
        return None


def download_file(url: str, destination: str, api_key: str) -> bool:
    """Downloads a file from a URL to a destination."""
    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_API_VERSION,
    }
    try:
        with requests.get(url, headers=headers, timeout=300, stream=True) as response:
            response.raise_for_status()
            with open(destination, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
        logger.info(f"Successfully downloaded file to {destination}")
        return True
    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading file: {e}")
        return False


# ------------------------- Enrichment context helpers ------------------------

def prepare_enrichment_mapping(enrichment_file: Optional[str]) -> Dict[int, Dict[str, List[dict]]]:
    """Load enrichment CSV and build a mapping: program_id -> {category -> [rows sorted by FDR]}.

    Expected columns: program_id, category, description, fdr, inputGenes
    """
    if not enrichment_file:
        return {}
    if not os.path.exists(enrichment_file):
        logger.error(f"Enrichment file not found: {enrichment_file}")
        return {}

    df = pd.read_csv(enrichment_file)
    required_cols = {"program_id", "category", "description", "fdr", "inputGenes"}
    missing = required_cols - set(df.columns)
    if missing:
        logger.error(f"Enrichment file missing required columns: {missing}")
        return {}

    df["category"] = df["category"].astype(str).str.strip()
    df_sorted = df.sort_values(["program_id", "category", "fdr"], ascending=[True, True, True])

    enrichment_by_program: Dict[int, Dict[str, List[dict]]] = {}
    for (pid, cat), sub in df_sorted.groupby(["program_id", "category"], sort=False):
        if cat not in ("KEGG", "Process"):
            continue
        program_map = enrichment_by_program.setdefault(int(pid), {})
        program_map[cat] = sub.to_dict(orient="records")

    logger.info(f"Prepared enrichment mapping for {len(enrichment_by_program)} programs.")
    return enrichment_by_program


def build_enrichment_context(
    enrichment_by_program: Dict[int, Dict[str, List[dict]]], topic_id: int, top_enrichment: int
) -> str:
    """Builds a concise textual enrichment context for a topic."""
    program_context = enrichment_by_program.get(topic_id, {})
    if not program_context:
        return ""

    lines: List[str] = []
    for category in ("KEGG", "Process"):
        rows = program_context.get(category, [])
        if not rows:
            continue
        lines.append(f"- {category} (top {min(top_enrichment, len(rows))}):")
        for row in rows[:top_enrichment]:
            desc = row.get("description") or row.get("term") or "NA"
            fdr = row.get("fdr")
            fdr_str = f"{float(fdr):.2e}" if isinstance(fdr, (float, int)) else str(fdr)
            genes = row.get("inputGenes", "")
            gene_list = genes.split("|") if isinstance(genes, str) else []
            if len(gene_list) > 25:
                shown = "|".join(gene_list[:25])
                genes_str = f"{shown}|… (+{len(gene_list) - 25} more)"
            else:
                genes_str = genes if isinstance(genes, str) else ""
            lines.append(f"  - {desc} (FDR={fdr_str}) — genes: {genes_str}")

    return "\n".join(["Optional enrichment context (STRING):", *lines]) if lines else ""


def generate_prompt(
    topic_id: int,
    gene_df: pd.DataFrame,
    prompt_template: str,
    max_genes: Optional[int],
    enrichment_by_program: Dict[int, Dict[str, List[dict]]],
    top_enrichment: int,
) -> Optional[str]:
    """Generates a prompt for a given topic ID."""
    topic_genes = (
        gene_df[gene_df["RowID"] == topic_id]
        .sort_values("Score", ascending=False)["Name"]
        .tolist()
    )
    if not topic_genes:
        logger.warning(f"No genes found for Topic {topic_id}")
        return None

    if max_genes is not None and max_genes > 0:
        topic_genes = topic_genes[:max_genes]

    gene_list_str = ", ".join(topic_genes)
    enrichment_context = build_enrichment_context(
        enrichment_by_program=enrichment_by_program,
        topic_id=topic_id,
        top_enrichment=top_enrichment,
    )

    return (
        prompt_template.replace("{gene_list}", gene_list_str)
        .replace("{topic_id}", str(topic_id))
        .replace("{enrichment_context}", enrichment_context)
    )


PROMPT_TEMPLATE = """
You are a vascular-biology specialist. Analyse **Topic {topic_id}**, a gene program extracted from single-cell Perturb-seq of brain endothelial cells (ECs).
Your goal is to help understand this program and map it onto EC fate trajectories, states and cellular phenotypes.

Gene list (top-loading): {gene_list}

{enrichment_context}
 
 Guidance: Treat the gene list as primary evidence. Optionally consult the enrichment context to cross-check themes, but do not rely on it as the main source. Avoid repeating enrichment terms verbatim; focus on synthesis grounded in the genes, and note enrichment only when it strengthens or challenges an interpretation.

────────────────────────────────────────────────────────
### Output requirements
Return the entire answer in GitHub-flavoured Markdown. Use standard HUGO gene symbols. Keep language precise and concise (avoid over-speculation).

CRITICALLY, include the following two lines near the top, exactly with these bold labels:

**Brief Summary:** <one to two sentences, concise>
**Three Key Words:** <word1>, <word2>, <word3>

Then provide the following sections:

1. **High-level overview (≤ 120 words)**
   - Main theme(s) and biological tendency of the set.

2. **Key functional clusters**
   Group the genes (give the group a short name) and briefly summarise their shared roles.
   - Cluster name: [1-sentence summary] — genes: ...
   - Highlight any hallmark genes or distinctive markers that stand out.
   - Summarize the major states or processes this program might represent.

3. **Key Pathways**
   List major pathways represented and example genes, including:
   - Pathway names (signaling, metabolic etc.)
   - Transcriptional regulations
   - Other plausible gene interactions
   - Positive vs negative regulators of the highlighted processes.

4. **Clinical Relevance**
   List potential clinical implications or disease associations, including:
   - Neurological / vascular diseases, BBB integrity, stroke, tumour angiogenesis, etc.
   - Potential therapeutic targets (rank 1-3).

5. **Statistical Summary**
   Provide quantitative analysis:
   - Total number of genes
   - Distribution across functional categories (percentages)
   - Notable patterns or biases

Start the document with the H2 heading: `## Topic {topic_id} annotation`.
"""


# --------------------------------- Commands ---------------------------------

def cmd_prepare(args: argparse.Namespace) -> int:
    """Prepare a batch request JSON file for the given gene CSV."""
    try:
        gene_df = pd.read_csv(args.gene_file)
    except FileNotFoundError as e:
        logger.error(f"Error reading gene file: {e}")
        return 2

    topic_ids = sorted(gene_df["RowID"].unique())
    if args.num_topics:
        topic_ids = topic_ids[: args.num_topics]
        logger.info(f"Limiting to first {args.num_topics} topics for testing.")

    enrichment_by_program = prepare_enrichment_mapping(args.enrichment_file)

    batch_requests: List[dict] = []
    for topic_id in topic_ids:
        prompt = generate_prompt(
            topic_id=topic_id,
            gene_df=gene_df,
            prompt_template=PROMPT_TEMPLATE,
            max_genes=args.max_genes,
            enrichment_by_program=enrichment_by_program,
            top_enrichment=args.top_enrichment,
        )
        if prompt:
            request = {
                "custom_id": f"topic_{topic_id}_annotation",
                "params": {
                    "model": MODEL,
                    "max_tokens": 3500,
                    "messages": [{"role": "user", "content": prompt}],
                },
            }
            batch_requests.append(request)

    if not batch_requests:
        logger.error("No requests were generated. Aborting.")
        return 3

    payload = {"requests": batch_requests}

    try:
        with open(args.output_file, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        logger.info(f"Successfully created batch request file at {args.output_file}")
        return 0
    except OSError as e:
        logger.error(f"Error writing to file: {e}")
        return 4


def cmd_submit(args: argparse.Namespace) -> int:
    """Build and submit the batch; print Batch ID."""
    rc = cmd_prepare(args)
    if rc != 0:
        return rc

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        logger.error("ANTHROPIC_API_KEY environment variable not set.")
        return 5

    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_API_VERSION,
        "content-type": "application/json",
    }
    try:
        with open(args.output_file, "r", encoding="utf-8") as f:
            payload = json.load(f)
        logger.info(f"Submitting batch request for {len(payload.get('requests', []))} topics...")
        response = requests.post(ANTHROPIC_API_URL, headers=headers, data=json.dumps(payload), timeout=120)
        response.raise_for_status()
        response_data = response.json()
        batch_id = response_data.get("id")
        status = response_data.get("processing_status")
        logger.info("Successfully submitted batch job!")
        logger.info(f"Batch ID: {batch_id}")
        logger.info(f"Status: {status}")
        print(batch_id)
        return 0
    except requests.exceptions.RequestException as e:
        logger.error(f"An error occurred while calling the Anthropic API: {e}")
        if getattr(e, "response", None) is not None:
            logger.error(f"Response Body: {e.response.text}")
        return 6
    except Exception as e:
        logger.error(f"Unexpected error during submission: {e}")
        return 7


def _monitor_once(batch_id: str, api_key: str) -> Optional[dict]:
    status_data = check_batch_status(batch_id, api_key)
    if not status_data:
        return None
    status = status_data.get("processing_status")
    logger.info(f"Current job status: {status}")
    return status_data


def cmd_monitor(args: argparse.Namespace) -> int:
    """Poll a Batch ID and optionally download the results when complete."""
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        logger.error("ANTHROPIC_API_KEY environment variable not set.")
        return 5

    logger.info(f"Starting to poll for batch job: {args.batch_id}")

    while True:
        status_data = _monitor_once(args.batch_id, api_key)
        if not status_data:
            logger.error("Failed to get batch status. Aborting.")
            return 8

        status = status_data.get("processing_status")
        if status in ["completed", "ended"]:
            logger.info(f"Batch job {status}!")

            # Check for both possible output file URL fields
            output_file_url = status_data.get("output_file_url") or status_data.get("results_url")

            if output_file_url and args.output_file:
                logger.info(f"Downloading results from: {output_file_url}")
                ok = download_file(output_file_url, args.output_file, api_key)
                return 0 if ok else 9
            elif output_file_url:
                logger.info("Results available at URL (no --output-file provided):")
                print(output_file_url)
                return 0
            else:
                logger.error(f"Job {status} but no output file URL was found.")
                logger.info("Full status response:")
                logger.info(status_data)
                return 10

        if status in ["failed", "cancelled"]:
            logger.error(f"Batch job {status}. Details: {status_data.get('error')}")
            return 11

        logger.info(f"Waiting for {args.poll_interval} seconds before next poll...")
        time.sleep(args.poll_interval)


def cmd_submit_monitor(args: argparse.Namespace) -> int:
    """End-to-end: prepare → submit → monitor → download."""
    rc = cmd_prepare(args)
    if rc != 0:
        return rc
    # Submit and capture batch id from stdout requires a direct call here
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        logger.error("ANTHROPIC_API_KEY environment variable not set.")
        return 5

    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_API_VERSION,
        "content-type": "application/json",
    }
    try:
        with open(args.output_file, "r", encoding="utf-8") as f:
            payload = json.load(f)
        logger.info(f"Submitting batch request for {len(payload.get('requests', []))} topics...")
        response = requests.post(ANTHROPIC_API_URL, headers=headers, data=json.dumps(payload), timeout=120)
        response.raise_for_status()
        response_data = response.json()
        batch_id = response_data.get("id")
        logger.info(f"Batch ID: {batch_id}")
    except Exception as e:
        logger.error(f"Failed to submit: {e}")
        return 6

    # Monitor until completion
    monitor_args = argparse.Namespace(batch_id=batch_id, output_file=args.download_to, poll_interval=args.poll_interval)
    return cmd_monitor(monitor_args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare, submit, monitor, and download Anthropic batch annotations.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Shared options for prepare/submit/submit-monitor
    def add_shared_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("--gene-file", required=True, help="Path to the gene CSV file (columns: Name, Score, RowID)")
        p.add_argument("--enrichment-file", default=None, help="Optional STRING enrichment CSV to include in prompts")
        p.add_argument("--top-enrichment", type=int, default=10, help="Top-N rows per category (KEGG, Process) to include")
        p.add_argument("--max-genes", type=int, default=120, help="Max number of top-loading genes per topic")
        p.add_argument("--output-file", default="anthropic_batch_request.json", help="Path to save the generated batch request JSON")
        p.add_argument("--num-topics", type=int, help="Limit the number of topics to process (testing)")

    # prepare
    p_prepare = subparsers.add_parser("prepare", help="Build the batch request JSON only")
    add_shared_options(p_prepare)
    p_prepare.set_defaults(func=cmd_prepare)

    # submit
    p_submit = subparsers.add_parser("submit", help="Build and submit the batch; prints Batch ID")
    add_shared_options(p_submit)
    p_submit.set_defaults(func=cmd_submit)

    # monitor
    p_monitor = subparsers.add_parser("monitor", help="Monitor a Batch ID and optionally download results")
    p_monitor.add_argument("--batch-id", required=True, help="The ID of the batch job to monitor")
    p_monitor.add_argument("--output-file", help="Where to save the results .jsonl if available")
    p_monitor.add_argument("--poll-interval", type=int, default=30, help="Polling interval in seconds (default: 30)")
    p_monitor.set_defaults(func=cmd_monitor)

    # submit-monitor
    p_submit_monitor = subparsers.add_parser("submit-monitor", help="End-to-end: prepare → submit → monitor → download")
    add_shared_options(p_submit_monitor)
    p_submit_monitor.add_argument("--download-to", required=True, help="Path to save the downloaded results .jsonl")
    p_submit_monitor.add_argument("--poll-interval", type=int, default=30, help="Polling interval in seconds (default: 30)")
    p_submit_monitor.set_defaults(func=cmd_submit_monitor)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    rc = args.func(args)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()


