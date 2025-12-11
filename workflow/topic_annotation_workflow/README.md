# Topic Annotation Workflow

This workflow automates building program gene lists, running STRING enrichment, submitting topics for LLM-based annotation via Anthropic Batch, and generating final reports.

---

## Prerequisites

- Conda env with `pandas`, `requests`, `python-dotenv`, `markdown`, `matplotlib`
- Environment variable `ANTHROPIC_API_KEY` available (e.g., via `.env` in project root)

---

## Unified CLIs

Two consolidated tools replace the older multiple scripts:

- `tools/topic_annotation_workflow/01_genes_to_string_enrichment.py` — extract top genes and run STRING enrichment
- `tools/topic_annotation_workflow/02_submit_and_monitor_batch.py` — prepare, submit, monitor, and download Anthropic batch jobs

The legacy scripts `annotate_cnmf_programs_string.py`, `run_string_enrichment.py`, `01_submit_batch_annotation.py`, and `check_and_download_batch.py` are now superseded by the above.

---

## 1) Extract genes and run STRING enrichment

### Extract top-N genes per program
```bash
python tools/topic_annotation_workflow/01_genes_to_string_enrichment.py extract \
  --input data/Hypoxia/p2_HnN_loading_gene_k30.csv \
  --n-top 100 \
  --json-out output/topic_annotations/hnn_k30_top100_genes.json \
  --csv-out output/topic_annotations/hnn_k30_top100_genes_overview.csv
```

### Run STRING enrichment (and optional figures)
```bash
python tools/topic_annotation_workflow/01_genes_to_string_enrichment.py enrich \
  --genes-json output/topic_annotations/hnn_k30_top100_genes.json \
  --species 10090 \
  --out-csv-full output/topic_annotations/hnn_k30_string_full.csv \
  --out-csv-filtered output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
  --figures-dir output/topic_annotations/enrichment_figures
```

### End-to-end extract → enrich
```bash
python tools/topic_annotation_workflow/01_genes_to_string_enrichment.py all \
  --input data/Hypoxia/p2_HnN_loading_gene_k30.csv \
  --n-top 100 \
  --json-out output/topic_annotations/hnn_k30_top100_genes.json \
  --csv-out output/topic_annotations/hnn_k30_top100_genes_overview.csv \
  --species 10090 \
  --out-csv-full output/topic_annotations/hnn_k30_string_full.csv \
  --out-csv-filtered output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
  --figures-dir output/topic_annotations/enrichment_figures
```

---

## 2) Prepare, submit, monitor, and download Anthropic batch

### Prepare batch request only
```bash
python tools/topic_annotation_workflow/02_submit_and_monitor_batch.py prepare \
  --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
  --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
  --top-enrichment 10 \
  --max-genes 120 \
  --num-topics 2 \
  --output-file output/topic_annotations/hypoxia_k30/anthropic_batch_request_with_enrichment_test.json
```

### Submit now (prints Batch ID)
```bash
python tools/topic_annotation_workflow/02_submit_and_monitor_batch.py submit \
  --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
  --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv
```

### Monitor an existing Batch ID and download results
```bash
python tools/topic_annotation_workflow/02_submit_and_monitor_batch.py monitor \
  --batch-id <BATCH_ID> \
  --output-file output/topic_annotations/hypoxia_k30/batch_results.jsonl \
  --poll-interval 30
```

### End-to-end: submit → monitor → download
```bash
python tools/topic_annotation_workflow/02_submit_and_monitor_batch.py submit-monitor \
  --gene-file data/Hypoxia/p2_HnN_loading_gene_k30.csv \
  --enrichment-file output/topic_annotations/hnn_k30_string_process_kegg_filtered_bg500.csv \
  --download-to output/topic_annotations/hypoxia_k30/batch_results.jsonl \
  --poll-interval 30
```

Notes:
- Enrichment file expected columns: `program_id, category (KEGG|Process), description, fdr, inputGenes`.
- Defaults: `top_enrichment=10`, `max_genes=120`, model is `claude-3-7-sonnet-20250219` with `max_tokens=3500`.

---

## 3) Parse results and build summary (merged)

### Parse API results to markdown files AND generate summary CSV
```bash
python tools/topic_annotation_workflow/03_parse_and_summarize.py \
  --results-jsonl <results.jsonl> \
  --markdown-dir <markdown_dir> \
  --summary-csv <output_summary.csv> \
  --gene-loading-file <genes.csv>
```

Notes:
- Use `--no-summary` to only parse JSONL → markdown.
- Use `--no-parse` with `--markdown-dir` and `--summary-csv` to only build the summary CSV.

### Create final HTML report
```bash
python tools/topic_annotation_workflow/04_generate_html_report.py <summary.csv> <markdown_dir> <final_report.html>
```
