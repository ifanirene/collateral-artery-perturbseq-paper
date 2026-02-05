# Repository Guidelines
## Instructions
These guidelines apply to the entire repository. Keep changes focused, tested, and easy to review.
When applying rules, reference them in brackets [rule-name] within responses for clarity and traceability.
During you interaction with the user, if you find anything generalizable, especially about a fix to a mistake you made or a correction you received, you should take note in the `learning and memory` section in the `AGENTS.md` file so you will not make the same mistake again. 

## Learning and Memory
- Document reusable project information (library versions, model names, fixes) in this file
- Track recurrent errors and their solutions to prevent repetition
- Token usage logs live under `results/output/logs/token_logs/` (default for `scripts/tools/token_tracker.py`). [Learning and Memory]
- `results/low_usage_program_summary.csv` flags low-usage programs via `is_low_usage`; map `response_id` like `X20` to `program_id` by extracting digits. [Learning and Memory]
- Avoid leaving compatibility symlinks after reorganizations; use real folders only (e.g., `scripts/`, `results/`). [Learning and Memory]
- Installed `hdbscan` (conda-forge) for UMAP clustering; downgraded `scikit-learn` to 1.7.2 to resolve `umap`/`check_array` incompatibility. [Learning and Memory]
- Scripts under `scripts/tools/pipelines/` should use `Path(__file__).resolve().parents[3]` to reach repo root; `parents[1]` only reaches `scripts/tools`. [Learning and Memory]
- For program clustermaps, filter downstream **genes (rows)** by p-value counts when requested, not target columns. [Learning and Memory]
- Harmonizome gene definitions for the union program set are written to `results/perturbation_analysis/harmonizome_gene_definitions_union_filtered.csv`. [Learning and Memory]
- Harmonizome metadata has boilerplate (dates/resources); clean text and remove months/years before keyword clustering. [Learning and Memory]
- When editing files, use the `apply_patch` tool directly instead of invoking it via a shell command. [Learning and Memory]
- Some conda envs here do not include `bin/activate`; call the env's Python directly if activation script is missing. [Learning and Memory]
- Pytest may not resolve `tools.*` imports unless `scripts/` is on `sys.path`; tests can append the repo `scripts/` directory. [Learning and Memory]
- For report search, avoid embedding full annotation text in HTML data attributes; build the search index from `PROGRAMS` data to prevent quote truncation. [Learning and Memory]
- TODO: No new workflows or commands discovered in this run; update when new ones appear. [Learning and Memory]

## Environment Management
- conda activate /Volumes/IF_PHAGE/conda_envs/perturb2 and use it for python tasks.
- conda activate /Volumes/IF_PHAGE/conda_envs/cellchat-r44 and use it for R tasks.
- build conda environments in: `/Volumes/IF_PHAGE/conda_envs/`
- Verify environment with `conda list` before installing packages
- Use mamba for package management when possible
- R env: `/Volumes/IF_PHAGE/conda_envs/cellchat-r44` (R 4.4.x, Seurat 5.1.0, CellChat 1.6.1). [Environment Management]
- Run with env R directly to pick correct libs: `/Volumes/IF_PHAGE/conda_envs/cellchat-r44/bin/Rscript scripts/src/create_filtered_bubble_plot.R`. Avoid `conda run` PATH issues. [Environment Management]

## Project Structure & Module Organization
- `scripts/` — single home for code and CLIs
  - `src/` (core Python source)
    - `analysis/` (reusable analysis modules; prefer entrypoints in `scripts/tools/pipelines/`)
    - `data_processing/` (preprocessing, gRNA comparisons, filtering)
    - `utils/` (I/O helpers, AnnData utilities)
    - `visualization/` (plotting utilities)
  - `tools/` (small CLIs/workflows; e.g., `search_engine.py`, `llm_api.py`, `cluster_annotation/`)
    - `pipelines/` (pipelines and one-off utilities)
      - `analysis/` (analysis driver scripts)
      - `env/` (environment helpers)
      - `scratch/` (temporary scripts)
- `docs/` — documentation and project narratives
  - `proposals/`, `reports/`, `third_party/`
- `notebooks/` — exploratory work (avoid committing large outputs)
- `data/` — inputs and reference datasets
  - `Topic_analysis_K100/`, `Topic_analysis_k100_new/`
  - `reference/`, `external/`, `archive/`, `tmp/`
- `results/` — analysis outputs
  - `perturbation_analysis/`, `markers/`, `programs/`, `metrics/`
- `results/output/` — runtime artifacts and intermediate outputs
  - `logs/`, `llm_batches/`
- `results/figures/` — final figures (PDFs only)
  - `venn/`
- No compatibility symlinks are used; put code under `scripts/` and artifacts under `results/`. [Project Structure & Module Organization]
- `tests/` — pytest suite (see `pytest.ini` for patterns)

## Build, Test, and Development Commands
- conda activate /Volumes/IF_PHAGE/conda_envs/perturb2 and use it for most tasks.
- build conda environments in: `/Volumes/IF_PHAGE/conda_envs/`
- Verify environment with `conda list` before installing packages
- Use mamba for package management when possible
- Run a tool locally:
  - `python scripts/tools/search_engine.py "artery enrichment" --max-results 3`
- Run tests:
  - `pytest -q`  (entire suite)
  - `pytest -k search_engine -q`  (subset)

## Coding Style & Naming Conventions
- PEP 8; 4‑space indents; target ~88‑char lines.
- `snake_case` for modules/functions, `PascalCase` for classes, `UPPER_CASE` for constants.
- Prefer `pathlib`, type hints, and docstrings (Google or NumPy style).
- Avoid hard‑coded absolute paths; accept inputs via args/flags.

## Testing Guidelines
- Framework: pytest (`tests/test_*.py`, classes `Test*`, functions `test_*`).
- Use `tmp_path` and mocks; do not hit network or external APIs in tests.
- Add unit tests with new logic; update tests when behavior changes.

## Commit & Pull Request Guidelines
- Use Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`).
- PRs include: problem/goal, summary of changes, repro command(s), before/after artifacts, and data/config paths touched.
- Ensure `pytest` passes; keep diffs minimal; avoid reformat‑only PRs.

## Security & Configuration
- Secrets via env vars; never commit keys. Required for `scripts/tools/llm_api.py`: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `DEEPSEEK_API_KEY`, `GOOGLE_API_KEY` (optional: `AZURE_OPENAI_API_KEY`, `AZURE_OPENAI_ENDPOINT`).
- Add `.env` to `.gitignore`; when adding variables, provide a minimal `.env.example`.

## Code Annotations
Use standardized annotation blocks for new/modified code:
```python
"""
@description 
This component handles [specific functionality].
It is responsible for [specific responsibilities].

Key features:
- Feature 1: Description
- Feature 2: Description

@dependencies
- DependencyA: Used for X
- DependencyB: Used for Y

@examples
- Example usage quick
- Example usage detailed parameter setting
"""
```

## Agent-Specific Instructions
- Always execute scripts after writing/editing them!
- Do not rename or move files without justification.
- Update docs/tests when changing public APIs.
- Always save plots as pdf only, unless specified otherwise. For editable text in PDF figures, set mpl.rcParams['pdf.fonttype'] = 42 to use TrueType fonts instead of Type 3 fonts.
