# Collateral Artery Perturb-seq: Code Repository

This repository contains the computational analysis code accompanying the manuscript:

**"A cross-species Perturb-seq platform identifies artery repressors that govern collateral artery formation"**

All computational analyses were performed using the custom scripts provided here. Raw sequencing data are deposited in GEO (accession pending).

> **Note on completeness**: Several analyses described in the manuscript (CellOracle GRN construction, scFEA, CellRank trajectory analysis, hypoxia scRNA-seq preprocessing, within-screen attenuation analysis, and endothelial subtype program enrichment) are not yet represented by a released script. A full inventory of missing scripts, missing software version information, and other details required for complete reproduction is provided in [Additional_Requested_Information.md](Additional_Requested_Information.md).

---

## Study Overview

Collateral arteries are natural arterial bypasses that limit ischemic tissue damage, yet no therapies exist to promote their growth. We exploited an evolutionary outlier—the guinea pig, which develops exceptionally dense collateral networks across multiple organs—to discover the genes that govern collateral artery formation. Using integrated single-cell RNA sequencing (scRNA-seq) and ATAC-seq across matched embryonic mouse and guinea pig hearts and brains, we identified endothelial gene expression programs that diverge between species. Rather than finding that canonical pro-arterial signals (VEGF, NOTCH) are globally enhanced in guinea pigs, we discovered that pathways normally active in mice—hypoxia response, TCA-cycle metabolism, unfolded protein response, and WNT signaling—are selectively reduced in guinea pig endothelial cells.

To test whether these cross-species differences causally influence arterial identity, we developed an in vivo CRISPRi Perturb-seq platform in the neonatal mouse brain. A pooled library of 162 guide RNAs targeting candidate genes was delivered to brain endothelial cells via AAV, and perturbation-resolved transcriptomes were captured at postnatal day 9. The screen identified 45 artery repressors, including multiple genes in the WNT receptor, HIF transcription factor, TCA-cycle enzyme, and protein chaperone families, all of which were downregulated in guinea pig endothelium relative to mouse. Targeted in vivo inhibition of selected repressors—Frizzled receptors, TCA-cycle enzymes, and ER chaperones—increased pial collateral artery abundance in neonatal mice. Consensus non-negative matrix factorization (cNMF) of the Perturb-seq transcriptomes resolved 100 co-regulated gene programs and revealed that these repressors converge on VEGF-dependent and neuronal-guidance programs in angiogenic tip and pre-artery tip endothelial cell states, representing a previously underappreciated bottleneck in the tip-to-artery transition. Environmental hypoxia independently recapitulated the same program shifts and reduced collateral formation, and the same programs were differentially active between guinea pig and mouse brain endothelial cells, providing orthogonal support for the proposed mechanism.

---

## Repository Structure

```
analysis/
├── 01_cross_species_atlas/          # Cross-species scRNA-seq/ATAC-seq integration and spatial validation
│   ├── 01_heart_ec_crossspecies_integration.R
│   ├── 02_brain_ec_crossspecies_integration.R
│   ├── 03_heart_ec_cellchat_ligand_receptor.R
│   ├── 04_heart_xenium_ec_pseudobulk_validation.R
│   ├── 11_heart_ec_atac_signac_cicero_preprocessing.R
│   └── 12_heart_brain_ec_hallmark_pathway_comparison.py
│
├── 02_perturbseq_screen/            # In vivo Perturb-seq preprocessing and SCEPTRE artery-score screen
│   ├── 05_perturbseq_endothelial_preprocessing.R
│   ├── 13_endothelial_state_ucell_scoring.R
│   ├── 14_endothelial_state_sceptre_artery_score.R
│   └── 15_perturbseq_gene_level_sceptre_transcriptome.R
│
├── 03_gene_program_analysis/        # cNMF program discovery, SCEPTRE testing, and cross-dataset projection
│   ├── 06_perturbseq_cnmf_k100_factorization.py
│   ├── 07_perturbseq_cnmf_discovery_notebook.ipynb
│   ├── 08_perturbseq_cnmf_program_sceptre.R
│   ├── 09_project_programs_onto_hypoxia_normoxia.py
│   └── 10_project_programs_onto_crossspecies_brain_ec.py
│
genome_reference/                    # Guinea pig genome reference build scripts and configuration
│   ├── Genome_reference.md
│   ├── CavPor4.config
│   ├── BSgenome.CPorcellus.NCBI.cavPor4-seed
│   └── reproduce_guinea_pig_genome_reference_build.R
```

### Script Overview

| # | Script | Analysis section | What it does |
|---|---|---|---|
| 01 | `01_heart_ec_crossspecies_integration.R` | Cross-species analysis | Heart scRNA-seq QC, ortholog remapping, integration, EC annotation (Fig 2a–d) |
| 02 | `02_brain_ec_crossspecies_integration.R` | Cross-species analysis | Brain scRNA-seq integration and EC annotation (Fig 2h, 2k) |
| 03 | `03_heart_ec_cellchat_ligand_receptor.R` | Cross-species analysis | CellChat ligand-receptor analysis, receptor candidate nomination (Fig 2g) |
| 04 | `04_heart_xenium_ec_pseudobulk_validation.R` | Cross-species analysis | Xenium spatial pseudobulk DE for spatial validation (Fig 2e–f, 2j) |
| 05 | `05_perturbseq_endothelial_preprocessing.R` | In vivo Perturb-seq | Perturb-seq QC, cluster pruning, guide filtering, SCEPTRE/cNMF export (Fig 3) |
| 06 | `06_perturbseq_cnmf_k100_factorization.py` | Gene program analysis | cNMF k=100 factorization (command-line, parallelizable) |
| 07 | `07_perturbseq_cnmf_discovery_notebook.ipynb` | Gene program analysis | Interactive cNMF exploration notebook (rank sweep, parameter rationale) |
| 08 | `08_perturbseq_cnmf_program_sceptre.R` | Gene program analysis | SCEPTRE testing of cNMF program usage by perturbation (Fig 5e–h) |
| 09 | `09_project_programs_onto_hypoxia_normoxia.py` | Gene program analysis | STARcat projection onto hypoxia/normoxia dataset (Fig 5m–q) |
| 10 | `10_project_programs_onto_crossspecies_brain_ec.py` | Gene program analysis | STARcat projection onto cross-species brain EC atlas (Fig 5q bottom) |
| 11 | `11_heart_ec_atac_signac_cicero_preprocessing.R` | Cross-species analysis | scATAC-seq QC, depth-matching, Cicero coaccessibility (Fig S1 ATAC panels) |
| 12 | `12_heart_brain_ec_hallmark_pathway_comparison.py` | Cross-species analysis | decoupler AUCell on MSigDB Hallmark gene sets, cross-organ concordance (Fig 2h) |
| 13 | `13_endothelial_state_ucell_scoring.R` | In vivo Perturb-seq | UCell scoring of endothelial-state marker genes; exports score matrix for SCEPTRE |
| 14 | `14_endothelial_state_sceptre_artery_score.R` | In vivo Perturb-seq | Primary artery-score SCEPTRE screen; identifies 45 artery repressors (Fig 3g–h, 5a–b) |
| 15 | `15_perturbseq_gene_level_sceptre_transcriptome.R` | In vivo Perturb-seq | Transcriptome-wide gene-level SCEPTRE testing (all cells and per-cell-type) |

---

## Analysis Pipeline

The analyses are organized into three sequential modules. Within each module, scripts are numbered in execution order. Data produced by earlier scripts serve as inputs to later ones.

### Module 1 — Cross-Species Endothelial Cell Atlas (`01_cross_species_atlas/`)

This module builds the cross-species single-cell atlas used to nominate candidate genes and to validate spatial transcriptomics signals. The two integration scripts (01 and 02) can be run independently of one another; the CellChat (03) and Xenium (04) scripts depend on the labeled endothelial objects produced by the heart integration script.

**Script 01 — Heart endothelial cell cross-species integration** (`01_heart_ec_crossspecies_integration.R`)

This script reproduces the preprocessing, ortholog remapping, integration, and subtype annotation of embryonic heart endothelial cells from three mouse stages (E13.5, E15.5, E17.5) and three matched guinea pig stages (GD25, GD32, GD35). Starting from filtered Cell Ranger gene-expression matrices, it performs sample-level quality control and doublet scoring with scDblFinder, remaps guinea pig gene symbols to mouse ortholog names using a BioMart-derived one-to-one ortholog table, and integrates the six samples into a shared endothelial space using Seurat v4 canonical correlation analysis anchoring. After whole-heart integration, endothelial clusters are isolated and reintegrated at higher resolution. The final labeled object distinguishes large artery, arteriole, and multiple capillary and venous states using canonical markers (*Gja5*, *Gja4*, *Cxcr4*, *Nr2f2*, *Apln*) and exports the annotated metadata alongside QC and UMAP summary plots. This script produces the heart endothelial atlas underlying Figure 2a–d and Figures S2–S3.

**Script 02 — Brain endothelial cell cross-species integration** (`02_brain_ec_crossspecies_integration.R`)

This script applies the same preprocessing and integration strategy to late-stage brain endothelial cells: two mouse samples (E18 and P0) and two guinea pig GD35 replicates. Because brain samples were FACS-sorted for CD31+ cells before sequencing, the filtering criteria differ from the heart workflow—higher gene and UMI thresholds are applied, and guinea pig cells are additionally screened by ribosomal content to exclude low-quality nuclei. After whole-brain integration and endothelial cluster isolation, a higher-stringency subset is reembedded and annotated into arterial, pre-artery, capillary, cycling, and venous states. This brain atlas supports the cross-species differential expression analyses shown in Figure 2h, 2k, and Figure S3a–c.

**Script 03 — CellChat ligand-receptor comparison of heart endothelial cells** (`03_heart_ec_cellchat_ligand_receptor.R`)

This script uses CellChat v2.1.2 to infer and compare intercellular signaling into endothelial cells in mouse versus guinea pig embryonic hearts. One CellChat object is constructed per species from the annotated, integrated heart atlas, using annotated cardiac cell classes as sender and receiver groups. Incoming pathway strengths directed toward endothelial cells are aggregated across all sender populations and compared between species, testing whether canonical pro-arterial pathways (VEGF, NOTCH, TGFb, WNT, CXCL) are globally elevated in guinea pigs. A second analysis maps cross-species differential expression onto the merged CellChat network to rank receptors by the strength and directionality of their species-biased incoming interactions, identifying guinea pig-enriched and mouse-enriched receptor candidates for the Perturb-seq library. Outputs include species-specific and merged CellChat objects, pathway-level summary tables, and receptor candidate ranking tables (Figure 2g, Figure S3d–e). The top 50 differentially expressed receptors contributed to the Perturb-seq screen library.

**Script 04 — Xenium spatial transcriptomics endothelial pseudobulk validation** (`04_heart_xenium_ec_pseudobulk_validation.R`)

This script processes spatial gene expression data from a homology-aware Xenium panel applied to one embryonic heart per species (mouse E18, guinea pig GD40). It discovers and loads the four vendor Xenium region outputs, applies cell-level quality filters (100–2,000 transcripts, at least 50 detected genes), and assigns spatial replicate units by partitioning cell centroids into four spatial clusters per region using k-means. The four region outputs are then integrated using Seurat v5 sketch-based CCA integration, followed by cluster annotation and endothelial cell isolation using either manual cluster labels or an endothelial marker fallback mask. Raw Xenium counts within each spatial replicate unit are aggregated into endothelial pseudobulk matrices and tested for differential expression between species using DESeq2 with a region-adjusted design. This spatial validation confirms cross-species differences in arterial markers (*Gja5*), pro-arterial ligands (*Vegfa*, *Notch1*), TCA-cycle enzymes (*Idh2*, *Mdh2*, *Ogdh*), and protein chaperones (*Hsp90ab1*, *Hspa5*) identified in the single-cell analyses (Figure 2e–f, 2j).

**Script 11 — Heart endothelial cell scATAC-seq QC and Cicero coaccessibility** (`11_heart_ec_atac_signac_cicero_preprocessing.R`)

This script runs the complete scATAC-seq preprocessing pipeline for the cross-species heart Multiome data. Starting from Cell Ranger ARC outputs, it constructs Signac ChromatinAssay objects for both mouse (mm10) and guinea pig (mCavPor4.1) using species-appropriate gene annotations derived from Ensembl and the custom NCBI GFF3 respectively. Per-cell quality metrics—nucleosome signal and TSS enrichment—are computed and cells are filtered to retain those with 2,000–25,000 in-peak fragments, nucleosome signal < 2.5, and TSS enrichment > 2. To enable a balanced cross-species comparison at the chromatin level, mouse cells are depth-matched to guinea pig by binomial thinning scaled to the ratio of per-species median in-peak fragment depths. Both depth-matched and unmatched mouse ATAC data as well as guinea pig ATAC data are then converted to Monocle3 CDS objects with LSI dimensionality reduction and UMAP embedding for Cicero input. Cicero is run independently in each species (k=50 for guinea pig, k=100 for depth-matched mouse) and the resulting peak-to-peak coaccessibility connections are exported in underscore-delimited CellOracle format (peak coordinates converted from dot notation). Random seeds are set explicitly for Cicero (2013) and binomial thinning (42). This script depends on the cell barcodes and ATAC metadata produced by Script 01 and produces the coaccessibility network used as the base GRN prior for CellOracle (Figure S1 ATAC panels). Run in the CellOracle conda environment (Python 3.8.18 / R 4.2.0).

**Script 12 — Cross-species endothelial Hallmark pathway comparison** (`12_heart_brain_ec_hallmark_pathway_comparison.py`)

This Python script computes and compares pathway activity between mouse and guinea pig endothelial cells in both heart and brain using the decoupler AUCell method. It loads preprocessed heart and brain endothelial AnnData objects, re-normalizes counts to 10,000 per cell, and runs AUCell scoring against all 50 MSigDB Hallmark gene sets (v2024.1). Per-pathway activity scores are compared between species using a Wilcoxon rank-sum test per organ, and only pathways passing FDR < 0.01 and |log2FC| ≥ 0.25 in both heart and brain are reported as concordantly cross-organ differential. Summary tables of per-organ DE results and the merged concordant pathway list are exported as CSV/TSV files, and a combined dot-plot figure is saved as PDF. Environment: Python 3.12.9, scanpy 1.11.3, decoupler 2.1.1, anndata 0.11.4, pandas 2.2.3, numpy 2.2.4, scipy 1.15.2, matplotlib 3.10.1, seaborn 0.13.2, adjustText 1.3.0 (Figure 2h).

---

### Module 2 — In Vivo Perturb-seq Screen (`02_perturbseq_screen/`)

**Script 05 — Perturb-seq endothelial preprocessing** (`05_perturbseq_endothelial_preprocessing.R`)

This script reconstructs the manuscript-facing preprocessing path for the neonatal mouse brain CRISPRi Perturb-seq screen. Starting from Cell Ranger outputs for ten sequencing batches, it combines resequenced RNA counts with original CRISPR guide capture matrices by intersecting shared barcodes. After RNA-level quality control (2,000–7,500 genes, 3,000–30,000 UMIs, less than 10% mitochondrial reads), the merged object is clustered and non-endothelial populations are removed through three rounds of manual cluster pruning to retain endothelial, perivascular, and screen-relevant populations. Doublet detection is then run in batch-aware mode using scDblFinder. Guide multiplicity is estimated from UMI thresholding, and cells with more than 15 assigned guides or zero guide UMIs are excluded. The final filtered RNA and gRNA count matrices, together with per-cell covariates (batch, cell type, mitochondrial fraction), are exported as an RData bundle consumed by all downstream SCEPTRE analyses (Script 08) and an h5ad file consumed by cNMF (Scripts 06–07). This script supports Figure 3b–e and Figure S4.

**Script 13 — Endothelial-state UCell marker scoring** (`13_endothelial_state_ucell_scoring.R`)

This script derives per-cell UCell scores for each endothelial subtype from the filtered Perturb-seq Seurat object produced by Script 05. The six endothelial states (artery, pre-artery, tip, capillary, venous, and arterio-venous) are grouped by combining artery_1 and artery_2 clusters from the upstream annotation. For each state, marker genes are identified using a ROC-based `FindAllMarkers` test and ranked by AUC. The top 50 markers are retained for the three artery-proximal states (large artery, arteriole, and pre-artery); the top 30 markers are retained for the remaining states. UCell `AddModuleScore_UCell` is then run on the full filtered object using these state-specific marker sets. The resulting per-cell score table and sparse score matrix are exported for input to Script 14. This script is a prerequisite for the primary artery-score screen.

**Script 14 — Primary artery-score SCEPTRE screen** (`14_endothelial_state_sceptre_artery_score.R`)

This script runs the central SCEPTRE screen that identifies 45 artery repressors from the Perturb-seq data. It loads the guide count matrix, screen inputs, singlet covariate table, and UCell score matrix from Scripts 05 and 13, and tests each perturbation against each of the six endothelial-state UCell scores using SCEPTRE v0.1.0 in high-MOI mode with mixture-model guide assignment. The regression formula includes log-transformed per-cell guide count and UMI covariates together with batch and broad cell-type indicators (`~ log(grna_n_nonzero) + log(grna_n_umis) + batch + cell_type_l1`). Output discovery tables record which perturbations significantly elevate or reduce each endothelial-state score and are used to define artery repressors as those that increase the artery or pre-artery score after CRISPRi knockdown (Figures 3g–h, 5a–b, S4k). This is the primary screen analysis.

**Script 15 — Transcriptome-wide gene-level SCEPTRE testing** (`15_perturbseq_gene_level_sceptre_transcriptome.R`)

This command-line R script runs a transcriptome-wide perturbation association test using SCEPTRE, testing each guide RNA against each expressed gene across the full Perturb-seq endothelial dataset. It accepts `--gene-guide-input`, `--technical-covariates`, `--guide-target-annotation`, and `--output-root` arguments for flexible submission on a compute cluster. Guide assignment uses a UMI threshold of 10 (thresholding rather than mixture-model assignment). Both an all-cells analysis and per-cell-type analyses (for subtypes with at least 200 cells) are run. A `--dry-run` flag supports validation without executing SCEPTRE, and an `--all-cells-only` flag restricts the run to the full-dataset analysis. Outputs are written to separate subdirectories for all-cells and per-cell-type results.

---

### Module 3 — Gene Program Discovery and Cross-Dataset Projection (`03_gene_program_analysis/`)

This module discovers co-regulated gene expression programs from the Perturb-seq transcriptomes using consensus NMF (cNMF), tests how each CRISPRi perturbation affects program usage using SCEPTRE, and projects the resulting program reference onto two independent datasets—a neonatal hypoxia experiment and the cross-species brain EC atlas. Scripts 06 and 07 are alternative entry points for the same cNMF factorization workflow (command-line versus interactive notebook); Script 08 depends on the outputs of Script 05 and Scripts 06/07; Scripts 09 and 10 both depend on the STARcat reference generated by Scripts 06/07.

**Script 06 — cNMF factorization of Perturb-seq endothelial cells** (`06_perturbseq_cnmf_k100_factorization.py`)

This command-line Python script runs the retained 100-program cNMF factorization on the preprocessed Perturb-seq endothelial AnnData (raw UMI counts). It provides four sequential stages: `prepare` (restores raw counts from the AnnData `raw` slot and writes the cNMF job ledger), `factorize-worker` (runs one factorization shard, intended to be parallelized across 12 workers), `combine` (merges per-worker spectra), and `consensus` (clusters replicate spectra at a local-density threshold of 0.2 and writes the final consensus outputs). The retained settings—k=100, 200 factorization replicates, 10,000 overdispersed genes—were chosen after a rank sweep documented in Script 07. An `all-sequential` fallback stage is available for single-machine runs. This script produces the consensus spectra, per-cell usage matrix, and STARcat reference file that serve as inputs to Scripts 08–10.

```bash
# Prepare the cNMF job ledger
python analysis/03_gene_program_analysis/06_perturbseq_cnmf_k100_factorization.py prepare \
  --source-h5ad data/final_pool_endothelial_cells.h5ad

# Factorize in parallel (recommended; adjust --total-workers to match available CPUs)
parallel python analysis/03_gene_program_analysis/06_perturbseq_cnmf_k100_factorization.py \
  factorize-worker --total-workers 12 --worker-index {} ::: $(seq 0 11)

# Merge and build consensus
python analysis/03_gene_program_analysis/06_perturbseq_cnmf_k100_factorization.py combine
python analysis/03_gene_program_analysis/06_perturbseq_cnmf_k100_factorization.py consensus
```

**Script 07 — cNMF exploratory and production workflow notebook** (`07_perturbseq_cnmf_discovery_notebook.ipynb`)

This Jupyter notebook documents the interactive cNMF exploration that preceded the production run. It illustrates the exploratory run at k=15 used to verify that batch-associated signal was absent without correction, the rank sweep across k=10–60 used to select the final program number, and the production run at k=100 using the same settings as Script 06. The notebook is provided as a companion to Script 06 to clarify the exploratory reasoning behind parameter choices. The final production outputs loaded at the end of the notebook are the same objects generated by Script 06.

**Gene program annotation** (not a script in this repository)

Before program-level perturbation testing, each of the 100 cNMF programs was assigned a biologically interpretable label using the ProgExplorer annotation tool ([https://github.com/ifanirene/ProgExplorer](https://github.com/ifanirene/ProgExplorer)). ProgExplorer integrates top gene loadings, STRING Gene Ontology and KEGG enrichment against the Mus musculus background (STRING v12.0), endothelial subtype-enrichment summaries, and perturbation-response results, and applies a structured LLM-backed annotation prompt to assign concise program names, functional module descriptions, and regulator hypotheses. Refer to the ProgExplorer repository for the annotation pipeline code and documentation. Annotated program outputs are referenced in the manuscript supplementary tables.

**Script 08 — SCEPTRE program-level perturbation testing** (`08_perturbseq_cnmf_program_sceptre.R`)

This script tests whether each CRISPRi perturbation significantly alters the usage of each of the 100 cNMF gene programs discovered in Scripts 06–07. Per-cell program usage values are converted into count-like responses by row-normalizing each cell's usage vector to sum to one and multiplying by the cell's total UMI count before rounding, making the SCEPTRE conditional-resampling model applicable. The script aligns cells across the program usage matrix, the gene-level SCEPTRE input bundle from Script 05, and the guide annotation table, then runs SCEPTRE in high-MOI mode with batch and library-size covariates. Both all-endothelial-cell and per-cell-type analyses are supported through a configuration flag. The output discovery tables identify which perturbations increase or decrease each program and are used to characterize how artery repressors act through WNT, HIF-glycolysis, VEGF-tip-cell, and neuronal-guidance programs (Figure 5e–h, Figure S5b–j).

**Script 09 — STARcat projection onto hypoxia/normoxia brain endothelial cells** (`09_project_programs_onto_hypoxia_normoxia.py`)

This script uses STARcat to project the 100 cNMF gene programs from the Perturb-seq screen onto an independently collected dataset of neonatal mouse brain endothelial cells exposed to either normoxia (21% O₂) or mild hypoxia (11% O₂) for seven days. The projection holds the cNMF consensus spectra fixed and estimates per-cell usage scores by non-negative least squares (NNLS), column-standardizing the query counts before fitting. Row-normalized usage scores (summing to one per cell) are saved alongside batch and cell-type annotations for downstream differential-usage analysis. This projection tests whether the HIF-glycolysis, canonical WNT, and tip-cell programs identified as targets of artery repressors in the Perturb-seq screen are also induced by environmental oxygen reduction, and whether artery-associated programs are correspondingly diminished (Figure 5m–q).

**Script 10 — STARcat projection onto cross-species brain endothelial cells** (`10_project_programs_onto_crossspecies_brain_ec.py`)

This script applies the same STARcat projection framework to the integrated cross-species brain endothelial cell dataset from Module 1, enabling direct comparison of program activity between mouse and guinea pig brain endothelial cells. Because the query includes guinea pig cells whose gene symbols were mapped to mouse orthologs during upstream integration, only the 6,704 reference genes present in the shared ortholog space are used for fitting (compared with approximately 9,755 for a pure mouse query). The resulting per-cell usage matrix is saved alongside species and cell-type metadata for stratified program comparisons. This analysis tests whether the artery-associated programs elevated in guinea pig ECs and the hypoxia, WNT, and tip-cell programs reduced in guinea pig ECs mirror the program shifts observed in the Perturb-seq screen, closing the loop between the cross-species atlas and the in vivo perturbation results (Figure 5q, bottom row).

---

## Guinea Pig Genome Reference (`genome_reference/`)

The cross-species atlas and scATAC-seq analyses require two custom genome reference resources for *Cavia porcellus* that are not available from 10x Genomics or Bioconductor. Both are derived from the NCBI mCavPor4.1 assembly (GCF_034190915.1, RS_2024_02, February 2023). Build scripts and configuration files are in `genome_reference/`; large binary outputs (genome FASTA, Cell Ranger ARC reference tarball, R package tarball) are available upon request.

> **Version note**: All references to the guinea pig genome in this study use the NCBI assembly name **mCavPor4.1**. Internal build-iteration labels (CavPor4.2, CavPor4.3) that appear in R package names are implementation details only; the underlying genome sequence is unchanged from GCF_034190915.1.

### Custom Cell Ranger ARC Multiome Reference

The Cell Ranger ARC reference was built with `cellranger-arc mkref` v2.0.2. The source GTF was filtered to retain only protein-coding and lncRNA genes (28,680 of 37,419 total) to reduce spurious multi-mapping; pseudogenes and small RNA genes were excluded. A custom motif file derived from JASPAR2024 Core vertebrates (non-redundant) replaced the standard Cell Ranger ARC motif library. The 86 largest NCBI scaffold accessions (NW_026947484.1 through NW_026948059.1) were designated as primary contigs; there is no mitochondrial contig in this assembly release.

```bash
cellranger-arc mkref --config=CavPor4.config --nthreads=32 --memgb=256
```

- Config file: `genome_reference/CavPor4.config`
- Full build documentation (input files, biotype filter rationale, primary contig list): `genome_reference/Genome_reference.md`

### Custom BSgenome Data Package

The BSgenome package (`BSgenome.CPorcellus.NCBI.CavPor4.3`, v4.0) provides guinea pig genome sequences for Signac chromatin assay construction and motif scanning in R/Bioconductor. It was built by converting the NCBI FASTA to 2-bit format with `faToTwoBit` (UCSC utilities, April 2021) and packaging with `BSgenome::forgeBSgenomeDataPkg()` in R 4.2.0. The package exposes 86 sequences with NCBI accession identifiers; in Signac and Cicero, underscores in seqlevels are replaced with hyphens before peak-matrix construction.

- Seed file: `genome_reference/BSgenome.CPorcellus.NCBI.cavPor4-seed`
- Build script (also covers TxDb and OrgDb construction): `genome_reference/reproduce_guinea_pig_genome_reference_build.R`
- Install from provided tarball: `R CMD INSTALL BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz`

---

## Data Availability

Raw sequencing data (FASTQ files) and processed count matrices for all experiments will be deposited in the NCBI Gene Expression Omnibus (GEO) upon publication. The following datasets are included:

| Dataset | Species | Assay | Stages / Conditions |
|---|---|---|---|
| Cross-species heart atlas | Mouse + guinea pig | scRNA-seq + scATAC-seq (10x Multiome) | Mouse E13.5, E15.5, E17.5; guinea pig GD25, GD32, GD35 |
| Cross-species brain atlas | Mouse + guinea pig | scRNA-seq (10x 3' v3.1), CD31+ FACS-sorted | Mouse E18, P0; guinea pig GD35 |
| Xenium spatial transcriptomics | Mouse + guinea pig | Xenium (299-gene homology-aware panel) | Mouse E18; guinea pig GD40 |
| In vivo Perturb-seq screen | Mouse | scRNA-seq + direct guide capture (10x GEM-X 3' v4) | Neonatal P1 AAV injection, P9 harvest; 10 replicates |
| Hypoxia vs. normoxia | Mouse | scRNA-seq (10x 3'), CD31+ FACS-sorted | P1–P8 exposure, 21% O₂ vs. 11% O₂ |

---

## Software Requirements

### Prerequisite — Cell Ranger (primary sequencing data processing)

| Software | Version | Purpose |
|---|---|---|
| Cell Ranger | 7.1.0 | Demultiplexing, alignment, and UMI counting for all scRNA-seq libraries (Perturb-seq, brain atlases, hypoxia) |
| Cell Ranger ARC | 2.0.2 | Joint RNA + ATAC counting for cross-species heart Multiome libraries |

### R Environment — Seurat 4 (heart atlas, brain atlas, CellChat, Perturb-seq preprocessing, artery-score screen)

Scripts 01–03 and 05 were developed and validated in R 4.2.0 on a Linux (CentOS 7) system using the packages listed below. An environment snapshot is recorded in the header of each script.

| Package | Version | Purpose |
|---|---|---|
| Seurat | 4.4.0 | Single-cell preprocessing, integration, clustering, and visualization |
| SeuratObject | 5.0.1 | Seurat object infrastructure |
| scDblFinder | 1.17.2 / 1.14.0 | Cluster-aware doublet detection (v1.17.2 for cross-species atlas; v1.14.0 for Perturb-seq preprocessing) |
| UCell | 2.2.0 | Rank-based gene signature scoring (artery-state scores) |
| BiocParallel | 1.32.5 | Parallel execution for doublet scoring |
| gprofiler2 | 0.2.3 | Ortholog mapping for cell-cycle gene sets |
| CellChat | 2.1.2 | Ligand-receptor communication inference (Script 03 only) |
| sceptre | 0.1.0 / 0.10.0 | CRISPRi perturbation testing (v0.1.0 for gene-level artery-score screen; v0.10.0 for program-level testing in Script 08) |
| Signac | 1.10.0 | scATAC-seq QC and chromatin assay construction (Script 11) |
| ggplot2 | 3.5.0 | Visualization |
| patchwork | 1.2.0 | Plot composition |
| dplyr | 1.1.4 | Table manipulation |

### R Environment — Seurat 5 (Xenium spatial workflow)

Script 04 requires Seurat v5 for sketch-based integration of Xenium region outputs:

| Package | Version | Purpose |
|---|---|---|
| Seurat | 5.3.0 | Xenium loading, sketch integration, and cluster projection |
| DESeq2 | 1.42.0 | Pseudobulk differential expression |
| jsonlite | 2.0.0 | Parsing Xenium experiment metadata |
| data.table | 1.17.8 | Reading `cells.csv.gz` per-cell metadata |

### Python Environment — cNMF and STARcat (Scripts 06, 07, 09, 10)

| Package | Version | Purpose |
|---|---|---|
| Python | ≥ 3.10 | Runtime |
| scanpy | 1.10.1 | AnnData handling and preprocessing (used for cNMF raw-count export) |
| anndata | 0.10.7 | AnnData I/O |
| cnmf | 1.7.0 | Consensus NMF factorization |
| starcatpy | 1.0.9 | NNLS-based projection onto cNMF reference spectra |
| numpy | | Numerical operations |
| pandas | | DataFrame handling |
| scikit-learn | 1.0.2 | NNLS solver (used internally by STARcat) |

### Python Environment — Pathway analysis (Script 12; cross-species Hallmark scoring)

Script 12 runs in a separate Python 3.12.9 environment (the `perturb2` conda environment used for this project):

| Package | Version | Purpose |
|---|---|---|
| Python | 3.12.9 | Runtime |
| scanpy | 1.11.3 | Renormalization of count matrices before scoring |
| decoupler | 2.1.1 | AUCell scoring of MSigDB Hallmark gene sets (v2024.1) |
| anndata | 0.11.4 | AnnData I/O |
| pandas | 2.2.3 | DataFrame handling |
| numpy | 2.2.4 | Numerical operations |
| scipy | 1.15.2 | Statistical testing |
| matplotlib | 3.10.1 | Figure rendering |
| seaborn | 0.13.2 | Plot styling |
| adjustText | 1.3.0 | Label de-overlapping |

### Python Environment — CellOracle (GRN construction and TF perturbation scoring)

Used for the CellOracle analysis (missing script; see [Additional_Requested_Information.md](Additional_Requested_Information.md)). This is a **separate environment** from the cNMF/STARcat environment:

| Package | Version | Purpose |
|---|---|---|
| Python | 3.8.18 | Runtime |
| CellOracle | 0.18.0 | GRN construction and in silico TF knockout simulations |
| scanpy | 1.9.6 | Single-cell operations within CellOracle pipeline |
| pybedtools | 0.9.1 | Guinea pig TSS annotation construction |
| gimmemotifs | 0.17.0 | Motif scanning for base GRN (motif database: gimme.vertebrate.v5.0) |
| Monocle3 | 1.3.1 | Latent semantic indexing backend for Cicero coaccessibility |
| Cicero | 1.3.9 | Peak-to-gene coaccessibility (cole-trapnell-lab/monocle3 branch) |

### Python Environment — CellRank (trajectory and fate probability analysis)

Used for CellRank trajectory analysis (script not yet released; see [Additional_Requested_Information.md](Additional_Requested_Information.md)):

| Package | Version | Purpose |
|---|---|---|
| CellRank | 2.0.7 | Fate probability computation and trajectory inference |
| scikit-learn | 1.6.1 | Internal ML utilities (separate environment from cNMF/STARcat) |

---

## Citation

If you use this code, please cite:

> [Citation will be added upon publication]

---

## Contact

For questions about the analysis code, please open an issue on this GitHub repository. For questions about the biological findings or raw data, please contact the corresponding author.
