# Collateral Artery Perturb-seq: Code Repository

This repository contains the computational analysis code accompanying the manuscript:

**"A cross-species Perturb-seq platform identifies artery repressors that govern collateral artery formation"**

All computational analyses were performed using the custom scripts provided here. Raw sequencing data are deposited in GEO (accession pending).

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
│   └── 04_heart_xenium_ec_pseudobulk_validation.R
│
├── 02_perturbseq_screen/            # In vivo Perturb-seq preprocessing
│   └── 05_perturbseq_endothelial_preprocessing.R
│
└── 03_gene_program_analysis/        # cNMF program discovery, SCEPTRE testing, and cross-dataset projection
    ├── 06_perturbseq_cnmf_k100_factorization.py
    ├── 07_perturbseq_cnmf_discovery_notebook.ipynb
    ├── 08_perturbseq_cnmf_program_sceptre.R
    ├── 09_project_programs_onto_hypoxia_normoxia.py
    └── 10_project_programs_onto_crossspecies_brain_ec.py
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

---

### Module 2 — In Vivo Perturb-seq Screen (`02_perturbseq_screen/`)

**Script 05 — Perturb-seq endothelial preprocessing** (`05_perturbseq_endothelial_preprocessing.R`)

This script reconstructs the manuscript-facing preprocessing path for the neonatal mouse brain CRISPRi Perturb-seq screen. Starting from Cell Ranger outputs for ten sequencing batches, it combines resequenced RNA counts with original CRISPR guide capture matrices by intersecting shared barcodes. After RNA-level quality control (2,000–7,500 genes, 3,000–30,000 UMIs, less than 10% mitochondrial reads), the merged object is clustered and non-endothelial populations are removed through three rounds of manual cluster pruning to retain endothelial, perivascular, and screen-relevant populations. Doublet detection is then run in batch-aware mode using scDblFinder. Guide multiplicity is estimated from UMI thresholding, and cells with more than 15 assigned guides or zero guide UMIs are excluded. The final filtered RNA and gRNA count matrices, together with per-cell covariates (batch, cell type, mitochondrial fraction), are exported as an RData bundle consumed by all downstream SCEPTRE analyses (Script 08) and an h5ad file consumed by cNMF (Scripts 06–07). This script supports Figure 3b–e and Figure S4.

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

**Script 08 — SCEPTRE program-level perturbation testing** (`08_perturbseq_cnmf_program_sceptre.R`)

This script tests whether each CRISPRi perturbation significantly alters the usage of each of the 100 cNMF gene programs discovered in Scripts 06–07. Per-cell program usage values are converted into count-like responses by row-normalizing each cell's usage vector to sum to one and multiplying by the cell's total UMI count before rounding, making the SCEPTRE conditional-resampling model applicable. The script aligns cells across the program usage matrix, the gene-level SCEPTRE input bundle from Script 05, and the guide annotation table, then runs SCEPTRE in high-MOI mode with batch and library-size covariates. Both all-endothelial-cell and per-cell-type analyses are supported through a configuration flag. The output discovery tables identify which perturbations increase or decrease each program and are used to characterize how artery repressors act through WNT, HIF-glycolysis, VEGF-tip-cell, and neuronal-guidance programs (Figure 5e–h, Figure S5b–j).

**Script 09 — STARcat projection onto hypoxia/normoxia brain endothelial cells** (`09_project_programs_onto_hypoxia_normoxia.py`)

This script uses STARcat to project the 100 cNMF gene programs from the Perturb-seq screen onto an independently collected dataset of neonatal mouse brain endothelial cells exposed to either normoxia (21% O₂) or mild hypoxia (11% O₂) for seven days. The projection holds the cNMF consensus spectra fixed and estimates per-cell usage scores by non-negative least squares (NNLS), column-standardizing the query counts before fitting. Row-normalized usage scores (summing to one per cell) are saved alongside batch and cell-type annotations for downstream differential-usage analysis. This projection tests whether the HIF-glycolysis, canonical WNT, and tip-cell programs identified as targets of artery repressors in the Perturb-seq screen are also induced by environmental oxygen reduction, and whether artery-associated programs are correspondingly diminished (Figure 5m–q).

**Script 10 — STARcat projection onto cross-species brain endothelial cells** (`10_project_programs_onto_crossspecies_brain_ec.py`)

This script applies the same STARcat projection framework to the integrated cross-species brain endothelial cell dataset from Module 1, enabling direct comparison of program activity between mouse and guinea pig brain endothelial cells. Because the query includes guinea pig cells whose gene symbols were mapped to mouse orthologs during upstream integration, only the 6,704 reference genes present in the shared ortholog space are used for fitting (compared with approximately 9,755 for a pure mouse query). The resulting per-cell usage matrix is saved alongside species and cell-type metadata for stratified program comparisons. This analysis tests whether the artery-associated programs elevated in guinea pig ECs and the hypoxia, WNT, and tip-cell programs reduced in guinea pig ECs mirror the program shifts observed in the Perturb-seq screen, closing the loop between the cross-species atlas and the in vivo perturbation results (Figure 5q, bottom row).

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

### R Environment — Seurat 4 (heart atlas, brain atlas, CellChat, Perturb-seq preprocessing)

Scripts 01–03 and 05 were developed and validated in R 4.2.0 on a Linux (CentOS 7) system using the packages listed below. An environment snapshot is recorded in the header of each script.

| Package | Version | Purpose |
|---|---|---|
| Seurat | 4.4.0 | Single-cell preprocessing, integration, clustering, and visualization |
| SeuratObject | 5.0.1 | Seurat object infrastructure |
| scDblFinder | 1.17.2 | Cluster-aware doublet detection |
| BiocParallel | 1.32.5 | Parallel execution for doublet scoring |
| gprofiler2 | 0.2.3 | Ortholog mapping for cell-cycle gene sets |
| CellChat | 2.1.2 | Ligand-receptor communication inference (Script 03 only) |
| sceptre | 0.1.0 / 0.10.0 | CRISPRi perturbation testing (Script 08) |
| ggplot2 | 3.5.0 | Visualization |
| patchwork | 1.2.0 | Plot composition |
| dplyr | 1.1.4 | Table manipulation |

### R Environment — Seurat 5 (Xenium spatial workflow)

Script 04 requires Seurat v5 for sketch-based integration of Xenium region outputs:

| Package | Version | Purpose |
|---|---|---|
| Seurat | ≥ 5.0 | Xenium loading, sketch integration, and cluster projection |
| DESeq2 | | Pseudobulk differential expression |
| jsonlite | | Parsing Xenium experiment metadata |
| data.table | | Reading `cells.csv.gz` per-cell metadata |

### Python Environment — cNMF and STARcat (Scripts 06, 07, 09, 10)

| Package | Version | Purpose |
|---|---|---|
| Python | ≥ 3.10 | Runtime |
| scanpy | | AnnData handling and preprocessing |
| anndata | 0.10.7 | AnnData I/O |
| cnmf | 1.7.0 | Consensus NMF factorization |
| starcatpy | 1.0.9 | NNLS-based projection onto cNMF reference spectra |
| numpy | | Numerical operations |
| pandas | | DataFrame handling |
| scikit-learn | | NNLS solver (used internally by STARcat) |

---

## Citation

If you use this code, please cite:

> [Citation will be added upon publication]

---

## Contact

For questions about the analysis code, please open an issue on this GitHub repository. For questions about the biological findings or raw data, please contact the corresponding author.
