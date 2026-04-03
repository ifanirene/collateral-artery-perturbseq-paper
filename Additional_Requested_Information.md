# Additional Requested Information

This is a checklist of information still needed before the repository is ready for final code deposition. Items are grouped by category. Please fill in or provide each item.

---

## 1. Missing Scripts

The following analyses are described in the Methods and appear in manuscript figures but have no corresponding script in `analysis/`. Please provide a script (or confirm it will not be released) for each.

- [x] **scATAC-seq QC and coaccessibility** — Added as Script 11 (`11_heart_ec_atac_signac_cicero_preprocessing.R`). *(Methods: "scATAC-seq-Specific Processing, Quality Control, and Integration"; Figure S1 ATAC panels)*
- [ ] **CellOracle GRN construction and TF perturbation scoring** — custom guinea pig TSS annotation, motif scanning, GRN edge fitting, pseudotime, ps_sum scoring, Bonferroni correction. *(Methods: "CellOracle GRN Construction"; Figures 3a, S3f–h)*
- [x] **Cross-species pathway enrichment** — Added as Script 12 (`12_heart_brain_ec_hallmark_pathway_comparison.py`). *(Methods: "Differential Expression, Pathway Enrichment"; Figure 2h)*
- [ ] **scFEA metabolic flux estimation** — applied to the heart atlas across all cardiac cell types. *(Methods: "Metabolic-Flux Inference"; Figure 2i)*
- [x] **Gene-level SCEPTRE artery-score screen** — Split into Script 13 (UCell scoring: `13_endothelial_state_ucell_scoring.R`) and Script 14 (SCEPTRE screen: `14_endothelial_state_sceptre_artery_score.R`). *(Methods: "Primary Artery-State Score Derivation"; Figures 3g–h, 5a–b, S4k)*
- [x] **Transcriptome-wide gene-level SCEPTRE testing** — Added as Script 15 (`15_perturbseq_gene_level_sceptre_transcriptome.R`).
- [ ] **CellRank trajectory and artery fate probability** — fate probability computation, pseudotime-binned program usage profiles for Fzd4/Fzd6 KD vs. NTC. *(Figure 5c, Figure S6a–g)*
- [ ] **Within-screen attenuation/permutation analysis** — testing whether Program 20 effects are attenuated after accounting for Programs 6 or 18; 2,000 restricted permutations. *(Figure S6g)*
- [ ] **Neonatal hypoxia scRNA-seq preprocessing** — Seurat loading of normoxia and two hypoxia batches, SCT integration, cluster annotation, export as h5ad for Script 09. *(Methods: "scRNA-seq Library Preparation and Endothelial Annotation" in hypoxia section; Figures 5o, S7)*
- [ ] **Endothelial subtype gene program enrichment** — per-cell usage normalization, one-sided Mann-Whitney U per (subtype, program) pair, BH correction per subtype. *(Methods: "Endothelial-Subtype Enrichment of Gene Programs"; Figures 5h, S5b)*

---

## 2. Missing Software Versions

The following packages are referenced in the Methods but had **no version number stated**. Items resolved by the user or recovered from script headers are checked off; the remaining item is still needed.

- [x] **Cell Ranger ARC** — v2.0.2 (confirmed by user)
- [x] **Monocle3** — v1.3.1 (recovered from environment snapshot in Script 11 header)
- [x] **CellRank** — v2.0.7 (confirmed by user; `traj-env` conda environment)
- [x] **scikit-learn (cNMF/STARcat env)** — v1.0.2 (confirmed by user; `cnmf_env`)
- [x] **scikit-learn (CellRank/traj env)** — v1.6.1 (confirmed by user; separate `traj-env`)
- [x] **Seurat v5 exact version** — v5.3.0 (confirmed by user; `seurat5_r` conda environment)
- [x] **DESeq2** — v1.42.0 (confirmed by user; `seurat5_r`)
- [x] **jsonlite** — v2.0.0 (confirmed by user; `seurat5_r`)
- [x] **data.table** — v1.17.8 (confirmed by user; `seurat5_r`)
- [ ] **scFEA** — version still unknown; no conda environment identified

---

## 3. Reference Genomes and Databases

- [x] **Guinea pig Cell Ranger reference**: Build script and configuration are now in `genome_reference/` (`CavPor4.config`, `reproduce_guinea_pig_genome_reference_build.R`, `Genome_reference.md`). Built with `cellranger-arc mkref` v2.0.2; 86 primary scaffolds, protein-coding + lncRNA genes only, custom JASPAR2024 motif file.
- [ ] **Ensembl BioMart access date**: When was the guinea pig–mouse ortholog table downloaded? (Needed for reproducibility.)
- [ ] **MSigDB Hallmark gene sets version**: Which version / access date was used for the decoupler pathway analysis? v2024.1

---

## 4. Supplementary Files Needed in the Repository or GEO Deposit

- [ ] **feature_reference.csv** (Table S4): The gRNA library sequence file needed to run the Cell Ranger guide-capture step. Should be deposited in GEO or added to the repository.
- [ ] **Guinea pig–mouse ortholog table**: The BioMart-derived one-to-one ortholog file used in Scripts 01 and 02. Should be provided as a data file so users can reproduce the integration without re-querying BioMart.
- [ ] **Xenium panel gene list**: The 299-gene panel design (gene names, target transcript regions) used for the spatial validation experiment.

---

## 5. Other Open Items

- [ ] **Random seeds**: Seeds are documented for cNMF (14, 123). Are seeds set and recorded for scDblFinder, UMAP layout, k-means spatial replicates (Xenium), or CellOracle pseudotime? If so, please provide.
- [ ] **CellChat database version**: Confirm the exact CellChat database version bundled with CellChat v2.1.2 used for the analysis.
- [ ] **scFEA metabolic module set**: Which curated metabolic module set was used for scFEA (the package ships with several options)?
