# Guinea Pig Genome Reference Resources for Cross-Species Multiome Analysis

This document describes the guinea pig (*Cavia porcellus*) genome reference resources
built for this study. Both resources are derived from the NCBI mCavPor4.1 assembly
(GCF_034190915.1): a custom Cell Ranger ARC reference for primary read alignment and
feature counting, and a custom BSgenome R data package for downstream Signac chromatin
assay construction. Neither is available from 10x Genomics or Bioconductor. Large
binary files (genome FASTA, indexed reference, R package tarballs) can be provided
upon request; the build scripts and configuration files in this directory are sufficient
to reconstruct all resources from the NCBI source files.

> **Version note:** All references to genome version throughout this study use the NCBI
> assembly name **mCavPor4.1**. Internal build iteration labels (CavPor4.2, CavPor4.3)
> that appear in R package names and output directory names are implementation details;
> the underlying genome sequence and annotation are unchanged from the NCBI
> GCF_034190915.1 release.

---

## 1. Source Genome Assembly

| Field | Value |
|---|---|
| Species | *Cavia porcellus* (guinea pig) |
| Assembly name | mCavPor4.1 |
| NCBI accession | GCF_034190915.1 |
| Assembly level | Scaffold |
| Annotation release | RS_2024_02 (February 2023) |
| Source URL | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034190915.1/ |

**Files downloaded from NCBI:**
- `GCF_034190915.1_mCavPor4.1_genomic.fna` — whole-genome FASTA (3.0 GB)
- `GCF_034190915.1_mCavPor4.1_genomic.gtf` — gene annotation in GTF format (599 MB)
- `GCF_034190915.1_mCavPor4.1_genomic.gff.gz` — gene annotation in GFF3 format (21 MB)

The assembly uses NCBI scaffold accession identifiers (e.g., NW_026947484.1) rather
than chromosome-style names. There is no mitochondrial genome contig in this assembly
release. The 86 largest scaffolds (NW_026947484.1 through NW_026948059.1), which
correspond to the chromosome-scale assignments in the mCavPor4.1 annotation, were
designated as primary contigs in both the Cell Ranger ARC reference and the BSgenome
package.

---

## 2. Custom Cell Ranger ARC Multiome Reference for mCavPor4.1

The Cell Ranger ARC reference was built with `cellranger-arc mkref` v2.0.2 following
the 10x Genomics non-standard genome guide.

### 2a. GTF preparation

The original NCBI GTF contains 37,419 annotated genes across 13 biotype categories.
Before running `mkref`, the GTF was filtered to retain only gene biotypes compatible
with single-cell RNA-seq quantification and that reduce spurious multi-mapping:

| Biotype | Genes in original GTF | Retained |
|---|---|---|
| protein_coding | 21,479 | Yes |
| lncRNA | 7,201 | Yes |
| pseudogene | 5,226 | No |
| snoRNA | 1,541 | No |
| snRNA | 1,000 | No |
| tRNA | 463 | No |
| V_segment | 313 | No |
| rRNA | 110 | No |
| ncRNA | 56 | No |
| C_region | 15 | No |
| misc_RNA | 10 | No |
| transcribed_pseudogene | 3 | No |
| other | 2 | No |

**28,680 genes retained** (protein_coding + lncRNA) in `GCF_034190915.1_mCavPor4.1_filtered.gtf`.

Pseudogenes and small RNA genes were excluded because they inflate apparent complexity
and increase multi-mapping; immunoglobulin gene segments (V_segment, C_region) were
excluded because they are absent from the heart endothelial cell types being profiled.

### 2b. Custom JASPAR2024 motif file

The standard Cell Ranger ARC motif library was replaced with a custom motif file
derived from JASPAR2024 Core vertebrates (non-redundant). The file
`JASPAR2024_modified.txt` contains PFM entries reformatted for compatibility with
`cellranger-arc mkref`. This file can be provided upon request.

### 2c. Reference build command

The build was driven by the config file `CavPor4.config` (copied to this directory):

```json
{
    organism: "Cavia_porcellus"
    genome: ["mCavPor4"]
    input_fasta: ["GCF_034190915.1_mCavPor4.1_genomic.fna"]
    input_gtf: ["GCF_034190915.1_mCavPor4.1_filtered.gtf"]
    input_motifs: "JASPAR2024_modified.txt"
}
```

Build command (run from the directory containing the config and input files):

```bash
cellranger-arc mkref \
    --config=CavPor4.config \
    --nthreads=32 \
    --memgb=256
```

Build parameters: 32 threads, 256 GB RAM. STAR v2.7.2a was used internally for
genome index generation.

### 2d. Reference output

The completed reference directory contains:

```
mCavPor4.1_cellranger_arc_ref/
  fasta/
    genome.fa           (3.0 GB — genomic FASTA restricted to primary contigs)
    genome.fa.{amb,ann,bwt,fai,pac,sa}  (BWA index for ATAC alignment)
  genes/
    genes.gtf.gz        (compressed GTF restricted to primary contigs)
  star/                 (STAR index for RNA alignment)
  regions/
    transcripts.bed     (transcript coordinate BED)
    tss.bed             (transcription start site BED)
    motifs.pfm          (custom JASPAR2024 motif file)
  reference.json        (build metadata; records mkref version and primary contigs)
```

The final reference indexes 1,552 protein-coding and 563 lncRNA genes on the 86
primary scaffolds. The reduction from 28,680 to 2,115 genes reflects that only genes
located on the 86 primary contigs are indexed; the remaining ~26,500 genes reside on
smaller unplaced scaffolds not included as primary contigs.

### 2e. Primary contigs

All 86 primary contig accessions are listed in `reference.json`. The full list:
NW_026947484.1, NW_026947485.1, NW_026947486.1, NW_026947487.1, NW_026947488.1,
NW_026947489.1, NW_026947490.1, NW_026947491.1, NW_026947492.1, NW_026947493.1,
NW_026947494.1, NW_026947495.1, NW_026947496.1, NW_026947497.1, NW_026947498.1,
NW_026947499.1, NW_026947500.1, NW_026947501.1, NW_026947502.1, NW_026947503.1,
NW_026947504.1, NW_026947505.1, NW_026947506.1, NW_026947507.1, NW_026947508.1,
NW_026947509.1, NW_026947510.1, NW_026947511.1, NW_026947512.1, NW_026947513.1,
NW_026947514.1, NW_026947515.1, NW_026947516.1, NW_026947517.1, NW_026947518.1,
NW_026947519.1, NW_026947520.1, NW_026947521.1, NW_026947522.1, NW_026947523.1,
NW_026947524.1, NW_026947525.1, NW_026947526.1, NW_026947527.1, NW_026947528.1,
NW_026947529.1, NW_026947530.1, NW_026947532.1, NW_026947533.1, NW_026947536.1,
NW_026947540.1, NW_026947541.1, NW_026947543.1, NW_026947544.1, NW_026947546.1,
NW_026947547.1, NW_026947549.1, NW_026947564.1, NW_026947565.1, NW_026947567.1,
NW_026947571.1, NW_026947574.1, NW_026947575.1, NW_026947579.1, NW_026947603.1,
NW_026947616.1, NW_026947617.1, NW_026947619.1, NW_026947638.1, NW_026947660.1,
NW_026947688.1, NW_026947713.1, NW_026947716.1, NW_026947760.1, NW_026947767.1,
NW_026947772.1, NW_026947778.1, NW_026947784.1, NW_026947794.1, NW_026947825.1,
NW_026947855.1, NW_026947908.1, NW_026947911.1, NW_026947943.1, NW_026947957.1,
NW_026947974.1, NW_026948011.1, NW_026948059.1

Note: the assembly does not include a mitochondrial contig; `non_nuclear_contigs` is
set to `null` in the config and reference.json.

---

## 3. Custom BSgenome Data Package for mCavPor4.1

The BSgenome package provides genome sequences for use in R/Bioconductor workflows
(Signac chromatin assay construction, motif scanning). It is **not available from
Bioconductor** and must be installed from the provided tarball. The R package name
is `BSgenome.CPorcellus.NCBI.CavPor4.3` (an internal build label); the underlying
genome is mCavPor4.1 (GCF_034190915.1).

### 3a. Genome FASTA to 2-bit conversion

The UCSC `faToTwoBit` utility converts the NCBI genomic FASTA to 2-bit format,
the storage format required by the BSgenome package infrastructure:

```bash
faToTwoBit GCF_034190915.1_mCavPor4.1_genomic.fna src_dir/cavpor4.2bit
```

The `faToTwoBit` binary (UCSC utilities, April 2021) is included in the build
directory alongside the seed file.

### 3b. BSgenome seed file

The seed file `BSgenome.CPorcellus.NCBI.cavPor4-seed` (in this directory) specifies
the package metadata and links to the 2-bit genome file:

```
Package: BSgenome.CPorcellus.NCBI.CavPor4.3
Title: Full genome sequences for Cavia porcellus (mCavPor4.1)
organism: Cavia porcellus
common_name: Guinea Pig
provider: NCBI
genome: GCF_034190915.1
release_date: Feb. 2023
source_url: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034190915.1/
BSgenomeObjname: CPorcellus
seqfile_name: cavpor4.2bit
seqs_srcdir: ./src_dir        (directory containing cavpor4.2bit)
circ_seqs: character(0)       (no circular sequences; no mitochondrial contig in assembly)
seqnames: [86 NCBI scaffold accessions from mCavPor4.1 GTF seqlevels]
```

### 3c. Package build command

```r
library(BSgenome)
# Run from the directory containing the seed file and src_dir/cavpor4.2bit
forgeBSgenomeDataPkg("BSgenome.CPorcellus.NCBI.cavPor4-seed")
```

This produces the installable package tarball (636 MB).
The complete build script is `reproduce_guinea_pig_genome_reference_build.R` in this
directory; it also covers TxDb and OrgDb construction.

### 3d. Installation

```r
# Install from the provided tarball
devtools::install_local("BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz")
# or from the shell:
# R CMD INSTALL BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz
```

### 3e. Usage in R

```r
library(BSgenome.CPorcellus.NCBI.CavPor4.3)
# The genome object is accessible as:
BSgenome.CPorcellus.NCBI.CavPor4.3
# or via the short alias registered at load time:
CPorcellus
```

### 3f. Sequence names and chromosome notation

The package exposes 86 sequences with NCBI accession identifiers (e.g.,
`NW_026947484.1`). These identifiers contain underscores and dots, which conflict
with downstream tools that use underscores or hyphens as coordinate-field delimiters.
In Signac and Cicero, all underscores in seqlevels are replaced with hyphens before
peak-matrix construction (see `reproduce_heart_ec_atac_signac_cicero_preprocessing.R`).

---

## 4. Additional Annotation Resources

The build script `reproduce_guinea_pig_genome_reference_build.R` also documents:

- **TxDb object** — built from the NCBI GFF3 annotation using
  `GenomicFeatures::makeTxDbFromGFF()` with `dbxrefTag = "GeneID"` to resolve NCBI
  gene identifiers; saved as `cavpor4.TxDb.sqlite`. Used for transcript-level
  annotation lookups in CellOracle TSS window construction.

- **OrgDb package (org.Cporcellus.eg.db v4.1)** — built using
  `AnnotationForge::makeOrgPackageFromNCBI()` for NCBI taxonomy ID 10141 (*Cavia
  porcellus*); provides gene symbol, Entrez ID, GO, and RefSeq mappings. The package
  tarball can be provided upon request.

---

## 5. Files in This Directory

| File | Description |
|---|---|
| `Genome_reference.md` | This document |
| `BSgenome.CPorcellus.NCBI.cavPor4-seed` | BSgenome seed file (package metadata and seqnames) |
| `reproduce_guinea_pig_genome_reference_build.R` | Annotated build script for BSgenome, TxDb, and OrgDb |
| `CavPor4.config` | `cellranger-arc mkref` configuration for the mCavPor4.1 reference |

## 6. Files Available Upon Request

The following large binary files were produced by the build and are available
upon request:

| File | Size | Contents |
|---|---|---|
| `BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz` | 636 MB | Installable BSgenome R package for mCavPor4.1 |
| `mCavPor4.1_cellranger_arc_ref.tar.gz` | ~31 GB | Complete Cell Ranger ARC reference directory |
| `GCF_034190915.1_mCavPor4.1_filtered.gtf` | 597 MB | Biotype-filtered GTF used for `mkref` |
| `JASPAR2024_modified.txt` | 278 KB | Custom JASPAR2024 Core vertebrate motif file |
| `org.Cporcellus.eg.db_4.1.tar.gz` | 42 MB | Guinea pig OrgDb annotation package |
