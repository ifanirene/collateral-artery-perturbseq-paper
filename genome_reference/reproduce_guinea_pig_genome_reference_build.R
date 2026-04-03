#!/usr/bin/env Rscript

# @theme: Custom genome annotation resources for Cavia porcellus (mCavPor4.1 assembly)
# @task_id: genome-ref-01
# @description:
#   Builds three R/Bioconductor annotation resources for guinea pig (*Cavia porcellus*)
#   from the NCBI mCavPor4.1 genome assembly (GCF_034190915.1, release RS_2024_02).
#   None of these resources is available from Bioconductor or CRAN.
#
#   Resources produced:
#     1. BSgenome data package (BSgenome.CPorcellus.NCBI.CavPor4.3 v4.0) — provides
#        the genome sequence for Signac chromatin assay construction and motif scanning.
#        The NCBI genomic FASTA is first converted to 2-bit format with faToTwoBit
#        (UCSC utility), then packaged with BSgenome::forgeBSgenomeDataPkg() using the
#        accompanying seed file (BSgenome.CPorcellus.NCBI.cavPor4-seed).
#
#     2. TxDb SQLite database (cavpor4.TxDb.sqlite) — transcript-level annotation
#        built from the NCBI GFF3 file using GenomicFeatures::makeTxDbFromGFF() with
#        dbxrefTag = "GeneID" to resolve NCBI gene IDs. Used for CellOracle TSS
#        window construction and general transcript annotation.
#
#     3. OrgDb package (org.Cporcellus.eg.db v4.1) — gene-level annotation (Entrez ID,
#        gene symbol, GO terms, RefSeq IDs) for Cavia porcellus (NCBI taxonomy ID 10141),
#        built from NCBI gene database flat files using AnnotationForge::makeOrgPackageFromNCBI().
#
# @inputs:
#   PRE-STEP (shell, not R): Convert genome FASTA to 2-bit format with faToTwoBit.
#   - `GCF_034190915.1_mCavPor4.1_genomic.fna` : NCBI whole-genome FASTA (3.0 GB);
#       downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034190915.1/
#   - `GCF_034190915.1_mCavPor4.1_genomic.gff.gz` : NCBI GFF3 annotation (21 MB);
#       downloaded from the same NCBI assembly page.
#   - `BSgenome.CPorcellus.NCBI.cavPor4-seed` : BSgenome seed file specifying the
#       package name, genome metadata, and the 86 primary scaffold accession names.
#   - `src_dir/cavpor4.2bit` : 2-bit genome file produced by the faToTwoBit pre-step.
#   - NCBI gene database flat files for Cavia porcellus (taxonomy ID 10141):
#       gene_info.gz, gene2accession.gz, gene2go.gz, gene2pubmed.gz, gene2refseq.gz;
#       downloaded from https://ftp.ncbi.nlm.nih.gov/gene/DATA/ and placed in the
#       NCBI files directory specified by ncbi_files_dir below.
#
# @outputs:
#   - `BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz` : installable BSgenome R package
#       tarball (636 MB); install with devtools::install_local() or R CMD INSTALL.
#   - `cavpor4.TxDb.sqlite` : TxDb SQLite database for guinea pig transcript annotation.
#   - `org.Cporcellus.eg.db_4.1.tar.gz` : OrgDb package tarball (42 MB) providing
#       gene-level ID mappings; install with devtools::install_local() or R CMD INSTALL.
#
# @key_params:
#   - genome_accession = "GCF_034190915.1"  : NCBI genome accession for mCavPor4.1.
#   - taxonomy_id = "10141"                  : NCBI taxonomy ID for Cavia porcellus.
#   - orgdb_version = "4.1"                  : version label for the OrgDb package.
#   - txdb_biotypes = c("mRNA","lnc_RNA","rRNA","tRNA","ncRNA") : GFF3 feature types
#       retained in the exported transcript GRanges (NCBI GFF3 biotype names differ
#       from GTF names; "mRNA" corresponds to protein_coding, "lnc_RNA" to lncRNA).
#
# @dependencies:
#   - R 4.2.0
#   - BSgenome (>= 1.66.3)      : forgeBSgenomeDataPkg(), getChromInfoFromNCBI().
#   - GenomicFeatures            : makeTxDbFromGFF(), transcripts(), saveDb(), loadDb().
#   - GenomicRanges              : GRanges manipulation.
#   - rtracklayer                : GFF3 import (used for verification; not required for
#                                  makeTxDbFromGFF itself).
#   - AnnotationForge            : makeOrgPackageFromNCBI().
#   - devtools                   : install_local() for the built package tarballs.
#   - httr                       : SSL configuration for NCBI download calls.
#   External tools (shell):
#   - faToTwoBit (UCSC utilities, April 2021) : genome FASTA to 2-bit conversion.
#
# @environment:
#   Run in the seurat4 conda environment (R 4.2.0) or any R 4.2.x environment with
#   the above Bioconductor packages installed. The BSgenome and AnnotationForge builds
#   require internet access to NCBI during the OrgDb step.
#
# @examples:
#   Rscript reproduce_guinea_pig_genome_reference_build.R
#
#   Run interactively block by block. Each STEP either produces a file on disk or
#   installs a package into the R library. The pre-step (faToTwoBit) must be run from
#   the shell before starting this script.
#
#   Shell pre-step (run once before this script):
#     ./faToTwoBit GCF_034190915.1_mCavPor4.1_genomic.fna src_dir/cavpor4.2bit
#
# Recorded environment (build completed May 2024):
#   R: 4.2.0 | platform: x86_64-conda-linux-gnu | OS: CentOS Linux 7
#   BSgenome: >= 1.66.3 | GenomicFeatures: >= 1.50 | AnnotationForge: >= 1.40
#   faToTwoBit: UCSC utilities April 2021

options(stringsAsFactors = FALSE)
options(timeout = 100000)  # Needed for large NCBI downloads during OrgDb build

suppressPackageStartupMessages({
  library(BSgenome)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(AnnotationForge)
  library(devtools)
  library(httr)
})

sessionInfo()

###############################################################################
# File paths. Replace placeholder paths with actual locations before running.
###############################################################################

# Directory containing the seed file, src_dir/cavpor4.2bit, and NCBI input files.
build_dir <- "path/to/genome_reference_build_directory"

# Path to the NCBI GFF3 annotation file.
gff3_path <- file.path(build_dir, "GCF_034190915.1_mCavPor4.1_genomic.gff.gz")

# Directory containing NCBI gene database flat files (gene_info.gz, gene2go.gz, etc.).
ncbi_files_dir <- "path/to/ncbi_gene_flatfiles_directory"

# Output directory for the OrgDb package source before tarballing.
output_dir <- build_dir

# Analysis metadata.
genome_accession <- "GCF_034190915.1"
taxonomy_id      <- "10141"            # NCBI taxonomy ID for Cavia porcellus
orgdb_version    <- "4.1"
author_name      <- "Your Name <your@email.edu>"

###############################################################################
# PRE-STEP (shell command — run before this script, not in R).
#
# The NCBI genomic FASTA must be converted to 2-bit format before the BSgenome
# package can be forged. The faToTwoBit binary (UCSC utilities) is provided in
# the build directory. Run this once from the shell:
#
#   cd path/to/genome_reference_build_directory
#   mkdir -p src_dir
#   ./faToTwoBit GCF_034190915.1_mCavPor4.1_genomic.fna src_dir/cavpor4.2bit
#
# The resulting cavpor4.2bit file (~500 MB) is the genome sequence store referenced
# by the BSgenome seed file (seqfile_name: cavpor4.2bit, seqs_srcdir: ./src_dir).
###############################################################################

###############################################################################
# STEP 1. Retrieve chromosome information from NCBI to verify the mCavPor4.1
#         scaffold accessions before forging the BSgenome package.
#
# getChromInfoFromNCBI() queries the NCBI assembly record and returns a data frame
# with scaffold accessions (RefSeqAccn), sequence lengths, and circularity flags.
# The mCavPor4.1 assembly has 86 primary chromosomal scaffolds (NW_026947484.1
# through NW_026948059.1) plus many smaller unplaced scaffolds. The BSgenome seed
# file restricts the package to the 86 primary scaffolds, which are confirmed here
# against the NCBI record.
#
# Output of this step: chrominfo data frame (used in STEP 3 for TxDb construction).
###############################################################################

cat("STEP 1: Retrieving chromosome information from NCBI for", genome_accession, "\n")
chrominfo <- getChromInfoFromNCBI(genome_accession)
cat(sprintf("  Retrieved %d scaffold entries from NCBI.\n", nrow(chrominfo)))
cat(sprintf("  First scaffold: %s (%d bp)\n",
            chrominfo$RefSeqAccn[1], chrominfo$SequenceLength[1]))

###############################################################################
# STEP 2. Forge the BSgenome data package from the seed file.
#
# forgeBSgenomeDataPkg() reads the seed file, locates the 2-bit genome file
# referenced by seqs_srcdir and seqfile_name, packages the 86 primary scaffold
# sequences, and writes an installable R package tarball named
# BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz.
#
# The seed file (BSgenome.CPorcellus.NCBI.cavPor4-seed) must be present in the
# working directory, and src_dir/cavpor4.2bit must have been created in the shell
# pre-step above.
#
# The R package name "BSgenome.CPorcellus.NCBI.CavPor4.3" is an internal build
# label; the underlying genome is mCavPor4.1 (GCF_034190915.1). The package
# exposes 86 genome sequences accessible by NCBI scaffold accession.
#
# Output: BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz (~636 MB) in the
# current working directory.
###############################################################################

cat("STEP 2: Forging BSgenome package from seed file...\n")
setwd(build_dir)
forgeBSgenomeDataPkg("BSgenome.CPorcellus.NCBI.cavPor4-seed")
cat("  Done. Package tarball written to:", build_dir, "\n")

###############################################################################
# STEP 3. Install the BSgenome package into the R library.
#
# devtools::install_local() installs from the tarball produced by STEP 2.
# Alternatively, run from the shell: R CMD INSTALL BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz
#
# After installation, load the package and verify that the genome object
# is accessible and returns the expected number of sequences.
###############################################################################

cat("STEP 3: Installing BSgenome package...\n")
devtools::install_local(
  file.path(build_dir, "BSgenome.CPorcellus.NCBI.CavPor4.3_4.0.tar.gz")
)

# Verify installation.
library(BSgenome.CPorcellus.NCBI.CavPor4.3)
gp_genome <- BSgenome.CPorcellus.NCBI.CavPor4.3
cat(sprintf("  BSgenome object loaded. Sequences: %d\n", length(seqnames(gp_genome))))
cat("  First five sequence names:", paste(head(seqnames(gp_genome), 5), collapse = ", "), "\n")

###############################################################################
# STEP 4. Build the TxDb transcript annotation database from the NCBI GFF3 file.
#
# makeTxDbFromGFF() parses the GFF3 annotation and constructs a SQLite-backed
# TxDb object that maps transcript identifiers to genomic coordinates. The
# dbxrefTag = "GeneID" argument instructs the function to use NCBI GeneID
# cross-references to name genes, which is necessary because NCBI GFF3 files
# use internal locus_tag identifiers rather than gene symbols as primary names.
#
# The chrominfo data frame from STEP 1 is used to supply scaffold lengths and
# circularity flags. The first entry in chrominfo (a non-canonical contig) is
# removed because it lacks a valid scaffold accession that matches the GFF3.
#
# After construction, the TxDb is queried to extract a GRanges of all transcripts
# with their gene IDs and biotype annotations. Transcripts are filtered to biotypes
# that correspond to protein-coding and functional RNA genes (NCBI GFF3 biotype
# names differ from Ensembl GTF names: "mRNA" = protein_coding, "lnc_RNA" = lncRNA).
# The genome field is set to the assembly accession for provenance tracking.
# The TxDb is then saved as a SQLite file for fast reloading.
#
# Output: cavpor4.TxDb.sqlite in the build directory.
###############################################################################

cat("STEP 4: Building TxDb from NCBI GFF3 annotation...\n")

# Prepare chrominfo: remove the first row, which corresponds to a non-standard
# contig that does not appear in the GFF3, and reformat column names to match
# the makeTxDbFromGFF() expectation.
chrominfo_txdb <- chrominfo[-1, ]
chrominfo_txdb$chrom       <- chrominfo_txdb$RefSeqAccn
chrominfo_txdb$length      <- chrominfo_txdb$SequenceLength
chrominfo_txdb$is_circular <- chrominfo_txdb$circular

txdb <- makeTxDbFromGFF(
  file       = gff3_path,
  format     = "gff3",
  chrominfo  = chrominfo_txdb,
  dataSource = "NCBI GFF3 for Cavia porcellus (GCF_034190915.1)",
  organism   = "Cavia porcellus",
  dbxrefTag  = "GeneID"   # Use NCBI GeneID cross-references to name genes
)
genome(txdb) <- genome_accession

txdb_path <- file.path(build_dir, "cavpor4.TxDb.sqlite")
saveDb(txdb, txdb_path)
cat(sprintf("  TxDb saved to: %s\n", txdb_path))

###############################################################################
# STEP 4b. Extract and filter the transcript GRanges for downstream use.
#
# This step exports a GRanges of filtered transcripts from the TxDb, suitable
# for use in CellOracle TSS annotation and general transcript-level lookups.
# Column names are standardized to the format expected by CellOracle's genome
# registration utilities (tx_id, gene_name, gene_id, gene_biotype).
#
# Biotypes retained (NCBI GFF3 nomenclature):
#   - mRNA       : protein-coding transcripts
#   - lnc_RNA    : long non-coding RNA
#   - rRNA       : ribosomal RNA
#   - tRNA       : transfer RNA
#   - ncRNA      : other non-coding RNA
#
# Output: guinea_pig_transcript_granges_filtered.rds in the build directory.
###############################################################################

txdb_biotypes <- c("mRNA", "lnc_RNA", "rRNA", "tRNA", "ncRNA")

txdb <- loadDb(txdb_path)
tx <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id", "tx_type"))

# Fill missing transcript names with gene IDs so every record has a non-NA name.
unnamed_idx <- which(is.na(mcols(tx)$tx_name))
mcols(tx)$tx_name[unnamed_idx] <- as.character(mcols(tx)$gene_id[unnamed_idx])

# Standardize column names to the CellOracle-expected schema.
mcols(tx)$tx_id      <- as.character(mcols(tx)$tx_name)
mcols(tx)$gene_name  <- as.character(mcols(tx)$gene_id)
mcols(tx)$gene_id    <- as.character(mcols(tx)$gene_id)
colnames(mcols(tx))  <- c("tx_id", "gene_name", "gene_id", "gene_biotype")

tx_filtered <- tx[mcols(tx)$gene_biotype %in% txdb_biotypes]
genome(tx_filtered) <- genome_accession

cat(sprintf("  Transcripts retained after biotype filter: %d of %d\n",
            length(tx_filtered), length(tx)))

saveRDS(tx_filtered,
        file.path(build_dir, "guinea_pig_transcript_granges_filtered.rds"))

###############################################################################
# STEP 5. Build the OrgDb annotation package from NCBI gene database flat files.
#
# makeOrgPackageFromNCBI() constructs a Bioconductor OrgDb package providing
# gene-level mappings (Entrez ID, gene symbol, GO terms, RefSeq accessions) for
# Cavia porcellus (NCBI taxonomy ID 10141). The function requires NCBI gene
# database flat files (gene_info.gz, gene2accession.gz, gene2go.gz, gene2pubmed.gz,
# gene2refseq.gz) to have been downloaded to the directory given by NCBIFilesDir.
# These files are available from https://ftp.ncbi.nlm.nih.gov/gene/DATA/
#
# SSL verification is disabled here because NCBI connections may reject certificate
# validation in some HPC network configurations; remove the httr::set_config() call
# if running in a standard network environment.
#
# Output: org.Cporcellus.eg.db source directory and _4.1.tar.gz tarball in output_dir.
###############################################################################

cat("STEP 5: Building OrgDb from NCBI gene database files...\n")

# Disable SSL verification if running behind an HPC proxy.
httr::set_config(httr::config(ssl_verifypeer = FALSE))

makeOrgPackageFromNCBI(
  version    = orgdb_version,
  author     = author_name,
  maintainer = author_name,
  outputDir  = output_dir,
  tax_id     = taxonomy_id,
  genus      = "Cavia",
  species    = "porcellus",
  NCBIFilesDir = ncbi_files_dir
)

# Install the OrgDb package from the tarball produced above.
orgdb_tarball <- file.path(output_dir,
                           paste0("org.Cporcellus.eg.db_", orgdb_version, ".tar.gz"))
devtools::install_local(orgdb_tarball)
cat(sprintf("  OrgDb package installed from: %s\n", orgdb_tarball))

# Quick verification: check that key columns are accessible.
library(org.Cporcellus.eg.db)
cat("  OrgDb key types available:", paste(keytypes(org.Cporcellus.eg.db), collapse = ", "), "\n")

cat("\nAll genome annotation resources for mCavPor4.1 built successfully.\n")
