#!/usr/bin/env Rscript

# @theme: scATAC-seq quality control, Signac chromatin assay construction, and Cicero
#         peak-coaccessibility scoring for CellOracle base GRN input
# @task_id: atac-preproc-01
# @description:
#   Parallel preprocessing of mouse (mm10) and guinea pig (mCavPor4.1) single-nucleus
#   ATAC-seq profiles for heart endothelial cells, culminating in Cicero coaccessibility
#   connections used as the chromatin-accessibility constraint for CellOracle gene-
#   regulatory network (GRN) construction. For each species, the workflow:
#     1. Extracts EC-assigned nuclei barcodes from the matched RNA-integrated Seurat object.
#     2. Builds a Signac ChromatinAssay from the aggregated multiome ATAC peak matrix.
#     3. Computes and applies ATAC quality-control filters (fragment depth, nucleosome
#        signal, and TSS enrichment) that select subnucleosomal open-chromatin profiles
#        from high-quality, transcriptionally active nuclei.
#     4. Converts the QC-filtered Signac object to a Monocle3 CellDataSet, applies
#        latent semantic indexing (LSI) and UMAP, and runs Cicero to estimate the
#        probability that pairs of regulatory elements co-occur in the same nucleus.
#     5. For mouse, performs binomial depth matching to align median in-peak fragment
#        depth to the guinea pig level before re-running Cicero, so that coaccessibility
#        score thresholds are applied equivalently across species.
#     6. Exports Cicero connections and peak lists in the underscore-delimited format
#        required by CellOracle.
#
# @inputs:
#   MOUSE
#   - [sample]_per_barcode_metrics.csv : per-nucleus ATAC QC metrics from Cell Ranger ARC
#       for three mouse embryonic heart samples (E13, E15, E17); barcodes are re-suffixed
#       to match the aggregated library index (E13=-1, E15=-2, E17=-3).
#   - mouse_heart_aggr_filtered_feature_bc_matrix.h5 : aggregated peak-barcode count
#       matrix from the three-sample Cell Ranger ARC aggregation run.
#   - mouse_heart_atac_fragments.tsv.gz : aggregated ATAC fragment file (tabix-indexed).
#   - mouse_heart_ec_integrated_seurat.rds : RNA-integrated Seurat object containing
#       EC-annotated nuclei from the three mouse embryonic heart samples; used only to
#       extract EC barcode lists.
#   - mm10_chromosome_length.txt : two-column tab-delimited file of mm10 chromosome
#       names and lengths required by Cicero.
#
#   GUINEA PIG
#   - [sample]_per_barcode_metrics.csv : per-nucleus ATAC QC metrics from Cell Ranger ARC
#       for three guinea pig embryonic heart samples (GD241, GD245, GD247); barcodes are
#       re-suffixed to match the aggregated library index (GD241=-2, GD245=-1, GD247=-3).
#   - guinea_pig_heart_aggr_filtered_feature_bc_matrix.h5 : aggregated peak-barcode count
#       matrix from the three-sample Cell Ranger ARC aggregation run.
#   - guinea_pig_heart_atac_fragments.tsv.gz : aggregated ATAC fragment file (tabix-indexed).
#   - guinea_pig_heart_ec_integrated_seurat.rds : RNA-integrated Seurat object containing
#       EC-annotated nuclei from the three guinea pig embryonic heart samples; used only to
#       extract EC barcode lists.
#   - GCF_034190915.1_mCavPor4.1_genomic.gtf : NCBI GTF annotation for the guinea pig
#       reference genome (mCavPor4.1 / GCF_034190915.1, release RS_2024_02); used to build
#       the Signac ChromatinAssay gene annotation.
#
# @outputs:
#   MOUSE
#   - mouse_heart_ec_atac_signac.rds : QC-filtered Signac/Seurat object for mouse EC nuclei
#       (peaks assay; mm10 genome; EnsDb.Mmusculus.v79 annotations).
#   - mouse_heart_ec_atac_monocle_cds.rds : Monocle3 CellDataSet after LSI preprocessing
#       and UMAP embedding; used as input to make_cicero_cds.
#   - mouse_heart_ec_cicero_connections.rds : Cicero coaccessibility connection table
#       (depth-matched run; k = 100).
#   - mouse_heart_ec_cicero_connections.csv : same connections, peak identifiers in
#       underscore-delimited format required by CellOracle.
#   - mouse_heart_ec_all_peaks.csv : all peak identifiers present in the QC-filtered
#       CellDataSet, underscore-delimited, for CellOracle peak-annotation step.
#
#   GUINEA PIG
#   - guinea_pig_heart_ec_atac_signac.rds : QC-filtered Signac/Seurat object for guinea
#       pig EC nuclei (peaks assay; mCavPor4.1 genome; NCBI GTF annotation).
#   - guinea_pig_heart_ec_atac_monocle_cds.rds : Monocle3 CellDataSet after LSI
#       preprocessing and UMAP embedding.
#   - guinea_pig_heart_ec_cicero_connections.rds : Cicero coaccessibility connection table.
#   - guinea_pig_heart_ec_cicero_connections.csv : same connections, underscore-delimited.
#   - guinea_pig_heart_ec_all_peaks.csv : all peak identifiers, underscore-delimited.
#
# @key_params:
#   - atac_min_fragments = 2000    : lower fragment-depth cutoff; removes low-complexity
#       profiles that are likely empty droplets or poorly captured nuclei.
#   - atac_max_fragments = 25000   : upper fragment-depth cutoff; removes likely doublets
#       or nuclei with anomalously high ATAC capture.
#   - max_nucleosome_signal = 2.5  : nucleosome-signal upper bound; retains sub-
#       nucleosomal open-chromatin profiles and excludes over-digested or closed chromatin.
#   - min_tss_enrichment = 2       : TSS enrichment lower bound; confirms that reads
#       are concentrated at actively regulated loci rather than distributed across
#       background regions.
#   - cicero_k_downsampled = 100   : number of nearest neighbors for make_cicero_cds in
#       the depth-matched mouse run; larger k compensates for reduced library complexity
#       after thinning.
#   - coaccessibility_threshold = 0.54 : minimum Cicero coaccessibility score retained
#       for CellOracle base GRN construction (applied downstream in CellOracle).
#
# @dependencies:
#   - R 4.2.0 (x86_64-conda-linux-gnu, CentOS 7)
#   - Signac 1.10.0          : Chromatin assay construction, QC metrics.
#   - Seurat 4.4.0           : Seurat object infrastructure.
#   - SeuratWrappers         : as.cell_data_set() for Monocle3 conversion.
#   - monocle3               : LSI, UMAP, and CDS infrastructure.
#   - cicero (monocle3 branch, cole-trapnell-lab/cicero-release) : coaccessibility scoring.
#   - Matrix / irlba         : sparse matrix and singular-value decomposition support.
#   - EnsDb.Mmusculus.v79    : mm10 gene annotation (mouse only).
#   - BSgenome.Mmusculus.UCSC.mm10 : mm10 genome reference (mouse only).
#   - BSgenome.CPorcellus.NCBI.CavPor4.2 : mCavPor4.1 genome reference (guinea pig only).
#   - rtracklayer            : GTF import for guinea pig annotation.
#   - GenomicRanges / GenomicFeatures : GRanges construction and manipulation.
#   - future                 : multisession parallelism for heavy Signac steps.
#   - ggplot2 / patchwork / cowplot : QC visualization.
#   - dplyr                  : data manipulation.
#
# @environment:
#   Run in the custom CellOracle R environment at:
#     /oak/stanford/groups/kredhors/irenefan/celloracle_env2.2/lib/r.4.2
#   This environment takes priority over the user default R library to ensure the
#   correct Cicero/Monocle3 versions are loaded.
#   Activate with: conda activate celloracle_env2.2 (or equivalent)
#
# @examples:
#   Rscript reproduce_heart_ec_atac_signac_cicero_preprocessing.R
#
#   Run interactively block by block. Each STEP saves a checkpoint object or CSV
#   that is loaded by the next STEP. Replace placeholder paths below with the
#   actual locations used in your release package.
#
# Recorded environment snapshot (celloracle conda environment, 2026-04-02):
#   R: 4.2.0 | platform: x86_64-conda-linux-gnu | OS: CentOS Linux 7
#   Signac: 1.10.0 | Seurat: 4.4.0 | monocle3: 1.3.1 | cicero: 1.3.9 (monocle3 branch)
#   Matrix: 1.5.1 | irlba: 2.3.5 | EnsDb.Mmusculus.v79: 79 | future: 1.33.0

options(stringsAsFactors = FALSE)

# Prepend the CellOracle-specific R library so Cicero/Monocle3 take priority
# over the user-level library, which may contain conflicting package versions.
.libPaths(c(
  "path/to/celloracle_env/lib/r.4.2",
  setdiff(.libPaths(), "path/to/user_default_r_library/4.2")
))

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(cicero)
  library(monocle3)
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(BSgenome.CPorcellus.NCBI.CavPor4.2)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(BSgenome)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(future)
})

sessionInfo()

###############################################################################
# SECTION 0. File paths and analysis parameters.
#
# Replace the placeholder paths below with the locations used in your release
# package. All output files are written to the same directory by default.
###############################################################################

# --- Input paths ---

# Mouse inputs
mouse_ec_rna_seurat_path  <- "path/to/mouse_heart_ec_integrated_seurat.rds"
mouse_e13_barcode_metrics <- "path/to/E13_per_barcode_metrics.csv"
mouse_e15_barcode_metrics <- "path/to/E15_per_barcode_metrics.csv"
mouse_e17_barcode_metrics <- "path/to/E17_per_barcode_metrics.csv"
mouse_atac_h5_path        <- "path/to/mouse_heart_aggr_filtered_feature_bc_matrix.h5"
mouse_fragment_file       <- "path/to/mouse_heart_atac_fragments.tsv.gz"
mouse_chrom_lengths_path  <- "path/to/mm10_chromosome_length.txt"

# Guinea pig inputs
gp_ec_rna_seurat_path  <- "path/to/guinea_pig_heart_ec_integrated_seurat.rds"
gp_gd241_barcode_metrics <- "path/to/GD241_per_barcode_metrics.csv"
gp_gd245_barcode_metrics <- "path/to/GD245_per_barcode_metrics.csv"
gp_gd247_barcode_metrics <- "path/to/GD247_per_barcode_metrics.csv"
gp_atac_h5_path          <- "path/to/guinea_pig_heart_aggr_filtered_feature_bc_matrix.h5"
gp_fragment_file         <- "path/to/guinea_pig_heart_atac_fragments.tsv.gz"
gp_gtf_path              <- "path/to/GCF_034190915.1_mCavPor4.1_genomic.gtf"

# Output directory
output_dir <- "path/to/output_directory"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Shared ATAC QC thresholds ---
# These are applied identically to both species so that cell-quality inclusion
# criteria are equivalent across the cross-species comparison.

atac_min_fragments   <- 2000   # Exclude low-complexity / empty-droplet profiles
atac_max_fragments   <- 25000  # Exclude likely doublets or over-captured nuclei
max_nucleosome_signal <- 2.5   # Retain subnucleosomal open-chromatin profiles
min_tss_enrichment    <- 2     # Confirm enrichment at actively regulated loci

# --- Cicero parameters ---
set.seed(2013)  # Applied before make_cicero_cds in all runs for reproducibility

cicero_k_default     <- 50    # Default nearest-neighbor count for make_cicero_cds
cicero_k_downsampled <- 100   # Larger k used after mouse depth matching to compensate
                               # for the reduced library complexity after thinning

# Gene biotypes retained for the guinea pig annotation.
gp_desired_biotypes <- c("protein_coding", "lncRNA", "rRNA")

###############################################################################
# STEP 1. Build genome annotations.
#
# Mouse annotations are retrieved directly from EnsDb.Mmusculus.v79 and
# reformatted to UCSC-style chromosome names (e.g., "chr1") to match the
# mm10-aligned ATAC peaks.
#
# Guinea pig annotations are imported from the NCBI GTF file for mCavPor4.1
# (GCF_034190915.1, release RS_2024_02). Because NCBI chromosome identifiers
# contain underscores (e.g., "NC_008769.1"), underscores are converted to
# hyphens to prevent Cicero from misinterpreting the identifier delimiter.
# Only protein-coding, lncRNA, and rRNA genes are retained.
###############################################################################

# ---- 1a. Mouse: mm10 gene annotations from EnsDb.Mmusculus.v79 ----

mouse_gene_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(mouse_gene_annotations) <- "UCSC"   # Convert to "chr1" style
genome(mouse_gene_annotations) <- "mm10"
saveRDS(mouse_gene_annotations,
        file.path(output_dir, "mouse_mm10_gene_annotations.rds"))

# ---- 1b. Guinea pig: gene annotations from NCBI GTF (mCavPor4.1) ----

gp_gene_annotations <- rtracklayer::import(gp_gtf_path, format = "gtf")

# NCBI chromosome names contain underscores (e.g., "NC_008769.1"). Replace
# underscores with hyphens so that Signac and Cicero can parse peak identifiers
# correctly; the genomic coordinates themselves are unaffected.
seqlevels(gp_gene_annotations) <- gsub("_", "-", seqlevels(gp_gene_annotations))

gp_gene_annotations <- gp_gene_annotations[
  mcols(gp_gene_annotations)$gene_biotype %in% gp_desired_biotypes
]

# The GTF 'gene_name' field is stored in column 10 of the mcols; rename it
# to the standard 'gene_name' slot expected by Signac's annotation functions.
names(mcols(gp_gene_annotations))[10] <- "gene_name"

genome(gp_gene_annotations) <- "GCF_034190915.1"
saveRDS(gp_gene_annotations,
        file.path(output_dir, "guinea_pig_cavpor4_gene_annotations.rds"))

###############################################################################
# STEP 2. Load RNA-integrated EC Seurat objects and extract EC barcode lists.
#
# The RNA-integration step (performed upstream in the heart EC preprocessing
# pipeline) annotated each nucleus with an endothelial cell-type label. Here,
# those labeled objects are loaded solely to extract EC barcode identifiers.
# ATAC peak counts for those barcodes are then subset from the aggregated
# multiome peak matrix.
#
# In the aggregated Cell Ranger ARC output, each sample is assigned a numeric
# barcode suffix (-1, -2, -3). The RNA Seurat object stores barcodes with
# sample-specific suffix tags appended during integration. The suffix is
# stripped and then replaced with the correct aggregation-library suffix so
# that EC barcodes can be matched against the aggregated ATAC matrix.
#
# Mouse samples and their aggregation suffixes:
#   E13 → -1  |  E15 → -2  |  E17 → -3
#
# Guinea pig samples and their aggregation suffixes:
#   GD241 → -2  |  GD245 → -1  |  GD247 → -3
###############################################################################

# ---- 2a. Mouse EC barcodes ----

mouse_ec_rna <- readRDS(mouse_ec_rna_seurat_path)
Idents(mouse_ec_rna) <- "sample"

# Extract barcodes per sample, strip the integration suffix, and apply the
# correct aggregation-library suffix.
ec_barcodes_e13 <- WhichCells(mouse_ec_rna, idents = "E13")
ec_barcodes_e13 <- sub("_.*$", "", ec_barcodes_e13)        # Strip integration tag
# E13 barcodes keep the default -1 aggregation suffix; no change needed.

ec_barcodes_e15 <- WhichCells(mouse_ec_rna, idents = "E15")
ec_barcodes_e15 <- sub("_.*$", "", ec_barcodes_e15)
ec_barcodes_e15 <- gsub("-1$", "-2", ec_barcodes_e15)      # E15 is library 2

ec_barcodes_e17 <- WhichCells(mouse_ec_rna, idents = "E17")
ec_barcodes_e17 <- sub("_.*$", "", ec_barcodes_e17)
ec_barcodes_e17 <- gsub("-1$", "-3", ec_barcodes_e17)      # E17 is library 3

mouse_ec_barcodes_all <- c(ec_barcodes_e13, ec_barcodes_e15, ec_barcodes_e17)

# ---- 2b. Guinea pig EC barcodes ----

gp_ec_rna <- readRDS(gp_ec_rna_seurat_path)
Idents(gp_ec_rna) <- "sample"

ec_barcodes_gd241 <- WhichCells(gp_ec_rna, idents = "G241")
ec_barcodes_gd241 <- sub("_.*$", "", ec_barcodes_gd241)
ec_barcodes_gd241 <- gsub("-1$", "-2", ec_barcodes_gd241)  # GD241 is library 2

ec_barcodes_gd245 <- WhichCells(gp_ec_rna, idents = "G245")
ec_barcodes_gd245 <- sub("_.*$", "", ec_barcodes_gd245)
# GD245 barcodes keep the default -1 aggregation suffix; no change needed.

ec_barcodes_gd247 <- WhichCells(gp_ec_rna, idents = "G247")
ec_barcodes_gd247 <- sub("_.*$", "", ec_barcodes_gd247)
ec_barcodes_gd247 <- gsub("-1$", "-3", ec_barcodes_gd247)  # GD247 is library 3

gp_ec_barcodes_all <- c(ec_barcodes_gd241, ec_barcodes_gd245, ec_barcodes_gd247)

###############################################################################
# STEP 3. Load aggregated ATAC peak counts and subset to EC barcodes.
#
# The aggregated Cell Ranger ARC filtered_feature_bc_matrix.h5 contains peak
# counts for all cell-associated barcodes across all three samples. Only
# barcodes that were annotated as endothelial in the RNA integration step are
# retained, and the peak matrix is restricted to the intersection with the EC
# barcode list.
#
# For mouse, peaks on non-standard chromosomes (contigs, unplaced scaffolds)
# are removed because Cicero requires chromosome lengths from the standard
# assembly and CellOracle motif scanning operates on canonical chromosomes only.
#
# For guinea pig, no standard-chromosome filter is applied because the mCavPor4.1
# assembly uses NCBI accession identifiers rather than "chr"-style names;
# the chromosome scope is controlled downstream by matching against the
# BSgenome seqInfo object.
###############################################################################

# ---- 3a. Mouse ----

mouse_atac_counts_all <- Read10X_h5(mouse_atac_h5_path)$Peaks
mouse_atac_ec_overlap  <- intersect(mouse_ec_barcodes_all, colnames(mouse_atac_counts_all))
mouse_atac_counts      <- mouse_atac_counts_all[, mouse_atac_ec_overlap]

# Restrict to peaks on standard (canonical) chromosomes.
grange_mouse <- StringToGRanges(rownames(mouse_atac_counts), sep = c(":", "-"))
standard_chr_mask <- seqnames(grange_mouse) %in% standardChromosomes(grange_mouse)
mouse_atac_counts <- mouse_atac_counts[as.vector(standard_chr_mask), ]

cat(sprintf("Mouse: retained %d peaks across %d EC nuclei\n",
            nrow(mouse_atac_counts), ncol(mouse_atac_counts)))

# ---- 3b. Guinea pig ----

gp_atac_counts_all <- Read10X_h5(gp_atac_h5_path)$Peaks
gp_atac_ec_overlap  <- intersect(gp_ec_barcodes_all, colnames(gp_atac_counts_all))
gp_atac_counts      <- gp_atac_counts_all[, gp_atac_ec_overlap]

cat(sprintf("Guinea pig: retained %d peaks across %d EC nuclei\n",
            nrow(gp_atac_counts), ncol(gp_atac_counts)))

###############################################################################
# STEP 4. Build Signac ChromatinAssay and Seurat ATAC objects.
#
# A Signac ChromatinAssay is constructed separately for each species, linking
# the EC peak count matrix to the species-specific genome and fragment file.
# Per-barcode ATAC QC metrics (atac_peak_region_fragments, passed_filters, etc.)
# are loaded from the Cell Ranger ARC per_barcode_metrics.csv files and merged
# into the Seurat object metadata so that QC filters can be applied in STEP 5.
#
# Guinea pig chromosome identifiers: because the NCBI assembly uses underscores
# in accession names (e.g., "NC_008769.1"), the seqInfo levels are converted to
# hyphen-delimited names (e.g., "NC-008769.1") before creating the assay,
# consistent with the GTF annotation reformatting in STEP 1.
###############################################################################

# ---- 4a. Mouse: per-barcode QC metadata ----

mouse_meta_e13 <- read.csv(mouse_e13_barcode_metrics, header = TRUE, row.names = 1)
mouse_meta_e15 <- read.csv(mouse_e15_barcode_metrics, header = TRUE, row.names = 1)
mouse_meta_e17 <- read.csv(mouse_e17_barcode_metrics, header = TRUE, row.names = 1)

# Apply the same suffix remapping used for EC barcodes so metadata row names
# align with the aggregated ATAC peak matrix column names.
row.names(mouse_meta_e15) <- gsub("-1$", "-2", row.names(mouse_meta_e15))
row.names(mouse_meta_e17) <- gsub("-1$", "-3", row.names(mouse_meta_e17))

mouse_meta <- rbind(mouse_meta_e13, mouse_meta_e15, mouse_meta_e17)

# ---- 4b. Mouse: Signac ChromatinAssay and Seurat object ----

mouse_chrom_assay <- CreateChromatinAssay(
  counts    = mouse_atac_counts,
  sep       = c(":", "-"),
  genome    = "mm10",
  fragments = mouse_fragment_file,
  min.cells = 1
)

mouse_atac <- CreateSeuratObject(
  counts    = mouse_chrom_assay,
  assay     = "peaks",
  meta.data = mouse_meta
)
Annotation(mouse_atac) <- mouse_gene_annotations

# ---- 4c. Guinea pig: per-barcode QC metadata ----

gp_meta_gd241 <- read.csv(gp_gd241_barcode_metrics, header = TRUE, row.names = 1)
gp_meta_gd245 <- read.csv(gp_gd245_barcode_metrics, header = TRUE, row.names = 1)
gp_meta_gd247 <- read.csv(gp_gd247_barcode_metrics, header = TRUE, row.names = 1)

row.names(gp_meta_gd241) <- gsub("-1$", "-2", row.names(gp_meta_gd241))
# GD245 retains the -1 suffix; no change needed.
row.names(gp_meta_gd247) <- gsub("-1$", "-3", row.names(gp_meta_gd247))

gp_meta <- rbind(gp_meta_gd241, gp_meta_gd245, gp_meta_gd247)

# ---- 4d. Guinea pig: Signac ChromatinAssay and Seurat object ----

# Use the BSgenome seqInfo object as the genome reference so that Signac
# receives full contig-length information for all mCavPor4.1 scaffolds.
# Underscores in NCBI chromosome names are replaced with hyphens to avoid
# conflicts with Signac's peak-identifier parsing (which uses "-" as delimiter).
gp_seq_info <- seqinfo(BSgenome.CPorcellus.NCBI.CavPor4.2)
seqlevels(gp_seq_info) <- gsub("_", "-", seqlevels(gp_seq_info))

gp_chrom_assay <- CreateChromatinAssay(
  counts    = gp_atac_counts,
  sep       = c(":", "-"),
  genome    = gp_seq_info,
  fragments = gp_fragment_file,
  min.cells = 1
)

guinea_pig_atac <- CreateSeuratObject(
  counts    = gp_chrom_assay,
  assay     = "peaks",
  meta.data = gp_meta
)

# Attach the NCBI GTF annotation; chromosome names in the GTF were already
# reformatted to hyphen-delimited in STEP 1.
Annotation(guinea_pig_atac) <- gp_gene_annotations

###############################################################################
# STEP 5. Compute ATAC quality-control metrics.
#
# Two per-nucleus chromatin-quality metrics are computed for each species:
#   - Nucleosome signal: ratio of mono-nucleosomal to sub-nucleosomal fragment
#       lengths. Low values indicate that the ATAC digestion captured open
#       chromatin accessible regions rather than wrapped nucleosomal DNA.
#   - TSS enrichment: ratio of fragment coverage at annotated transcription
#       start sites to background flanking regions. High values confirm that
#       fragments are concentrated at actively regulated genomic loci.
#
# Additional per-nucleus metrics (fraction of reads in peaks, blacklist ratio)
# are pre-computed by Cell Ranger ARC and are available in the per-barcode
# metadata loaded in STEP 4; they are not recomputed here but can be inspected
# in the object metadata for exploratory QC.
###############################################################################

# Enable multi-session parallelism for the computationally intensive QC steps.
plan("multisession", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)  # 100 GB for large peak matrices

# ---- 5a. Mouse ----
mouse_atac <- NucleosomeSignal(object = mouse_atac)
mouse_atac <- TSSEnrichment(object = mouse_atac, fast = TRUE)

# ---- 5b. Guinea pig ----
guinea_pig_atac <- NucleosomeSignal(object = guinea_pig_atac)
guinea_pig_atac <- TSSEnrichment(object = guinea_pig_atac, fast = TRUE)

plan("sequential")
gc()

###############################################################################
# STEP 6. Apply ATAC QC filters to select high-quality EC nuclei.
#
# The same three thresholds are applied to both species:
#   - atac_peak_region_fragments 2,000–25,000: retains nuclei with adequate
#       sequencing depth for coaccessibility scoring while excluding empty
#       droplets (< 2,000) and likely doublets or over-captured nuclei (> 25,000).
#   - nucleosome_signal < 2.5: selects sub-nucleosomal open-chromatin profiles;
#       nuclei with high nucleosome signal contain predominantly wrapped DNA
#       that would not reflect regulatory-element accessibility.
#   - TSS.enrichment > 2: confirms that the ATAC signal is enriched at
#       transcriptionally active regions rather than distributed as background noise.
#
# These thresholds match those reported in the scATAC-seq Methods subsection
# and are applied equivalently across species so that retained nuclei
# represent the same chromatin-quality criteria in each.
###############################################################################

# ---- 6a. Mouse ----

mouse_atac <- subset(
  x      = mouse_atac,
  subset = atac_peak_region_fragments > atac_min_fragments  &
           atac_peak_region_fragments < atac_max_fragments  &
           nucleosome_signal          < max_nucleosome_signal &
           TSS.enrichment             > min_tss_enrichment
)

cat(sprintf("Mouse: %d EC nuclei retained after ATAC QC\n", ncol(mouse_atac)))
saveRDS(mouse_atac,
        file.path(output_dir, "mouse_heart_ec_atac_signac.rds"))

# ---- 6b. Guinea pig ----

guinea_pig_atac <- subset(
  x      = guinea_pig_atac,
  subset = atac_peak_region_fragments > atac_min_fragments  &
           atac_peak_region_fragments < atac_max_fragments  &
           nucleosome_signal          < max_nucleosome_signal &
           TSS.enrichment             > min_tss_enrichment
)

cat(sprintf("Guinea pig: %d EC nuclei retained after ATAC QC\n", ncol(guinea_pig_atac)))
saveRDS(guinea_pig_atac,
        file.path(output_dir, "guinea_pig_heart_ec_atac_signac.rds"))

###############################################################################
# STEP 7. Convert QC-filtered Signac objects to Monocle3 CellDataSets.
#
# Cicero requires a Monocle3 CellDataSet as input. The Signac object is
# converted with SeuratWrappers::as.cell_data_set(). Peaks with zero total
# counts across all retained nuclei are removed before dimensionality reduction.
#
# For guinea pig, the NCBI chromosome accession format (e.g., "NW_026947484.1")
# is incompatible with Monocle3's feature-name parsing because the underscore
# separator is also used as the peak-component delimiter. To resolve this,
# the first hyphen in each peak identifier (which separates the chromosome
# accession from the start coordinate) is converted to an underscore. After
# this substitution, peak names take the form "NW_026947484.1-start-end",
# consistent with Monocle3 and Cicero expectations. A new CellDataSet is
# constructed from the reformatted count matrix so that all internal indices
# are consistent.
#
# Latent semantic indexing (LSI) followed by UMAP embedding is computed for
# each species to provide the low-dimensional coordinates used by Cicero's
# k-nearest-neighbor coaccessibility model.
###############################################################################

# ---- 7a. Mouse: Monocle3 CDS and LSI-UMAP ----

mouse_atac_cds <- as.cell_data_set(mouse_atac)
mouse_atac_cds <- mouse_atac_cds[Matrix::rowSums(exprs(mouse_atac_cds)) != 0, ]
mouse_atac_cds <- monocle3::detect_genes(mouse_atac_cds)
mouse_atac_cds <- estimate_size_factors(mouse_atac_cds)
mouse_atac_cds <- preprocess_cds(mouse_atac_cds, method = "LSI")
mouse_atac_cds <- reduce_dimension(mouse_atac_cds,
                                   reduction_method  = "UMAP",
                                   preprocess_method = "LSI")

saveRDS(mouse_atac_cds,
        file.path(output_dir, "mouse_heart_ec_atac_monocle_cds.rds"))

# ---- 7b. Guinea pig: peak-name reformatting, Monocle3 CDS, and LSI-UMAP ----

gp_atac_cds_raw <- as.cell_data_set(guinea_pig_atac)

# Extract the count matrix and metadata to reconstruct the CDS with corrected
# peak identifiers. The first hyphen in each NCBI peak name is replaced with an
# underscore (e.g., "NW-026947484.1-start-end" → "NW_026947484.1-start-end").
gp_counts_matrix <- counts(gp_atac_cds_raw)
gp_row_data      <- rowData(gp_atac_cds_raw)
gp_col_data      <- colData(gp_atac_cds_raw)

gp_peak_names_reformatted <- sub("-", "_", rownames(gp_counts_matrix), fixed = TRUE)
rownames(gp_counts_matrix) <- gp_peak_names_reformatted
rownames(gp_row_data)      <- gp_peak_names_reformatted

guinea_pig_atac_cds <- new_cell_data_set(
  gp_counts_matrix,
  cell_metadata = gp_col_data,
  gene_metadata = gp_row_data
)

guinea_pig_atac_cds <- guinea_pig_atac_cds[
  Matrix::rowSums(exprs(guinea_pig_atac_cds)) != 0,
]
guinea_pig_atac_cds <- monocle3::detect_genes(guinea_pig_atac_cds)
guinea_pig_atac_cds <- estimate_size_factors(guinea_pig_atac_cds)
guinea_pig_atac_cds <- preprocess_cds(guinea_pig_atac_cds, method = "LSI")
guinea_pig_atac_cds <- reduce_dimension(guinea_pig_atac_cds,
                                        reduction_method  = "UMAP",
                                        preprocess_method = "LSI")

saveRDS(guinea_pig_atac_cds,
        file.path(output_dir, "guinea_pig_heart_ec_atac_monocle_cds.rds"))

###############################################################################
# STEP 8. Mouse depth matching: binomial thinning to equalize in-peak fragment
#         depth before running Cicero.
#
# Because mouse multiome libraries were sequenced more deeply than guinea pig
# libraries, the median in-peak fragment count per nucleus differs between
# species. Running Cicero at the native mouse depth and then applying the same
# coaccessibility threshold (>= 0.54) as guinea pig would give unequal
# stringency across species. To address this, mouse ATAC fragment counts are
# down-sampled by binomial thinning: each observed count is replaced by a
# Binomial(count, p) draw, where p = median(guinea_pig_in_peak_depth) /
# median(mouse_in_peak_depth). This preserves the relative distribution of
# counts across peaks within each nucleus while aligning the per-nucleus
# sequencing depth to the guinea pig median.
#
# The depth-matched mouse Signac object is used only for Cicero; it is not
# saved as a final QC-filtered object and does not replace the full-depth
# mouse_heart_ec_atac_signac.rds checkpoint.
###############################################################################

# Compute global thinning probability from median in-peak depths.
mouse_in_peak_depth    <- Matrix::colSums(
  GetAssayData(mouse_atac, assay = "peaks", slot = "counts")
)
gp_in_peak_depth       <- Matrix::colSums(
  GetAssayData(guinea_pig_atac, assay = "peaks", slot = "counts")
)
p_depth_match <- min(1, median(gp_in_peak_depth) / median(mouse_in_peak_depth))
cat(sprintf(
  "Depth matching: median in-peak fragments mouse=%.0f, guinea pig=%.0f; p=%.3f\n",
  median(mouse_in_peak_depth), median(gp_in_peak_depth), p_depth_match
))

# Binomial thinning applied column-wise to the sparse count matrix.
thin_sparse_counts_by_column <- function(sparse_mat, p) {
  sparse_mat <- as(sparse_mat, "dgCMatrix")
  if (p < 1) {
    for (j in seq_len(ncol(sparse_mat))) {
      idx <- (sparse_mat@p[j] + 1L):sparse_mat@p[j + 1L]
      if (length(idx)) {
        sparse_mat@x[idx] <- rbinom(
          n    = length(idx),
          size = sparse_mat@x[idx],
          prob = p
        )
      }
    }
    sparse_mat <- drop0(sparse_mat)
  }
  sparse_mat
}

mouse_atac_for_depth_matched_cicero <- mouse_atac
set.seed(42)
mouse_counts_thinned <- thin_sparse_counts_by_column(
  GetAssayData(mouse_atac, assay = "peaks", slot = "counts"),
  p_depth_match
)

# Update the peak assay counts and recompute total in-peak depth metadata.
mouse_atac_for_depth_matched_cicero[["peaks"]]@counts        <- mouse_counts_thinned
mouse_atac_for_depth_matched_cicero$atac_peak_region_fragments <- Matrix::colSums(
  mouse_counts_thinned
)

# Rebuild the Monocle3 CDS from the depth-matched Signac object and re-run
# LSI-UMAP. A larger k (100 neighbors) is used for make_cicero_cds to
# compensate for the reduced per-nucleus count depth after thinning.
mouse_atac_cds_depth_matched <- as.cell_data_set(mouse_atac_for_depth_matched_cicero)
mouse_atac_cds_depth_matched <- mouse_atac_cds_depth_matched[
  Matrix::rowSums(exprs(mouse_atac_cds_depth_matched)) != 0,
]
mouse_atac_cds_depth_matched <- monocle3::detect_genes(mouse_atac_cds_depth_matched)
mouse_atac_cds_depth_matched <- estimate_size_factors(mouse_atac_cds_depth_matched)
mouse_atac_cds_depth_matched <- preprocess_cds(mouse_atac_cds_depth_matched, method = "LSI")
mouse_atac_cds_depth_matched <- reduce_dimension(mouse_atac_cds_depth_matched,
                                                  reduction_method  = "UMAP",
                                                  preprocess_method = "LSI")

###############################################################################
# STEP 9. Run Cicero to estimate peak-to-peak coaccessibility scores.
#
# Cicero groups nearby nuclei by k-nearest-neighbor (k = cicero_k_downsampled
# for the depth-matched mouse run; k = cicero_k_default for guinea pig) to
# form "metacells", then fits a graphical LASSO model across each chromosome
# to estimate the probability that pairs of peaks are co-accessible in the same
# nucleus. Coaccessibility scores range from -1 to 1; pairs with a score >= 0.54
# are retained downstream in CellOracle as candidate regulatory links.
#
# For mouse, the depth-matched CDS is used so that the coaccessibility model
# operates at a sequencing depth comparable to guinea pig. For guinea pig, the
# CDS prepared in STEP 7 is used directly.
#
# Chromosome lengths are required by run_cicero to define the genomic windows
# for the graphical LASSO fits. For mouse, lengths are loaded from a pre-
# downloaded mm10 chromosome length table. For guinea pig, lengths are derived
# from the BSgenome seqInfo object and restricted to chromosomes present in the
# QC-filtered peak set.
###############################################################################

# ---- 9a. Mouse: Cicero with depth-matched counts and k = 100 ----

mouse_umap_coords_depth_matched <- reducedDims(mouse_atac_cds_depth_matched)$UMAP
mouse_cicero_cds_depth_matched  <- make_cicero_cds(
  mouse_atac_cds_depth_matched,
  reduced_coordinates = mouse_umap_coords_depth_matched,
  k = cicero_k_downsampled
)

mouse_chrom_lengths <- read.table(mouse_chrom_lengths_path)

mouse_cicero_connections <- run_cicero(
  mouse_cicero_cds_depth_matched,
  mouse_chrom_lengths
)
saveRDS(mouse_cicero_connections,
        file.path(output_dir, "mouse_heart_ec_cicero_connections.rds"))

# ---- 9b. Guinea pig: Cicero with standard k ----

gp_umap_coords <- reducedDims(guinea_pig_atac_cds)$UMAP
guinea_pig_cicero_cds <- make_cicero_cds(
  guinea_pig_atac_cds,
  reduced_coordinates = gp_umap_coords,
  k = cicero_k_default
)

# Derive chromosome lengths from BSgenome; restrict to chromosomes present in
# the QC-filtered peak set to avoid errors on unrepresented scaffolds.
gp_chrom_lengths_full <- seqlengths(BSgenome.CPorcellus.NCBI.CavPor4.2)
gp_chrom_lengths_df   <- data.frame(
  V1 = names(gp_chrom_lengths_full),
  V2 = as.numeric(gp_chrom_lengths_full)
)

gp_chromosomes_in_peaks <- unique(
  sapply(strsplit(rownames(guinea_pig_atac_cds), split = "-"), `[`, 1)
)
gp_chrom_lengths_df <- gp_chrom_lengths_df[
  gp_chrom_lengths_df$V1 %in% gp_chromosomes_in_peaks,
]

guinea_pig_cicero_connections <- run_cicero(
  guinea_pig_cicero_cds,
  gp_chrom_lengths_df
)
saveRDS(guinea_pig_cicero_connections,
        file.path(output_dir, "guinea_pig_heart_ec_cicero_connections.rds"))

###############################################################################
# STEP 10. Export Cicero connections and peak lists in CellOracle format.
#
# CellOracle requires peak identifiers with underscore delimiters throughout
# (e.g., "chr1_100_500" for mouse; "NW_026947484.1_start_end" for guinea pig).
# The default Cicero output uses mixed delimiters inherited from the input CDS
# (which uses hyphens between coordinates). All hyphens in Peak1 and Peak2
# columns are replaced with underscores before export.
#
# The all_peaks list contains every peak present in the final QC-filtered CDS
# and is used by CellOracle to define the universe of regulatory elements
# against which TF motif scanning is performed.
###############################################################################

# ---- 10a. Mouse ----

mouse_cicero_connections_export         <- mouse_cicero_connections
mouse_cicero_connections_export$Peak1   <- gsub("-", "_", mouse_cicero_connections_export$Peak1)
mouse_cicero_connections_export$Peak2   <- gsub("-", "_", mouse_cicero_connections_export$Peak2)

mouse_all_peaks_export <- gsub("-", "_", row.names(exprs(mouse_atac_cds_depth_matched)))

write.csv(
  mouse_cicero_connections_export,
  file = file.path(output_dir, "mouse_heart_ec_cicero_connections.csv")
)
write.csv(
  mouse_all_peaks_export,
  file = file.path(output_dir, "mouse_heart_ec_all_peaks.csv")
)

cat(sprintf(
  "Mouse: exported %d Cicero connections and %d peaks for CellOracle.\n",
  nrow(mouse_cicero_connections_export),
  length(mouse_all_peaks_export)
))

# ---- 10b. Guinea pig ----

gp_cicero_connections_export         <- guinea_pig_cicero_connections
gp_cicero_connections_export$Peak1   <- gsub("-", "_", gp_cicero_connections_export$Peak1)
gp_cicero_connections_export$Peak2   <- gsub("-", "_", gp_cicero_connections_export$Peak2)

gp_all_peaks_export <- gsub("-", "_", row.names(exprs(guinea_pig_atac_cds)))

write.csv(
  gp_cicero_connections_export,
  file = file.path(output_dir, "guinea_pig_heart_ec_cicero_connections.csv")
)
write.csv(
  gp_all_peaks_export,
  file = file.path(output_dir, "guinea_pig_heart_ec_all_peaks.csv")
)

cat(sprintf(
  "Guinea pig: exported %d Cicero connections and %d peaks for CellOracle.\n",
  nrow(gp_cicero_connections_export),
  length(gp_all_peaks_export)
))

cat("\nATAC Signac/Cicero preprocessing complete.\n")
cat(sprintf("Output directory: %s\n", output_dir))
