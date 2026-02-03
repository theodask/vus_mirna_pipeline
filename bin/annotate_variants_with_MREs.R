#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(readr)
  library(rtracklayer)
  library(optparse)
})

option_list <- list(
  make_option(c("--cancer_type"), type="character", default=NULL,
              help="Cancer type identifier", metavar="character"),
  make_option(c("--input_bed"), type="character", default=NULL,
              help="Input BED file with ClinVar-annotated variants", metavar="character"),
  make_option(c("--mres_file"), type="character", default=NULL,
              help="MREs annotation file", metavar="character"),
  make_option(c("--output_annotated"), type="character", default=NULL,
              help="Output BED file with MRE annotations", metavar="character"),
  make_option(c("--output_microt"), type="character", default=NULL,
              help="Output BED file for microT-CNN (no header)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cancer_type) || is.null(opt$input_bed) || 
    is.null(opt$mres_file) || is.null(opt$output_annotated) || 
    is.null(opt$output_microt)) {
  print_help(opt_parser)
  stop("All required arguments must be provided.", call.=FALSE)
}

cancer_type <- opt$cancer_type
input_bed_file <- opt$input_bed
mres_file <- opt$mres_file
output_annotated_file <- opt$output_annotated
output_microt_file <- opt$output_microt

cat("Annotating Variants with MREs\n")
cat("Cancer type         : ", cancer_type, "\n")
cat("Input BED           : ", input_bed_file, "\n")
cat("MREs file           : ", mres_file, "\n")
cat("Output annotated    : ", output_annotated_file, "\n")
cat("Output microT-CNN   : ", output_microt_file, "\n\n")

# Check if input files exist
if (!file.exists(input_bed_file)) {
  stop("Input BED file not found: ", input_bed_file, call.=FALSE)
}
if (!file.exists(mres_file)) {
  stop("MREs file not found: ", mres_file, call.=FALSE)
}

# Load already-filtered variants
cat("Loading variants...\n")
variants_df <- fread(input_bed_file)
cat("Loaded ", nrow(variants_df), " variants\n\n")

# Convert to GRanges (BED is 0-based, GRanges is 1-based)
cat("Converting variants to GRanges...\n")
variants_gr <- GRanges(
  seqnames = variants_df$CHROM,
  ranges = IRanges(start = variants_df$START + 1, end = variants_df$END),
  REF = variants_df$REF,
  ALT = variants_df$ALT,
  Patient_Hosts = variants_df$Patient_Hosts,
  DUMMY = variants_df$DUMMY,
  ALLELE = variants_df$ALLELE,
  CLINVAR_ID = variants_df$CLINVAR_ID
)
seqlevelsStyle(variants_gr) <- "UCSC"  # Match with MREs

# Load MREs
cat("Loading MREs...\n")
mres_filtered <- fread(mres_file)
cat("Loaded ", nrow(mres_filtered), " MREs\n\n")

# Convert MREs to GRanges
cat("Converting MREs to GRanges...\n")
mres_gr <- GRanges(
  seqnames = mres_filtered$chromosome,
  ranges = IRanges(start = mres_filtered$start, end = mres_filtered$end),
  mirna_name = mres_filtered$mirna_name,
  cds_utr = mres_filtered$cds_utr
)
seqlevelsStyle(mres_gr) <- "UCSC"

# Find overlaps
cat("Finding overlaps between variants and MREs...\n")
overlaps <- findOverlaps(variants_gr, mres_gr, ignore.strand = TRUE)
cat("Found ", length(overlaps), " overlaps\n\n")

if (length(overlaps) == 0) {
  stop("No overlaps found between variants and MREs. Check chromosome naming conventions.", call.=FALSE)
}

# Extract matching rows
variant_hits <- variants_df[queryHits(overlaps), ]
mres_hits <- mres_filtered[subjectHits(overlaps), ]

# Annotate variants with MRE information
cat("Annotating variants with MRE information...\n")
seed_hits <- variant_hits %>%
  mutate(
    miRNA = mres_hits$mirna_name,
    MRE_START = mres_hits$start,
    MRE_END = mres_hits$end,
    POS_IN_MRE = START - mres_hits$start + 1,
    CDS_UTR = mres_hits$cds_utr
  )

# Reorder columns for output 1
seed_hits <- seed_hits[, .(
  CHROM, START, END, DUMMY, ALLELE, REF, ALT, CLINVAR_ID,
  miRNA, MRE_START, MRE_END, POS_IN_MRE, CDS_UTR,
  Patient_Hosts
)]

# Create output directories if needed
for (output_file in c(output_annotated_file, output_microt_file)) {
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory: ", output_dir, "\n")
  }
}

# Full annotated results (with header)
cat("\nWriting annotated variants BED...\n")
write_tsv(seed_hits, output_annotated_file)
cat("Annotated variant BED saved to: ", output_annotated_file, "\n")

# microT-CNN input (distinct core columns, no header)
cat("Preparing microT-CNN input...\n")
microt_df <- dplyr::distinct(
  seed_hits, 
  CHROM, START, END, DUMMY, ALLELE, REF, ALT, CLINVAR_ID, 
  .keep_all = FALSE
)
write_tsv(microt_df, output_microt_file, col_names = FALSE)
cat("microT-CNN input saved to: ", output_microt_file, "\n")

cat("Completed.\n")
cat("Total input variants:        ", nrow(variants_df), "\n")
cat("Variants with MRE overlaps:  ", length(unique(queryHits(overlaps))), "\n")
cat("Total variant-MRE pairs:     ", nrow(seed_hits), "\n")
cat("Unique miRNAs:               ", length(unique(seed_hits$miRNA)), "\n")
cat("Distinct variants (microT):  ", nrow(microt_df), "\n")

