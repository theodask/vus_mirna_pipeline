#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(VariantAnnotation)
  library(data.table)
  library(readr)
  library(optparse)
})

option_list <- list(
  make_option(c("--cancer_type"), type="character", default=NULL,
              help="Cancer type identifier", metavar="character"),
  make_option(c("--input_bed"), type="character", default=NULL,
              help="Input BED file with filtered variants", metavar="character"),
  make_option(c("--clinvar_vcf"), type="character", default=NULL,
              help="ClinVar VCF file path", metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output BED file with ClinVar annotations", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cancer_type) || is.null(opt$input_bed) || 
    is.null(opt$clinvar_vcf) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All required arguments must be provided.", call.=FALSE)
}

cancer_type <- opt$cancer_type
input_bed_file <- opt$input_bed
clinvar_vcf_file <- opt$clinvar_vcf
output_file <- opt$output

cat("ClinVar Intersection Analysis\n")
cat("Cancer type    : ", cancer_type, "\n")
cat("Input BED      : ", input_bed_file, "\n")
cat("ClinVar VCF    : ", clinvar_vcf_file, "\n")
cat("Output file    : ", output_file, "\n\n")

# Check if input files exist
if (!file.exists(input_bed_file)) {
  stop("Input BED file not found: ", input_bed_file, call.=FALSE)
}
if (!file.exists(clinvar_vcf_file)) {
  stop("ClinVar VCF file not found: ", clinvar_vcf_file, call.=FALSE)
}

cat("Reading input variants...\n")
my_variants <- fread(input_bed_file, header = TRUE)
total_input_variants <- nrow(my_variants)
cat("Total input variants: ", total_input_variants, "\n\n")

# Strip 'chr' from CHROM if needed (ClinVar uses 1, 2, MT instead of chr1, chr2, chrM)
my_variants$CHROM <- sub("^chr", "", my_variants$CHROM)
my_variants$CHROM[my_variants$CHROM == "M"] <- "MT"

# Convert BED to GRanges (GRanges is 1-based, BED is 0-based)
cat("Converting variants to GRanges...\n")
gr_my <- GRanges(
  seqnames = my_variants$CHROM,
  ranges = IRanges(start = my_variants$START + 1, end = my_variants$END),
  REF = my_variants$REF,
  ALT = my_variants$ALT,
  Patient_Hosts = my_variants$Patient_Hosts
)

# Read ClinVar VCF
cat("Reading ClinVar VCF file...\n")
vcf <- readVcf(clinvar_vcf_file, "hg38")
total_clinvar_variants <- nrow(vcf)
cat("Total variants in ClinVar file: ", total_clinvar_variants, "\n")

# Keep only single ALT allele variants
is_single_alt <- elementNROWS(alt(vcf)) == 1
vcf <- vcf[is_single_alt]

# Extract REF and ALT
ref_alleles <- as.character(ref(vcf))
alt_alleles <- as.character(unlist(alt(vcf)))

# Keep SNVs only (exclude indels, structural variants, etc.)
is_snv <- nchar(ref_alleles) == 1 & nchar(alt_alleles) == 1
snv_count <- sum(is_snv)
filtered_out <- total_clinvar_variants - snv_count

cat("SNVs kept: ", snv_count, "\n")
cat("Non-SNV variants filtered out: ", filtered_out, "\n\n")

vcf <- vcf[is_snv]

# Get GRanges from ClinVar VCF
gr_clinvar <- rowRanges(vcf)
mcols(gr_clinvar)$REF <- ref_alleles[is_snv]
mcols(gr_clinvar)$ALT <- alt_alleles[is_snv]

# Intersect by genomic coordinates
cat("Finding overlapping variants...\n")
hits <- findOverlaps(gr_my, gr_clinvar)

# Filter by matching REF and ALT alleles
cat("Filtering by matching REF and ALT alleles...\n")
matches <- hits[
  mcols(gr_my)[queryHits(hits), "REF"] == mcols(gr_clinvar)[subjectHits(hits), "REF"] &
  mcols(gr_my)[queryHits(hits), "ALT"] == mcols(gr_clinvar)[subjectHits(hits), "ALT"]
]

# Get original rows that matched
filtered_indices <- queryHits(matches)
filtered_df <- my_variants[filtered_indices, ]

# Extract ClinVar IDs
clinvar_ids <- names(gr_clinvar)
matched_ids <- clinvar_ids[subjectHits(matches)]
filtered_df$CLINVAR_ID <- matched_ids

# Restore 'chr' prefix in CHROM
filtered_df$CHROM <- paste0("chr", filtered_df$CHROM)
filtered_df$CHROM <- sub("^chrMT$", "chrM", filtered_df$CHROM)

# Reorder columns
filtered_df <- filtered_df[, c("CHROM", "START", "END", "DUMMY", "ALLELE",
                               "REF", "ALT", "CLINVAR_ID", "Patient_Hosts")]

# Report statistics
matched_variants <- length(unique(queryHits(matches)))
cat("RESULTS\n")
cat("Total input variants:           ", total_input_variants, "\n")
cat("Variants matched in ClinVar:    ", matched_variants, "\n")
cat("Percentage matched:             ", 
    sprintf("%.2f%%", (matched_variants/total_input_variants)*100), "\n")
cat("Output variants with ClinVar ID:", nrow(filtered_df), "\n")

# Create output directory if needed
output_dir <- dirname(output_file)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory: ", output_dir, "\n")
}

write_tsv(filtered_df, output_file)

cat("Completed.\n")
cat("Output saved to: ", output_file, "\n")

