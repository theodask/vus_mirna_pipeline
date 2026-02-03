#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(data.table)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("--cancer_type"), type="character", default=NULL,
              help="Cancer type identifier", metavar="character"),
  make_option(c("--mited_file"), type="character", default=NULL,
              help="Path to miTED expression file", metavar="character"),
  make_option(c("--mres_file"), type="character", default=NULL,
              help="Path to MREs annotation file", metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output BED file path", metavar="character"),
  make_option(c("--rpm_threshold"), type="numeric", default=150,
              help="RPM threshold for filtering expressed miRNAs [default=%default]", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cancer_type) || is.null(opt$mited_file) || is.null(opt$mres_file) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All required arguments must be provided.", call.=FALSE)
}

cancer_type <- opt$cancer_type
mited_file <- opt$mited_file
mres_file <- opt$mres_file
output_file <- opt$output
rpm_threshold <- opt$rpm_threshold

cat("Selecting MREs from Top Expressed miRNAs\n")
cat("Cancer type:", cancer_type, "\n")
cat("miTED file:", mited_file, "\n")
cat("MREs file:", mres_file, "\n")
cat("RPM threshold:", rpm_threshold, "\n")
cat("Output file:", output_file, "\n\n")

# Check if input files exist
if (!file.exists(mited_file)) {
  stop("miTED file not found: ", mited_file, call.=FALSE)
}
if (!file.exists(mres_file)) {
  stop("MREs file not found: ", mres_file, call.=FALSE)
}

# Read miTED data and filter for highly expressed miRNAs
cat("Reading miTED expression data...\n")
mited_df <- fread(mited_file)
cat("Total miRNAs in miTED file:", nrow(mited_df), "\n")

expressed_mirnas <- mited_df %>%
  filter(expr > rpm_threshold) %>%
  pull(mirna)

cat("Found", length(expressed_mirnas), "expressed miRNAs in", cancer_type, 
    "with RPM >", rpm_threshold, "\n\n")

if (length(expressed_mirnas) == 0) {
  stop("No RPM threshold was provided.", call.=FALSE)
}

# Read and filter MREs based on highly expressed miRNAs
cat("Reading MREs annotation file...\n")
mres_df <- fread(mres_file)
cat("Total MREs in annotation file:", nrow(mres_df), "\n")

cat("Filtering MREs for expressed miRNAs...\n")
mres_filtered <- mres_df %>%
  filter(mirna_name %in% expressed_mirnas)

cat("Filtered MREs:", nrow(mres_filtered), "\n\n")

if (nrow(mres_filtered) == 0) {
  stop("No MREs found for the expressed miRNAs. Check miRNA name matching.", 
       call.=FALSE)
}

# Select MRE coordinates for output (BED format: chr, start, end)
mre_coordinates <- mres_filtered %>%
  dplyr::select(
    chrom = chromosome,
    start = start,
    end = end
  )

# Create output directory if it doesn't exist
output_dir <- dirname(output_file)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

write_delim(mre_coordinates, output_file, delim = "\t", col_names = FALSE)

cat("MRE coordinates for top expressed miRNAs saved to:\n")
cat(output_file, "\n")
