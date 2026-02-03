#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# Function to get patient IDs who are carriers
get_hosting_patients_robust <- function(row_of_genotypes_with_names) {
  
  # Initialize a character vector to store patient IDs who are carriers
  carriers <- character(0)
  
  # Get the patient IDs from the names of the input vector
  patient_ids <- names(row_of_genotypes_with_names)
  
  # Loop through each patient's genotype string in the current row
  for (i in seq_along(row_of_genotypes_with_names)) {
    full_gt_string <- row_of_genotypes_with_names[i]
    current_patient_id <- patient_ids[i]
    
    # Extract the GT field (the first part before the first colon)
    gt_field_parts <- strsplit(full_gt_string, ":")[[1]]
    gt <- gt_field_parts[1]
    
    # Define common missing genotypes
    is_missing <- gt %in% c("./.", ".|.", ".")
    
    # Check if it's a carrier of the alternate allele
    is_carrier <- !is_missing && grepl("1", gt)
    
    if (is_carrier) {
      carriers <- c(carriers, current_patient_id)
    }
  }
  
  # Return comma-separated patient IDs or NA if no carriers
  if (length(carriers) > 0) {
    return(paste(carriers, collapse = ","))
  } else {
    return(NA_character_)
  }
}

# Define a function to process a single VCF file
process_single_vcf <- function(input_vcf_path, output_bed_dir, cancer_type) {
  # Infer a unique identifier for the chunk from the input_vcf_path
  chunk_id <- tools::file_path_sans_ext(basename(input_vcf_path))
  
  cat(paste0("Processing chunk: ", input_vcf_path, "\n"))
  
  # Read the VCF file, skipping header lines until #CHROM
  vcf_df <- fread(input_vcf_path, skip = "#CHROM")
  
  # Extract INFO fields
  vcf_df <- vcf_df %>%
    mutate(
      DP = as.numeric(str_extract(INFO, "DP=[0-9]+") %>% str_replace("DP=", "")),
      GERMQ = as.numeric(str_extract(INFO, "GERMQ=[0-9]+") %>% str_replace("GERMQ=", "")),
      TLOD = as.numeric(str_extract(INFO, "TLOD=[0-9.eE+-]+") %>% str_replace("TLOD=", "")),
      AF = as.numeric(str_extract(INFO, "AF=[0-9.]+") %>% str_replace("AF=", "")),
      FILTER_PASS = FILTER == "PASS"
    )
  
  # Apply filtering criteria
  vcf_filtered <- vcf_df %>%
    filter(
      DP >= 10,
      GERMQ >= 20,
      AF >= 0.05,
      FILTER_PASS
    )
  
  # Identify patient columns (assumes all patient columns start with "TCGA")
  patient_cols <- grep("^TCGA", colnames(vcf_filtered), value = TRUE)
  
  # Apply function to get hosting patients
  if (length(patient_cols) > 0 && nrow(vcf_filtered) > 0) {
    vcf_filtered[, Patient_Hosts := apply(.SD, 1, get_hosting_patients_robust), .SDcols = patient_cols]
  } else {
    vcf_filtered[, Patient_Hosts := NA_character_]
  }
  
  # Remove the original patient columns
  vcf_filtered[, (patient_cols) := NULL]
  
  # Rename "#CHROM" to "CHROM"
  setnames(vcf_filtered, old = "#CHROM", new = "CHROM")
  
  # Choose columns and replace POS to START
  bed_df <- vcf_filtered[, .(CHROM, START = POS, REF, ALT, Patient_Hosts)]
  
  # Remove multi-allelic ALT values
  bed_df <- bed_df[!grepl(",", ALT)]
  
  # Add extra columns
  bed_df[, END := START]
  bed_df[, DUMMY := "*"]
  bed_df[, ALLELE := "+"]
  
  # Reorder columns and convert bed file to 0-based
  bed_df <- bed_df[, .(CHROM, START = START - 1, END, DUMMY, ALLELE, REF, ALT, Patient_Hosts)]
  
  output_bed_file <- file.path(output_bed_dir, paste0("temp_chunk_", chunk_id, "_", cancer_type, ".bed"))
  fwrite(bed_df, output_bed_file, sep = "\t", col.names = FALSE)
  
  cat(paste0("Processed chunk saved to: ", output_bed_file, "\n"))
}

# Command line interface
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  input_vcf_path <- args[1]
  output_bed_dir <- args[2]
  cancer_type <- args[3]
  process_single_vcf(input_vcf_path, output_bed_dir, cancer_type)
} else {
  cat("Usage: Rscript process_vcf_chunk.R <input_vcf_path> <output_bed_dir> <cancer_type>\n")
  quit(status = 1)
}

