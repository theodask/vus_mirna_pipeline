#!/bin/bash

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Catch errors in pipes

# Parse command-line arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_vcf> <output_bed> <cancer_type> <lines_per_chunk> <r_script>"
    echo ""
    echo "Arguments:"
    echo "  input_vcf       : Input VCF file (gzipped)"
    echo "  output_bed      : Output BED file path"
    echo "  cancer_type     : Cancer type identifier"
    echo "  lines_per_chunk : Number of variants per chunk (e.g., 200000)"
    echo "  r_script        : Path to process_vcf_chunk.R script"
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_BED="$2"
CANCER_TYPE="$3"
LINES_PER_CHUNK="$4"
R_SCRIPT="$5"

echo "Processing VCF in Chunks"
echo "Cancer type      : ${CANCER_TYPE}"
echo "Input VCF        : ${INPUT_VCF}"
echo "Output BED       : ${OUTPUT_BED}"
echo "Lines per chunk  : ${LINES_PER_CHUNK}"
echo "R script         : ${R_SCRIPT}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input VCF file not found: ${INPUT_VCF}" >&2
    exit 1
fi

# Check if R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R script not found: ${R_SCRIPT}" >&2
    exit 1
fi

# Create temporary directories in current working directory
# Nextflow will manage these in the work directory
TEMP_VCF_DIR="./temp_vcf_chunks"
TEMP_BED_DIR="./temp_bed_chunks"

mkdir -p "${TEMP_VCF_DIR}"
mkdir -p "${TEMP_BED_DIR}"

# Define the BED header
BED_HEADER="CHROM\tSTART\tEND\tDUMMY\tALLELE\tREF\tALT\tPatient_Hosts"

echo "Created temporary directories:"
echo "  VCF chunks: ${TEMP_VCF_DIR}"
echo "  BED chunks: ${TEMP_BED_DIR}"
echo ""

# Extract VCF header
echo "Extracting VCF header."
VCF_HEADER=$(zcat "${INPUT_VCF}" | grep -E "^##|^#CHROM")
HEADER_LINES=$(echo "${VCF_HEADER}" | wc -l)
echo "Header has ${HEADER_LINES} lines."

# Count total variants
TOTAL_VARIANTS=$(zcat "${INPUT_VCF}" | grep -v "^#" | wc -l)
echo "Total variants in VCF: ${TOTAL_VARIANTS}"
echo ""

# Split the VCF file into chunks
echo "Splitting VCF file into chunks of ${LINES_PER_CHUNK} lines."
zcat "${INPUT_VCF}" | tail -n +$((HEADER_LINES + 1)) | split -l "${LINES_PER_CHUNK}" - "${TEMP_VCF_DIR}/${CANCER_TYPE}_chunk_"

# Count chunks created
CHUNK_COUNT=$(ls -1 "${TEMP_VCF_DIR}"/${CANCER_TYPE}_chunk_* 2>/dev/null | wc -l)
echo "Created ${CHUNK_COUNT} chunks."
echo ""

# Process each chunk
echo "Processing each VCF chunk with R."
CHUNK_NUM=0
for chunk_file in "${TEMP_VCF_DIR}"/${CANCER_TYPE}_chunk_*; do
    CHUNK_NUM=$((CHUNK_NUM + 1))
    echo "Processing chunk ${CHUNK_NUM}/${CHUNK_COUNT}: $(basename ${chunk_file})..."
    
    # Re-add VCF header to each chunk
    TEMP_CHUNK_WITH_HEADER="${chunk_file}.tmp.vcf"
    echo "${VCF_HEADER}" > "${TEMP_CHUNK_WITH_HEADER}"
    cat "${chunk_file}" >> "${TEMP_CHUNK_WITH_HEADER}"
    
    # Run the R script for this chunk
    Rscript "${R_SCRIPT}" "${TEMP_CHUNK_WITH_HEADER}" "${TEMP_BED_DIR}" "${CANCER_TYPE}"
    
    # Clean up the temporary chunk with header
    rm "${TEMP_CHUNK_WITH_HEADER}"
done

echo ""
echo "Concatenating all temporary BED files."

# Create temporary file for data without header
TEMP_FINAL_BED_NO_HEADER="${OUTPUT_BED}.tmp_no_header"

# Concatenate all BED chunks (header-less)
cat "${TEMP_BED_DIR}"/temp_chunk_*.bed > "${TEMP_FINAL_BED_NO_HEADER}"

# Count lines in concatenated file
LINES_IN_CONCAT=$(wc -l < "${TEMP_FINAL_BED_NO_HEADER}")
echo "Lines in concatenated data: ${LINES_IN_CONCAT}"

# Prepend the header to create final BED file
echo "Creating final BED file with header..."
echo -e "${BED_HEADER}" > "${OUTPUT_BED}"
cat "${TEMP_FINAL_BED_NO_HEADER}" >> "${OUTPUT_BED}"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "${TEMP_VCF_DIR}"
rm -rf "${TEMP_BED_DIR}"
rm -f "${TEMP_FINAL_BED_NO_HEADER}"

# Final count
TOTAL_LINES=$(wc -l < "${OUTPUT_BED}")
TOTAL_VARIANTS_BED=$((TOTAL_LINES - 1))  # Subtract header

echo ""
echo "Completed."
echo "Total lines in final BED (with header): ${TOTAL_LINES}"
echo "Total variants in BED: ${TOTAL_VARIANTS_BED}"
echo "Output saved to: ${OUTPUT_BED}"
echo ""
echo "First 5 lines of output:"
head -n 5 "${OUTPUT_BED}"

