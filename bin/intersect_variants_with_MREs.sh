#!/bin/bash

# Script to intersect VCF variants with MRE seed regions using bedtools
# Preserves VCF header and outputs only variants that overlap seed regions

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Catch errors in pipes

# Parse command-line arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_vcf> <seed_bed> <output_vcf> <cancer_type>"
    echo ""
    echo "Arguments:"
    echo "  input_vcf   : Input VCF file (can be gzipped)"
    echo "  seed_bed    : BED file with seed MRE coordinates"
    echo "  output_vcf  : Output VCF file (will be gzipped)"
    echo "  cancer_type : Cancer type identifier (for logging)"
    exit 1
fi

INPUT_VCF="$1"
SEED_BED="$2"
OUTPUT_VCF="$3"
CANCER_TYPE="$4"

echo "========================================="
echo "Intersecting Variants with MRE Seed Regions"
echo "========================================="
echo "Cancer type     : ${CANCER_TYPE}"
echo "Input VCF       : ${INPUT_VCF}"
echo "Seed regions BED: ${SEED_BED}"
echo "Output VCF      : ${OUTPUT_VCF}"
echo ""

# Check if input files exist
if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input VCF file not found: ${INPUT_VCF}" >&2
    exit 1
fi

if [ ! -f "${SEED_BED}" ]; then
    echo "ERROR: Seed regions BED file not found: ${SEED_BED}" >&2
    exit 1
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools is not installed or not in PATH" >&2
    exit 1
fi

# Count total variants in input VCF
total_variants=$(zcat "${INPUT_VCF}" | grep -v "^#" | wc -l)
echo "Total variants in input VCF: ${total_variants}"

# Count seed regions
total_seeds=$(wc -l < "${SEED_BED}")
echo "Total seed regions: ${total_seeds}"
echo ""
echo "Running bedtools intersect..."

# Perform intersection:
# 1. Extract VCF header using zcat and grep
# 2. Use bedtools intersect to find overlapping variants
#    -a: VCF file (variants)
#    -b: BED file (seed regions)
#    -wa: Write original entry from A (variant) for each overlap
#    -u: Write unique entries from A (report each variant only once)
# 3. Compress output with gzip
(
    zcat "${INPUT_VCF}" | grep '^#'
    bedtools intersect \
        -a "${INPUT_VCF}" \
        -b "${SEED_BED}" \
        -wa -u
) | gzip > "${OUTPUT_VCF}"

# Count variants in output
intersecting_variants=$(zcat "${OUTPUT_VCF}" | grep -v "^#" | wc -l)

echo ""
echo "Completed."
echo "Total variants:        ${total_variants}"
echo "Intersecting variants: ${intersecting_variants}"
echo "Percentage:            $(awk "BEGIN {printf \"%.2f\", ($intersecting_variants/$total_variants)*100}")%"
echo "Output saved to:       ${OUTPUT_VCF}"

