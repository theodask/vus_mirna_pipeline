#!/bin/bash

# Script to calculate seed regions from MRE coordinates
# Adds a 150bp window (±75bp) around the seed region (positions 2-7 of the MRE)

set -e  # Exit on error
set -u  # Exit on undefined variable

# Parse command-line arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bed> <output_bed> <cancer_type>"
    echo ""
    echo "Arguments:"
    echo "  input_bed    : Input BED file with MRE coordinates"
    echo "  output_bed   : Output BED file with seed region windows"
    echo "  cancer_type  : Cancer type identifier (for logging)"
    exit 1
fi

INPUT_BED="$1"
OUTPUT_BED="$2"
CANCER_TYPE="$3"

echo "Calculating Seed Regions"
echo "Cancer type  : ${CANCER_TYPE}"
echo "Input BED    : ${INPUT_BED}"
echo "Output BED   : ${OUTPUT_BED}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_BED}" ]; then
    echo "ERROR: Input BED file not found: ${INPUT_BED}" >&2
    exit 1
fi

# Count input lines
input_lines=$(wc -l < "${INPUT_BED}")
echo "Processing ${input_lines} MRE regions..."

# Calculate seed regions with 150bp window (±75bp from center)
awk 'BEGIN{OFS="\t"} {
    chr = $1;
    start = $2;
    
    # Seed region is positions 2-7 of the MRE (1-indexed in biology, 0-indexed in BED)
    seed_start = start + 1;
    seed_end = seed_start + 6;
    
    # Calculate center of seed region
    center = int((seed_start + seed_end) / 2);
    
    # Create 150bp window around center (±75bp)
    window_start = center - 75;
    window_end = center + 75;
    
    # Ensure we don'\''t go below 0
    if (window_start < 0) window_start = 0;
    
    print chr, window_start, window_end;
}' "${INPUT_BED}" > "${OUTPUT_BED}"

# Count output lines
output_lines=$(wc -l < "${OUTPUT_BED}")

echo ""
echo "Completed."
echo "Processed: ${input_lines} regions"
echo "Output:    ${output_lines} seed windows"
echo "Saved to:  ${OUTPUT_BED}"

