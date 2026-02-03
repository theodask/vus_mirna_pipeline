#!/bin/bash

# Script to help users set up the VUS-miRNA pipeline directory structure

set -e  # Exit on error

PIPELINE_DIR="${1:-.}"

echo "VUS-miRNA Pipeline Setup"
echo "Setting up directory structure in: $PIPELINE_DIR"
echo ""

# Create directory structure
mkdir -p "$PIPELINE_DIR"/{modules,bin,resources/{aggregated_variants,chip_seq_atlas,genomic_annotation,microRNA_annotation,miTED,variants,metadata},results,microt_cnn}

echo "Directory structure created"
echo ""
echo "Next Steps:"
echo ""
echo "1. Download required reference files:"
echo ""
echo "   Genomic Annotation Files:"
echo "   - GTF: Homo_sapiens.GRCh38.113.gtf"
echo "   - FAI: Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
echo "   Place in: $PIPELINE_DIR/resources/genomic_annotation/"
echo ""
echo "   microRNA Annotation:"
echo "   - TarBase MREs file"
echo "   Place in: $PIPELINE_DIR/resources/microRNA_annotation/"
echo ""
echo "   miTED Expression Data:"
echo "   - Cancer-specific miRNA expression files"
echo "   Place in: $PIPELINE_DIR/resources/miTED/[cancer_type]/"
echo ""
echo "   ClinVar Database:"
echo "   - clinvar.vcf"
echo "   Place in: $PIPELINE_DIR/resources/variants/"
echo ""
echo "   Clinical Metadata:"
echo "   - TCGA clinical data (tar.gz files)"
echo "   Place in: $PIPELINE_DIR/resources/metadata/[cancer_type]/"
echo ""
echo "   microT-CNN References:"
echo "   - hsa_miRBase_22v1_0based.gff3"
echo "   - Homo_sapiens.GRCh38.113.gtf"
echo "   Place in: $PIPELINE_DIR/microt_cnn/"
echo ""
echo "2. Test the pipeline:"
echo "   nextflow run main.nf -profile docker --help"
echo ""
echo "3. Run with your data:"
echo "   nextflow run main.nf -profile docker \\"
echo "       --cancer_type BRCA \\"
echo "       --vcf /path/to/merged_BRCA.vcf.gz \\"
echo "       --skip_genomic_regions true \\"
echo "       --rpm_threshold 100"
echo ""
echo "For detailed documentation, see README.md"
