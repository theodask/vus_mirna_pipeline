process FIND_GENOMIC_REGIONS {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/variant_regions", mode: 'copy'
    
    input:
    tuple path(script_file), val(cancer_type), path(vcf_file), path(gtf_file), path(fai_file)
    
    output:
    tuple val(cancer_type), path("variant_regions_${cancer_type}.tsv"), emit: regions
    tuple val(cancer_type), path("variant_regions_${cancer_type}_summary.tsv"), emit: summary
    tuple val(cancer_type), path("variant_density_barplot_log_scaled.png"), emit: plot
    
    script:
    """
    echo "Processing cancer type: ${cancer_type}"
    echo "VCF file: ${vcf_file}"
    echo "GTF file: ${gtf_file}"
    echo "FAI file: ${fai_file}"
    
    python3 ${script_file} \\
        --vcf ${vcf_file} \\
        --gtf ${gtf_file} \\
        --fai ${fai_file} \\
        --cancer_type ${cancer_type} \\
        --output variant_regions_${cancer_type}.tsv \\
        --summary variant_regions_${cancer_type}_summary.tsv \\
        --plot variant_density_barplot_log_scaled.png
    """
}

