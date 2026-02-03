process PROCESS_VCF_CHUNKS {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/filtered_variants", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(intersected_vcf), path(bash_script), path(r_script)
    
    output:
    tuple val(cancer_type), path("TCGA_filtered_variants_${cancer_type}.bed"), emit: filtered_bed
    
    script:
    """
    echo "Processing and filtering VCF for cancer type: ${cancer_type}"
    echo "Input VCF: ${intersected_vcf}"
    echo "Lines per chunk: ${params.lines_per_chunk}"
    
    bash ${bash_script} \\
        ${intersected_vcf} \\
        TCGA_filtered_variants_${cancer_type}.bed \\
        ${cancer_type} \\
        ${params.lines_per_chunk} \\
        ${r_script}
    """
}

