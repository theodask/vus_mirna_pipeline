process CALCULATE_SEED_REGIONS {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/seed_regions", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(mre_bed), path(script_file)
    
    output:
    tuple val(cancer_type), path("seed_mre_coordinates_top_miRNA_expressed_in_mited_${cancer_type}.bed"), emit: seed_bed
    
    script:
    """
    echo "Calculating seed regions for cancer type: ${cancer_type}"
    echo "Input MRE BED: ${mre_bed}"
    
    bash ${script_file} \\
        ${mre_bed} \\
        seed_mre_coordinates_top_miRNA_expressed_in_mited_${cancer_type}.bed \\
        ${cancer_type}
    """
}
