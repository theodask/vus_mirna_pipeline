process SELECT_MRES {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/mre_selection", mode: 'copy'
    
    input:
    tuple path(script_file), val(cancer_type), path(mited_file), path(mres_file)
    
    output:
    tuple val(cancer_type), path("mre_coordinates_top_miRNA_expressed_in_mited_${cancer_type}.bed"), emit: mre_bed
    
    script:
    """
    echo "Selecting MREs for cancer type: ${cancer_type}"
    echo "miTED file: ${mited_file}"
    echo "MREs file: ${mres_file}"
    
    Rscript ${script_file} \\
        --cancer_type ${cancer_type} \\
        --mited_file ${mited_file} \\
        --mres_file ${mres_file} \\
        --output mre_coordinates_top_miRNA_expressed_in_mited_${cancer_type}.bed \\
        --rpm_threshold ${params.rpm_threshold}
    """
}
