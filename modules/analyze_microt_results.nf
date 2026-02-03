process ANALYZE_MICROT_RESULTS {
    tag "${cancer_type}"
    publishDir "${params.outdir}/${cancer_type}/microt_analysis", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(wildtype_results), path(mutated_results), path(analysis_script)
    
    output:
    tuple val(cancer_type), path("only_mutated_MREs_${cancer_type}.tsv"), emit: only_mutated
    tuple val(cancer_type), path("only_wildtype_MREs_${cancer_type}.tsv"), emit: only_wildtype
    tuple val(cancer_type), path("microT_CNN_wt_mut_results_${cancer_type}.tsv"), emit: disrupted_mres
    tuple val(cancer_type), path("histogram_MREs_disruption_scores.tiff"), emit: hist_disruption
    tuple val(cancer_type), path("histogram_gained_MREs_scores.tiff"), emit: hist_gained
    tuple val(cancer_type), path("histogram_lost_MREs_scores.tiff"), emit: hist_lost
    
    script:
    """
    Rscript ${analysis_script} \\
        --cancer_type ${cancer_type} \\
        --wildtype_file ${wildtype_results} \\
        --mutated_file ${mutated_results} \\
        --output_dir .
    """
}
