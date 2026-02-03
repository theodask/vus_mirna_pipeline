process ANNOTATE_VARIANTS_MRES {
    tag "$cancer_type"
    
    // Multiple publishDir directives for different outputs
    publishDir "${params.outdir}/${cancer_type}/mre_annotated", mode: 'copy', pattern: "variants_in_seed_${cancer_type}.bed"
    publishDir "${params.microt_input_dir}", mode: 'copy', pattern: "microt_input_${cancer_type}.bed"
    
    input:
    tuple val(cancer_type), path(clinvar_bed), path(mres_file), path(script_file)
    
    output:
    tuple val(cancer_type), path("variants_in_seed_${cancer_type}.bed"), emit: annotated_bed
    tuple val(cancer_type), path("microt_input_${cancer_type}.bed"), emit: microt_input
    
    script:
    """
    echo "Annotating variants with MREs for cancer type: ${cancer_type}"
    echo "Input BED: ${clinvar_bed}"
    echo "MREs file: ${mres_file}"
    
    Rscript ${script_file} \\
        --cancer_type ${cancer_type} \\
        --input_bed ${clinvar_bed} \\
        --mres_file ${mres_file} \\
        --output_annotated variants_in_seed_${cancer_type}.bed \\
        --output_microt microt_input_${cancer_type}.bed
    """
}
