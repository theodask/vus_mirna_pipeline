process CLINVAR_INTERSECTION {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/clinvar_annotated", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(filtered_bed), path(clinvar_vcf), path(script_file)
    
    output:
    tuple val(cancer_type), path("filtered_variants_intersected_with_clinvar_plus_clinvar_ids_${cancer_type}.bed"), emit: clinvar_bed
    
    script:
    """
    echo "Annotating variants with ClinVar for cancer type: ${cancer_type}"
    echo "Input BED: ${filtered_bed}"
    echo "ClinVar VCF: ${clinvar_vcf}"
    
    Rscript ${script_file} \\
        --cancer_type ${cancer_type} \\
        --input_bed ${filtered_bed} \\
        --clinvar_vcf ${clinvar_vcf} \\
        --output filtered_variants_intersected_with_clinvar_plus_clinvar_ids_${cancer_type}.bed
    """
}






