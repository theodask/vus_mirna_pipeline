process PREPARE_MICROT_INPUTS {
    tag "$cancer_type"
    
    // Publish all outputs to the microT-CNN shared directory
    publishDir "${params.microt_input_dir}", mode: 'copy', pattern: "*.{fa,tab,fasta,txt,bed}"
    
    input:
    tuple val(cancer_type), path(microt_input_bed), path(mirbase_gff), path(ensembl_gtf), path(mre_interactions), path(script_file)
    
    output:
    tuple val(cancer_type), path("wildtype.fa"), emit: wildtype_fa, optional: true
    tuple val(cancer_type), path("mutated.fa"), emit: mutated_fa, optional: true
    tuple val(cancer_type), path("wildtype_mirna_sequences.fasta"), emit: wildtype_mirna, optional: true
    tuple val(cancer_type), path("exons.wildtype.tab"), emit: exons_wildtype, optional: true
    tuple val(cancer_type), path("exons.mutated.tab"), emit: exons_mutated, optional: true
    tuple val(cancer_type), path("transcript_id_list.txt"), emit: transcript_list, optional: true
    tuple val(cancer_type), path("ssm_intersection.bed"), emit: ssm_intersection, optional: true
    
    script:
    """
    echo "Preparing microT-CNN Inputs"
    echo "Cancer type      : ${cancer_type}"
    echo "Input BED        : ${microt_input_bed}"
    echo "miRBase GFF3     : ${mirbase_gff}"
    echo "Ensembl GTF      : ${ensembl_gtf}"
    echo "MRE interactions : ${mre_interactions}"
    echo ""
    
    # Run the R script
    Rscript ${script_file} \\
        --cancer_type ${cancer_type} \\
        --input ${microt_input_bed} \\
        --mirbase_input ${mirbase_gff} \\
        --ensembl_input ${ensembl_gtf} \\
        --analysis mre \\
        --mre_interactions ${mre_interactions}
    
    echo ""
    echo "microT-CNN input preparation completed"
    """
}





