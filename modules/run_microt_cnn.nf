process RUN_MICROT_CNN_WILDTYPE {
    tag "$cancer_type"

    container 'penny0lane/microt_cnn:v2.2'

    publishDir "${params.outdir}/${cancer_type}/microt_cnn_results",
        mode: 'copy',
        pattern: "Results_wildtype/**"

    input:
    tuple val(cancer_type), path(wildtype_fa), path(exons_wildtype), path(wildtype_mirna)

    output:
    tuple val(cancer_type),
          path("Results_wildtype/microTCNN_output_wildtype.txt"),
          emit: wildtype_results

    script:
    """
    echo "Running microT-CNN WILDTYPE for ${cancer_type}"

    # Snowfall permissions
    export HOME="\$PWD"
    mkdir -p .sfCluster Results_wildtype
    
    mkdir -p Results_wildtype
    mkdir -p .sfCluster

    cp ${params.microt_input_dir}/config.wildtype.yml .

    Rscript /R/main.R config.wildtype.yml

    echo "WILDTYPE completed"
    """
}


process RUN_MICROT_CNN_MUTATED {
    tag "$cancer_type"

    container 'penny0lane/microt_cnn:v2.2'

    publishDir "${params.outdir}/${cancer_type}/microt_cnn_results",
        mode: 'copy',
        pattern: "Results_mutated/**"

    input:
    tuple val(cancer_type), path(mutated_fa), path(exons_mutated),
          path(wildtype_mirna), path(wildtype_results)

    output:
    tuple val(cancer_type),
          path("Results_mutated/microTCNN_output_mutated.txt"),
          emit: mutated_results

    script:
    """
    echo "Running microT-CNN MUTATED for ${cancer_type}"

    # Snowfall permissions
    export HOME="\$PWD"
    mkdir -p .sfCluster Results_mutated
    
    mkdir -p Results_mutated
    mkdir -p .sfCluster

    cp ${params.microt_input_dir}/config.mutated.yml .

    Rscript /R/main.R config.mutated.yml

    echo "MUTATED completed"
    """
}

