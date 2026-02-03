process SURVIVAL_ANALYSIS {
    tag "${cancer_type}"
    publishDir "${params.outdir}/${cancer_type}/survival_analysis", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(microt_results), path(variant_bed), path(clinical_tar), path(survival_script), path(lasso_script)
    
    output:
    // Captures the first R script outputs
    tuple val(cancer_type), path("survival_plots/survival_pvalues.tsv"), emit: survival_pvalues, optional: true
    tuple val(cancer_type), path("survival_plots/*.tiff"), emit: survival_plots, optional: true
    
    // Capture ALL TSV and PNG files inside the lasso_analysis folder
    tuple val(cancer_type), path("lasso_analysis/*.tsv"), emit: lasso_tables
    tuple val(cancer_type), path("lasso_analysis/*.png"), emit: lasso_plots, optional: true

    
    script:
    """
    # Run survival analysis with age adjustment
    echo "=== Running Survival Analysis ==="
    Rscript ${survival_script} \
        --cancer_type ${cancer_type} \
        --microt_results ${microt_results} \
        --variant_bed ${variant_bed} \
        --clinical_tar ${clinical_tar} \
        --output_dir .
    
    # Run LASSO Cox regression
    echo ""
    echo "=== Running LASSO Cox Regression ==="
    Rscript ${lasso_script} \
        --cancer_type ${cancer_type} \
        --microt_results ${microt_results} \
        --variant_bed ${variant_bed} \
        --clinical_tar ${clinical_tar} \
        --output_dir ./lasso_analysis
    
    echo ""
    echo "=== All Survival Analyses Complete ==="
    """
}


