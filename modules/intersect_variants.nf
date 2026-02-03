process INTERSECT_VARIANTS {
    tag "$cancer_type"
    
    publishDir "${params.outdir}/${cancer_type}/intersected_variants", mode: 'copy'
    
    input:
    tuple val(cancer_type), path(vcf_file), path(seed_bed), path(script_file)
    
    output:
    tuple val(cancer_type), path("merged_${cancer_type}_binding_sites_only.vcf.gz"), emit: intersected_vcf
    
    script:
    """
    echo "Intersecting variants with seed regions for cancer type: ${cancer_type}"
    echo "VCF file: ${vcf_file}"
    echo "Seed BED: ${seed_bed}"
    
    bash ${script_file} \\
        ${vcf_file} \\
        ${seed_bed} \\
        merged_${cancer_type}_binding_sites_only.vcf.gz \\
        ${cancer_type}
    """
}
