#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { FIND_GENOMIC_REGIONS } from './modules/find_genomic_regions'
include { SELECT_MRES } from './modules/select_mres'
include { CALCULATE_SEED_REGIONS } from './modules/calculate_seed_regions'
include { INTERSECT_VARIANTS } from './modules/intersect_variants'
include { PROCESS_VCF_CHUNKS } from './modules/process_vcf_chunks'
include { CLINVAR_INTERSECTION } from './modules/clinvar_intersection'
include { ANNOTATE_VARIANTS_MRES } from './modules/annotate_variants_mres'
include { PREPARE_MICROT_INPUTS } from './modules/prepare_microt_inputs'
include { RUN_MICROT_CNN_WILDTYPE } from './modules/run_microt_cnn'
include { RUN_MICROT_CNN_MUTATED } from './modules/run_microt_cnn'
include { ANALYZE_MICROT_RESULTS } from './modules/analyze_microt_results'
include { SURVIVAL_ANALYSIS } from './modules/survival_analysis'

// Print pipeline header
log.info """\
    VUS-miRNA ANALYSIS PIPELINE
    Cancer type    : ${params.cancer_type ?: 'Not specified'}
    VCF file       : ${params.vcf ?: 'Not specified'}
    GTF file       : ${params.gtf}
    FAI file       : ${params.fai}
    MREs file      : ${params.mres_file}
    miTED dir      : ${params.mited_dir}
    ClinVar VCF    : ${params.clinvar_vcf}
    microT dir     : ${params.microt_input_dir}
    miRBase GFF3   : ${params.mirbase_gff}
    Ensembl GTF    : ${params.ensembl_gtf_full}
    RPM threshold  : ${params.rpm_threshold}
    Lines/chunk    : ${params.lines_per_chunk}
    Output dir     : ${params.outdir}
    Skip genomic   : ${params.skip_genomic_regions ?: false}
    """
    .stripIndent()

workflow {
    // Create reference file channels (these are the same for all samples)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    fai_ch = Channel.fromPath(params.fai, checkIfExists: true)
    genomic_regions_script_ch = Channel.fromPath('bin/find_genomic_regions.py', checkIfExists: true)
    select_mres_script_ch = Channel.fromPath('bin/select_MREs_from_top_mirnas_expressed_in_mited.R', checkIfExists: true)
    seed_regions_script_ch = Channel.fromPath('bin/calculate_seed_regions_from_mited_cords.sh', checkIfExists: true)
    intersect_variants_script_ch = Channel.fromPath('bin/intersect_variants_with_MREs.sh', checkIfExists: true)
    vcf_chunks_bash_script_ch = Channel.fromPath('bin/handle_vcf_chunks.sh', checkIfExists: true)
    vcf_chunks_r_script_ch = Channel.fromPath('bin/process_and_filter_vcf.R', checkIfExists: true)
    clinvar_script_ch = Channel.fromPath('bin/clinvar_intersection.R', checkIfExists: true)
    annotate_mres_script_ch = Channel.fromPath('bin/annotate_variants_with_MREs.R', checkIfExists: true)
    prepare_microt_script_ch = Channel.fromPath('bin/prepare_microt_cnn_inputs.R', checkIfExists: true)
    clinvar_vcf_ch = Channel.fromPath(params.clinvar_vcf, checkIfExists: true)
    mirbase_gff_ch = Channel.fromPath(params.mirbase_gff, checkIfExists: true)
    ensembl_gtf_full_ch = Channel.fromPath(params.ensembl_gtf_full, checkIfExists: true)
    mres_file_ch = Channel.fromPath(params.mres_file, checkIfExists: true)
    analyze_microt_script_ch = Channel.fromPath('bin/microt_cnn_output_analysis.R', checkIfExists: true)

    // Option 1: Single cancer type
    if (params.cancer_type && params.vcf) {
        vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
        cancer_type_ch = Channel.value(params.cancer_type)
        
        // Store VCF info for later use in intersection
        vcf_for_intersect_ch = Channel.value(params.cancer_type).combine(vcf_ch)
        
        // Build miTED file path
        mited_file = file("${params.mited_dir}/${params.cancer_type}/miTED_top_expressed_miRNA_normal_and_cancer_RPM_${params.cancer_type}.tsv")
        if (!mited_file.exists()) {
            error "miTED file not found: ${mited_file}"
        }
        mited_ch = Channel.fromPath(mited_file)
        
        // Combine inputs for FIND_GENOMIC_REGIONS
        genomic_input_ch = genomic_regions_script_ch
            .combine(cancer_type_ch)
            .combine(vcf_ch)
            .combine(gtf_ch)
            .combine(fai_ch)
        
        // Combine inputs for SELECT_MRES
        mres_input_ch = select_mres_script_ch
            .combine(cancer_type_ch)
            .combine(mited_ch)
            .combine(mres_file_ch)
    }
    // Option 2: Multiple cancer types from CSV
    else if (params.sample_sheet) {
        def sample_data = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map { row -> 
                def mited_file = file("${params.mited_dir}/${row.cancer_type}/miTED_top_expressed_miRNA_normal_and_cancer_RPM_${row.cancer_type}.tsv")
                if (!mited_file.exists()) {
                    error "miTED file not found for ${row.cancer_type}: ${mited_file}"
                }
                tuple(row.cancer_type, file(row.vcf_path, checkIfExists: true), mited_file)
            }
        
        // Separate channels for each process
        cancer_vcf_ch = sample_data.map { cancer_type, vcf, mited -> tuple(cancer_type, vcf) }
        cancer_mited_ch = sample_data.map { cancer_type, vcf, mited -> tuple(cancer_type, mited) }
        
        // Store VCF info for intersection
        vcf_for_intersect_ch = cancer_vcf_ch
        
        // Combine for FIND_GENOMIC_REGIONS
        genomic_input_ch = genomic_regions_script_ch
            .combine(cancer_vcf_ch)
            .combine(gtf_ch)
            .combine(fai_ch)
        
        // Combine for SELECT_MRES
        mres_input_ch = select_mres_script_ch
            .combine(cancer_mited_ch)
            .combine(mres_file_ch)
    }
    // Option 3: Auto-discover VCF files
    else if (params.vcf_dir) {
        def sample_data = Channel
            .fromPath("${params.vcf_dir}/**/merged_*.vcf.gz")
            .map { vcf_file -> 
                def cancer_type = vcf_file.name.replaceAll(/merged_(.+)\.vcf\.gz/, '$1')
                def mited_file = file("${params.mited_dir}/${cancer_type}/miTED_top_expressed_miRNA_normal_and_cancer_RPM_${cancer_type}.tsv")
                if (!mited_file.exists()) {
                    log.warn "miTED file not found for ${cancer_type}, skipping: ${mited_file}"
                    return null
                }
                tuple(cancer_type, vcf_file, mited_file)
            }
            .filter { it != null }
        
        // Separate channels for each process
        cancer_vcf_ch = sample_data.map { cancer_type, vcf, mited -> tuple(cancer_type, vcf) }
        cancer_mited_ch = sample_data.map { cancer_type, vcf, mited -> tuple(cancer_type, mited) }
        
        // Store VCF info for intersection
        vcf_for_intersect_ch = cancer_vcf_ch
        
        // Combine for FIND_GENOMIC_REGIONS
        genomic_input_ch = genomic_regions_script_ch
            .combine(cancer_vcf_ch)
            .combine(gtf_ch)
            .combine(fai_ch)
        
        // Combine for SELECT_MRES
        mres_input_ch = select_mres_script_ch
            .combine(cancer_mited_ch)
            .combine(mres_file_ch)
    }
    else {
        error "Please provide either --cancer_type and --vcf, or --sample_sheet, or --vcf_dir"
    }
    
    // --- Conditional Execution of FIND_GENOMIC_REGIONS ---
    if (!params.skip_genomic_regions) {
        log.info "--- Executing FIND_GENOMIC_REGIONS process. This may take time. ---"
        // Run the process
        FIND_GENOMIC_REGIONS(genomic_input_ch)
        
        // Emit outputs for logging/downstream visibility
        FIND_GENOMIC_REGIONS.out.regions.view { cancer_type, file -> 
            "Generated regions file for ${cancer_type}: ${file}" 
        }
        FIND_GENOMIC_REGIONS.out.summary.view { cancer_type, file -> 
            "Generated summary file for ${cancer_type}: ${file}" 
        }
        FIND_GENOMIC_REGIONS.out.plot.view { cancer_type, file -> 
            "Generated plot for ${cancer_type}: ${file}" 
        }
    } else {
        log.info "--- FIND_GENOMIC_REGIONS process skipped as requested (--skip_genomic_regions true). ---"
    }
    
    // --- Execute SELECT_MRES process ---
    log.info "--- Executing SELECT_MRES process ---"
    SELECT_MRES(mres_input_ch)
    
    // Emit outputs for logging
    SELECT_MRES.out.mre_bed.view { cancer_type, file ->
        "Generated MRE coordinates file for ${cancer_type}: ${file}"
    }
    
    // --- Execute CALCULATE_SEED_REGIONS process ---
    log.info "--- Executing CALCULATE_SEED_REGIONS process ---"
    
    // Combine SELECT_MRES output with the script
    seed_input_ch = SELECT_MRES.out.mre_bed
        .combine(seed_regions_script_ch)
    
    CALCULATE_SEED_REGIONS(seed_input_ch)
    
    // Emit outputs for logging
    CALCULATE_SEED_REGIONS.out.seed_bed.view { cancer_type, file ->
        "Generated seed region coordinates file for ${cancer_type}: ${file}"
    }
    
    // --- Execute INTERSECT_VARIANTS process ---
    log.info "--- Executing INTERSECT_VARIANTS process ---"
    
    // Combine VCF, seed regions, and script for intersection
    intersect_input_ch = vcf_for_intersect_ch
        .join(CALCULATE_SEED_REGIONS.out.seed_bed)
        .combine(intersect_variants_script_ch)
    
    INTERSECT_VARIANTS(intersect_input_ch)
    
    // Emit outputs for logging
    INTERSECT_VARIANTS.out.intersected_vcf.view { cancer_type, file ->
        "Generated intersected VCF for ${cancer_type}: ${file}"
    }
    
    // --- Execute PROCESS_VCF_CHUNKS process ---
    log.info "--- Executing PROCESS_VCF_CHUNKS process ---"
    
    // Combine intersected VCF with both scripts
    vcf_chunks_input_ch = INTERSECT_VARIANTS.out.intersected_vcf
        .combine(vcf_chunks_bash_script_ch)
        .combine(vcf_chunks_r_script_ch)
    
    PROCESS_VCF_CHUNKS(vcf_chunks_input_ch)
    
    // Emit outputs for logging
    PROCESS_VCF_CHUNKS.out.filtered_bed.view { cancer_type, file ->
        "Generated filtered variants BED for ${cancer_type}: ${file}"
    }
    
    // --- Execute CLINVAR_INTERSECTION process ---
    log.info "--- Executing CLINVAR_INTERSECTION process ---"
    
    // Combine filtered BED with ClinVar VCF and script
    clinvar_input_ch = PROCESS_VCF_CHUNKS.out.filtered_bed
        .combine(clinvar_vcf_ch)
        .combine(clinvar_script_ch)
    
    CLINVAR_INTERSECTION(clinvar_input_ch)
    
    // Emit outputs for logging
    CLINVAR_INTERSECTION.out.clinvar_bed.view { cancer_type, file ->
        "Generated ClinVar-annotated variants for ${cancer_type}: ${file}"
    }
    
    // --- Execute ANNOTATE_VARIANTS_MRES process ---
    log.info "--- Executing ANNOTATE_VARIANTS_MRES process ---"
    
    // Combine ClinVar-annotated BED with MREs file and script
    annotate_mres_input_ch = CLINVAR_INTERSECTION.out.clinvar_bed
        .combine(mres_file_ch)
        .combine(annotate_mres_script_ch)
    
    ANNOTATE_VARIANTS_MRES(annotate_mres_input_ch)
    
    // Emit outputs for logging
    ANNOTATE_VARIANTS_MRES.out.annotated_bed.view { cancer_type, file ->
        "Generated MRE-annotated variants for ${cancer_type}: ${file}"
    }
    ANNOTATE_VARIANTS_MRES.out.microt_input.view { cancer_type, file ->
        "Generated microT-CNN input for ${cancer_type}: ${file}"
    }
    
    // --- Execute PREPARE_MICROT_INPUTS process ---
    log.info "--- Executing PREPARE_MICROT_INPUTS process ---"
    
    // Combine microT input BED with resource files and script
    prepare_microt_input_ch = ANNOTATE_VARIANTS_MRES.out.microt_input
        .combine(mirbase_gff_ch)
        .combine(ensembl_gtf_full_ch)
        .combine(mres_file_ch)
        .combine(prepare_microt_script_ch)
    
    PREPARE_MICROT_INPUTS(prepare_microt_input_ch)
    
    // Emit outputs for logging
    PREPARE_MICROT_INPUTS.out.wildtype_fa.view { cancer_type, file ->
        "Generated wildtype sequences for ${cancer_type}: ${file}"
    }
    PREPARE_MICROT_INPUTS.out.mutated_fa.view { cancer_type, file ->
        "Generated mutated sequences for ${cancer_type}: ${file}"
    }
    
    // --- Execute RUN_MICROT_CNN_WILDTYPE process ---
    log.info "--- Executing microT-CNN WILDTYPE analysis ---"
    
    // Combine wildtype inputs for microT-CNN
    microt_wildtype_input_ch = PREPARE_MICROT_INPUTS.out.wildtype_fa
        .join(PREPARE_MICROT_INPUTS.out.exons_wildtype)
        .join(PREPARE_MICROT_INPUTS.out.wildtype_mirna)
    
    RUN_MICROT_CNN_WILDTYPE(microt_wildtype_input_ch)
    
    RUN_MICROT_CNN_WILDTYPE.out.wildtype_results.view { cancer_type, file ->
        "microT-CNN WILDTYPE analysis completed for ${cancer_type}"
    }
    
    // --- Execute RUN_MICROT_CNN_MUTATED process ---
    log.info "--- Executing microT-CNN MUTATED analysis ---"
    
    // Combine mutated inputs with wildtype results
    microt_mutated_input_ch = PREPARE_MICROT_INPUTS.out.mutated_fa
        .join(PREPARE_MICROT_INPUTS.out.exons_mutated)
        .join(PREPARE_MICROT_INPUTS.out.wildtype_mirna)
        .join(RUN_MICROT_CNN_WILDTYPE.out.wildtype_results)
    
    RUN_MICROT_CNN_MUTATED(microt_mutated_input_ch)
    
    RUN_MICROT_CNN_MUTATED.out.mutated_results.view { cancer_type, file ->
        "microT-CNN MUTATED analysis completed for ${cancer_type}"
    }
    
    // --- Execute ANALYZE_MICROT_RESULTS process ---
    log.info "--- Executing microT-CNN results analysis ---"
    
    // Combine wildtype and mutated results with analysis script
    analyze_microt_input_ch = RUN_MICROT_CNN_WILDTYPE.out.wildtype_results
        .join(RUN_MICROT_CNN_MUTATED.out.mutated_results)
        .combine(analyze_microt_script_ch)
    
    ANALYZE_MICROT_RESULTS(analyze_microt_input_ch)
    
    // Emit outputs for logging
    ANALYZE_MICROT_RESULTS.out.only_mutated.view { cancer_type, file ->
        "Generated gained MREs file for ${cancer_type}: ${file}"
    }
    ANALYZE_MICROT_RESULTS.out.only_wildtype.view { cancer_type, file ->
        "Generated lost MREs file for ${cancer_type}: ${file}"
    }
    ANALYZE_MICROT_RESULTS.out.disrupted_mres.view { cancer_type, file ->
        "Generated disrupted MREs analysis for ${cancer_type}: ${file}"
    }
    ANALYZE_MICROT_RESULTS.out.hist_disruption.view { cancer_type, file ->
        "Generated disruption scores histogram for ${cancer_type}: ${file}"
    }
    ANALYZE_MICROT_RESULTS.out.hist_gained.view { cancer_type, file ->
        "Generated gained MREs histogram for ${cancer_type}: ${file}"
    }
    ANALYZE_MICROT_RESULTS.out.hist_lost.view { cancer_type, file ->
        "Generated lost MREs histogram for ${cancer_type}: ${file}"
    }
    
    // --- Execute SURVIVAL_ANALYSIS process ---
    log.info "--- Executing survival analysis ---"
    
    // Create script channels
    survival_script_ch = Channel.fromPath('bin/survival_analysis_plus_age.R', checkIfExists: true)
    lasso_script_ch = Channel.fromPath('bin/lasso_cox_regression_plus_age.R', checkIfExists: true)
    
    // Handle clinical data based on input mode
    if (params.cancer_type && params.vcf) {
        // Single cancer type mode
        clinical_tar_file = file("${params.resources_dir}/metadata/${params.cancer_type}/clinical.cart.${params.cancer_type}.tar.gz")
        if (!clinical_tar_file.exists()) {
            error "Clinical data not found: ${clinical_tar_file}"
        }
        clinical_tar_ch = Channel.fromPath(clinical_tar_file)
        
        // Combine inputs for survival analysis
        survival_input_ch = ANALYZE_MICROT_RESULTS.out.disrupted_mres
            .join(ANNOTATE_VARIANTS_MRES.out.annotated_bed)
            .combine(clinical_tar_ch)
            .combine(survival_script_ch)
            .combine(lasso_script_ch)
    } else if (params.sample_sheet || params.vcf_dir) {
        // Multiple cancer types mode - need to map each cancer type to its clinical file
        clinical_tar_ch = ANALYZE_MICROT_RESULTS.out.disrupted_mres
            .map { cancer_type, file ->
                def clinical_file = file("${params.resources_dir}/metadata/${cancer_type}/clinical.cart.${cancer_type}.tar.gz")
                if (!clinical_file.exists()) {
                    log.warn "Clinical data not found for ${cancer_type}: ${clinical_file}"
                    return null
                }
                tuple(cancer_type, clinical_file)
            }
            .filter { it != null }
        
        // Combine inputs for survival analysis
        survival_input_ch = ANALYZE_MICROT_RESULTS.out.disrupted_mres
            .join(ANNOTATE_VARIANTS_MRES.out.annotated_bed)
            .join(clinical_tar_ch)
            .combine(survival_script_ch)
            .combine(lasso_script_ch)
    }
    
    SURVIVAL_ANALYSIS(survival_input_ch)
    
    // Emit outputs for logging
    SURVIVAL_ANALYSIS.out.survival_pvalues.view { cancer_type, file ->
        "Generated survival p-values for ${cancer_type}: ${file}"
    }
    
    SURVIVAL_ANALYSIS.out.lasso_tables.view { cancer_type, file ->
        "Generated LASSO Cox results for ${cancer_type}: ${file}"
    }
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Output dir: ${params.outdir}
        """
        .stripIndent()
}

