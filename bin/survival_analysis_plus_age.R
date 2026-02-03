#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(stringr)
    library(dplyr)
    library(survminer)
    library(optparse)
})

# Define command-line arguments
option_list <- list(
    make_option(c("--cancer_type"), type="character", default=NULL,
                help="Cancer type identifier", metavar="character"),
    make_option(c("--microt_results"), type="character", default=NULL,
                help="Path to microT-CNN results file", metavar="character"),
    make_option(c("--variant_bed"), type="character", default=NULL,
                help="Path to variant BED file with patient info", metavar="character"),
    make_option(c("--clinical_tar"), type="character", default=NULL,
                help="Path to clinical tar.gz file", metavar="character"),
    make_option(c("--output_dir"), type="character", default=NULL,
                help="Output directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cancer_type) || is.null(opt$microt_results) || 
    is.null(opt$variant_bed) || is.null(opt$clinical_tar) || 
    is.null(opt$output_dir)) {
    print_help(opt_parser)
    stop("All required arguments must be provided.", call.=FALSE)
}

# Assign variables
cancer_type <- opt$cancer_type
microt_results_file <- opt$microt_results
variant_bed_file <- opt$variant_bed
clinical_tar_file <- opt$clinical_tar
output_dir <- opt$output_dir

cat("========================================\n")
cat("Survival Analysis\n")
cat("========================================\n")
cat("Cancer type:      ", cancer_type, "\n")
cat("microT results:   ", microt_results_file, "\n")
cat("Variant BED:      ", variant_bed_file, "\n")
cat("Clinical data:    ", clinical_tar_file, "\n")
cat("Output directory: ", output_dir, "\n\n")

# Create output directories
survival_plots_dir <- file.path(output_dir, "survival_plots")
dir.create(survival_plots_dir, recursive = TRUE, showWarnings = FALSE)

extract_dir <- file.path(output_dir, "extracted_clinical")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
cat("Loading data files...\n")
microt_results <- fread(microt_results_file)
variant_patient <- fread(variant_bed_file, header = TRUE)
cat("Loaded", nrow(microt_results), "microT-CNN results\n")
cat("Loaded", nrow(variant_patient), "variant-patient mappings\n\n")

# --- Extract and Process Clinical Data ---
cat("Extracting clinical data from tar.gz...\n")
untar(clinical_tar_file, files = "clinical.tsv", exdir = extract_dir)
clinical <- fread(file.path(extract_dir, "clinical.tsv"))

untar(clinical_tar_file, files = "follow_up.tsv", exdir = extract_dir)
follow_up <- fread(file.path(extract_dir, "follow_up.tsv"))

clinical <- as.data.frame(clinical)
follow_up <- as.data.frame(follow_up)

# Rename and select relevant clinical columns
cat("Processing clinical data...\n")
clinical <- clinical %>%
    dplyr::rename(
        cases_submitter_id = `cases.submitter_id`,
        vital_status = `demographic.vital_status`,
        days_to_death = `demographic.days_to_death`,
        age_at_index = `demographic.age_at_index`
    ) %>%
    dplyr::select(cases_submitter_id, vital_status, days_to_death, age_at_index)

# Convert to character first to replace "'--", then convert to numeric after cleaning
clinical$days_to_death <- as.character(clinical$days_to_death)
clinical$days_to_death[clinical$days_to_death %in% c("--", "'--")] <- NA
clinical$days_to_death <- as.numeric(clinical$days_to_death)

clinical$age_at_index <- as.character(clinical$age_at_index)
clinical$age_at_index[clinical$age_at_index %in% c("--", "'--")] <- NA
clinical$age_at_index <- as.numeric(clinical$age_at_index)

# Normalize patient IDs in clinical data
clinical$cases_submitter_id <- toupper(trimws(clinical$cases_submitter_id))
clinical$cases_submitter_id <- substr(clinical$cases_submitter_id, 1, 12)

# Process follow-up data
follow_up_clean <- follow_up %>%
    transmute(
        cases_submitter_id = toupper(substr(trimws(`cases.submitter_id`), 1, 12)),
        days_to_follow_up = na_if(as.character(`follow_ups.days_to_follow_up`), "'--")
    ) %>%
    mutate(days_to_follow_up = as.numeric(days_to_follow_up)) %>%
    group_by(cases_submitter_id) %>%
    summarise(days_to_follow_up = max(days_to_follow_up, na.rm = TRUE), .groups = "drop")

# Build event/time columns
clinical <- clinical %>%
    left_join(follow_up_clean, by = "cases_submitter_id") %>%
    mutate(
        event = ifelse(vital_status == "Dead", 1, 0),
        time = ifelse(!is.na(days_to_death), days_to_death, days_to_follow_up)
    ) %>%
    dplyr::select(cases_submitter_id, time, event, age_at_index) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::distinct(cases_submitter_id, .keep_all = TRUE)

cat(sprintf("Total unique clinical patients after cleaning: %d\n\n", nrow(clinical)))

# --- Process Variant-Patient Mappings ---
cat("Processing variant-patient mappings...\n")
variant_patient[, Patient_Hosts := lapply(Patient_Hosts, function(x) {
    patients <- unlist(strsplit(x, ",|;"))
    short_ids <- toupper(substr(trimws(patients), 1, 12))
    return(short_ids)
})]

# Map CLINVAR_ID to patient list
clinvar_to_patients <- variant_patient[, .(Patient_Hosts = unique(unlist(Patient_Hosts))), by = CLINVAR_ID]
clinvar_to_patients <- split(clinvar_to_patients$Patient_Hosts, clinvar_to_patients$CLINVAR_ID)
cat(sprintf("Mapped %d unique variants to patients\n\n", length(clinvar_to_patients)))

# --- Survival Analysis Function ---
run_survival_analysis <- function(variant, patients_with_variant, clinical_df, output_plots_dir) {
    clinical_df$mutated <- ifelse(clinical_df$cases_submitter_id %in% patients_with_variant, 1, 0)
    
    # Remove patients with missing age for this multivariate analysis
    clinical_df <- clinical_df %>% dplyr::filter(!is.na(age_at_index))
    
    group_sizes <- table(clinical_df$mutated)
    
    # Check if both groups exist and meet the minimum size requirement
    if (length(group_sizes) < 2) {
        cat(sprintf("  Skipping variant %s: Only one group found\n", variant))
        return(list(pval = NA, hr = NA, hr_lci = NA, hr_uci = NA))
    }
    
    if (any(group_sizes < 2)) {
        cat(sprintf("  Skipping variant %s: Small group sizes (WT: %d, MUT: %d)\n", 
                    variant, group_sizes["0"], group_sizes["1"]))
        return(list(pval = NA, hr = NA, hr_lci = NA, hr_uci = NA))
    }
    
    # Perform Cox Proportional Hazards Model, including age_at_index
    cox_model <- coxph(Surv(time, event) ~ mutated + age_at_index, data = clinical_df)
    cox_summary <- summary(cox_model)
    
    # Extract p-value for the 'mutated' variable
    mutated_row_index <- which(rownames(cox_summary$coefficients) == "mutated")
    if (length(mutated_row_index) == 1) {
        pval_mutated <- cox_summary$coefficients[mutated_row_index, "Pr(>|z|)"]
        hr <- cox_summary$conf.int[mutated_row_index, "exp(coef)"]
        hr_lci <- cox_summary$conf.int[mutated_row_index, "lower .95"]
        hr_uci <- cox_summary$conf.int[mutated_row_index, "upper .95"]
    } else {
        pval_mutated <- NA
        hr <- NA
        hr_lci <- NA
        hr_uci <- NA
    }
    
    # For plotting, use Kaplan-Meier and log-rank test
    surv_diff <- survdiff(Surv(time, event) ~ mutated, data = clinical_df)
    if (!is.null(surv_diff$chisq) && !is.na(surv_diff$chisq)) { 
        pval_logrank <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)
    } else {
        pval_logrank <- NA
    }
    
    fit <- survfit(Surv(time, event) ~ mutated, data = clinical_df)
    plot_file <- file.path(output_plots_dir, paste0("survival_variant_", variant, ".tiff"))
    
    tiff(plot_file, width = 2400, height = 2000, res = 300)
    ggsurv <- ggsurvplot(
        fit,
        data = clinical_df,
        conf.int = TRUE,
        pval = paste0("Log-rank P = ", signif(pval_logrank, 3)),
        legend.labs = c("Wildtype", "Mutated"),
        legend.title = "Group",
        palette = c("blue", "red"),
        xlab = "Days",
        ylab = "Survival Probability",
        title = paste("Survival for variant", variant),
        subtitle = paste0(
            "Cox Model (Adj. for Age): Mutated HR = ", signif(hr, 3), 
            " (95% CI: ", signif(hr_lci, 3), "-", signif(hr_uci, 3), 
            "), P = ", signif(pval_mutated, 3)
        ),
        risk.table = TRUE,
        risk.table.col = "strata",
        risk.table.height = 0.25
    )
    print(ggsurv)
    dev.off()
    
    return(list(pval = pval_mutated, hr = hr, hr_lci = hr_lci, hr_uci = hr_uci))
}

# --- Run Survival Analysis for All Variants ---
cat("Running survival analysis for all variants...\n")
survival_results <- data.frame(
    variant_id = character(), 
    pvalue_mutated_cox = numeric(),
    hr_mutated = numeric(),
    hr_lci = numeric(),
    hr_uci = numeric(),
    stringsAsFactors = FALSE
)

# Extract unique variants
unique_variants <- unique(unlist(strsplit(microt_results$variant_id, ";")))
cat(sprintf("Total of %d unique variants disrupting MREs were found.\n\n", length(unique_variants)))

# Loop over variants
for (variant in unique_variants) {
    cat(sprintf(">>> Variant: %s\n", variant))
    
    if (variant %in% names(clinvar_to_patients)) {
        patients_with_variant <- unlist(clinvar_to_patients[[variant]])
        patients_with_variant <- toupper(trimws(patients_with_variant))
        patients_with_variant <- substr(patients_with_variant, 1, 12)
        
        intersected <- intersect(patients_with_variant, clinical$cases_submitter_id)
        cat(sprintf("  Mutated patients: %d in BED, %d in clinical data\n", 
                    length(patients_with_variant), length(intersected)))
        
        results_list <- tryCatch({
            run_survival_analysis(variant, intersected, clinical, survival_plots_dir)
        }, error = function(e) {
            message(sprintf("  Error for variant %s: %s", variant, e$message))
            return(NULL)
        })
        
        if (!is.null(results_list)) {
            survival_results <- rbind(survival_results, data.frame(
                variant_id = variant, 
                pvalue_mutated_cox = results_list$pval, 
                hr_mutated = results_list$hr,
                hr_lci = results_list$hr_lci,
                hr_uci = results_list$hr_uci,
                stringsAsFactors = FALSE
            ))
        }
    } else {
        cat(sprintf("  Variant %s not found in clinvar_to_patients\n", variant))
    }
}

# Sort and save results
survival_results_sorted <- survival_results %>%
    arrange(is.na(pvalue_mutated_cox), pvalue_mutated_cox)

fwrite(survival_results_sorted, 
       file.path(survival_plots_dir, "survival_pvalues.tsv"), 
       sep = "\t")

# Clean up extracted directory
unlink(extract_dir, recursive = TRUE)

cat("\n========================================\n")
cat("SUCCESS!\n")
cat("Survival analysis complete for:", cancer_type, "\n")
cat("Results saved to:", survival_plots_dir, "\n")
cat("Total variants analyzed:", nrow(survival_results_sorted), "\n")
cat("========================================\n")

