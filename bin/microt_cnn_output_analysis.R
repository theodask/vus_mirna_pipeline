#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("--cancer_type"), type="character", default=NULL,
              help="Cancer type identifier", metavar="character"),
  make_option(c("--wildtype_file"), type="character", default=NULL,
              help="Path to wildtype MRE file", metavar="character"),
  make_option(c("--mutated_file"), type="character", default=NULL,
              help="Path to mutated MRE file", metavar="character"),
  make_option(c("--output_dir"), type="character", default=NULL,
              help="Output directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cancer_type) || is.null(opt$wildtype_file) || 
    is.null(opt$mutated_file) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("All required arguments (--cancer_type, --wildtype_file, --mutated_file, --output_dir) must be provided.", call.=FALSE)
}

cancer_type   <- opt$cancer_type
wildtype_file <- opt$wildtype_file
mutated_file  <- opt$mutated_file
output_dir    <- opt$output_dir

cat("microT-CNN Output Analysis\n")
cat("Cancer type:     ", cancer_type, "\n")
cat("Wildtype file:   ", wildtype_file, "\n")
cat("Mutated file:    ", mutated_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Check if input files exist
if (!file.exists(wildtype_file)) {
  stop("Wildtype file not found: ", wildtype_file, call.=FALSE)
}
if (!file.exists(mutated_file)) {
  stop("Mutated file not found: ", mutated_file, call.=FALSE)
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

cat("Reading input files...\n")
wildtype <- fread(wildtype_file)
mutated  <- fread(mutated_file)

cat("Read", nrow(wildtype), "wildtype MREs\n")
cat("Read", nrow(mutated), "mutated MREs\n\n")

# ID Extraction
cat("Extracting transcript IDs and variant IDs...\n")
wildtype <- wildtype %>%
  mutate(transcript_id = str_extract(trRegion_name, "ENST\\d+"))

mutated <- mutated %>%
  mutate(
    transcript_id = str_extract(trRegion_name, "ENST\\d+"),
    variant_id = str_extract_all(trRegion_name, "_\\d+") %>%
      lapply(function(x) gsub("_", "", x)) %>%
      sapply(function(x) paste(x, collapse = ";"))
  )

# Clean trRegion_name in mutated set
mutated$trRegion_name <- str_remove_all(mutated$trRegion_name, "_\\d+")

# Unique MRE Identification
cat("Identifying unique MREs per scenario...\n")

# Only Mutated
only.mutated.df <- mutated[!(mutated$mre.rel.tr.start %in% wildtype$mre.rel.tr.start &
                               mutated$mre.rel.tr.end %in% wildtype$mre.rel.tr.end),]

# Only Wildtype
only.wildtype.df <- wildtype[!(wildtype$mre.rel.tr.start %in% mutated$mre.rel.tr.start &
                                 wildtype$mre.rel.tr.end %in% mutated$mre.rel.tr.end),]

fwrite(only.mutated.df, 
       file.path(output_dir, paste0("only_mutated_MREs_", cancer_type, ".tsv")),
       sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

fwrite(only.wildtype.df, 
       file.path(output_dir, paste0("only_wildtype_MREs_", cancer_type, ".tsv")),
       sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Keep unique mutated rows
only.mutated.df <- unique(only.mutated.df)

# Disruption Score Analysis 
cat("Finding disrupted MREs (common to both)...\n")
disrupted.MREs <- merge(mutated, wildtype,
                        by = c("trRegion_name", "mirnaname", "chr",
                               "collapsed.mre.abs.coordinates", "strand",
                               "mre.rel.tr.start", "mre.rel.tr.end",
                               "region", "region_training"),
                        suffixes = c(".mut", ".wt"))

cat("Computing disruption scores for", nrow(disrupted.MREs), "MREs...\n")

# Compute disruption score and process columns
disrupted.MREs <- disrupted.MREs %>%
  mutate(
    disruption_score = ((MRE_score.mut - MRE_score.wt) * log2(1 + ((MRE_score.mut + MRE_score.wt) / 2))) / log2(1.5),
    Gene = str_extract(trRegion_name, "ENSG\\d+"),
    binding_coordinates = collapsed.mre.abs.coordinates,
    binding_start = as.numeric(gsub(":", "", str_extract(collapsed.mre.abs.coordinates, ":\\d+")))
  )

# Rename columns and select
disrupted.MREs <- as.data.frame(disrupted.MREs) %>%
  dplyr::rename(
    miRNA = mirnaname,
    interaction_score_mut = MRE_score.mut,
    interaction_score_wt = MRE_score.wt,
    binding_type_mut = binding_type.mut,
    binding_type_wt = binding_type.wt,
    Strand = strand
  ) %>%
  dplyr::select(miRNA, Gene, binding_coordinates, Strand, variant_id,
                interaction_score_mut, interaction_score_wt, disruption_score,
                binding_type_mut, binding_type_wt)

# Sort by absolute disruption score
disrupted.MREs <- disrupted.MREs[order(abs(disrupted.MREs$disruption_score), decreasing = TRUE), ]

fwrite(disrupted.MREs, 
       file.path(output_dir, paste0("microT_CNN_wt_mut_results_", cancer_type, ".tsv")),
       sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Visualization
cat("\nGenerating histograms...\n")

# Disruption Scores
cat("- Plotting disruption scores...\n")
histogram_data <- disrupted.MREs$disruption_score
breaks_disruption <- seq(floor(min(histogram_data)), ceiling(max(histogram_data)), by = 0.1)

tiff(file.path(output_dir, "histogram_MREs_disruption_scores.tiff"),
     res = 300, compression = "lzw", units = "px", width = 2400, height = 2000)
hist(histogram_data, freq = TRUE,
     xlab = 'Disruption score',
     ylab = 'Number of predicted MREs',
     main = 'Distribution of MRE disruption scores',
     labels = TRUE,
     col = "#41b3a3",
     breaks = breaks_disruption)
dev.off()

# Gained MREs
cat("- Plotting gained MREs...\n")
gained_scores <- only.mutated.df$MRE_score
breaks_gained <- seq(floor(min(gained_scores)), ceiling(max(gained_scores)), by = 0.1)

tiff(file.path(output_dir, "histogram_gained_MREs_scores.tiff"),
     res = 300, compression = "lzw", units = "px", width = 2400, height = 2000)
hist(gained_scores, freq = TRUE,
     xlab = 'MRE score',
     ylab = 'Number of gained MREs',
     main = 'Distribution of gained MRE scores',
     labels = TRUE,
     col = "#1f78b4",
     breaks = breaks_gained)
dev.off()

# Abolished MREs
cat("- Plotting abolished MREs...\n")
lost_scores <- only.wildtype.df$MRE_score
breaks_lost <- seq(floor(min(lost_scores)), ceiling(max(lost_scores)), by = 0.1)

tiff(file.path(output_dir, "histogram_lost_MREs_scores.tiff"),
     res = 300, compression = "lzw", units = "px", width = 2400, height = 2000)
hist(lost_scores, freq = TRUE,
     xlab = 'MRE score',
     ylab = 'Number of abolished MREs',
     main = 'Distribution of abolished MRE scores',
     labels = TRUE,
     col = "#e6550d",
     breaks = breaks_lost)
dev.off()

cat("Complete. Files saved to:", output_dir, "\n")

