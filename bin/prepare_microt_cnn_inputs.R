suppressPackageStartupMessages({
  library(dbplyr)
  library(data.table)
  library(dplyr)
  library(optparse)
  library(seqinr)
  library(stringr)

  library(biomaRt)
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(BiocFileCache)
})

# Create cache directory in current working directory
cache_dir <- file.path(getwd(), ".BiocFileCache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Set cache location environment variables
Sys.setenv(BIOCFILECACHE_CACHE = cache_dir)
Sys.setenv(BIOCFILECACHE_ASK = "FALSE")
Sys.setenv(R_USER_CACHE_DIR = cache_dir)

cat("BiocFileCache configured to:", cache_dir, "\n")

##### Read input arguments #####
option_list = list(
  make_option(c("--cancer_type"), type="character", default=NULL,
              help="Required! Cancer type identifier", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Required! File containing Simple Somatic Mutations in BED format. (e.g. chr6 100 100 *  +  G  A  sample1)", metavar="character"),
  make_option(c("-m", "--mirbase_input"), type="character", default="hsa_miRBase_22v1_0based.gff3",
              help="Required! The miRBase gff3 file.", metavar="character"),
  make_option(c("-e", "--ensembl_input"), type="character", default="Homo_sapiens.GRCh38.113.gtf",
              help="Required! The Ensembl gtf file.", metavar="character"),
  make_option("--analysis", type="character", default="mre",
              help="Required! The analysis option between 'splicing', 'orf' and 'mre' that determines the pipeline workflow accordingly", metavar="character"),
  make_option("--orf_min_length", type="integer", default=1000,
              help="Used in 'orf' analysis. Threshold input for NCBI ORF finder, determines the minimum required length for ORFs to be found. [default= %default]", metavar="number"),
  make_option("--splicing_probability_threshold", type = "double", default = 0.5,
              help="Used in 'splicing' analysis. Threshold input for SpliceAI, determines the probability threshold for a mutation variant to be determined as splice-altering in either donor or acceptor role. [default= %default]", metavar="number"),
  make_option("--splicing_annotation", type="character", default=NULL,
              help="Required in 'splicing' analysis. The Gencode annotation in the specific format required by SpliceAI. Find more info in the SpliceAI and SpliceAI-lookup Github pages.", metavar="character"),
  make_option("--splicing_genome", type="character", default=NULL,
              help="Required in 'splicing' analysis. One of the GRCh38 or GRCh37 genomes provided by the SpliceAI Github page.", metavar="character"),
  make_option("--mre_interactions", type="character", default="TarBasev9_with_microT_binding_sites_plus_mre_type_and_cds_utr.tsv",
              help="Required in 'mre' analysis. Homo sapiens TarBase interactions provided by https://dianalab.e-ce.uth.gr/tarbasev9/downloads.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$cancer_type)){
  print_help(opt_parser)
  stop("--cancer_type argument must be supplied.", call.=FALSE)
}

# Set cancer_type variable
cancer_type <- opt$cancer_type

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("--input argument must be supplied.", call.=FALSE)
} else if (!file.exists(opt$input)) {
  stop(paste0("\"", opt$input, "\" does not exist."), call.=FALSE)
} else if (is.null(opt$mirbase_input)){
  print_help(opt_parser)
  stop("--mirbase_input argument must be supplied.", call.=FALSE)
} else if (!file.exists(opt$mirbase_input)) {
  stop(paste0("\"", opt$mirbase_input, "\" does not exist."), call.=FALSE)
} else if (is.null(opt$ensembl_input)){
  print_help(opt_parser)
  stop("--ensembl_input argument must be supplied.", call.=FALSE)
} else if (!file.exists(opt$ensembl_input)) {
  stop(paste0("\"", opt$ensembl_input, "\" does not exist."), call.=FALSE)
} else if (is.null(opt$analysis)){
  print_help(opt_parser)
  stop("--analysis argument must be supplied.", call.=FALSE)
} else if (!(opt$analysis %in% c('splicing','orf','mre'))){
  print_help(opt_parser)
  stop("--analysis argument must be specified as 'splicing', 'orf' or 'mre'.", call.=FALSE)
}

##### Define functions #####

bed_intersect_data_frames <- function(df1, df2) {
  
  temp_dir <- tempdir()
  file1 <- file.path(temp_dir, "file1.bed")
  file2 <- file.path(temp_dir, "file2.bed")
  fwrite(df1, file1, col.names = FALSE, quote = FALSE, sep = "\t")
  fwrite(df2, file2, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Run bedtools intersect using system()
  output_file <- file.path(temp_dir, "output.bed")
  bedtools_cmd <- paste0("bedtools intersect -a ", file1, " -b ", file2, " -wa -wb > ", output_file)
  system(bedtools_cmd)
  
  # Read the output BED file back into R as a data frame
  output_df <- fread(output_file)
  
  # Clean up temporary files if needed
  file.remove(file1, file2, output_file)
  return(output_df)
}

getSequence <- function(gr=merge.bbmap.results.gr[as.character(strand(merge.bbmap.results.gr))=="+"], genome=Hsapiens, ucscToEnsembl=1, upflank=0, strand="+", addMetadata=TRUE)
{
  if(addMetadata){
    values(gr)$seq<-rep( "null", length(gr) )
  } else{
    mcols(gr) <- DataFrame(seq=rep( "null", length(gr) ))  
  }
  
  for(chr in (sort(unique(seqnames(gr))))){
    # message(chr)
    gr.select <- gr[seqnames(gr)==chr]
    if(ucscToEnsembl){
      chr.ucsc<-noquote(paste0("chr",chr))
    } else{ chr.ucsc<-chr }
    if(chr=="MT" || chr=="chrMT"){chr.ucsc="chrM"}
    
    chr.seq <- genome[[chr.ucsc]]
    
    if(upflank>0){
      
      if(strand=="+"){
        
        start <- start(gr.select) - upflank
        
        start[start<1] <- 1
        seq <- DNAStringSet(Views(chr.seq, start=start, end=end(gr.select)))
        
      } else{
        
        end <- end(gr.select) + upflank
        
        end[end>length(chr.seq)] <- length(chr.seq)
        seq <- DNAStringSet(Views(chr.seq, start=start(gr.select), end=end))
        
      }
      
    } else {
      seq <- DNAStringSet(Views(chr.seq, start=start(gr.select), end=end(gr.select)))  
    }
    
    if(strand=="-"){
      seq <- reverseComplement(DNAStringSet(seq))  
    }
    
    gr[seqnames(gr)==chr]$seq <- seq
    
  }
  
  return (gr);
}

# Function to read transcript names from a file
readTranscriptFile <- function(file_path) {
  return(readLines(file_path))
}

# Function to retrieve information for specified transcripts
getTranscriptInfo <- function(list = transcript_file, version = version, dir = dir) {
  
  # download data from ensembl
  #ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = version)
  
  cat("Connecting to Ensembl via biomaRt...\n")
  ensembl <- tryCatch({
    useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        version = version
    )
  }, error = function(e) {
    cat("Note: Using default Ensembl version (connection issue with specified version)\n")
    useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = version)
  })
   
  # Specify the transcript names you are interested in
  target_transcripts <- readTranscriptFile(list)
  types <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'ensembl_transcript_id', 'transcript_biotype'),
                 filters = 'ensembl_transcript_id',
                 values = unique(target_transcripts),
                 mart = ensembl)
  
  UTR3 <- getBM(attributes = c('chromosome_name', '3_utr_start', '3_utr_end', 'ensembl_exon_id', 'rank', 'strand', 'ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'gene_biotype'),
                filters = 'ensembl_transcript_id',
                values = unique(target_transcripts), 
                mart = ensembl)
  
  UTR5 <- getBM(attributes = c('chromosome_name', '5_utr_start', '5_utr_end', 'ensembl_exon_id', 'rank', 'strand', 'ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'gene_biotype'),
                filters = 'ensembl_transcript_id',
                values = unique(target_transcripts), 
                mart = ensembl)
  
  CDS <- getBM(attributes = c('chromosome_name', 'genomic_coding_start', 'genomic_coding_end', 'ensembl_exon_id', 'rank', 'strand', 'ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'gene_biotype'),
               filters = 'ensembl_transcript_id',
               values = unique(target_transcripts), 
               mart = ensembl)
  
  # Fix exon file
  cols <- c("chromosome", "start", "end", "exon.id", "exon.rank", "strand", "gene.id", "gene.name", "transcript", "gene.type")
  
  clean_bed_convert <- function(annotation = UTR5, type = "UTR5", cols, biotypes = types) {
    annotation <- annotation[complete.cases(annotation[, 2:3]), ]
    colnames(annotation) <- cols
    annotation$type = type
    annotation$strand <- as.character(annotation$strand)
    annotation[strand == "-1", strand := "-"]
    annotation[strand == "1", strand := "+"]
    annotation$transcript.id = paste0(annotation$gene.id, "@", annotation$gene.name, "@", annotation$transcript)
    
    setnames(biotypes, c(3, 4), c("transcript", "transcript.type"))
    annotation <- merge(annotation, biotypes[, c("transcript", "transcript.type")], by = "transcript", all.x = TRUE)
    annotation$biotype = paste0(annotation$gene.type, "|", annotation$transcript.type)
    return(annotation[, c("chromosome", "start", "end", "exon.rank", "strand", "transcript.id", "exon.id", "biotype", "type")])
  }
  
  ensembl.UTR5 <- clean_bed_convert(as.data.table(UTR5), "UTR5", cols, types)
  ensembl.UTR3 <- clean_bed_convert(as.data.table(UTR3), "UTR3", cols, types)
  ensembl.CDS <- clean_bed_convert(as.data.table(CDS), "CDS", cols, types)
  
  exons <- as.data.table(rbind(ensembl.CDS, ensembl.UTR3, ensembl.UTR5))
  exons$chromosome <- paste0("chr",exons$chromosome)
  exons[chromosome == "chrMT", chromosome := "chrM"]
  exons$set <- '1'
  exons <- exons[, c('chromosome', 'start', 'end', 'transcript.id', 'set', 'strand', 'type')]
  write.table(exons, file = paste0(dir, "/exons.wildtype.tab"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  # Create fasta for transcripts
  gets_spliced_seq <- function(exons.table = exons, output.file = fi2) {
    exons.table <- exons
    exons.plus <- exons.table[strand == "+"]
    exons.minus <- exons.table[strand == "-"]
    exons.plus.gr <- makeGRangesFromDataFrame(as.data.frame(exons.plus), seqnames.field = "chromosome", start.field = "start",
                                              end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
    exons.minus.gr <- makeGRangesFromDataFrame(as.data.frame(exons.minus), seqnames.field = "chromosome", start.field = "start",
                                               end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
    gen <- "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"
    exons.plus.gr <- getSequence(exons.plus.gr, genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, ucscToEnsembl = 0, upflank = 0, strand = "+")
    exons.minus.gr <- getSequence(exons.minus.gr, genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, ucscToEnsembl = 0, upflank = 0, strand = "-")
    
    exons.plus.dt <- setDT(as.data.frame(exons.plus.gr))
    exons.minus.dt <- setDT(as.data.frame(exons.minus.gr))
    
    exons.grouped.plus <- exons.plus.dt %>% group_by(transcript.id, seqnames, strand) %>%
      arrange(start) %>%
      summarise(transcript.seq = paste(seq, collapse = "")) %>%
      ungroup() %>% as.data.table()
    
    exons.grouped.minus <- exons.minus.dt %>% group_by(transcript.id, seqnames, strand) %>%
      arrange(desc(start)) %>%
      summarise(transcript.seq = paste(seq, collapse = "")) %>%
      ungroup() %>% as.data.table()
    
    transcripts <- rbind(exons.grouped.minus, exons.grouped.plus)
    transcripts$header = paste0(">", transcripts$seqnames, ":", transcripts$transcript.id, "|", transcripts$strand)
    write.table(transcripts[, c("header", "transcript.seq")], file = output.file, sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
    return(transcripts)
  }
  
  fi2 <- paste0(dir, "/wildtype.fa")
  transcripts.seq <- gets_spliced_seq(exons, output.file = fi2)
}

getmutated <- function(exons = exons.wildtype.tab, snps = ssm.gtf.bed, dir = dir) {
  setnames(exons, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7'), c('chromosome', 'start', 'end', 'transcript.id', 'set', 'strand', 'type'))
  exons.plus <- exons[strand == "+"]
  exons.minus <- exons[strand == "-"]
  exons.plus.gr <- makeGRangesFromDataFrame(as.data.frame(exons.plus), seqnames.field = "chromosome", start.field = "start",
                                            end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
  exons.minus.gr <- makeGRangesFromDataFrame(as.data.frame(exons.minus), seqnames.field = "chromosome", start.field = "start",
                                             end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
  
  exons.plus.gr <- getSequence(exons.plus.gr, genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, ucscToEnsembl = 0, upflank = 0, strand = "+")
  exons.minus.gr <- getSequence(exons.minus.gr, genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, ucscToEnsembl = 0, upflank = 0, strand = "-")
  
  exons.plus.dt <- setDT(as.data.frame(exons.plus.gr))
  exons.minus.dt <- setDT(as.data.frame(exons.minus.gr))
  
  exons <- rbind(exons.plus.dt, exons.minus.dt)
  rm(exons.plus.dt, exons.minus.dt, exons.plus.gr, exons.minus.gr, exons.plus, exons.minus)
  
  # Process mutations grouped by gene
  new_coords <- data.table()
  snps[, unique_id := V8]  
  # Store unique mutation IDs
  snps_grouped <- snps[, .(mut_ids = paste(unique(V8), collapse = "_")), by = .(V12, V15, V16)]
  
  for (i in 1:nrow(snps_grouped)) {
    gene_id <- snps_grouped[i, V12]
    transcript <- paste0(snps_grouped[i, V15], '@', gene_id, '@', snps_grouped[i, V16])
    mutation_ids <- snps_grouped[i, mut_ids]
    
    coords <- exons[transcript.id == transcript]
    
    if (nrow(coords) > 0) {
      snp_subset <- snps[V12 == gene_id]
      
      for (j in 1:nrow(snp_subset)) {
        ssm.pos <- snp_subset[j, V2]
        genomic.strand <- snp_subset[j, V14]
        mutation <- snp_subset[j, V7]
        
        result_row <- coords[which(ssm.pos >= coords$start & ssm.pos <= coords$end), ]
        if (nrow(result_row) > 0) {
          if (genomic.strand == "+") {
            substr(result_row$seq, (ssm.pos - result_row$start)+1, (ssm.pos - result_row$start)+1) <- mutation
          } else {
            reversed_sequence <- rev(strsplit(result_row$seq, NULL)[[1]])
            reversed_sequence <- paste(reversed_sequence, collapse = "")
            reversed_sequence <- chartr("ATGC","TACG",reversed_sequence)
            
            substr(reversed_sequence, (ssm.pos - result_row$start)+1, (ssm.pos - result_row$start)+1) <- mutation
            reversed_sequence <- rev(strsplit(reversed_sequence, NULL)[[1]])
            reversed_sequence <- paste(reversed_sequence, collapse = "")
            reversed_sequence <- chartr("ATGC","TACG",reversed_sequence)
            
            result_row$seq <- reversed_sequence
          }
          
          coords[start == result_row$start, seq := result_row$seq]
        }
      }
      
      # Append all mutation IDs to the transcript ID
      coords[, transcript.id := paste0(transcript.id, '_', mutation_ids)]
      new_coords <- rbind(new_coords, coords)
    }
  }
  
  exons.plus.2 <- new_coords[strand == "+"]
  exons.minus.2 <- new_coords[strand == "-"]
  
  exons.grouped.plus <- exons.plus.2 %>% group_by(transcript.id, seqnames, strand) %>%
    arrange(start) %>%
    summarise(transcript.seq = paste(seq, collapse = "")) %>%
    ungroup() %>% as.data.table()
  
  exons.grouped.minus <- exons.minus.2 %>% group_by(transcript.id, seqnames, strand) %>%
    arrange(desc(start)) %>%
    summarise(transcript.seq = paste(seq, collapse = "")) %>%
    ungroup() %>% as.data.table()
  
  transcripts <- rbind(exons.grouped.minus, exons.grouped.plus)
  transcripts$header = paste0(">", transcripts$seqnames, ":", transcripts$transcript.id, "|", transcripts$strand)
  write.table(transcripts[, c("header", "transcript.seq")], file = 'mutated.fa', sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  new_coords <- new_coords[, c('seqnames', 'start', 'end', 'transcript.id', 'set', 'strand', 'type')]
  write.table(new_coords, file = 'exons.mutated.tab', sep='\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
}

cat("\n>>> Analysis Start Time: ")
system("date")

##### Read necessary resources #####

ssm.bed <- fread(opt$input, sep = "\t", col.names = c("chromosome","start","end","status","strand","reference","allele","sample"))

ssm.bed$chromosome <- gsub("chrMT", "chrM", ssm.bed$chromosome)

# miRBase miRNA annotation (v22.1)
#mirbase <- as.data.frame(import("/mnt/raid1/thanos/metallaso/resources/hsa_miRBase_22v1_0based.gff3"))

mirbase <- as.data.frame(import(opt$mirbase_input))
mirbase <- mirbase[mirbase$type == "miRNA",]

mirbase.bed <- mirbase[,c("seqnames","start","end","Name","score","strand")]

if (is.factor(mirbase.bed$seqnames)) {
  levels(mirbase.bed$seqnames) <- gsub("chrMT", "chrM", levels(mirbase.bed$seqnames))
} else {
  mirbase.bed$seqnames <- gsub("chrMT", "chrM", mirbase.bed$seqnames)
}
# Ensembl Genes annotation (v110)
#ensembl.homo.sapiens.gtf <- as.data.frame(import("/mnt/raid1/thanos/metallaso/resources/Homo_sapiens.GRCh38.110.chr.gtf"))
ensembl.homo.sapiens.gtf <- as.data.frame(import(opt$ensembl_input))
ensembl.homo.sapiens.gtf <- ensembl.homo.sapiens.gtf[ensembl.homo.sapiens.gtf$gene_biotype %in% c("protein_coding", "lncRNA", "pseudogene"),]
ensembl.homo.sapiens.gtf <- ensembl.homo.sapiens.gtf[ensembl.homo.sapiens.gtf$type %in% c("transcript", "three_prime_utr"), ]
ensembl.homo.sapiens.gtf <- ensembl.homo.sapiens.gtf[ensembl.homo.sapiens.gtf$transcript_biotype %in% c("protein_coding", "lncRNA", "pseudogene"),]

ensembl.homo.sapiens.bed <- ensembl.homo.sapiens.gtf[,c("seqnames","start","end","gene_name","score","strand","gene_id","transcript_id")]

if (is.factor(ensembl.homo.sapiens.bed$seqnames)) { # Assuming seqnames is the chromosome col name
  levels(ensembl.homo.sapiens.bed$seqnames) <- gsub("chrMT", "chrM", levels(ensembl.homo.sapiens.bed$seqnames))
} else {
  ensembl.homo.sapiens.bed$seqnames <- gsub("chrMT", "chrM", ensembl.homo.sapiens.bed$seqnames)
}

rm("mirbase", "ensembl.homo.sapiens.gtf")
cat("\n>>> Resources loaded successfully.\n")

##### Correct all input chromosome fields to be in the same form #####

if (all(grepl("chr", ssm.bed$chromosome))) {
  if (!all(grepl("chr", mirbase.bed$seqnames)))
    mirbase.bed$seqnames <- paste0("chr", mirbase.bed$seqnames)
  else if (!all(grepl("chr", ensembl.homo.sapiens.bed$seqnames)))
    ensembl.homo.sapiens.bed$seqnames <- paste0("chr", ensembl.homo.sapiens.bed$seqnames)
} else {
  if (any(grepl("chr", mirbase.bed$seqnames)))
    mirbase.bed$seqnames <- gsub("chr", "", mirbase.bed$seqnames)
  else if (any(grepl("chr", ensembl.homo.sapiens.bed$seqnames)))
    ensembl.homo.sapiens.bed$seqnames <- gsub("chr", "", ensembl.homo.sapiens.bed$seqnames)
}

##### Convert all input SSMs in "+" strand #####

for (ssm in 1:nrow(ssm.bed)) {
  if (ssm.bed[ssm,]$strand == "-") {
    ssm.bed[ssm,]$reference <- reverse(chartr("ATCG", "TAGC", ssm.bed[ssm,]$reference))
    ssm.bed[ssm,]$allele <- reverse(chartr("ATCG", "TAGC", ssm.bed[ssm,]$allele))
    ssm.bed[ssm,]$strand <- "+"
  }
}

print("Checking for chrMT in ssm.bed:")
print(unique(ssm.bed$chromosome[grepl("chrMT", ssm.bed$chromosome)]))
print("Checking for chrMT in mirbase.bed:")
print(unique(mirbase.bed$seqnames[grepl("chrMT", mirbase.bed$seqnames)]))
print("Checking for chrMT in ensembl.homo.sapiens.bed:")
print(unique(ensembl.homo.sapiens.bed$seqnames[grepl("chrMT", ensembl.homo.sapiens.bed$seqnames)]))

##### Intersect all inputs #####

if (opt$analysis == 'mre' ) {
  ssm.mirbase.bed <- bed_intersect_data_frames(ssm.bed, mirbase.bed)
  ssm.gtf.bed <- bed_intersect_data_frames(ssm.bed, ensembl.homo.sapiens.bed)
  
  ssm.intersection <- rbind(ssm.mirbase.bed, ssm.gtf.bed, fill = TRUE)
} else {
  ssm.intersection <- bed_intersect_data_frames(ssm.bed, ensembl.homo.sapiens.bed)
}
fwrite(ssm.intersection, "ssm_intersection.bed", col.names = FALSE, quote = FALSE, sep = "\t")
cat("\n>>> Inputs intersected successfully.\n")


##### Split execution based on input #####

if (opt$analysis == 'splicing') {
  
  ##### Convert SSM intersection to SpliceAI input vcf #####
  
  ssm.intersection.vcf = ssm.intersection[,c(1,2,6,7)]
  ssm.intersection.vcf = cbind(ssm.intersection.vcf,'.','.','.','.')
  ssm.intersection.vcf = ssm.intersection.vcf[,c(1,2,5,3,4,7,8,6)]
  colnames(ssm.intersection.vcf) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  ssm.intersection.vcf$ID = paste0(ssm.intersection$V8, ",", ssm.intersection$V4)
  ssm.intersection.vcf$INFO = paste0("Gene=",ssm.intersection$V12, ",", ssm.intersection$V10, "-", ssm.intersection$V11)
  ssm.intersection.vcf$CHROM <- gsub("chr", "", ssm.intersection.vcf$CHROM)
  
  ## Create, format and fill the SpliceAI vcf input
  write("##fileformat=VCFv4.2", "ssm_intersection.vcf")
  write(paste0("##fileDate=",Sys.Date()), "ssm_intersection.vcf", append = TRUE)
  write("##reference=GRCh38/hg38", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=1,length=248956422>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=2,length=242193529>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=3,length=198295559>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=4,length=190214555>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=5,length=181538259>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=6,length=170805979>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=7,length=159345973>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=8,length=145138636>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=9,length=138394717>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=10,length=133797422>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=11,length=135086622>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=12,length=133275309>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=13,length=114364328>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=14,length=107043718>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=15,length=101991189>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=16,length=90338345>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=17,length=83257441>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=18,length=80373285>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=19,length=58617616>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=20,length=64444167>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=21,length=46709983>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=22,length=50818468>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=X,length=156040895>", "ssm_intersection.vcf", append = TRUE)
  write("##contig=<ID=Y,length=57227415>", "ssm_intersection.vcf", append = TRUE)
  write("##INFO=<ID=Gene,Number=.,Type=String,Description='Gene Info'>", "ssm_intersection.vcf", append = TRUE)
  write(paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",sep = "\t"), "ssm_intersection.vcf", append = TRUE)
  fwrite(ssm.intersection.vcf, "ssm_intersection.vcf", col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  cat("\n>>> SSM Intersection VCF generated successfully.\n")
  
  ## Execute SpliceAI
  system(paste0("spliceai -I ssm_intersection.vcf -O spliceAI_out.vcf -R ", opt$splicing_genome, " -A ", opt$splicing_annotation, " -D 50"), intern = TRUE)
  cat("\n>>> SpliceAI executed successfully.\n")
  
  ## Filter and store results using the provided splicing_probability_threshold
  splicing_threshold = opt$splicing_probability_threshold
  spliceAI_results = read.table("spliceAI_out_partial.vcf", sep = "\t", col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
  spliceAI_results[c("GENE_INFO","SPLICE_INFO")] = str_split_fixed(spliceAI_results$INFO, ';', 2)
  spliceAI_results[c("GENE","GENE_COORDS")] = str_split_fixed(spliceAI_results$GENE_INFO, ',', 2)
  spliceAI_results["GENE"] = str_split_fixed(spliceAI_results$GENE, '=', 2)[,2]
  spliceAI_results[c("GENE_START","GENE_END")] = str_split_fixed(spliceAI_results$GENE_COORDS, '-', 2)
  spliceAI_results = spliceAI_results[spliceAI_results$SPLICE_INFO != "",]
  spliceAI_results["SPLICE_INFO"] = str_split_fixed(spliceAI_results$SPLICE_INFO, '=', 2)[,2]
  spliceAI_results = spliceAI_results[c("CHROM","POS","ID","REF","ALT","GENE","GENE_START","GENE_END","SPLICE_INFO")]
  
  splicing_filtered_results = data.frame(matrix(ncol = 17, nrow = 0))
  colnames(splicing_filtered_results) = c("CHROM","POS","ID","REF","ALT","GENE","GENE_START","GENE_END","TRANSCRIPT","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
  splicing_significant_mutations = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(splicing_filtered_results) = c("CHROM","POS","ID","REF","ALT","GENE","GENE_START","GENE_END","SPLICE_INFO")
  
  for (mutation_idx in 1:nrow(spliceAI_results)) {
    mutation_splicing_events = str_split(spliceAI_results$SPLICE_INFO[mutation_idx],',', simplify = TRUE)
    for (mutation_splicing_event_idx in 1:ncol(mutation_splicing_events)) {
      tmp = subset(spliceAI_results[mutation_idx,], select = -SPLICE_INFO)
      significant_events = 0
      tmp[c("ALLELE","TRANSCRIPT","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")] = str_split_fixed(mutation_splicing_events[mutation_splicing_event_idx], "\\|", 10)
      tmp = tmp[c("CHROM","POS","ID","REF","ALT","GENE","GENE_START","GENE_END","TRANSCRIPT","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")]
      if ((tmp$DS_AG > splicing_threshold) | (tmp$DS_AL > splicing_threshold) | (tmp$DS_DG > splicing_threshold) | (tmp$DS_DL > splicing_threshold)) {
        splicing_filtered_results = rbind(splicing_filtered_results,tmp)
        significant_events = significant_events + 1
      }
    }
    if (significant_events > 0) {
      spliceAI_results$SPLICE_INFO[mutation_idx] = significant_events
      splicing_significant_mutations = rbind(splicing_significant_mutations,spliceAI_results[mutation_idx,])
    }
  }
  rm(tmp)
  colnames(splicing_significant_mutations)[9] = "SIGNIFICANT_SPLICING_EVENTS"
  splicing_significant_mutations <- splicing_significant_mutations[order(splicing_significant_mutations$SIGNIFICANT_SPLICING_EVENTS,decreasing=TRUE),]
  fwrite(splicing_filtered_results, "significant_splicing_events.txt", col.names = TRUE, quote = FALSE, sep = "\t")
  fwrite(splicing_significant_mutations, "significant_splicing_mutations.txt", col.names = TRUE, quote = FALSE, sep = "\t")
  
} else if (opt$analysis == 'orf') {
  
  ##### Generate the wildtype and mutated genomic sequences and write them to separate FASTA files #####
  
  file.create("wildtype_sequences.fasta")
  file.create("mutated_sequences.fasta")
  
  for (bed.line.index in 1:nrow(ssm.intersection)) {
    if (!grepl("chr", ssm.intersection[bed.line.index,]$V1)) {
      chromosome <- paste0("chr", ssm.intersection[bed.line.index,]$V1)
    } else {
      chromosome <- ssm.intersection[bed.line.index,]$V1
    }
    
    ssm.pos <- ssm.intersection[bed.line.index,]$V2
    ssm.strand <- ssm.intersection[bed.line.index,]$V5
    
    genomic.start <- ssm.intersection[bed.line.index,]$V10
    genomic.end <- ssm.intersection[bed.line.index,]$V11
    genomic.strand <- ssm.intersection[bed.line.index,]$V14
    
    genomic.sequence <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, chromosome, genomic.start, genomic.end))
    
    if (genomic.strand == "+") {
      mutated.genomic.sequence = genomic.sequence
      substr(mutated.genomic.sequence, (genomic.end - ssm.pos) + 1, (genomic.end - ssm.pos) + 1) <- ssm.intersection[bed.line.index,]$V7
    } else {
      genomic.sequence <- reverse(chartr("ATCG", "TAGC", genomic.sequence))
      mutated.genomic.sequence = genomic.sequence
      substr(mutated.genomic.sequence, (genomic.end - ssm.pos) + 1, (genomic.end - ssm.pos) + 1) <- ssm.intersection[bed.line.index,]$V7
    }
    
    fasta_name_line = paste0(ssm.intersection[bed.line.index,]$V12,":",ssm.intersection[bed.line.index,]$V16," - ",ssm.intersection[bed.line.index,]$V9,", ",ssm.intersection[bed.line.index,]$V10,"-",ssm.intersection[bed.line.index,]$V11," - ",ssm.intersection[bed.line.index,]$V8,", ",ssm.intersection[bed.line.index,]$V4)
    write.fasta(genomic.sequence, fasta_name_line, "wildtype_sequences.fasta", open = "a", nbchar = 100, as.string = TRUE)
    
    fasta_mutated_name_line = paste0(fasta_name_line," - ",ssm.intersection[bed.line.index,]$V2,", ",ssm.intersection[bed.line.index,]$V6,"-",ssm.intersection[bed.line.index,]$V7)
    write.fasta(mutated.genomic.sequence, fasta_mutated_name_line, "mutated_sequences.fasta", open = "a", nbchar = 100, as.string = TRUE)
  }
  cat("\n>>> Wildtype and Mutated sequences FASTA files generated successfully.\n")
  
  ##### Execute NCBI ORFfinder and filter results #####
  
  write(paste0("FROM WILDTYPE FASTA: \n"), "transcripts_without_orfs.txt")
  system(paste0("./resources/ORFfinder -in wildtype_sequences.fasta -ml ", opt$orf_min_length," -out orffinder_wildtype_results.fasta 2>> transcripts_without_orfs.txt"))
  write(paste0("\nFROM MUTATED FASTA: \n\n"), "transcripts_without_orfs.txt", append = TRUE)
  system(paste0("./resources/ORFfinder -in mutated_sequences.fasta -ml ", opt$orf_min_length," -out orffinder_mutated_results.fasta 2>> transcripts_without_orfs.txt"))
  cat("\n>>> NCBI ORF Finder executed successfully.\n")
  
  ## Aggregate the ORFs found in both conditions
  system("grep '>' orffinder_wildtype_results.fasta > tmp.txt")
  wildtype_orfs = read.table("tmp.txt",sep = "|")
  wildtype_orfs$V2 = lapply(wildtype_orfs$V2, function(x) sub(" .*", "", x))
  wildtype_orfs = unique(wildtype_orfs$V2)
  
  system("grep '>' orffinder_mutated_results.fasta > tmp.txt")
  mutated_orfs = read.table("tmp.txt",sep = "|")
  mutated_orfs$V2 = lapply(mutated_orfs$V2, function(x) sub(" .*", "", x))
  mutated_orfs = unique(mutated_orfs$V2)
  invisible(file.remove("tmp.txt"))
  
  ## Iterate over the wildtype ORFs list and remove any common ORFs from both lists. What remains are uniques.
  for (orf.index in length(wildtype_orfs):1) {
    if (wildtype_orfs[[orf.index]] %in% mutated_orfs) {
      mutated_orfs = mutated_orfs[mutated_orfs != wildtype_orfs[[orf.index]]]
      wildtype_orfs = wildtype_orfs[-orf.index]
    }
  }
  
  ## Report findings
  if ((length(wildtype_orfs) == 0) & (length(mutated_orfs) == 0)){
    cat("\n>>> Analysis Complete. NO UNIQUE ORFs FOUND!\n")
  } else {
    file.create("unique_ORFs.txt")
    if (length(wildtype_orfs) > 0) {
      write(paste0("\nWILDTYPE UNIQUE ORFs: \n"), "unique_ORFs.txt", append = TRUE)
      for (orf in wildtype_orfs) {
        write(paste0(orf), "unique_ORFs.txt", append = TRUE)
      }
    }
    if (length(mutated_orfs) > 0) {
      write(paste0("\nMUTATED UNIQUE ORFs: \n"), "unique_ORFs.txt", append = TRUE)
      for (orf in mutated_orfs) {
        write(paste0(orf), "unique_ORFs.txt", append = TRUE)
      }
    }
    cat("\n>>> Analysis Complete. Unique ORFs can be found in the 'unique_ORFs.txt' file generated.\n")
  }
  
} else if (opt$analysis == 'mre') {
  tarBase.v9 <- fread(opt$mre_interactions)
  tarBase.v9$dummy <- "*"
  tarBase.v9 <- tarBase.v9[, c("chromosome", "start", "end", "dummy", "strand", "mirna_name", "transcript_id")]
  
  tarbase.ssm.intersect <- bed_intersect_data_frames(tarBase.v9, ssm.bed)
  tarbase.ssm.intersect <- tarbase.ssm.intersect[tarbase.ssm.intersect$V1 != "chrMT",]
  #tarbase.ssm.intersect <- tarbase.ssm.intersect[tarbase.ssm.intersect$V3 - tarbase.ssm.intersect$V2 == 12,]
  if (nrow(tarbase.ssm.intersect) == 0) {
    stop("0 transcripts overlapped by input SSMs. Exiting.")
  }
  filtered_data <- unique(tarbase.ssm.intersect[,c(6,7)])
  
  ##### If SSMs for miRNAs exist in input, Generate the wildtype and mutated genomic sequences and write them to separate FASTA files #####
  
  # If miRNAs are included in the SSMs:
  if (length(ssm.mirbase.bed) > 0) {
    file.create("wildtype_mirna_sequences.fasta")
    file.create("mutated_mirna_sequences.fasta")
    
    for (bed.line.index in 1:nrow(ssm.mirbase.bed)) {
      if (!grepl("chr", ssm.mirbase.bed[bed.line.index,]$V1)) {
        chromosome <- paste0("chr", ssm.mirbase.bed[bed.line.index,]$V1)
      } else {
        chromosome <- ssm.mirbase.bed[bed.line.index,]$V1
      }
      
      ssm.pos <- ssm.mirbase.bed[bed.line.index,]$V2
      ssm.strand <- ssm.mirbase.bed[bed.line.index,]$V5
      
      genomic.start <- ssm.mirbase.bed[bed.line.index,]$V10
      genomic.end <- ssm.mirbase.bed[bed.line.index,]$V11
      genomic.strand <- ssm.mirbase.bed[bed.line.index,]$V14
      mirna.name <- ssm.mirbase.bed[bed.line.index,]$V12
      
      if (mirna.name %in% filtered_mirnas$V1) {
        genomic.sequence <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, chromosome, genomic.start, genomic.end))
        
        if (genomic.strand == "+") {
          mutated.genomic.sequence = genomic.sequence
          substr(mutated.genomic.sequence, (genomic.end - ssm.pos) + 1, (genomic.end - ssm.pos) + 1) <- ssm.mirbase.bed[bed.line.index,]$V7
        } else {
          genomic.sequence <- reverse(chartr("ATCG", "TAGC", genomic.sequence))
          mutated.genomic.sequence = genomic.sequence
          substr(mutated.genomic.sequence, (genomic.end - ssm.pos) + 1, (genomic.end - ssm.pos) + 1) <- ssm.mirbase.bed[bed.line.index,]$V7
        }
        
        fasta_name_line = paste0(ssm.mirbase.bed[bed.line.index,]$V12)
        write.fasta(genomic.sequence, fasta_name_line, "wildtype_mirna_sequences.fasta", open = "a", nbchar = 100, as.string = TRUE)
        
        fasta_mutated_name_line = paste0(fasta_name_line)
        write.fasta(mutated.genomic.sequence, fasta_mutated_name_line, "mutated_mirna_sequences.fasta", open = "a", nbchar = 100, as.string = TRUE)
      }
    }
  } else {
    ## If there are no SSMs for miRNAs, Get the unique miRNAs from the interacting pairs of Tarbase9 and gather their wildtype sequences for the microT input
    # Get unique miRNAs and process them
    filtered_mirnas <- as.data.frame(filtered_data[, 1, with = FALSE])
    colnames(filtered_mirnas) <- "mirna"
    unique_filtered_mirnas <- unique(filtered_mirnas)
    file.create("wildtype_mirna_sequences.fasta")
    
    cat("Number of unique miRNAs to process after filtering:", nrow(unique_filtered_mirnas), "\n")
    
    # Filter out miRNAs not annotated in miRBase v22.1
    unique_filtered_mirnas <- unique_filtered_mirnas[unique_filtered_mirnas$mirna %in% mirbase.bed$Name, , drop = FALSE]
    cat("Number of miRNAs after filtering against miRBase:", nrow(unique_filtered_mirnas), "\n")
    
    if (nrow(unique_filtered_mirnas) == 0) {
      cat("Warning: No miRNAs from filtered data found in miRBase v22.1\n")
    } else {
      for (i in 1:nrow(unique_filtered_mirnas)) {
        current_mirna <- unique_filtered_mirnas[i, "mirna"]
        #cat("Processing miRNA:", current_mirna, "\n")
        
        mirIndex <- mirbase.bed[mirbase.bed$Name == current_mirna, ]
        
        if (nrow(mirIndex) == 0) {
          cat("Warning: miRNA", current_mirna, "not found in mirbase.bed. Skipping.\n")
          next
        }
        
        #cat("  Found", nrow(mirIndex), "match(es)\n")
        
        genomic.start <- mirIndex[1, ]$start
        genomic.end <- mirIndex[1, ]$end
        chromosome <- mirIndex[1, ]$seqnames
        
        if (is.na(chromosome) || is.na(genomic.start) || is.na(genomic.end)) {
          cat("Warning: NA values found for miRNA", current_mirna, ". Skipping.\n")
          next
        }
        
        #cat("  Chromosome:", chromosome, "Start:", genomic.start, "End:", genomic.end, "\n")
        
        genomic.sequence <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, chromosome, genomic.start, genomic.end))
        fasta_name_line <- paste0(mirIndex[1, ]$Name)
        write.fasta(genomic.sequence, fasta_name_line, "wildtype_mirna_sequences.fasta", open = "a", nbchar = 100, as.string = TRUE)
        #cat("  Successfully processed miRNA:", current_mirna, "\n")
      }
    }
  }
  
  file.create("transcript_id_list.txt")
  
  colnames(filtered_data) <- c("mirna", "transcript_id")
  filtered_transcripts <- filtered_data$transcript_id
  for (bed.line.index in 1:nrow(ssm.gtf.bed)) {
    transcript.id <- ssm.gtf.bed[bed.line.index,]$V16
    
    if (transcript.id %in% filtered_transcripts) {
      write(transcript.id, "transcript_id_list.txt", sep="\n", append = TRUE)
    }
  }
  
  # set inputs
  version <- '113'
  #dir <- file.path(getwd(), "microt_temp")
  dir <- getwd()
  # if (!dir.exists(dir)) {
  #   dir.create(dir, recursive = TRUE)
  # }
  
  transcript_file <- paste0(dir, '/transcript_id_list.txt')
  
  # for wildtype 
  # two files created (wildtype.fa, exons.wildtype.tab)
  getTranscriptInfo(transcript_file, version, dir)
  
  # for mutated
  exons.wildtype.tab <- fread(paste0(dir, '/exons.wildtype.tab'))
  ssm.gtf.bed <- fread(paste0(dir, "/ssm_intersection.bed"))
  # two files created (mutated.fa, exons.mutated.tab)
  getmutated(exons.wildtype.tab, ssm.gtf.bed, dir)
  
  cat("\n>>> MicroT input files generated successfully.\n")
}

#system(paste0("docker run -v /mnt/raid1/theo/Documents/docker/omics:/microt_temp penny0lane/microt_cnn Rscript main.R /microt_temp/config.yml"))

cat("\n>>> Analysis End Time: ")
system("date")
