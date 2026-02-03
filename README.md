# VUS-miRNA Analysis Pipeline

A comprehensive Nextflow pipeline for identifying and characterizing Variants of Uncertain Significance (VUS) in microRNA binding regions, with integrated survival analysis.

# Overview

This pipeline analyzes genomic variants to identify those that disrupt microRNA (miRNA) binding sites, predicts their functional impact using microT-CNN, and performs survival analysis to assess clinical significance in cancer patients.

# Key Features

- Variant Analysis: Identifies variants in miRNA seed regions
- Variants Genomic localization (optional)
- MRE Disruption Prediction: Uses microT-CNN to predict binding affinity changes
- Survival Analysis: Cox regression and LASSO for clinical significance
- Containerization: Docker support for reproducibility
- Reporting: Automated generation of plots and statistics

# Prerequisites

- Nextflow (≥23.04.0)
- Docker
- 16GB+ RAM recommended
- 50GB+ disk space

# Installation

```
# Clone the repository
git clone https://github.com/theodask/vus_mirna_pipeline.git
cd vus_mirna_pipeline

# Setup directory structure
chmod +x setup_pipeline.sh
./setup_pipeline.sh

# Pull Docker image
docker pull theodask/vus_mirna_pipeline:1.0
```


# Running the Pipeline

## Single Cancer Type

```
nextflow run main.nf \
    -profile docker \
    --cancer_type BRCA \
    --vcf /path/to/merged_BRCA.vcf.gz \
    --rpm_threshold 100 \
    --skip_genomic_regions true \
    -resume
```

## Multiple Cancer Types (Sample Sheet)

Create `samples.csv`:
```csv
cancer_type,vcf_path
BRCA,/path/to/merged_BRCA.vcf.gz
BLCA,/path/to/merged_BLCA.vcf.gz
COAD,/path/to/merged_COAD.vcf.gz
```

Run:
```
nextflow run main.nf \
    -profile docker \
    --sample_sheet samples.csv \
    --rpm_threshold 100
    --skip_genomic_regions true \
```


# Required Files

## Download and Place in resources/ directory:

**Genomic Annotation** (`resources/genomic_annotation/`):
- `Homo_sapiens.GRCh38.113.gtf` - [Ensembl GTF](http://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/)
- `Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai` - [Ensembl FAI](http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/)

**microRNA Annotation** (`resources/microRNA_annotation/`):
- TarBase v9.0 MRE file with microT binding sites

**ClinVar Database** (`resources/variants/`):
- `clinvar.vcf` - [NCBI ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)

**miTED Expression Data** (`resources/miTED/[cancer_type]/`):
- Cancer-specific miRNA expression files - [miTED Database](http://www.mited.org/)
- Format: `miTED_top_expressed_miRNA_normal_and_cancer_RPM_[cancer_type].tsv`

**Clinical Metadata** (`resources/metadata/[cancer_type]/`):
- `clinical.cart.[cancer_type].tar.gz` - [GDC Portal](https://portal.gdc.cancer.gov/)

**microT-CNN References** (`microt_cnn/`):
- `hsa_miRBase_22v1_0based.gff3`
- `Homo_sapiens.GRCh38.113.gtf`


# Configuration

## Parameters

 Parameter | Description | Default |
 
 `--cancer_type` | Cancer type (BRCA, LUAD, etc.) | null |
 
 `--vcf` | Path to input VCF file | null |
 
 `--vcf_dir` | Directory for auto-discovery | null |
 
 `--sample_sheet` | CSV with cancer types and VCFs | null |
 
 `--rpm_threshold` | RPM threshold for miRNA expression | null |
 
 `--skip_genomic_regions` | Skip genomic region analysis | false |
 
 `--lines_per_chunk` | Variants per chunk | 200000 |
 
 `--outdir` | Output directory | results/ |
 
 `--base_dir` | Base directory for resources | pipeline dir |

# Resource Configuration

Modify `nextflow.config` to adjust resources:
```
process {
    withName: PROCESS_NAME {
        cpus = 8
        memory = '16 GB'
        time = '12 h'
    }
}
```

# Output Structure
```
results/
└── [cancer_type]/
    ├── genomic_regions/
    ├── mre_coordinates/
    ├── seed_regions/
    ├── intersected_variants/
    ├── filtered_variants/
    ├── clinvar_annotated/
    ├── mre_annotated/
    ├── microt_inputs/
    ├── microt_cnn_results/
    ├── microt_analysis/
    │   ├── microT_CNN_wt_mut_results_*.tsv
    │   ├── only_mutated_MREs_*.tsv
    │   ├── only_wildtype_MREs_*.tsv
    │   └── histogram_*.tiff
    └── survival_analysis/
        ├── survival_plots/
        │   └── survival_pvalues.tsv
        │   └── survival_variant_*.tiff
        └── lasso_analysis/
            ├── lasso_cox_results_*.tsv
            ├── cox_results_with_mutation_stats.tsv
            ├── lasso_selected_predictors.tsv
            ├── lasso_top_variants_barplot.png
            ├── mutation_frquency_distribution.png
            └── score_logrank_global_pvalue.tsv
```

# Docker Image
```
# Pull from Docker Hub
docker pull theodask/vus_mirna_pipeline:1.0

# Or build locally
docker build -t theodask/vus_mirna_pipeline:1.0 .
```

# Included Software

- **R 4.0.3**: tidyverse, Bioconductor, survival, survminer, glmnet
- **Python 3**: numpy, matplotlib, gffutils, intervaltree
- **Tools**: bedtools, samtools, bcftools, tabix


# Troubleshooting

## Common Issues

**Out of memory**
```
# Increase memory in nextflow.config
process.memory = '32 GB'
```

**Process timeout**
```
# Increase time limit
process.time = '48 h'
```

# Citation

If you use this pipeline in your research, please cite it as:

Thodoris Daskalopoulos.
**VUS-miRNA-Pipeline**: A Nextflow pipeline for identifying and characterizing variants of uncertain significance in miRNA binding regions.
GitHub repository: https://github.com/theodask/vus_mirna_pipeline
Version 1.0.0 (2025).

# License

MIT License (see LICENSE file)

# Author

Thodoris Daskalopoulos, MSc in Bioinformatics

# Acknowledgments

- [microT-CNN] (https://github.com/dianalabgr/microT-CNN)
Zacharopoulou E, Paraskevopoulou MD, Tastsoglou S, Alexiou A, Karavangeli A, Pierros V, Digenis S, Mavromati G, Hatzigeorgiou AG#, Karagkouni D#, microT-CNN: An avant-garde Deep Convolutional Neural Network unravels functional miRΝΑ targets beyond canonical sites.
- [TarBase v9] (https://dianalab.e-ce.uth.gr/tarbasev9)
Giorgos Skoufos, Panos Kakoulidis, Spyros Tastsoglou, Elissavet Zacharopoulou, Vasiliki Kotsira, Marios Miliotis, Galatea Mavromati, Dimitris Grigoriadis, Maria Zioga, Angeliki Velli, Ioanna Koutou, Dimitra Karagkouni, Steve Stavropoulos, Filippos S Kardaras, Anna Lifousi, Eustathia Vavalou, Armen Ovsepian, Anargyros Skoulakis, Sotiris K Tasoulis, Spiros V Georgakopoulos, Vassilis P Plagianakos, Artemis G Hatzigeorgiou, TarBase-v9.0 extends experimentally supported miRNAgene interactions to cell-types and virally encoded miRNAs, Nucleic Acids Research, 2023, DOI:  https://doi.org/10.1093/nar/gkad1071
- [miTED] (http://www.mited.org/)
Ioannis Kavakiotis, Athanasios Alexiou, Spyros Tastsoglou, Ioannis S Vlachos and Artemis G Hatzigeorgiou, DIANA-miTED: a microRNA tissue expression database, Nucleic Acids Research, 2021;, gkab733, https://doi.org/10.1093/nar/gkab733
- [ClinVar] (https://www.ncbi.nlm.nih.gov/clinvar/)
- [TCGA] (https://www.cancer.gov/tcga.)
- [Nextflow] (https://www.nextflow.io/) 
