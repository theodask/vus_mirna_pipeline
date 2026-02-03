FROM rocker/tidyverse:4.0.3

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies + compiler
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    bedtools \
    samtools \
    bcftools \
    tabix \
    wget \
    curl \
    build-essential \
    g++ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip3 install --no-cache-dir \
    numpy==1.24.4 \
    matplotlib==3.7.1 \
    gffutils==0.13 \
    intervaltree==3.1.0

# CRAN packages
RUN R -e "install.packages(c('data.table','seqinr','optparse','glmnet','survival','survminer','dbplyr'))"

# Bioconductor packages
RUN R -e "BiocManager::install(version='3.12', update=FALSE, ask=FALSE)"
RUN R -e "BiocManager::install(c('Rhtslib','Rsamtools','Biostrings','GenomicRanges','rtracklayer','VariantAnnotation','BiocFileCache'), update=FALSE, ask=FALSE)"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38', update=FALSE, ask=FALSE)"

WORKDIR /pipeline

