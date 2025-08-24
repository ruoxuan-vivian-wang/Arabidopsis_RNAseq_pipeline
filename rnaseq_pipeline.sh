#!/bin/bash
# RNA-seq Analysis Pipeline
# Author: Ruoxuan
# Date: 2025-07-03
# Description: Fastq QC, genome indexing, alignment, and gene counting using STAR and related tools

# -----------------------------
# 1. Create and activate environment
# -----------------------------
# Install necessary software packages
micromamba create -n rnaseq_env python=3.10 -y --override-channels -c conda-forge -c bioconda fastp star subread samtools multiqc fastqc htseq

# Initialize micromamba shell integration
eval "$(micromamba shell hook --shell bash)"

# Activate the environment
micromamba activate rnaseq_env

# -----------------------------
# 2. Define input and output directories (user configurable)
# -----------------------------
RAW_DIR="raw_fastq"          # folder containing raw fastq files
CLEAN_DIR="clean_fastq"      # folder to save cleaned fastq files
STAR_INDEX_DIR="star_index"  # folder to store STAR genome index
STAR_OUT_DIR="star_out"      # folder for STAR alignment output
REF_GENOME="reference/genome.fasta"    # reference genome fasta
GTF_FILE="reference/annotation.gtf"   # genome annotation GTF

mkdir -p $CLEAN_DIR $STAR_INDEX_DIR $STAR_OUT_DIR

# -----------------------------
# 3. Fastq QC and adapter trimming using fastp
# -----------------------------
fastp \
  -i ${RAW_DIR}/sample_R1.fastq.gz \   # forward reads
  -I ${RAW_DIR}/sample_R2.fastq.gz \   # reverse reads
  -o ${CLEAN_DIR}/sample_clean_R1.fastq.gz \  # output cleaned forward reads
  -O ${CLEAN_DIR}/sample_clean_R2.fastq.gz \  # output cleaned reverse reads
  --detect_adapter_for_pe \                     # auto detect paired-end adapters
  -w 8 \                                       # threads (adjustable)
  -h ${CLEAN_DIR}/sample_fastp.html \          # HTML QC report
  -j ${CLEAN_DIR}/sample_fastp.json            # JSON QC report

# -----------------------------
# 4. Generate STAR genome index
# -----------------------------
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $STAR_INDEX_DIR \
     --genomeFastaFiles $REF_GENOME \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 99   # for 100bp reads

# -----------------------------
# 5. Align reads with STAR and count genes
# -----------------------------
STAR --runThreadN 8 \
     --genomeDir $STAR_INDEX_DIR \
     --readFilesIn ${CLEAN_DIR}/sample_clean_R1.fastq.gz ${CLEAN_DIR}/sample_clean_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${STAR_OUT_DIR}/sample_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile $GTF_FILE

# End of pipeline
