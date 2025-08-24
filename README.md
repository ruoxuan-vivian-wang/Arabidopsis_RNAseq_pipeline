# 🧬 Arabidopsis RNA-seq Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.10-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

A reproducible and easy-to-use RNA-seq analysis pipeline for **quality control, alignment, and gene counting**. Designed for Linux/macOS environments with common bioinformatics tools.

## ⚡️ Features

- Quality control and adapter trimming with **fastp**
- Genome indexing and read alignment with **STAR**
- Gene quantification with **HTSeq / STAR GeneCounts**
- Multi-threaded support for large datasets
- Easy to adapt for different species or reference genomes

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone git clone https://github.com/ruoxuan-vivian-wang/Arabidopsis_RNAseq_pipeline.git
cd Arabidopsis_RNAseq_pipeline
```

### 2. Setup environment (Micromamba / Conda)

```bash
micromamba create -n rnaseq_env python=3.10 -y -c conda-forge -c bioconda fastp star subread samtools multiqc fastqc htseq gffread
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq_env
```

### 3. Prepare reference genome & annotation

```bash
gffread reference.gff -T -o reference.gtf
```

### 4. Generate STAR genome index

```bash
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ./STAR_index \
  --genomeFastaFiles reference.fasta \
  --sjdbGTFfile reference.gtf \
  --sjdbOverhang 99
```

### 5. Run Fastp for QC

```bash
fastp -i sample_R1.fq.gz -I sample_R2.fq.gz \
  -o sample_clean_R1.fq.gz -O sample_clean_R2.fq.gz \
  --detect_adapter_for_pe \
  -w 8 -h sample_fastp.html -j sample_fastp.json
```

### 6. Align reads with STAR & generate gene counts

```bash
STAR --runThreadN 8 \
  --genomeDir ./STAR_index \
  --readFilesIn sample_clean_R1.fq.gz sample_clean_R2.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix ./STAR_out/sample_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --sjdbGTFfile reference.gtf
```

## 📦 Output

Cleaned FASTQ files

QC reports (HTML + JSON)

Sorted BAM files

Gene count tables (ReadsPerGene.out.tab)

## ⚠️ Notes

Adjust --runThreadN based on your system.

Replace reference genome and annotation files with your species-specific files.

Modify file paths as needed.

## 📝 License

Released under the MIT License. Feel free to use, modify, and share.

Designed & maintained by Ruoxuan Wang.
