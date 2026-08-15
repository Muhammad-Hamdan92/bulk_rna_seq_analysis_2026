#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

PROJECT_DIR="/home/muhammad/Documents/Transcriptomics/transcriptomics_2026"
RAW_READS_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/00_raw_reads"
REF_GENOME_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/reference_genome"
RAW_READS_DIR="${PROJECT_DIR}/00_raw_reads"
REFERENCE_GENOME_DIR="${PROJECT_DIR}/00_reference_genome"

echo "Create project directories."
echo "---------------------------------------------------------------------------------------"

# Create project directories.
mkdir -p \
    "${RAW_READS_DIR}" \
    "${REFERENCE_GENOME_DIR}" \
    "${PROJECT_DIR}/01_qc_before_processing" \
    "${PROJECT_DIR}/02_process_reads" \
    "${PROJECT_DIR}/03_qc_after_processing" \
    "${PROJECT_DIR}/04_star_index" \
    "${PROJECT_DIR}/05_star_result"

echo "------------------------------------------------------------------------------------------"
echo "Directories have created"
echo "------------------------------------------------------------------------------------------"

echo "------------------------------------------------------------------------------------------"
echo "coping raw files"
echo "------------------------------------------------------------------------------------------"
# Copy raw read files into the raw reads directory.
cp -v "${RAW_READS_SOURCE_DIR}"/*.{fasta.gz,fa.gz,fastq.gz,fq.gz} "${RAW_READS_DIR}"/

echo "-------------------------------------------------------------------------------------------"
echo "Setup complete. Raw reads are in ${RAW_READS_DIR}"
echo "-------------------------------------------------------------------------------------------"

echo "Now, copy reference genome and annotation file is starting"

echo "-------------------------------------------------------------------------------------------"
# Copy reference genome and annotation files into the reference genome directory.
cp -v "${REF_GENOME_SOURCE_DIR}"/*.{fa,fasta,fa.gz,fasta.gz,gtf,gff,gff3} "${REFERENCE_GENOME_DIR}"/

echo "Setup complete. Raw reads are in ${RAW_READS_DIR}"
echo "Reference genome files are in ${REFERENCE_GENOME_DIR}"
exit 0

## fasqc report for quality control
