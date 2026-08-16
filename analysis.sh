#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

eval "$(conda shell.bash hook)"

PROJECT_DIR="/home/muhammad/Documents/Transcriptomics/transcriptomics_2026"
RAW_READS_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/00_raw_reads"
REF_GENOME_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/reference_genome"
RAW_READS_DIR="${PROJECT_DIR}/00_raw_reads"
REFERENCE_GENOME_DIR="${PROJECT_DIR}/00_reference_genome"
SAMPLES="SRR3734796 SRR3734797 SRR3734798 SRR3734816 SRR3734817 SRR3734818"  ##EDIT ACCORDING TO YOUR READS
THREADS=3  #EDIT ACCORDING TO CPU THREADS YOU HAVE

echo "Create project directories."
echo "---------------------------------------------------------------------------------------"

# Create project directories.
mkdir -p \
    "${RAW_READS_DIR}" \
    "${REFERENCE_GENOME_DIR}" \
    "${PROJECT_DIR}/01_qc_before_processing" \
    "${PROJECT_DIR}/02_processed_reads" \
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

## fasqc report for quality control
echo "-------------------------------------------------------------------------------------------"
echo "QC on raw reads (before processing)"
echo "-------------------------------------------------------------------------------------------"

cd "${PROJECT_DIR}/01_qc_before_processing"
mkdir -p reports
conda activate 01_short_reads_qc

# Run fastqc for all samples
echo "Processing fastqc for 6 samples..."
fastqc -O reports --extract --svg -t "${THREADS}" "${RAW_READS_DIR}"/*.{fasta.gz,fa.gz,fastq.gz,fq.gz}

echo "fastqc reports complete"
echo "-------------------------------------------------------------------------------------------"
echo "multiqc reporting start"
echo "-------------------------------------------------------------------------------------------"

conda activate 02_short_reads_multiqc

mkdir -p "${PROJECT_DIR}/01_qc_before_processing/multiqc"
multiqc -p -o "${PROJECT_DIR}/01_qc_before_processing/multiqc" "${PROJECT_DIR}/01_qc_before_processing/reports"

echo "-------------------------------------------------------------------------------------------"
echo "QC on raw reads (before processing) completed"
echo "-------------------------------------------------------------------------------------------"

echo "-------------------------------------------------------------------------------------------"
echo "fastp for trimming the bad adapter"
echo "-------------------------------------------------------------------------------------------"

cd "${PROJECT_DIR}/02_processed_reads"
conda activate 01_short_reads_qc

mkdir -p "${PROJECT_DIR}/02_processed_reads/fastp_reports"

for s in ${SAMPLES}; do
    echo "-> Trimming ${s}"
    fastp \
    -i "${PROJECT_DIR}/00_raw_reads/${s}.fastq.gz" \
    -o "${PROJECT_DIR}/02_processed_reads/${s}.fastq.gz" \
    --detect_adapter_for_pe \
    -q 25 \
    -l 36 \
    -h "${PROJECT_DIR}/02_processed_reads/fastp_reports/${s}.fastp.html" \
    -j "${PROJECT_DIR}/02_processed_reads/fastp_reports/${s}.fastp.json" \
    -w "${THREADS}"
done

echo "-------------------------------------------------------------------------------------------"
echo "raw reads have processed"
echo "-------------------------------------------------------------------------------------------"

echo "-------------------------------------------------------------------------------------------"
echo "multiqc and fastqc quality reports on processed reads"
echo "-------------------------------------------------------------------------------------------"

cd "${PROJECT_DIR}/03_qc_after_processing"
mkdir -p "${PROJECT_DIR}/03_qc_after_processing/processed_reads_reports"

conda activate 01_short_reads_qc
echo "Processing fastqc for 6 processed samples..."
fastqc -O processed_reads_reports --extract --svg -t "${THREADS}" \
    "${PROJECT_DIR}/02_processed_reads"/*.fastq.gz

echo "-------------------------------------------------------------------------------------------"
echo "multiqc start"
echo "-------------------------------------------------------------------------------------------"

conda activate 02_short_reads_multiqc
mkdir -p "${PROJECT_DIR}/03_qc_after_processing/processed_multiqc"

multiqc -p -o "${PROJECT_DIR}/03_qc_after_processing/processed_multiqc" \
"${PROJECT_DIR}/03_qc_after_processing/processed_reads_reports"

echo "-------------------------------------------------------------------------------------------"
echo "processed reads quality report completed"
echo "-------------------------------------------------------------------------------------------"

exit 0


