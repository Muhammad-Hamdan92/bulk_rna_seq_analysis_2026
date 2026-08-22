#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

#===========================================================================================
#CONFIGURATION
#===========================================================================================

eval "$(conda shell.bash hook)"

PROJECT_DIR="/home/muhammad/Documents/Transcriptomics/transcriptomics_2026"
RAW_READS_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/00_raw_reads"
REF_GENOME_SOURCE_DIR="/home/muhammad/Documents/Transcriptomics/reference_genome"
RAW_READS_DIR="${PROJECT_DIR}/00_raw_reads"
REFERENCE_GENOME_DIR="${PROJECT_DIR}/00_reference_genome"
SAMPLES="SRR3734796 SRR3734797 SRR3734798 SRR3734816 SRR3734817 SRR3734818"  ##EDIT ACCORDING TO YOUR READS
THREADS=3  #EDIT ACCORDING TO CPU THREADS YOU HAVE
SJDB_OVERHANG=99
RAM_USAGE=8000000000
#===========================================================================================
#CREATING DIRECTORIES
#===========================================================================================

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

#===========================================================================================
#COPING RAW FILES AND RFERENCE GENOME FOR STAR INDEXING
#===========================================================================================

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

#===========================================================================================
# FASTQC AND MULTIQC REPORTS ON RAW READS/BEFORE PROCESSING
#===========================================================================================

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

#===========================================================================================
#FASTP FOR QUALITY CONTROL/TRIMMING BAD READS
#===========================================================================================


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

#===========================================================================================
#FACTQC AND MULTIQC REPORTS ON PROCESSED READS
#===========================================================================================

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

#===========================================================================================
#CREATING STAR INDEX 
#===========================================================================================
cd "${PROJECT_DIR}/04_star_index"
conda activate 03_star


STAR --runMode genomeGenerate \
    --genomeDir "${PROJECT_DIR}/04_star_index" \
    --genomeFastaFiles "${PROJECT_DIR}/00_reference_genome/mouse_reference.fa" \
    --sjdbGTFfile "${PROJECT_DIR}/00_reference_genome/mouse_reference.gtf" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    --runThreadN "${THREADS}" \
    --genomeSAsparseD 3 \
    --genomeSAindexNbases 13 \
    --limitGenomeGenerateRAM "${RAM_USAGE}"
echo "-------------------------------------------------------------------------------------------"
echo "STAR genome index built at ${PROJECT_DIR}/04_star_index"
echo "-------------------------------------------------------------------------------------------"

#===========================================================================================
#ALIGNING THE READS STAR INDEX PER SAMPLESTAR \
echo "------------------------------------------------------------------------------------------"
echo "Aligning reads with STAR"
echo "------------------------------------------------------------------------------------------"

mkdir -p "${PROJECT_DIR}/05_star_result/samtool_qc"

for s in ${SAMPLES}; do
    echo "  -> Aligning ${s}"
    mkdir -p "${PROJECT_DIR}/05_star_result/${s}"

    STAR --runMode alignReads \
        --genomeDir "${PROJECT_DIR}/04_star_index" \
        --readFilesIn \
            "${PROJECT_DIR}/02_processed_reads/${s}.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --quantMode GeneCounts \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${PROJECT_DIR}/05_star_result/${s}/${s}_"

    echo "  -> Indexing + flagstat for ${s}"
    samtools index "${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam"
    samtools flagstat "${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam" \
        > "${PROJECT_DIR}/05_star_result/samtool_qc/${s}.flagstat.txt"
done

echo "STAR alignment complete. Results are in ${PROJECT_DIR}/05_star_result"
#===========================================================================================
exit 0


