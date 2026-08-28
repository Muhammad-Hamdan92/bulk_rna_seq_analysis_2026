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
PSICLASS_DIR="${PROJECT_DIR}/06_psiclass_result"
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
    "${PROJECT_DIR}/05_star_result" \
    "${PROJECT_DIR}/06_psiclass_result" \
    "${PROJECT_DIR}/07_featureCounts"

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
#===========================================================================================

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
        --outFileNamePrefix "${PROJECT_DIR}/05_star_result/${s}/${s}_" \
        --limitBAMsortRAM "${RAM_USAGE}"
done

#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "samtool indexing and flagstat is running"
echo "------------------------------------------------------------------------------------------"

echo "------------------------------------------------------------------------------------------"
echo "04_samtool is activated"
echo "------------------------------------------------------------------------------------------"

conda activate 04_samtools

for s in ${SAMPLES}; do
    echo "  -> Indexing + flagstat for ${s}"
    samtools index "${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam"
    samtools flagstat "${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam" \
        > "${PROJECT_DIR}/05_star_result/samtool_qc/${s}.flagstat.txt"
done

echo "------------------------------------------------------------------------------------------"
echo "STAR alignment complete. Results are in ${PROJECT_DIR}/05_star_result"
echo "------------------------------------------------------------------------------------------"
#===========================================================================================

#===========================================================================================
#ASSEMBLE THE TRANSCRIPT BY USING [MULTI ASSEMBLER]PSICALSS
#===========================================================================================

conda activate 05_psiclass

# STAR writes one sorted BAM and index per sample directory.
bam_list=""
for s in ${SAMPLES}; do
    bam_path="${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam"
    if [[ ! -f "${bam_path}" ]]; then
        echo "ERROR: STAR BAM not found: ${bam_path}" >&2
        exit 1
    fi
    bam_list+="${bam_path},"
done
bam_list="${bam_list%,}"

echo "BAM list: $bam_list"

# Run PsiCLASS
psiclass \
    -b "${bam_list}" \
    -o "${PSICLASS_DIR}/transcriptomics_2026" \
    -p "${THREADS}" \
    -s 0

# ADD THE GENE NAMES BY USING PSICLASS
echo "add gene name ........................."
mkdir -p "${PROJECT_DIR}/06_psiclass_result/06.1_with_gene_names"

ADDGENENAME="$(which add-genename)"
ANNOTATE="${REFERENCE_GENOME_DIR}/mouse_reference.gtf"
GENENAME_DIR="${PSICLASS_DIR}/06.1_with_gene_names"
PSICLASS_GTF_LIST="${PSICLASS_DIR}/transcriptomics_2026_gtf.list"

add-gene-name \
    "${ANNOTATE}" \
    "${PSICLASS_GTF_LIST}" \
    -o "${GENENAME_DIR}"
#===========================================================================================
./annotation.sh     # USING BEDTOOLS FOR ANNOTATION
#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "gene name annotation is completed"
echo "------------------------------------------------------------------------------------------"

#===========================================================================================
# GTF QUALITY CHECK BEFORE FEATURECOUNTS
#===========================================================================================

echo "============================================================"
echo "GTF QUALITY CHECK"
echo "============================================================"

GTF="${PSICLASS_DIR}/06.3_vote_with_strand/transcriptomics_2026_vote.withStrand.gtf"

echo
echo "Final GTF:"
echo "$GTF"

echo
echo "------------------------------------------------------------"
echo "1. Total records"
echo "------------------------------------------------------------"

grep -vc '^#' "$GTF"

echo
echo "------------------------------------------------------------"
echo "2. Number of genes"
echo "------------------------------------------------------------"

awk -F'\t' '$3=="transcript"' "$GTF" \
| grep -o 'gene_id "[^"]*"' \
| sort -u \
| wc -l

echo
echo "------------------------------------------------------------"
echo "3. Number of transcripts"
echo "------------------------------------------------------------"

awk -F'\t' '$3=="transcript"' "$GTF" \
| grep -o 'transcript_id "[^"]*"' \
| sort -u \
| wc -l

echo
echo "------------------------------------------------------------"
echo "4. Number of exons"
echo "------------------------------------------------------------"

awk -F'\t' '$3=="exon"' "$GTF" | wc -l

echo
echo "------------------------------------------------------------"
echo "5. Known annotated transcripts"
echo "------------------------------------------------------------"

grep -c 'reference_gene_id' "$GTF"

echo
echo "------------------------------------------------------------"
echo "6. Novel transcripts"
echo "------------------------------------------------------------"

grep -c 'annotation_status "novel"' "$GTF" || true

echo
echo "------------------------------------------------------------"
echo "7. Strand check"
echo "------------------------------------------------------------"

echo "Plus strand:"
awk -F'\t' '$7=="+"' "$GTF" | wc -l

echo "Minus strand:"
awk -F'\t' '$7=="-"' "$GTF" | wc -l

echo "Unknown strand:"
awk -F'\t' '$7=="."' "$GTF" | wc -l

echo
echo "------------------------------------------------------------"
echo "8. Example annotated genes"
echo "------------------------------------------------------------"

awk -F'\t' '$3=="transcript"' "$GTF" \
| grep 'gene_name' \
| head -10

echo
echo "============================================================"
echo "GTF QC COMPLETED"
echo "============================================================"


cd /home/muhammad/Documents/Transcriptomics/transcriptomics_2026
wc -l 06_psiclass_result/06.2_vote_with_gene_names/transcriptomics_2026_vote.annotated.gtf
wc -l 06_psiclass_result/06.3_vote_with_strand/transcriptomics_2026_vote.withStrand.gtf

echo "------------------------------------------------------------------------------------------"
echo "PsiCLASS assembly completed. Results are in ${PSICLASS_DIR}"
echo "------------------------------------------------------------------------------------------"

#===========================================================================================
#FEATURECOUNTS FOR GENE LEVEL QUANTIFICATION
#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "making directories"
echo "------------------------------------------------------------------------------------------"

mkdir -p "${PROJECT_DIR}/07_featureCounts"
cd "${PROJECT_DIR}/07_featureCounts"

mkdir -p "gene_counts"

echo "------------------------------------------------------------------------------------------"
echo "making path for gtf files and out dir of gene counts .txt"
echo "------------------------------------------------------------------------------------------"

GTF_PATH="${PROJECT_DIR}/06_psiclass_result/06.3_vote_with_strand/transcriptomics_2026_vote.withStrand.gtf"
GENECOUNT_DIR="${PROJECT_DIR}/07_featureCounts/gene_counts"

echo "------------------------------------------------------------------------------------------"

echo "------------------------------------------------------------------------------------------"
echo "listing BAM files"
echo "------------------------------------------------------------------------------------------"

BAM_FILES=""

for s in ${SAMPLES}
do
    BAM="${PROJECT_DIR}/05_star_result/${s}/${s}_Aligned.sortedByCoord.out.bam"
    echo "${s}: ${BAM}"
    BAM_FILES="${BAM_FILES} ${BAM}"
done

echo "------------------------------------------------------------------------------------------"
echo "activate 07_subread env"
echo "------------------------------------------------------------------------------------------"

conda activate 07_subread

featureCounts \
    -T 3 \
    -t exon \
    -g gene_id \
    -s 0 \
    -a ${GTF_PATH} \
    -o ${GENECOUNT_DIR}/gene_counts.txt \
    ${BAM_FILES} \
    #"${PROJECT_DIR}/05_star_result/SRR3734796/SRR3734796_Aligned.sortedByCoord.out.bam" \
    #"${PROJECT_DIR}/05_star_result/SRR3734797/SRR3734797_Aligned.sortedByCoord.out.bam" \
    #"${PROJECT_DIR}/05_star_result/SRR3734798/SRR3734798_Aligned.sortedByCoord.out.bam" \
    #"${PROJECT_DIR}/05_star_result/SRR3734816/SRR3734816_Aligned.sortedByCoord.out.bam" \
    #"${PROJECT_DIR}/05_star_result/SRR3734817/SRR3734817_Aligned.sortedByCoord.out.bam" \
    #"${PROJECT_DIR}/05_star_result/SRR3734818/SRR3734818_Aligned.sortedByCoord.out.bam"


#===========================================================================================
# RENAME THE COLUMNS TO ACCESSION NUMBER
#===========================================================================================

cd /home/muhammad/Documents/Transcriptomics/transcriptomics_2026

awk 'BEGIN{OFS="\t"}
NR==2 {
    for(i=7;i<=NF;i++){
        n=split($i,a,"/")
        sub(/_Aligned\.sortedByCoord\.out\.bam$/,"",a[n])
        $i=a[n]
    }
}
{print}
' \
${GENECOUNT_DIR}/gene_counts.txt \
> ${GENECOUNT_DIR}/gene_counts_clean.txt

echo "------------------------------------------------------------------------------------------"
echo "remove the first line"
echo "------------------------------------------------------------------------------------------"
sed -i '1d' \
${GENECOUNT_DIR}/gene_counts_clean.txt

echo "------------------------------------------------------------------------------------------"
echo "gene level qauntification have completed"
echo "------------------------------------------------------------------------------------------"


exit 0


