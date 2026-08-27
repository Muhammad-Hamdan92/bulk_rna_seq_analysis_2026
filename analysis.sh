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
    "${PROJECT_DIR}/06_psiclass_result"

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

#"${ADDGENENAME}" \
    #"${ANNOTATE}" \
    #"${PSICLASS_GTF_LIST}" \
    #-o "${GENENAME_DIR}"
#===========================================================================================
./annotation.sh     # USING BEDTOOLS FOR ANNOTATION
#===========================================================================================




#===========================================================================================
#STATISTICALS ANALYSIS BEFORE GOING EXPRESSION ANALYSIS
#[GENE COUNTS, TRANSCRIPT COUNTS, COMPARE THE ANNOTATED FILE WITH STRAND ANNOTATED]
#===========================================================================================

#COMPARE THE ANNOTATED FILE WITH STRAND ANNOTATED


# ============================================================
# GTF summary statistics: raw vs gene-named vs stranded
#
# Computes, for each GTF version:
#   - total lines
#   - unique genes   (unique gene_id values)
#   - unique transcripts (unique transcript_id values)
#   - genes with >1 transcript
#
# Same numbers as the manual `cut -d' ' -f2/-f4` pipeline, but
# using grep -o on the gene_id/transcript_id patterns directly,
# so it doesn't depend on attribute position or spacing (which
# differs between vote.gtf, the 06.2 gene-named file, and the
# 06.3 stranded file, since each adds different extra fields).
# ============================================================

# ---- Project directories ----
PROJECT_DIR="${PROJECT_DIR:-$HOME/Documents/Transcriptomics/transcriptomics_2026}"
PSICLASS_DIR="${PROJECT_DIR}/06_psiclass_result"

VOTE_GTF="${PSICLASS_DIR}/transcriptomics_2026_vote.gtf"
ANNOTATED_GTF="${PSICLASS_DIR}/06.2_vote_with_gene_names/transcriptomics_2026_vote.annotated.gtf"
STRANDED_GTF="${PSICLASS_DIR}/06.3_vote_with_strand/transcriptomics_2026_vote.withStrand.gtf"

OUT_DIR="${PSICLASS_DIR}/07_summary_stats"
SUMMARY_TSV="${OUT_DIR}/gtf_summary_stats.tsv"

mkdir -p "$OUT_DIR"

echo "============================================================"
echo "GTF summary statistics"
echo "============================================================"

printf "Version\tTotalLines\tNumGenes\tNumTranscripts\tNumGenes_gt1_txpt\n" > "$SUMMARY_TSV"

# Raw PsiCLASS vote GTF
raw_total_lines=$(wc -l < "$VOTE_GTF")
raw_tx_lines=$(awk -F'\t' '$3=="transcript"' "$VOTE_GTF")
raw_num_genes=$(echo "$raw_tx_lines" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
raw_num_transcripts=$(echo "$raw_tx_lines" | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)
raw_num_multi_tx_genes=$(echo "$raw_tx_lines" | grep -o 'gene_id "[^"]*"; transcript_id "[^"]*"' | sort -u | grep -o 'gene_id "[^"]*"' | sort | uniq -c | awk '$1 > 1' | wc -l)
printf "raw_vote\t%s\t%s\t%s\t%s\n" "$raw_total_lines" "$raw_num_genes" "$raw_num_transcripts" "$raw_num_multi_tx_genes" >> "$SUMMARY_TSV"

# Gene-named GTF
named_total_lines=$(wc -l < "$ANNOTATED_GTF")
named_tx_lines=$(awk -F'\t' '$3=="transcript"' "$ANNOTATED_GTF")
named_num_genes=$(echo "$named_tx_lines" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
named_num_transcripts=$(echo "$named_tx_lines" | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)
named_num_multi_tx_genes=$(echo "$named_tx_lines" | grep -o 'gene_id "[^"]*"; transcript_id "[^"]*"' | sort -u | grep -o 'gene_id "[^"]*"' | sort | uniq -c | awk '$1 > 1' | wc -l)
printf "gene_named_06.2\t%s\t%s\t%s\t%s\n" "$named_total_lines" "$named_num_genes" "$named_num_transcripts" "$named_num_multi_tx_genes" >> "$SUMMARY_TSV"

# Stranded GTF
stranded_total_lines=$(wc -l < "$STRANDED_GTF")
stranded_tx_lines=$(awk -F'\t' '$3=="transcript"' "$STRANDED_GTF")
stranded_num_genes=$(echo "$stranded_tx_lines" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
stranded_num_transcripts=$(echo "$stranded_tx_lines" | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)
stranded_num_multi_tx_genes=$(echo "$stranded_tx_lines" | grep -o 'gene_id "[^"]*"; transcript_id "[^"]*"' | sort -u | grep -o 'gene_id "[^"]*"' | sort | uniq -c | awk '$1 > 1' | wc -l)
printf "stranded_06.3\t%s\t%s\t%s\t%s\n" "$stranded_total_lines" "$stranded_num_genes" "$stranded_num_transcripts" "$stranded_num_multi_tx_genes" >> "$SUMMARY_TSV"

echo
echo "============================================================"
echo "Summary table:"
echo "============================================================"
column -t -s $'\t' "$SUMMARY_TSV"

echo
echo "Saved to: $SUMMARY_TSV"
echo
echo "To open in Excel/LibreOffice, copy this file over and import"
echo "as tab-delimited, or run:"
echo "  cp \"$SUMMARY_TSV\" ~/Desktop/gtf_summary_stats.tsv"

cd /home/muhammad/Documents/Transcriptomics/transcriptomics_2026
wc -l 06_psiclass_result/06.2_vote_with_gene_names/transcriptomics_2026_vote.annotated.gtf
wc -l 06_psiclass_result/06.3_vote_with_strand/transcriptomics_2026_vote.withStrand.gtf


echo "gene name annotation is completed"


echo "PsiCLASS assembly completed. Results are in ${PSICLASS_DIR}"
  
exit 0


