#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

#===========================================================================================
#CONFIGURATION
#===========================================================================================

eval "$(conda shell.bash hook)"

PROJECT_DIR="/home/muhammad/Documents/Transcriptomics/transcriptomics_2026"
REFERENCE_GENOME_DIR="${PROJECT_DIR}/00_reference_genome"
STAR_RESULT_DIR="${PROJECT_DIR}/05_star_result"
STRAND_DIR="${PROJECT_DIR}/07_strandedness"
SAMPLES="SRR3734796 SRR3734797 SRR3734798 SRR3734816 SRR3734817 SRR3734818"  ##EDIT ACCORDING TO YOUR READS
THREADS=3  #EDIT ACCORDING TO CPU THREADS YOU HAVE

REFERENCE_GTF="${REFERENCE_GENOME_DIR}/mouse_reference.gtf"
GENEPRED_FILE="${STRAND_DIR}/mouse_reference.genePred"
BED12_FILE="${STRAND_DIR}/mouse_reference.bed12"

#===========================================================================================
#CREATING DIRECTORIES
#===========================================================================================

echo "Create strandedness directories."
echo "---------------------------------------------------------------------------------------"

mkdir -p \
    "${STRAND_DIR}" \
    "${STRAND_DIR}/reports"

echo "------------------------------------------------------------------------------------------"
echo "Directories have created"
echo "------------------------------------------------------------------------------------------"

#===========================================================================================
#CREATING CONDA ENV FOR RSEQC + UCSC GTF-TO-BED12 TOOLS
#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "Creating 07_rseqc conda env (if not already present)"
echo "------------------------------------------------------------------------------------------"

if ! conda env list | grep -q '^07_rseqc '; then
    conda create -n 07_rseqc -y -c bioconda -c conda-forge \
        rseqc \
        ucsc-gtftogenepred \
        ucsc-genepredtobed
else
    echo "07_rseqc env already exists, skipping creation."
fi

echo "------------------------------------------------------------------------------------------"
echo "07_rseqc env is ready"
echo "------------------------------------------------------------------------------------------"

conda activate 07_rseqc

#===========================================================================================
#CONVERTING REFERENCE GTF -> BED12 (RSeQC needs BED12, not GTF)
#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "Converting reference GTF to BED12 for infer_experiment.py"
echo "------------------------------------------------------------------------------------------"

if [[ ! -s "${BED12_FILE}" ]]; then
    gtfToGenePred "${REFERENCE_GTF}" "${GENEPRED_FILE}"
    genePredToBed "${GENEPRED_FILE}" "${BED12_FILE}"
    echo "BED12 written to ${BED12_FILE}"
else
    echo "BED12 already exists at ${BED12_FILE}, skipping conversion."
fi

wc -l "${BED12_FILE}"

#===========================================================================================
#RUNNING infer_experiment.py PER SAMPLE
#===========================================================================================

echo "------------------------------------------------------------------------------------------"
echo "Running infer_experiment.py for ${SAMPLES}"
echo "------------------------------------------------------------------------------------------"

SUMMARY_FILE="${STRAND_DIR}/strandedness_summary.txt"
> "${SUMMARY_FILE}"

for s in ${SAMPLES}; do
    echo "  -> Checking strandedness for ${s}"

    BAM_FILE="${STAR_RESULT_DIR}/${s}/${s}_Aligned.sortedByCoord.out.bam"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "ERROR: BAM not found: ${BAM_FILE}" >&2
        exit 1
    fi

    REPORT_FILE="${STRAND_DIR}/reports/${s}.infer_experiment.txt"

    infer_experiment.py \
        -i "${BAM_FILE}" \
        -r "${BED12_FILE}" \
        > "${REPORT_FILE}" 2>&1

    echo "===== ${s} =====" >> "${SUMMARY_FILE}"
    cat "${REPORT_FILE}" >> "${SUMMARY_FILE}"
    echo "" >> "${SUMMARY_FILE}"
done

echo "------------------------------------------------------------------------------------------"
echo "Strandedness check completed for all samples"
echo "Per-sample reports: ${STRAND_DIR}/reports/"
echo "Combined summary:   ${SUMMARY_FILE}"
echo "------------------------------------------------------------------------------------------"

cat "${SUMMARY_FILE}"

exit 0
