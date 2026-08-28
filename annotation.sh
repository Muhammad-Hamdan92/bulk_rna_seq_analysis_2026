#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# PsiCLASS vote.gtf -> GRCm39 gene-name annotation
# ============================================================

# ---- Activate your PsiCLASS environment if needed ----
# conda activate 05_psiclass

# ---- Project directories ----
PROJECT_DIR="${HOME/Documents/Transcriptomics/transcriptomics_2026}"

REFERENCE_GENOME_DIR="${PROJECT_DIR}/00_reference_genome}"
PSICLASS_DIR="${PSICLASS_DIR:-${PROJECT_DIR}/06_psiclass_result}"

REFERENCE_GTF="${REFERENCE_GENOME_DIR}/mouse_reference.gtf"
VOTE_GTF="${PSICLASS_DIR}/transcriptomics_2026_vote.gtf"

OUT_DIR="${PSICLASS_DIR}/06.2_vote_with_gene_names"
OUT_GTF="${OUT_DIR}/transcriptomics_2026_vote.annotated.gtf"

mkdir -p "$OUT_DIR"

echo "============================================================"
echo "PsiCLASS vote.gtf -> GRCm39 gene-name annotation"
echo "============================================================"

echo
echo "Reference GTF:"
echo "$REFERENCE_GTF"

echo
echo "PsiCLASS vote GTF:"
echo "$VOTE_GTF"

echo
echo "Output:"
echo "$OUT_GTF"

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [[ ! -s "$REFERENCE_GTF" ]]; then
    echo "ERROR: Reference GTF not found:"
    echo "$REFERENCE_GTF"
    exit 1
fi

if [[ ! -s "$VOTE_GTF" ]]; then
    echo "ERROR: PsiCLASS vote GTF not found:"
    echo "$VOTE_GTF"
    exit 1
fi

# ------------------------------------------------------------
# Check required tools
# ------------------------------------------------------------

for tool in bedtools awk sort; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is not installed."
        echo
        echo "Install bedtools with:"
        echo "conda install -c bioconda bedtools"
        exit 1
    fi
done

# ------------------------------------------------------------
# Temporary files
# ------------------------------------------------------------

TMP_DIR="${OUT_DIR}/tmp_annotation"
mkdir -p "$TMP_DIR"

REF_BED="${TMP_DIR}/grcm39_transcripts.bed"
PSI_BED="${TMP_DIR}/psiclass_transcripts.bed"
BEST_MATCH="${TMP_DIR}/best_matches.tsv"

# ------------------------------------------------------------
# 1. Extract GRCm39 transcript coordinates
#
# We use transcript features because PsiCLASS output contains
# transcript structures.
#
# NOTE: -F'\t' is required here. Without it, awk's default
# whitespace splitting also breaks on the spaces inside the
# attributes column (col 9), so $9 ends up truncated to just
# "gene_id" and the inner split(...,";") loop never matches
# anything -> every record gets silently dropped.
# ------------------------------------------------------------

echo
echo "[1/6] Extracting GRCm39 transcripts..."

awk -F'\t' '
BEGIN { OFS="\t" }

$0 !~ /^#/ && $3 == "transcript" {

    gene_id=""
    gene_name=""
    transcript_id=""

    n=split($9,a,";")

    for(i=1;i<=n;i++){

        gsub(/^[ \t]+|[ \t]+$/,"",a[i])

        if(a[i] ~ /^gene_id /){
            sub(/^gene_id /,"",a[i])
            gsub(/"/,"",a[i])
            gene_id=a[i]
        }

        if(a[i] ~ /^gene_name /){
            sub(/^gene_name /,"",a[i])
            gsub(/"/,"",a[i])
            gene_name=a[i]
        }

        if(a[i] ~ /^transcript_id /){
            sub(/^transcript_id /,"",a[i])
            gsub(/"/,"",a[i])
            transcript_id=a[i]
        }
    }

    if(gene_id != "" && gene_name != "" && transcript_id != ""){
        print $1,$4-1,$5,gene_id,gene_name,transcript_id,$7
    }
}
' "$REFERENCE_GTF" > "$REF_BED"

echo "Reference transcripts:"
wc -l "$REF_BED"

# ------------------------------------------------------------
# 2. Extract PsiCLASS transcript coordinates
# ------------------------------------------------------------

echo
echo "[2/6] Extracting PsiCLASS transcripts..."

awk -F'\t' '
BEGIN { OFS="\t" }

$0 !~ /^#/ && $3 == "transcript" {

    gene_id=""
    transcript_id=""

    n=split($9,a,";")

    for(i=1;i<=n;i++){

        gsub(/^[ \t]+|[ \t]+$/,"",a[i])

        if(a[i] ~ /^gene_id /){
            sub(/^gene_id /,"",a[i])
            gsub(/"/,"",a[i])
            gene_id=a[i]
        }

        if(a[i] ~ /^transcript_id /){
            sub(/^transcript_id /,"",a[i])
            gsub(/"/,"",a[i])
            transcript_id=a[i]
        }
    }

    if(gene_id != "" && transcript_id != ""){
        print $1,$4-1,$5,gene_id,transcript_id,$7
    }
}
' "$VOTE_GTF" > "$PSI_BED"

echo "PsiCLASS transcripts:"
wc -l "$PSI_BED"

# ------------------------------------------------------------
# 3. Find overlaps
#
# -s requires same strand
# -wa -wb reports both records
# ------------------------------------------------------------

echo
echo "[3/6] Finding PsiCLASS/GRCm39 transcript overlaps..."

# NOTE: -s (same-strand) is intentionally omitted. PsiCLASS assembled
# these transcripts without strand information (column 6 is "." for
# every record), so requiring strand agreement against the stranded
# reference GTF would make every intersection fail (0 overlaps).
# We match on genomic position only.

bedtools intersect \
    -wa \
    -wb \
    -a "$PSI_BED" \
    -b "$REF_BED" \
    > "${TMP_DIR}/overlaps.tsv"

echo "Overlaps found:"
wc -l "${TMP_DIR}/overlaps.tsv"

# ------------------------------------------------------------
# 4. Select the best reference match
#
# Score = number of overlapping bases.
#
# For each PsiCLASS transcript, retain the reference transcript
# having the largest genomic overlap.
# ------------------------------------------------------------

echo
echo "[4/6] Selecting best GRCm39 match..."

bedtools intersect \
    -wo \
    -a "$PSI_BED" \
    -b "$REF_BED" \
    > "${TMP_DIR}/overlaps_with_length.tsv"

awk '
BEGIN { OFS="\t" }

{
    psi_gene=$4
    psi_tx=$5

    ref_gene=$10
    ref_name=$11
    ref_tx=$12
    ref_strand=$13

    overlap=$NF

    key=psi_gene SUBSEP psi_tx

    if (!(key in best) || overlap > best[key]) {
        best[key]=overlap
        gene[key]=ref_gene
        name[key]=ref_name
        tx[key]=ref_tx
        strand[key]=ref_strand
    }
}

END {
    for(key in best) {
        split(key,x,SUBSEP)

        print x[1],x[2],gene[key],name[key],tx[key],best[key],strand[key]
    }
}
' "${TMP_DIR}/overlaps_with_length.tsv" \
| sort -k1,1 -k2,2 \
> "$BEST_MATCH"

echo "Best matches:"
wc -l "$BEST_MATCH"

# ------------------------------------------------------------
# 5. Create lookup table and annotate the vote GTF
# ------------------------------------------------------------

echo
echo "[5/6] Creating annotation lookup..."

awk -F'\t' '
BEGIN { OFS="\t" }

FNR==NR {
    key=$1 "|" $2

    gene[key]=$3
    name[key]=$4
    tx[key]=$5
    overlap[key]=$6
    strand[key]=$7

    next
}

$0 ~ /^#/ {
    print
    next
}

{
    # Only modify transcript and exon records.
    if($3 != "transcript" && $3 != "exon"){
        print
        next
    }

    psi_gene=""
    psi_tx=""

    n=split($9,a,";")

    for(i=1;i<=n;i++){

        gsub(/^[ \t]+|[ \t]+$/,"",a[i])

        if(a[i] ~ /^gene_id /){
            sub(/^gene_id /,"",a[i])
            gsub(/"/,"",a[i])
            psi_gene=a[i]
        }

        if(a[i] ~ /^transcript_id /){
            sub(/^transcript_id /,"",a[i])
            gsub(/"/,"",a[i])
            psi_tx=a[i]
        }
    }

    key=psi_gene "|" psi_tx

    # Remove existing gene_name if present.
    gsub(/[ \t]*gene_name "[^"]*";/,"",$9)

    if(key in name){

        # Keep original PsiCLASS gene_id/transcript_id,
        # but add GRCm39 annotation.
        $9=$9 " gene_name \"" name[key] "\";"
        $9=$9 " reference_gene_id \"" gene[key] "\";"
        $9=$9 " reference_transcript_id \"" tx[key] "\";"
        $9=$9 " annotation_overlap_bp \"" overlap[key] "\";"

        # Import strand: PsiCLASS assembled this record with
        # strand "." (unstranded data), so overwrite the actual
        # GTF strand column (field 7) with the reference strand.
        $7=strand[key]

    } else {

        # Truly unmatched transcript.
        $9=$9 " gene_name \"novel_" psi_gene "\";"
        $9=$9 " annotation_status \"novel\";"

    }

    print
}
' "$BEST_MATCH" "$VOTE_GTF" > "$OUT_GTF"

# ------------------------------------------------------------
# 6. Validation
# ------------------------------------------------------------

echo
echo "[6/6] Validating output..."

echo
echo "Output file:"
ls -lh "$OUT_GTF"

echo
echo "First annotated records:"
grep -m 10 'gene_name' "$OUT_GTF"

echo
echo "Known GRCm39 gene names:"
{ grep 'reference_gene_id' "$OUT_GTF" \
    | grep -v 'reference_gene_id ""' \
    | head -10 ; } || true

echo
echo "Novel records:"
grep -c 'annotation_status "novel"' "$OUT_GTF" || true

echo
echo "Known records:"
grep -c 'reference_gene_id' "$OUT_GTF" || true

echo
echo "Unique gene names:"
{ grep -o 'gene_name "[^"]*"' "$OUT_GTF" \
    | sort -u \
    | head -20 ; } || true

# ------------------------------------------------------------
# 7. Import strand + remove no-strand entries
#
# PsiCLASS leaves strand "." for every transcript on unstranded
# data. Step 5 already overwrote column 7 with the matched
# reference strand wherever a match was found. Any record that
# STILL has "." in column 7 here means no reference transcript
# overlapped it at all -> we can't assign strand -> drop it,
# since most downstream tools (featureCounts, IGV, etc.) need a
# real +/- strand.
#
# Set REMOVE_NOVEL=true to ALSO drop genuinely novel/unannotated
# transcripts (annotation_status "novel"), same as chaining a
# `grep -v "novel"` onto the output.
# ------------------------------------------------------------

REMOVE_NOVEL="${REMOVE_NOVEL:-false}"

# Separate output folder (not 06.2) so the unstranded and stranded
# versions sit side by side and are easy to diff/compare directly.
STRAND_OUT_DIR="${PSICLASS_DIR}/06.3_vote_with_strand"
mkdir -p "$STRAND_OUT_DIR"

STRANDED_GTF="${STRAND_OUT_DIR}/transcriptomics_2026_vote.withStrand.gtf"

echo
echo "[7/7] Importing strand, removing no-strand entries..."

if [[ "$REMOVE_NOVEL" == "true" ]]; then
    echo "REMOVE_NOVEL=true -> also dropping novel/unannotated transcripts"
    awk -F'\t' 'BEGIN{OFS="\t"} $0 ~ /^#/ {print; next} $7 != "." && $9 !~ /annotation_status "novel"/ {print}' \
        "$OUT_GTF" > "$STRANDED_GTF"
else
    awk -F'\t' 'BEGIN{OFS="\t"} $0 ~ /^#/ {print; next} $7 != "." {print}' \
        "$OUT_GTF" > "$STRANDED_GTF"
fi

echo
echo "Stranded/filtered GTF:"
wc -l "$STRANDED_GTF"

echo
echo "Dropped (no-strand) records:"
echo $(( $(grep -vc '^#' "$OUT_GTF" || true) - $(grep -vc '^#' "$STRANDED_GTF" || true) ))

echo
echo "============================================================"
echo "DONE"
echo "============================================================"
echo
echo "Final annotated GTF (all records, incl. no-strand/novel):"
echo "$OUT_GTF"
echo
echo "Final stranded + filtered GTF (recommended for downstream use):"
echo "$STRANDED_GTF"
