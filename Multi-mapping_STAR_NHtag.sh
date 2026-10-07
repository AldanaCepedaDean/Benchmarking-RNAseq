#!/bin/bash

# ============================================================
# Split BAM by NH tag and quantify reads with featureCounts
#
# NH:i:1  -> uniquely mapped reads
# NH:i:>1 -> multimapping reads
# ============================================================

set -euo pipefail


# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

# Working/output directory
WORKDIR="/path/to/output_directory"

# Input STAR BAM file
BAM="/path/to/input/sample_Aligned.sortedByCoord.out.bam"

# Gene annotation
ANNOTATION="/path/to/annotation.gtf"

# Number of threads
THREADS=8

# Final output table
FINAL_CSV="${WORKDIR}/NH_unique_vs_multi.csv"


# ------------------------------------------------------------
# CHECK INPUT FILES
# ------------------------------------------------------------

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found:"
    echo "$BAM"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo "ERROR: Annotation file not found:"
    echo "$ANNOTATION"
    exit 1
fi


echo ""
echo "============================================================"
echo " BAM:        $BAM"
echo " ANNOTATION: $ANNOTATION"
echo "============================================================"
echo ""


# ============================================================
# 1. SPLIT BAM FILES BY NH TAG
# ============================================================

echo ">>> Extracting reads with NH=1..."

samtools view -@ "$THREADS" -h "$BAM" | \
awk 'BEGIN{OFS="\t"}
     /^@/ || $0 ~ /(^|\t)NH:i:1(\t|$)/' | \
samtools view -@ "$THREADS" -b -o "${WORKDIR}/NH1.bam" -

samtools index -@ "$THREADS" "${WORKDIR}/NH1.bam"


echo ">>> Extracting reads with NH>1..."

samtools view -@ "$THREADS" -h "$BAM" | \
awk 'BEGIN{OFS="\t"}
     /^@/ || $0 ~ /(^|\t)NH:i:[2-9][0-9]*(\t|$)/' | \
samtools view -@ "$THREADS" -b -o "${WORKDIR}/NHmulti.bam" -

samtools index -@ "$THREADS" "${WORKDIR}/NHmulti.bam"


echo ""
echo ">>> BAM files generated:"
echo "    ${WORKDIR}/NH1.bam"
echo "    ${WORKDIR}/NHmulti.bam"
echo ""


# ============================================================
# 2. FEATURECOUNTS
# ============================================================

echo ">>> Running featureCounts for NH=1 reads..."

featureCounts \
    -T "$THREADS" \
    -a "$ANNOTATION" \
    -p \
    -M \
    --countReadPairs \
    -o "${WORKDIR}/NH1_counts.txt" \
    "${WORKDIR}/NH1.bam"


echo ""
echo ">>> Running featureCounts for NH>1 reads..."

featureCounts \
    -T "$THREADS" \
    -a "$ANNOTATION" \
    -p \
    -M \
    --countReadPairs \
    -o "${WORKDIR}/NHmulti_counts.txt" \
    "${WORKDIR}/NHmulti.bam"


# ============================================================
# 3. EXTRACT COUNT TABLES
# ============================================================

echo ""
echo ">>> Extracting gene-level counts..."

awk 'BEGIN{FS=OFS="\t"} NR>2 {print $1, $7}' \
    "${WORKDIR}/NH1_counts.txt" \
    > "${WORKDIR}/NH1_clean.txt"


awk 'BEGIN{FS=OFS="\t"} NR>2 {print $1, $7}' \
    "${WORKDIR}/NHmulti_counts.txt" \
    > "${WORKDIR}/NHmulti_clean.txt"


# ============================================================
# 4. SORT TABLES BEFORE JOINING
# ============================================================

sort -k1,1 \
    "${WORKDIR}/NH1_clean.txt" \
    > "${WORKDIR}/NH1_clean_sorted.txt"

sort -k1,1 \
    "${WORKDIR}/NHmulti_clean.txt" \
    > "${WORKDIR}/NHmulti_clean_sorted.txt"


# ============================================================
# 5. CALCULATE MULTIMAPPING PERCENTAGE
# ============================================================

echo ">>> Calculating the percentage of multimapping reads..."

join -t $'\t' \
    "${WORKDIR}/NH1_clean_sorted.txt" \
    "${WORKDIR}/NHmulti_clean_sorted.txt" \
| awk 'BEGIN{
        OFS=",";
        print "gene,percent_multimapping"
    }
    {
        gene=$1;
        unique=$2;
        multi=$3;

        total=unique+multi;

        if (total == 0)
            percent=0;
        else
            percent=100*multi/total;

        print gene, percent
    }' \
> "$FINAL_CSV"


# ============================================================
# 6. CLEAN UP INTERMEDIATE FILES
# ============================================================

rm -f \
    "${WORKDIR}/NH1_clean.txt" \
    "${WORKDIR}/NHmulti_clean.txt" \
    "${WORKDIR}/NH1_clean_sorted.txt" \
    "${WORKDIR}/NHmulti_clean_sorted.txt"


# ============================================================
# FINISH
# ============================================================

echo ""
echo "============================================================"
echo "Process completed successfully."
echo "Final table:"
echo "$FINAL_CSV"
echo "============================================================"
