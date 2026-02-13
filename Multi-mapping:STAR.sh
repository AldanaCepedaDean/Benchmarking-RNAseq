#!/bin/bash

# ==========================
# Script: STAR_MultiMapping.sh
# Description: This script run STAR + samtools + featureCounts 
#              and calculate the percentage of multi-mapping per gene.
# ==========================

set -euo pipefail

# ---------- Función de ayuda ----------
usage() {
    echo "Uso: $0 -g GENOME_DIR -a ANNOTATION -1 READ1 -2 READ2 -o OUTPUT_DIR [-t THREADS]"
    echo
    echo "Obligatory parameters:"
    echo "  -g   Directory with the STAR genome index"
    echo "  -a   Annotation file (GTF/GFF)"
    echo "  -1   Read 1 FASTQ file"
    echo "  -2   Read 2 FASTQ file"
    echo "  -o   Directory for the output"
    echo
    echo "Optional parameters:"
    echo "  -t   Threads (default=8)"
    echo "  -h   Help"
    exit 1
}

# ---------- Valores por defecto ----------
THREADS=8
READ1=""
READ2=""

# ---------- Parsear opciones ----------
while getopts ":g:a:1:2:o:t:h" opt; do
  case $opt in
    g) GENOME_DIR="$OPTARG" ;;
    a) ANNOTATION="$OPTARG" ;;
    1) READ1="$OPTARG" ;;
    2) READ2="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    \?) echo " Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# ---------- Verificar argumentos ----------
if [[ -z "${GENOME_DIR:-}" || -z "${ANNOTATION:-}" || -z "${READ1:-}" || -z "${READ2:-}" || -z "${OUT_DIR:-}" ]]; then
    echo " Error: missing some obligatory parameters."
    usage
fi

mkdir -p "$OUT_DIR"
FINAL_CSV="${OUT_DIR}/multimapping_percentages.csv"
echo "GeneID,Multi_Percent" > "$FINAL_CSV"

# ---------- Nombre de muestra ----------
sample=$(basename "$READ1" | sed 's/_1.*//')
prefix="${OUT_DIR}/${sample}_STAR"

echo "🔹 Processing sample: $sample"

# ---------- STAR ----------
STAR \
    --runThreadN $THREADS \
    --readFilesIn "$READ1" "$READ2" \
    --genomeDir "$GENOME_DIR" \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 50 \
    --quantMode TranscriptomeSAM \
    --outSJfilterReads Unique \
    --outSJfilterOverhangMin 150 150 150 150 \
    --outFilterType BySJout \
    --outFilterMismatchNoverReadLmax 0.01 \
    --alignEndsType EndToEnd \
    --outSAMunmapped Within \
    --outFileNamePrefix "${prefix}_"

BAM="${prefix}_Aligned.sortedByCoord.out.bam"

# ---------- SAMTOOLS ----------
samtools view -b -F 0x100 "$BAM" > "${prefix}_unique.bam"
samtools view -b -f 0x100 "$BAM" > "${prefix}_multi.bam"

# ---------- FEATURECOUNTS ----------
featureCounts -T $THREADS -a "$ANNOTATION" -p -o "${prefix}_unique_counts.txt" "${prefix}_unique.bam" 
featureCounts -T $THREADS -a "$ANNOTATION" -p -M --countReadPairs --fraction --fracOverlap 0.1 -o "${prefix}_multi_counts.txt" "${prefix}_multi.bam" 

# ---------- Limpieza y cálculo de porcentaje ----------
awk 'BEGIN{FS=OFS="\t"} NR>2 {print $1, $7}' "${prefix}_unique_counts.txt" > "${prefix}_unique_clean.txt"
awk 'BEGIN{FS=OFS="\t"} NR>2 {print $1, $7}' "${prefix}_multi_counts.txt" > "${prefix}_multi_clean.txt"

join -t $'\t' "${prefix}_unique_clean.txt" "${prefix}_multi_clean.txt" | \
awk 'BEGIN{OFS=","}
    {
        gene=$1; unique=$2; multi=$3;
        total=unique+multi;
        if (total == 0) percent=0;
        else percent=100*multi/total;
        print gene, percent
    }' >> "$FINAL_CSV"

echo "Process finished. Final table: $FINAL_CSV"
