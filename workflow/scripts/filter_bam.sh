#!/usr/bin/env bash
set -euo pipefail

INPUT=$1
OUTPUT=$2
BLACKLIST=$3
LAYOUT=$4
TMP_PREFIX=${OUTPUT%.bam}.tmp

FLAGS=(-b -q 30 -F 3852)
if [[ "$LAYOUT" == "PE" ]]; then
  FLAGS+=(-f 2)
fi

samtools view "${FLAGS[@]}" -o "${TMP_PREFIX}.q30.bam" "$INPUT"
if [[ "$LAYOUT" == "PE" ]]; then
  # Keep pairing intact: remove both alignments when either mate intersects a
  # blacklist interval.  Filtering alignments independently can leave an
  # orphan that still carries the proper-pair flag and is then miscounted as a
  # fragment by downstream QC/downsampling.
  samtools sort -n -o "${TMP_PREFIX}.name.bam" "${TMP_PREFIX}.q30.bam"
  bedtools pairtobed -abam "${TMP_PREFIX}.name.bam" -b "$BLACKLIST" -type neither \
    > "${TMP_PREFIX}.noblacklist.bam"
  samtools view -h "${TMP_PREFIX}.noblacklist.bam" \
    | awk 'BEGIN{OFS="\t"} /^@/ || $3 != "chrM"' \
    | samtools sort -o "$OUTPUT"
else
  samtools view -h "${TMP_PREFIX}.q30.bam" \
    | awk 'BEGIN{OFS="\t"} /^@/ || $3 != "chrM"' \
    | samtools view -b -o "${TMP_PREFIX}.nochrM.bam"
  bedtools intersect -v -abam "${TMP_PREFIX}.nochrM.bam" -b "$BLACKLIST" > "$OUTPUT"
fi
samtools index "$OUTPUT"
rm -f "${TMP_PREFIX}.q30.bam" "${TMP_PREFIX}.q30.bam.bai" \
  "${TMP_PREFIX}.name.bam" "${TMP_PREFIX}.noblacklist.bam" \
  "${TMP_PREFIX}.nochrM.bam"
