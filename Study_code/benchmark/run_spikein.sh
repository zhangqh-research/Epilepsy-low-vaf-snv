#!/usr/bin/env bash
set -euo pipefail
# Arguments: varfile lesion.bam reference.fa picard.jar output_dir array_index
VARFILE=$1
LESION_BAM=$2
REF=$3
PICARD_JAR=$4
OUTDIR=$5
TASK_ID=$6
SEED=$((2026071500 + TASK_ID))
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)
RAW_BAM=${OUTDIR}/spikein.unsorted.bam
SORTED_BAM=${OUTDIR}/spikein.sorted.bam
TMPDIR=${OUTDIR}/tmp
ADD_STDOUT=${OUTDIR}/addsnv.stdout.log
ADD_STDERR=${OUTDIR}/addsnv.stderr.log
[[ ! -e "$RAW_BAM" && ! -e "$SORTED_BAM" && ! -e "$ADD_STDOUT" && ! -e "$ADD_STDERR" ]]
mkdir -p "$TMPDIR"
addsnv.py \
-v "${VARFILE}" \
-f "${LESION_BAM}" \
-r "${REF}" \
-o "${RAW_BAM}" \
-n 5000 \
-p 8 \
--picardjar "${PICARD_JAR}" \
--aligner mem \
--alignerthreads 1 \
--mindepth 1 \
--maxdepth 100000 \
--minmutreads 1 \
--ignoresnps \
--tagreads \
--maxopen 500 \
--seed "${SEED}" \
--tmpdir "${TMPDIR}" \
> "${ADD_STDOUT}" \
2> "${ADD_STDERR}"
samtools sort -@ 8 -m 4G -o "$SORTED_BAM" "$RAW_BAM"
samtools index -@ 8 "$SORTED_BAM"
samtools quickcheck -v "$SORTED_BAM"
