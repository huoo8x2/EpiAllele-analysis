#!/usr/bin/env bash

source ~/.bashrc
conda activate cuttag


OUT_DIR=$1
MAPPING_DIR="${OUT_DIR}/mapping"
ALLELE_DIR="${OUT_DIR}/allele"
REGION=$2

mkdir -p "${ALLELE_DIR}/bedgraph"

for DIR in "${MAPPING_DIR}"/*; do
    SAMPLE="$(basename "$DIR")"
    echo "Processing sample: $SAMPLE"

    BAM_FILE="${MAPPING_DIR}/${SAMPLE}/${SAMPLE}_bowtie2.mapped.sorted.bam"
    ## calculate scaling factor 
    seqDepthFile="${BAM_FILE/_bowtie2.mapped.sorted.bam/_bowtie2_mapped.seqDepth}"
    if [ -f "$seqDepthFile" ]; then
        seqDepth=$(cat "$seqDepthFile")
    else
        seqDepth=$(samtools view -F 0x04 "$BAM_FILE" | awk 'END {print int(NR/2)}')
        echo "$seqDepth" > "$seqDepthFile"
    fi

    scale_factor_seqdepth=$(echo "10000000 / $seqDepth" | bc -l)
    echo "Scale factor (seqDepth): $scale_factor_seqdepth"

    ## generate bedgraph
    types=("total" "dba" "c57")
    REGION_DPTOOLS="${REGION//-/:}"

    for type in "${types[@]}"; do
        BAM_PATH="${ALLELE_DIR}/bamfile/${SAMPLE}.${type}.bam"

        ## bin50
        BG_RAW="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin50.raw.bedgraph"
        BG_SCALED="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin50.seqdepth_scaled.bedgraph"
        bamCoverage --binSize 50 --region "$REGION_DPTOOLS" -of bedgraph -b "$BAM_PATH" -o "$BG_RAW"
        bamCoverage --binSize 50 --region "$REGION_DPTOOLS" -of bedgraph --scaleFactor "$scale_factor_seqdepth" -b "$BAM_PATH" -o "$BG_SCALED"

        ## bin100
        BG_RAW="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin100.raw.bedgraph"
        BG_SCALED="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin100.seqdepth_scaled.bedgraph"
        bamCoverage --binSize 100 --region "$REGION_DPTOOLS" -of bedgraph -b "$BAM_PATH" -o "$BG_RAW"
        bamCoverage --binSize 100 --region "$REGION_DPTOOLS" -of bedgraph --scaleFactor "$scale_factor_seqdepth" -b "$BAM_PATH" -o "$BG_SCALED"

        ## bin500
        BG_RAW="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin500.raw.bedgraph"
        BG_SCALED="${ALLELE_DIR}/bedgraph/${SAMPLE}.${type}.bin500.seqdepth_scaled.bedgraph"
        bamCoverage --binSize 500 --region "$REGION_DPTOOLS" -of bedgraph -b "$BAM_PATH" -o "$BG_RAW"
        bamCoverage --binSize 500 --region "$REGION_DPTOOLS" -of bedgraph --scaleFactor "$scale_factor_seqdepth" -b "$BAM_PATH" -o "$BG_SCALED"
    done
done





