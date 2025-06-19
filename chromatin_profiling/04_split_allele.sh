#!/bin/bash

source ~/.bashrc
conda activate cuttag

BAM_FILE=$1
SAMPLE=$2
ALLELE_DIR=$3
REGION=$4
SNP_FILE=$5
SPIKEIN=$6
MULTIPLIER=$7

mkdir -p "${ALLELE_DIR}/bamfile" "${ALLELE_DIR}/bigwig" "${ALLELE_DIR}/stat"

## extract target regions
samtools view -bS -h $BAM_FILE $REGION | samtools sort -o "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 
samtools index "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 

## split allele
python "~/project/EpiAllele/script/chromatin_profiling/04_split_allele.py" --samplename "$SAMPLE" --totalfile "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" --bamfilePath "${ALLELE_DIR}/bamfile/" --statPath "${ALLELE_DIR}/stat/" --snpfile "$SNP_FILE"

# calculate scaling factor, normalization and file format conversion

if [ "$SPIKEIN" = TRUE ]; then
    seqDepthFile="${BAM_FILE/_bowtie2.mapped.sorted.bam/_bowtie2_mapped.seqDepth}"
    spikeInSeqDepthFile="${BAM_FILE/_bowtie2.mapped.sorted.bam/_bowtie2_spikeIn.seqDepth}"
    spikeInSamFile="${BAM_FILE/_bowtie2.mapped.sorted.bam/_bowtie2_spikeIn.sam}"

    # calculate seqDepth from main BAM
    if [ -f "$seqDepthFile" ]; then
        seqDepth=$(cat "$seqDepthFile")
    else
        seqDepth=$(samtools view -F 0x04 "$BAM_FILE" | awk 'END {print int(NR/2)}')
        echo "$seqDepth" > "$seqDepthFile"
    fi

    # calculate spikeInDepth from spike-in SAM
    if [ -f "$spikeInSeqDepthFile" ]; then
        spikeInDepth=$(cat "$spikeInSeqDepthFile")
    else
        spikeInDepth=$(samtools view -F 0x04 "$spikeInSamFile" | awk 'END {print int(NR/2)}')
        echo "$spikeInDepth" > "$spikeInSeqDepthFile"
    fi

    scale_factor_seqdepth=$(echo "$MULTIPLIER / $seqDepth" | bc -l)
    scale_factor_spikein=$(echo "10000 / $spikeInDepth" | bc -l)

    echo "Scale factor (seqDepth): $scale_factor_seqdepth"
    echo "Scale factor (spikeIn): $scale_factor_spikein"

else
    seqDepthFile="${BAM_FILE/_bowtie2.mapped.sorted.bam/_bowtie2_mapped.seqDepth}"

    if [ -f "$seqDepthFile" ]; then
        seqDepth=$(cat "$seqDepthFile")
    else
        seqDepth=$(samtools view -F 0x04 "$BAM_FILE" | awk 'END {print int(NR/2)}')
        echo "$seqDepth" > "$seqDepthFile"
    fi

    scale_factor_seqdepth=$(echo "$MULTIPLIER / $seqDepth" | bc -l)
    echo "Scale factor (seqDepth): $scale_factor_seqdepth"
fi

types=("total" "dba" "c57")
REGION_DPTOOLS="${REGION//-/:}"

for type in "${types[@]}"; do
    BAM_PATH="${ALLELE_DIR}/bamfile/${SAMPLE}.${type}.bam"
    BW_RAW="${ALLELE_DIR}/bigwig/${SAMPLE}.${type}.raw.bw"

    samtools index -@ 15 "$BAM_PATH"
    bamCoverage --binSize 5 --region "$REGION_DPTOOLS" -b "$BAM_PATH" -o "$BW_RAW"

    if [ "$SPIKEIN" = TRUE ]; then
        BW_SCALED_SEQ="${ALLELE_DIR}/bigwig/${SAMPLE}.${type}.seqdepth_scaled.bw"
        BW_SCALED_SPIKEIN="${ALLELE_DIR}/bigwig/${SAMPLE}.${type}.spikein_scaled.bw"

        bamCoverage --binSize 5 --region "$REGION_DPTOOLS" --scaleFactor "$scale_factor_seqdepth" -b "$BAM_PATH" -o "$BW_SCALED_SEQ"
        bamCoverage --binSize 5 --region "$REGION_DPTOOLS" --scaleFactor "$scale_factor_spikein" -b "$BAM_PATH" -o "$BW_SCALED_SPIKEIN"
    else
        BW_SCALED="${ALLELE_DIR}/bigwig/${SAMPLE}.${type}.seqdepth_scaled.bw"
        bamCoverage --binSize 5 --region "$REGION_DPTOOLS" --scaleFactor "$scale_factor_seqdepth" -b "$BAM_PATH" -o "$BW_SCALED"
    fi
done



