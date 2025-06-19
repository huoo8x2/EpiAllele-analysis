#!/bin/bash
#SBATCH --mem 128G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -J splitAllele
#SBATCH -e /data02/hukaijie/EpiAllele/logs/log/errput_%x_%J
#SBATCH -o /data02/hukaijie/EpiAllele/logs/log/output_%x_%J


source ~/.bashrc
conda activate cuttag

BAM_FILE=$1
SAMPLE=$2
ALLELE_DIR=$3
REGION=$4
SNP_FILE=$5
SPECIES=$6

mkdir -p "${ALLELE_DIR}/bamfile" "${ALLELE_DIR}/stat"

## extract target regions
if [ ! -e "${BAM_FILE}.bai" ]; then
    samtools index -@ 2 "$BAM_FILE"
fi

samtools view -@ 2 -F 0x04 -bS -h $BAM_FILE $REGION | samtools sort -@ 2 -o "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 
samtools index -@ 2 "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 

## split allele
if [ "$SPECIES" = "mouse" ]; then
    python "/data02/hukaijie/EpiAllele/script/cuttag/pipeline/split_allele.py" \
        --samplename "$SAMPLE" \
        --totalfile "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" \
        --bamfilePath "${ALLELE_DIR}/bamfile/" \
        --statPath "${ALLELE_DIR}/stat/" \
        --snpfile "$SNP_FILE"
elif [ "$SPECIES" = "human" ]; then
    python "/data02/hukaijie/EpiAllele/script/RNAseq/pipeline/split_MYH7_allele.py" \
        --samplename "$SAMPLE" \
        --totalfile "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" \
        --bamfilePath "${ALLELE_DIR}/bamfile/" \
        --statPath "${ALLELE_DIR}/stat/" \
        --snpfile "$SNP_FILE"
fi




