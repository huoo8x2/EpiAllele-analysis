#!/bin/bash
#SBATCH --mem 128G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -J wgbs_merge_rep
#SBATCH -o /data02/hukaijie/EpiAllele/logs/log/output_%x_%J
#SBATCH -e /data02/hukaijie/EpiAllele/logs/log/errput_%x_%J

SAMPLE=$1
ALLELE_DIR=$2
BAM_DIR="${ALLELE_DIR}/bamfile"
MERGE_DIR="${ALLELE_DIR}/merge_bamfile"
METHY_DIR="${ALLELE_DIR}/merge_callmethy"
mkdir -p $MERGE_DIR $METHY_DIR
REFERENCE=$3

source ~/.bashrc
methylCtools="/data01/huangpinzheng/methylCtools/methylCtools_raw/methylCtools"

date
echo "----${SAMPLE} START MERGE REP----"


for type in c57 dba total; do
    ## ==============================merge bamfiles===================================
    conda activate cuttag
    echo $type
    samtools merge -@ 2 ${BAM_DIR}/${SAMPLE}*.${type}.bam -o "${MERGE_DIR}/${SAMPLE}.${type}.merge.bam"
    samtools index -@ 2 "${MERGE_DIR}/${SAMPLE}.${type}.merge.bam"

    ## ===============================call methy again================================
    conda activate methyl2
    $methylCtools bcall "$REFERENCE.Myh6.pos.gz" "${MERGE_DIR}/${SAMPLE}.${type}.merge.bam"  - | bgzip > "${METHY_DIR}/${SAMPLE}.${type}.merge.call.gz" 
    tabix -s 1 -b 2 -e 2 "${METHY_DIR}/${SAMPLE}.${type}.merge.call.gz"

    ## ================generate bedgraph (filter CpG region < 10)=====================
    zcat "${METHY_DIR}/${SAMPLE}.${type}.merge.call.gz" | grep CG | awk '$5+$6 > 10 {print $1, $2, $2+1, $5/($5+$6)*100}' > "${METHY_DIR}/${SAMPLE}.${type}.merge.bedgraph"
done

echo "----${SAMPLE} SPLIT MERGE REP----"
date
