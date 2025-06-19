#!/bin/bash
#SBATCH --mem 128G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -p batch
#SBATCH -J split_allele
#SBATCH -o /data02/hukaijie/EpiAllele/logs/log/output_%x_%J
#SBATCH -e /data02/hukaijie/EpiAllele/logs/log/errput_%x_%J

source ~/.bashrc
conda activate cuttag

methylCtools="/data01/huangpinzheng/methylCtools/methylCtools_raw/methylCtools"

BAM_FILE=$1
SAMPLE=$2
ALLELE_DIR=$3
REGION=$4
SNP_FILE=$5
REFERENCE=$6

mkdir -p "${ALLELE_DIR}/bamfile" "${ALLELE_DIR}/callmethy" "${ALLELE_DIR}/stat" "${ALLELE_DIR}/bedgraph"

## extract target regions
samtools view -bS -h $BAM_FILE $REGION | samtools sort -o "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 
samtools index "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" 

## split allele
python "/data02/hukaijie/EpiAllele/script/cuttag/pipeline/split_allele.py" --samplename "$SAMPLE" --totalfile "${ALLELE_DIR}/bamfile/${SAMPLE}.total.bam" --bamfilePath "${ALLELE_DIR}/bamfile/" --statPath "${ALLELE_DIR}/stat/" --snpfile "$SNP_FILE"

## create index
samtools index "${ALLELE_DIR}/bamfile/${SAMPLE}.c57.bam" 
samtools index "${ALLELE_DIR}/bamfile/${SAMPLE}.dba.bam" 

conda activate methyl2

## call methyl
for suffix in c57 dba total; do
    $methylCtools bcall "$REFERENCE.Myh6.pos.gz" "${ALLELE_DIR}/bamfile/${SAMPLE}.${suffix}.bam"  - | bgzip > "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz" 
    tabix -s 1 -b 2 -e 2 "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz"
    zcat "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz" | grep CG | awk '$5+$6 > 10 {print $1, $2, $2+1, $5/($5+$6)*100}' > "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.bedgraph"

    ## bin50 bedgraph
    zcat "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz" | grep CG | \
    awk '$2>=55199278 && $2<=55200278 {
        bin = int($2 / 50);
        start = bin * 50;
        end = start + 50;
        key = $1 "\t" start "\t" end;
        methy[key] += $5;
        total[key] += ($5 + $6);
    }
    END {
        for (k in methy) {
            if (total[k] > 0)
                print k "\t" (methy[k] / total[k]) * 100;
        }
    }' > "${ALLELE_DIR}/bedgraph/${SAMPLE}.${suffix}.bin50.bedgraph"

    ## bin100 bedgraph
    zcat "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz" | grep CG | \
    awk '$2>=55199278 && $2<=55200278 {
        bin = int($2 / 100);
        start = bin * 100;
        end = start + 100;
        key = $1 "\t" start "\t" end;
        methy[key] += $5;
        total[key] += ($5 + $6);
    }
    END {
        for (k in methy) {
            if (total[k] > 0)
                print k "\t" (methy[k] / total[k]) * 100;
        }
    }' > "${ALLELE_DIR}/bedgraph/${SAMPLE}.${suffix}.bin100.bedgraph"

    ## bin500 bedgraph
    zcat "${ALLELE_DIR}/callmethy/${SAMPLE}.${suffix}.call.gz" | grep CG | \
    awk '$2>=55199278 && $2<=55200278 {
        bin = int($2 / 500);
        start = bin * 500;
        end = start + 500;
        key = $1 "\t" start "\t" end;
        methy[key] += $5;
        total[key] += ($5 + $6);
    }
    END {
        for (k in methy) {
            if (total[k] > 0)
                print k "\t" (methy[k] / total[k]) * 100;
        }
    }' > "${ALLELE_DIR}/bedgraph/${SAMPLE}.${suffix}.bin500.bedgraph"
done



