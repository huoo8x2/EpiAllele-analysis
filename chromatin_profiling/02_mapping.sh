#!/usr/bin/env bash

source ~/.bashrc
conda activate cuttag

FQ1=$1
FQ2=$2
WD_H=$3
GENOME=$4
SPIKEIN=$5
DUP=$6

NAME=$(basename $FQ1 _val_1.fq.gz)

# date
# echo "----START MAPPING----"


if [[ "$GENOME" == "hg38" ]]; then
    ref="/data01/META/bowtie_ref_grch38/grch38_B2/grch38_1kgmaj_B2"
    chr_len="/data01/wuleilei/LYQ/GRCH38_chr_lens"
elif [[ "$GENOME" == "mm39" ]]; then
    ref="/data01/META/GRCm39/GRCm39_cuttag/GRCm39_cuttag"
    chr_len="/data01/META/GRCm39/GRCm39_cuttag_chr_lens"
else
    echo "Unknown GENOME: $GENOME"
    exit 1
fi

source ~/.bashrc
conda activate cuttag

#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
cores=10
#spikeInRef=/data01/META/GRCm39/Spike_in/Spike_in
#spikeInRef="/data01/META/ncbi-genomes-2023-05-10/Ecoil"
spikeInRef="/data01/META/GRCm39/Spike_in/Spike_in"

if [ ! -d ${WD_H}/${NAME} ];then
mkdir -p ${WD_H}/${NAME}
fi
WD=${WD_H}/${NAME}


## mapping reads to ref genome
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} --rg-id ${NAME} --rg SM:${NAME} --rg LB:CHIPSEQ --rg PL:ILLUMINA -x ${ref} -1 ${FQ1} -2 ${FQ2} -S ${WD}/${NAME}_bowtie2.sam &> ${WD}/${NAME}_bowtie2_map.txt

if [ $? -eq 0 ];then
echo "mapping genome S"
else
echo "mapping genome F" >&2
exit 1
fi

## mapping reads to spike-in
## for chipseq, no spike-in exp during library construction
if [ "$SPIKEIN" = TRUE ]; then
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-overlap --no-dovetail -p ${cores} --rg-id ${NAME}_spin --rg SM:${NAME}_spin --rg LB:CHIPSEQ --rg PL:ILLUMINA -x ${spikeInRef} -1 ${FQ1} -2 ${FQ2} -S ${WD}/${NAME}_bowtie2_spikeIn.sam &> ${WD}/${NAME}_bowtie2_spikeIn.txt

    if [ $? -eq 0 ];then
    echo "mapping Spikein S"
    else
    echo "mapping Spikein F" >&2
    exit 1
    fi
fi

## remove PCR duplicat4es
## this step necessary for chipseq data, not recommended for cuttag data
if [ "$DUP" = TRUE ]; then
    gatk SortSam -I ${WD}/${NAME}_bowtie2.sam -O ${WD}/${NAME}.sorted.sam --SORT_ORDER coordinate 
    if [ $? -eq 0 ];then
    echo "SortSam S"
    else
    echo "SortSam F" >&2
    exit 1
    fi


    gatk MarkDuplicates -I ${WD}/${NAME}.sorted.sam -O ${WD}/${NAME}.sorted.dupMarked.sam -M ${WD}/${NAME}.duplicate_metrics
    if [ $? -eq 0 ];then
    echo "MarkDuplicates S"
    else
    echo "MarkDuplicates F" >&2
    exit 1
    fi

    gatk MarkDuplicates -I ${WD}/${NAME}.sorted.sam -O ${WD}/${NAME}.sorted.rmDup.sam --REMOVE_DUPLICATES true --METRICS_FILE ${WD}/${NAME}.picard.rmDup.txt
    if [ $? -eq 0 ];then
    echo "MarkDuplicates S"
    else
    echo "MarkDuplicates F" >&2
    exit 1
    fi

    MAPPED_SAM=${WD}/${NAME}.sorted.rmDup.sam 
else
    MAPPED_SAM=${WD}/${NAME}.sorted.rmDup.sam 
fi


## get fragment length distribution
samtools view -F 0x04 ${MAPPED_SAM} | awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${WD}/${NAME}_fragmentLen.txt
if [ $? -eq 0 ];then
echo "lengthing S"
else
echo "lengthing F" >&2
exit 1
fi


## filter reads and file format conversion
### filter and convert samfile into bamfile format
samtools view -@ 15 -bS -F 0x04 ${MAPPED_SAM}  >${WD}/${NAME}_bowtie2.mapped.bam
if [ $? -eq 0 ];then
echo "samtobam S"
else
echo "samtobam F" >&2
exit 1
fi

### convert bamfile into bedfile format 
bedtools bamtobed -i ${WD}/${NAME}_bowtie2.mapped.bam -bedpe >${WD}/${NAME}_bowtie2.bed
if [ $? -eq 0 ];then
echo "bamtobed S"
else
echo "bamtobed F" >&2
exit 1
fi

### keep the read pairs that are on the same chromosome and fragment length less than 1000bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${WD}/${NAME}_bowtie2.bed >${WD}/${NAME}_bowtie2.clean.bed

### Only extract the fragment related columns
cut -f 1,2,6 ${WD}/${NAME}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >${WD}/${NAME}_bowtie2.fragments.bed


## Assess replicate reproducibility
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${WD}/${NAME}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ${WD}/${NAME}_bowtie2.fragmentsCount.bin${binLen}.bed

if [ $? -eq 0 ];then
echo "binLen 500 success!"
else
echo "binLen 500 Failure!"
exit 1
fi

## calculate scaling factor, normalization and file format conversion
### for chipseq without spike-in, using total mapped reads as seq depth to calculate scaling factor
if [ "$SPIKEIN" = TRUE ]; then
    seqDepthDouble=`samtools view -F 0x04 ${WD}/${NAME}_bowtie2_spikeIn.sam | wc -l`
    seqDepth=$((seqDepthDouble/2))
    echo $seqDepth >${WD}/${NAME}_bowtie2_spikeIn.seqDepth
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
else
    seqDepthDouble=`samtools view -F 0x04 ${WD}/${NAME}_bowtie2.mapped.bam | wc -l`
    seqDepth=$((seqDepthDouble/2))
    echo $seqDepth >${WD}/${NAME}_bowtie2_mapped.seqDepth
    scale_factor=`echo "10000000 / $seqDepth" | bc -l`
fi

echo "Scaling factor for $histName is: $scale_factor!"

bedtools genomecov -bg -scale $scale_factor -i ${WD}/${NAME}_bowtie2.fragments.bed -g ${chr_len} >${WD}/${NAME}_bowtie2.fragments.normalized.bedgraph
if [ $? -eq 0 ];then
    echo "scale_genomecov success!"
else
    echo "scale genomecov Failure!"
    exit 1
fi


### samfile sort and index
samtools sort -@ 15 -o ${WD}/${NAME}_bowtie2.mapped.sorted.bam  ${WD}/${NAME}_bowtie2.mapped.bam 
if [ $? -eq 0 ];then
echo "samtools sort success!"
else
echo "samtools sort Failure!"
exit 1
fi

samtools index -@ 15 ${WD}/${NAME}_bowtie2.mapped.sorted.bam 
if [ $? -eq 0 ];then
echo "samtools index success!"
else
echo "samtools index Failure!"
exit 1
fi

bamCoverage -b ${WD}/${NAME}_bowtie2.mapped.sorted.bam -o ${WD}/${NAME}_raw.bw
if [ $? -eq 0 ];then
echo "bamCoverage success!"
else
echo "bamCoverage Failure!"
exit 1
fi

bamCoverage --scaleFactor ${scale_factor} -b ${WD}/${NAME}_bowtie2.mapped.sorted.bam -o ${WD}/${NAME}_scaled.bw
if [ $? -eq 0 ];then
echo "scaled bamCoverage success!"
else
echo "scaled bamCoverage Failure!"
exit 1
fi


# echo "----MAPPING DONE----"
# date