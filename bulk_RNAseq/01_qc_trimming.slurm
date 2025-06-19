#!/bin/bash
#SBATCH --mem 128G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -J runRNAseq
#SBATCH -e /data02/hukaijie/EpiAllele/logs/log/errput_runRNAseq_%J
#SBATCH -o /data02/hukaijie/EpiAllele/logs/log/output_runRNAseq_%J

source ~/.bashrc
conda activate RNAseq

file_1=$1
file_2=$2
outDir=$3


## run fastqc
mkdir -p ${outDir}/fastqc/
fastqcOutput=${outDir}/fastqc/

echo "fastqc start"
fastqc -t 10 -o $fastqcOutput $file_1 $file_2
echo "fastqc done"


## run trimmomatic
filename=$(basename $file_1 _1.fq.gz)
mkdir -p ${outDir}/trimmomatic
trimOutput=${outDir}/trimmomatic

##rna.adapter.fa for trimmomatic 
#>adapter_1
#AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#>adapter_2
#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

echo "Trimmomatic start"
trimmomatic PE $file_1 $file_2 \
	${trimOutput}/${filename}.forward_paired.fq.gz ${trimOutput}/${filename}.forward_unpaired.fq.gz \
	${trimOutput}/${filename}.reverse_paired.fq.gz ${trimOutput}/${filename}.reverse_unpaired.fq.gz \
    ILLUMINACLIP:rna.adapter.fa:2:30:10:8:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:50 -threads 10

echo "Trimmomatic done"

## run fastqc again
mkdir -p ${outDir}/trim_fastqc
fastqcOutput2=${outDir}/trim_fastqc

trim_file_1=${trimOutput}/${filename}.forward_paired.fq.gz
trim_file_2=${trimOutput}/${filename}.reverse_paired.fq.gz

echo "fastqc start"
fastqc -t 10 -o $fastqcOutput2 $trim_file_1 $trim_file_2
echo "fastqc done"
