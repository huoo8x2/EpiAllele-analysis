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

forward_file=$1
reverse_file=$2
outDir=$3
index=$4
gtf=$5


## alignment + counts quantification
mkdir -p ${outDir}/hisat2
mkdir -p ${outDir}/featureCount
hisat2Path=${outDir}/hisat2
featureCountPath=${outDir}/featureCount
name=$(basename $forward_file .forward_unpaired.fq.gz)

echo "alignment start"
# alignment
hisat2  -p 10 -x $index -1 $forward_file -2 $reverse_file -S ${hisat2Path}/${name}.sam 
echo "******hisat2 sort done !******"
# samtools sorting and generate bam file
samtools sort -O bam -@ 10 -o ${hisat2Path}/${name}.bam ${hisat2Path}/${name}.sam     
echo "******Samtools sort done !******"
# counts quantification
featureCounts  -T 10 -p --countReadPairs -t exon -g gene_id -a $gtf \
	-o ${featureCountPath}/${name}.featurecounts.txt ${hisat2Path}/${name}.bam
echo "******featureCount sort done !******"

echo "alignment done"
