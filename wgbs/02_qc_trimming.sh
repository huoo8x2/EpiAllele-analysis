#!/usr/bin/env bash

source ~/.bashrc
conda activate cuttag

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <file1> <file2> <output_directory>"
    exit 1
fi


file_1=$1
file_2=$2
outDir=$3



## run fastqc
mkdir -p ${outDir}/fastqc
fastqcOutput=${outDir}/fastqc

echo "fastqc start"
fastqc -t 10 -o $fastqcOutput $file_1 $file_2
echo "fastqc done"


## run trim_galore
mkdir -p ${outDir}/trim_fq
trimOutput=${outDir}/trim_fq
samplename=$(basename $file_1 _1.fq.gz)


echo "trim_galore start"
trim_galore -q 30 --cores 4 --fastqc --length 40 --paired \
    -o "$trimOutput" --basename "$samplename" "$file_1" "$file_2"
echo "trim_galore done"s