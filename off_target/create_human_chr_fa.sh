#!/bin/bash
ref_genome="/data01/META/genome/human/GRCh38/GRCh38.p13.genome.fa"
output_dir="/data02/hukaijie/EpiAllele/ref/GRCh38_chr"
mkdir -p "$output_dir"

while read -r chr; do
    echo $chr
    output_file="${output_dir}/GRCh38_chr${chr}.fa"
    echo $output_file
    samtools faidx "$ref_genome" "chr${chr}" > "$output_file"
done < "/data02/hukaijie/EpiAllele/script/offinder/chromosomes.txt"
