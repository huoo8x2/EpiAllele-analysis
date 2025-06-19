#!/bin/bash
SNP_DIR="/data02/hukaijie/EpiAllele/data/snp_files"
mkdir -p $SNP_DIR

source ~/.bashrc
conda activate cuttag

## mybpc3 snp sort
mybpc3_file="/data02/hukaijie/EpiAllele/data/snp_files/Mybpc3_mpg_snp.csv"
out_mybpc3_file="/data02/hukaijie/EpiAllele/data/snp_files/Mybpc3_mpg_snp_sort.csv"
Rscript "/data02/hukaijie/EpiAllele/script/snp_sort/mpg_snp_file_sort.R" --snp_file $mybpc3_file --out_file $out_mybpc3_file


## mybpc3 snp sort
myh6_file="/data02/hukaijie/EpiAllele/data/snp_files/Myh6_mpg_snp.csv"
out_myh6_file="/data02/hukaijie/EpiAllele/data/snp_files/Myh6_mpg_snp_sort.csv"
Rscript "/data02/hukaijie/EpiAllele/script/snp_sort/mpg_snp_file_sort.R" --snp_file $myh6_file --out_file $out_myh6_file
