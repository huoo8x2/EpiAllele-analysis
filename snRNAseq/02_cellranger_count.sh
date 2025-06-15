#!/bin/bash
id=$1
transcriptome=$2
fastqs=$3
sample=$1
outdir=$4

echo "START----"
cellranger count --id=$id \
                --transcriptome=$transcriptome \
                --fastqs=$fastqs \
                --sample=$sample \
                --output-dir=$outdir

echo "cellranger done"
