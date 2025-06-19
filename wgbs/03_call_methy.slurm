#!/bin/bash
#SBATCH --mem 100G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -J WGBS
#SBATCH -p batch
#SBATCH -o /home/LYQ/hukj_logs/output_%x_%J
#SBATCH -e /home/LYQ/hukj_logs/errput_%x_%J

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <reference> <reads1> <reads2> <output_folder>"
    exit 1
fi


# Assign input parameters to variablesls
REFERENCE=$1
INPUT_FOLDER=$2
SAMPLE_NAME=$3
OUTPUT_FOLDER=$4


# Activate the 'methy' conda environment
echo "Activating the 'methyl' conda environment..."
conda activate methylctools

# Path to methylCtools executable
methylCtools_path="/data01/huangpinzheng/methylCtools/methylCtools_raw/methylCtools"


# Create the output folder if it doesn't exist
mkdir -p ${OUTPUT_FOLDER}/conv_fq
mkdir -p ${OUTPUT_FOLDER}/conv_bamfile
mkdir -p ${OUTPUT_FOLDER}/bamfile
mkdir -p ${OUTPUT_FOLDER}/callmethy


# Extract sample name from reads1 file path
#SAMPLE_NAME=$(basename "$READS1" | cut -d_ -f1)


# Define file names based on sample name and output folder
READS1="${INPUT_FOLDER}/${SAMPLE_NAME}_1.fq.gz"
READS2="${INPUT_FOLDER}/${SAMPLE_NAME}_2.fq.gz"
READS_CONV="${OUTPUT_FOLDER}/conv_fq/${SAMPLE_NAME}.conv.fq"
BAM_CONV="${OUTPUT_FOLDER}/conv_bamfile/${SAMPLE_NAME}.conv.bam"
BAM_FINAL="${OUTPUT_FOLDER}/bamfile/${SAMPLE_NAME}.bam"
CALL_GZ="${OUTPUT_FOLDER}/callmethy/${SAMPLE_NAME}.call.gz"


# Step 3: Convert FASTQ files using methylCtools
echo "Step 3: Converting FASTQ files using methylCtools"
/data01/huangpinzheng/methylCtools/methylCtools_raw/methylCtools fqconv -1 "$READS1" -2 "$READS2" "$READS_CONV"

# Step 4: Align reads to the converted reference using bwa and convert to BAM
echo "Step 4: Aligning reads and converting to BAM"
bwa mem -p -t 32 -M "$REFERENCE.conv.fa" "$READS_CONV" | samtools view -Sb - > "$BAM_CONV"
rm -f $READS_CONV


# Step 5a: Convert BAM file using methylCtools and sort it
echo "Step 5a: Converting BAM and sorting"
"$methylCtools_path" bconv "$BAM_CONV" - | samtools sort -@ 32 -o "$BAM_FINAL"

# Step 5b: Rename the sorted BAM file to final BAM name and index it
echo "Step 5b: Renaming sorted BAM and indexing"
samtools index -@ 32 "$BAM_FINAL"

# Step 6a: Call methylation using methylCtools and compress the output
echo "Step 6a: Calling methylation and compressing output"
"$methylCtools_path" bcall "$REFERENCE.pos.gz" "$BAM_FINAL" - | bgzip > "$CALL_GZ"

# Step 6b: Index the compressed methylation call file
echo "Step 6b: Indexing the compressed methylation call file"
tabix -s 1 -b 2 -e 2 "$CALL_GZ"

