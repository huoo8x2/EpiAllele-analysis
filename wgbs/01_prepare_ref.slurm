#!/bin/bash

#SBATCH --mem 100G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -J WGBS
#SBATCH -p batch
#SBATCH -o /home/LYQ/hukj_logs/output_%x_%J
#SBATCH -e /home/LYQ/hukj_logs/errput_%x_%J

# Script to prepare reference genome using methylCtools

# Activate the 'methy' conda environment
echo "Activating the 'methylctools' conda environment..."

# Initialize conda
source ~/.bashrc

# Activate the 'methy' environment
conda activate methylctools

# Check if the reference.fa is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 reference.fa input_dir"
    exit 1
fi

# Input reference file (full path)
reference_fa="$1"
input_dir="$2"

# Get the directory and base name of the input file
base=$(basename "$reference_fa" .fa)

# Path to methylCtools executable
methylCtools_path="/data01/huangpinzheng/methylCtools/methylCtools_raw/methylCtools"

# Step 1a: Generate positions and compress
echo "Running methylCtools fapos..."
"$methylCtools_path" fapos -3 "$reference_fa" - | bgzip > "$input_dir/$base.pos.gz"

# Step 1b: Index the positions
echo "Running tabix..."
tabix -s 1 -b 2 -e 2 "$input_dir/$base.pos.gz"

# Step 2a: Convert reference genome
echo "Running methylCtools faconv..."
"$methylCtools_path" faconv "$reference_fa" "$input_dir/$base.conv.fa"

# Step 2b: Index the converted reference with bwa
echo "Running bwa index..."
bwa index -a bwtsw "$input_dir/$base.conv.fa"

echo "All steps completed successfully. Output files are located in: $input_dir"

