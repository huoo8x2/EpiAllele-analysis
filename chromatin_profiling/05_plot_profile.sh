#!/usr/bin/env bash

source ~/.bashrc
conda activate cuttag

DATA_DIR=$1
RESULT_DIR=$2
GROUP=$3
TYPE=$4
START=$5
END=$6
CORES=$7

TSS=55204384

if [[ $START -lt $TSS ]]; then
    A=$(( ( (TSS - START + 9) / 10 ) * 10 )) 
    B=$(( ( (END - TSS + 9) / 10 ) * 10 ))  
else
    A=$(( ( (END - TSS + 9) / 10 ) * 10 ))  
    B=$(( ( (TSS - START + 9) / 10 ) * 10 )) 
fi

A_kb=$(echo "scale=1; ($A / 1000 + 0.05) / 1" | bc)
B_kb=$(echo "scale=1; ($B / 1000 + 0.05) / 1" | bc)

echo "Computed parameters: -a $A -b $B"


if [[ "$TYPE" == both.* ]]; then
    suffix="${TYPE#both.}"  # remove "both."
    INPUT_FILES="${DATA_DIR}/${GROUP}_*.c57.${suffix}.bw ${DATA_DIR}/${GROUP}_*.dba.${suffix}.bw"
    PLOT_COLORS="#8ECFC9 #FFBE7A #FA7F6F #82B0D2"
    PLOT_NAME="${GROUP}_${suffix}_TSS-${A_kb}kb_+${B_kb}kb"
else
    INPUT_FILES="${DATA_DIR}/${GROUP}_*.${TYPE}.bw"
    PLOT_COLORS="#94AED3 #F29391"
    PLOT_NAME="${GROUP}.${TYPE}_${A_kb}kb_TSS_${B_kb}kb"
fi

MATRIX_FILE="${RESULT_DIR}/${PLOT_NAME}.mat.gz"

echo $MATRIX_FILE
echo $PLOT_COLORS

computeMatrix reference-point -S $INPUT_FILES \
    -R "/data02/hukaijie/EpiAllele/ref/gene_bed/myh6.bed" \
    --referencePoint TSS -a $A -b $B \
    --skipZeros -o "$MATRIX_FILE" -p "$CORES"

plotProfile -m "$MATRIX_FILE" \
    -out "${RESULT_DIR}/${PLOT_NAME}_profile.pdf" \
    --numPlotsPerRow 2 \
    --perGroup \
    --yMin 0 \
    --color $PLOT_COLORS \
    --plotTitle "$PLOT_NAME"