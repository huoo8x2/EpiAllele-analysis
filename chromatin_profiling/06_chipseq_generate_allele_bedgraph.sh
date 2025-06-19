source ~/.bashrc
conda activate cuttag


OUT_DIR=$1
REGION=$2

REGION_DPTOOLS="${REGION//-/:}"
ALLELE_DIR="${OUT_DIR}/allele"
mkdir -p "${ALLELE_DIR}/bedgraph"

suffixes=("c57" "dba" "total")
while read -r target control; do
    for suffix in "${suffixes[@]}"; do
        TARGET_BW="${ALLELE_DIR}/bigwig/${target}.${suffix}.scaled.bw"
        CONTROL_BW="${ALLELE_DIR}/bigwig/${control}.${suffix}.scaled.bw"
        SAMPLE=$(basename "$TARGET_BW" .scaled.bw)
        echo $TARGET_BW
        echo $CONTROL_BW
        echo $SAMPLE

        ## bin50
        BG_FILE="${ALLELE_DIR}/bedgraph/${SAMPLE}.bin50.ratio.bedgraph"
        bigwigCompare --operation ratio --binSize 50 --region "$REGION_DPTOOLS" --fixedStep \
            -b1 "$TARGET_BW" -b2 "$CONTROL_BW" -of bedgraph -o "$BG_FILE"
        ## bin100
        BG_FILE="${ALLELE_DIR}/bedgraph/${SAMPLE}.bin100.ratio.bedgraph"
        bigwigCompare --operation ratio --binSize 100 --region "$REGION_DPTOOLS" --fixedStep \
            -b1 "$TARGET_BW" -b2 "$CONTROL_BW" -of bedgraph -o "$BG_FILE"
        ## bin500
        BG_FILE="${ALLELE_DIR}/bedgraph/${SAMPLE}.bin500.ratio.bedgraph"
        bigwigCompare --operation ratio --binSize 500 --region "$REGION_DPTOOLS" --fixedStep \
            -b1 "$TARGET_BW" -b2 "$CONTROL_BW" -of bedgraph -o "$BG_FILE"
    done
done < "/data02/hukaijie/EpiAllele/data/chipseq/target_control.txt"


