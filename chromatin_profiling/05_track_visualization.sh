#!/bin/bash

OUT_DIR=$1
ALLELE_DIR="${OUT_DIR}/allele"
GENOME=$2
TYPE=$3
SPIKEIN=$4

source ~/.bashrc
conda activate cuttag

dirs=( $(ls "$OUT_DIR/mapping") )
# extract groups
declare -A group_set
for dir in "${dirs[@]}"; do
    group=$(sed 's/_rep[0-9]\+$//' <<< "$dir")
    group_set["$group"]=1
done
sorted_groups=($(echo "${!group_set[@]}" | tr ' ' '\n' | sort | tr '\n' ' '))
MERGE_GROUPS=$(IFS=,; echo "${sorted_groups[*]}")
COMPARE_GROUPS=("${sorted_groups[@]}")

# extract histone part
declare -A factor_set
for group in "${sorted_groups[@]}"; do
    factor=${group%_*}
    factor_set["$factor"]=1
done
sorted_factors=($(echo "${!factor_set[@]}" | tr ' ' '\n' | sort | tr '\n' ' '))
PROFILE_GROUPS=("${sorted_factors[@]}")
REGION="chr14:55169378:55214384"
PROFILE_REGIONS=(
    "55194524_55208581 55194524 55208581"
    "55202000_55205000 55202000 55205000"
)

## suffix
if [ "$SPIKEIN" = TRUE ]; then
    types=("total.raw" "c57.raw" "dba.raw" "total.seqdepth_scaled" "c57.seqdepth_scaled" "dba.seqdepth_scaled" "total.spikein_scaled" "c57.spikein_scaled" "dba.spikein_scaled")
    plotProfile_types=("c57.seqdepth_scaled" "dba.seqdepth_scaled" "both.seqdepth_scaled" "c57.raw" "dba.raw" "both.raw" "c57.spikein_scaled" "dba.spikein_scaled" "both.spikein_scaled")
    pattern_sets=(
    "c57.spikein_scaled.bw,c57.seqdepth_scaled.bw,c57.raw.bw"
    "dba.spikein_scaled.bw,dba.seqdepth_scaled.bw,dba.raw.bw"
    "total.spikein_scaled.bw,total.seqdepth_scaled.bw,total.raw.bw"
    )
else
    types=("total.raw" "c57.raw" "dba.raw" "total.seqdepth_scaled" "c57.seqdepth_scaled" "dba.seqdepth_scaled")
    plotProfile_types=("c57.seqdepth_scaled" "dba.seqdepth_scaled" "both.seqdepth_scaled" "c57.raw" "dba.raw" "both.raw")
    pattern_sets=(
    "c57.seqdepth_scaled.bw,c57.raw.bw"
    "dba.seqdepth_scaled.bw,dba.raw.bw"
    "total.seqdepth_scaled.bw,total.raw.bw"
    )
fi

## merge replicates (bigwig)
BW_PATH="${ALLELE_DIR}/bigwig"
prefixes=$(ls ${BW_PATH} | sed -E 's/_rep[123].*//' | uniq)
MERGE_PATH="${ALLELE_DIR}/merge"
mkdir -p $MERGE_PATH

for group in $prefixes; do
    for type in "${types[@]}"; do
        bigwigAverage -b ${BW_PATH}/${group}_rep*.${type}.bw --binSize 5 --region $REGION -o ${MERGE_PATH}/${group}.${type}.bw
    done
done

## visualization
### merged
MERGE_DIR="${ALLELE_DIR}/merge"
PLOT_DIR="${ALLELE_DIR}/merge_trackplot"
mkdir -p $PLOT_DIR
for i in "${!types[@]}"; do
    PATTERN="${types[$i]}"
    trackplot_name="${PATTERN}.trackplot"
    echo "$PATTERN -> $trackplot_name"
    Rscript "/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/05_merge_compare_trackplot.R" --data_path $MERGE_DIR --results_path $PLOT_DIR --groups $MERGE_GROUPS --pattern $PATTERN --genome $GENOME --trackplot_name $trackplot_name --type $TYPE
done

### individual (compare replicates)
PLOT_DIR="${ALLELE_DIR}/rep_compare_trackplot"
mkdir -p $PLOT_DIR

for group in "${COMPARE_GROUPS[@]}"; do
    for patterns in "${pattern_sets[@]}"; do
        IFS=',' read -r file1 file2 <<< "$patterns"
        prefix="${file1%%.*}"
        trackplot_name="${group}.${prefix}.trackplot"
        echo "$group -> $trackplot_name with patterns: $patterns"
        Rscript "/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/05_rep_compare_trackplot.R" \
            --merge_path "$MERGE_DIR" --rep_path "$BW_PATH" --results_path "$PLOT_DIR" \
            --group "$group" --patterns "$patterns" --genome "$GENOME" --trackplot_name "$trackplot_name"
    done
done

## plotprofile
PROFILE_DIR="${ALLELE_DIR}/plotProfile"
mkdir -p $PROFILE_DIR
CORES=4

for group in "${PROFILE_GROUPS[@]}"; do
    for region_info in "${PROFILE_REGIONS[@]}"; do
        read -r region_name start end <<< "$region_info"

        for type in "${plotProfile_types[@]}"; do
            echo "Submitting: $group - $type - $region_name"
            bash "/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/05_plot_profile.sh" \
                "$MERGE_DIR" "$PROFILE_DIR" "$group" "$type" "$start" "$end" "$CORES"
        done
    done
done