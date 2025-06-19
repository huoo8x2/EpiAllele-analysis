#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_path output_path"
    exit 1
fi

input_path=$1
output_path=$2
mkdir -p "$output_path"

for file in "$input_path"/*.call.gz; do
    filename=$(basename "$file" .call.gz)
    zcat $file | awk '$4 == "CG"' > "$output_path/$filename.cg" &
done

wait

