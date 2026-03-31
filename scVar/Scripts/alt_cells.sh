#!/bin/bash

input_file=$1

if [ ! -f "$input_file" ]; then
    echo "Error: Input file not found!"
    exit 1
fi

declare -A map

while IFS=$'\t' read -r col1 col2 col3 col4 barcode _; do
    if [[ -n ${barcode} ]]; then
        map["${col1}_${col2}_${col3}_${col4}"]="${map["${col1}_${col2}_${col3}_${col4}"]:+${map["${col1}_${col2}_${col3}_${col4}"]},}${barcode}"
    fi
done < "$input_file"

for key in "${!map[@]}"; do
    IFS='_' read -r -a key_parts <<< "$key"
    echo -e "${key_parts[0]}\t${key_parts[1]}\t${key_parts[2]}\t${key_parts[3]}\t${map[$key]}"
done