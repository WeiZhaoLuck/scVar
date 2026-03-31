#!/bin/bash  
###  
 # @Author: WeiZhaoLuck zhaow@big.ac.cn  
 # @Date: 2025-03-06 14:29:01  
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-28 09:48:16
 # @FilePath: \scvar_docker\codes\split_maf_bycelltype.sh  
###   

if [ "$#" -ne 3 ]; then  
    echo "用法: $0 <maf_file> <tsv_file> <output_dir>"  
    exit 1  
fi  

maf_file=$1  
tsv_file=$2  
output_dir=$3  
min_count=100  

output_dir_path="$output_dir"  
mkdir -p "$output_dir_path"  

declare -A cell_type_counts  

while IFS=$'\t' read -r barcode cell_type; do  
    ((cell_type_counts["$cell_type"]++))  
done < <(tail -n +2 "$tsv_file")  

valid_cell_types=()  
for cell_type in "${!cell_type_counts[@]}"; do  
    if [ "${cell_type_counts[$cell_type]}" -gt "$min_count" ]; then  
        valid_cell_types+=("$cell_type")  
    fi  
done  

for cell_type in "${valid_cell_types[@]}"; do  
    barcodes=$(awk -v ct="$cell_type" -F'\t' '$2 == ct {print $1}' "$tsv_file")  

    if [ -n "$barcodes" ]; then  
        # 输出maf文件的第一行和符合条件的行  
        awk -v barcodes="$barcodes" -v maf_file="$maf_file" -F'\t' 'BEGIN {  
            # 输出maf_file的第一行  
            getline header < maf_file  
            print header  
            split(barcodes, b, " ")  
            for (i in b) {  
                barcode_set[b[i]] = 1  
            }  
        }  
        {  
            if ($2 in barcode_set) {  
                print  
            }  
        }' "$maf_file" > "$output_dir_path/$cell_type.maf"  

        line_count=$(wc -l < "$output_dir_path/$cell_type.maf")  
        echo "Output file: $output_dir_path/$cell_type.maf, Cell type: $cell_type, Number of lines: $line_count"  
    fi  
done  