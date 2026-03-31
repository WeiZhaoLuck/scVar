#!/bin/bash
###
# @Author: WeiZhaoLuck zhaow@big.ac.cn
# @Date: 2025-03-04 19:11:21
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2026-01-08 11:50:33
# @FilePath: \scvar_docker\codes\Input.sh
###

set -euo pipefail

# Check input parameters
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <sample_name>"
    exit 1
fi

# Get the input file path
input_dir="$1"
output_dir="$2"
sample_name="$3"

# Find input file (.h5ad or .rds)
input_file="$(find "$input_dir" -type f \( -name "*.h5ad" -o -name "*.rds" \) | head -n 1 || true)"

if [ -z "$input_file" ]; then
    echo "Error: No .h5ad or .rds file found in directory: $input_dir"
    exit 1
fi

echo "Input file: $input_file"
echo "Output dir: $output_dir"
echo "Sample name: $sample_name"

# Run preprocessing (must succeed)
case "$input_file" in
    *.rds)
        echo "Detected RDS file, running R script..."
        Rscript /codes/input.R "$input_file" "$output_dir"
        ;;
    *.h5ad)
        echo "Detected H5AD file, running Python script..."
        python /codes/input.py "$input_file" "$output_dir"
        # 注意：这里不再强行 conda deactivate（除非你确实在本脚本里 activate 过）
        ;;
    *)
        echo "Error: Unsupported file format. Please provide a .rds or .h5ad file."
        exit 1
        ;;
esac

echo "Preprocessing finished successfully. Continue to move BAM..."

# Prepare target directory
basepath="$(dirname "$output_dir")"
target_dir="$basepath/mapping/run_count_${sample_name}/outs"
mkdir -p "$target_dir"

# Find BAM file and check if exists
input_bam="$(find "$input_dir" -type f -name "*.bam" | head -n 1 || true)"
if [ -z "$input_bam" ]; then
    echo "Error: No BAM file found in directory $input_dir"
    exit 1
fi

# Check BAM index exists (optional but safer)
if [ ! -f "${input_bam}.bai" ]; then
    echo "Error: BAM index not found: ${input_bam}.bai"
    exit 1
fi

# Move BAM and BAI
mv "$input_bam" "$target_dir/possorted_genome_bam.bam"
mv "${input_bam}.bai" "$target_dir/possorted_genome_bam.bam.bai"

echo "Moved BAM and BAI to: $target_dir"
