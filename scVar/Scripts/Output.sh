#!/bin/bash 
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-28 14:03:10
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-31 10:30:04
 # @FilePath: \scvar_docker\codes\Output.sh
### 

# 目标文件夹（可替换为你的路径）
TARGET_DIR=$1
RESULT_DIR=$2
THREAD=$(( $3 ))

# # 检查是否存在 .h5ad 文件
# if compgen -G "${TARGET_DIR}/*.h5ad" > /dev/null; then
#     echo "Found .h5ad files. Running command for h5ad..."
#     # 替换为你的命令（例如读取 h5ad 文件的 Python 脚本）
#     nohup python /codes/output.py -r ${RESULT_DIR}/genotype/genotype_ref_cell_type.tsv -a ${RESULT_DIR}/genotype/genotype_alt_cell_type.tsv -m ${RESULT_DIR}/genotype/all_final.tsv -d /data/combined_samples_final.h5ad -o ${RESULT_DIR} -t 4
#     exit 0  # 执行后退出脚本
# fi

if h5ad_files=$(compgen -G "${TARGET_DIR}/*.h5ad") && [ -n "$h5ad_files" ]; then
    echo "Found .h5ad files: $h5ad_files"
    # 取第一个文件路径（空格分隔的列表）
    first_h5ad=$(echo "$h5ad_files" | head -n 1)
    echo "Processing first file: $first_h5ad"
    python /codes/output.py -r ${RESULT_DIR}/genotype/genotype_ref_cell_type.tsv -a ${RESULT_DIR}/genotype/genotype_alt_cell_type.tsv -m ${RESULT_DIR}/genotype/all_final.tsv -d $first_h5ad -o ${RESULT_DIR} -t $THREAD
    # exit 0

else
    python /codes/output.py -r ${RESULT_DIR}/genotype/genotype_ref_cell_type.tsv -a ${RESULT_DIR}/genotype/genotype_alt_cell_type.tsv -m ${RESULT_DIR}/genotype/all_final.tsv -d ${RESULT_DIR}/seurat/rna.h5ad -o ${RESULT_DIR} -t $THREAD
    # exit 0
fi

# if rds_files=$(compgen -G "${TARGET_DIR}/*.rds") && [ -n "$h5ad_files" ]; then
#     echo "Found .rds files: $rds_files"
#     # 取第一个文件路径（空格分隔的列表）
#     first_rds=$(echo "$rds_files" | head -n 1)
#     echo "Processing first file: $first_rds"

#     nohup python /codes/output.py -r ${RESULT_DIR}/genotype/genotype_ref_cell_type.tsv -a ${RESULT_DIR}/genotype/genotype_alt_cell_type.tsv -m ${RESULT_DIR}/genotype/all_final.tsv -d $first_h5ad -o ${RESULT_DIR} -t 4
#     exit 0
# fi

# # 检查是否存在 .rds 文件
# if compgen -G "${TARGET_DIR}/*.rds" > /dev/null; then
#     echo "Found .rds files. Running command for RDS..."
#     # 替换为你的命令（例如读取 RDS 文件的 R 脚本）
#     Rscript process_rds.R "$TARGET_DIR"
#     exit 0
# fi

# # 检查是否存在 .fastq 文件
# if compgen -G "${TARGET_DIR}/*.fastq" > /dev/null || compgen -G "${TARGET_DIR}/*.fastq.gz" > /dev/null; then
#     echo "Found FASTQ files. Running command for FASTQ..."
#     # 替换为你的命令（例如运行 FastQC）
#     fastqc "${TARGET_DIR}"/*.fastq*
#     exit 0
# fi

# 如果没有任何匹配的文件
# echo "No relevant files found (.h5ad, .rds, .fastq)."
# mv ${RESULT_DIR}/genotype ${RESULT_DIR}/Genotype
# mv ${RESULT_DIR}/seurat ${RESULT_DIR}/Seurat
# mv ${RESULT_DIR}/mapping ${RESULT_DIR}/Mapping
# mv ${RESULT_DIR}/mutation ${RESULT_DIR}/Mutation
# exit 1
