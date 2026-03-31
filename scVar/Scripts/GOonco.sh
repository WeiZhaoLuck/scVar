#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-31 14:21:51
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-31 14:27:59
 # @FilePath: \scvar_docker\codes\GOonco.sh
### 

### 
# 初始化变量  
all_path=""  
sample=""  
pCutoff=""  
qCutoff=""  

# 解析命令行参数  
while [[ $# -gt 0 ]]; do  
    case "$1" in  
        --path)  
            all_path="$2"  
            shift 2  
            ;;  
        --sample)  
            sample="$2"  
            shift 2  
            ;;  
        --pCutoff)
            pCutoff="$2"
            # 检查是否为合法浮点数
            if ! [[ "$pCutoff" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                echo "错误：--pCutoff 必须是浮点数（如 0.05），实际输入为 '$pCutoff'"
                exit 1
            fi
            shift 2
            ;;
        --qCutoff)
            qCutoff="$2"
            # 检查是否为合法浮点数
            if ! [[ "$qCutoff" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                echo "错误：--qCutoff 必须是浮点数（如 0.05），实际输入为 '$qCutoff'"
                exit 1
            fi
            shift 2
            ;;
        *)  
            echo "未知选项: $1"  
            exit 1  
            ;;  
    esac  
done 

work_path=$all_path/$sample
errorlog=$work_path/logs/${sample}_analysis_stderr.log
outlog=$work_path/logs/${sample}_analysis_stdout.log


mkdir -p $work_path/analysis/pathway
bash /codes/split_maf_bycelltype.sh $work_path/genotype/final.maf $work_path/seurat/output_cell_barcodes.tsv $work_path/analysis/pathway > ${outlog} 2>> ${errorlog} 
Rscript /codes/GO.R $work_path/analysis/pathway $work_path/analysis/pathway $pCutoff $qCutoff > ${outlog} 2>> ${errorlog} 