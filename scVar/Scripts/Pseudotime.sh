#!/bin/bash

# 初始化变量
all_path=""
sample=""
mutation=""

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
        --mutation)
            mutation="$2"
            shift 2
            ;;
        *)
            echo "未知选项: $1"
            exit 1
            ;;
    esac
done

work_dir=$all_path/$sample
errorlog=$work_dir/logs/${sample}_analysis_stderr.log
outlog=$work_dir/logs/${sample}_analysis_stdout.log

mkdir -p $work_dir/analysis/pseudotime

if [ -d "$work_dir/mapping/run_count_${sample}/outs/soupX_matrix" ]; then
    mkdir -p $work_dir/analysis/pseudotime/matrix_files
    cp $work_dir/mapping/run_count_${sample}/outs/soupX_matrix/* $work_dir/analysis/pseudotime/matrix_files/
    python /codes/pre_pseudotime2.py -i $work_dir/output_data.h5mu -m $mutation -o $work_dir/analysis/pseudotime  >> ${outlog} 2>> ${errorlog}
else
    python /codes/pre_pseudotime.py -i $work_dir/output_data.h5mu -m $mutation -o $work_dir/analysis/pseudotime  >> ${outlog} 2>> ${errorlog}
fi

# python /codes/pre_pseudotime.py -i $work_dir/output_data.h5mu -m $mutation -o $work_dir/analysis/pseudotime  >> ${outlog} 2>> ${errorlog}
Rscript /codes/pseudotime.R $work_dir/analysis/pseudotime $mutation >> ${outlog} 2>> ${errorlog}