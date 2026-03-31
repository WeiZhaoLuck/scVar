#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-04 20:14:47
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-25 14:35:23
 # @FilePath: \scvar_docker\codes\extract_bam_parallel.sh

### 
bam_path=$1
barcodes_path=$2
out_path=`dirname $bam_path`

mkdir $out_path/bam_tmp
python /codes/multi_extract_parallel3.py $barcodes_path $bam_path $out_path/bam_tmp &&\
find $out_path/bam_tmp/ -name "*.extract.bam" > $out_path/bam_tmp/all_bam_path  &&\
samtools merge -f -b $out_path/bam_tmp/all_bam_path -@ 6 ${out_path}/extracted.bam &&\
samtools index ${out_path}/extracted.bam
rm $out_path/bam_tmp/*.extract.bam

cd $out_path/bam_tmp/
uniques=()
while IFS= read -r -d '' f; do
    base=$(basename "$f" .bam)
    extracted="${base#*_}"
    if [[ ! "$extracted" == *.extract ]]; then
        uniques+=("$extracted")
    fi
done < <(find . -maxdepth 1 -name "*.bam" -print0)

# 去重排序
sorted_uniqs=($(printf "%s\n" "${uniques[@]}" | sort -u))

# 设置最大并发数（根据需求调整）
max_jobs=4
current_jobs=0

for xxx in "${sorted_uniqs[@]}"; do
    (
        # 处理逻辑
        files=(*_"${xxx}".bam)
        
        if [ ${#files[@]} -ge 1 ]; then
            if [ ${#files[@]} -gt 1 ]; then
                echo "Merging files for ${xxx}..."
                samtools merge -@ 4 "${xxx}.bam" "${files[@]}"  # 添加线程参数
            else
                echo "Copying single file for ${xxx}..."
                cp "${files[0]}" "${xxx}.bam"
            fi
            samtools index -@ 4 "${xxx}.bam"  # 添加线程参数
        else
            echo "No files found for ${xxx}, something is wrong." >&2
        fi
        rm "${files[@]}"
    ) &

    # 并发控制
    ((current_jobs++))
    if (( current_jobs % max_jobs == 0 )); then
        wait
    fi
done
wait
