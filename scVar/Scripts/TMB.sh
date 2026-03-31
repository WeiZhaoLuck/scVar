#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-05 16:52:26
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-04-08 15:04:24
### 


results_path=$1
sample=$2
all_bam=$results_path/mapping/run_count_$sample/outs/extracted.bam
barcode_path=$results_path/seurat/output_cell_barcodes.tsv
maf_path=$results_path/mutation/results/variants/filter.maf
tsv_path=$results_path/genotype/all_final.tsv
out_path=$results_path/analysis/TMB/TMB.txt
# Region=BT1301
# cell_type="Alveolar cell"  # 定义变量
result_dir=$results_path/analysis/TMB

directory=$(dirname "$all_bam")

# if [ ! -d "$directory/bam_tmp" ]; then
#     # 如果不存在，则创建文件夹
#     mkdir "$directory/bam_tmp"
#     echo "Directory created: $directory/bam_tmp"
# else
#     echo "Directory already exists: $directory/bam_tmp"
# fi

if [ ! -d "$result_dir" ]; then
    # 如果不存在，则创建文件夹
    mkdir "$result_dir"
    echo "Directory created: $result_dir"
else
    echo "Directory already exists: $result_dir"
fi
if [ -f "$out_path" ]; then  
    > "$out_path"  # 清空文件  
else  
    touch "$out_path"  # 创建文件  
fi  
cd $directory/bam_tmp
barcode_dir=$(dirname "$barcode_path")
# bamtools split -in "$all_bam" -reference &&\
# mv $directory/extracted.REF_* $directory/bam_tmp/ &&\
#
## 处理过滤后的条形码
##awk -F, '{count[$NF]++} END {for (str in count) if (count[str] > 200) print str}' "$barcode_path" | \
##while read -r str; do
##  grep ",$str$" "$barcode_path"
##done > "$barcode_dir/barcodes_type_filtered.csv" &&\
##
##{
##  head -n 1 $barcode_path
##  awk -F, '{count[$NF]++} END {for (str in count) if (count[str] > 200) print str}' $barcode_path | \
##   while read -r str; do grep ",$str$" $barcode_path; done
##} > "$barcode_dir/barcodes_type_filtered.csv" &&\
#
#
#{
#  head -n 1 "$barcode_path"  # 保留第一行
#  awk -F'\t' -v type="$cell_type" '$2 == type {print}' "$barcode_path"  # 仅保留第二列等于变量值的行
#} > "$barcode_dir/${Region}_filtered.tsv"

{
 head -n 1 $barcode_path
 awk -F'\t' '{count[$NF]++} END {for (str in count) if (count[str] > 100) print str}' $barcode_path | \
  while read -r str; do grep "$str$" $barcode_path; done
} > "$barcode_dir/barcodes_type_filtered.csv" &&\

# python /codes/extract_bam_bycelltype.py "$barcode_dir/${Region}_filtered.tsv" "$directory/bam_tmp" &&\

mapfile -t cell_types < <(awk -F'\t' 'NR>1 {print $NF}' "$barcode_dir/barcodes_type_filtered.csv" | sort -u)

maf_dir=$(dirname "$maf_path")
tsv_dir=$(dirname "$tsv_path")

# 定义处理每个细胞类型的函数
process_cell_type() {
  local cell_type=${1// /_}
#   samtools merge -@ 6 "${cell_type}.bam" *_${cell_type}.extract.bam && samtools index "${cell_type}.bam"
#  rm *_${cell_type}.extract.bam # 清理中间处理文件

  grep "${cell_type}" "$tsv_path" > "$tsv_dir/${cell_type}.tsv"

  awk 'NR==FNR{a[$1]; next} $1 in a' "$tsv_dir/${cell_type}.tsv" "$maf_path" > "$maf_dir/${cell_type}.maf"

  local bam_path="${directory}/bam_tmp/${cell_type}.bam"
  local sample="${cell_type}"
  local vcf_path="$maf_dir/${cell_type}.maf"

  name=$(basename "$bam_path" .bam)
  cd "$directory/bam_tmp"

  /opt/mosdepth -t 6 "$name" "$bam_path"

  gunzip -c "$name.per-base.bed.gz" > "$name.per-base.bed" # 避免不必要解压再压缩

  total_length=0
  while IFS=$'\t' read -r chromosome start end reads; do
      if [[ $chromosome =~ ^[0-9]+$ && $chromosome -le 22 || $chromosome == "X" || $chromosome == "Y" ]] && [ $reads -gt 0 ]; then
          region_length=$((end - start))
          total_length=$((total_length + region_length))
      fi
  done < "$name.per-base.bed"

  echo "Total length of $sample qualifying regions: $total_length" >> "${out_path}"

  count=$(awk -F'\t' '$8 == "Missense_Mutation" || $8 == "Splice_Site" || $8 == "Nonstop_Mutation" || $8 == "Nonsense_Mutation" || $8 == "Translation_Start_Site" {count++} END {print count}' "$vcf_path")
  new_variable=$(echo "$count * 1000000 / $total_length" | bc -l)
  echo "Mutations count of $sample: $count" >> "${out_path}"
  echo "TMB of $sample qualifying regions: $new_variable" >> "${out_path}"

  rm "$name.mosdepth.*"
  rm "$name.per-base.*"
  rm $maf_dir/${cell_type}.maf
}


# 并行处理每个细胞类型
for cell_type in "${cell_types[@]}"; do
  process_cell_type "$cell_type" &
done

# 等待所有后台任务完成
wait


sample="All"
cd ${directory}
/opt/mosdepth -t 6 "All" "$all_bam" 

gunzip -c "All.per-base.bed.gz" > "All.per-base.bed" # 避免不必要解压再压缩
total_length=0
while IFS=$'\t' read -r chromosome start end reads; do
    if [[ $chromosome =~ ^[0-9]+$ && $chromosome -le 22 || $chromosome == "X" || $chromosome == "Y" ]] && [ $reads -gt 0 ]; then
        region_length=$((end - start))
        total_length=$((total_length + region_length))
    fi
done < "$sample.per-base.bed"
rm "$name.mosdepth.*"
rm "$name.per-base.*"

echo "Total length of $sample qualifying regions: $total_length" >> "${out_path}"

count=$(awk -F'\t' '$8 == "Missense_Mutation" || $8 == "Splice_Site" || $8 == "Nonstop_Mutation" || $8 == "Nonsense_Mutation" || $8 == "Translation_Start_Site" {count++} END {print count}' "$maf_path")
new_variable=$(echo "$count * 1000000 / $total_length" | bc -l)
echo "Mutations count of $sample: $count" >> "${out_path}"
echo "TMB of $sample qualifying regions: $new_variable" >> "${out_path}"
