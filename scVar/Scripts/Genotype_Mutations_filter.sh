#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-07 15:20:20
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-04-08 15:26:38
 # @FilePath: \scvar_docker\codes\Genotype_Mutations_filter.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

genotype_path=$1
vcf_path=$2
vaf=$3

genotype_all=`dirname $genotype_path`
cd $genotype_all &&\
mutation_path=`dirname $vcf_path`

##genotype each cell
#awk -v FS='\t' -v OFS='\t' '{ sums[$1"\t"$2"\t"$4"\t"$5"\t"$9] += $10; sum2[$1"\t"$2"\t"$4"\t"$5"\t"$9] += $11; sum3[$1"\t"$2"\t"$4"\t"$5"\t"$9] += $12 ; } END { for (key in sums) { split(key,arr,"\t"); print key,sums[key],sum2[key],sum3[key]} }' genotype_all.tsv | sort -k 5 > genotype_all_each_cell_sorted_bybarcode.tsv &&\

##define ref cells and alt cells
awk -v FS='\t' '{if($7>=1){print }}'  genotype_all_each_cell_sorted_bybarcode.tsv > genotype_all_alt.tsv
awk -v FS='\t' '{if($7==0){print }}'  genotype_all_each_cell_sorted_bybarcode.tsv > genotype_all_ref.tsv
awk -v FS='\t' -v OFS='\t' '{key=$1"\t"$2"\t"$3"\t"$4; count[key]++} END {for (key in count) {if (count[key] >= 2) print key}}'  genotype_all_alt.tsv > Filtered_mutation.tsv
#awk -v FS='\t' -v OFS='\t' '{key=$1"\t"$2"\t"$3"\t"$4; count[key]++} END {for (key in count) {if (count[key] >= 2) print key}}' genotype_all_alt.tsv > Filtered_mutation.tsv

##filter vcf with filtered mutations
awk -v FS='\t' -v OFS='\t' 'FNR==NR {key2=$1"\t"$2"\t"$3"\t"$4; map[key2]=1; next} {key1=$1"\t"$2"\t"$4"\t"$5; if (map[key1]) print}' Filtered_mutation.tsv $vcf_path > $mutation_path/vcf_filtered.tsv

##filter genotype result with filtered mutations
#awk -v FS='\t' -v OFS='\t' 'FNR==NR {key2=$1"\t"$2"\t"$3"\t"$4; map[key2]=1; next} {key1=$1"\t"$2"\t"$4"\t"$5; if (map[key1]) print}' Filtered_mutation.tsv genotype_all.tsv > genotype_all_filtered.tsv

##final vcf file
grep '^#CHROM' $vcf_path | awk -v FS='\t' -v OFS='\t' '{print $0"\tVAF\tbarcode"}' > $mutation_path/header.vcf
#cat $mutation_path/header.vcf $mutation_path/vcf_filtered.tsv > $mutation_path/vcf_filtered_final.vcf

### filter vcf with vaf
perl_vaf=/codes/filter_vaf.pl
perl ${perl_vaf} $mutation_path/vcf_filtered.tsv ${vaf} > $mutation_path/snp_vaf_before.vcf
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1$2$3$4];next} !($1$2$4$5 in a){print $_}' /reference/filter_all_001.vcf $mutation_path/snp_vaf_before.vcf > $mutation_path/snp_vaf.vcf
bash /codes/alt_cells.sh genotype_all_alt.tsv > genotype_all_alt_cells.tsv
#awk -v FS='\t' -v OFS='\t' '{key2=$1"\t"$2"\t"$3"\t"$4; barcode=$5; if (key2 in data) {if (match(","data[key2]",", ","barcode",")) {next} else {data[key2]=data[key2]","barcode}} else {data[key2]=barcode}} END {for (key in data) {print key"\t"substr(data[key], 2)}}' genotype_all_alt.tsv > genotype_all_alt_cells.tsv
awk -v FS='\t' -v OFS='\t' 'FNR==NR {key2=$1"\t"$2"\t"$3"\t"$4; map[key2]=$5; next} {key1=$1"\t"$2"\t"$4"\t"$5; if (map[key1]) print $_,map[key1]}' genotype_all_alt_cells.tsv  $mutation_path/snp_vaf.vcf > $mutation_path/snp_vaf_barcodes.vcf

cat $mutation_path/header.vcf $mutation_path/snp_vaf_barcodes.vcf > $mutation_path/snp_vaf_final.vcf
grep '^##'  $vcf_path | cat  - $mutation_path/snp_vaf_final.vcf > $mutation_path/Mutations.vcf

grep -v '^#' $mutation_path/snp_vaf.vcf | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$4,$5}' > $mutation_path/Filtered_mutation_final.tsv
awk -v FS='\t' -v OFS='\t' 'FNR==NR {key2=$1"\t"$2"\t"$3"\t"$4; map[key2]=1; next} {key1=$1"\t"$2"\t"$4"\t"$5; if (map[key1]) print}' $mutation_path/Filtered_mutation_final.tsv ${genotype_path} > genotype_all_filtered.tsv

rm $mutation_path/snp_vaf.vcf
rm $mutation_path/snp_vaf_barcodes.vcf
rm $mutation_path/snp_vaf_before.vcf

rm $mutation_path/vcf_filtered.tsv
rm $genotype_all/genotype_all_each_cell_sorted_bybarcode.tsv


