#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-17 14:54:11
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-18 10:20:18
 # @FilePath: \scvar_docker\codes\SNV_Annotation_Celltype2.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

genotype_path=$1
barcode_celltypes_path=$2
annotation_vcf_path=$3

genotype_all=`dirname $genotype_path`
cd $genotype_all &&\

##new
##barcode_normal:
barcode_path_out=`dirname $barcode_celltypes_path`
sed 's/"//g' $barcode_celltypes_path | tail -n +2 | awk -v FS='\t' -v OFS="\t" '{print $1,$2}' > $barcode_path_out/barcodes_cell_type.tsv
sed 's/"//g' $barcode_path_out/output_cell_cluster.tsv | tail -n +2 | awk -v FS='\t' -v OFS="\t" '{print $1,$2}' > $barcode_path_out/barcodes_cluster.tsv &&\


vcf_path=`dirname $annotation_vcf_path`
awk -v FS='\t' -v OFS='\t' 'FNR==NR {key=$1"\t"$2"\t"$4"\t"$5; map[key]; next} {key=$1"\t"$2"\t"$3"\t"$4; if (key in map) print}'  $vcf_path/snp_vaf_final.vcf  genotype_all_alt.tsv > genotype_all_alt2.tsv &&\
awk -v FS='\t' -v OFS='\t' 'FNR==NR {key=$1"\t"$2"\t"$4"\t"$5; map[key]; next} {key=$1"\t"$2"\t"$3"\t"$4; if (key in map) print}'  $vcf_path/snp_vaf_final.vcf  genotype_all_ref.tsv > genotype_all_ref2.tsv &&\
##join cell type with barcode
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1]=$2;next} {print $1"_"$2"_"$3"_"$4,$5,$1,$2,$3,$4,$6,$7,$8,a[$5]}' $barcode_path_out/barcodes_cell_type.tsv genotype_all_alt2.tsv > genotype_alt_cell_type.tsv
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1]=$2;next} {print $1"_"$2"_"$3"_"$4,$5,$1,$2,$3,$4,$6,$7,$8,a[$5]}' $barcode_path_out/barcodes_cell_type.tsv genotype_all_ref2.tsv > genotype_ref_cell_type.tsv

##join cluster with barcode
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1]=$2;next} {print $1"_"$2"_"$3"_"$4,$5,$1,$2,$3,$4,$6,$7,$8,a[$5]}' $barcode_path_out/barcodes_cluster.tsv genotype_all_alt2.tsv > genotype_alt_cluster.tsv
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1]=$2;next} {print $1"_"$2"_"$3"_"$4,$5,$1,$2,$3,$4,$6,$7,$8,a[$5]}' $barcode_path_out/barcodes_cluster.tsv genotype_all_ref2.tsv > genotype_ref_cluster.tsv
# {print $1"_"$2"_"$4"_"$5,$9,$1,$2,$4,$5,$10,$11,$12,a[$9]}

python /codes/add_celltype_barcodes.py genotype_alt_cell_type.tsv genotype_all_alt_sorted_each_cell_celltypes.tsv &&\
python /codes/add_celltype_barcodes.py genotype_ref_cell_type.tsv genotype_all_ref_sorted_each_cell_celltypes.tsv &&\

python /codes/add_celltype_barcodes.py genotype_alt_cluster.tsv genotype_all_alt_sorted_each_cell_clusters.tsv &&\
python /codes/add_celltype_barcodes.py genotype_ref_cluster.tsv genotype_all_ref_sorted_each_cell_clusters.tsv &&\

tail -n +2  genotype_all_alt_sorted_each_cell_celltypes.tsv | awk -v FS="\t" -v OFS="\t" '{print $1,$3,$4}' | sort -k 1  > key_celltypes_barcodes.tsv
## mutation cells count
awk -v OFS="\t" '{ count[$1]++ } END { for (key in count) print key,count[key] }' genotype_alt_cell_type.tsv | sort -k 1 > genotype_mutation_count.tsv

###annotation with mutation cell counts

awk 'BEGIN {FS=OFS="\t"} {print $1"_"$2"_"$3"_"$4,$0}' $annotation_vcf_path | awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0,"mutation_cells_count","cell_types","barcodes"} NR!=1 {print}' > vcf_normaled_final.tsv &&\

tail -n +2 vcf_normaled_final.tsv | sort -k 1 > vcf_normaled_final_sorted.tsv &&\
awk -v OFS="\t" 'FNR==NR{a[$1];next} $1 in a'  genotype_mutation_count.tsv vcf_normaled_final_sorted.tsv > annotation_final.tsv
join -1 1 -2 1 -t $'\t' annotation_final.tsv genotype_mutation_count.tsv > vcf_normaled_final_count.tsv

### vcf2maf
# vcf_path=`dirname $annotation_vcf_path`
bash /codes/vcf2maf.sh $vcf_path $genotype_all/genotype_alt_cell_type.tsv

###annotation with mutation cell type
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1]=$2"\t"$3;next} {print $_,a[$1]}' key_celltypes_barcodes.tsv vcf_normaled_final_count.tsv > all.tsv
head -n 1 vcf_normaled_final.tsv > header.tsv
cat header.tsv all.tsv > all_final.tsv

rm $genotype_all/annotation_final.tsv
rm $genotype_all/genotype_all_alt2.tsv
rm $genotype_all/genotype_all_ref2.tsv

rm $genotype_all/genotype_all_alt_sorted_each_cell_celltypes.tsv
rm $genotype_all/genotype_all_ref_sorted_each_cell_celltypes.tsv

rm $genotype_all/vcf_normaled_final.tsv
rm $genotype_all/vcf_normaled_final_count.tsv
rm $genotype_all/vcf_normaled_final_sorted.tsv
rm $genotype_all/key_celltypes_barcodes.tsv