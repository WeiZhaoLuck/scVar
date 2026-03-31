#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-21 16:27:51
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-31 10:00:58
 # @FilePath: \scvar_docker\codes\mutation_signature.sh
### 
work_path=$1
mkdir -p $work_path/analysis/mutation_signature
signature_path=$work_path/analysis/mutation_signature
mkdir -p $signature_path/cell_type
mutation_path=$work_path/genotype
barcode_cell_type_path=$work_path/seurat
cd $mutation_path && perl -lane '$raw=$_;my @result=split("\t",$raw);if($result[7]>0){print $_}' genotype_all_alt.tsv | sort -k 5 > genotype_all_alt_sorted.tsv &&\
sed 's/"//g' $barcode_cell_type_path/output_cell_barcodes.tsv | tail -n +2 | awk -v FS='\t' -v OFS="\t" '{print $1,$2}' | sed 's/ /_/g' > $barcode_cell_type_path/barcodes_cell_type.tsv &&\

#sort -k 8  $mutation_path/genotype_all_alt.tsv >  $mutation_path/genotype_all_alt_sorted.tsv &&\

join -t $'\t' -1 5 -2 1 $mutation_path/genotype_all_alt_sorted.tsv $barcode_cell_type_path/barcodes_cell_type.tsv > $signature_path/genotype_cell_type.tsv &&\
awk -v FS=' ' -v OFS='\t' '{print $2,$3,$4,$5,"Sample_all"}' $signature_path/genotype_cell_type.tsv > $signature_path/Mut.table.for.signatures_all.tsv &&\
awk -v FS='\t' -v OFS='\t' '{print $2,$3,$4,$5,$9}' $signature_path/genotype_cell_type.tsv > $signature_path/Mut.table.for.signatures_cell_type.tsv &&\
cd $signature_path &&\
awk -v FS='\t' -v OFS='\t' '$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)$/{print $5,"chr"$1,$2,$3,$4}' Mut.table.for.signatures_all.tsv > Mut.table.for.signatures_all_normal.tsv
awk -v FS='\t' -v OFS='\t' '$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)$/{print $5,"chr"$1,$2,$3,$4}' Mut.table.for.signatures_cell_type.tsv > Mut.table.for.signatures_cell_type_normal.tsv
###mutation signature
echo "###Step8: Mutation signaturecd  generating........" &&\
MutSigRun=/codes/MutSigRunV2.R
#importantvcf_builder=${cellranger_results}/final_annovar.hg38
#cd ${cellranger_results} && python ${ConvertVcf2Maf} --invcf ${invcf} --outtable ${outtable} && echo "[succeed]:<vcf2maf>" >> ${logfile} || echo "[failed]:<vcf2maf>" >> ${errorlog} &&\
/usr/bin/Rscript ${MutSigRun} && echo "[succeed]:<Generate mutation signature>"
echo "Complete!" > $signature_path/log.txt
