#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-21 16:32:30
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-04-08 14:51:09
 # @FilePath: \scvar_docker\codes\Statistic.sh

### 
work_path=$1/$2/analysis
bam_path=$1/$2/mapping/run_count_${2}/outs/extracted.bam
SNV_annotation=$1/$2/genotype/all_final.tsv
out_path=$work_path/report/data

#bam  statistics
mkdir $work_path/bam_statistics
/opt/bamdst-master/bamdst -p /codes/GRCh38.KnownGene.cds.bed -o $work_path/bam_statistics $bam_path &&\
tail -n +2 $work_path/bam_statistics/chromosomes.report | sort -k1,1n -k1,1V |  sed '1d' | awk -v OFS="," '{print $1,$3,$5}' |  awk -F',' '$1 !~ /_/'  > $work_path/bam_statistics/chromosomes_results_sorted.csv

##top mutations
cols=`awk -F'\t' '{print NF;exit}' $SNV_annotation`
mutation_count=$(($cols -2))
genotype_path=`dirname $SNV_annotation`
tail -n +2 $SNV_annotation | sort -t $'\t' -k $mutation_count,$mutation_count -n -r | head -n 30 > $genotype_path/top30.tsv &&\
cat $genotype_path/header.tsv $genotype_path/top30.tsv > $genotype_path/top30_final.tsv &&\
sed 's/ /_/g' $1/$2/seurat/output_cell_barcodes.tsv| tail -n +2 | cut -f 2 | sort | uniq -c | awk '{print $2","$1}' > $1/$2/seurat/result.csv


### data to json
mkdir -p $work_path/report/data
python /codes/csvtojson_new.py $1/$2

## top cell
awk -F'\t' '{print $1 "\t" $NF}' $genotype_path/top30.tsv  | tr -d "['\[\]]" | awk '{print NR-1 "\t" $0}' >  $genotype_path/top30_cell_final.tsv
# tr -d "['\[\]]" < $genotype_path/top30_cell.tsv > $out_path/top30_cell_final.tsv


### top cell highlight
mkdir $out_path/top_mutation_tsne
mapping_path=`dirname $bam_path`
python /codes/Topumap.py $1/$2/output_data.h5mu $genotype_path/top30_cell_final.tsv $work_path/report/data

##report
js_rescorce=/codes/sources
scripts=/codes/scripts
html_path=/codes/test.html
cp -r $js_rescorce $work_path/report
cp -r $scripts $work_path/report
# cp $1/seurat/cell_type.png $work_path/report/data
cp $html_path $work_path/report

rm $genotype_path/top30*
rm $1/$2/seurat/result.csv