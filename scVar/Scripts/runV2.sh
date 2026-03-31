#!/bin/bash

work_path=$1
bam_path=$2
SNV_annotation=$work_path/genotype/all_final.tsv
out_path=$work_path/report/data

###bam  statistics
mkdir $work_path/bam_statistics
/opt/bamdst-master/bamdst -p /codes/GRCh38.KnownGene.cds.bed -o $work_path/bam_statistics $bam_path &&\
tail -n +2 $work_path/bam_statistics/chromosomes.report | sort -k1,1n -k1,1V |  sed '1d' | awk -v OFS="," '{print $1,$3,$5}' |  awk -F',' '$1 !~ /_/'  > $work_path/bam_statistics/chromosomes_results_sorted.csv

##top mutations
cols=`awk -F'\t' '{print NF;exit}' $SNV_annotation`
mutation_count=$(($cols -2))
genotype_path=`dirname $SNV_annotation`
tail -n +2 $SNV_annotation | sort -t $'\t' -k $mutation_count,$mutation_count -n -r | head -n 30 > $genotype_path/top30.tsv &&\
cat $genotype_path/header.tsv $genotype_path/top30.tsv > $genotype_path/top30_final.tsv &&\
sed 's/ /_/g' $work_path/seurat/barcodes_cell_type.tsv | cut -f 2 | sort | uniq -c | awk '{print $2","$1}' > $work_path/seurat/result.csv
### data to json
mkdir -p $work_path/report/data
python /codes/csvtojson_new.py $work_path

### top cell
tr -d "['\[\]]" < $genotype_path/top30_cell.tsv > $out_path/top30_cell_final.tsv


### top cell highlight
mkdir $out_path/top_mutation_tsne
mkdir $out_path/top_mutation_umap
mapping_path=`dirname $bam_path`
Rscript /codes/top_mutation.R $out_path/top30_cell_final.tsv  $mapping_path/soupX_matrix $out_path

### Statistics
mkdir $work_path/statistics
awk -v FS='\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4,$9}' $work_path/mutation/results/variants/vcf_normaled.tsv | awk -F '\t' '{if($2!=""){print $_}}'> $work_path/genotype/mutation.info
python /codes/MutationAndGene_Matrix.py --type1 gene --type2 cluster --ref $work_path/genotype/genotype_ref_cluster.tsv --alt $work_path/genotype/genotype_alt_cluster.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cluster.tsv --out $work_path/statistics/ &
P1=$!
python /codes/MutationAndGene_Matrix.py --type1 mutation --type2 cluster --ref $work_path/genotype/genotype_ref_cluster.tsv --alt $work_path/genotype/genotype_alt_cluster.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cluster.tsv --out $work_path/statistics/ &
P2=$!
python /codes/MutationAndGene_Matrix.py --type1 gene --type2 cell_type --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cell_type.tsv --out $work_path/statistics/ &
P3=$!
python /codes/MutationAndGene_Matrix.py --type1 mutation --type2 cell_type --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cell_type.tsv --out $work_path/statistics/ &
P4=$!

wait $P1 $P2 $P3 $P4
python /codes/Calculate_p.py --type mutation --ref $work_path/statistics/all_table_ref_mutation_cluster.csv --alt $work_path/statistics/all_table_alt_mutation_cluster.csv --out $work_path/statistics/ &
python /codes/Calculate_p.py --type gene --ref $work_path/statistics/all_table_ref_gene_cluster.csv --alt $work_path/statistics/all_table_alt_gene_cluster.csv --out $work_path/statistics/ &
python /codes/Calculate_p.py --type mutation --ref $work_path/statistics/all_table_ref_mutation_cell_type.csv --alt $work_path/statistics/all_table_alt_mutation_cell_type.csv --out $work_path/statistics/ &
python /codes/Calculate_p.py --type gene --ref $work_path/statistics/all_table_ref_gene_cell_type.csv --alt $work_path/statistics/all_table_alt_gene_cell_type.csv --out $work_path/statistics/ &
python /codes/Binary_Matrix.py --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv -m $work_path/genotype/mutation.info  --out $work_path/statistics/



##report
js_rescorce=/codes/sources
scripts=/codes/scripts
html_path=/codes/test.html
cp -r $js_rescorce $work_path/report
cp -r $scripts $work_path/report
cp $work_path/seurat/cell_type.png $work_path/report/data
cp $html_path $work_path/report
