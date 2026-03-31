#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-07 09:28:47
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-04-08 11:03:47
 # @FilePath: \scvar_docker\codes\Statistics_Report.sh
 # @Description: calculates the number of mutations and genes per cell type and cluster, calculates the cell_type specific mutations, and generates a binary matrix.
### 
work_path=$1/$2
bam_path=$1/$2/mapping/run_count_${2}/outs/extracted.bam
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' $work_path/seurat/output_cell_cluster.tsv | tail -n +2 > $work_path/seurat/barcodes_cluster.tsv
### Statistics
mkdir -p  $work_path/analysis/celltype_mutations
awk -v FS='\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4,$7}' $work_path/mutation/results/variants/vcf_normaled.tsv | awk -F '\t' '{if($2!=""){print $_}}'> $work_path/genotype/mutation.info
python /codes/MutationAndGene_Matrix.py --type1 gene --type2 cluster --ref $work_path/genotype/genotype_ref_cluster.tsv --alt $work_path/genotype/genotype_alt_cluster.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cluster.tsv --out $work_path/analysis/celltype_mutations/ &
P1=$!
python /codes/MutationAndGene_Matrix.py --type1 mutation --type2 cluster --ref $work_path/genotype/genotype_ref_cluster.tsv --alt $work_path/genotype/genotype_alt_cluster.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cluster.tsv --out $work_path/analysis/celltype_mutations/ &
P2=$!
python /codes/MutationAndGene_Matrix.py --type1 gene --type2 cell_type --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cell_type.tsv --out $work_path/analysis/celltype_mutations/ &
P3=$!
python /codes/MutationAndGene_Matrix.py --type1 mutation --type2 cell_type --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv --mutation $work_path/genotype/mutation.info --barcode $work_path/seurat/barcodes_cell_type.tsv --out $work_path/analysis/celltype_mutations/ &
P4=$!

wait $P1 $P2 $P3 $P4
python /codes/Calculate_p.py --type mutation --ref $work_path/analysis/celltype_mutations/all_table_ref_mutation_cluster.csv --alt $work_path/analysis/celltype_mutations/all_table_alt_mutation_cluster.csv --out $work_path/analysis/celltype_mutations/ &
P1=$!
python /codes/Calculate_p.py --type gene --ref $work_path/analysis/celltype_mutations/all_table_ref_gene_cluster.csv --alt $work_path/analysis/celltype_mutations/all_table_alt_gene_cluster.csv --out $work_path/analysis/celltype_mutations/ &
P2=$!
python /codes/Calculate_p.py --type mutation --ref $work_path/analysis/celltype_mutations/all_table_ref_mutation_cell_type.csv --alt $work_path/analysis/celltype_mutations/all_table_alt_mutation_cell_type.csv --out $work_path/analysis/celltype_mutations/ &
P3=$!
python /codes/Calculate_p.py --type gene --ref $work_path/analysis/celltype_mutations/all_table_ref_gene_cell_type.csv --alt $work_path/analysis/celltype_mutations/all_table_alt_gene_cell_type.csv --out $work_path/analysis/celltype_mutations/ &
P4=$!

wait $P1 $P2 $P3 $P4
rm $work_path/analysis/celltype_mutations/all_table_*

# ### Binary Matrix
# python /codes/Binary_Matrix.py --ref $work_path/genotype/genotype_ref_cell_type.tsv --alt $work_path/genotype/genotype_alt_cell_type.tsv -m $work_path/genotype/mutation.info  --out $work_path/CellType_Mutations/

