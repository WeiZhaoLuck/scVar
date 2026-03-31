#!/bin/bash

mutation_path=$1
genotype_alt_celltypes=$2

python /codes/Vcftomaf.py $mutation_path/test_anno.vep.annovar.vcf.hg38_multianno.txt 1 all | awk -v FS='\t' -v OFS='\t' '{print $2"_"$3"_"$5"_"$6,$_}'  > $mutation_path/test_anno.vep.annovar.vcf.hg38_multianno_final.maf &&\
awk '{print $1"_"$2"_"$4"_"$5}' $mutation_path/snp_vaf_final.vcf | awk 'NR==FNR{a[$1]; next} $1 in a' - $mutation_path/test_anno.vep.annovar.vcf.hg38_multianno_final.maf > $mutation_path/filter.maf &&\

awk -v FS='\t' -v OFS='\t' '{print $1,$2}' $genotype_alt_celltypes | awk 'NR==FNR{a[$1]; next} $1 in a' $mutation_path/filter.maf - > $mutation_path/alt_mutations_celltype_filter.tsv &&\
awk -v FS='\t' -v OFS='\t'  'NR==FNR{a[$1]=$0; next} {print $0, a[$1]}' $mutation_path/filter.maf $mutation_path/alt_mutations_celltype_filter.tsv | awk -v FS='\t' -v OFS='\t'  '{$12 = $2}1' > final_test.maf &&\
head -n 1 $mutation_path/test_anno.vep.annovar.vcf.hg38_multianno_final.maf | awk -v FS='\t' -v OFS='\t' '{print "ID","barcode",$_}' | cat - final_test.maf > final.maf

rm $mutation_path/alt_mutations_celltype_filter.tsv
