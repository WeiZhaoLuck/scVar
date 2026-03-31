#!/bin/bash

bam_path=$1
file_snp=$2
tbl_file=$3
upn=$4
filter=$5
mapq=$6
baseq=$7
parallel=$8
tbl_all_path=`dirname $file_snp`
# tbl_all_path=$(cd "$dirname "$file_snp"";pwd) 
if [ ! -d "$tbl_file" ]; then
    mkdir -p $tbl_file
else
    rm $tbl_file/*
fi
echo ${tbl_all_path}
echo ${tbl_file}
echo ${upn}
grep -v '#'  ${file_snp} | awk -v FS="\t" -v OFS="\t" '{if($4~/^[A|T|G|C]$/){if($5~/^[A|T|G|C]$/){print $1,$2,$2,$4,$5,"gene","SNV"}}}' | awk '{if($0!="") print}' >  ${tbl_all_path}/test_final.tbl &&\
bash /codes/split.sh ${tbl_all_path}/test_final.tbl $tbl_file &&\
result_name=`dirname ${upn}`
cd ${result_name} && python /codes/genotype.py $bam_path $tbl_file $upn --filter $filter --mapq $mapq --baseq $baseq --parallel $parallel
cd $result_name && cat genotype+*_All.tsv  | grep -v '^MT' > genotype_all.tsv
cd $result_name && cat ${upn}+*_CB.tsv | grep -v '^MT' > ${upn}_all_each_cell_sorted_bybarcode.tsv &&\
rm ${upn}+*_All.tsv
rm ${upn}+*_CB.tsv
rm ${tbl_all_path}/test_final.tbl
