#!/bin/sh
bam_path=$1
barcodes_path=$2
parallel=$3
outdir_path=`dirname $bam_path`
tmp_path=${outdir_path}/bam_tmp
bam_base=`basename ${bam_path}`
index_bam=`echo "${bam_base%%.*}"`

#samtools sort -o ${outdir_path}/${index_bam}_sorted.bam -@ ${parallel} ${bam_path} &&\
#samtools index ${outdir_path}/${index_bam}_sorted.bam
#mkdir ${outdir_path}/bam_tmp
#bamtools split -in ${outdir_path}/${index_bam}_sorted.bam -reference &&\
#mv ${outdir_path}/*.REF_* ${tmp_path}/ &&\
#python /codes/multi_extract.py ${barcodes_path} ${tmp_path} &&\
find ${tmp_path}/ -name "*.extract.bam" > ${tmp_path}/all_bam_path  &&\
samtools merge -f -b ${tmp_path}/all_bam_path -@ ${parallel} ${outdir_path}/extracted.bam &&\
samtools index ${outdir_path}/extracted.bam



