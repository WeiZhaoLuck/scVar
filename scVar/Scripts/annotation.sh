#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-07 15:18:36
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-08 10:51:43
 # @FilePath: \scvar_docker\codes\annotation.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 
vcf_path=$1
out_path=`dirname $vcf_path`
vcfanno_coml=$2
genotype_path=$3
barcode_path=$4
vcfanno_thre=$5


#vep annotation
/opt/ensembl-vep/vep -i ${vcf_path} --fork 6 --sift b -filter "SIFT is deleterious" --keep_csq  --vcf --everything --pick --vcf_info_field CSQ -o ${out_path}/test_anno.vep.vcf --assembly GRCh38 --cache --cache_version 107 --dir /reference/vep_data --offline --fasta /reference/fasta/genome.fa --force_overwrite &&\


##annovar annotation

/codes/table_annovar.pl ${out_path}/test_anno.vep.vcf  /reference/humandb -buildver hg38 -out ${out_path}/test_anno.vep.annovar.vcf -remove -protocol refGene,ensGene,knownGene,cytoBand,avsnp150,dbnsfp41a,clinvar_20210123,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,exac03  -operation g,g,g,r,f,f,f,f,f,f,f,f,f,f  -nastring . -vcfinput -polish &&\
gzip ${out_path}/test_anno.vep.annovar.vcf.hg38_multianno.vcf &&\

##vcfanno annotation
vcfanno  -p ${vcfanno_thre} ${vcfanno_coml} ${out_path}/test_anno.vep.annovar.vcf.hg38_multianno.vcf.gz > ${out_path}/scvar_anno.vep.annovar.vcfanno.vcf

grep -v '^##' ${out_path}/scvar_anno.vep.annovar.vcfanno.vcf > ${out_path}/tmp.vcf &&\
python /codes/normalize_new_v2.py  ${out_path}/tmp.vcf ${out_path}

bash /codes/SNV_Annotation_Celltype2.sh $genotype_path $barcode_path ${out_path}/vcf_normaled.tsv


rm ${out_path}/scvar_anno.vep.annovar.vcfanno.vcf
