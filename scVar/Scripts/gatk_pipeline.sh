#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-05 09:32:09
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-27 17:40:57
 # @FilePath: \scvar_docker\codes\gatk_pipeline.sh

### 

bam_path=$1
out_path=$2
tmp_path=`dirname $bam_path`

java -jar "-Xmx20g" /opt/picard.jar SortSam \
    INPUT=${bam_path} \
    OUTPUT=${tmp_path}/SortSam2.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true &&\
# java -jar /opt/picard.jar BuildBamIndex \
#       I=${tmp_path}/SortSam2.bam &&\
# samtools  sort -@ 5 -o ${tmp_path}/SortSam2.bam ${bam_path} &&\
# samtools index ${tmp_path}/SortSam2.bam &&\

java -jar /opt/picard.jar MarkDuplicates \
    INPUT=${tmp_path}/SortSam2.bam \
    OUTPUT=${tmp_path}/possorted_genome_bam_merge_sorted_duplicates.bam \
    METRICS_FILE=${tmp_path}/output_markduplicates_metrics.txt \
    CREATE_INDEX=true &&\
/opt/gatk-4.2.2.0/gatk --java-options "-Xmx60g" SplitNCigarReads \
    -R /reference/fasta/genome.fa \
    -I ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates.bam \
    -O ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates_splited.bam &&\

/opt/gatk-4.2.2.0/gatk --java-options "-Xmx60g" BaseRecalibrator \
-R /reference/fasta/genome.fa \
-I ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates_splited.bam \
--known-sites /reference/1000G_phase1.snps.high_confidence.hg38_nochr_sort.vcf \
--known-sites /reference/Mills_and_1000G_gold_standard.indels.hg38_nochr_sort.vcf \
--known-sites /reference/dbsnp_138.hg38.nonchr.vcf \
-O ${tmp_path}/test.sorted.markdup.recal_data.table &&\

/opt/gatk-4.2.2.0/gatk --java-options "-Xmx40g"  ApplyBQSR \
--bqsr-recal-file ${tmp_path}/test.sorted.markdup.recal_data.table \
-R /reference/fasta/genome.fa \
-I ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates_splited.bam \
-O ${tmp_path}/test.sorted.markdup.BQSR.bam &&\


/opt/gatk-4.2.2.0/gatk --java-options "-Xmx40g"  HaplotypeCaller \
-R /reference/fasta/genome.fa \
-I ${tmp_path}/test.sorted.markdup.BQSR.bam \
-ploidy 40 \
-O ${out_path}/output_gatk_BQSR.vcf.gz &&\

gunzip ${out_path}/output_gatk_BQSR.vcf.gz && mv ${out_path}/output_gatk_BQSR.vcf ${out_path}/variants.vcf

rm ${tmp_path}/SortSam2.bam
rm ${tmp_path}/SortSam2.bai
rm ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates.bam
rm ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates.bai
rm ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates_splited.bam
rm ${tmp_path}/possorted_genome_bam_merge_sorted_duplicates_splited.bai
