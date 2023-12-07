#!/bin/bash
sample1=$1
sample2=$2
path_sample1=$3
path_sample2=$4
result_path=$5
file_path=$6
pre1=$7
pre2=$8

### Reference
genomefa_path=genome.fa

### Scripts
regenotyping=regenotyping.py

### Bam file to regions
bedtools bamtobed -i ${path_sample2}/possorted_genome_bam_sorted.bam > ${path_sample2}/possorted_genome_bam_sorted.bed &&\
bedtools bamtobed -i ${path_sample1}/possorted_genome_bam_sorted.bam > ${path_sample1}/possorted_genome_bam_sorted.bed &&\
awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3}' ${path_sample2}/possorted_genome_bam_sorted.bed | uniq > ${path_sample2}/possorted_genome_bam_sorted_u.bed &
pid_task1=$!
awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3}' ${path_sample1}/possorted_genome_bam_sorted.bed | uniq > ${path_sample1}/possorted_genome_bam_sorted_u.bed &
pid_task2=$!
wait $pid_task1
wait $pid_task2
bedtools merge -i ${path_sample1}/possorted_genome_bam_sorted_u.bed > ${path_sample1}/possorted_genome_bam_sorted_final.bed &
pid_task1=$!
bedtools merge -i ${path_sample2}/possorted_genome_bam_sorted_u.bed > ${path_sample2}/possorted_genome_bam_sorted_final.bed &
pid_task2=$!
wait $pid_task1
wait $pid_task2

### Get common regions between sample1 and sample2
bedtools intersect -a ${path_sample1}/possorted_genome_bam_sorted_final.bed -b ${path_sample2}/possorted_genome_bam_sorted_final.bed > ${result_path}/comm_s1_s2_final.bed &&\

### Extract bam files with common regions
samtools view -bh -L ${result_path}/comm_s1_s2_final.bed ${path_sample1}/possorted_genome_bam_sorted.bam > ${result_path}/possorted_genome_bam_${sample1}.bam &
pid_task1=$!
samtools view -bh -L ${result_path}/comm_s1_s2_final.bed ${path_sample2}/possorted_genome_bam_sorted.bam > ${result_path}/possorted_genome_bam_${sample2}.bam &
pid_task2=$!
wait $pid_task1
wait $pid_task2
samtools index -@ 6 ${result_path}/possorted_genome_bam_${sample1}.bam &
pid_task1=$!
samtools index -@ 6 ${result_path}/possorted_genome_bam_${sample2}.bam &
pid_task2=$!
wait $pid_task1
wait $pid_task2

### Extract 10 percentage and 90 percentage separately of above extracted bam files as new bam file
java -Xmx16g -jar /software/biosoft/software/picard2.22.8/picard/build/libs/picard.jar DownsampleSam I=${result_path}/possorted_genome_bam_${sample1}.bam O=${result_path}/possorted_genome_bam_${sample1}_${pre1}.bam p=${pre1} R=100 &
pid_task1=$!
java -Xmx16g -jar /software/biosoft/software/picard2.22.8/picard/build/libs/picard.jar DownsampleSam I=${result_path}/possorted_genome_bam_${sample2}.bam O=${result_path}/possorted_genome_bam_${sample2}_${pre2}.bam p=${pre2} R=100 &
pid_task2=$!
wait $pid_task1
wait $pid_task2
echo ${result_path}/possorted_genome_bam_${sample1}_${pre1}.bam  > ${result_path}/all_path
echo ${result_path}/possorted_genome_bam_${sample2}_${pre2}.bam  >> ${result_path}/all_path

### Merge 10% bam file and 90% bam file to merge.bam
samtools merge -p -c -f -b ${result_path}/all_path  -@ 8 ${result_path}/possorted_genome_bam_merge.bam &&\
samtools index -@ 4 ${result_path}/possorted_genome_bam_merge.bam &&\

### Get germline mutations from 10% bam file and 90% bam file
mkdir ${result_path}/sample1
mkdir ${result_path}/sample2
cd ${result_path}/sample1 && python strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam ${result_path}/possorted_genome_bam_${sample1}_${pre2}.bam \
--referenceFasta ${genomefa_path} \
--rna --useAllDataForSequenceErrorEstimation --exome  --runDir ${result_path}/sample1 && python ${result_path}/8 sample1/runWorkflow.py -m local -j 30 &
pid_task1=$!
cd ${result_path}/sample2 && python strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam ${result_path}/possorted_genome_bam_${sample2}_${pre2}.bam \
--referenceFasta ${genomefa_path} \
--rna --useAllDataForSequenceErrorEstimation --exome --runDir ${result_path}/sample2 && python ${result_path}/sample2/runWorkflow.py -m local -j 30 &
pid_task2=$!
wait $pid_task1
wait $pid_task2


### Get variants with vaf == 100 and depth greater than 5 
grep -v '^#' ${result_path}/sample1/results/variants/variants.vcf | awk -F "\t" '{if($4~/^[A|T|G|C]$/){if($5~/^[A|T|G|C]$/){print $_}}}' |  perl -lane '$raw=$_;if(/#CHROM/){print $_."\t"."VAF"}else{ if(/:(\d+\,\d*):/){my @result=split /\,/,$1;if(($result[0]+$result[1])>5){my $vaf=$result[1]*100/($result[0]+$result[1]);if($vaf==100){print $raw."\t".$vaf}}}}' > ${result_path}/sample1/results/variants/germline.vcf &
pid_task1=$!
grep -v '^#' ${result_path}/sample2/results/variants/variants.vcf | awk -F "\t" '{if($4~/^[A|T|G|C]$/){if($5~/^[A|T|G|C]$/){print $_}}}' |  perl -lane '$raw=$_;if(/#CHROM/){print $_."\t"."VAF"}else{ if(/:(\d+\,\d*):/){my @result=split /\,/,$1;if(($result[0]+$result[1])>5){my $vaf=$result[1]*100/($result[0]+$result[1]);if($vaf==100){print $raw."\t".$vaf}}}}' > ${result_path}/sample2/results/variants/germline.vcf &
pid_task2=$!
grep -v '^#' ${result_path}/sample2/results/variants/variants.vcf | awk -F "\t" '{if($4~/^[A|T|G|C]$/){if($5~/^[A|T|G|C]$/){print $_}}}'  > ${result_path}/sample2/results/variants/snp.vcf &
pid_task3=$!
wait $pid_task1
wait $pid_task2
wait $pid_task3

awk -v FS='\t' -v OFS='\t' '{print $1,$2}' ${result_path}/sample1/results/variants/germline.vcf >${result_path}/sample1/results/variants/sample1_germ_position &&\
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' ${result_path}/sample2/results/variants/germline.vcf >${result_path}/sample2/results/variants/sample2_germ_position &&\
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' ${result_path}/sample2/results/variants/snp.vcf > ${result_path}/sample2/results/variants/sample2_snp_position &&\

### Make whitelist
samtools depth ${result_path}/possorted_genome_bam_${sample2}_${pre2}.bam > ${file_path}/${pre2}.depth &&\
grep -vFf ${result_path}/sample2/results/variants/sample2_snp_position ${result_path}/sample1/results/variants/sample1_germ_position > ${file_path}/whitelist.pre &&\
awk -v FS='\t' -v OFS='\t' '{if($3>0){print $1,$2}}' ${file_path}/${pre2}.depth > ${file_path}/${pre2}.dep &&\
cat ${file_path}/90xbam.dep ${file_path}/whitelist.pre | sort | uniq -d > ${file_path}/whitelist &&\
awk 'FNR==NR{a[$1$2]; next} ($1$2 in a)' whitelist ${result_path}/sample1/results/variants/germline.vcf  > whitelist.vcf &&\
touch whitelist_${pre1}.vcf
python ${regenotyping} $bam_f1 whitelist.vcf whitelist_${pre1}.vcf
touch whitelist_${pre2}.vcf
python ${regenotyping} $bam_f2 whitelist.vcf whitelist_${pre2}.vcf
awk -F '\t' '{if($10!=0){print }}' whitelist_${pre1}.vcf | awk -v FS='\t' -v OFS='\t' '{if($8==0){print $1,$2}}' | grep '^[0-9XY]'  > whitelist_${pre1}_filter
awk -F '\t' '{if($10!=0){print }}' whitelist_${pre2}.vcf | awk -v FS='\t' -v OFS='\t' '{if($9==0){print $1,$2}}' | grep '^[0-9XY]'  > whitelist_${pre2}_filter
cat whitelist_${pre1}_filter whitelist_${pre2}_filter | sort | uniq -d > whitelist_filter_final

### Filter Whitelist
awk 'FNR==NR{a[$1$2]; next} ($1$2 in a)' whitelist_filter_final ${result_path}/sample1/results/variants/germline.vcf  > whitelist_filter_final.vcf &&\
grep '^[0-9XY]' whitelist_filter_final.vcf | perl -lane '$raw=$_;if(/#CHROM/){print $_."\t"."VAF"}else{ if(/:(\d+\,\d*):/){my @result=split /\,/,$1;if(($result[0]+$result[1])>5){print $_}}}' | perl -lane 'if(/:(\d.*):/){my @result=split /:/,$1;my @adf=split /,/,$result[5];my @adr=split /,/,$result[6];if($adf[1]!=0 && $adr[1]!=0){print $_}}' |perl -lane '$raw=$_;if(/MQ=(\d.*)/){@MQ=split /\t/,$1;if($MQ[0]==255){print $raw}}' | awk '{print $1"_"$2}'> Whitelist