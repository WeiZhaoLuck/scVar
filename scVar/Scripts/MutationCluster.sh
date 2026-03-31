#!/bin/bash  
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-31 10:50:33
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-31 14:17:27
 # @FilePath: \scvar_docker\codes\MutationCluster.sh
 
### 
# 初始化变量  
all_path=""  
sample=""  
method=""  
flag=""  
number=""  
clustermethod=""  

# 解析命令行参数  
while [[ $# -gt 0 ]]; do  
    case "$1" in  
        --path)  
            all_path="$2"  
            shift 2  
            ;;  
        --sample)  
            sample="$2"  
            shift 2  
            ;;  
        --method)  
            method="$2"  
            shift 2  
            ;;  
        --flag)  
            flag="$2"  
            shift 2  
            ;;  
        --number)  
            number="$2"  
            shift 2  
            ;;  
        --clustermethod)  
            clustermethod="$2"  
            shift 2  
            ;;  
        *)  
            echo "未知选项: $1"  
            exit 1  
            ;;  
    esac  
done  

# 检查参数是否提供  
if [[ -z "$all_path" || -z "$sample" || -z "$method" || -z "$flag" || -z "$number" || -z "$clustermethod" ]]; then  
    echo "用法: $0 --path <工作路径> --sample <样本> --method <方法(可选: var1|var2|TF_IDF)> --flag <0|1> --number <整数> --clustermethod <聚类方法>"  
    exit 1  
fi  

# 验证method参数  
if [[ "$method" != "var1" && "$method" != "var2" && "$method" != "TF_IDF" ]]; then  
    echo "错误: method 参数必须是 var1、var2 或 TF_IDF。"  
    exit 1  
fi  

# 验证flag参数，确保它是0或1  
if [[ ! "$flag" =~ ^[01]$ ]]; then  
    echo "错误: flag 参数必须是 0 或 1。"  
    exit 1  
fi  

# 将flag转换为整数  
# declare -i flag_int=$flag  

# 验证number参数，确保它是一个整数  
if [[ ! "$number" =~ ^[0-9]+$ ]]; then  
    echo "错误: number 参数必须是一个整数。"  
    exit 1  
fi  
valid_methods=("ward.D" "ward.D2" "single" "complete" "average" "mcquitty" "median" "centroid")  
if [[ ! " ${valid_methods[@]} " =~ " ${clustermethod} " ]]; then  
    echo "错误: clustermethod 参数必须是 ward.D、ward.D2、single、complete、average、mcquitty、median 或 centroid。"  
    exit 1  
fi  

# 赋值给有意义的变量  
work_path=$all_path/$sample  
# sample=$2  
# method=$3  
# flag=$4
# number=$5
# clustermethod=$6
declare -i number_int=$number 
declare -i number_int=$flag 

# 验证method参数  
# if [[ "$method" != "var1" && "$method" != "var2" && "$method" != "TF_IDF" ]]; then  
#     echo "错误: method 参数必须是 var1、var2 或 TF_IDF。"  
#     exit 1  
# fi  

# # 验证第四个参数（flag）  
# if [[ ! "$flag" =~ ^[01]$ ]]; then  
#     echo "错误: flag 参数必须是 0 或 1。"  
#     exit 1  
# fi  


# 验证method参数  
# valid_methods=("ward.D" "ward.D2" "single" "complete" "average" "mcquitty" "median" "centroid")  
# if [[ ! " ${valid_methods[@]} " =~ " ${clustermethod} " ]]; then  
#     echo "错误: clustermethod 参数必须是 ward.D、ward.D2、single、complete、average、mcquitty、median 或 centroid。"  
#     exit 1  
# fi  

errorlog=$work_path/logs/${sample}_analysis_stderr.log
outlog=$work_path/logs/${sample}_analysis_stdout.log
mkdir -p $work_path/analysis/mutation_cluster
awk -v FS='\t' -v OFS=',' '{print $1,$2}' $work_path/seurat/output_cell_cluster.tsv | sed 's/cluster/Cluster/g'> $work_path/analysis/mutation_cluster/df_annotation.csv &&\
nohup python /codes/MutationCluster.py -i $work_path/output_data.h5mu  -o $work_path/analysis/mutation_cluster -m $method -z $flag -n $number >> ${outlog} 2>> ${errorlog} &&\
Rscript /codes/MutationCluster.R $work_path/analysis/mutation_cluster $number $flag $clustermethod >> ${outlog} 2>> ${errorlog}