#!/bin/bash
###
 # @Author: WeiZhaoLuck zhaow@big.ac.cn
 # @Date: 2025-03-31 10:04:56
 # @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
 # @LastEditTime: 2025-03-31 10:42:18
 # @FilePath: \scvar_docker\codes\Analysis_SigTMBEntroy.sh
### 

work_path=$1/$2
sample=$2
errorlog=$work_path/logs/${sample}_analysis_stderr.log
outlog=$work_path/logs/${sample}_analysis_stdout.log

echo "###Step1: Mutation signature generating........" > ${outlog} 2>> ${errorlog}
bash /codes/mutation_signature.sh  $work_path > ${outlog} 2>> ${errorlog} 

echo "###Step2: TMB generating........" >> ${outlog} 2>> ${errorlog}
bash /codes/TMB.sh $work_path $sample  > ${outlog} 2>> ${errorlog} 

echo "###Step3: Entroy Simpson generating........" >> ${outlog} 2>> ${errorlog}
mkdir $work_path/analysis/entroy_simpson
python /codes/EntroySimpson.py --input $work_path/output_data.h5mu --out  $work_path/analysis/entroy_simpson > ${outlog} 2>> ${errorlog} 