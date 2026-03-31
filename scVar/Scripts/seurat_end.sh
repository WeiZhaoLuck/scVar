#!/bin/bash
input=$1
out1=$2
out2=$3
sed 's/"//g' $input | tail -n +2 | awk -v FS=',' -v OFS="\t" '{print $1,$9}' > $out1
sed 's/"//g' $input | tail -n +2 | awk -v FS=',' -v OFS="\t" '{print $1,$8}' > $out2
