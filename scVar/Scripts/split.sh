#!/bin/bash
function split_tbl(){
  threads=7
  tbl="$1"
  chr=`basename $tbl | sed 's/\.[^.]*$//'`
  echo ${chr}
  total_wc=$(wc -l $tbl | awk '{print $1}')
  line_per_file=$(echo "scale=0; ($total_wc+$((threads-1)))/$threads" | bc)
  cd $split_out && split --numeric-suffixes=1 --suffix-length=1 --additional-suffix=.tbl -l "$line_per_file" "$tbl" "$chr."
}

tbl_path_all=$1
split_out=$2
cd ${split_out} && awk -F '\t' '{print > ($1".tbl")}' ${tbl_path_all} &&\
array_name=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
for var in `find $split_out -name "*.tbl"`
do
  chr_name=`basename $var`
  chr=(${chr_name//.tbl/ })
  if [[ "${array_name[@]}"  =~ ${chr[0]} ]]; then
      echo "${chr[0]} exists"
      split_tbl $var
      rm $var
#  else
#      cp $tbl_path_all/${chr}.tbl $split_out/.
  fi
#  split_tbl $var
done
