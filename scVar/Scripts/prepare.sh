#!/bin/bash

file_path=$1
gunzip ${file_path}/filtered_feature_bc_matrix/features.tsv.gz && gunzip ${file_path}/soupX_matrix/features.tsv.gz && awk -F '\t' '{print $2}' ${file_path}/filtered_feature_bc_matrix/features.tsv | paste - ${file_path}/soupX_matrix/features.tsv | awk '{$(2+1) = $1; print $0}' | sed 's/[^ ]* //' > ${file_path}/soupX_matrix/test.tsv && mv ${file_path}/soupX_matrix/test.tsv ${file_path}/soupX_matrix/features.tsv && gzip ${file_path}/soupX_matrix/features.tsv && gzip ${file_path}/filtered_feature_bc_matrix/features.tsv