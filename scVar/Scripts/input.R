# /*
#  * @Author: WeiZhaoLuck zhaow@big.ac.cn
#  * @Date: 2025-03-04 16:32:51
#  * @LastEditors: WeiZhaoLuck zhaow@big.ac.cn
#  * @LastEditTime: 2025-03-04 16:39:02
#  * @FilePath: \scvar_docker\codes\input.R
#  */

library(Seurat)
# library(dplyr)
# library(patchwork)
# library(SoupX)
# library(celldex) ##http://bioconductor.org/packages/release/data/experiment/html/celldex.html
# library(scRNAseq) ##https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html
# library(SingleR)
# library(scater) ##https://bioconductor.org/packages/release/bioc/html/scater.html
library(sceasy)

args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
output_path <- args[2]


# /**
#   * @description: extract_celltype_from_rds
#   * @param {string} rds_file_path
#   * @param {string} raw_data_path
#   * @return {*}
#  */
extract_celltype_from_rds <- function(rds_file_path, raw_data_path) {  
    # 读取RDS文件  
    object_test <- readRDS(rds_file_path)  
    out_h5ad <- file.path(output_path, "rna.h5ad")
    sceasy::convertFormat(object_test, from = "seurat", to = "anndata", outFile = out_h5ad)
    # 获取meta.data  
    meta_data <- object_test@meta.data  
    
    # 定义要匹配的关键字  
    keywords <- c("celltype", "celltypes", "label")  
    
    # 使用grep进行模糊匹配，查找包含关键字的列名  
    matched_columns <- unlist(lapply(keywords, function(keyword) {  
      grep(keyword, colnames(meta_data), value = TRUE)  
    }))  
    
    # 选择第一个匹配的列名  
    selected_column <- if (length(matched_columns) > 0) matched_columns[1] else NULL  
    
    # 如果找到有效的列名，提取细胞barcode和对应的列  
    if (!is.null(selected_column)) {  
      output_data <- data.frame(  
        cell = rownames(meta_data),  
        celltype = gsub(" ", "_", meta_data[[selected_column]])  
      )  
      # seurat_path <- file.path(raw_data_path, "seurat")

      # 判断文件或目录是否存在
      if (file.exists(raw_data_path)) {
        cat("File or directory exists:", raw_data_path, "\n")
      } else {
        dir.create(raw_data_path, recursive = TRUE)
        cat("File or directory does not exist:", raw_data_path, "\n")
      }
      # 输出结果的文件路径  
      celltype_out <- file.path(raw_data_path, "output_cell_barcodes.tsv")
      
      # 输出结果为TSV文件  
      write.table(output_data, celltype_out, sep = "\t", row.names = FALSE, quote = FALSE)  
      cat("Output successful, file path is:", celltype_out, "\n")
    } else {  
      cat("No valid column names found.\n")
    }  


    keywords <- c("leiden", "cluster")  
    
    # 使用grep进行模糊匹配，查找包含关键字的列名  
    matched_columns <- unlist(lapply(keywords, function(keyword) {  
      grep(keyword, colnames(meta_data), value = TRUE)  
    }))  
    
    # 选择第一个匹配的列名  
    selected_column <- if (length(matched_columns) > 0) matched_columns[1] else NULL  
    
    # 如果找到有效的列名，提取细胞barcode和对应的列  
    if (!is.null(selected_column)) {  
      output_data <- data.frame(  
        cell = rownames(meta_data),  
        cluster = gsub(" ", "_", meta_data[[selected_column]])
      )  
      # seurat_path <- file.path(raw_data_path, "seurat")

      # 判断文件或目录是否存在
      if (file.exists(raw_data_path)) {
        cat("File or directory exists:", raw_data_path, "\n")
      } else {
        dir.create(raw_data_path, recursive = TRUE)
        cat("File or directory does not exist:", raw_data_path, "\n")
      }
      # 输出结果的文件路径  
      celltype_out <- file.path(raw_data_path, "output_cell_cluster.tsv")
      
      # 输出结果为TSV文件  
      write.table(output_data, celltype_out, sep = "\t", row.names = FALSE, quote = FALSE)  
      cat("Output successful, file path is:", celltype_out, "\n")
    } else {  
      cat("No valid column names found.\n")
    }  
}  


extract_celltype_from_rds(rds_path, output_path)

