
library(Seurat)
library(monocle)
library(dplyr)
library(cowplot) 

args <- commandArgs(trailingOnly = TRUE)  
matrix_path <- args[1]
mutation_input <- args[2]  # 获取用户输入的第一个参数  
mutation_list <- unlist(strsplit(mutation_input, ","))  # 拆分成向量  

setwd(matrix_path)
# mutation_list <- c('10_100208621_T_C','10_100208676_T_C')
test <- Read10X(data.dir = "matrix_files/")
metadata <- read.csv("metadata.csv")
rownames(metadata) <- metadata$X

test_obj <- CreateSeuratObject(counts = test, meta.data = metadata)
rm(test, metadata)
gc()

for (mutation in mutation_list) {
  # 动态生成列名称
  column_name <- mutation # 假设列名与突变名称相同

  # 检查 mutation 是否以数字开头
  if (grepl("^[0-9]", mutation)) {
    # 如果以数字开头，列名前加 "X"
    col <- paste0("X", column_name)
  } else {
    # 否则直接使用列名
    col <- column_name
  }

  # 使用 ifelse 更新指定列
  test_obj@meta.data[[col]] <- ifelse(test_obj@meta.data[[col]] == 1, "Mutated", "NA")
}

data <- as(as.matrix(test_obj@assays$RNA@counts), "sparseMatrix")
featureData <- data.frame(gene_short_name = row.names(test_obj), row.names = row.names(test_obj))
fd <- new("AnnotatedDataFrame", data = featureData)
phenoData <- test_obj@meta.data
pd <- new("AnnotatedDataFrame", data = phenoData)

cds <- newCellDataSet(data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)
rm(test_obj, data, fd, pd, phenoData, featureData)
gc()
##
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(
  fData(cds),
  num_cells_expressed >= 10
))
### strategy1

# diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~leiden_Ave",cores=4)
# 找出那些在不同细胞类型之间有显著差异表达的基因，以用于细胞排序。
# ordergene <- rownames(subset(diff, qval < 0.01))
# 标记为用于细胞排序的基因


### strategy2
gene_Disp <- dispersionTable(cds)
ordergene <- gene_Disp %>%
  dplyr::filter(mean_expression >= 0.1, dispersion_empirical >= dispersion_fit) %>%
  pull(gene_id) %>%
  unique()

cds <- setOrderingFilter(cds, ordergene)
# 其中黑点就是后面构建发育轨迹需要的基因
# plot_ordering_genes(cds)
cds <- reduceDimension(cds,
  max_components = 2, num_dim = 20,
  method = "DDRTree"
)

cds <- orderCells(cds)
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(tibble)
for (mutation in mutation_list) {
  if (grepl("^[0-9]", mutation)) {
    # 如果以数字开头，列名前加 "X"
    col <- paste0("X", mutation)
  } else {
    # 否则直接使用列名
    col <- mutation
  }

  data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% # 提取坐标
    select_(Component_1 = 1, Component_2 = 2) %>% # 重命名
    rownames_to_column("cells")
  p_data <- pData(cds) %>%
    rownames_to_column("cells") %>%
    select(cells, State, Pseudotime, col)

  # 将 data_df 与 p_data 按照 cells 列横向合并
  merged_df <- left_join(data_df, p_data, by = "cells")
  colnames(merged_df)[colnames(merged_df) == col] <- "mutation"
  merged_df_sorted <- merged_df[order(merged_df$mutation != "Mutated", decreasing = TRUE), ]
  g <- ggplot() +
    geom_point_rast(data = merged_df_sorted, aes(
      x = Component_1,
      y = Component_2,
      color = Pseudotime
    ))
  g2 <- ggplot() +
    geom_point_rast(data = merged_df_sorted, aes(
      x = Component_1,
      y = Component_2,
      color = mutation
    )) # 散点图


  final_plot_g <- plot_grid(g, g2, ncol = 2)
  name <- paste(mutation, ".png", sep = "")
  ggsave(name, plot = final_plot_g, width = 10, height = 3)
}

