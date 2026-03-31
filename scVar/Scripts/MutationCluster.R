library(ComplexHeatmap)  
library(circlize)  # 用于颜色设置 
library(dplyr)
# library(Cairo)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)  

condition <- as.numeric(args[3])
# 条件判断  
if (condition == 0) {  
  # 如果 condition 为 1，执行以下操作  
  color_gradient <- c("gray","white", "red")
  AT <- c(-1,0,1)
  label <- c("Absent","No","Yes")
  # 这里可以加入具体的操作代码  
} else if (condition == 1) {  
  # 如果 condition 为 1，执行以下操作  
  color_gradient <- colorRampPalette(c("white", "red"))(100) 
  AT <- c(0,1)
  label <- c("No","Yes")
} 
allowed_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
selected_method <- args[4]

# 验证用户输入的方法是否在允许的列表中  
if (!(selected_method %in% allowed_methods)) {  
  stop("无效的方法: ", selected_method, "\n可用的选项: ", paste(allowed_methods, collapse = ", "))  
}  
setwd(args[1])
##read file
feature_matrix <- read.csv("./mutation_matrix_features.csv",row.names = 1,  check.names = FALSE)  
cell_info <- read.csv("./df_annotation.csv")  
filtered_data_frame <- feature_matrix[, 1:as.numeric(args[2])]  

##annotation dataframe
annotation_col_all <- data.frame(  
  Cluster = factor(cell_info$Cluster)  # 使用 Cluster 列作为 CellType  
)  
rownames(annotation_col_all) <- cell_info$cell  # 设置行名为细胞名  


annotation_col <- subset(annotation_col_all, rownames(annotation_col_all) %in% rownames(filtered_data_frame)) 
annotation_col$Cluster <- as.character(annotation_col$Cluster) 

 
jaccard_distance <- dist(filtered_data_frame, method = "binary") 


# 绘制热图，使用 Jaccard 距离进行行聚类  

transposed_matrix <- as.data.frame(t(filtered_data_frame))
# 设置SVG输出  



cellwidth <- 0.2
cellheight <- 1
pdf("./heatmap_final.pdf", width = cellwidth * nrow(cell_info) / 25 + 1, height = as.numeric(args[2]) / 50 * 2)
print(pheatmap(transposed_matrix,  
               color = color_gradient,
               clustering_distance_cols = jaccard_distance,  # 使用 Jaccard 距离  
               clustering_method = selected_method,  # 聚类方法  
               use_raster = FALSE,
               cluster_rows = TRUE,  # 行聚类  
               cluster_cols = TRUE,  # 列聚类  
               show_rownames = FALSE, 
               show_colnames = FALSE,
               row_title = "Mutations",  
               column_title = "Cells",  
               annotation_col = annotation_col,  # 添加列注释  色  
               fontsize_row = 10,  # 设置行字体大小  
               fontsize_col = 10,
               name="Mutated status",
               cellwidth = cellwidth,
               cellheight = cellheight,
               heatmap_legend_param=list(
                 at=AT,
                 labels=label,
                 title="Mutation",
                 legend_height=unit(4,"cm"),
                 title_position="topcenter",
                 color_bar = "discrete",
                 border = gpar(col = "black", lwd = 2)
               )) )  # 设置列字体大小  )+ 
# Modify the legend  
dev.off()  

# Add a custom legend below the existing legend  
# Create a new viewport for the custom legend  
#pushViewport(viewport(x = 0.82, y = 0.13, width = 0.1, height = 0.2, just = c("left", "bottom")))  

# Draw the custom legend squares with borders  
# Red square for Mutated  
#grid.rect(gp = gpar(fill = "red", col = "black"), width = unit(0.7, "lines"), height = unit(0.7, "lines"))  
#grid.text("Mutated", x = 0.7, y = 0.5, gp = gpar(col = "black",fontsize=10),just = "left")  

# White square for NA  
#grid.rect(y = 0.35, gp = gpar(fill = "white", col = "black"), width = unit(0.7, "lines"), height = unit(0.7, "lines"))  
#grid.text("NA", x = 0.7, y = 0.35, gp = gpar(col = "black",fontsize=10),just = "left")  

# Pop the viewport  
#popViewport()  


