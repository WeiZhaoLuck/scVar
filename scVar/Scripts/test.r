library(Seurat)
library(dplyr)
library(patchwork)
library(SoupX)
library(celldex) ##http://bioconductor.org/packages/release/data/experiment/html/celldex.html
library(scRNAseq) ##https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html
library(SingleR)
library(scater) ##https://bioconductor.org/packages/release/bioc/html/scater.html
library(sceasy)
### export PATH=/p300s/biosoft/app/python/python3.9/bin:$PATH
### source activate R4.2.2
### awk -F '\t' '{print $2}' filtered_feature_bc_matrix/features.tsv | paste - soupX_matrix/features.tsv | awk '{$(2+1) = $1; print $0}' | sed 's/[^ ]* //' > test.tsv
library(anndata)
py_require("anndata")

args <- commandArgs(trailingOnly = TRUE)
results_path=args[1]
raw_data_path=args[2]

test.data<-Read10X(data.dir = results_path)
object_test <- CreateSeuratObject(counts = test.data, project = "test", min.cells = 3, min.features = 10)
object_test[["percent.mt"]] <- PercentageFeatureSet(object_test, pattern = "^MT-")
VlnPlot(object_test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object_test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object_test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_out1=file.path(raw_data_path, "percentmt_nFeature_RNA.png")
ggsave(plot_out1,plot1 + plot2,dpi = 300,width = 20,height = 10)
rb.genes <- rownames(object_test)[grep("^RP[SL]",rownames(object_test))]
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
object_test[["percent.mito"]] <- PercentageFeatureSet(object_test, pattern = "^RP[SL]")
object_test <- subset(object_test, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & percent.mt < 10 & percent.mito <50)
object_test <- NormalizeData(object = object_test, normalization.method = "LogNormalize", scale.factor = 1e4)
object_test <- FindVariableFeatures(object_test, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(object_test), 10)
object_seruat<-object_test@assays[["RNA"]]@data
csv_out1=file.path(raw_data_path, "result_test2.csv")
#write.csv(object_test@assays[["RNA"]]@data,csv_out1,qmethod="double")

len<-length(colnames(object_test@assays[["RNA"]]@data))
test<-matrix(colnames(object_test@assays[["RNA"]]@data),len,1)
names(test)<-'barcodes'
csv_out2=file.path(raw_data_path, "result_barcodes.csv")
write.csv(test,csv_out2)
barcodes_all<-colnames(object_test@assays[["RNA"]]@data)
barcode_out=file.path(raw_data_path, "barcodes.csv")
write.table(barcodes_all,barcode_out, quote = FALSE, row.names = FALSE,col.names = FALSE)

all.genes <- rownames(object_test)
object_test <- ScaleData(object_test, features = all.genes)
object_test  <- RunPCA(object_test , features = VariableFeatures(object = object_test))
VizDimLoadings(object_test , dims = 1:2,nfeatures = 20,reduction = "pca")
object_test <- FindNeighbors(object_test, dims = 1:30)
object_test <- FindClusters(object_test, resolution = 0.5)
object_test <- RunUMAP(object_test, dims = 1:30)
DimPlot(object_test , reduction = "umap")
object_test  <- RunTSNE(object_test , dims = 1:30)
plot1<-DimPlot(object_test , reduction = "umap",label = TRUE)
plot2<-DimPlot(object_test , reduction = "tsne",label = TRUE)
plot_out2=file.path(raw_data_path, "umap_tsne.png")
ggsave(plot_out2,plot1 + plot2,dpi = 300,width = 20,height = 10)
meta=object_test@meta.data
hpca.se <-  get(load("/reference/ref_Human_all.RData"))
object_test_for_SingleR <- GetAssayData(object_test, slot="data")
object_test.hesc <- SingleR(test = object_test_for_SingleR , ref = hpca.se, labels = hpca.se$label.main)
table(object_test.hesc$labels,meta$seurat_clusters)
object_test@meta.data$labels <-object_test.hesc$labels
plot3<-DimPlot(object_test, group.by = c("seurat_clusters", "labels"),reduction = "umap")
#print(DimPlot(object_test, group.by = c("seurat_clusters", "labels"),reduction = "umap"))
plot_out3=file.path(raw_data_path, "cell_type.png")
ggsave(plot_out3,plot3,dpi = 300,width = 20,height = 10)
csv_out3=file.path(raw_data_path, "barcodes_type.csv")
write.csv(object_test@meta.data,csv_out3)
# 提取metadata的label列
label_data <- object_test@meta.data$label

# 创建新的数据框, 将行名作为cell列
extracted_data <- data.frame(cell = rownames(object_test@meta.data), celltype = label_data)

csv_out3 <- file.path(raw_data_path, "output_cell_barcodes.tsv")
write.table(extracted_data, csv_out3, sep = "\t", quote = F, row.names = F)

cluster_data <- object_test@meta.data$seurat_clusters
extracted_data_cluster <- data.frame(cell = rownames(object_test@meta.data), cluster = cluster_data)
csv_out3 <- file.path(raw_data_path, "output_cell_cluster.tsv")
write.table(extracted_data_cluster, csv_out3, sep = "\t", quote = F, row.names = F)

out_rds <- file.path(raw_data_path, "rna.rds")
saveRDS(object_seruat, out_rds)

out_h5ad <- file.path(raw_data_path, "rna.h5ad")
sceasy::convertFormat(object_test, from = "seurat", to = "anndata", outFile = out_h5ad)

