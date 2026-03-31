library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(scRNAseq)
library(SingleR)
library(scater)

args <- commandArgs(T)

top20_cell=read.csv(args[1],sep="\t",header=FALSE)
print(args[1])
print(args[2])
print(args[3])
test.data<-Read10X(data.dir = args[2])
object_test <- CreateSeuratObject(counts = test.data, project = "test", min.cells = 3, min.features = 10)
object_test[["percent.mt"]] <- PercentageFeatureSet(object_test, pattern = "^MT-")
rb.genes <- rownames(object_test)[grep("^RP[SL]",rownames(object_test))]
object_test[["percent.mito"]] <- PercentageFeatureSet(object_test, pattern = "^RP[SL]")
object_test <- subset(object_test, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & percent.mt < 10 & percent.mito <50)
object_test <- NormalizeData(object = object_test, normalization.method = "LogNormalize", scale.factor = 1e4)
object_test <- FindVariableFeatures(object_test, selection.method = "vst", nfeatures = 2000)
object_seruat<-object_test@assays[["RNA"]]@data

len<-length(colnames(object_test@assays[["RNA"]]@data))
test<-matrix(colnames(object_test@assays[["RNA"]]@data),len,1)
names(test)<-'barcodes'

all.genes <- rownames(object_test)
object_test <- ScaleData(object_test, features = all.genes)
object_test  <- RunPCA(object_test , features = VariableFeatures(object = object_test))
VizDimLoadings(object_test , dims = 1:2,nfeatures = 20,reduction = "pca")
object_test <- FindNeighbors(object_test, dims = 1:30)
object_test <- FindClusters(object_test, resolution = 0.5)
object_test <- RunUMAP(object_test, dims = 1:30)
DimPlot(object_test , reduction = "umap")
object_test  <- RunTSNE(object_test , dims = 1:30)


meta=object_test@meta.data
hpca.se <- get(load("/codes/ref_Human_all.RData"))
object_test_for_SingleR <- GetAssayData(object_test, slot="data")
object_test.hesc <- SingleR(test = object_test_for_SingleR , ref = hpca.se, labels = hpca.se$label.main)
table(object_test.hesc$labels,meta$seurat_clusters)
object_test@meta.data$labels <-object_test.hesc$labels
plot5<-print(DimPlot(object_test, group.by = c("labels"),reduction = "tsne"))
for (i in 1:length(top20_cell$V2)){
  h_cells<-strsplit(gsub(" ","",top20_cell$V3[i]),",")[[1]]
  name<-top20_cell$V2[i]
  plot6<-DimPlot(object_test, cells.highlight = h_cells, cols.highlight = "red",reduction = "tsne",label = TRUE)
  plot_all<-plot5+plot6
  # raw_data_path=paste("/results/E-MTAB-8410/ERS15977681/report/data","top_mutation_tsne",sep="/")
  raw_data_path=paste(args[3],"top_mutation_tsne",sep="/")
  # raw_data_path="D:/E_extension/scvar/pipeline/top20_cell_highlight/top20_mutation_tsne"
  name_out<-paste(name,"png",sep=".")
  plot_out=file.path(raw_data_path,name_out)
  ggsave(plot_out,plot_all,dpi = 300,width = 20,height = 10)
}
plot5<-print(DimPlot(object_test, group.by = c( "labels"),reduction = "umap"))
for (i in 1:length(top20_cell$V2)){
  h_cells<-strsplit(gsub(" ","",top20_cell$V3[i]),",")[[1]]
  name<-top20_cell$V2[i]
  plot6<-DimPlot(object_test, cells.highlight = h_cells, cols.highlight = "red",reduction = "umap",label = TRUE)
  plot_all<-plot5+plot6
  raw_data_path=paste(args[3],"top_mutation_umap",sep="/")
  name_out<-paste(name,"png",sep=".")
  plot_out=file.path(raw_data_path,name_out)
  ggsave(plot_out,plot_all,dpi = 300,width = 20,height = 10)
}




