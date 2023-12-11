library(Seurat)
library(dplyr)
library(patchwork)
library(SoupX)
library(celldex) ##http://bioconductor.org/packages/release/data/experiment/html/celldex.html
library(scRNAseq) ##https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html
library(SingleR)
library(scater) 
library(scCustomize)
test.data<-Read10X(data.dir = "F:/scvar_docker/example3/E-MTAB-6129/BT1292/mapping/run_count_BT1292/outs/soupX_matrix")
object_test <- CreateSeuratObject(counts = test.data, project = "test", min.cells = 3, min.features = 10)
object_test[["percent.mt"]] <- PercentageFeatureSet(object_test, pattern = "^MT-")
VlnPlot(object_test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object_test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object_test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
rb.genes <- rownames(object_test)[grep("^RP[SL]",rownames(object_test))]
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
object_test[["percent.mito"]] <- PercentageFeatureSet(object_test, pattern = "^RP[SL]")
object_test <- subset(object_test, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & percent.mt < 10 & percent.mito <50)
object_test <- NormalizeData(object = object_test, normalization.method = "LogNormalize", scale.factor = 1e4)
object_test <- FindVariableFeatures(object_test, selection.method = "vst", nfeatures = 2000)
top30 <- head(VariableFeatures(object_test), 30)
FindAllMarkers(object = object_test)
object_seruat<-object_test@assays[["RNA"]]@data
len<-length(colnames(object_test@assays[["RNA"]]@data))
test<-matrix(colnames(object_test@assays[["RNA"]]@data),len,1)
names(test)<-'barcodes'
barcodes_all<-colnames(object_test@assays[["RNA"]]@data)
all.genes <- rownames(object_test)
object_test <- ScaleData(object_test, features = all.genes)
object_test  <- RunPCA(object_test , features = VariableFeatures(object = object_test))
VizDimLoadings(object_test , dims = 1:4,nfeatures = 20,reduction = "pca")
object_test <- FindNeighbors(object_test, dims = 1:30)
object_test <- FindClusters(object_test, resolution = 0.5)
object_test <- RunUMAP(object_test, dims = 1:30)
DimPlot(object_test , reduction = "umap")
object_test  <- RunTSNE(object_test , dims = 1:30)
plot1<-DimPlot(object_test , reduction = "umap",label = TRUE)
plot2<-DimPlot(object_test , reduction = "tsne",label = TRUE)
cluster_mutation=read.csv("F:/github/scvar/mutation_cells_matrix/cluster_result.tsv",sep="\t")
h_cells1=strsplit(gsub(" ","",cluster_mutation$barcode[11]),",")[[1]]
h_cells2=strsplit(gsub(" ","",cluster_mutation$barcode[4]),",")[[1]]
re_1=rep('red', times = length(h_cells1))
re_2=rep('green', times = length(h_cells2))
DimPlot(object_test, cells.highlight = c(h_cells1,h_cells2), cols.highlight = c(re_1, re_2),reduction = "umap",label = TRUE)
DimPlot(object_test, cells.highlight = h_cells1, cols.highlight = "red",reduction = "umap",label = TRUE)
###并列输出多个图
dimplot_list <- list()
for (i in cluster_mutation$cluster){
  h_cells1=strsplit(gsub(" ","",cluster_mutation[cluster_mutation$cluster==i,"barcode"]),",")[[1]]
  p <- DimPlot(object_test, cells.highlight =h_cells1 , cols.highlight = "red",reduction = "umap",label = TRUE)
  dimplot_list[[as.character(i)]] <- p
}
library(cowplot)
cowplot::plot_grid(plotlist = dimplot_list)


cells <- list("cluster6" = h_cells1,
              "cluster4" = h_cells2)

p2 <- Cell_Highlight_Plot(seurat_object = object_test, cells_highlight = cells, highlight_color = c("dodgerblue", "forestgreen"),reduction = "umap",label = TRUE,pt.size="0.2")


colors_test<-c("Unselected"="grey50","Group_2"="red")
p<-DimPlot(object_test, cells.highlight = h_cells1, cols.highlight = "red",reduction = "umap",label = TRUE,cols = colors)
p <- p + guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
saveRDS(object_test,file="test.rds")



##细胞周期
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(object_test))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(object_test))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(object_test))
object_test <- CellCycleScoring(object=object_test,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAa <- RunPCA(object_test, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
p


###特异表达基因
#BiocManager::install("MAST")
markers <- FindAllMarkers(object = object_test, test.use="MAST" ,
                          only.pos = FALSE)   
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
all.markers_neg=all.markers[all.markers$pct.1 < all.markers$pct.2, ]
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10_neg=all.markers_neg %>% group_by(cluster) %>% top_n(n = 10, wt = -avg_log2FC)
all.markers_neg_out <- data.frame(
  gene = all.markers_neg$gene,
  cluster = all.markers_neg$cluster
)
all.markers_neg_out
all.markers_out<- data.frame(
  gene = all.markers$gene,
  cluster = all.markers$cluster
)
write.table(all.markers_neg_out, file = "F:/github/scvar/mutation_cells_matrix/diff_gene.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write.table(all.markers_out, file = "F:/github/scvar/mutation_cells_matrix/diff_gene_all.csv", sep = ",", col.names = TRUE, row.names = FALSE)

###系统发育树
object_test<-BuildClusterTree(object_test)
PlotClusterTree(object_test)

library(maftools)
#laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml = read.maf(maf = "F:/scvar_docker/example3/E-MTAB-6129/BT1298/genotype/merge_final_2.maf",vc_nonSyn=c("RNA","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","5'UTR","3'UTR","Intron","5'Flank","3'Flank","IGR","Frameshift_INDEL","Inframe_INDEL","Silent","Unknown","Translation_Start_Site"))
pdf('F:/scvar_docker/example3/E-MTAB-6129/BT1298/genotype/test.pdf')s
plotmafSummary(maf=laml, rmOutlier=FALSE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
result<-oncoplot(maf=laml,top=10,fontSize = 0.5,showTumorSampleBarcodes = T,legend_height = 1,legendFontSize = 1,writeMatrix=TRUE,SampleNamefontSize = 0.6)
###转换和颠换
titv(maf=laml, plot=TRUE, useSyn=TRUE) s
oncostrip(maf=laml,genes=top30,showTumorSampleBarcodes = T)
###基因互斥与共现
Interact <-somaticInteractions(maf=laml,top=30,pvalue = c(0.05,0.01))
###基因云
#geneCloud(input=laml,top=10)
#rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.6)
###肿瘤突变负荷
tcgaCompare(maf=laml, cohortName="laml")
###VAF
plotVaf(maf=laml,vafCol='VAF',top=20) 
###检测驱动癌症突变
luad.sig <- oncodrive(maf=laml,AACol='aachange',minMut=10, pvalMethod="zscore")
plot1<-plotOncodrive(res = luad.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
plot1
###信号富集分析
OncogenicPathways(maf = laml)
###肿瘤异质性分析
vafclust <-inferHeterogeneity(maf=laml,vafCol='VAF',tsb="Epithelial_cells")
plotClusters(clusters=vafclust)
dev.off()
###肿瘤突变负荷计算
maf_tmb = tmb(maf = laml,captureSize = 50, logScale = TRUE)
