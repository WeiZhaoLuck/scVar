
library(SoupX)
library(Seurat)
library(DropletUtils)
run_soupx <- function(toc_path, tod_path, all_path, rho = NULL) {
  toc <- Read10X(toc_path, gene.column = 1)
  tod <- Read10X(tod_path, gene.column = 1)

  tod <- tod[rownames(toc), ]

  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 1)
  all <- RunUMAP(all, dims = 1:30)

  matx <- all@meta.data

  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))

  if (is.null(rho)) {
    tryCatch(
      {
        sc <- autoEstCont(sc)
        if (!is.null(sc@contaminationFraction)) {
          # 如果污染比例已设置，则进行调整并输出
          out <- adjustCounts(sc)
          out_file <- file.path(all_path, "soupX_matrix")
          DropletUtils:::write10xCounts(out_file, out, version = "3")
        } else {
          # 污染比例未成功计算，直接输出原始数据
          out_file <- file.path(all_path, "soupX_matrix")
          DropletUtils:::write10xCounts(out_file, tod, version = "3")
        }
      },
      error = function(e) {
        print("autoEstCont Error! Retrying with forceAccept...")
        sc <- autoEstCont(sc, tfidfMin = 0.5, forceAccept = TRUE)
        out <- adjustCounts(sc)
        out_file <- file.path(all_path, "soupX_matrix")
        DropletUtils:::write10xCounts(out_file, out, version = "3")
        # if (!is.null(sc@contaminationFraction)) {
        #   out <- adjustCounts(sc)
        #   out_file <- file.path(all_path, "soupX_matrix")
        #   DropletUtils:::write10xCounts(out_file, out, version = "3")
        # } else {
        #   # PB: 如果污染比例仍然未成功计算，直接输出原始数据
        #   out_file <- file.path(all_path, "soupX_matrix")
        #   DropletUtils:::write10xCounts(out_file, tod, version = "3")
        # }
      }
    )
  } else {
    # 设置污染比例
    sc <- setContaminationFraction(sc, rho)
    out <- adjustCounts(sc)
    out_file <- file.path(all_path, "soupX_matrix")
    DropletUtils:::write10xCounts(out_file, out, version = "3")
  }
}

args <- commandArgs(trailingOnly = TRUE)
all_path=args[1]
print(all_path)
toc_path<-file.path(all_path, "filtered_feature_bc_matrix")
tod_path<-file.path(all_path, "raw_feature_bc_matrix")
run_soupx(toc_path,tod_path,all_path)