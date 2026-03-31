library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)  
library(stringr)
library(maftools)

args <- commandArgs(trailingOnly = TRUE)
main_dir <- args[1]
out_path <- args[2]
pCutoff = as.numeric(args[3])
qCutoff = as.numeric(args[4])

# for (folder_path in list.dirs(main_dir, full.names = TRUE, recursive = FALSE)) { 
#   maf_files <- list.files(folder_path, pattern = "*.maf", full.names = TRUE)
#   # 创建一个空的数据框来存储合并的结果  
#   combined_results_bp <- data.frame() 
#   combined_results_mf <- data.frame()  
#   combined_results_cc <- data.frame()  
  
#   # 循环每个maf文件  
#   for (maf_file in maf_files) {  
#     # 读取maf文件  
#     MAF_info <- read.csv(maf_file, sep = "\t", header = TRUE)  
#     variant_classes <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site","Translation_Start_Site", "Nonsense_Mutation","Nonstop_Mutation", "Missense_Mutation","IGR", "Frameshift_INDEL", "5'UTR", "3'UTR")  
#     filtered_MAF_info <- MAF_info[MAF_info$Variant_Classification %in% variant_classes, ]
#     # 转换基因ID  
#     gene_up_entrez <- tryCatch({  
#       as.character(na.omit(bitr(filtered_MAF_info$Gene_vep,   
#                                 fromType = "ENSEMBL",   
#                                 toType = "ENTREZID",   
#                                 OrgDb = "org.Hs.eg.db")[, 2]))  
#     }, error = function(e) {  
#       message(paste("Error in converting gene IDs for file:", maf_file, ":", e$message))  
#       return(NULL)  # Return NULL if there's an error  
#     })  
    
#     # Check if gene_up_entrez is NULL, if so, skip the rest of the loop  
#     if (is.null(gene_up_entrez) || length(gene_up_entrez) == 20) {  
#       next  # Skip to the next maf_file  
#     }  
    
#     # GO富集分析  
#     enrich_result_mutect <- enrichGO(gene = gene_up_entrez,   
#                                      OrgDb = org.Hs.eg.db,   
#                                      ont = "ALL",   
#                                      pAdjustMethod = "BH",   
#                                      pvalueCutoff = pCutoff,   
#                                      qvalueCutoff = qCutoff)  
#     if (is.null(enrich_result_mutect)) {  
#       next  # Skip to the next maf_file  
#     }
    
#     # 计算GeneRatio_Value  
#     enrich_result_mutect <- enrich_result_mutect %>%   
#       mutate(GeneRatio_Value = as.numeric(str_extract(GeneRatio, "^[^/]+")) /   
#                as.numeric(str_extract(GeneRatio, "(?<=/)[^/]+"))) %>%   
#       mutate(Sample = basename(maf_file)) # 使用文件名作为样本名  
    
#     # 获取top10数据  
#     top10_dat_bp <- enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())
#     top10_dat_cc <- enrich_result_mutect %>% filter(ONTOLOGY == "CC") %>% head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())
#     top10_dat_mf <- enrich_result_mutect %>% filter(ONTOLOGY == "MF") %>% head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())
    
#     # 合并结果  
#     combined_results_bp <- rbind(combined_results_bp, top10_dat_bp)  
#     combined_results_cc <- rbind(combined_results_cc, top10_dat_cc)   
#     combined_results_mf <- rbind(combined_results_mf, top10_dat_mf)   
#   }
#   # 去掉 .maf 后缀  
#   combined_results_bp$Sample <- gsub("\\.maf$", "", combined_results_bp$Sample)  
#   combined_results_cc$Sample <- gsub("\\.maf$", "", combined_results_cc$Sample)  
#   combined_results_mf$Sample <- gsub("\\.maf$", "", combined_results_mf$Sample)
#   file_name_all <- paste0(out_path,basename(folder_path),"/BP.csv") 
#   write.csv(combined_results_bp, file = file_name_all, row.names = FALSE) 
#   file_name_all <- paste0(out_path,basename(folder_path),"/CC.csv")
#   write.csv(combined_results_cc, file =file_name_all, row.names = FALSE) 
#   file_name_all <- paste0(out_path,basename(folder_path),"/MF.csv") 
#   write.csv(combined_results_mf, file = file_name_all, row.names = FALSE)
# }

maf_files <- list.files(main_dir, pattern = "*.maf", full.names = TRUE)
# 创建一个空的数据框来存储合并的结果
combined_results_bp <- data.frame()
combined_results_mf <- data.frame()
combined_results_cc <- data.frame()

combined_onco <- data.frame()
combined_cancer <- data.frame()


# 循环每个maf文件
for (maf_file in maf_files) {
    # 读取maf文件
    MAF_info <- read.csv(maf_file, sep = "\t", header = TRUE)
    variant_classes <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "Missense_Mutation", "IGR", "Frameshift_INDEL", "5'UTR", "3'UTR")
    filtered_MAF_info <- MAF_info[MAF_info$Variant_Classification %in% variant_classes, ]
    # 转换基因ID
    gene_up_entrez <- tryCatch(
        {
            as.character(na.omit(bitr(filtered_MAF_info$Gene_vep,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db"
            )[, 2]))
        },
        error = function(e) {
            message(paste("Error in converting gene IDs for file:", maf_file, ":", e$message))
            return(NULL) # Return NULL if there's an error
        }
    )

    # Check if gene_up_entrez is NULL, if so, skip the rest of the loop
    if (is.null(gene_up_entrez) || length(gene_up_entrez) == 20) {
        next # Skip to the next maf_file
    }

    # GO富集分析
    enrich_result_mutect <- enrichGO(
        gene = gene_up_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = pCutoff,
        qvalueCutoff = qCutoff
    )
    if (is.null(enrich_result_mutect)) {
        next # Skip to the next maf_file
    }

    # 计算GeneRatio_Value
    enrich_result_mutect <- enrich_result_mutect %>%
        mutate(GeneRatio_Value = as.numeric(str_extract(GeneRatio, "^[^/]+")) /
            as.numeric(str_extract(GeneRatio, "(?<=/)[^/]+"))) %>%
        mutate(Sample = basename(maf_file)) # 使用文件名作为样本名

    # 获取top10数据
    top10_dat_bp <- enrich_result_mutect %>%
        filter(ONTOLOGY == "BP") %>%
        head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())
    top10_dat_cc <- enrich_result_mutect %>%
        filter(ONTOLOGY == "CC") %>%
        head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())
    top10_dat_mf <- enrich_result_mutect %>%
        filter(ONTOLOGY == "MF") %>%
        head(enrich_result_mutect %>% filter(ONTOLOGY == "BP") %>% nrow())

    # 合并结果
    combined_results_bp <- rbind(combined_results_bp, top10_dat_bp)
    combined_results_cc <- rbind(combined_results_cc, top10_dat_cc)
    combined_results_mf <- rbind(combined_results_mf, top10_dat_mf)
    laml <- tryCatch(
        {
            read.maf(maf = maf_file, vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "Missense_Mutation", "IGR", "Frameshift_INDEL", "5'UTR", "3'UTR"))
        },
        error = function(e) {
            # 输出错误信息并返回 NULL
            message(paste("Error reading MAF file:", maf_file, " - ", e$message))
            return(NULL)
        }
    )

    pathways <- tryCatch(
        {
            pathways(
                laml,
                pathdb = "sigpw",
                pathways = NULL,
                fontSize = 1,
                panelWidths = c(2, 4, 4),
                plotType = NA,
                col = "#f39c12"
            )
        },
        error = function(e) {
            message(paste("Error processing pathways for file:", maf_file, " - ", e$message))
            return(NULL)
        }
    )
    pathways_cancer <- tryCatch(
        {
            pathways(
                laml,
                pathdb = "smgbp",
                pathways = NULL,
                fontSize = 1,
                panelWidths = c(2, 4, 4),
                plotType = NA,
                col = "#f39c12"
            )
        },
        error = function(e) {
            message(paste("Error processing pathways for file:", maf_file, " - ", e$message))
            return(NULL)
        }
    )

    # 仅在成功获取 pathways 后继续
    if (is.null(pathways)) {
        next # 跳过当前文件，继续下一个文件
    }
    sample_name <- tools::file_path_sans_ext(basename(maf_file))
    pathways_df <- as.data.frame(pathways) %>%
        dplyr::select(Pathway, N, n_affected_genes, fraction_affected) %>%
        mutate(sample = sample_name)
    combined_onco <- rbind(combined_onco, pathways_df)

    # cancer pathways
    if (is.null(pathways_cancer)) {
        next # 跳过当前文件，继续下一个文件
    }
    pathways_cancer_df <- as.data.frame(pathways_cancer) %>%
        mutate(sample = sample_name)
    combined_cancer <- rbind(combined_cancer, pathways_cancer_df)
    

}
# 去掉 .maf 后缀
combined_results_bp$Sample <- gsub("\\.maf$", "", combined_results_bp$Sample)
combined_results_cc$Sample <- gsub("\\.maf$", "", combined_results_cc$Sample)
combined_results_mf$Sample <- gsub("\\.maf$", "", combined_results_mf$Sample)
file_name_all <- paste0(out_path, "/BP.csv")
write.csv(combined_results_bp, file = file_name_all, row.names = FALSE)
file_name_all <- paste0(out_path, "/CC.csv")
write.csv(combined_results_cc, file = file_name_all, row.names = FALSE)
file_name_all <- paste0(out_path, "/MF.csv")
write.csv(combined_results_mf, file = file_name_all, row.names = FALSE)

file_name_all <- paste0(out_path, "/OncogenicPathway.csv")
write.csv(combined_onco, file = file_name_all, row.names = FALSE)

file_name_all <- paste0(out_path, "/PanCancer.csv")
write.csv(combined_cancer, file = file_name_all, row.names = FALSE)
# 'aplot', 'ggnewscale', 'ggtangle', 'scatterpie', 'ggtree' are not available for package 'enrichplot'
# BiocManager::install("clusterProfiler")