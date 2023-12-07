import pandas as pd
import os
import sys
import numpy as np
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
# from scipy.stats.contingency import odds_ratio
import sys


def each_group(labels1,labels2):
#
# labels1=sys.argv[1]
# labels2=sys.argv[2]
# barcode_path="/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/snakemake/test/E-MTAB-8410/ERS3858524/seurat/barcodes_type.csv"
# barcode_dall=pd.read_csv(barcode_path,sep=",")
# barcodes=barcode_dall.iloc[:,0]
#
# genotype_results1="/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/snakemake/test/E-MTAB-8410/ERS3858524/SNV_celltype_annotation.tsv"
# genotype_results=pd.read_csv(genotype_results1,sep="\t")

    # genotype_results2="/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/SNV_celltype/files/barcode_dall.csv"
    # barcode_dall=pd.read_csv(genotype_results2,sep=",")
    # df=barcode_dall[['gene','SNV_label']]
    # gene_SNV=df.drop_duplicates(subset=['gene','SNV_label'],keep='first')


    data_labels_alt=all_alt
    data_labels_ref=all_ref
    resultsP=[]
    results_odd=[]
    for i in range(len(data_labels_ref)):
        data_each_ref=[data_labels_ref.loc[i,labels1],data_labels_ref.loc[i,labels2]]
        data_each_alt=[data_labels_alt.loc[i,labels1],data_labels_alt.loc[i,labels2]]
        res,result_p=fisher_exact([data_each_ref,data_each_alt])
        # res = odds_ratio([data_each_ref,data_each_alt]).statistic
    #     result_all[result_all['SNV_label']==data_labels_ref.iloc[i,0]]['Pvalue']=result_p
        resultsP.append(result_p)
        results_odd.append(res)

    reject,p_adjusted,_,_=smm.multipletests(resultsP,method='fdr_bh')

    SNV=data_labels_ref.iloc[:,0].tolist()
    result_all=pd.DataFrame()
    result_all['SNV_label']=SNV
    result_all['p-adjusted']=p_adjusted
    result_all['pvalue']=resultsP
    result_all['odd_ratio']=results_odd
    result_all['ref_label1']=data_labels_ref.loc[:,labels1]
    result_all['ref_label2']=data_labels_ref.loc[:,labels2]
    result_all['alt_label1']=data_labels_alt.loc[:,labels1]
    result_all['alt_label2']=data_labels_alt.loc[:,labels2]
    results_gene_SNV=pd.merge(result_all,snv_info,on=['SNV_label'],how="right")
    # results_gene_SNV_sorted=results_gene_SNV.sort_values('p-adjusted')
    # results_gene_SNV_sorted_out=results_gene_SNV_sorted[results_gene_SNV_sorted['p-adjusted']<0.05]
    # results_gene_SNV=results_gene_SNV[results_gene_SNV['p-adjusted']<0.05]
    name="padjusted"+"_"+labels1+"_"+labels2+".csv"
    # name2="pvalue"+labels+".csv"
    path_out=os.path.join(sys.argv[4],name)
    # path_out2=os.path.join("/mnt/f/github/scvar/cell_specify_mutations_analysis/results","summary.report")
    # path_out2=os.path.join("/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/SNV_celltype/files",name2)
    print(results_gene_SNV)
    # if(len(result_all)>0):
    #     result_all.to_csv(path_out)        
    results_gene_SNV.to_csv(path_out)
    # with open(path_out2, 'w+') as f:
    #     title=labels1+" VS "+labels2
    #     f.write(title)        
    #     f.write(result_all)

if __name__ == '__main__':
    all_ref=pd.read_csv(sys.argv[1],sep=",")
    all_alt=pd.read_csv(sys.argv[2],sep=",")
    snv_info=pd.read_csv(sys.argv[3],header=0,sep="\t")
    snv_info.columns=['SNV_label','SYMBOL_vep','Gene_vep']
    # each_group(all_ref,all_alt)
    labels_all_ref=all_ref.columns.tolist()[1:]
    # print(len(labels_all_ref))
    labels_all_alt = all_alt.columns.tolist()[1:]
    print (len(labels_all_alt))
    for i in range(len(labels_all_ref)):
        for l in range(len(labels_all_alt)):
            if labels_all_ref[i]==labels_all_alt[l]:
                continue
            else:
                each_group(labels_all_ref[i],labels_all_alt[l])