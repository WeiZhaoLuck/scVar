import pandas as pd
import os
import sys
import numpy as np
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
# from scipy.stats.contingency import odds_ratio
import sys


def each_group(labels1,labels2):
    data_labels_alt=all_alt
    data_labels_ref=all_ref
    resultsP=[]
    results_odd=[]
    data_labels_ref_out=[]
    SNV=data_labels_ref.index.values.tolist()
    for a in range(len(data_labels_ref)):
        print(a)
        if(data_labels_alt.iloc[a,labels2]>0 and data_labels_ref.iloc[a,labels1]>0):
            data_each_ref=[data_labels_ref.iloc[a,labels1],data_labels_ref.iloc[a,labels2]]
            data_each_alt=[data_labels_alt.iloc[a,labels1],data_labels_alt.iloc[a,labels2]]
            res,result_p=fisher_exact([data_each_alt,data_each_ref])
            # res = odds_ratio([data_each_ref,data_each_alt]).statistic
        #     result_all[result_all['SNV_label']==data_labels_ref.iloc[i,0]]['Pvalue']=result_p
            if(result_p<0.05):
                resultsP.append(result_p)
                results_odd.append(res)
                data_labels_ref_out.append(SNV[a])
    if(len(resultsP)>0):
        reject,p_adjusted,_,_=smm.multipletests(resultsP,method='fdr_bh')

        # SNV=data_labels_ref.iloc[:,0].tolist()
        result_all=pd.DataFrame()
        result_all['SNV_label']=data_labels_ref_out
        result_all['p-adjusted']=p_adjusted
        result_all['pvalue']=resultsP
        result_all['odd_ratio']=results_odd
        result_all['ref_label1']=data_labels_ref.loc[data_labels_ref_out].iloc[:,labels1].tolist()
        result_all['ref_label2']=data_labels_ref.loc[data_labels_ref_out].iloc[:,labels2].tolist()
        result_all['alt_label1']=data_labels_alt.loc[data_labels_ref_out].iloc[:,labels1].tolist()
        result_all['alt_label2']=data_labels_alt.loc[data_labels_ref_out].iloc[:,labels2].tolist()
        # results_gene_SNV=pd.merge(result_all,snv_info,on=['SNV_label'],how="right")
        # results_gene_SNV_sorted=results_gene_SNV.sort_values('p-adjusted')
        # results_gene_SNV_sorted_out=results_gene_SNV_sorted[results_gene_SNV_sorted['p-adjusted']<0.05]
        results_gene_SNV=result_all[result_all['pvalue']<0.05]
        name="padjusted"+"_"+labels_all_ref[labels1]+"_"+labels_all_alt[labels2]+".csv"
        # name2="pvalue"+labels+".csv"
        path_out=os.path.join("results",name)
        # path_out2=os.path.join("/mnt/f/github/scvar/cell_specify_mutations_analysis/results","summary.report")
        # path_out2=os.path.join("/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/SNV_celltype/files",name2)
        print(result_all)
        # if(len(result_all)>0):
        #     result_all.to_csv(path_out)        
        results_gene_SNV.to_csv(path_out)

def cluster_other(label):
    data_labels_alt=all_alt
    data_labels_ref=all_ref
    resultsP=[]
    results_odd=[]
    data_labels_ref_out=[]
    SNV=data_labels_ref.index.values.tolist()
    result_list = [element for element in labels_all_ref if element != label]
    data_labels_ref['RowSum'] = data_labels_ref[result_list].sum(axis=1)
    data_labels_alt['RowSum'] = data_labels_alt[result_list].sum(axis=1)
    for a in range(len(data_labels_ref)):
        name=SNV[a]
        if(data_labels_alt.loc[name,'RowSum']>0 and data_labels_ref.loc[name,'RowSum']>0):
            data_each_ref=[data_labels_ref.loc[name,label],data_labels_ref.loc[name,'RowSum']]
            data_each_alt=[data_labels_alt.loc[name,label],data_labels_alt.loc[name,'RowSum']]
            res,result_p=fisher_exact([data_each_alt,data_each_ref])
            # res = odds_ratio([data_each_ref,data_each_alt]).statistic
        #     result_all[result_all['SNV_label']==data_labels_ref.iloc[i,0]]['Pvalue']=result_p
            resultsP.append(result_p)
            results_odd.append(res)
            data_labels_ref_out.append(SNV[a])
    reject,p_adjusted,_,_=smm.multipletests(resultsP,method='fdr_bh')

    # SNV=data_labels_ref.iloc[:,0].tolist()
    result_all=pd.DataFrame()
    result_all['SNV_label']=data_labels_ref_out
    result_all['p-adjusted']=p_adjusted
    result_all['pvalue']=resultsP
    result_all['odd_ratio']=results_odd
    result_all['ref_label1']=data_labels_ref.loc[data_labels_ref_out].loc[:,label].tolist()
    result_all['ref_label2']=data_labels_ref.loc[data_labels_ref_out].loc[:,'RowSum'].tolist()
    result_all['alt_label1']=data_labels_alt.loc[data_labels_ref_out].loc[:,label].tolist()
    result_all['alt_label2']=data_labels_alt.loc[data_labels_ref_out].loc[:,'RowSum'].tolist()
    # results_gene_SNV=pd.merge(result_all,snv_info,on=['SNV_label'],how="right")
    # results_gene_SNV_sorted=results_gene_SNV.sort_values('p-adjusted')
    # results_gene_SNV_sorted_out=results_gene_SNV_sorted[results_gene_SNV_sorted['p-adjusted']<0.05]
    results_gene_SNV=result_all[result_all['pvalue']<0.05]
    name="padjusted"+"_"+label+"_"+'other'+".csv"
    # name2="pvalue"+labels+".csv"
    path_out=os.path.join("results",name)
    # path_out2=os.path.join("/mnt/f/github/scvar/cell_specify_mutations_analysis/results","summary.report")
    # path_out2=os.path.join("/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/SNV_celltype/files",name2)
    print(result_all)
    # if(len(result_all)>0):
    #     result_all.to_csv(path_out)        
    result_all.to_csv(path_out)
    
    


if __name__ == '__main__':
    all_ref=pd.read_csv("all_table_ref_gene_cluster.csv",sep=",",index_col=0)
    all_alt=pd.read_csv("all_table_alt_gene_cluster.csv",sep=",",index_col=0)
    # snv_info=pd.read_csv(sys.argv[3],header=0,sep="\t")
    # snv_info.columns=['SNV_label','SYMBOL_vep','Gene_vep']
    # each_group(all_ref,all_alt)
    labels_all_ref=all_ref.columns.tolist()[:]
    # print(len(labels_all_ref))
    labels_all_alt = all_alt.columns.tolist()[:]
    print (len(labels_all_alt))
    # for i in range(len(labels_all_ref)):
    #     for l in range(len(labels_all_alt)):
    #         if labels_all_ref[i]==labels_all_alt[l]:
    #             continue
    #         else:
    #             each_group(i,l)
    for i in range(len(labels_all_ref)):
        cluster_other(labels_all_ref[i])