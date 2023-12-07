import pandas as pd
import sys
from collections import Counter

def merge_df1_df2(df1,df2,key):
    merged_df = pd.merge(df1, df2, on=key, how='left')
    return merged_df

genotype_path="BT1292_1/genotype_alt_cell_type.tsv"
ref_path="BT1292_1/genotype_ref_cell_type.tsv"
genotype_results=pd.read_table(genotype_path,sep="\t",header=None,names=['key','barcode','chrom','pos','ref','alt','ref_num','alt_num','all_num','cell_type'])
ref_results=pd.read_table(ref_path,sep="\t",header=None,names=['key','barcode','chrom','pos','ref','alt','ref_num','alt_num','all_num','cell_type'])
result_final_alt = genotype_results.groupby('key')['barcode'].apply(list).reset_index()
result_final_ref = ref_results.groupby('key')['barcode'].apply(list).reset_index()
merge_result=merge_df1_df2(result_final_alt,result_final_ref,'key')

mutation_Gene=pd.read_csv("BT1292_1/mutaion.info",sep="\t",header=0,names=['mutation','gene'],index_col=False)
barcode_cluster_info=pd.read_csv("BT1292_1/barcodes_cell_type.tsv",sep="\t",header=None,names=['barcode','cluster'],index_col=False)
unique_df = mutation_Gene.drop_duplicates()
all_gene=list(set(unique_df['gene']))
cell_type=list(set(barcode_cluster_info['cluster'].values.tolist()))
out_data_alt=pd.DataFrame(0,index=all_gene,columns=cell_type)
out_data_ref=pd.DataFrame(0,index=all_gene,columns=cell_type)
for i in range(len(all_gene)):
    print(i)
    mutations_each_gene=unique_df[unique_df['gene']==all_gene[i]]['mutation'].tolist()
    intersection = list(set(merge_result['key'].tolist()).intersection(set(mutations_each_gene)))
    if(len(intersection)>0):
        mutations_each_gene_data=merge_result[merge_result['key'].isin(intersection)]
        mutations_each_gene_data_sum_alt=list(set(mutations_each_gene_data.explode('barcode_x')['barcode_x'].tolist()))
        sum_alt=Counter(barcode_cluster_info[barcode_cluster_info['barcode'].isin(mutations_each_gene_data_sum_alt)]['cluster'])
        for column_name, column_data in sum_alt.items():
            out_data_alt.loc[all_gene[i],column_name] = column_data
        mutations_each_gene_data_sum_ref=list(set(mutations_each_gene_data.explode('barcode_y')['barcode_y'].tolist()))
        sum_ref=Counter(barcode_cluster_info[barcode_cluster_info['barcode'].isin(mutations_each_gene_data_sum_ref)]['cluster'])
        for column_name, column_data in sum_ref.items():
            out_data_ref.loc[all_gene[i],column_name] = column_data
out_data_alt.to_csv("BT1292_1/all_table_alt_gene_cell_type.csv",sep=",")
out_data_ref.to_csv("BT1292_1/all_table_ref_gene_cell_type.csv",sep=",")
# SNVs=set(genotype_results['key'].tolist())
# SNVs_celltype=pd.DataFrame(data=SNVs,columns=['SNV'])
# result_final = genotype_results.groupby('key')['barcode'].apply(list).reset_index()
# SNVs_celltype['mutation_barcodes']=result_final['barcode']

