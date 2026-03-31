import pandas as pd 
import numpy as np
import argparse
import copy

def ToMatrix(ref_path,alt_path):
    all_ref=pd.read_csv(ref_path,sep="\t",header=None,names=['SNV', 'barcode', 'chrom','position','ref','alt','ref_count','alt_count','all_count','celltype'])
    all_alt=pd.read_csv(alt_path,sep="\t",header=None,names=['SNV', 'barcode', 'chrom','position','ref','alt','ref_count','alt_count','all_count','celltype'])
    SNV=sorted(list(set(all_alt['SNV'].tolist())))
    barcodes =sorted(list(set(all_alt['barcode']).union(set(all_ref['barcode']))))
    results=pd.DataFrame(data=None,index=SNV, columns=barcodes)
    for i in range(len(barcodes)):
        bar=all_ref[all_ref['barcode']==barcodes[i]]['SNV'].values.tolist()
        bar_o=list(set(bar).intersection(set(SNV)))
        results.loc[bar_o,barcodes[i]]=0
        bar_alt=all_alt[all_alt['barcode']==barcodes[i]]['SNV'].values.tolist()
        bar_o=list(set(bar_alt).intersection(set(SNV)))
        results.loc[bar_o,barcodes[i]]=1
    results=results.fillna(0)
    final=args.out+"mutations_matrix.csv"
    results.to_csv(final,sep=",")
    mutation_info=pd.read_csv(args.mutation_information,sep="\t")
    all_gene=list(set(mutation_info['Gene_vep']))
    all_cells=results.columns.tolist()
    gene_cells=pd.DataFrame(0,index=all_gene,columns=all_cells)
    for i in range(len(all_gene)):
        mutations_each_gene=mutation_info[mutation_info['Gene_vep']==all_gene[i]]['CHROM_POS_REF_ALT'].tolist()
        intersection=list(set(results.index.tolist()).intersection(set(mutations_each_gene)))
        sum_result = results.loc[intersection].sum()
        gene_cells.loc[all_gene[i]]=sum_result
    gene_cells[gene_cells > 0] = 1
    results_gene=copy.deepcopy(gene_cells)
    final_gene=args.out+"gene_matrix.csv"
    results_gene.to_csv(final_gene,sep=",")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ToMatrix')
    parser.add_argument('-r', '--ref', required=True, help='ref_path')
    parser.add_argument('-a', '--alt', required=True, help='alt_path')
    parser.add_argument('-m', '--mutation_information', required=True, help='mutations information path')
    parser.add_argument('-o', '--out', required=True, help='out_path')
    
    args = parser.parse_args()
    ToMatrix(args.ref,args.alt)     
        

    