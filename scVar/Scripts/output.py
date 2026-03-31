'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-19 14:47:04
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-25 18:19:36
FilePath: \scvar_docker\codes\output.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import pandas as pd 
import numpy as np
import argparse
import copy 
from scipy import sparse
from scipy.io import mmwrite  
import scanpy as sc  
import mudata as md
import anndata as ad
import pandas as pd  
from concurrent.futures import ProcessPoolExecutor  
from collections import defaultdict  

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
    results=results.fillna(-1)
    sparse_matrix = sparse.csr_matrix(results.values)
    col = results.columns.tolist()
    del results
    return col,sparse_matrix

def read_metadata(chunk):  
    """读取一个数据块的 SNVs 和 barcodes"""  
    return set(chunk['SNV']), set(chunk['barcode'])  

def process_chunk(chunk, value):  
    """处理单个数据块并返回结果字典"""  
    result = defaultdict(dict)  
    for _, row in chunk.iterrows():  
        result[row['SNV']][row['barcode']] = value  
    return result  

def merge_results(acc, new):  
    """合并并行处理结果"""  
    for snv, barcodes in new.items():  
        for barcode, value in barcodes.items():  
            # 保持alt数据优先（值1覆盖值0）  
            if acc.at[snv, barcode] != 1:  
                acc.at[snv, barcode] = value  
    return acc  

def ToMatrix_optimized(ref_path, alt_path, chunk_size=5000, workers=4):  
    # 第一阶段：收集所有SNV和barcode  
    snvs = set()  
    barcodes = set()  
    
    # 并行收集元数据  
    with ProcessPoolExecutor(max_workers=workers) as executor:  
        
        # 处理ref文件元数据  
        futures = []  
        for chunk in pd.read_csv(ref_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):  
            futures.append(executor.submit(read_metadata, chunk))  
        
        # 处理alt文件元数据  
        for chunk in pd.read_csv(alt_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):  
            futures.append(executor.submit(read_metadata, chunk))  
        
        # 合并元数据结果  
        for future in futures:  
            chunk_snvs, chunk_barcodes = future.result()  
            snvs.update(chunk_snvs)  
            barcodes.update(chunk_barcodes)  

    # 创建结果矩阵框架  
    SNV = sorted(snvs)  
    barcodes = sorted(barcodes)  
    results = pd.DataFrame(-1, index=SNV, columns=barcodes)  

    # 第二阶段：并行处理数据  
    # 并行处理ref文件  
    with ProcessPoolExecutor(max_workers=workers) as executor:  
        futures = []  
        for chunk in pd.read_csv(ref_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):  
            futures.append(executor.submit(process_chunk, chunk, 0))  
        
        for future in futures:  
            results = merge_results(results, future.result())  

    # 并行处理alt文件  
    with ProcessPoolExecutor(max_workers=workers) as executor:  
        futures = []  
        for chunk in pd.read_csv(alt_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):  
            futures.append(executor.submit(process_chunk, chunk, 1))  
        
        for future in futures:  
            results = merge_results(results, future.result())  

    return results.columns.tolist(), sparse.csr_matrix(results.values)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ToMatrix')
    parser.add_argument('-r', '--ref', required=True, help='ref_path')
    parser.add_argument('-a', '--alt', required=True, help='alt_path')
    parser.add_argument('-m', '--mutation_information', required=True, help='mutations information path')
    parser.add_argument('-d', '--rna_h5ad', required=True, help='rna information path')
    parser.add_argument('-o', '--out', required=True, help='out_path')
    parser.add_argument('-t', '--threads', required=True, help='threads')
    args = parser.parse_args()
    COLUMN_NAMES = ['SNV', 'barcode', 'chrom', 'position', 'ref', 'alt',   
                'ref_count', 'alt_count', 'all_count', 'celltype']  
    # A,MutationMatrix = ToMatrix(args.ref,args.alt)
    A,MutationMatrix =ToMatrix_optimized(args.ref,args.alt,   
                                 chunk_size=5000, workers=int(args.threads)) 
    adata_combined = sc.read_h5ad(args.rna_h5ad)
    Feature = pd.read_csv(args.mutation_information,sep='\t' ,header=0,index_col=0)
    Feature = Feature.astype(str)
    mutation_matrix = ad.AnnData(X=MutationMatrix.T,obs=pd.DataFrame(index=A),var=Feature)
    del MutationMatrix
    common_cells = adata_combined.obs_names.intersection(mutation_matrix.obs_names) 
    rna_comm = adata_combined[common_cells].copy()  
    mut_comm = mutation_matrix[common_cells].copy() 
    mdata = md.MuData({"rna": rna_comm, "mutation": mut_comm})
    mdata.write(args.out+'/output_data.h5mu')
    
    