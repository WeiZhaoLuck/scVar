'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-26 18:02:51
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-26 18:02:59
FilePath: \scvar_docker\codes\output2.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import pandas as pd
import numpy as np
import argparse
from scipy import sparse
import scanpy as sc
import mudata as md
import anndata as ad
from collections import defaultdict

COLUMN_NAMES = ['SNV', 'barcode', 'chrom', 'position', 'ref', 'alt',
                'ref_count', 'alt_count', 'all_count', 'celltype']

def collect_metadata(file_path, chunk_size=5000):
    """收集文件中的所有SNV和barcode"""
    snvs = set()
    barcodes = set()
    for chunk in pd.read_csv(file_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):
        snvs.update(chunk['SNV'].unique())
        barcodes.update(chunk['barcode'].unique())
    return snvs, barcodes

def process_chunk(chunk, snv_to_idx, barcode_to_idx, value, entries):
    """处理数据块并更新entries字典"""
    for _, row in chunk.iterrows():
        snv = row['SNV']
        barcode = row['barcode']
        if snv in snv_to_idx and barcode in barcode_to_idx:
            r = snv_to_idx[snv]
            c = barcode_to_idx[barcode]
            entries[(r, c)] = value

def ToMatrix_optimized(ref_path, alt_path, chunk_size=5000):
    # 阶段一：收集所有SNV和barcode
    ref_snvs, ref_barcodes = collect_metadata(ref_path, chunk_size)
    alt_snvs, alt_barcodes = collect_metadata(alt_path, chunk_size)
    
    all_snvs = ref_snvs.union(alt_snvs)
    all_barcodes = ref_barcodes.union(alt_barcodes)
    
    # 创建索引映射
    snv_list = sorted(all_snvs)
    barcode_list = sorted(all_barcodes)
    snv_to_idx = {snv: i for i, snv in enumerate(snv_list)}
    barcode_to_idx = {barcode: j for j, barcode in enumerate(barcode_list)}
    
    # 阶段二：处理数据
    entries = defaultdict(int)
    
    # 处理ref文件
    for chunk in pd.read_csv(ref_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):
        process_chunk(chunk, snv_to_idx, barcode_to_idx, 0, entries)
    
    # 处理alt文件（覆盖ref的值）
    for chunk in pd.read_csv(alt_path, sep="\t", header=None, names=COLUMN_NAMES, chunksize=chunk_size):
        process_chunk(chunk, snv_to_idx, barcode_to_idx, 1, entries)
    
    # 生成稀疏矩阵
    rows, cols, data = [], [], []
    for (r, c), val in entries.items():
        rows.append(r)
        cols.append(c)
        data.append(val)
    
    coo_matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(snv_list), len(barcode_list)))
    csr_matrix = coo_matrix.tocsr()
    
    return barcode_list, csr_matrix, snv_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate mutation matrix')
    parser.add_argument('-r', '--ref', required=True, help='Path to ref counts file')
    parser.add_argument('-a', '--alt', required=True, help='Path to alt counts file')
    parser.add_argument('-m', '--mutation_info', required=True, help='Path to mutation information')
    parser.add_argument('-d', '--rna_h5ad', required=True, help='Path to RNA h5ad file')
    parser.add_argument('-o', '--out', required=True, help='Output directory')
    args = parser.parse_args()

    # 生成突变矩阵
    barcodes, mutation_matrix, snvs = ToMatrix_optimized(args.ref, args.alt)
    
    # 读取RNA数据和突变信息
    adata_combined = sc.read_h5ad(args.rna_h5ad)
    feature_df = pd.read_csv(args.mutation_info, sep='\t', index_col=0).astype(str)
    
    # 创建AnnData对象并对齐细胞
    mutation_adata = ad.AnnData(
        X=mutation_matrix.T,
        obs=pd.DataFrame(index=barcodes),
        var=feature_df.loc[snvs]  # 确保顺序一致
    )
    
    # 对齐细胞
    common_cells = adata_combined.obs_names.intersection(mutation_adata.obs_names)
    rna_subset = adata_combined[common_cells].copy()
    mutation_subset = mutation_adata[common_cells].copy()
    
    # 创建MuData对象并保存
    mdata = md.MuData({"rna": rna_subset, "mutation": mutation_subset})
    mdata.write(f"{args.out}/output_data.h5mu")