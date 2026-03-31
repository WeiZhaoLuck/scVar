import pandas as pd 
import numpy as np
import argparse
import copy 
from scipy import sparse
from scipy.io import mmwrite  
import scanpy as sc  
import mudata as md
import anndata as ad
import math

def filter_rows(df):
    filtered_rows = []  
    for index, row in df.iterrows():  
        total_cells = len(row)  # 细胞总数  
        count_neg1 = (row == -1).sum()  # 计算 -1 的数量  
        
        # 计算 -1 占总数的比例  
        negative_ratio = count_neg1 / total_cells  
        
        # 判断是否小于 95%  
        if negative_ratio <= 0.95:  
            filtered_rows.append(row)  
    
    # 创建过滤后的数据框  
    filtered_df = pd.DataFrame(filtered_rows, index=df.index[:len(filtered_rows)])
    return filtered_df

def calculate_entropy_simpson(df):  
    results = []  
    # 获取突变ID（行名）  
    mutation_ids = df.index.tolist()  
    
    for idx, row in df.iterrows():  
        # 过滤掉未覆盖的数据（-1）  
        effective = [x for x in row if x != -1]  
        n0 = effective.count(0)  
        n1 = effective.count(1)  
        N = len(effective)  
        
        # 计算信息熵  
        if N == 0:  
            entropy = math.nan  
        else:  
            p0 = n0 / N  
            p1 = n1 / N  
            entropy_val = 0.0  
            if p0 > 0:  
                entropy_val += p0 * math.log2(p0)  
            if p1 > 0:  
                entropy_val += p1 * math.log2(p1)  
            entropy = -entropy_val  
        
        # 计算辛普森指数  
        if N < 2:  
            simpson = math.nan  
        else:  
            numerator = n0 * (n0 - 1) + n1 * (n1 - 1)  
            denominator = N * (N - 1)  
            simpson = numerator / denominator  
        
        # 将突变ID、信息熵和辛普森指数添加到结果中   
        print(idx, entropy, simpson)  
        results.append({'mutation_id': idx, 'entropy': entropy, 'simpson': simpson})  
    
    # 创建结果数据框  
    results_df = pd.DataFrame(results)  
    return results_df 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculating Entropy and Simpsons Index')
    parser.add_argument('-i', '--input', required=True, help='h5mu_path')
    parser.add_argument('-o', '--out', required=True, help='out_path')
    args = parser.parse_args()
    Mudata = md.read(args.input)
    df = pd.DataFrame.sparse.from_spmatrix(Mudata['mutation'].X, index=Mudata['mutation'].obs.index, columns=Mudata['mutation'].var.index).T
    results_df = calculate_entropy_simpson(filter_rows(df))
    results_df.to_csv(args.out+'./entropy_simpson_results.tsv', sep='\t', index=False)
    