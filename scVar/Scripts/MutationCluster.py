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

import numpy as np  
import pandas as pd  
from scipy.sparse import csr_matrix, csc_matrix  
from sklearn.feature_extraction.text import TfidfTransformer  
import subprocess
import os
from scipy.sparse import csr_matrix, issparse

def set_workspace(path):  
    # 如果文件夹不存在，则创建它  
    if not os.path.exists(path):  
        os.makedirs(path)  
        print(f"文件夹已创建: {path}")  
    else:  
        print(f"文件夹已存在: {path}")  

    # 设置工作空间  
    os.chdir(path)  
    print(f"工作空间已设置为: {os.getcwd()}") 
    

def filter_mutations_chunked(sparse_matrix,mutations, cell_threshold=5, chunk_size=1000):
    
    n_mutations, n_cells = sparse_matrix.shape  
    min_cells = cell_threshold  
    print(f"筛选条件：每个突变至少需要在 {min_cells} 个细胞中不为-1")  

 
    valid_rows = []  
    index_rows= []

    # 分块处理  
    for chunk_start in range(0, n_mutations, chunk_size):  
        chunk_end = min(chunk_start + chunk_size, n_mutations)  
        chunk = sparse_matrix[chunk_start:chunk_end]  
        
        # 转换为COO格式便于逐行处理  
        coo_chunk = chunk.tocoo()  
        
        # 计算每行非-1的细胞数  
        row_counts = np.bincount(  coo_chunk.row,  weights=(coo_chunk.data != -1).astype(int),  minlength=chunk.shape[0])  
        
        # 筛选符合条件的行（局部索引）  
        chunk_valid = np.where(row_counts >= min_cells)[0]  
        
        # 转换为全局索引  
        global_valid = chunk_start + chunk_valid  
        index_rows.extend(global_valid)
        result = [mutations[i] for i in global_valid.tolist()] 
        valid_rows.extend(result)  
        
        print(f"处理块 [{chunk_start}-{chunk_end}]，发现 {len(chunk_valid)} 个有效突变")  

    # 提取最终结果  
    print(f"共筛选出 {len(valid_rows)} 个符合要求的突变")  
    return valid_rows,sparse_matrix[index_rows]

# def process_data(option, sparse_matrix, row_names, col_names, treat_neg_one_as_zero):  
#     """  
#     处理稀疏矩阵并生成特征矩阵，支持 pandas.DataFrame 或 scipy.sparse 矩阵。  

#     参数:  
#     - option: 选择处理方法 ('var1', 'var2', 'TF_IDF')  
#     - sparse_matrix: 输入的稀疏矩阵 (pandas.DataFrame 或 scipy.sparse 矩阵)  
#     - treat_neg_one_as_zero: 是否将 -1 视为 0 (0: 不视为 0, 1: 视为 0)  
#     """  

#     # # 确保输入为 scipy.sparse 矩阵  
#     # if isinstance(sparse_matrix, pd.DataFrame): 
#     #     row_names = df.index
#     #     col_names = df.columns
#     #     sparse_matrix = csr_matrix(sparse_matrix.sparse.to_coo())  # 转为 Scipy 稀疏矩阵格式  
#     # elif not isinstance(sparse_matrix, (csr_matrix, csc_matrix)):  
#     #     raise ValueError("Input sparse_matrix must be a pandas.DataFrame with sparse format or a scipy.sparse matrix.")  

    
    
#     # 如果 treat_neg_one_as_zero 为 1，将 -1 替换为 0  
#     if treat_neg_one_as_zero == 1:  
#         sparse_matrix = sparse_matrix.copy()  
#         sparse_matrix.data[sparse_matrix.data == -1] = 0  # 仅替换非零值中的 -1  

#     # 方法1: 处理 var1  
#     if option == 'var1':  
#         # 使用稀疏矩阵计算每行方差  
#         mean = sparse_matrix.mean(axis=1).A1  # 稀疏矩阵的行均值  
#         mean_squared = sparse_matrix.multiply(sparse_matrix).mean(axis=1).A1  # 每行平方均值  
#         variances = mean_squared - np.square(mean)  # 方差  

#         # 计算方差的 85% 分位数  
#         threshold = np.quantile(variances, 0.85)  
#         high_variance_indices = np.where(variances >= threshold)[0]  
#         sorted_indices = high_variance_indices[np.argsort(variances[high_variance_indices])[::-1]]
#         # 提取对应的行  
#         high_variance_rows = sparse_matrix[sorted_indices]  
        
#         # 提取对应的行名  
#         high_variance_row_names = [row_names[i] for i in sorted_indices]  
        
#         # 构建 DataFrame，包含行名和列名  
#         high_variance_df = pd.DataFrame(  
#             high_variance_rows.toarray(),  
#             index=high_variance_row_names,  
#             columns=col_names  
#         )   
        
#         # 保存为 CSV  
#         high_variance_df.T.to_csv('./mutation_matrix_features.csv')  

#     # 方法2: 处理 var2  
#     elif option == 'var2':  
#         # 使用稀疏矩阵优化计算  
#         mean = sparse_matrix.mean(axis=1).A1  
#         mean_squared = sparse_matrix.multiply(sparse_matrix).mean(axis=1).A1  
#         variances = mean_squared - np.square(mean)  

#         # 获取方差最大的 100 行索引  
#         top_variance_indices = np.argsort(variances)[-int(args.number):]
#         high_variance_rows = sparse_matrix[top_variance_indices]  
        
#         # 提取对应的行名  
#         high_variance_row_names = [row_names[i] for i in top_variance_indices]  
        
#         # 构建 DataFrame，包含行名和列名  
#         high_variance_df = pd.DataFrame(  
#             high_variance_rows.toarray(),  
#             index=high_variance_row_names,  
#             columns=col_names  
#         ) 
        
#         # 保存为 CSV  
#         high_variance_rows.T.to_csv('./mutation_matrix_features.csv')  

#     # 方法3: 处理 TF-IDF  
#     elif option == 'TF_IDF':  
#         # 计算 TF-IDF  
#         tfidf_transformer = TfidfTransformer(smooth_idf=True)  
#         tfidf_matrix = tfidf_transformer.fit_transform(sparse_matrix.T)  

#         # 计算每个列的方差  
#         tfidf_variances = tfidf_matrix.var(axis=0).A1  

#         # 获取方差最大的 100 个特征  
#         top_variance_features = np.argsort(tfidf_variances)[-int(args.number):]
        
#         binary_matrix_top_features = tfidf_matrix[:, top_variance_features]  
#         high_variance_row_names = [row_names[i] for i in top_variance_indices]
#         high_variance_df = pd.DataFrame(  
#             binary_matrix_top_features.toarray(),  
#             index=col_names,  
#             columns=high_variance_row_names  
#         )
#         # 保存为 CSV  
#         binary_matrix_top_features.to_csv('./mutation_matrix_features.csv')  

#     else:  
#         raise ValueError("Invalid option. Please choose 'var1', 'var2', or 'TF_IDF'.")  
    
def compute_variance_chunked(sparse_matrix, chunk_size=1000):  
    """分块计算行方差（适用于var1/var2方法）"""  
    n_rows = sparse_matrix.shape[0]  
    n_cols = sparse_matrix.shape[1]  
    
    # 初始化累加器  
    sum_x = np.zeros(n_rows)  
    sum_x2 = np.zeros(n_rows)  
    
    # 分块处理  
    for chunk_start in range(0, n_rows, chunk_size):  
        chunk_end = min(chunk_start + chunk_size, n_rows)  
        chunk = sparse_matrix[chunk_start:chunk_end, :]  
        
        # 转换为CSR以便逐行处理  
        if not issparse(chunk):  
            chunk = csr_matrix(chunk)  
            
        # 逐行计算  
        for i in range(chunk.shape[0]):  
            row_idx = chunk_start + i  
            row_data = chunk[i].data  
            sum_x[row_idx] = np.sum(row_data)  
            sum_x2[row_idx] = np.sum(row_data ** 2)  
    
    # 计算全局统计量  
    mean = sum_x / n_cols  
    mean_squared = sum_x2 / n_cols  
    variances = mean_squared - (mean ** 2)  
    return variances  

def process_data(option, sparse_matrix, row_names, col_names, treat_neg_one_as_zero, chunk_size=1000):  
    # 预处理：将-1替换为0（如果启用）  
    if treat_neg_one_as_zero == 1:  
        sparse_matrix = sparse_matrix.copy()  
        sparse_matrix.data[sparse_matrix.data == -1] = 0  

    # 方法选择分支  
    if option in ['var1', 'var2']:  
        # 分块计算方差  
        variances = compute_variance_chunked(sparse_matrix, chunk_size)  
        
        if option == 'var1':  
            # 取方差前85%分位数  
            threshold = np.quantile(variances, 0.85)  
            selected_indices = np.where(variances >= threshold)[0]  
            sorted_indices = selected_indices[np.argsort(variances[selected_indices])[::-1]]  
        elif option == 'var2':  
            # 取方差最大的N个特征  
            selected_indices = np.argsort(variances)[-int(args.number):]  
            sorted_indices = selected_indices[::-1]  # 降序排列  
            
        # 提取选中行并保存  
        high_var_matrix = sparse_matrix[sorted_indices]  
        high_var_df = pd.DataFrame(  
            high_var_matrix.toarray(),  
            index=[row_names[i] for i in sorted_indices],  
            columns=col_names  
        )  
        high_var_df.T.to_csv('./mutation_matrix_features.csv')  

    elif option == 'TF_IDF':  
        # 分块计算TF-IDF方差  
        def tfidf_variance_chunked(matrix, chunk_size=1000):  
            n_cols = matrix.shape[1]  
            sum_x = np.zeros(n_cols)  
            sum_x2 = np.zeros(n_cols)  
            
            # 分块处理列  
            for chunk_start in range(0, n_cols, chunk_size):  
                chunk_end = min(chunk_start + chunk_size, n_cols)  
                chunk = matrix[:, chunk_start:chunk_end].toarray()  
                sum_x[chunk_start:chunk_end] = chunk.sum(axis=0)  
                sum_x2[chunk_start:chunk_end] = (chunk ** 2).sum(axis=0)  
                
            mean = sum_x / matrix.shape[0]  
            return (sum_x2 / matrix.shape[0]) - (mean ** 2)  
        
        # 分块计算TF-IDF  
        tfidf = TfidfTransformer().fit_transform(sparse_matrix.T)  
        variances = tfidf_variance_chunked(tfidf)  
        
        # 选择top特征  
        top_indices = np.argsort(variances)[-int(args.number):]  
        selected_features = sparse_matrix.T[:, top_indices]  
        
        # 保存结果  
        pd.DataFrame(  
            selected_features.toarray(),   
            index=col_names,  
            columns=[row_names[i] for i in top_indices]  
        ).to_csv('./mutation_matrix_features.csv')  
        
    else:  
        raise ValueError("Invalid method option")    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculating Entropy and Simpsons Index')
    parser.add_argument('-i', '--input', required=True, help='h5mu_path')
    parser.add_argument('-o', '--output', required=True, help='out_path')
    parser.add_argument( '-m', '--method', choices=['var1', 'var2', 'TF_IDF'], default='var1',help='method for cho (options: var1, var2, TF_IDF, default: var1)') 
    parser.add_argument('-z', '--zetro',choices=['0', '1'], default='1',help='whether to treat -1 as 0 (0:No;1:Yes)')
    parser.add_argument('-n', '--number', default='100',help='the number of features for clustering')
    args = parser.parse_args()
    set_workspace(args.output) 
    Mudata = md.read(args.input)
    # df = pd.DataFrame.sparse.from_spmatrix(Mudata['mutation'].X, index=Mudata['mutation'].obs.index, columns=Mudata['mutation'].var.index).T
    A,B = filter_mutations_chunked(Mudata['mutation'].X.T,Mudata['mutation'].var.index.tolist(),chunk_size=1000)
    print(B)
    process_data(args.method,B,A,Mudata['mutation'].obs.index,int(args.zetro),chunk_size=1000)