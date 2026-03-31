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
import os

import numpy as np  
import pandas as pd  
from scipy.sparse import csr_matrix, csc_matrix  
from sklearn.feature_extraction.text import TfidfTransformer  
import subprocess  

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

def process_data(option, sparse_matrix, treat_neg_one_as_zero):  
    """  
    处理稀疏矩阵并生成特征矩阵，支持 pandas.DataFrame 或 scipy.sparse 矩阵。  

    参数:  
    - option: 选择处理方法 ('var1', 'var2', 'TF_IDF')  
    - sparse_matrix: 输入的稀疏矩阵 (pandas.DataFrame 或 scipy.sparse 矩阵)  
    - treat_neg_one_as_zero: 是否将 -1 视为 0 (0: 不视为 0, 1: 视为 0)  
    """  

    # 确保输入为 scipy.sparse 矩阵  
    if isinstance(sparse_matrix, pd.DataFrame): 
        row_names = df.index
        col_names = df.columns
        sparse_matrix = csr_matrix(sparse_matrix.sparse.to_coo())  # 转为 Scipy 稀疏矩阵格式  
    elif not isinstance(sparse_matrix, (csr_matrix, csc_matrix)):  
        raise ValueError("Input sparse_matrix must be a pandas.DataFrame with sparse format or a scipy.sparse matrix.")  

    # 如果 treat_neg_one_as_zero 为 1，将 -1 替换为 0  
    if treat_neg_one_as_zero == 1:  
        sparse_matrix = sparse_matrix.copy()  
        sparse_matrix.data[sparse_matrix.data == -1] = 0  # 仅替换非零值中的 -1  

    # 方法1: 处理 var1  
    if option == 'var1':  
        # 使用稀疏矩阵计算每行方差  
        mean = sparse_matrix.mean(axis=1).A1  # 稀疏矩阵的行均值  
        mean_squared = sparse_matrix.multiply(sparse_matrix).mean(axis=1).A1  # 每行平方均值  
        variances = mean_squared - np.square(mean)  # 方差  

        # 计算方差的 85% 分位数  
        threshold = np.quantile(variances, 0.85)  
        high_variance_indices = np.where(variances >= threshold)[0]  
        sorted_indices = high_variance_indices[np.argsort(variances[high_variance_indices])[::-1]]
        # 提取对应的行  
        high_variance_rows = sparse_matrix[sorted_indices]  
        
        # 提取对应的行名  
        high_variance_row_names = [row_names[i] for i in sorted_indices]  
        
        # 构建 DataFrame，包含行名和列名  
        high_variance_df = pd.DataFrame(  
            high_variance_rows.toarray(),  
            index=high_variance_row_names,  
            columns=col_names  
        )   
        
        # 保存为 CSV  
        high_variance_df.T.to_csv('./mutation_matrix_features.csv')  

    # 方法2: 处理 var2  
    elif option == 'var2':  
        # 使用稀疏矩阵优化计算  
        mean = sparse_matrix.mean(axis=1).A1  
        mean_squared = sparse_matrix.multiply(sparse_matrix).mean(axis=1).A1  
        variances = mean_squared - np.square(mean)  

        # 获取方差最大的 100 行索引  
        top_variance_indices = np.argsort(variances)[-100:]
        high_variance_rows = sparse_matrix[top_variance_indices]  
        
        # 提取对应的行名  
        high_variance_row_names = [row_names[i] for i in top_variance_indices]  
        
        # 构建 DataFrame，包含行名和列名  
        high_variance_df = pd.DataFrame(  
            high_variance_rows.toarray(),  
            index=high_variance_row_names,  
            columns=col_names  
        ) 
        
        # 保存为 CSV  
        high_variance_rows.T.to_csv('./mutation_matrix_features.csv')  

    # 方法3: 处理 TF-IDF  
    elif option == 'TF_IDF':  
        # 计算 TF-IDF  
        tfidf_transformer = TfidfTransformer(smooth_idf=True)  
        tfidf_matrix = tfidf_transformer.fit_transform(sparse_matrix.T)  

        # 计算每个列的方差  
        tfidf_variances = tfidf_matrix.var(axis=0).A1  

        # 获取方差最大的 100 个特征  
        top_variance_features = np.argsort(tfidf_variances)[-100:]
        
        binary_matrix_top_features = tfidf_matrix[:, top_variance_features]  
        high_variance_row_names = [row_names[i] for i in top_variance_indices]
        high_variance_df = pd.DataFrame(  
            binary_matrix_top_features.toarray(),  
            index=col_names,  
            columns=high_variance_row_names  
        )
        # 保存为 CSV  
        binary_matrix_top_features.to_csv('./mutation_matrix_features.csv')  

    else:  
        raise ValueError("Invalid option. Please choose 'var1', 'var2', or 'TF_IDF'.")  
    
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
    df = pd.DataFrame.sparse.from_spmatrix(Mudata['mutation'].X, index=Mudata['mutation'].obs.index, columns=Mudata['mutation'].var.index).T
    process_data(args.method,df,int(args.zetro))
    result = subprocess.run(["Rscript", "/codes/MutationCluster.R", args.output,int(args.number),int(args.zetro)], capture_output=True, text=True)  
    # 打印输出  
    print("STDOUT:", result.stdout)  
    print("STDERR:", result.stderr) 