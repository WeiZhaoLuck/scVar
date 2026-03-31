'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-21 16:15:44
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-21 16:21:15
FilePath: \scvar_docker\codes\pre_pseudotime.py
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
import math 
from scipy import io
import os

def pre(Mudata,MutationList):
    df = pd.DataFrame(index=Mudata['mutation'].obs.index, columns=MutationList,data=0)
    positions = {}  
    for m in MutationList:  
        try:  
            position = Mudata['mutation'].var.index.tolist().index(m)  
            positions[m] = position  
        except ValueError:  
            positions[m] = None  # 如果突变不在参考列表中 
    result_indices = {}  

    # 遍历字典  
    for mutation, column_index in positions.items():  
        # 获取第 X 列（column_index）并查找值为 1 的行索引  
        column_data = Mudata['mutation'].X[:, column_index].todense()  # 使用 todense() 方法  
        # 将稠密数组转换为一维数组  
        column_data = np.array(column_data).flatten()  
        # column_data = Mudata['mutation'].X[:, column_index].A1  # 转换为一维数组  
        indices = np.where(column_data == 1)[0]  # 查找值为 1 的行索引  
    
        # 将找到的行索引保存在结果字典中  
        result_indices[mutation] = indices.tolist()  # 将 NumPy 数组转换为列表 
    for index, (mutation, indices) in enumerate(result_indices.items()):  
        df.iloc[indices, index] = 1  
    return df


def ouput(df):
    
    os.makedirs('matrix_files', exist_ok=True)
    adata = Mudata['rna'].raw.to_adata()
    with open('matrix_files/barcodes.tsv','w') as f:
        for item in adata.obs_names:
            f.write(item + '\n')
    with open('matrix_files/features.tsv','w') as f:
        for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
            f.write(item + '\n')
    io.mmwrite('matrix_files/matrix.mtx',adata.X.T)
    os.system('ls matrix_files/')  # 使用 os.system 运行 shell 命令  

    # 压缩文件  
    os.system('gzip matrix_files/*')  # 压缩文件  
    merged_df = pd.concat([Mudata['rna'].obs, df], axis=1)  
    merged_df.to_csv('metadata.csv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculating Entropy and Simpsons Index')
    parser.add_argument('-i', '--input', required=True, help='h5mu_path')
    parser.add_argument('-m', '--mutation', required=True, help='h5mu_path')
    parser.add_argument('-o', '--out', required=True, help='out_path')
    args = parser.parse_args()
    os.chdir(args.out)
    Mutations = args.mutation.split(',')
    Mudata = md.read(args.input)
    df = pre(Mudata,Mutations)
    ouput(df)