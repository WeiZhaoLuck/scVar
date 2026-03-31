'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-04 17:29:25
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-27 11:12:32
'''
import scanpy as sc
import os
import pandas as pd
import sys


'''
description: extract celltype from h5ad
param {*} h5ad_file_path
param {*} output_folder
return {*}
'''
def extract_celltype_from_h5ad(h5ad_file_path, output_folder):  
    # 读取h5ad文件  
    adata_combined = sc.read_h5ad(h5ad_file_path)  
    
    # 获取obs数据框  
    obs_data = adata_combined.obs  
    
    # 定义要匹配的关键字  
    keywords = ['celltype', 'cell_type']  
    
    # 使用模糊匹配查找包含关键字的列名  
    matched_columns = [col for col in obs_data.columns if any(keyword in col for keyword in keywords)]  
    
    # 选择第一个匹配的列名  
    selected_column = matched_columns[0] if matched_columns else None  
    
    # 判断输出文件夹是否存在，不存在则创建  
    if not os.path.exists(output_folder):  
        os.makedirs(output_folder)  
        print(f"File or directory exists: {output_folder}")  
    else:  
        print(f"File or directory does not exist: {output_folder}")  
    
    # 如果找到有效的列名，提取细胞和对应的细胞类型  
    if selected_column:  
        data = pd.DataFrame({  
            'cell': adata_combined.obs.index.tolist(),  # 使用obs的索引作为细胞名称  
            'celltype': obs_data[selected_column].str.replace(" ", "_")  # 替换空格为下划线  
        })
        output_data =  data[data['celltype'].notna() & (data['celltype'] != '')] 
        # output_data[output_data['celltype'] ]
        # 输出TSV文件的路径  
        output_file_path = os.path.join(output_folder, "output_cell_barcodes.tsv")  
        
        # 输出TSV文件  
        output_data.to_csv(output_file_path, sep='\t', index=False)  
        print(f"Output successful, file path is: {output_file_path}")  
    else:  
        print("No valid column names found.\n") 
    
    
    keywords = ['leiden', 'cluster']  
    
    # 使用模糊匹配查找包含关键字的列名  
    matched_columns = [col for col in obs_data.columns if any(keyword in col for keyword in keywords)]  
    
    # 选择第一个匹配的列名  
    selected_column = matched_columns[0] if matched_columns else None  
    
    # 判断输出文件夹是否存在，不存在则创建  
    if not os.path.exists(output_folder):  
        os.makedirs(output_folder)  
        print(f"File or directory exists: {output_folder}")  
    else:  
        print(f"File or directory does not exist: {output_folder}")  
    
    # 如果找到有效的列名，提取细胞和对应的细胞类型  
    if selected_column:  
        data = pd.DataFrame({  
            'cell': adata_combined.obs.index.tolist(),  # 使用obs的索引作为细胞名称  
            'cluster': obs_data[selected_column].str.replace(" ", "_")  # 替换空格为下划线  
        })
        output_data =  data[data['cluster'].notna() & (data['cluster'] != '')] 
        # output_data[output_data['celltype'] ]
        # 输出TSV文件的路径  
        output_file_path = os.path.join(output_folder, "output_cell_cluster.tsv")  
        
        # 输出TSV文件  
        output_data.to_csv(output_file_path, sep='\t', index=False)  
        print(f"Output successful, file path is: {output_file_path}")  
    else:  
        print("No valid column names found.\n") 


if __name__=='__main__':
    input_h5 = sys.argv[1]
    output_folder = sys.argv[2]
    extract_celltype_from_h5ad(input_h5, output_folder) 