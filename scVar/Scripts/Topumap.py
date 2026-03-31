'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-28 16:34:43
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-31 16:03:32
FilePath: \scvar_docker\codes\Topumap.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import anndata as ad
import networkx as nx
import scanpy as sc
from matplotlib import rcParams
import pandas as pd
import scipy.sparse as sp
import numpy as np
import pandas as pd
import scanpy as sc
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import os
from scipy.stats import hypergeom  
from statsmodels.stats.multitest import multipletests 
import mudata as md
import re
import matplotlib.pyplot as plt
import sys

def plot_mutation_umap(
    mut_df,
    celltype_cols: list,
    rna,
):
    for idx, row in mut_df.iterrows():
        mut_id = row['mut_id']
        barcodes = row['barcodes'].split(',')  # 分割逗号连接的barcode
        
        # 筛选存在于RNA数据中的有效barcode
        valid_bcs = [bc.strip() for bc in barcodes if bc.strip() in rna.obs.index]
        if not valid_bcs:
            print(f"跳过突变 {mut_id}: 无有效barcode")
            continue
        # 创建画布：左子图显示所有细胞类型，右子图高亮突变细胞
        n_cols = len(celltype_cols) + 1
        fig, axes = plt.subplots(1, n_cols, figsize=(5*n_cols, 4), dpi=150)
        fig.suptitle(f"Mutation: {mut_id}", y=1.05)
        
        # 绘制每个细胞类型UMAP
        for i, col in enumerate(celltype_cols):
            sc.pl.umap(rna, color=col, ax=axes[i], show=False, title=col, frameon=False)
        
        # 高亮突变相关细胞：红色 vs 灰色
        highlight = rna.obs.index.isin(valid_bcs)
        rna.obs['highlight'] = ['mutated' if x else 'others' for x in highlight]
        sc.pl.umap(
            rna,
            color='highlight',
            palette={'mutated': 'red', 'others': 'lightgrey'},
            ax=axes[-1],
            show=False,
            title=f"Mutated Cells (n={len(valid_bcs)})",
            frameon=False
        )
        
        # 保存图片
        plt.tight_layout()
        os.makedirs('top_mutation_umap', exist_ok=True)
        plt.savefig(f"top_mutation_umap/{mut_id}.png", bbox_inches='tight')
        plt.close()
        
        print(f"已生成: {mut_id}.png")

if __name__ == "__main__":
    
    Md = md.read(sys.argv[1])
    rna = Md.mod['rna'] 
    mut_df = pd.read_csv(sys.argv[2], sep='\t', header=None, names=['chr_pos', 'mut_id', 'barcodes'])
    os.chdir(sys.argv[3])
    CELLTYPE_PATTERN = re.compile(r'cell_?type|leiden|cluster|labels', re.IGNORECASE)
    celltype_cols = [col for col in rna.obs.columns if CELLTYPE_PATTERN.search(col)]
    if not celltype_cols:
        print(f"警告: 无细胞类型列 (celltype/cell_type/leiden/cluster)")
    else:
        plot_mutation_umap(mut_df,celltype_cols,rna)
    n_cols = len(celltype_cols)
    # fig, axes = plt.subplots(1, n_cols, figsize=(5*n_cols, 4), dpi=150)
    # fig.suptitle(f"Mutation: {mut_id}", y=1.05)

    # 绘制每个细胞类型UMAP
    if len(celltype_cols) > 1:
        fig, axes = plt.subplots(1, len(celltype_cols), figsize=(5*len(celltype_cols), 4), dpi=150)
        for col, ax in zip(celltype_cols, axes):
            sc.pl.umap(rna, color=col, ax=ax, show=False, frameon=False)
    else:
        fig, ax = plt.subplots(figsize=(5, 4), dpi=150)
        sc.pl.umap(rna, color=celltype_cols, ax=ax, show=False, frameon=False)
    # fig, axes = plt.subplots(1, len(celltype_cols), figsize=(5*len(celltype_cols), 4), dpi=150)
    # for col, ax in zip(celltype_cols, axes):
    #     sc.pl.umap(rna, color=col, ax=ax, show=False, frameon=False)
    # sc.pl.umap(rna, color=celltype_cols, ax=ax, show=False, frameon=False)  # 一次性绘制所有细胞
    plt.savefig("cell_type.png", bbox_inches='tight')