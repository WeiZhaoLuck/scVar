import pandas as pd
import sys

genotype_path=sys.argv[1]
# cell_type_path=sys.argv[2]
out_path=sys.argv[2]

# cell_type=pd.read_table(cell_type_path,sep="\t",header=None,names=['barcode','cell_type'])
genotype_results=pd.read_table(genotype_path,sep="\t",header=None,names=['key','barcode','chrom','pos','ref','alt','ref_num','alt_num','all_num','cell_type'])

SNVs=set(genotype_results['key'].tolist())
SNVs_celltype=pd.DataFrame(data=SNVs,columns=['SNV'])
result = genotype_results.groupby('key')['cell_type'].apply(list).reset_index()
result['Unique_List_celltypes'] = result['cell_type'].apply(lambda x: list(set(x)))

result_final = genotype_results.groupby('key')['barcode'].apply(list).reset_index()
result['mutation_barcodes']=result_final['barcode']
result.to_csv(out_path,sep='\t', index=False)