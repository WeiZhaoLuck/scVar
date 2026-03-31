import requests
import os
import csv
import json
import pandas as pd
import re
import yaml
import subprocess
# from selenium import webdriver
# from bs4 import BeautifulSoup
from pandas.core.frame import DataFrame
import sys

path_all=sys.argv[1]
out_path=path_all+'/analysis/report/data'
# 文件夹中的CSV文件列表

target_char_signature = "_signature"
pattern_signature = re.compile(target_char_signature)
csv_folder = path_all+'/analysis/mutation_signature/cell_type'
csv_files = [file for file in os.listdir(csv_folder) if file.endswith('.csv')]
matching_files_signature = [file for file in csv_files if re.search(pattern_signature, file)]
# 创建一个用于存储结果的列表
data_list = []

#signature_all
signature_all = pd.read_csv(
    path_all+"/analysis/mutation_signature/Signature_celltype.csv", index_col=0)
target_char_cosmic = "_cosmic"
pattern_cosmic = re.compile(target_char_cosmic)
matching_files_cosmic = [file for file in csv_files if re.search(pattern_cosmic, file)]

# 逐个处理CSV文件
for csv_file in sorted(matching_files_cosmic, key=lambda x: (x.lower(), x)):
    file_path = os.path.join(csv_folder, csv_file)
    # 获取文件名（去除扩展名）
    file_name = os.path.splitext(csv_file)[0].split('_cosmic')[0]
    signature = file_name+'_signature.csv'
    file_path_signature = os.path.join(csv_folder, signature)
    # 读取CSV文件的数据
    data = []
    with open(file_path, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            data.append(row)
    bar_data=[]
    with open(file_path_signature, 'r') as signaturecsv:
        csv_reader = csv.reader(signaturecsv)
        for row_bar in csv_reader:
            bar_data.append(row_bar)
    # bar_data = signature_all.loc[file_name].values.tolist()
    # 将数据添加到结果列表中
    data_list.append({'id': file_name, 'name': file_name,
                     'bardata': bar_data, 'piedata': data})

# 将数据列表保存为JSON文件
json_text = 'let all=' + json.dumps(data_list, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)
# json_file_path = 'output.json'
# with open(json_file_path, 'w', encoding='utf-8') as json_file:
#     json.dump(data_list, json_file, ensure_ascii=False, indent=4)

# print(f'Data has been saved to {json_file_path}')


##table_data
# print(os.path.dirname(os.path.dirname(path_all))+"/config.yaml")
f1 = open(os.path.dirname(os.path.dirname(path_all))+"/config.yaml")

test = yaml.load(f1, Loader=yaml.FullLoader)
overview=[]
for key_index in range(len(list(test.keys()))):
    data_table_overview = {}
    data_table_overview['Parameters'] = list(test.keys())[key_index]
    data_table_overview['Value'] = list(test.values())[key_index]
    # if key_index == 2:
    #     data_table_overview['Value'] = list(
    #         list(test.values())[key_index].keys())[0]
    # else:
    #     data_table_overview['Value'] = list(test.values())[key_index]
    overview.append(data_table_overview)

# 将数据列表保存为JSON文件
json_text = 'let table1=' + \
    json.dumps(overview, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output_table_data.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)

# json_file_path_table = 'output_table_data.json'
# with open(json_file_path_table, 'w', encoding='utf-8') as json_file:
#     json.dump(overview, json_file, ensure_ascii=False, indent=4)

os.chdir(path_all)
result = subprocess.run(['wc', '-l', './seurat/barcodes_cell_type.tsv'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
cell_number = int(result.stdout.strip().split()[0])
result = subprocess.run(['wc', '-l', './seurat/result.csv'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
cell_type = int(result.stdout.strip().split()[0])

data_table_cell=[]
for i in range(2):
    data_table_cell_each={}
    if i==0:
       data_table_cell_each['Name'] = 'Number of cells after filtration'
       data_table_cell_each['Value'] = cell_number
    else:
        data_table_cell_each['Name'] = 'Number of cell types'
        data_table_cell_each['Value'] = cell_type
    data_table_cell.append(data_table_cell_each)

# 将数据列表保存为JSON文件
json_text = 'let cell=' + \
    json.dumps(data_table_cell, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output_table_data_cell.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)
# json_file_path_table_cell = 'output_table_data_cell.json'
# with open(json_file_path_table_cell, 'w', encoding='utf-8') as json_file:
#     json.dump(data_table_cell, json_file, ensure_ascii=False, indent=4)


##pie_cell_type_data
pie_cell_type_data = pd.read_csv(path_all+"/seurat/result.csv",header=None,names=['celltype','value'])
cell_type_total=[]
cell_type_count=[]
for i in range(len(pie_cell_type_data)):
    each_cell={}
    cell_count={}
    each_cell['name'] = pie_cell_type_data.loc[i,'celltype']
    cell_count['name'] = pie_cell_type_data.loc[i, 'celltype']
    each_cell['value'] = str(pie_cell_type_data.loc[i, 'value']/sum(pie_cell_type_data.loc[:, 'value'])*100)
    cell_count['value'] = str(pie_cell_type_data.loc[i, 'value'])
    cell_type_total.append(each_cell)
    cell_type_count.append(cell_count)
# json_file_path_pie = 'output_pie_cell_type.json'
# with open(json_file_path_pie, 'w', encoding='utf-8') as json_file:
#     json.dump(cell_type_total, json_file, ensure_ascii=False, indent=4)
# 将数据列表保存为JSON文件
cell_count = {}
cell_count['name'] ="all"
cell_count['value'] = cell_number
cell_type_count.append(cell_count)
json_text = 'let pie_data=' + \
    json.dumps(cell_type_total, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output_pie_cell_type.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)

json_text = 'let cell_data_count=' + \
    json.dumps(cell_type_count, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output_cell_data_count.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)


##bar_line_bam
bam_data = pd.read_csv(path_all+"/analysis/bam_statistics/chromosomes_results_sorted.csv", header=None, names=['Chromosome', 'Avgdepth', 'Coverage'])
bam_data_final = bam_data.loc[1:]
bam_data_total=[]
for i in range(len(bam_data_final.columns.tolist())):
    each={}
    if i==0:
        each['Chromosome'] = [h.strip()
                              for h in bam_data_final.loc[:, bam_data_final.columns.tolist()[i]]]
        
    else:
        each[bam_data_final.columns.tolist()[i].replace(' ', '')] = [
            h for h in bam_data_final.loc[:, bam_data_final.columns.tolist()[i]]]
    bam_data_total.append(each)
# json_file_path_bam_data = 'output_bam_data.json'
# with open(json_file_path_bam_data, 'w', encoding='utf-8') as json_file:
#     json.dump(bam_data_total, json_file, ensure_ascii=False, indent=4)
json_text = 'let bam_data=' + \
    json.dumps(bam_data_total, ensure_ascii=False, indent=4)
# 将JSON文本写入文件
with open(out_path+'/output_bam_data.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)


# ###TOP20
# top20_vcf = pd.read_csv(path_all+"/genotype/top20_final.tsv",sep='\t', header=0)
# col_name_vcf = ['CHROM_POS','CHROM', 'POS', 'REF', 'ALT',
#                 'Consequence_vep', 'IMPACT_vep', 'SYMBOL_vep', 'Gene_vep', 'HGNC_ID_vep', 'CCDS_vep', 'AF_vep', 'gnomADg_AF_vep', 'MAX_AF_vep', 'DamagePredCount', 'VEST4_score', 'REVEL_score', 'MPC_score', 'CADD_phred', 'DANN_score', 'CLNDISDB', '1000g2015aug_all', 'PMID', 'cancerhotsports_Hugo_Symbol', 'cancerhotsports_Variant_Classification', 'cancerhotsports_Variant_Type','cancerhotsports_Gene','cancerhotsports_IMPACT','cancerhotsports_TUMORTYPE','cancerhotsports_oncotree_organtype','CIVIC_CSQ','COS_CODING_ID','COS_CODING_Gene_name','COS_CODING_GENOMIC_ID','COS_CODING_HGVSC','cgi_mutation_type','cgi_cancer_types','brca_variant_nomenclature_clinical_significance','brca_variant_nomenclature_variant_haplotype','brca_ENIGMA_clinical_significance_citation','brca_ClinVar_SCV_accession','clinivar_ALLELEID','clinivar_CLNHGVS','clinivar_CLNVCSO']
#
# top20_num = pd.read_csv(path_all+"/genotype/top20.tsv", sep='\t', header=None,names=['Chrom_pos','chrom','pos','count'])
# order = top20_num.loc[:, 'Chrom_pos'].tolist()
# top20_vcf['CHROM_POS'] = pd.Categorical(
#     top20_vcf['CHROM_POS'], categories=order, ordered=True)
# sorted_df = top20_vcf.sort_values('CHROM_POS')
# sorted_df['mutation_cell_count'] = top20_num.loc[:, 'count'].tolist()
# col_name_vcf.append('mutation_cell_count')
# out = pd.DataFrame(sorted_df, columns=col_name_vcf).astype(str).reset_index(drop=True)
#
# # data={'CHROM_POS':order}
# # df_output=pd.DataFrame(data)
# # for i in col_name_vcf:
# #     # print(i)
# #     df_output[i] = top20_vcf.loc[:, col_name_vcf[i]].tolist()
# output_all_top20=[]
# for each in range(len(out)):
#     each_mutation = out.loc[each, :].to_dict()
#     output_all_top20.append(each_mutation)
# # json_file_path_top20 = 'output_table_data_top20.json'
# # with open(json_file_path_top20, 'w', encoding='utf-8') as json_file:
# #     json.dump(output_all_top20, json_file, ensure_ascii=False, indent=4)
# json_text = 'let top20=' + \
#     json.dumps(output_all_top20, ensure_ascii=False, indent=4)
# # 将JSON文本写入文件
# with open(out_path+'/output_table_data_top20.json', 'w', encoding='utf-8') as json_file:
#     json_file.write(json_text)


##top20_new
top20_vcf = pd.read_csv(path_all+"/genotype/top30_final.tsv",sep='\t', header=0)
col_name_vcf = ['CHROM_POS_REF_ALT','CHROM', 'POS', 'REF', 'ALT',
                'Consequence_vep', 'IMPACT_vep', 'SYMBOL_vep', 'Gene_vep', 'HGNC_ID_vep', 'CCDS_vep', 'AF_vep', 'gnomADg_AF_vep', 'MAX_AF_vep', 'DamagePredCount', 'VEST4_score', 'REVEL_score', 'MPC_score', 'CADD_phred', 'DANN_score', 'CLNDISDB', '1000g2015aug_all', 'PMID', 'cancerhotsports_Hugo_Symbol', 'cancerhotsports_Variant_Classification', 'cancerhotsports_Variant_Type','cancerhotsports_Gene','cancerhotsports_IMPACT','cancerhotsports_TUMORTYPE','cancerhotsports_oncotree_organtype','CIVIC_CSQ','COS_CODING_ID','COS_CODING_Gene_name','COS_CODING_GENOMIC_ID','COS_CODING_HGVSC','cgi_mutation_type','cgi_cancer_types','brca_variant_nomenclature_clinical_significance','brca_variant_nomenclature_variant_haplotype','brca_ENIGMA_clinical_significance_citation','brca_ClinVar_SCV_accession','clinivar_ALLELEID','clinivar_CLNHGVS','clinivar_CLNVCSO','mutation_cells_count','cell_types']
out = pd.DataFrame(top20_vcf, columns=col_name_vcf).astype(str).reset_index(drop=True)
output_all_top20=[]
for each in range(len(out)):
    each_mutation = out.loc[each, :].to_dict()
    output_all_top20.append(each_mutation)
json_text = 'let top20=' + \
    json.dumps(output_all_top20, ensure_ascii=False, indent=4)

with open(out_path+'/output_table_data_top20.json', 'w', encoding='utf-8') as json_file:
    json_file.write(json_text)

# with open("/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/snakemake/test/other/ERS15977678/genotype/output_table_data_top20.json", 'w', encoding='utf-8') as json_file:
#     json_file.write(json_text)

##output top mutation cell barcodes
out_top_barcode=pd.DataFrame(top20_vcf,columns=['CHROM_POS_REF_ALT','barcodes'])
out_top_barcode.to_csv(path_all+"/genotype/top30_cell.tsv",sep="\t",header=False)

###cosmic information
# all_cosmic = ['sbs'+str(i) for i in range(1, 7)] + \
#     ['sbs'+str(i) for i in range(8, 10)]+['sbs'+str(i) for i in range(11, 17)]+['sbs'+str(i) for i in range(18, 61)]+['sbs'+str(i)
#                                                                                                                       for i in range(84, 96)]+['sbs7a', 'sbs7b', 'sbs7c', 'sbs7d', 'sbs10a', 'sbs10b', 'sbs10c', 'sbs10d', 'sbs17a', 'sbs17b']
# data_all = {'cosmicid': all_cosmic}
# cosmic_data =DataFrame(data_all)
# pa_all=[]
# comm_all=[]
# for i in all_cosmic:
#     url_sub = "https://cancer.sanger.ac.uk/signatures/sbs/"
#     url = url_sub+i+'/'
#     response = requests.get(url)
#     if response.status_code == 200:
#         page_content = response.text
#         soup = BeautifulSoup(page_content, 'html.parser')
#         div_element = soup.find('section', id='proposed-aetiology')
#         if div_element:
#             # p_content = div_element.findAll('p')
#             pa = div_element.findAll('p')[0].get_text()
#             if len(div_element.findAll('p')) == 2:
#                 comm = div_element.findAll('p')[1].get_text()
#             else:
#                 comm = ""
#             # p_content = div_element.find('p').get_text()
#             pa_all.append(pa)
#             comm_all.append(comm)
#             # print(comm)
#         else:
#             print("Div element not found.")
#     else:
#         print("Failed to retrieve the page content.")
# cosmic_data['commens'] = comm_all
# cosmic_data['pa_data'] = pa_all
#
# output_all_cosmic_info = []
# for each in range(len(cosmic_data)):
#     each_mutation = cosmic_data.loc[each, :].to_dict()
#     output_all_cosmic_info.append(each_mutation)
# # json_file_path_top20 = 'output_table_data_top20.json'
# # with open(json_file_path_top20, 'w', encoding='utf-8') as json_file:
# #     json.dump(output_all_top20, json_file, ensure_ascii=False, indent=4)
# json_text = 'let cismic_info=' + \
#     json.dumps(output_all_cosmic_info, ensure_ascii=False, indent=4)
# # 将JSON文本写入文件
# with open(out_path+'/output_table_data_cosmic_info.json', 'w', encoding='utf-8') as json_file:
#     json_file.write(json_text)
