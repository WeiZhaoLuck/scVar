'''
Author: WeiZhaoLuck zhaow@big.ac.cn
Date: 2025-03-19 16:08:54
LastEditors: WeiZhaoLuck zhaow@big.ac.cn
LastEditTime: 2025-03-20 09:47:03
FilePath: \scvar_docker\codes\extract_bam_bycelltypes.py
Description: 
'''
import pysam
from collections import defaultdict
import sys

work_path = sys.argv[1]
sample = sys.argv[2]
barcodes_path = work_path+"/seurat/output_cell_barcodes.tsv"
out_path = "{}/mapping/run_count_{}/outs/bam_tmp".format(work_path, sample)  
input_bam = "{}/mapping/run_count_{}/outs/extracted.bam".format(work_path, sample) 
# 构建barcode到细胞类型的映射
barcode_map = {}
ALL_CELL_TYPES = [] 
with open(barcodes_path) as f:
    next(f)  # 跳过标题行
    for line in f:
        barcode, cell_type = line.strip().split("\t")
        barcode_map[barcode] = cell_type
        if cell_type not in ALL_CELL_TYPES:
            ALL_CELL_TYPES.append(cell_type)


infile=pysam.Samfile(input_bam, "rb")
DICT_files = {}
for cell_type in ALL_CELL_TYPES:
    outfile="{}/{}.bam".format(out_path, cell_type)
    outfile_wb = pysam.AlignmentFile(outfile, "wb",template=infile)
    DICT_files[cell_type] = outfile_wb

for read in infile.fetch():
    # if read.has_tag("CB"):
    #     cb = read.get_tag("CB")
    #     if cb in barcode_map:
    #         cell_type = barcode_map[cb]
    #         DICT_files[cell_type].write(read)
    try:
        cb = read.opt("CB")
        print(cb)
        if cb in barcode_map:
            cell_type = barcode_map[cb]
            DICT_files[cell_type].write(read)
    except:
        FILTER2 = 'CB_not_found'
    
for i in DICT_files.values():
    i.close()


# 初始化输出文件
# in_bam = pysam.AlignmentFile(input_bam, "rb")
# writers = defaultdict(lambda: None)

# for read in in_bam:
#     if not read.has_tag("CB"):
#         continue
#     cb = read.get_tag("CB")
#     if cb not in barcode_map:
#         continue
#     cell_type = barcode_map[cb]
#     if not writers[cell_type]:
#         writers[cell_type] = pysam.AlignmentFile(
#             f"{out_path}/{cell_type}.bam", "wb", template=in_bam
#         )
#     writers[cell_type].write(read)

# # 关闭所有文件
# in_bam.close()
# for writer in writers.values():
#     if writer:
#         writer.close()
