import os
import pandas as pd
import pysam
import multiprocessing as mp
from tqdm import tqdm
import sys


def extract_bam(barcodes_path, chr, out_path):
    # 临时 BAM 文件路径
    """
    Extract the reads from a BAM file which have the barcodes in barcodes_path.

    Parameters
    ----------
    barcodes_path : str
        The path of the barcodes file.
    chr : str
        The name of the chromosome.
    out_path : str
        The path of the output directory.

    Returns
    -------
    None
    """
    tmpfilename = f"{out_path}/{chr}.extract.bam"

    # 读取条形码列表
    barcodes_all = pd.read_csv(barcodes_path, sep="\t")
    barcodes = set(barcodes_all.iloc[:, 0])  # 使用集合提高查找效率

    # 打开对应染色体的 BAM 文件
    bam_path = f"{out_path}/{chr}.bam"
    bam_all = pysam.AlignmentFile(bam_path, "rb")

    # 创建输出 BAM 文件
    with pysam.AlignmentFile(tmpfilename, "wb", header=bam_all.header) as outf:
        for item in bam_all:
            if item.has_tag('CB'):
                a = item.get_tag('CB')
                if a in barcodes:
                    outf.write(item)


def split_bam_by_chr(bam_file, chr, out_path):
    # 打开 BAM 文件并获取染色体信息
    """
    Split a BAM file into multiple files, each containing the reads from a particular chromosome.

    Parameters
    ----------
    bam_file : str
        The path of the BAM file.
    chr : str
        The name of the chromosome.
    out_path : str
        The path of the output directory.

    Returns
    -------
    None
    """
   
    bam_all = pysam.AlignmentFile(bam_file, "rb")

    # 创建临时 BAM 文件
    tmp_bam_path = f"{out_path}/{chr}.bam"
    with pysam.AlignmentFile(tmp_bam_path, "wb", header=bam_all.header) as tmp_bam:
        for read in bam_all.fetch(chr):
            tmp_bam.write(read)


def process_chromosome(barcodes_path, bam_file, chr, out_path):
    # 将 BAM 文件按照染色体分割
    """
    Process a single chromosome.

    Parameters
    ----------
    barcodes_path : str
        The path of the barcodes file.
    bam_file : str
        The path of the BAM file.
    chr : str
        The name of the chromosome.
    out_path : str
        The path of the output directory.

    Returns
    -------
    None
    """
    
    split_bam_by_chr(bam_file, chr, out_path)
    # 提取与条形码匹配的 reads
    extract_bam(barcodes_path, chr, out_path)
    tmp_bam_path = f"{out_path}/{chr}.bam"
    os.remove(tmp_bam_path)



def multicore_process_extract_bam(barcodes_path, bam_file, out_path):
    # 打开 BAM 文件并获取染色体信息
    """
    Process a BAM file by extracting reads that match the barcodes in parallel.

    Parameters
    ----------
    barcodes_path : str
        The path of the barcodes file.
    bam_file : str
        The path of the BAM file.
    out_path : str
        The path of the output directory.

    Returns
    -------
    None
    """
    
    bam_all = pysam.AlignmentFile(bam_file, "rb")

    # 获取所有染色体
    chromosomes = set()
    for read in bam_all:
        chromosomes.add(read.reference_name)
    chromosomes_cleaned = {chrom for chrom in chromosomes if chrom is not None}

        # 创建进程池
    pool = mp.Pool()
    multi_res = []

    # 为每个染色体调用 process_chromosome 函数
    for chr in chromosomes_cleaned:
        multi_res.append(pool.apply_async(process_chromosome, (barcodes_path, bam_file, chr, out_path)))

        # 显示进度条
    with tqdm(total=len(chromosomes_cleaned), desc="Processing chromosomes") as pbar:
        for res, chr_name in zip(multi_res, chromosomes_cleaned):
            res.get()  # 等待结果
            pbar.update(1)  # 更新进度条
            pbar.set_postfix({"Completed": chr_name})  # 设置进度条后缀信息

    pool.close()  # 关闭进程池
    pool.join()  # 等待所有进程完成


if __name__ == "__main__":
    # 输入参数

    barcodes_path = sys.argv[1]
    bam_file = sys.argv[2]
    out_path = sys.argv[3]

    # 确保输出目录存在
    os.makedirs(out_path, exist_ok=True)

    # 并行处理提取 BAM 文件
    multicore_process_extract_bam(barcodes_path, bam_file, out_path)