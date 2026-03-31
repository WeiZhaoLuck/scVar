import sys
import pysam
import pandas as pd
import os
from multiprocessing import Process, Queue, Pool
import multiprocessing as mp
from tqdm import tqdm  

### Extract bam file based on barcodes
def extract_bam(barcodes_path, chr):  
    bam_path = "%s/%s.bam" % (out_path, chr)  
    tmpfilename = "%s/%s.extract.bam" % (out_path, chr)  
    barcodes_all = pd.read_csv(barcodes_path, sep=",")  
    barcodes = barcodes_all.iloc[:, 0]  
    bam_all = pysam.AlignmentFile(bam_path)  
    with pysam.AlignmentFile(tmpfilename, "wb", header=bam_all.header) as outf:  
        for item in bam_all:  
            if item.has_tag('CB'):  
                a = item.get_tag('CB')  
                if a in set(barcodes):  
                    outf.write(item)  

def multicore_process_extract_bam(barcodes_path, dir):  
    pool = mp.Pool()   
    all_bam = os.listdir(dir)  
    all_chr = [o.replace(".bam", "") for o in all_bam]  
 
    multi_res = []  
    for n in range(len(all_chr)):   
        multi_res.append(pool.apply_async(extract_bam, (barcodes_path, all_chr[n],)))  
 
    result = []  
    with tqdm(total=len(all_chr), desc="Processing chromosomes") as pbar:  
        for res, chr_name in zip(multi_res, all_chr):  
            result.append(res.get()) 
            pbar.update(1)  
            pbar.set_postfix({"Completed": chr_name})  

    pool.close()  
    pool.join() 
    return result


if __name__=='__main__':
    barcodes_path = sys.argv[1]
    out_path = sys.argv[2]
    multicore_process_extract_bam(barcodes_path,out_path)