import sys
import pysam
import pandas as pd
import os
from multiprocessing import Process, Queue, Pool
import multiprocessing as mp


###bamtools split -in tmp.bam -reference
###samtools merge [options] -o out.bam [options] in1.bam ... inN.bam
# barcodes_path=sys.argv[1]
# out_path=sys.argv[2]

# barcodes_path="/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/snakemake/test/E-MTAB-8410/ERS3858523/seurat/barcodes_type.csv"
# out_path="/p300s/baoym_group/zhaow/projects/scVar/PipelineV2/snakemake/test/E-MTAB-8410/ERS3858523/run_count_ERS3858523/outs/bam_right"

def extract_bam(chr,label,barcode_list):
    bam_path="%s/%s.bam"%(out_path,chr)
    tmpfilename="%s/%s_%s.extract.bam"%(out_path,chr,label)
    bam_all = pysam.AlignmentFile(bam_path)
    with pysam.AlignmentFile(tmpfilename, "wb", header=bam_all.header) as outf:
        for item in bam_all:
            if(item.has_tag('CB')):
                a=item.get_tag('CB')
                if(a in set(barcode_list)):
                    outf.write(item)

def multicore_process_extract_bam(barcode_list,all_chr,label):
    pool = mp.Pool(processes=8)
    multi_res = [pool.apply_async(extract_bam, (all_chr[n],label,barcode_list)) for n in range(0, len(all_chr))]
    result = [res.get() for res in multi_res]


if __name__=='__main__':
    barcodes_path = sys.argv[1]
    out_path = sys.argv[2]
    barcodes_all = pd.read_csv(barcodes_path, sep="\t",header=0)
    celltypes=list(set(barcodes_all.iloc[:,1]))
    all_bam = os.listdir(out_path)
    all_chr = [o.replace(".bam", "") for o in all_bam if o.startswith("extracted.REF_")]
    for i in celltypes:
        each_celltype=barcodes_all[barcodes_all['major_celltype']==i]
        each_celltype_barcodes=each_celltype.iloc[:,0]
        multicore_process_extract_bam(each_celltype_barcodes,all_chr,i.replace(' ', '_'))