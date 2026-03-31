from collections import defaultdict
import pysam
import logging as logger
import sys
import argparse
import os
import multiprocessing as mp
from collections import Counter
import time
from datetime import datetime
import cProfile, pstats, io
import subprocess
import tempfile
from tqdm import tqdm


### Call indels
def indel(alt_infor, bam_fil):
    """
    :param alt_infor: list of variant information (chr, start, end, ref, alt, gene, annotation)
    :param bam_fil: bam file
    :return: barcodes and barcode counts (bar_count includes ref and alt reads labels [CB:UB])
    """
    c, s, e, ref, alt, gene, anno = alt_infor
    barcodes = defaultdict(list)  # CB
    barUcodes = defaultdict(list)  # UB
    bar_count = {'ref': [], 'alt': []}  # UB
    with pysam.AlignmentFile(bam_fil, 'rb') as pile:
        logger.info("processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            c, e, e, ref, alt, gene, anno))

        for pileupcolumn in pile.pileup(reference=c, start=(int(s) - 1), end=int(s), truncate=True, stepper="nofilter",
                                        max_depth=100000000):
            for read in pileupcolumn.pileups:

                if read.is_refskip or not read.alignment.has_tag('UB'):
                    continue
                q_name = read.alignment.query_name
                q_tag = read.alignment.get_tag('CB')
                qU_tag = read.alignment.get_tag('UB')
                qstring = q_tag + ':' + qU_tag
                barcodes['{}'.format(q_name)].append(q_tag)
                barUcodes['{}'.format(q_name)].append(qstring)

                if alt == '-':
                    if read.query_position == None:
                        alt_tag = read.alignment.get_tag('CB')
                        altU_tag = read.alignment.get_tag('UB')
                        ustringa = alt_tag + ':' + altU_tag
                        bar_count['alt'].append(ustringa)
                    else:
                        ref_tag = read.alignment.get_tag('CB')
                        refU_tag = read.alignment.get_tag('UB')
                        ustring = ref_tag + ':' + refU_tag
                        bar_count['ref'].append(ustring)
                else:
                    if read.indel == indel and \
                            read.alignment.query_alignment_sequence[read.query_position +
                                                                    1:read.query_position + read.indel + 1] == alt:
                        alt_tag = read.alignment.get_tag('CB')
                        # UB
                        altU_tag = read.alignment.get_tag('UB')
                        ustringa = alt_tag + ':' + altU_tag
                        bar_count['alt'].append(ustringa)
                    else:
                        ref_tag = read.alignment.get_tag('CB')
                        # UB
                        refU_tag = read.alignment.get_tag('UB')
                        ustring = ref_tag + ':' + refU_tag
                        bar_count['ref'].append(ustring)
        return barUcodes, bar_count


### Call snps
def snp(alt_infor, bam_fil, baseq, mapq):  ##right
    """
    :param alt_infor: list of variant information (chr, start, end, ref, alt, gene, annotation)
    :param bam_fil: bam file
    :param baseq: base quality
    :param mapq: map quality
    :return: barcodes and barcode counts (bar_count includes ref and alt reads labels [CB:UB])
    """
    c, s, e, ref, alt, gene, anno = alt_infor
    barcodes = defaultdict(list)  # CB
    barUcodes = defaultdict(list)  # UB
    bar_count = {'ref': [], 'alt': [], 'total': []}  # UB
    with pysam.AlignmentFile(bam_fil, 'rb') as pile:
        logger.info("processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            c, e, e, ref, alt, gene, anno))
        for pileupcolumn in pile.pileup(reference=c, start=(int(s) - 1), end=int(s),
                                        truncate=True, stepper="nofilter", max_depth=100000000, min_base_quality=baseq,
                                        min_mapping_quality=mapq):
            for read in pileupcolumn.pileups:
                #                 print(read.alignment.query_name)
                if read.is_del or read.is_refskip or not read.alignment.has_tag('UB') or \
                        int(read.alignment.query_qualities[read.query_position]) < 0 or \
                        int(read.alignment.mapping_quality) < 0 or not read.alignment.has_tag('CB') or \
                        read.alignment.query_sequence[read.query_position] == 'N':
                    continue
                q_name = read.alignment.query_name
                q_tag = read.alignment.get_tag('CB')
                qU_tag = read.alignment.get_tag('UB')
                qstring = q_tag + ':' + qU_tag
                barcodes['{}'.format(q_name)].append(q_tag)
                barUcodes['{}'.format(q_name)].append(qstring)
                ustring = read.alignment.get_tag('CB') + ':' + read.alignment.get_tag('UB')
                bar_count['total'].append(ustring)
                if read.alignment.query_sequence[read.query_position] == ref:
                    bar_count['ref'].append(ustring)
                elif read.alignment.query_sequence[read.query_position] == alt:
                    bar_count['alt'].append(ustring)
    return barUcodes, bar_count


### Summarize each variant
def sum_each(each_var, var_info, wu):
    """
    :param each_var: barcode counts (bar_count includes ref and alt reads labels [CB:UB], results from snp or indel)
    :param var_info: variant information (chr, start, end, ref, alt, gene, annotation)
    :param wu: output file
    :return: None
    """
    ref = each_var['ref']
    alt = each_var['alt']
    ref_barcode_ub = Counter([i for i in ref])
    # print(ref_barcode_ub)
    alt_barcode_ub = Counter([i for i in alt])
    all_barcode_ub = list(set(list(ref_barcode_ub.keys()) + list(alt_barcode_ub.keys())))
    c, s, e, ref, alt, gene, anno = var_info
    logger.info("Consensus calc variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(c, s, e, ref, alt, gene, anno))
    for ub in all_barcode_ub:
        cb = ub.split(':')[0]
        ref_counts = ref_barcode_ub.get(ub, 0)
        alt_counts = alt_barcode_ub.get(ub, 0)
        cbtager = '{chrm}\t{st}\t{en}\t{ref}\t{alt}\t' \
                  '{type}\t{gene}\t{bar}\t{bar_cb}\t{ref1}\t{alt1}\t{tot}\n'.format(
            chrm=c,
            st=s,
            en=e,
            ref=ref,
            alt=alt,
            type=anno,
            gene=gene,
            bar=ub,
            bar_cb=cb,
            alt1=alt_counts,
            ref1=ref_counts,
            tot=(ref_counts + alt_counts))
        wu.write(cbtager)


### Summarize each cell
def sum_cell(each_var, var_info, wu, rf):
    """
    :param each_var: barcode counts (bar_count includes ref and alt reads labels [CB:UB], results from snp or indel)
    :param var_info: variant information (chr, start, end, ref, alt, gene, annotation)
    :param wu: output file
    :param rf: minimum number of reads covering the locus in the cell
    :return: None
    """
    ref_cells = Counter([item.split(":")[0] for item in each_var['ref']])
    alt_cells = Counter([item.split(":")[0] for item in each_var['alt']])
    total_cb = list(set([item.split(":")[0] for item in each_var['total']]))
    c, s, e, ref, alt, gene, anno = var_info
    for cb in total_cb:
        ref_counts = ref_cells.get(cb, 0)
        alt_counts = alt_cells.get(cb, 0)
        if (ref_counts + alt_counts >= rf):
            cbtager = '{chrm}\t{st}\t{ref}\t{alt}\t' \
                      '{bar}\t{ref1}\t{alt1}\t{tot}\n'.format(
                chrm=c,
                st=s,
                en=e,
                ref=ref,
                alt=alt,
                type=anno,
                gene=gene,
                bar=cb,
                alt1=alt_counts,
                ref1=ref_counts,
                tot=(ref_counts + alt_counts))
            wu.write(cbtager)


def each_work(bam_path, tbl_path, outfile, indel_count, bq, mq, rf):
    command = "head -n 1 " + tbl_path + " | awk -F " + "\'\\t\'" + " '{print $2}'"
    start_each = subprocess.getoutput(command)
    command = "tail -n 1 " + tbl_path + " | awk -F " + "\'\\t\'" + " '{print $2}'"
    end_each = subprocess.getoutput(command)
    command = "head -n 1 " + tbl_path + " | awk -F " + "\'\\t\'" + " '{print $1}'"
    chr = subprocess.getoutput(command)
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    fp = tempfile.NamedTemporaryFile(mode='w+b', suffix=".bam", delete=False)
    pairedreads = pysam.AlignmentFile(fp.name, "wb", header=bamfile.header)
    for read in bamfile.fetch(contig=str(chr), start=int(start_each) - 10, stop=int(end_each) + 10):
        pairedreads.write(read)
    pairedreads.close()
    tempfile_path = fp.name
    pysam.index(tempfile_path)

    num_lines_variants = sum(1 for line in open(tbl_path))  # minus header
    logger.info("Number of variants:\t{}".format(num_lines_variants))
    name = os.path.basename(tbl_path)
    name_out = name.split(".tbl")[0]
    varfp = outfile + '+' + name_out + '_All.tsv'
    CBfp = outfile + '+' + name_out + '_counts_CB.tsv'
    with open(varfp, 'w+') as varfile, \
            open(CBfp, 'w+') as CB, \
            open(tbl_path, 'r') as regions:
        if (indel_count == 0):
            for lines in regions:
                var = lines.rstrip('\n').split('\t')
                a, b = snp(var, tempfile_path, bq, mq)
                if (len(b['total']) == len(b['ref']) + len(b['alt'])):
                    sum_each(b, var, varfile)
                    sum_cell(b, var, CB, rf)
                # else:
                #     print(len(b['total']))
                #     print(len(b['ref']))
                #     print(len(b['alt']))
        else:
            for lines in regions:
                var = lines.rstrip('\n').split('\t')
                a, b = indel(var, bam_path)
                sum_each(b, var, CB, varfile)
    os.remove(tempfile_path)
    os.remove(tempfile_path + ".bai")


def multi_process2(bampath, tbl_path, outfile, bq, indel_count, threads):
    pool = mp.Pool(processes=threads)
    result = []
    chr_vcf = [f.split('.tbl')[0] for f in os.listdir(tbl_path) if f.endswith('.tbl')]
    # print(chr_vcf)
    for i in range(len(chr_vcf)):
        chr = chr_vcf[i].split('.')[0]
        if (chr[0].isdigit()):
            bam_path_each = bampath + "/" + chr + ".bam"
        else:
            bam_path_each = bampath + "/" + chr_vcf[i] + ".bam"
        tbl_path_each = tbl_path + "/" + chr_vcf[i] + ".tbl"
        result.append(pool.apply_async(each_work, args=(bam_path_each, tbl_path_each, outfile, indel_count, bq,)))
    pool.close()
    pool.join()
    for i in result:
        i.get()


### Multi process
def multi_process(bampath, tbl_path, outfile, bq, mq, rf, indel_count, threads):
    """_summary_

    Args:
        bampath (string): bam file path
        tbl_path (string): variant file path
        outfile (string): output file path
        bq (int): base quality
        mq (int): map quality
        rf (int): minimum number of reads covering the locus in the cell
        indel_count (int): indel is or no(is:1;no:0)
        threads (int): number of threads
    """
    pool = mp.Pool(processes=threads)
    result = []
    chr_vcf = [f.split('.tbl')[0] for f in os.listdir(tbl_path) if f.endswith('.tbl')]
    # print(chr_vcf)
    for chr_name in tqdm(chr_vcf, desc="Processing chromosomes"):
        tbl_path_each = os.path.join(tbl_path, f"{chr_name}.tbl")
        result.append(pool.apply_async(each_work, args=(bampath, tbl_path_each, outfile, indel_count, bq, mq, rf)))

    for i in tqdm(result, desc="Waiting for results"):
        i.get()

    pool.close()
    pool.join()

    print("All done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse CB barcodes from Single cell rna seq data')
    parser.add_argument('bam_file', help='BAM file')
    parser.add_argument('variant_file', help='variants file with header')
    # parser.add_argument('barcodes', help='list of good barcodes file')
    parser.add_argument('upn', help='upn/sample name: will be used as prefix for out_file')
    parser.add_argument('-f', "--filter", type=int, default=0,
                        help='number of reads required per barcode default: 0')
    parser.add_argument('-mq', "--mapq", type=int, default=0,
                        help='Skip read with mapq smaller than default : 0')
    parser.add_argument('-bq', "--baseq", type=int, default=1,
                        help='Skip bases with base quality less than default : 1')
    parser.add_argument('-p', "--parallel", type=int, default=1,
                        help='Parallel, default : 1')
    args = parser.parse_args()

    bam_file = args.bam_file
    variants = args.variant_file
    # barcodes_good = args.barcodes
    outfile = args.upn
    bq = args.baseq
    mq = args.mapq
    rf = args.filter
    parallel = args.parallel
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    indel_count = 0
    pr = cProfile.Profile()
    pr.enable()
    logger.basicConfig(filename=outfile + '.log', filemode='w+',
                       level=logger.DEBUG,
                       format='%(asctime)s %(levelname)s %(message)s')
    logger.info("Start process")
    multi_process(bam_file, variants, outfile, bq, mq, rf, indel_count, parallel)
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logger.info("end process")

