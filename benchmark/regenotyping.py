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
import pandas as pd


def snp(alt_infor, bam_fil, baseq=0):  ##right
    c, s, e, ref, alt, gene, anno = alt_infor
    print(alt_infor)
    barcodes = defaultdict(list)  # CB
    # ======================================
    barUcodes = defaultdict(list)  # UB
    bar_count = {'ref': [], 'alt': []}  # UB
    with pysam.AlignmentFile(bam_fil, 'rb') as pile:
        # logger.info("processing variant: {}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        #     c,e,e, ref,alt, gene, anno))
        for pileupcolumn in pile.pileup(reference=c, start=(int(s) - 1), end=int(s),
                                        truncate=True, stepper="nofilter", max_depth=100000000, min_base_quality=baseq):
            for read in pileupcolumn.pileups:
                if read.is_del or read.is_refskip:
                    continue

                # print(read.query_position)
                ustring = read.alignment.query_name
                # print(ustring)
                # print(read.alignment.query_sequence[read.query_position])
                if read.alignment.query_sequence[read.query_position] != alt and read.alignment.is_duplicate == False:
                    bar_count['ref'].append(ustring)
                elif read.alignment.query_sequence[read.query_position] == alt and read.alignment.is_duplicate == False:
                    bar_count['alt'].append(ustring)
    return bar_count

if __name__ == '__main__':
    bam_file=sys.argv[1]
    variant_vcf=sys.argv[2]
    out_path=sys.argv[3]
    with open(out_path, 'w+') as varfile:
        for line in open(variant_vcf):
            result = line.rstrip('\n').split('\t')
            out = result[0:2]
            out.extend([result[1], result[3], result[4], 'test', 'anno'])
            ab = snp(out, bam_file, 0)
            ref_count = len(ab['ref'])
            alt_count = len(ab['alt'])
            out.extend([ref_count, alt_count])
            # print(out)
            outline = '{chrm}\t{st}\t{end}\t{ref}\t{alt}\t{gene}\t{anno}\t{ref_count}\t{alt_count}\t{tot}\n'.format(
                chrm=out[0],
                st=out[1],
                end=out[2],
                ref=out[3],
                alt=out[4],
                gene=out[5],
                anno=out[6],
                ref_count=out[7],
                alt_count=out[8],
                tot=(out[7] + out[8]))
            # print(outline)
            varfile.write(outline)


