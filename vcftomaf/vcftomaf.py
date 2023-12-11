#!/usr/bin/env python3
'''
convert annovar format (multianno.txt) to maftools required basic format
maf fields:
Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode
non maf fields:
VAF (Variant Allele Frequecy) and amino acid change information
usage: python convert_vcf_to_maf.py XXX.txt 1 barcode > XXX.maf
'''

import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.DEBUG)


def DictAnnoHgvs():
    convdict = {"exonic":"RNA",
                "splicing":"Splice_Site",
                "ncRNA":"RNA",
                "UTR5":"5'UTR",
                "UTR3":"3'UTR",
                "intronic":"Intron",
                "upstream":"5'Flank",
                "downstream":"3'Flank",
                "intergenic":"IGR",
                "frameshift insertion":"Frame_Shift_Ins",
                "frameshift deletion":"Frame_Shift_Del",
                "frameshift block substitution":"Frameshift_INDEL",
                "frameshift substitution":"Frameshift_INDEL",
                "stopgain":"Nonsense_Mutation",
                "stoploss":"Nonstop_Mutation",
                "nonframeshift insertion":"In_Frame_Ins",
                "nonframeshift deletion":"In_Frame_Del",
                "nonframeshift block substitution":"Inframe_INDEL",
                "nonframeshift substitution":"Inframe_INDEL",
                "nonsynonymous SNV":"Missense_Mutation",
                "synonymous SNV":"Silent",
                "unknown":"Unknown",
                "ncRNA_exonic":"RNA",
                "ncRNA_intronic":"RNA",
                "ncRNA_UTR3":"RNA",
                "ncRNA_UTR5":"RNA",
                "ncRNA":"RNA",
                "ncRNA_splicing":"RNA",
                "startloss":"Translation_Start_Site",
                "startgain":"Unknown"}
    return (convdict)


def VariantTypeHGVS(ref, alt):
    '''SNP, DNP, TNP, ONP, INS, DEL,'''
    l_ref = [i.upper() for i in ref.strip()]
    l_alt = [i.upper() for i in alt.strip()]
    for i in l_ref:
        if not i in ['A', 'T', 'C', 'G', '-']:
            raise Exception("invalid ref base : " + str(ref))
    for i in l_alt:
        if not i in ['A', 'T', 'C', 'G', '-']:
            raise Exception("invalid alt base : " + str(alt))
    if ref == '-' and alt != '-':
        return ('INS')
    if alt == '-' and ref != '-':
        return ('DEL')
    if ref != '-' and alt != '-':
        len_ref = len(ref)
        len_alt = len(alt)
        if len_alt == 1 and len_ref == 1:
            return ('SNP')
        if len_alt > len_ref:
            return ('INS')
        if len_alt < len_ref:
            return ('DEL')


def ConvVariantClassAnnovar2HGVS(a, b):
    d = DictAnnoHgvs()
    c = None
    if a != ".":
        for k, v in d.items():
            if k == a.strip():
                c = v
    else:
        for k, v in d.items():
            if k == b.strip():
                c = v
    if c:
        return (c)
    else:
        print(a + "\t" + b)
        raise Exception("conversion failed")


def ConvVariantClassHGVS2Annovar(a):
    d = DictAnnoHgvs()
    c = []
    for k, v in d.items():
        if v == a.strip():
            c.append(k)
    return (c)


def vep_col(info_vep):
    vep_result = info_vep.split("|")
    return (vep_result)


def convert_annovar_2_maftools(annovar=None, skip=0, hugo=None, chro=None, start=None, end=None, ref=None, alt=None,
                               func=None, vc=None, barcode=None, vaf=None, aachage=None):
    hugo = int(hugo) - 1
    chro = int(chro) - 1
    start = int(start) - 1
    end = int(end) - 1
    ref = int(ref) - 1
    alt = int(alt) - 1
    func = int(func) - 1
    vc = int(vc) - 1
    # barcode = int(barcode) - 1
    barcode = str(barcode)
    vaf = int(vaf) - 1
    aachage = int(aachage) - 1

    csq_normal = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
    all_csq_col = [u + "_vep" for u in csq_normal.split("|")]

    data_all = pd.read_csv(annovar, sep="\t")
    colnames = data_all.columns
    result = []
    name_list = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                 "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "VAF", "aachange"]
    other_col = [colnames[y] for y in range(len(colnames)) if
                 y not in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 102, 103, 104,
                           105, 106]]
    name_list_all = name_list + other_col + all_csq_col
    # result.append("\t".join(["Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","VAF","aachange"]
    # name_list.insert(colnames[y] for y in range(len(colnames)) if y not in [0,1,2,3,4,5,6,7,8,9,92,93,94,95,96,97,98,99,100,101,102,102,103,104,105,106])
    # name_list.append([colnames[y] for y in range(len(colnames)) if y not in [0,1,2,3,4,5,6,7,8,9,92,93,94,95,96,97,98,99,100,101,102,102,103,104,105,106]])
    # result.append("\t".join(["Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","VAF","aachange"]))
    result.append("\t".join(item for item in name_list_all))
    # with open(annovar,'r') as fh :
    with open(annovar, 'r') as fh:
        if skip > 0:
            for a in range(int(skip), 0, -1):
                next(fh)
        for r in fh.readlines():
            l_r = r.strip().split("\t")
            ahugo = l_r[hugo]
            achro = l_r[chro]
            astart = l_r[start]
            aend = l_r[end]
            aref = l_r[ref]
            aalt = l_r[alt]
            avt = VariantTypeHGVS(aref, aalt)
            abarcode = barcode
            avaf = l_r[vaf]
            aaachage = l_r[aachage]
            other_info = [l_r[y] for y in range(len(l_r)) if
                          y not in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 102,
                                    103, 104, 105, 106]]
            # print(str(l_r[vc]) +" "+str(avc))
            # print (r)
            count = len(l_r[func].strip().split(";"))
            # each_row=[]

            if count > 1:
                l_hugo = [i.strip() for i in ahugo.strip().split(";")]
                l_func = [i.strip() for i in l_r[func].strip().split(";")]
                for i in range(count):
                    # print(astart+" "+l_r[vc]+" "+l_func[i])
                    avc = ConvVariantClassAnnovar2HGVS(l_r[vc], l_func[i])
                    each_row_1 = [str(l_hugo[i]), str(achro), str(astart), str(aend), str(aref), str(aalt), str(avc),
                                  str(avt), str(abarcode), str(avaf), str(aaachage)]
                    # result.append("\t".join([str(l_hugo[i]),str(achro),str(astart),str(aend),str(aref),str(aalt),str(avc),str(avt),str(abarcode),str(avaf),str(aaachage)]))
            else:
                avc = ConvVariantClassAnnovar2HGVS(l_r[vc], l_r[func])
                each_row_1 = [str(ahugo), str(achro), str(astart), str(aend), str(aref), str(aalt), str(avc), str(avt),
                              str(abarcode), str(avaf), str(aaachage)]
            vep_info = vep_col(l_r[102])
            # result.append("\t".join([str(ahugo),str(achro),str(astart),str(aend),str(aref),str(aalt),str(avc),str(avt),str(abarcode),str(avaf),str(aaachage)]))
            each_row = each_row_1 + other_info + vep_info
            result.append("\t".join(each_row))
            # result.append("\t".join([l_r[u] for u in range(len(l_r)) if u not in [0,1,2,3,4,5,6,7,8,9,92,93,94,95,96,97,98,99,100,101,102,102,103,104,105,106]]))
    return (result)


if __name__ == "__main__":
    if len(sys.argv[1:]) < 2:
        print("python convert_annovar_2_maftools.py annvarfile(table) skiplines")
        sys.exit(1)
    # int(sys.argv[2])
    a = convert_annovar_2_maftools(sys.argv[1], int(sys.argv[2]),
                                   hugo=7,
                                   chro=1,
                                   start=2,
                                   end=3,
                                   ref=4,
                                   alt=5,
                                   func=6,
                                   vc=9,  # variant class
                                   barcode=sys.argv[3],
                                   vaf=106,
                                   aachage=10)
    for i in a:
        print(i)







