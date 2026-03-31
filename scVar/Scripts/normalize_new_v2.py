import logging
import os.path
import sys
import pandas as pd
import re
import codecs

def vcfanno(string):
    """
    This function takes a vcfanno annotation string as input and processes it to extract specific information. The input string is expected to be in a specific format, with each piece of information separated by a semicolon (;). The function splits the input string by semicolons and iterates over the resulting list. For each item in the list, it checks if it contains an equal sign (=). If it does, it splits the item into a key-value pair. The function then checks if the key is present in a predefined list of fields. If the key is present, it updates the corresponding value in the `all_result` DataFrame. Finally, the function returns the first row of the `all_result` DataFrame as a list.

    Parameters:
    - string: A string containing the information to be processed.

    Returns:
    - A list of values representing the extracted information.

    Example Usage:
    ```
    string = "field1=value1;field2=value2;field3=value3"
    result = vcfanno(string)
    print(result)
    ```
    Output:
    ```
    ['value1', 'value2', 'value3']
    ```
    """
    all_fields=["DOCM_DISEASE","PMID","cancerhotsports_Hugo_Symbol","cancerhotsports_Variant_Classification","cancerhotsports_Variant_Type","cancerhotsports_dbSNP_RS","cancerhotsports_HGVSc","cancerhotsports_HGVSp","cancerhotsports_HGVSp_Short","cancerhotsports_Exon_Number","cancerhotsports_Gene","cancerhotsports_Feature","cancerhotsports_Feature_type","cancerhotsports_Consequence","cancerhotsports_IMPACT","cancerhotsports_VARIANT_CLASS","cancerhotsports_TUMORTYPE","cancerhotsports_Amino_Acid_Change","cancerhotsports_oncotree_organtype","cancerhotsports_oncotree_parent","cancerhotsports_oncotree_detailed","cancerhotspot_SNV_Tumer_TotalCounts","cancerhotspot_SNV_Tumer_Counts","CIVIC_CSQ","COS_CODING_ID","COS_CODING_Gene_name","COS_CODING_STRAND","COS_CODING_GENOMIC_ID","COS_CODING_LEGACY_ID","COS_CODING_CDS","COS_CODING_AA","COS_CODING_HGVSC","COS_CODING_HGVSP","COS_CODING_HGVSG","COS_CODING_CNT","COS_CODING_ID","COS_CODING_Gene_name","COS_CODING_STRAND","COS_CODING_LEGACY_ID","COS_CODING_CNT","in_cgi", "cgi_mutation_type", "cgi_cancer_types", "cgi_abbreviations", "cgi_reported_by","cgi_gene", "cgi_gene_tumorigenesis", "cgi_gene_alterations", "cgi_gene_translocations", "cgi_gene_cancer_types", "cgi_gene_abbreviations", "cgi_gene_sources","brca_geneSymbol","brca_variant_nomenclature_reference_cDNA_sequence","brca_variant_nomenclature_HGVS_nucleotide","brca_variant_nomenclature_BIC_designation","brca_variant_nomenclature_genome","brca_variant_nomenclature_RNA_LOVD","brca_variant_nomenclature_synonyms","brca_variant_nomenclature_pathogenicity_all","brca_variant_nomenclature_clinical_significance","brca_variant_nomenclature_variant_haplotype","brca_ENIGMA_clinical_significance_citation","brca_ENIGMA_allel_origin","brca_ClinVar_SCV_accession","brca_ClinVar_DBID_LOVD","brca_ClinVar_BIC_Nomenclature","clinivar_ALLELEID","clinivar_CLNDN","clinivar_CLNDNINCL","clinivar_CLNDISDBINCL","clinivar_CLNHGVS","clinivar_CLNREVSTAT","clinivar_CLNSIG","clinivar_CLNSIGCONF","clinivar_CLNSIGINCL","clinivar_CLNVC","clinivar_CLNVCSO","clinivar_CLNVI","clinivar_DBVARID","clinivar_GENEINFO","clinivar_MC","clinivar_ORIGIN","clinivar_RS","clinivar_SSR"]
    all_result=pd.DataFrame(data=[["."]*len(all_fields)],columns=all_fields)
    info_list=string.split(";")
    for item in info_list:
        if("=" in item):
            [k,y]=item.split('=',1)
            if k in all_fields:
                all_result.iloc[0][k]=y
        else:
            continue
    return all_result.iloc[0,:].tolist()
def annovar_split(string):
    """
    Splits a annovar annotation string into multiple substrings based on a delimiter and returns a list of the substrings.

    Parameters:
    - string (str): The string to be split.

    Returns:
    - all_out (list): A list of the substrings.

    Example:
    >>> annovar_split("abc=x;def=y;ghi=z")
    ['x', 'y', 'z']
    """
    all_out=[]
    info_split=string.split(";")
    for item in range(len(info_split)-1):
        [k,y]=info_split[item].split('=',1)
        all_out.append(y)
    return all_out
def header():
    """
    Generate the header of the annotation file. 

    :return: header.
    """
    position_info=['CHROM','POS','REF','ALT']
    csq_normal = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
    vcfanno_fields=["DOCM_DISEASE","PMID","cancerhotsports_Hugo_Symbol","cancerhotsports_Variant_Classification","cancerhotsports_Variant_Type","cancerhotsports_dbSNP_RS","cancerhotsports_HGVSc","cancerhotsports_HGVSp","cancerhotsports_HGVSp_Short","cancerhotsports_Exon_Number","cancerhotsports_Gene","cancerhotsports_Feature","cancerhotsports_Feature_type","cancerhotsports_Consequence","cancerhotsports_IMPACT","cancerhotsports_VARIANT_CLASS","cancerhotsports_TUMORTYPE","cancerhotsports_Amino_Acid_Change","cancerhotsports_oncotree_organtype","cancerhotsports_oncotree_parent","cancerhotsports_oncotree_detailed","cancerhotspot_SNV_Tumer_TotalCounts","cancerhotspot_SNV_Tumer_Counts","CIVIC_CSQ","COS_CODING_ID","COS_CODING_Gene_name","COS_CODING_STRAND","COS_CODING_GENOMIC_ID","COS_CODING_LEGACY_ID","COS_CODING_CDS","COS_CODING_AA","COS_CODING_HGVSC","COS_CODING_HGVSP","COS_CODING_HGVSG","COS_CODING_CNT","COS_CODING_ID","COS_CODING_Gene_name","COS_CODING_STRAND","COS_CODING_LEGACY_ID","COS_CODING_CNT","in_cgi", "cgi_mutation_type", "cgi_cancer_types", "cgi_abbreviations", "cgi_reported_by","cgi_gene", "cgi_gene_tumorigenesis", "cgi_gene_alterations", "cgi_gene_translocations", "cgi_gene_cancer_types", "cgi_gene_abbreviations", "cgi_gene_sources","brca_geneSymbol","brca_variant_nomenclature_reference_cDNA_sequence","brca_variant_nomenclature_HGVS_nucleotide","brca_variant_nomenclature_BIC_designation","brca_variant_nomenclature_genome","brca_variant_nomenclature_RNA_LOVD","brca_variant_nomenclature_synonyms","brca_variant_nomenclature_pathogenicity_all","brca_variant_nomenclature_clinical_significance","brca_variant_nomenclature_variant_haplotype","brca_ENIGMA_clinical_significance_citation","brca_ENIGMA_allel_origin","brca_ClinVar_SCV_accession","brca_ClinVar_DBID_LOVD","brca_ClinVar_BIC_Nomenclature","clinivar_ALLELEID","clinivar_CLNDN","clinivar_CLNDNINCL","clinivar_CLNDISDBINCL","clinivar_CLNHGVS","clinivar_CLNREVSTAT","clinivar_CLNSIG","clinivar_CLNSIGCONF","clinivar_CLNSIGINCL","clinivar_CLNVC","clinivar_CLNVCSO","clinivar_CLNVI","clinivar_DBVARID","clinivar_GENEINFO","clinivar_MC","clinivar_ORIGIN","clinivar_RS","clinivar_SSR"]
    vep_header = [u + "_vep" for u in csq_normal.split("|")]
    annovar_header=['ANNOVAR_DATE','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','Func.ensGene','Gene.ensGene','GeneDetail.ensGene','ExonicFunc.ensGene','AAChange.ensGene','Func.knownGene','Gene.knownGene','GeneDetail.knownGene','ExonicFunc.knownGene','AAChange.knownGene','cytoBand','avsnp150','DamagePredCount','SIFT_pred','SIFT4G_pred','Polyphen2_HDIV_pred','Polyphen2_HVAR_pred','LRT_pred','MutationTaster_pred','MutationAssessor_pred','FATHMM_pred','PROVEAN_pred','VEST4_score','MetaSVM_pred','MetaLR_pred','M-CAP_pred','REVEL_score','MutPred_score','MVP_score','MPC_score','PrimateAI_pred','DEOGEN2_pred','BayesDel_addAF_pred','BayesDel_noAF_pred','ClinPred_pred','LIST-S2_pred','CADD_raw','CADD_phred','DANN_score','fathmm-MKL_coding_pred','fathmm-XF_coding_pred','Eigen-raw_coding','Eigen-phred_coding','Eigen-PC-raw_coding','Eigen-PC-phred_coding','GenoCanyon_score','integrated_fitCons_score','GM12878_fitCons_score','H1-hESC_fitCons_score','HUVEC_fitCons_score','LINSIGHT','GERP++_NR','GERP++_RS','phyloP100way_vertebrate','phyloP30way_mammalian','phyloP17way_primate','phastCons100way_vertebrate','phastCons30way_mammalian','phastCons17way_primate','bStatistic','Interpro_domain','GTEx_V8_gene','GTEx_V8_tissue','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','1000g2015aug_all','1000g2015aug_afr','1000g2015aug_eas','1000g2015aug_eur','1000g2015aug_sas','1000g2015aug_amr','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS']
    all_header=position_info+annovar_header+vep_header+vcfanno_fields
    print(len(all_header))
    return all_header
def each(each_row):
    """
    Generate the each annotation information.

    Parameters:
    - each_row: The input each_row dataframe.

    Returns:
    - all_out: A string representing the concatenated values of the input dataframe.

    """
    test=each_row.loc['INFO']
    pattern_vep=r"(CSQ=[A-Z]\|.*\|.*;ANNOVAR_DATE)"
    pattern_annovar=r"(ANNOVAR_DATE=.*;ALLELE_END)"
    pattern_vcfanno=r"(ALLELE_END.*)"
    m = re.search(pattern_vep, test)
    print(each_row)
    n = re.search(pattern_annovar, test)
    l= re.search(pattern_vcfanno, test)
    test_re=re.match(r'(CSQ=C|.*;)',test)
    l.group()=="ALLELE_END"
    if(l.group()=="ALLELE_END"):
        b1=[" " for m in range(0,85)]

    #         vcfanno_info=[" " for m in range(len(list(range(0,79)))]
    else:
        vcfanno_info=l.group().split("ALLELE_END;")[1]
        b1=vcfanno(vcfanno_info)
    #     vcfanno_info=l.group().split("ALLELE_END;")[1]

    vep_info_1=m.group().split("CSQ=")[1]
    vep_info=vep_info_1.split(";ANNOVER_DATE")[0]
    annovar_info=n.group()

    vep_all=vep_info.split("|")[0:]
    vep_all[-1]=vep_all[-1].split(";")[0]
    position=each_row[['#CHROM','POS','REF','ALT']].tolist()
    all=[]
    a1=annovar_split(annovar_info)
    # b1=vcfanno(vcfanno_info)
    all=position+a1+vep_all+b1
    len(all)
    all_out="\t".join(map(str, all))
    return all_out




if __name__ == '__main__':
    sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())
    vcf_path=sys.argv[1]
    out_path=sys.argv[2]
    out_file=os.path.join(out_path,"vcf_normaled.tsv")
    data=pd.read_csv(vcf_path,sep="\t")
    sep="\t"
    file_write_obj = open(out_file, 'w',encoding='utf-8')
    header_first=sep.join(header())
    file_write_obj.writelines(header_first)
    file_write_obj.write('\n')
    for i in range(len(data)):
        # print(i) 
        i_each=each(data.loc[i,:])
        # print(i_each)
        file_write_obj.writelines(i_each)
        file_write_obj.write('\n')

    file_write_obj.close()