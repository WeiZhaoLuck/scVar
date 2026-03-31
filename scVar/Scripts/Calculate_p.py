import pandas as pd
import os
import numpy as np
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import argparse


def each_group_task(args):
    """
    Task for processing each group in parallel.
    """
    labels1, labels2, all_alt, all_ref, labels_all_ref, labels_all_alt = args
    data_labels_alt = all_alt
    data_labels_ref = all_ref
    resultsP = []
    results_odd = []
    data_labels_ref_out = []
    SNV = data_labels_ref.index.values.tolist()
    for a in range(len(data_labels_ref)):
        if (data_labels_alt.iloc[a, labels2] > 0 and data_labels_ref.iloc[a, labels1] > 0):
            data_each_ref = [data_labels_ref.iloc[a, labels1], data_labels_ref.iloc[a, labels2]]
            data_each_alt = [data_labels_alt.iloc[a, labels1], data_labels_alt.iloc[a, labels2]]
            res, result_p = fisher_exact([data_each_alt, data_each_ref])
            resultsP.append(result_p)
            results_odd.append(res)
            data_labels_ref_out.append(SNV[a])
    if len(resultsP) > 0:
        reject, p_adjusted, _, _ = smm.multipletests(resultsP, method='fdr_bh')

        result_all = pd.DataFrame()
        result_all['SNV_label'] = data_labels_ref_out
        result_all['p-adjusted'] = p_adjusted
        result_all['pvalue'] = resultsP
        result_all['odd_ratio'] = results_odd
        result_all['ref_label1'] = data_labels_ref.loc[data_labels_ref_out].iloc[:, labels1].tolist()
        result_all['ref_label2'] = data_labels_ref.loc[data_labels_ref_out].iloc[:, labels2].tolist()
        result_all['alt_label1'] = data_labels_alt.loc[data_labels_ref_out].iloc[:, labels1].tolist()
        result_all['alt_label2'] = data_labels_alt.loc[data_labels_ref_out].iloc[:, labels2].tolist()
        results_gene_SNV = result_all[result_all['pvalue'] < 0.05]
        name = "padjusted" + "_" + labels_all_ref[labels1] + "_" + labels_all_alt[labels2] + ".csv"
        path_out = os.path.join("results", name)
        results_gene_SNV.to_csv(path_out)


def cluster_other_task(args):
    """
    Task for processing cluster_other in parallel.
    """
    label, all_alt, all_ref, labels_all_ref, args_out, args_type = args
    data_labels_alt = all_alt.copy()
    data_labels_ref = all_ref.copy()
    resultsP = []
    results_odd = []
    data_labels_ref_out = []
    SNV = data_labels_ref.index.values.tolist()
    result_list = [element for element in labels_all_ref if element != label]
    data_labels_ref['RowSum'] = data_labels_ref[result_list].sum(axis=1)
    data_labels_alt['RowSum'] = data_labels_alt[result_list].sum(axis=1)
    for a in range(len(data_labels_ref)):
        name = SNV[a]
        if data_labels_ref.loc[name, 'RowSum'] > 0:
            data_each_ref = [data_labels_ref.loc[name, label], data_labels_ref.loc[name, 'RowSum']]
            data_each_alt = [data_labels_alt.loc[name, label], data_labels_alt.loc[name, 'RowSum']]
            res, result_p = fisher_exact([data_each_alt, data_each_ref])
            resultsP.append(result_p)
            results_odd.append(res)
            data_labels_ref_out.append(SNV[a])
    reject, p_adjusted, _, _ = smm.multipletests(resultsP, method='fdr_bh')

    result_all = pd.DataFrame()
    result_all['SNV_label'] = data_labels_ref_out
    result_all['p-adjusted'] = p_adjusted
    result_all['pvalue'] = resultsP
    result_all['odd_ratio'] = results_odd
    result_all['ref_label1'] = data_labels_ref.loc[data_labels_ref_out].loc[:, label].tolist()
    result_all['ref_label2'] = data_labels_ref.loc[data_labels_ref_out].loc[:, 'RowSum'].tolist()
    result_all['alt_label1'] = data_labels_alt.loc[data_labels_ref_out].loc[:, label].tolist()
    result_all['alt_label2'] = data_labels_alt.loc[data_labels_ref_out].loc[:, 'RowSum'].tolist()
    results_gene_SNV = result_all[result_all['pvalue'] < 0.05]
    name = "padjusted" + "_" + label + "_" + 'other' + ".csv"
    out_path = os.path.join(args_out, args_type)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    path_out = os.path.join(out_path, name)
    result_all.to_csv(path_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Example of using parameters')
    parser.add_argument('--type', type=str, default='gene', help='gene or mutation')
    parser.add_argument('--ref', type=str, default='BT1292_1/genotype_ref_cell_type.tsv', help='ref path')
    parser.add_argument('--alt', type=str, default='BT1292_1/genotype_alt_cell_type.tsv', help='alt path')
    parser.add_argument('--out', type=str, default='BT1292_1/', help='out path')
    args = parser.parse_args()

    # Load data
    all_ref = pd.read_csv(args.ref, sep=",", index_col=0)
    all_alt = pd.read_csv(args.alt, sep=",", index_col=0)
    labels_all_ref = all_ref.columns.tolist()[:]
    labels_all_alt = all_alt.columns.tolist()[:]

    # Parallel processing for cluster_other
    tasks = [(label, all_alt, all_ref, labels_all_ref, args.out, args.type) for label in labels_all_ref]

    with Pool(cpu_count() - 1) as pool:
        list(tqdm(pool.imap(cluster_other_task, tasks), total=len(tasks), desc="Processing cluster_other"))

    # Uncomment the following block for parallel processing of each_group
    # tasks = [(i, l, all_alt, all_ref, labels_all_ref, labels_all_alt)
    #          for i in range(len(labels_all_ref)) for l in range(len(labels_all_alt)) if labels_all_ref[i] != labels_all_alt[l]]
    #
    # with Pool(cpu_count() - 1) as pool:
    #     list(tqdm(pool.imap(each_group_task, tasks), total=len(tasks), desc="Processing each_group"))