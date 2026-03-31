import pandas as pd
import argparse
from collections import Counter
from multiprocessing import Pool, cpu_count
from tqdm import tqdm


def merge_df1_df2(df1, df2, key):
    return pd.merge(df1, df2, on=key, how='left')


def process_gene_batch(args):
    """
    Process a batch of genes in parallel.
    :param args: Tuple containing (gene_batch, merge_result, barcode_cluster_info, unique_df, available_keys_set)
    :return: List of tuples: (gene_name, alt_counts(Counter), ref_counts(Counter))
    """
    gene_batch, merge_result, barcode_cluster_info, unique_df, available_keys = args
    results = []

    for gene_name in gene_batch:
        out_data_alt = Counter()
        out_data_ref = Counter()

        # Mutations for this gene (already filtered upstream, but keep safe)
        mutations_each_gene = unique_df.loc[unique_df["gene"] == gene_name, "mutation"].tolist()
        intersection = list(available_keys.intersection(mutations_each_gene))

        if intersection:
            mutations_each_gene_data = merge_result[merge_result["key"].isin(intersection)]

            # ---- alt barcodes ----
            # barcode_x is a list per key; explode to barcodes
            alt_barcodes = mutations_each_gene_data.explode("barcode_x")["barcode_x"].dropna().tolist()
            alt_barcodes = list(set(alt_barcodes))
            if alt_barcodes:
                sum_alt = Counter(
                    barcode_cluster_info.loc[
                        barcode_cluster_info["barcode"].isin(alt_barcodes), "cluster"
                    ]
                )
                out_data_alt.update(sum_alt)

            # ---- ref barcodes ----
            ref_barcodes = mutations_each_gene_data.explode("barcode_y")["barcode_y"].dropna().tolist()
            ref_barcodes = list(set(ref_barcodes))
            if ref_barcodes:
                sum_ref = Counter(
                    barcode_cluster_info.loc[
                        barcode_cluster_info["barcode"].isin(ref_barcodes), "cluster"
                    ]
                )
                out_data_ref.update(sum_ref)

        results.append((gene_name, out_data_alt, out_data_ref))

    return results


def gene(ref_path, alt_path, mutation_path, barcode_path, out_path, type_tag):
    """
    Process genes and generate alt/ref count tables.

    MODIFIED BEHAVIOR:
    - Only output genes whose mutations exist in args.ref/args.alt (i.e., in merge_result['key']).
    - NOT all genes appearing in mutation_path.
    """
    # Load data
    genotype_results = pd.read_table(
        alt_path, sep="\t", header=None,
        names=["key", "barcode", "chrom", "pos", "ref", "alt", "ref_num", "alt_num", "all_num", "cell_type"]
    )
    ref_results = pd.read_table(
        ref_path, sep="\t", header=None,
        names=["key", "barcode", "chrom", "pos", "ref", "alt", "ref_num", "alt_num", "all_num", "cell_type"]
    )

    result_final_alt = genotype_results.groupby("key")["barcode"].apply(list).reset_index()
    result_final_ref = ref_results.groupby("key")["barcode"].apply(list).reset_index()
    merge_result = merge_df1_df2(result_final_alt, result_final_ref, "key")

    # Available mutations (keys) that actually exist in alt/ref tables
    available_keys = set(merge_result["key"].astype(str).tolist())

    # Load mutation->gene mapping, then FILTER to only mutations present in alt/ref
    mutation_Gene = pd.read_csv(mutation_path, sep="\t", header=0, names=["mutation", "gene"], index_col=False)
    mutation_Gene["mutation"] = mutation_Gene["mutation"].astype(str)
    mutation_Gene = mutation_Gene[mutation_Gene["mutation"].isin(available_keys)].copy()

    # Load barcode->cluster mapping
    barcode_cluster_info = pd.read_csv(barcode_path, sep="\t", header=None, names=["barcode", "cluster"], index_col=False)

    # After filtering, only keep genes that have >=1 mutation in alt/ref
    unique_df = mutation_Gene.drop_duplicates()
    all_gene = sorted(unique_df["gene"].dropna().unique().tolist())

    # If nothing remains after filtering, still write empty files (optional) or raise
    if len(all_gene) == 0:
        # create empty tables with cluster columns
        clusters = sorted(barcode_cluster_info["cluster"].dropna().unique().tolist())
        out_data_alt = pd.DataFrame(index=[], columns=clusters, dtype="int32")
        out_data_ref = pd.DataFrame(index=[], columns=clusters, dtype="int32")
        out_data_alt.to_csv(out_path + f"all_table_alt_gene_{type_tag}.csv", sep=",")
        out_data_ref.to_csv(out_path + f"all_table_ref_gene_{type_tag}.csv", sep=",")
        return

    # Split genes into batches
    batch_size = 100
    gene_batches = [all_gene[i:i + batch_size] for i in range(0, len(all_gene), batch_size)]

    # Prepare tasks
    tasks = [(gene_batch, merge_result, barcode_cluster_info, unique_df, available_keys) for gene_batch in gene_batches]

    # Parallel processing
    nproc = max(1, cpu_count() - 1)
    with Pool(nproc) as pool:
        results = list(
            tqdm(pool.imap(process_gene_batch, tasks), total=len(tasks), desc="Processing gene batches")
        )

    flattened_results = [item for batch_result in results for item in batch_result]

    # Build output tables (only for filtered genes)
    clusters = sorted(barcode_cluster_info["cluster"].dropna().unique().tolist())
    out_data_alt = pd.DataFrame(0, index=all_gene, columns=clusters, dtype="int32")
    out_data_ref = pd.DataFrame(0, index=all_gene, columns=clusters, dtype="int32")

    for gene_name, alt_counts, ref_counts in flattened_results:
        for cluster, count in alt_counts.items():
            out_data_alt.loc[gene_name, cluster] = int(count)
        for cluster, count in ref_counts.items():
            out_data_ref.loc[gene_name, cluster] = int(count)

    # Save results
    out_data_alt.to_csv(out_path + f"all_table_alt_gene_{type_tag}.csv", sep=",")
    out_data_ref.to_csv(out_path + f"all_table_ref_gene_{type_tag}.csv", sep=",")


def process_mutation_batch(args):
    """
    Process a batch of mutations in parallel.
    :param args: Tuple containing the mutation batch, merge_result, barcode_cluster_info.
    :return: List of tuples containing each mutation and its alt/ref counts.
    """
    mutation_batch, merge_result, barcode_cluster_info = args
    results = []

    for mutation in mutation_batch:
        out_data_alt = Counter()
        out_data_ref = Counter()

        # Process alt barcodes
        bx = merge_result.loc[merge_result["key"] == mutation, "barcode_x"]
        if len(bx) > 0:
            mutations_each_gene_data_sum_alt = bx.values[0]
            if isinstance(mutations_each_gene_data_sum_alt, list):
                sum_alt = Counter(
                    barcode_cluster_info.loc[
                        barcode_cluster_info["barcode"].isin(mutations_each_gene_data_sum_alt), "cluster"
                    ]
                )
                out_data_alt.update(sum_alt)

        # Process ref barcodes
        by = merge_result.loc[merge_result["key"] == mutation, "barcode_y"]
        if len(by) > 0:
            mutations_each_gene_data_sum_ref = by.values[0]
            if isinstance(mutations_each_gene_data_sum_ref, list):
                sum_ref = Counter(
                    barcode_cluster_info.loc[
                        barcode_cluster_info["barcode"].isin(mutations_each_gene_data_sum_ref), "cluster"
                    ]
                )
                out_data_ref.update(sum_ref)

        results.append((mutation, out_data_alt, out_data_ref))

    return results


def mutations(ref_path, alt_path, barcode_path, out_path, type_tag):
    """
    Process mutations and generate alt/ref count tables.
    """
    genotype_results = pd.read_table(
        alt_path, sep="\t", header=None,
        names=["key", "barcode", "chrom", "pos", "ref", "alt", "ref_num", "alt_num", "all_num", "cell_type"]
    )
    ref_results = pd.read_table(
        ref_path, sep="\t", header=None,
        names=["key", "barcode", "chrom", "pos", "ref", "alt", "ref_num", "alt_num", "all_num", "cell_type"]
    )

    result_final_alt = genotype_results.groupby("key")["barcode"].apply(list).reset_index()
    result_final_ref = ref_results.groupby("key")["barcode"].apply(list).reset_index()
    merge_result = merge_df1_df2(result_final_alt, result_final_ref, "key")

    barcode_cluster_info = pd.read_csv(barcode_path, sep="\t", header=None, names=["barcode", "cluster"], index_col=False)
    cell_type = list(set(barcode_cluster_info["cluster"].values.tolist()))
    SNV_all = merge_result["key"].tolist()

    batch_size = 300
    mutation_batches = [SNV_all[i:i + batch_size] for i in range(0, len(SNV_all), batch_size)]
    tasks = [(mutation_batch, merge_result, barcode_cluster_info) for mutation_batch in mutation_batches]

    out_data_alt = pd.DataFrame(0, index=SNV_all, columns=cell_type)
    out_data_ref = pd.DataFrame(0, index=SNV_all, columns=cell_type)

    nproc = max(1, cpu_count() - 1)
    with Pool(nproc) as pool:
        results = list(
            tqdm(pool.imap(process_mutation_batch, tasks), total=len(tasks), desc="Processing mutation batches")
        )

    flattened_results = [item for batch_result in results for item in batch_result]

    for mutation, alt_counts, ref_counts in flattened_results:
        for cluster, count in alt_counts.items():
            out_data_alt.loc[mutation, cluster] = count
        for cluster, count in ref_counts.items():
            out_data_ref.loc[mutation, cluster] = count

    out_data_alt.to_csv(out_path + f"all_table_alt_mutation_{type_tag}.csv", sep=",")
    out_data_ref.to_csv(out_path + f"all_table_ref_mutation_{type_tag}.csv", sep=",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example of using parameters")
    parser.add_argument("--type1", type=str, default="gene", help="gene or mutation")
    parser.add_argument("--type2", type=str, default="BT1292_1/", help="cell_types or clusters")
    parser.add_argument("--ref", type=str, default="BT1292_1/genotype_ref_cell_type.tsv", help="ref path")
    parser.add_argument("--alt", type=str, default="BT1292_1/genotype_alt_cell_type.tsv", help="alt path")
    parser.add_argument("--mutation", type=str, default="BT1292_1/mutaion.info", help="mutation path")
    parser.add_argument("--barcode", type=str, default="BT1292_1/barcodes_cell_type.tsv", help="barcode path")
    parser.add_argument("--out", type=str, default="BT1292_1/", help="out path")
    args = parser.parse_args()

    if args.type1 == "gene":
        gene(args.ref, args.alt, args.mutation, args.barcode, args.out, args.type2)
    elif args.type1 == "mutation":
        mutations(args.ref, args.alt, args.barcode, args.out, args.type2)
    else:
        print("please input gene or mutations")
