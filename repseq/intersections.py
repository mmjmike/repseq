import math
import pandas as pd
import numpy as np
import itertools

from statsmodels.stats.multitest import multipletests
from scipy.stats import binom, poisson



from .common_functions import (print_progress_bar, run_parallel_calculation, overlap_type_to_flags,
                               jaccard_index, bray_curtis_dissimilarity, jensen_shannon_divergence)
from .io import read_clonoset
from .clonosets import filter_nonfunctional_clones, recount_fractions_for_clonoset, get_column_names_from_clonoset
from repseq.clone_filter import Filter
from .clustering import pool_clonotypes_from_clonosets_df



def intersect_clones_in_samples_batch(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=True, clonosets_df2=None, cl_filter2=None):
    """
    Calculating frequencies of intersecting clonotypes between multiple repseq samples.
    The result of this function may be used for scatterplots of frequencies/counts of 
    overlapping clonotypes
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - `sample_id` and `filename` columns,
            `filename` - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        by_umi (bool): set `=True` for MiXCR4 clonosets to select count/frequency of clonotypes 
            in UMI's if they exist in implemented protocol
        by_freq (bool): default is `True` - this means that the intersect metric is frequency of clonotype, 
            but not its count
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
    
    Important: when using particular overlap type, similar clonotypes in one particular clonoset are
    combined into one with summation of counts/frequencies.

    Returns:
        df (pd.DataFrame): dataframe with following columns: `clone`, `sample1_count`, `sample2_count`, `sample1`, `sample2`, `pair`
            clone - is tuple, containing sequence (aa or nt), plus V or J if they are required by the metric
            count columns contain freq/count of the clone in sample
            pair column is made for easy separation of possibly huge DataFrame into overlapping pairs
    """


    print("Intersecting clones in clonosets\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    
    clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2 = prepare_clonotypes_dfs_for_intersections(clonosets_df, clonosets_df2,
                                                                                                                        cl_filter, cl_filter2,
                                                                                                                        overlap_type, by_freq=by_freq,
                                                                                                                        strict=True)
    # generating a set of tasks
    
    tasks = []
    
    if two_dataframes:
        for sample1 in sample_list:
            for sample2 in sample_list2:
                tasks.append((sample1, sample2, clonoset_lists))
    else:
        for i in range(samples_total):
            sample1 = sample_list[i]
            for j in range(samples_total-i-1):
                sample2 = sample_list[j+i+1]
                tasks.append((sample1, sample2, clonoset_lists))
    
    results = run_parallel_calculation(intersect_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")

    # df = pd.concat(results).index.set_names()
    df = pd.concat(results).reset_index(drop=True)
    df = split_tuple_clone_column(df, overlap_type)

    return df


def count_table(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0, strict_presense=False, by_freq=False):
    
    print("Creating clonotypes count table\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=by_freq, strict=not bool(mismatches))
    unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts)
    
    tasks = []
    for sample_id in clonoset_dicts:
        task = [unique_clonotypes, sample_id, clonoset_dicts[sample_id], mismatches, strict_presense]
        tasks.append(task)
    
    results = run_parallel_calculation(count_table_mp, tasks, "Counting features", object_name="clonosets")
    result_dict = dict()
    for result in results:
        result_dict.update(result)
    count_table = pd.DataFrame(result_dict)
    count_table.index = unique_clonotypes
    return count_table



def count_table_mp(args):
    (features, sample_id, clonoset_dict, mismatches, strict_presense) = args
    result = []
    for feature in features:
        count = 0
        len_feature = len(feature[0])
        
        if mismatches:
            feature_to_check = (len_feature, *feature[1:])
            if feature_to_check in clonoset_dict:
                clonotype_present = False
                for clonotype in clonoset_dict[feature_to_check]:
                    if clonotype[0] == feature[0]:
                        clonotype_present = True
                    if sum([a != b for a,b in zip(feature[0],clonotype[0])]) <= mismatches:
                        count += clonotype[-1]
                if strict_presense and not clonotype_present:
                    result.append(0)
                    continue
        else:
            if feature in clonoset_dict:
                count += clonoset_dict[feature]

        result.append(count)
    return {sample_id: result}



def count_table_by_cluster(clonosets_df, clusters_list, cl_filter=None, overlap_type="aaV", mismatches=0, by_freq=True):
    
    print("Creating clonotypes count table\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=by_freq, strict=not bool(mismatches))
    
    clonotypes_by_cluster = convert_clusters_to_clonotype_list(clusters_list, aa, check_v, check_j, mismatches)

    tasks = []
    for sample_id in clonoset_dicts:
        task = [clonotypes_by_cluster, sample_id, clonoset_dicts[sample_id], mismatches]
        tasks.append(task)
    
    results = run_parallel_calculation(count_table_by_cluster_mp, tasks, "Counting cluster presence", object_name="clonosets")
    result_dict = dict()
    for result in results:
        result_dict.update(result)
    count_table = pd.DataFrame(result_dict).reset_index().rename(columns = {"index":"feature_id"})
    count_table["feature_id"] = count_table["feature_id"].apply(lambda x: f"cluster_{x}")

    return count_table

def convert_clusters_to_clonotype_list(clusters_list, aa, check_v, check_j, mismatches):
    
    clonotypes_by_cluster = dict()
    
    for cluster in clusters_list:
        for node in cluster:
            break
        cluster_no = node.additional_properties["cluster_no"]
        cluster_clonotypes_dict = dict()
        for node in cluster:
            if aa:
                seq = node.seq_aa
            else:
                seq = node.seq_nt

            if mismatches:            
                clone = [len(seq)]
            else:
                clone = [seq]
                
            if check_v:
                clone.append(node.v)
            if check_j:
                clone.append(node.j)
            clone = tuple(clone)
            
            if mismatches:
                if clone in cluster_clonotypes_dict:
                    cluster_clonotypes_dict[clone].add(seq)
                else:
                    cluster_clonotypes_dict[clone] = {seq}
            else:
                cluster_clonotypes_dict.update({clone:1})
        clonotypes_by_cluster[cluster_no] = cluster_clonotypes_dict
        
    return clonotypes_by_cluster

def count_table_by_cluster_mp(args):
    (clonotypes_by_cluster, sample_id, clonoset_dict, mismatches) = args
    result = dict()
    for cluster_no in clonotypes_by_cluster:
        cluster_clonotypes_dict = clonotypes_by_cluster[cluster_no]
        count = 0
        for clone in clonoset_dict:
            if mismatches:
                if clone in cluster_clonotypes_dict:
                    for seq_count in clonoset_dict[clone]:
                        (seq, clone_count) = seq_count
                        for seq2 in cluster_clonotypes_dict[clone]:
                            if sum([a != b for a,b in zip(seq,seq2)]) <= mismatches:
                                count += clone_count
                                break
            else:
                if clone in cluster_clonotypes_dict:
                    count += clonoset_dict[clone]
        result[cluster_no] = count
    return {sample_id: result}



def tcrnet(clonosets_df_exp, clonosets_df_control, cl_filter=None, cl_filter_c=None, overlap_type="aaVJ", mismatches=1):
    
    print("Running TCRnet neighbour count\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")

    clonoset_exp = pool_clonotypes_from_clonosets_df(clonosets_df_exp, cl_filter=cl_filter)
    clonoset_exp_dict = prepare_clonoset_for_intersection(clonoset_exp, overlap_type=overlap_type, by_freq=False, len_vj_format=True)
    
    unique_clonotypes = [(seq_count[0], *len_vj[1:]) for len_vj, seq_counts in clonoset_exp_dict.items() for seq_count in seq_counts]
    
    clonoset_control = pool_clonotypes_from_clonosets_df(clonosets_df_control, cl_filter=cl_filter_c)
    clonoset_control_dict = prepare_clonoset_for_intersection(clonoset_control, overlap_type=overlap_type, by_freq=False, len_vj_format=True)
    

    tasks = []
    chunks = 40
    chunk_size = len(unique_clonotypes)//chunks+1
    for i in range(chunks):
        first = i*chunk_size
        last = (i+1)*chunk_size
        task = (unique_clonotypes[first:last], clonoset_exp_dict, clonoset_control_dict, mismatches)
        tasks.append(task)
        
    results = run_parallel_calculation(tcrnet_mp, tasks, "Calc neighbours (TCRnet)", object_name="parts")
    results = list(itertools.chain.from_iterable(results)) # unpack results from several workers

    df = pd.DataFrame(results, columns=["clone", "count_exp", "count_control", "group_count_exp", "group_count_control"])
    df = tcrnet_stats_calc(df)
    return df


def tcrnet_stats_calc(df):
    print("Calculating TCRnet statistics...")
    result_df = df.copy()
    result_df["fold"] = result_df.apply(lambda x: (x["count_exp"]+1)/x["group_count_exp"]/(x["count_control"]+1)*x["group_count_control"],axis=1)
    result_df["p_value_b"] = result_df.apply(lambda x: 1-binom.cdf(x["count_exp"]-1, x["group_count_exp"], x["count_control"]/(x["group_count_control"]+1)), axis=1)
    result_df["p_value_p"] = result_df.apply(lambda x: 1-poisson.cdf(x["count_exp"]-1, x["group_count_exp"]*x["count_control"]/(x["group_count_control"]+1)), axis=1)
    result_df["p_value_b_adj"] = multipletests(result_df["p_value_b"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
    result_df["p_value_p_adj"] = multipletests(result_df["p_value_p"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
    result_df["log10_b_adj"] = -np.log(result_df["p_value_b_adj"])/np.log(10)
    result_df["log10_p_adj"] = -np.log(result_df["p_value_p_adj"])/np.log(10)
    result_df["log2_fc"] = np.log(result_df["fold"])/np.log(2)
    return result_df


def tcrnet_mp(args):
    (unique_clonotypes, clonoset_exp_dict, clonoset_control_dict, mismatches) = args
    results = []

    for unique_clone in unique_clonotypes:
        compact_clone = (len(unique_clone[0]), *unique_clone[1:])
        seq1 = unique_clone[0]
        count_exp = 0
        count_control = 0
        group_count_exp = 0
        group_count_control = 0

        if compact_clone in clonoset_exp_dict:
            for seq_count2 in clonoset_exp_dict[compact_clone]:
                if sum([a != b for a,b in zip(seq1,seq_count2[0])]) <= mismatches:
                    count_exp += 1
            group_count_exp = len(clonoset_exp_dict[compact_clone])

        if compact_clone in clonoset_control_dict:
            for seq_count2 in clonoset_control_dict[compact_clone]:
                if sum([a != b for a,b in zip(seq1,seq_count2[0])]) <= mismatches:
                    count_control += 1
            group_count_control = len(clonoset_control_dict[compact_clone])
        results.append([unique_clone, count_exp, count_control, group_count_exp, group_count_control])
    return results
        


def find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts):
    unique_clonotypes = set()
    # first = True
    for sample_id, clone_groups in clonoset_dicts.items():
        # if first:
        #     first= False
        #     i = 10
        #     for clone_group, clone_subgroup in clone_groups.items():    
        #         print(clone_group, clone_subgroup)
        #         i-=1
        #         if i <0:
        #             break
        for clone_group, clone_subgroup in clone_groups.items():
            # if clone_group is int - it is cdr3len value
            # clone_subgroup is dict with keys - tuples of cdr3seq,(v),(j) - and freq as values
            if isinstance(clone_group, int):        
                for clone in clone_subgroup:
                    unique_clonotypes.add(tuple(clone[:-1]))
            # if clone_subgroup is numeric, it means strict comparison
            # and clone_subgroup is itself a clone
            elif isinstance(clone_subgroup, (int, float, complex)):
                unique_clonotypes.add(clone_group)
            # if clone_group is not int - it is (cdr3len,(v),(j)) tuple
            # clone_subgroup is list of (cdr3seq,freq) tuples
            else:
                for clone in clone_subgroup:
                    unique_clonotypes.add(tuple([clone[0]] + list(clone_group[1:])))

        # for cdr3len in clonoset_dicts[sample_id]:
        #     for clonotype in clonoset_dicts[sample_id][cdr3len]:
        #         clone_len = 1
        #         if check_v:
        #             clone_len += 1
        #         if check_j:
        #             clone_len += 1
        #         unique_clonotypes.add(tuple(clonotype[:clone_len]))
    return list(unique_clonotypes)
    

def overlap_distances(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0, metric="F2", clonosets_df2=None, cl_filter2=None):
    """
    Calculating overlap distances between multiple repseq samples using F2 of F metrics
    The result of this function may be used for heatmap+clusterization of samples or for MDS plots
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - `sample_id` and `filename` columns,
            filename - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): The permissible number of single-letter mismatches in clonotypes sequences 
            for them to be treated similar, i.e. hamming distance.
        by_umi (bool): set =True for MiXCR4 clonosets to select count/frequency of clonotypes 
            in UMI's if they exist in implemented protocol
        metric (str): possible values - `F`, `F2` or `C`. Default `F2`. `F2` - sum of sqrt of product of 
            similar clonotype frequencies in two clonosets. `F` - sqrt of the sum of frequency products.
            `C` - total frequency of clonotypes in `sample1`, that are similar to clonotypes in `sample2`
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
    
    Important: similar clonotypes by `overlap_type` in one particular clonoset will be combined into one
        clonotype with sum for count.

    Returns:
        df (pd.DataFrame): dataframe with following columns: `clone`, `sample1_count`, `sample2_count`, `sample1`, `sample2`, `pair`. 
            `clone` - is tuple, containing sequence (aa or nt), plus V or J if they are required by the metric
            count columns contain freq/count of the clone in sample
            pair column is made for easy separation of possibly huge DataFrame into overlapping pairs
    """
    
    
    print("Intersecting clones in clonosets\n"+"-"*50)
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    print(f"Overlap type: {overlap_type}")
    
    metric = metric.upper()
    metrics = ["F", "F2", "C", "BC", "J", "JSD"]
    mismatch_metrics = ["F", "C"]
    non_symmetry_metrics = ["C"]
    frequency_metrics = ["F", "F2", "C"]
    

    if metric not in metrics:
        raise ValueError(f"Metric {metric} is not supported. Possible values: {', '.join(metrics)}")
    
    if mismatches and metric not in mismatch_metrics:
        raise ValueError(f"Metric {metric} does not allow mismatches. Mismatches only possible for: {', '.join(mismatch_metrics)}")

    by_freq = metric in frequency_metrics

    clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2 = prepare_clonotypes_dfs_for_intersections(clonosets_df, clonosets_df2,
                                                                                                                        cl_filter, cl_filter2,
                                                                                                                        overlap_type, by_freq=by_freq)
    
    # generating a set of tasks
    
    tasks = []
    
    if two_dataframes:
        for sample1 in sample_list:
            for sample2 in sample_list2:
                tasks.append((sample1, sample2, clonoset_lists, mismatches, metric))        
    else:
        if metric not in non_symmetry_metrics and not two_dataframes:
            for i in range(samples_total):
                sample1 = sample_list[i]
                for j in range(samples_total-i-1):
                    sample2 = sample_list[j+i+1]
                    tasks.append((sample1, sample2, clonoset_lists, mismatches, metric))
                if metric == "F2":
                    tasks.append((sample1, sample1, clonoset_lists, mismatches, metric))
        else:
            for i in range(samples_total):
                for j in range(samples_total):
                    sample1 = sample_list[i]
                    sample2 = sample_list[j]
                    if sample1 != sample2:
                        tasks.append((sample1, sample2, clonoset_lists, mismatches, metric))
    
    
    # run calculation in parallel
    result_list = run_parallel_calculation(overlap_metric_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")
    
    if not two_dataframes and metric != "C":
        result_list = result_list + [(result[1], result[0], result[2]) for result in result_list]
    overlap_df = pd.DataFrame(result_list, columns=["sample1", "sample2", metric.lower()]).pivot_table(index="sample1", columns=["sample2"], values=metric.lower()).reset_index().set_index("sample1").fillna(1)
    return overlap_df


def find_intersecting_clonotypes(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0, metric="F2", clonosets_df2=None, cl_filter2=None):
    """
    Calculating overlap distances between multiple repseq samples using F2 of F metrics
    The result of this function may be used for heatmap+clusterization of samples or for MDS plots
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - `sample_id` and `filename` columns,
            filename - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): The permissible number of single-letter mismatches in clonotypes sequences 
            for them to be treated similar, i.e. hamming distance.
        by_umi (bool): set =True for MiXCR4 clonosets to select count/frequency of clonotypes 
            in UMI's if they exist in implemented protocol
        metric (str): possible values - `F`, `F2` or `C`. Default `F2`. `F2` - sum of sqrt of product of 
            similar clonotype frequencies in two clonosets. `F` - sqrt of the sum of frequency products.
            `C` - total frequency of clonotypes in `sample1`, that are similar to clonotypes in `sample2`
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
    
    Important: similar clonotypes by `overlap_type` in one particular clonoset are NOT combined into one
    and are treated as different clonotypes.

    Returns:
        df (pd.DataFrame): dataframe with following columns: `clone`, `sample1_count`, `sample2_count`, `sample1`, `sample2`, `pair`. 
            `clone` - is tuple, containing sequence (aa or nt), plus V or J if they are required by the metric
            count columns contain freq/count of the clone in sample
            pair column is made for easy separation of possibly huge DataFrame into overlapping pairs
    """
    
    
    print("Intersecting clones in clonosets\n"+"-"*50)
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    print(f"Overlap type: {overlap_type}")
        
    clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2 = prepare_clonotypes_dfs_for_intersections(clonosets_df, clonosets_df2,
                                                                                                                        cl_filter, cl_filter2,
                                                                                                                        overlap_type, strict=not bool(mismatches))
    
    # generating a set of tasks
    
    tasks = []
    
    if two_dataframes:
        for sample1 in sample_list:
            for sample2 in sample_list2:
                tasks.append((sample1, sample2, clonoset_lists, check_v, check_j, mismatches))        
    else:
        for i in range(samples_total):
            for j in range(samples_total):
                sample1 = sample_list[i]
                sample2 = sample_list[j]
                if sample1 != sample2:
                    tasks.append((sample1, sample2, clonoset_lists, check_v, check_j, mismatches))
    
    
    # run calculation in parallel
    result_list = run_parallel_calculation(find_overlapping_clones_in_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")
    
    return pd.concat(result_list).reset_index(drop=True)


### Supporting functions


def prepare_clonotypes_dfs_for_intersections(clonosets_df, clonosets_df2, cl_filter, cl_filter2, overlap_type, by_freq=True, strict=False):
    """
    Args:
        clonosets_df (pd.DataFrame): _description_
        clonosets_df2 (pd.DataFrame): _description_
        cl_filter (Filter): _description_
        cl_filter2 (Filter): _description_
        overlap_type (str): _description_
        by_freq (bool, optional): _description_. Defaults to True.

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        clonoset_lists (dict): dict of 
        samples_total (int): 
        two_dataframes (bool):
        sample_list (list):
        sample_list2 (list):
    """
    # output:
    ### clonoset_lists
    
    if len(clonosets_df.sample_id.unique()) < len(clonosets_df):
        raise ValueError("Input clonosets in DataFrame have non-unique sample_id's")
    clonosets_df_1 = clonosets_df[["sample_id", "filename"]]
    two_dataframes = False
    if isinstance(clonosets_df2, pd.DataFrame):
        two_dataframes = True
        if len(clonosets_df2.sample_id.unique()) < len(clonosets_df2):
            raise ValueError("Input clonosets in DataFrame2 have non-unique sample_id's")
        clonosets_df_2 = clonosets_df2[["sample_id", "filename"]]
        intersecting_sample_ids = set(clonosets_df2.sample_id.unique()).intersection(set(clonosets_df.sample_id.unique()))
        if len(intersecting_sample_ids) > 0 and cl_filter2 is not None:
            print("WARNING! Some samples have the same sample_id in two sample_df's. The second filter will be applied to common samples")
    

    # converting clonosets to compact lists of clonotypes separated by CDR3 lengths to dictionary based on overlap type and count/freq/umi
    clonoset_lists = convert_clonosets_to_compact_dicts(clonosets_df_1, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=by_freq, strict=strict)
    if two_dataframes:
        if cl_filter2 is None:
            cl_filter2 = cl_filter
        clonoset_lists_2 = convert_clonosets_to_compact_dicts(clonosets_df_2, cl_filter=cl_filter2,
                                                        overlap_type=overlap_type, by_freq=by_freq, strict=strict)
        clonoset_lists.update(clonoset_lists_2)
    
    samples_total = len(clonosets_df_1)
    if two_dataframes:
        samples_total = len(pd.concat([clonosets_df_1, clonosets_df_2]))

    sample_list = list(clonosets_df_1.sort_values(by="sample_id").sample_id)
    sample_list2 = None
    if two_dataframes:
        sample_list2 = list(clonosets_df_2.sort_values(by="sample_id").sample_id)

    return clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2


def convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=True, strict=False):
    clonoset_dicts = {}
    len_vj_format=not strict
    
    if cl_filter is None:
        cl_filter = Filter()

    samples_total = len(clonosets_df)
    samples_read = 0
    print_progress_bar(samples_read, samples_total, "Reading clonosets")
    for i, r in clonosets_df.sort_values(by="sample_id").iterrows():
        filename = r["filename"]
        sample_id = r["sample_id"]
        clonoset = read_clonoset(filename)
        clonoset = cl_filter.apply(clonoset)
        cl_dict = prepare_clonoset_for_intersection(clonoset, overlap_type=overlap_type,
                                                    by_freq=by_freq, len_vj_format=len_vj_format)
        samples_read += 1
        print_progress_bar(samples_read, samples_total, "Reading clonosets")
        clonoset_dicts[sample_id] = cl_dict
    return clonoset_dicts


def overlap_metric_two_clone_dicts(args):
    (sample_id_1, sample_id_2, clonoset_dicts, mismatches, metric) = args


    f_metric = False
    c_metric = False
    if metric == "F":
        f_metric = True
    if metric == "C":
        c_metric = True
    

    cl1_dict = clonoset_dicts[sample_id_1]
    cl2_dict = clonoset_dicts[sample_id_2]

    if metric == "J":
        return (sample_id_1, sample_id_2, jaccard_index(cl1_dict, cl2_dict))
    
    if metric == "BC" or metric == "JSD":
        clonoset_dicts_for_pair = {sample_id_1: cl1_dict,
                                   sample_id_2: cl2_dict}
        unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts_for_pair)
        counts_dict = dict()
        for sample_id in clonoset_dicts_for_pair:
            args = (unique_clonotypes, sample_id, clonoset_dicts_for_pair, mismatches, False)
            result = count_table_mp(args)
            counts_dict.update(result)
        count_table = pd.DataFrame(counts_dict)
        count_table.index = unique_clonotypes
        if metric == "BC":
            metric_value = bray_curtis_dissimilarity(count_table[sample_id_1], count_table[sample_id_2])
        if metric == "JSD":
            metric_value = jensen_shannon_divergence(count_table[sample_id_1], count_table[sample_id_2])
        return (sample_id_1, sample_id_2, metric_value)


    frequency = 0
    for c1_key, c1_seq_freq in cl1_dict.items():
        if c1_key in cl2_dict:
            for c1 in c1_seq_freq:
                for c2 in cl2_dict[c1_key]:
                    if clonotypes_equal(c1, c2, False, False, mismatches=mismatches):
                        if f_metric:
                            frequency += c1[-1]*c2[-1]
                        elif c_metric:
                            frequency += c1[-1]
                            break
                        else:
                            frequency += math.sqrt(c1[-1]*c2[-1])
    if f_metric:
        frequency = math.sqrt(frequency)

    # for c1_len, c1_clones in cl1_dict.items():
    #     if c1_len in cl2_dict:
    #         for c1 in c1_clones:
    #             for c2 in cl2_dict[c1_len]:
    #                 if clonotypes_equal(c1, c2, check_v, check_j, mismatches=mismatches):
    #                     if f_metric:
    #                         frequency += c1[-1]*c2[-1]
    #                     elif c_metric:
    #                         frequency += c1[-1]
    #                         break
    #                     else:
    #                         frequency += math.sqrt(c1[-1]*c2[-1])
    # if f_metric:
    #     frequency = math.sqrt(frequency)


    return (sample_id_1, sample_id_2, frequency)

def clonotypes_equal(clonotype_1, clonotype_2, check_v, check_j, mismatches=0):
    seq1 = clonotype_1[0]
    seq2 = clonotype_2[0]
    if len(seq1) != len(seq2):
        return False
    if check_v and clonotype_1[1] != clonotype_2[1]:
        return False
    if check_j and clonotype_1[2] != clonotype_2[2]:
        return False
    if mismatches==0:
        return seq1 == seq2
    return sum([a != b for a,b in zip(seq1,seq2)]) <= mismatches

def prepare_clonoset_for_intersection(clonoset, overlap_type="aaV", by_freq=True, len_vj_format=False, pool_clonotypes=True):
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    
    colnames = get_column_names_from_clonoset(clonoset)
    # if only_functional:
    #     clonoset = filter_nonfunctional_clones(clonoset, colnames=colnames)
    #     clonoset = recount_fractions_for_clonoset(clonoset, colnames=colnames)
    
    cl_seq_col = colnames["cdr3aa_column"]
    if not aa:
        cl_seq_col = colnames["cdr3nt_column"]
    clonoset["seq"] = clonoset[cl_seq_col]
    
    weight_column = colnames["count_column"]
    if by_freq:
        weight_column = colnames["fraction_column"]
    
    result_colnames = ["seq"]
    
    if check_v:
        cl_v_col = colnames["v_column"]
        result_colnames.append(cl_v_col)
    
    if check_j:
        cl_j_col = colnames["j_column"]
        result_colnames.append(cl_j_col)
    
    
    if not len_vj_format or pool_clonotypes:
        clonoset["clone"] = clonoset.apply(lambda x: tuple(x[col] for col in result_colnames), axis=1)
        clonoset_dict = clonoset[["clone", weight_column]].groupby("clone").sum().sort_values(by=weight_column,ascending=False).to_dict()[weight_column]
        if len_vj_format:
            clonoset_dict = clone_dict_to_len_vj_format(clonoset_dict)
    else:
        result_colnames.append(weight_column)
        clonoset_compact = clonoset[result_colnames]
        clonoset_dict = dict()
        for i,r in clonoset_compact.iterrows():
            clone = [v for v in r]
            clone_value = [clone[0], clone[-1]]
            clone_key = tuple([len(clone[0])] + clone[1:-1])
            if clone_key not in clonoset_dict:
                clonoset_dict[clone_key] = [clone_value]
            else:
                clonoset_dict[clone_key].append(clone_value)

        
    # else:
        
    #     clonoset_list = list(clonoset.apply(lambda x: tuple(x[col] for col in result_colnames), axis=1))
    #     clonoset_dict = dict()
    #     for clone in clonoset_list:
    #         seq_len = len(clone[0])
    #         if seq_len in clonoset_dict:
    #             clonoset_dict[seq_len].append(clone)
    #         else:
    #             clonoset_dict[seq_len] = [clone]
    
    # if len_vj_format:
    #     clonoset_dict = clone_dict_to_len_vj_format(clonoset_dict)
    return clonoset_dict

def clone_dict_to_len_vj_format(clone_dict):
    len_vj_dict = dict()
    for clone, count in clone_dict.items():
        new_key = (len(clone[0]), *clone[1:])
        seq = clone[0]
        clone_value = (seq, count)
        if new_key not in len_vj_dict:
            len_vj_dict[new_key] = [clone_value]
        else:
            len_vj_dict[new_key].append(clone_value)
    return len_vj_dict

def intersect_two_clone_dicts(args):
    (sample_id_1, sample_id_2, clonoset_dicts) = args
    cl1_dict = clonoset_dicts[sample_id_1]
    cl2_dict = clonoset_dicts[sample_id_2]
    all_clones = set(cl1_dict.keys()).union(set(cl2_dict.keys()))
    results = []
    for clone in all_clones:
        freq1, freq2 = 0, 0
        if clone in cl1_dict:
            freq1 = cl1_dict[clone]
        if clone in cl2_dict:
            freq2 = cl2_dict[clone]
        results.append([clone, freq1, freq2])
    clones_intersect = pd.DataFrame(results, columns = ["clone", "sample1_count", "sample2_count"])
    clones_intersect["sample1"] = sample_id_1
    clones_intersect["sample2"] = sample_id_2
    clones_intersect["pair"] = f"{sample_id_1}_vs_{sample_id_2}"
    return clones_intersect

def split_tuple_clone_column(df, overlap_type):
    clone_column = "clone"
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    
    seq_column = "cdr3nt"
    if aa:
        seq_column = "cdr3aa"
    df[seq_column] = df[clone_column].apply(lambda x: x[0])

    if check_j:
        df["j"] = df[clone_column].apply(lambda x: x[2])
        df.insert(0, "j", df.pop("j"))

    if check_v:
        df["v"] = df[clone_column].apply(lambda x: x[1])
        df.insert(0, "v", df.pop("v"))

    df.insert(0, seq_column, df.pop(seq_column))
    df = df.drop(columns=[clone_column])
    return df


def find_overlapping_clones_in_two_clone_dicts(args):
    (sample_id_1, sample_id_2, clonoset_dicts, check_v, check_j, mismatches) = args
    
    cl1_dict = clonoset_dicts[sample_id_1]
    cl2_dict = clonoset_dicts[sample_id_2]

    results = []

    if mismatches:
        for c1_len_vj, c1_clone_list in cl1_dict.items():
            if c1_len_vj in cl2_dict:
                for c1 in c1_clone_list:
                    clone_1 = (c1[0], *c1_len_vj[1:])
                    for c2 in cl2_dict[c1_len_vj]:
                        if sum([a != b for a,b in zip(c1[0],c2[0])]) <= mismatches:
                            clone_2 = (c2[0], *c1_len_vj[1:])
                            results.append([clone_1, clone_2, c1[-1], c2[-1]])    
    else:
        for c1_len, c1_clones in cl1_dict.items():
            if c1_len in cl2_dict:
                for c1 in c1_clones:
                    for c2 in cl2_dict[c1_len]:
                        if clonotypes_equal(c1, c2, check_v, check_j, mismatches=mismatches):
                            results.append([c1, c1, c1[-1], c2[-1]])
    clones_intersect = pd.DataFrame(results, columns = ["clone1", "clone2", "sample1_count", "sample2_count"])
    clones_intersect["sample1"] = sample_id_1
    clones_intersect["sample2"] = sample_id_2
    clones_intersect["pair"] = f"{sample_id_1}_vs_{sample_id_2}"
    return clones_intersect