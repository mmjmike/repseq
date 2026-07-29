import math
import pandas as pd
import numpy as np
import itertools
import warnings

from statsmodels.stats.multitest import multipletests
from scipy.stats import binom, poisson



from .common_functions import (print_progress_bar, run_parallel_calculation, overlap_type_to_flags,
                               jaccard_index, bray_curtis_dissimilarity, jensen_shannon_divergence,
                               overlap_type_uses_sequence)
from .io import read_clonoset
from .clonosets import get_column_names_from_clonoset, pool_clonotypes_from_clonosets_df
from repseq.clone_filter import Filter



def intersect_clones_in_samples_batch(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=None,
                                      clonosets_df2=None, cl_filter2=None, cpu=None):
    """Build a count-and-frequency table for every requested sample pair.

    Clonotypes are always pooled and intersected using raw counts. Frequencies
    are calculated after each pairwise union is known, so they sum to one for
    each sample within each pair.

    Args:
        clonosets_df (pd.DataFrame): Table with unique ``sample_id`` values and
            clonoset ``filename`` paths.
        cl_filter (Filter, optional): Filter applied to the first sample table.
        overlap_type (str): One of ``aa``, ``aaV``, ``aaVJ``, ``nt``, ``ntV``,
            ``ntVJ``, ``VJ``, or ``VJlen``.
        by_freq (bool, optional): Deprecated compatibility argument. Its value
            is ignored because the function always intersects raw counts.
        clonosets_df2 (pd.DataFrame, optional): Optional second sample table for
            rectangular pairwise comparisons.
        cl_filter2 (Filter, optional): Filter applied to the second sample table.
        cpu (int, optional): Number of worker processes. ``None`` uses the
            executor default.

    Returns:
        pd.DataFrame: Full pairwise union table. Feature columns are followed by
        ``sample1_count``, ``sample2_count``, ``sample1_freq``, ``sample2_freq``,
        ``sample1``, ``sample2``, and ``pair``.
    """
    if by_freq is not None:
        warnings.warn(
            "`by_freq` is deprecated and ignored; "
            "intersect_clones_in_samples_batch always returns counts and derived frequencies.",
            DeprecationWarning,
            stacklevel=2,
        )
    print("Intersecting clones in clonosets\n" + "-" * 50)
    print(f"Overlap type: {overlap_type}")

    clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2 = prepare_clonotypes_dfs_for_intersections(
        clonosets_df,
        clonosets_df2,
        cl_filter,
        cl_filter2,
        overlap_type,
        by_freq=False,
        strict=True,
    )

    tasks = []
    if two_dataframes:
        for sample1 in sample_list:
            for sample2 in sample_list2:
                tasks.append((sample1, sample2, clonoset_lists))
    else:
        for i in range(samples_total):
            sample1 = sample_list[i]
            for j in range(samples_total - i - 1):
                sample2 = sample_list[j + i + 1]
                tasks.append((sample1, sample2, clonoset_lists))

    results = run_parallel_calculation(
        intersect_two_clone_dicts,
        tasks,
        "Intersecting clonosets",
        object_name="pairs",
        cpu=cpu,
    )
    df = pd.concat(results).reset_index(drop=True)
    df = split_tuple_clone_column(df, overlap_type)
    df.attrs["sample_list"] = sample_list
    df.attrs["sample_list2"] = sample_list2
    return df


def similarity(
    clonosets_df,
    cl_filter=None,
    overlap_type="aaV",
    by_freq=None,
    clonosets_df2=None,
    cl_filter2=None,
    mismatches=1,
    result="freq",
    cpu=None,
):
    """Calculate directional clonotype similarity between repertoire samples.

    The S-metric is the total count or frequency of clonotypes in a target
    sample that are similar to at least one clonotype in a comparison sample.
    Matrix rows are target samples and columns are comparison samples. A target
    clonotype contributes once even when it matches several comparison
    clonotypes; ``result="table"`` retains every matching clonotype pair.

    Similarity definitions are selected with ``overlap_type``: ``aa`` and
    ``nt`` compare CDR3 sequences, ``aaV``/``ntV`` additionally require equal V
    segments, ``aaVJ``/``ntVJ`` require equal V and J segments, ``VJ`` ignores
    CDR3 sequence, and ``VJlen`` additionally requires equal amino-acid CDR3
    length. ``mismatches`` is the maximum Hamming distance for sequence-based
    overlap types.

    Parameters
    ----------
    clonosets_df : pandas.DataFrame
        Target sample table with unique ``sample_id`` and ``filename`` columns.
    cl_filter : Filter, optional
        Filter applied to target samples.
    overlap_type : str, default "aaV"
        One of ``aa``, ``aaV``, ``aaVJ``, ``nt``, ``ntV``, ``ntVJ``, ``VJ``,
        or ``VJlen``.
    by_freq : bool, optional
        Deprecated compatibility argument. The value is ignored.
    clonosets_df2 : pandas.DataFrame, optional
        Comparison sample table. When omitted, all target samples are compared
        directionally with each other.
    cl_filter2 : Filter, optional
        Filter applied to comparison samples. Defaults to ``cl_filter``.
    mismatches : int, default 1
        Maximum CDR3 Hamming distance for sequence-based comparisons.
    result : {"freq", "count", "number", "table"}, default "freq"
        Frequency S-metric, count S-metric, number of matched target
        clonotypes, or the full table of matching clonotype pairs.
    cpu : int, optional
        Number of worker processes.

    Returns
    -------
    pandas.DataFrame
        Directional sample matrix or a matching-clonotype pair table.
    """
    if by_freq is not None:
        warnings.warn(
            "`by_freq` is deprecated and ignored; use result='freq' or "
            "result='count' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    overlap_type_to_flags(overlap_type)
    if not isinstance(mismatches, (int, np.integer)) or isinstance(
        mismatches,
        bool,
    ):
        raise TypeError("mismatches must be a non-negative integer")
    if mismatches < 0:
        raise ValueError("mismatches must be a non-negative integer")
    result = str(result).casefold()
    possible_results = {"freq", "count", "number", "table"}
    if result not in possible_results:
        raise ValueError(
            "result must be one of: " + ", ".join(sorted(possible_results))
        )

    uses_sequence = overlap_type_uses_sequence(overlap_type)
    effective_mismatches = int(mismatches) if uses_sequence else 0
    print("Calculating clonotype similarity\n" + "-" * 50)
    print(f"Overlap type: {overlap_type}")
    if uses_sequence:
        print(f"Maximum mismatches: {effective_mismatches}")

    clonoset_lists, _, two_dataframes, sample_list, sample_list2 = (
        prepare_clonotypes_dfs_for_intersections(
            clonosets_df,
            clonosets_df2,
            cl_filter,
            cl_filter2,
            overlap_type,
            by_freq=False,
            strict=not uses_sequence,
        )
    )
    column_samples = sample_list2 if two_dataframes else sample_list

    if result == "table":
        tasks = [
            (
                target_sample,
                comparison_sample,
                clonoset_lists,
                overlap_type,
                effective_mismatches,
                "table",
            )
            for target_sample in sample_list
            for comparison_sample in column_samples
            if two_dataframes or target_sample != comparison_sample
        ]
        pair_tables = (
            run_parallel_calculation(
                _similarity_pair_worker,
                tasks,
                "Calculating similarity",
                object_name="pairs",
                cpu=cpu,
            )
            if tasks
            else []
        )
        table = (
            pd.concat(pair_tables, ignore_index=True)
            if pair_tables
            else pd.DataFrame(columns=_similarity_table_columns())
        )
        table.attrs["sample_list"] = sample_list
        table.attrs["sample_list2"] = sample_list2
        table.attrs["result"] = "table"
        table.attrs["overlap_type"] = overlap_type
        table.attrs["mismatches"] = effective_mismatches
        return table

    tasks = [
        (
            target_sample,
            comparison_sample,
            clonoset_lists,
            overlap_type,
            effective_mismatches,
            result,
        )
        for target_sample in sample_list
        for comparison_sample in column_samples
    ]
    values = run_parallel_calculation(
        _similarity_pair_worker,
        tasks,
        "Calculating similarity",
        object_name="pairs",
        cpu=cpu,
    )
    matrix = pd.DataFrame(0.0, index=sample_list, columns=column_samples)
    for target_sample, comparison_sample, value in values:
        matrix.loc[target_sample, comparison_sample] = value
    if result in {"count", "number"}:
        matrix = matrix.astype(int)
    matrix.index.name = "sample1"
    matrix.columns.name = "sample2"
    matrix.attrs["sample_list"] = sample_list
    matrix.attrs["sample_list2"] = sample_list2
    matrix.attrs["result"] = result
    matrix.attrs["overlap_type"] = overlap_type
    matrix.attrs["mismatches"] = effective_mismatches
    return matrix


def count_table(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0, strict_presence=False, by_freq=False):
    """
    Creates a table that shows how many times each unique clonotype appears across different clonosets. It processes a given dataset of clonotypes (clonosets_df) 
    and generates a frequency/count table based on a specified overlap type.
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - `sample_id` and `filename` columns,
            `filename` - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`, `VJ`, and `VJlen`.
            aa/nt define which CDR3 sequence to use (amino acid or nucleotide). V/J in the overlap_type define whether
            to check V or J segments to decide if clonotypes are equal. `VJ` and `VJlen` do not compare sequences.
        mismatches (int): Max number of single-letter mismatches in clonotypes sequences 
            for them to be treated similar, i.e. hamming distance.
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
        strict_presence (bool, default:): if set to `True` and the clonotype is not found in the clonoset, it will not be counted, even when the `mismatches` option is not set to 0.
            If `False`, mismatched sequences are counted even if the exact match does not exist.
                by_freq (bool): default is `True` - this means that the intersect metric is frequency of clonotype, 
            but not its count
        by_freq (bool): default is `True` - this means that the intersect metric is frequency of clonotype, 
            but not its count
    
    Returns:
        df (pd.DataFrame): dataframe containing a pipe-delimited `clonotype` column, its component columns,
            and one count or frequency column per sample.
    """

    
    print("Creating clonotypes count table\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    effective_mismatches = mismatches if overlap_type_uses_sequence(overlap_type) else 0
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=by_freq,
                                                        strict=not bool(effective_mismatches))
    unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts)
    
    tasks = []
    for sample_id in clonoset_dicts:
        task = [unique_clonotypes, sample_id, clonoset_dicts[sample_id], effective_mismatches, strict_presence]
        tasks.append(task)
    
    results = run_parallel_calculation(count_table_mp, tasks, "Counting features", object_name="clonosets")
    result_dict = dict()
    for result in results:
        result_dict.update(result)
    count_table = pd.DataFrame(result_dict)
    count_table.insert(0, "clone", unique_clonotypes)
    return format_clonotype_columns(count_table, overlap_type)



def count_table_mp(args):
    (features, sample_id, clonoset_dict, mismatches, strict_presence) = args
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
                if strict_presence and not clonotype_present:
                    result.append(0)
                    continue
        else:
            if feature in clonoset_dict:
                count += clonoset_dict[feature]

        result.append(count)
    return {sample_id: result}



def count_table_by_cluster(clonosets_df, clusters_list, cl_filter=None, overlap_type="aaV", mismatches=0, by_freq=True):
    
    """
    This function creates a table that shows the presence of clonotypes grouped into user-provided clusters across different clonosets. Instead of counting individual clonotypes, it 
    calculates how many clonotypes from each cluster appear in each clonoset.
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - `sample_id` and `filename` columns,
            filename - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        cluster_list (?): description
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): Max number of single-letter mismatches in clonotypes sequences 
            for them to be treated similar, i.e. hamming distance.
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
    
    Returns:
        df (pd.DataFrame): dataframe with the following columns: description
    """
    
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
    
    """
    This is an implementation of TCRnet (TCR neighbour enrichment test) algorithm.  It identifies similar clonotypes for the experimental dataset and the control one based on sequence 
    similarity (allowing up to `mismatches` differences).    
    Args:
        clonosets_df_exp (pd.DataFrame): a DataFrame with experimental clonosets containing three columns - `sample_id` and `filename` columns,
            filename - full path to a clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
        clonosets_df_control (pd.DataFrame): a DataFrame with control clonosets containing three columns - `sample_id` and `filename` columns,
            filename - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
        cl_filter (repseq.clone_filter.Filter): A filter applied to the experimental dataset before processing
        cl_filter_c (repseq.clone_filter.Filter): A filter applied to the control dataset before processing
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): Max number of single-letter mismatches in clonotype sequences 
            for them to be treated similar, i.e. hamming distance.
    
    Returns:
        df (pd.DataFrame): dataframe beginning with a pipe-delimited `clonotype` column and its component columns,
            followed by neighbour counts and statistics. `p` in `p_value` stands for `poisson`, `b` for `binomial`,
            `adj` for multiple testing correction, and `log2_fc` for log2 fold change.
    """
    
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
    df = format_clonotype_columns(df, overlap_type)
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
        mismatches (int): Max number of single-letter mismatches in clonotypes sequences 
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
        cl_filter (Filter): clonoset filter - object from `clone_filter.py` module. It is applied to clonosets in `clonosets_df`.
        overlap_type (str): possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): Max number of single-letter mismatches in clonotypes sequences 
            for them to be treated similar, i.e. hamming distance.
        by_umi (bool): set =True for MiXCR4 clonosets to select count/frequency of clonotypes 
            in UMI's if they exist in implemented protocol
        metric (str): possible values - `F`, `F2` or `C`. Default `F2`. `F2` - sum of sqrt of product of 
            similar clonotype frequencies in two clonosets. `F` - sqrt of the sum of frequency products.
            `C` - total frequency of clonotypes in `sample1`, that are similar to clonotypes in `sample2`
        clonosets_df2 (pd.DataFrame): If `clonosets_df2` is None (default), samples within `clonosets_df` are compared with each other; 
            Otherwise, the comparison is performed exclusively between samples from `clonosets_df` and `clonosets_df2`. 
        cl_filter2 (Filter): clonoset filter - object from `clone_filter.py` module. It is applied to clonosets in `clonosets_df2`. 
            If there are samples with non-unique sample_ids between the two dataframes, both filters will be applied to those samples.


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


#
#
#
# In development
#
# 
# 

def intersect_clonosets_with_clonotype_list(clonosets_df, clones_df, cl_filter=None, overlap_type="aaV", mismatches=0, mode='counts', strict_presence=False, by_freq=False):
    if mode == 'counts':
        return count_table_with_custom_clonotypes(clonosets_df, 
                                           clones_df, 
                                           cl_filter=cl_filter, 
                                           overlap_type=overlap_type, 
                                           mismatches=mismatches, 
                                           strict_presence=strict_presence, 
                                           by_freq=by_freq)
    elif mode == 'clonotypes':
        return find_intersecting_clonotypes_with_custom_clonotypes(clonosets_df, 
                                                  clones_df, 
                                                  cl_filter=cl_filter, 
                                                  overlap_type=overlap_type, 
                                                  mismatches=mismatches, 
                                                  strict_presence=strict_presence)
      

def find_intersecting_clonotypes_with_custom_clonotypes(clonosets_df, clones_df, cl_filter=None, 
                                                        clonosets_df2=None, cl_filter2=None, overlap_type="aaV", 
                                                        mismatches=0, mode="clonotypes", strict_presence=True):
    
    print("Intersecting clones in clonosets\n"+"-"*50)
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    print(f"Overlap type: {overlap_type}")

    # generating a set of tasks
    tasks = []

    clonoset_lists, samples_total, two_dataframes, sample_list, sample_list2 = prepare_clonotypes_dfs_for_intersections(clonosets_df, clonosets_df2=clonosets_df2,
                                                                                                                    cl_filter=cl_filter, cl_filter2=cl_filter2,
                                                                                                                    overlap_type=overlap_type, strict=not bool(mismatches))
    # removing clonosets not present in clonoset_lists
    
    # ! intersection between clones_df and clonosets_df
    if mode == 'clonotypes':
        # creating dummy counts
        clones_df_copy = clones_df.copy()
        clones_df_copy['count'] = [0 for _ in range(clones_df.shape[0])]
        clones_df_copy['freq'] = [0 for _ in range(clones_df.shape[0])]
        clonoset_lists_with_custom_clonotypes = clonoset_lists

        clonoset_lists_with_custom_clonotypes['clonotype_list'] = prepare_clonoset_for_intersection(clones_df_copy, 
                                                                                                    overlap_type=overlap_type,
                                                                                               len_vj_format=bool(mismatches))        

        sample_list2 = ['clonotype_list']
        for sample1 in sample_list:
            for sample2 in sample_list2:
                tasks.append((sample2, sample1, clonoset_lists_with_custom_clonotypes, check_v, check_j, mismatches))    
    
    # run calculation in parallel
    result_list = run_parallel_calculation(find_overlapping_clones_in_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")
    
    # mismatches for the final dataframe
    result_df = pd.concat([df for df in result_list if not df.empty]).reset_index(drop=True)
    result_df["mismatches"] = [sum(a != b for a, b in zip(cl1[0], cl2[0]))
                                for cl1, cl2 in zip(result_df["clone1"], result_df["clone2"])]
    
    # removing extra columns
    result_df = result_df.rename(columns={'clone1': 'clonotype',
                      'sample2_count': 'count', # or freq depending on parameter
                      'sample2': 'sample_id',
                      'clone2': 'sample_clonotype'}).drop(columns=['pair', 'sample1', 'sample1_count'])
    return result_df


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


#
#
# in development
#
#


def count_table_with_custom_clonotypes(clonosets_df, clones_df=None, cl_filter=None, overlap_type="aaV", mismatches=0, strict_presence=False, by_freq=False):
    """
    Args:
        clones_df (pd.DataFrame): a dataframe with clonotypes to look for in `clonosets_df`. Should have the following columns:
        `cdr3aa`, `cdr3nt` (at least one of these two), `v`, `j` (could be none of those). For instance, your input might have only `cdr3aa` column.
        If `overlap_type` isn't specified, it is determined individually for each clonotype in `clones_df`. In this case, if a row corresponding to
        a particular clonotype has non-empty values in both `cdr3aa` and `cdr3nt` columns, `cdr3aa` is chosen.
    """

    print("Creating clonotypes count table\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=by_freq, strict=not bool(mismatches))

    if clones_df is not None:
        # creating dummy counts
        clones_df_copy = clones_df.copy()
        clones_df_copy['count'] = [0 for _ in range(clones_df.shape[0])]
        clones_df_copy['freq'] = [0 for _ in range(clones_df.shape[0])]
        clonotype_list_dict = {'clonotype_list': None}

        clonotype_list_dict['clonotype_list'] = prepare_clonoset_for_intersection(clones_df_copy, 
                                                                                overlap_type=overlap_type,
                                                                                len_vj_format=bool(mismatches))  
                           
        unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonotype_list_dict)
    else: 
        unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts)

    tasks = []
    for sample_id in clonoset_dicts:
        task = [unique_clonotypes, sample_id, clonoset_dicts[sample_id], mismatches, strict_presence]
        tasks.append(task)
    
    results = run_parallel_calculation(count_table_mp, tasks, "Counting features", object_name="clonosets")
    result_dict = dict()
    for result in results:
        result_dict.update(result)
    count_table = pd.DataFrame(result_dict)
    if aa:
        count_table.insert(0, 'cdr3aa', [ct[0] for ct in unique_clonotypes])
    else:
        count_table.insert(0, 'cdr3nt', [ct[0] for ct in unique_clonotypes])
    if check_v:
        count_table.insert(1, 'v', [ct[1] for ct in unique_clonotypes])
    if check_j:
        count_table.insert(1, 'j', [ct[2] for ct in unique_clonotypes])
    count_table.index = unique_clonotypes
    return count_table



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
    uses_sequence = overlap_type_uses_sequence(overlap_type)

    colnames = get_column_names_from_clonoset(clonoset)
    clonoset = clonoset.copy()
    weight_column = colnames["fraction_column"] if by_freq else colnames["count_column"]

    result_colnames = []
    if uses_sequence:
        cl_seq_col = colnames["cdr3aa_column"] if aa else colnames["cdr3nt_column"]
        clonoset["seq"] = clonoset[cl_seq_col]
        result_colnames.append("seq")
    if check_v:
        result_colnames.append(colnames["v_column"])
    if check_j:
        result_colnames.append(colnames["j_column"])
    if overlap_type == "VJlen":
        clonoset["cdr3_len"] = clonoset[colnames["cdr3aa_column"]].str.len()
        result_colnames.append("cdr3_len")

    clonoset["clone"] = clonoset.apply(lambda row: tuple(row[col] for col in result_colnames), axis=1)
    clonoset_dict = (clonoset[["clone", weight_column]].groupby("clone").sum()
                     .sort_values(by=weight_column, ascending=False).to_dict()[weight_column])

    if not uses_sequence:
        return clonoset_dict
    if len_vj_format and pool_clonotypes:
        return clone_dict_to_len_vj_format(clonoset_dict)
    if not len_vj_format:
        return clonoset_dict

    compact_columns = result_colnames + [weight_column]
    clonoset_dict = {}
    for _, row in clonoset[compact_columns].iterrows():
        clone = [value for value in row]
        clone_value = [clone[0], clone[-1]]
        clone_key = tuple([len(clone[0])] + clone[1:-1])
        clonoset_dict.setdefault(clone_key, []).append(clone_value)
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

def _similarity_table_columns():
    return [
        "clone1",
        "clone2",
        "sample1_count",
        "sample2_count",
        "sample1_freq",
        "sample2_freq",
        "mismatches",
        "sample1",
        "sample2",
        "pair",
    ]


def _similarity_total_count(clonoset_dict, uses_sequence):
    if not uses_sequence:
        return float(sum(clonoset_dict.values()))
    return float(
        sum(clone[-1] for clones in clonoset_dict.values() for clone in clones)
    )


def _similarity_sequence_matches(target_dict, comparison_dict, mismatches):
    for target_key, target_clones in target_dict.items():
        comparison_clones = comparison_dict.get(target_key, [])
        if not comparison_clones:
            continue
        clone_suffix = target_key[1:]
        for target_clone in target_clones:
            target_sequence = target_clone[0]
            target_identity = (target_sequence, *clone_suffix)
            target_count = target_clone[-1]
            for comparison_clone in comparison_clones:
                comparison_sequence = comparison_clone[0]
                distance = sum(
                    first != second
                    for first, second in zip(target_sequence, comparison_sequence)
                )
                if distance <= mismatches:
                    yield (
                        target_identity,
                        (comparison_sequence, *clone_suffix),
                        target_count,
                        comparison_clone[-1],
                        distance,
                    )


def _similarity_exact_matches(target_dict, comparison_dict):
    for target_clone, target_count in target_dict.items():
        if target_clone in comparison_dict:
            yield (
                target_clone,
                target_clone,
                target_count,
                comparison_dict[target_clone],
                0,
            )


def _similarity_pair_matches(
    target_dict,
    comparison_dict,
    overlap_type,
    mismatches,
):
    if overlap_type_uses_sequence(overlap_type):
        return _similarity_sequence_matches(
            target_dict,
            comparison_dict,
            mismatches,
        )
    return _similarity_exact_matches(target_dict, comparison_dict)


def _similarity_pair_worker(args):
    (
        target_sample,
        comparison_sample,
        clonoset_dicts,
        overlap_type,
        mismatches,
        result,
    ) = args
    target_dict = clonoset_dicts[target_sample]
    comparison_dict = clonoset_dicts[comparison_sample]
    uses_sequence = overlap_type_uses_sequence(overlap_type)
    target_total = _similarity_total_count(target_dict, uses_sequence)
    comparison_total = _similarity_total_count(comparison_dict, uses_sequence)
    matches = _similarity_pair_matches(
        target_dict,
        comparison_dict,
        overlap_type,
        mismatches,
    )

    if result == "table":
        rows = [
            [
                target_clone,
                comparison_clone,
                target_count,
                comparison_count,
                target_count / target_total if target_total else 0,
                comparison_count / comparison_total if comparison_total else 0,
                distance,
                target_sample,
                comparison_sample,
                f"{target_sample}_vs_{comparison_sample}",
            ]
            for (
                target_clone,
                comparison_clone,
                target_count,
                comparison_count,
                distance,
            ) in matches
        ]
        return pd.DataFrame(rows, columns=_similarity_table_columns())

    matched_target_counts = {}
    for target_clone, _, target_count, _, _ in matches:
        matched_target_counts.setdefault(target_clone, target_count)
    if result == "number":
        value = len(matched_target_counts)
    else:
        value = float(sum(matched_target_counts.values()))
        if result == "freq":
            value = value / target_total if target_total else 0.0
    return target_sample, comparison_sample, value


def intersect_two_clone_dicts(args):
    (sample_id_1, sample_id_2, clonoset_dicts) = args
    cl1_dict = clonoset_dicts[sample_id_1]
    cl2_dict = clonoset_dicts[sample_id_2]
    total1 = sum(cl1_dict.values())
    total2 = sum(cl2_dict.values())
    all_clones = set(cl1_dict.keys()).union(set(cl2_dict.keys()))
    results = []
    for clone in all_clones:
        count1 = cl1_dict.get(clone, 0)
        count2 = cl2_dict.get(clone, 0)
        freq1 = count1 / total1 if total1 else 0
        freq2 = count2 / total2 if total2 else 0
        results.append([clone, count1, count2, freq1, freq2])
    clones_intersect = pd.DataFrame(
        results,
        columns=["clone", "sample1_count", "sample2_count", "sample1_freq", "sample2_freq"],
    )
    clones_intersect["sample1"] = sample_id_1
    clones_intersect["sample2"] = sample_id_2
    clones_intersect["pair"] = f"{sample_id_1}_vs_{sample_id_2}"
    return clones_intersect

def split_tuple_clone_column(df, overlap_type):
    clone_column = "clone"
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)

    if overlap_type in {"VJ", "VJlen"}:
        df["v"] = df[clone_column].apply(lambda clone: clone[0])
        df["j"] = df[clone_column].apply(lambda clone: clone[1])
        if overlap_type == "VJlen":
            df["len"] = df[clone_column].apply(lambda clone: clone[2]).astype(int)
        feature_columns = ["v", "j"] + (["len"] if overlap_type == "VJlen" else [])
        return df[feature_columns + [column for column in df.columns if column not in feature_columns + [clone_column]]]

    seq_column = "cdr3aa" if aa else "cdr3nt"
    df[seq_column] = df[clone_column].apply(lambda clone: clone[0])
    if check_v:
        df["v"] = df[clone_column].apply(lambda clone: clone[1])
    if check_j:
        df["j"] = df[clone_column].apply(lambda clone: clone[2])
    feature_columns = [seq_column] + (["v"] if check_v else []) + (["j"] if check_j else [])
    return df[feature_columns + [column for column in df.columns if column not in feature_columns + [clone_column]]]


def format_clonotype_columns(df, overlap_type, clone_column="clone"):
    """Replace a tuple clone column with string and component columns."""
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    clonotypes = df[clone_column]

    if overlap_type in {"VJ", "VJlen"}:
        component_columns = ["v", "j"]
        if overlap_type == "VJlen":
            component_columns.append("len")
    else:
        component_columns = ["cdr3aa" if aa else "cdr3nt"]
        if check_v:
            component_columns.append("v")
        if check_j:
            component_columns.append("j")

    clone_position = df.columns.get_loc(clone_column)
    result = df.drop(columns=clone_column).copy()
    result.insert(
        clone_position,
        "clonotype",
        clonotypes.apply(lambda clone: "|".join(str(value) for value in clone)),
    )
    for offset, component_column in enumerate(component_columns, start=1):
        values = clonotypes.apply(lambda clone, index=offset - 1: clone[index])
        if component_column == "len":
            values = values.astype(int)
        result.insert(clone_position + offset, component_column, values)
    return result


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
                    # feature_to_check = (c1_len_vj, *c1[1:])
                    for c2 in cl2_dict[c1_len_vj]:
                        # check if present in cl2_dict (corresponds to those extracted from clonosets_df)
                        # cl1_dict clonotypes are from clones_df
                        # if c1[0] == c2[0]:
                        # clonotype_present = True
                        if sum([a != b for a,b in zip(c1[0],c2[0])]) <= mismatches:
                            # if strict_presence and not clonotype_present:
                                 # continue
                            # else
                            clone_2 = (c2[0], *c1_len_vj[1:])
                            results.append([clone_1, clone_2, c1[-1], c2[-1]])    
    else:
        for c1, c1_count in cl1_dict.items():
            if c1 in cl2_dict:
                results.append([c1, c1, c1_count, cl2_dict[c1]])
    clones_intersect = pd.DataFrame(results, columns = ["clone1", "clone2", "sample1_count", "sample2_count"])
    clones_intersect["sample1"] = sample_id_1
    clones_intersect["sample2"] = sample_id_2
    clones_intersect["pair"] = f"{sample_id_1}_vs_{sample_id_2}"
    return clones_intersect