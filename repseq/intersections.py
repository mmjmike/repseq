import math
import pandas as pd
from .common_functions import print_progress_bar, run_parallel_calculation
from .io import read_clonoset
from .clonosets import filter_nonfunctional_clones, recount_fractions_for_clonoset, get_column_names_from_clonoset
from repseq.clone_filter import Filter


def intersect_clones_in_samples_batch(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=True):
    """
    Calculating frequencies of intersecting clonotypes between multiple repseq samples.
    The result of this function may be used for scatterplots of frequencies/counts of 
    overlapping clonotypes
    
    Args:
        clonosets_df (pd.DataFrame): contains three columns - 'sample_id' and 'filename' columns,
            filename - full path to clonoset file. Clonoset file may be of MiXCR3/MiXCR4 or VDJtools format
            sample_id's should be all unique in this DF
        overlap_type (str): possible values are aa, aaV, aaVJ, nt, ntV, ntVJ. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        by_umi (bool): set =True for MiXCR4 clonosets to select count/frequency of clonotypes 
            in UMI's if they exist in implemented protocol
        by_freq (bool): default is True - this means that the intersect metric is frequency of clonotype, 
            but not its count
        only_functional (bool): use only functional clonotypes (do not contain stop codons or
            frameshifts in CDR3 sequences: * or _ symbol in CDR3aa sequence). The frequences are recounted to
            1 after filtering of non-functional clonotypes
    
    Important: when using particular overlap type, similar clonotypes in one particular clonoset are
    combined into one with summation of counts/frequencies.

    Returns:
        df (pd.DataFrame): dataframe with following columns: "clone", "sample1_count", "sample2_count", "sample1", "sample2", "pair"
            clone - is tuple, containing sequence (aa or nt), plus V or J if they are required by the metric
            count columns contain freq/count of the clone in sample
            pair column is made for easy separation of possibly huge DataFrame into overlapping pairs
    """


    print("Intersecting clones in clonosets\n"+"-"*50)
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    print(f"Overlap type: {overlap_type}")
    
    if len(clonosets_df.sample_id.unique()) < len(clonosets_df):
        raise ValueError("Input clonosets DataFrame contains non-unique sample ID's")
    
    # converting clonosets to compact dicts (clone: freq) based on overlap type and count/freq/umi
    samples_total = len(clonosets_df)
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=True,
                                                        retain_counts=False)
        
    sample_list = list(clonosets_df.sort_values(by="sample_id").sample_id)
    tasks = []
    clonosets_df.sort_values(by="sample_id")["sample_id"]
    
    for i in range(samples_total):
        for j in range(samples_total-i-1):
            sample1 = sample_list[i]
            sample2 = sample_list[j+i+1]
            tasks.append((sample1, sample2, clonoset_dicts))
    
    results = run_parallel_calculation(intersect_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")
    return pd.concat(results).reset_index(drop=True)


def count_table(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0):
    print("Creating clonotypes count table\n"+"-"*50)
    print(f"Overlap type: {overlap_type}")
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    clonoset_dicts = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=False)
    unique_clonotypes = find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts, check_v, check_j)
    
    tasks = []
    for sample_id in clonoset_dicts:
        task = [unique_clonotypes, sample_id, clonoset_dicts[sample_id], check_v, check_j, mismatches]
        tasks.append(task)
    
    results = run_parallel_calculation(count_table_mp, tasks, "Counting features", object_name="pairs")
    result_dict = dict()
    for result in results:
        result_dict.update(result)
    count_table = pd.DataFrame(result_dict)
    count_table.index = unique_clonotypes
    return count_table

def count_table_mp(args):
    (features, sample_id, clonoset_dict, check_v, check_j, mismatches) = args
    result = []
    for feature in features:
        count = 0
        len_feature = len(feature[0])
        if len_feature in clonoset_dict:
            clonotypes = clonoset_dict[len_feature]
            for clonotype in clonotypes:
                if clonotypes_equal(feature, clonotype, check_v, check_j, mismatches=mismatches):
                    count += clonotype[-1]
        result.append(count)
    return {sample_id: result}
        
        


def find_unique_clonotypes_in_clonoset_dicts(clonoset_dicts, check_v, check_j):
    unique_clonotypes = set()
    for sample_id in clonoset_dicts:
        for cdr3len in clonoset_dicts[sample_id]:
            for clonotype in clonoset_dicts[sample_id][cdr3len]:
                clone_len = 1
                if check_v:
                    clone_len += 1
                if check_j:
                    clone_len += 1
                unique_clonotypes.add(tuple(clonotype[:clone_len]))
    return unique_clonotypes
    

def overlap_distances(clonosets_df, cl_filter=None, overlap_type="aaV", mismatches=0, metric="F2"):
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
        metric (str): possible values - `F` or `F2`. Default `F2`. `F2` - sum of sqrt of product of 
            similar clonotype frequencies in two clonosets. F - sqrt of the sum of frequency products.
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
    
    metric = metric.upper()
    metrics = ["F", "F2"]
    if metric not in metrics:
        raise ValueError(f"Metric {metric} is not supported. Possible values: {', '.join(metrics)}")
    
    if len(clonosets_df.sample_id.unique()) < len(clonosets_df):
        raise ValueError("Input clonosets DataFrame contains non-unique sample ID's")
    
    # converting clonosets to compact lists of clonotypes separated by CDR3 lengths to dictionary based on overlap type and count/freq/umi
    samples_total = len(clonosets_df)
    clonoset_lists = convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=cl_filter,
                                                        overlap_type=overlap_type, by_freq=True)
    
    # generating a set of tasks
    
    sample_list = list(clonosets_df.sort_values(by="sample_id").sample_id)
    tasks = []
    
    for i in range(samples_total):
        for j in range(samples_total-i-1):
            sample1 = sample_list[i]
            sample2 = sample_list[j+i+1]
            tasks.append((sample1, sample2, clonoset_lists, check_v, check_j, mismatches, metric))
    
    # run calculation in parallel
    result_list = run_parallel_calculation(overlap_metric_two_clone_dicts, tasks, "Intersecting clonosets", object_name="pairs")
    
    result_list = result_list + [(result[1], result[0], result[2]) for result in result_list]
    overlap_df = pd.DataFrame(result_list, columns=["sample1", "sample2", "f2"]).pivot_table(index="sample1", columns=["sample2"], values="f2").reset_index().set_index("sample1").fillna(1)
    return overlap_df


### Supporting functions

def convert_clonosets_to_compact_dicts(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=True, retain_counts=True):
    clonoset_dicts = {}
    
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
                                                    by_freq=by_freq, retain_counts=retain_counts)
        samples_read += 1
        print_progress_bar(samples_read, samples_total, "Reading clonosets")
        clonoset_dicts[sample_id] = cl_dict
    return clonoset_dicts


def overlap_metric_two_clone_dicts(args):
    (sample_id_1, sample_id_2, clonoset_dicts, check_v, check_j, mismatches, metric) = args
    
    f_metric = False
    if metric == "F":
        f_metric = True
    
    cl1_dict = clonoset_dicts[sample_id_1]
    cl2_dict = clonoset_dicts[sample_id_2]

    frequency = 0
    for c1_len, c1_clones in cl1_dict.items():
        if c1_len in cl2_dict:
            for c1 in c1_clones:
                for c2 in cl2_dict[c1_len]:
                    if clonotypes_equal(c1, c2, check_v, check_j, mismatches=mismatches):
                        if f_metric:
                            frequency += c1[-1]*c2[-1]
                        else:
                            frequency += math.sqrt(c1[-1]*c2[-1])
    if f_metric:
        frequency = math.sqrt(frequency)
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

def overlap_type_to_flags(overlap_type):
    possible_overlap_types = ["aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ"]
    if overlap_type not in possible_overlap_types:
        raise ValueError("Incorrect overlap type. Possible values: {}".format(", ".join(possible_overlap_types)))    
    aa = False
    if overlap_type[0:2] == "aa":
        aa = True
    check_v = False
    if "V" in overlap_type:
        check_v = True
    check_j = False
    if "J" in overlap_type:
        check_j = True
    return aa, check_v, check_j

def prepare_clonoset_for_intersection(clonoset, overlap_type="aaV", by_freq=True, retain_counts=False):
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
    
    
    if not retain_counts:
        clonoset["clone"] = clonoset.apply(lambda x: tuple(x[col] for col in result_colnames), axis=1)
        return clonoset[["clone", weight_column]].groupby("clone").sum().sort_values(by=weight_column,ascending=False).to_dict()[weight_column]
    else:
        result_colnames.append(weight_column)
        clonoset_list = list(clonoset.apply(lambda x: tuple(x[col] for col in result_colnames), axis=1))
        clonoset_dict = {}
        for clone in clonoset_list:
            seq_len = len(clone[0])
            if seq_len in clonoset_dict:
                clonoset_dict[seq_len].append(clone)
            else:
                clonoset_dict[seq_len] = [clone]
        return clonoset_dict

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