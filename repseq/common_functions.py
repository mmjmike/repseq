import os
import pandas as pd
import concurrent.futures
import numpy as np
import math
from collections import OrderedDict

def print_progress_bar(samples_done, samples_total, program_name="", object_name="sample(s)"):
    if samples_total == 0:
        print(f"{program_name} |{'-' * 50}| 0/0 {object_name} processed")
        return
    total_steps = 50
    done = int(samples_done/samples_total*total_steps)
    bar = '#'*done + '-'*(total_steps-done)
    progress_bar = f"{program_name} |{bar}| {samples_done}/{samples_total} {object_name} processed"
    end = "\r"
    if samples_done == samples_total:
        end = "\n"
    print(progress_bar, end=end)
    
def run_parallel_calculation(function, tasks, program_name, object_name="tasks", verbose=True, cpu=None):
    result_list = []
    tasks_total = len(tasks)
    tasks_done = 0
    if verbose:
        if cpu == 1:
            print("Using 1 core")
        elif cpu is None:
            print("Using all available worker processes by default. Set cpu=1 or another integer to change this.")
        else:
            print(f"Using {cpu} cores")
        print_progress_bar(tasks_done, tasks_total, program_name, object_name=object_name)
    if cpu == 1:
        for task in tasks:
            result = function(task)
            result_list.append(result)
            tasks_done+=1
            if verbose:
                print_progress_bar(tasks_done, tasks_total, program_name, object_name=object_name)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
            for result in executor.map(function, tasks):
                result_list.append(result)
                tasks_done+=1
                if verbose:
                    print_progress_bar(tasks_done, tasks_total, program_name, object_name=object_name)
    return result_list

def diversity_metrics(list_of_numbers):
    counts = np.asarray(list(list_of_numbers), dtype=np.float64)
    counts = counts[counts > 0]
    total_size = counts.sum()
    diversity = int(len(counts))

    if diversity == 0 or total_size == 0:
        return {
            "diversity": 0,
            "norm_shannon_wiener": np.nan,
            "clonality": np.nan,
            "shannon_wiener": np.nan,
            "chao1": np.nan,
            "richness": 0,
            "ace": np.nan,
            "goods_coverage": np.nan,
            "d50": np.nan,
            "simpson": np.nan,
            "inverse_simpson": np.nan,
            "gini_simpson": np.nan,
            "berger_parker": np.nan,
            "gini_coefficient": np.nan,
        }

    freqs = counts / total_size
    sw = -np.sum(freqs * np.log(freqs))
    sw_norm = sw / np.log(diversity) if diversity > 1 else 0
    clonality = 1 - sw_norm
    chao1 = calc_chao1_index(counts)
    ace = calc_ace_index(counts)
    goods_coverage = calc_goods_coverage(counts)
    simpson = np.sum(freqs ** 2)
    inverse_simpson = 1 / simpson if simpson > 0 else np.nan
    gini_simpson = 1 - simpson
    berger_parker = np.max(freqs)
    gini_coefficient = calc_gini_coefficient(counts)
    d50 = calc_d50(counts)

    results = {"diversity": diversity,
               "norm_shannon_wiener": sw_norm,
               "clonality": clonality,
               "shannon_wiener": sw,
               "chao1": chao1,
               "richness": diversity,
               "ace": ace,
               "goods_coverage": goods_coverage,
               "d50": d50,
               "simpson": simpson,
               "inverse_simpson": inverse_simpson,
               "gini_simpson": gini_simpson,
               "berger_parker": berger_parker,
               "gini_coefficient": gini_coefficient}

    return results

def calc_chao1_index(counts):
    S_obs = np.sum(np.array(counts) > 0)

    f1 = np.sum(np.array(counts) == 1)
    f2 = np.sum(np.array(counts) == 2)
    
    if f2 == 0:  # To avoid division by zero
        chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    else:
        chao1 = S_obs + (f1 * (f1 - 1)) / (2 * f2)
    
    return chao1


def calc_ace_index(counts, rare_threshold=10):
    counts = np.asarray(counts, dtype=np.float64)
    counts = counts[counts > 0]
    abundant = counts > rare_threshold
    rare = counts <= rare_threshold
    s_abundant = np.sum(abundant)
    s_rare = np.sum(rare)
    n_rare = np.sum(counts[rare])
    if s_rare == 0:
        return float(s_abundant)
    if n_rare == 0:
        return np.nan
    f1 = np.sum(counts == 1)
    c_ace = 1 - f1 / n_rare
    if c_ace <= 0:
        return np.nan
    if n_rare <= 1:
        gamma_sq_ace = 0
    else:
        rare_counts = counts[rare]
        freqs = np.array([np.sum(rare_counts == i) for i in range(1, rare_threshold + 1)])
        i_values = np.arange(1, rare_threshold + 1)
        gamma_sq_ace = (
            s_rare
            / c_ace
            * np.sum(i_values * (i_values - 1) * freqs)
            / (n_rare * (n_rare - 1))
            - 1
        )
        gamma_sq_ace = max(gamma_sq_ace, 0)
    return s_abundant + s_rare / c_ace + f1 / c_ace * gamma_sq_ace


def calc_goods_coverage(counts):
    counts = np.asarray(counts, dtype=np.float64)
    counts = counts[counts > 0]
    total = counts.sum()
    if total == 0:
        return np.nan
    f1 = np.sum(counts == 1)
    return 1 - f1 / total


def calc_d50(counts):
    counts = np.asarray(counts, dtype=np.float64)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return np.nan
    ordered = np.sort(counts)[::-1]
    dominant_clones = np.searchsorted(np.cumsum(ordered), ordered.sum() * 0.5, side="left") + 1
    return dominant_clones / len(ordered)


def calc_gini_coefficient(counts):
    counts = np.sort(np.asarray(counts, dtype=np.float64))
    counts = counts[counts > 0]
    n = len(counts)
    if n == 0:
        return np.nan
    total = counts.sum()
    if total == 0:
        return np.nan
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * counts)) / (n * total) - (n + 1) / n


def extract_segment(s):
    segm = str(s).split("*")[0]
    segm = str(segm).split("(")[0]
    if segm == "nan":
        return "."
    else:
        return segm

def extract_refpoint_position(p, n, minus=False):
    pos = p.split(":")[n]
    if pos == "":
        return -1
    elif minus:
        return int(pos)-1
    else:
        return int(pos)
    
def round_down_to_2_significant(x):
    divisions = 0
    while x > 100:
        x = x/10
        divisions += 1
    return math.floor(x) * 10 ** divisions

def center_5(string):
    if len(string) <= 5:
        return string
    start = max(0, math.ceil(len(string) / 2) - 3)
    return string[start:start + 5]

def calc_insert_size(vend,dstart,dend,jstart):
    if dstart == -1:
        insert = jstart-vend-1
        if insert < 0:
            insert = 0
    else:
        vd = dstart-vend-1
        dj = jstart-dend-1
        if vd<0:
            vd = 0
        if dj<0:
            dj = 0
        insert = vd+dj
    return insert


def overlap_type_to_flags(overlap_type):
    possible_overlap_types = ["aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ", "VJ", "VJlen"]
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


def overlap_type_uses_sequence(overlap_type):
    overlap_type_to_flags(overlap_type)
    return overlap_type not in {"VJ", "VJlen"}


def jaccard_index(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    return len(intersection) / len(union)

def bray_curtis_dissimilarity(list1, list2):
    list1 = np.array(list1)
    list2 = np.array(list2)
    
    numerator = np.sum(np.abs(list1 - list2))
    denominator = np.sum(list1 + list2)
    
    return numerator / denominator

def kl_divergence(p, q, epsilon=1e-10):
    """
    Calculate the KL divergence between two probability distributions p and q.

    Args:
        p (list of float): First probability distribution.
        q (list of float): Second probability distribution.
        epsilon (float): Small value to avoid log(0) and division by zero.

    Returns:
        float: KL divergence D(P || Q)
    """
    
    p = np.array(p, dtype=np.float64, copy=True)
    q = np.array(q, dtype=np.float64, copy=True)

    # normalize
    p /= np.sum(p)
    q /= np.sum(q)

    # escape division by zero and log(0) with epsilon value
    p_safe = np.where(p == 0, epsilon, p)
    q_safe = np.where(q == 0, epsilon, q)

    result = np.sum(np.where(p != 0, p_safe * np.log(p_safe / q_safe), 0))

    return result

def jensen_shannon_divergence(p, q, epsilon=1e-10):
    p = np.array(p, dtype=np.float64, copy=True)
    q = np.array(q, dtype=np.float64, copy=True)
    
    # Normalize the distributions
    p /= np.sum(p)
    q /= np.sum(q)
    
    # Calculate the average distribution
    m = 0.5 * (p + q)
    
    # Calculate KLD for each distribution
    kld_p_m = kl_divergence(p, m, epsilon=epsilon)
    kld_q_m = kl_divergence(q, m, epsilon=epsilon)
    
    # Calculate JSD
    jsd = 0.5 * kld_p_m + 0.5 * kld_q_m
    return jsd



def _validate_metric_vectors(values1, values2):
    values1 = np.asarray(values1, dtype=np.float64)
    values2 = np.asarray(values2, dtype=np.float64)
    if values1.shape != values2.shape:
        raise ValueError("Metric vectors must have the same shape")
    if np.any(values1 < 0) or np.any(values2 < 0):
        raise ValueError("Metric vectors must contain non-negative values")
    return values1, values2


def _normalized_metric_vectors(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    sum1, sum2 = np.sum(values1), np.sum(values2)
    if sum1 == 0 or sum2 == 0:
        return values1, values2
    return values1 / sum1, values2 / sum2


def intersecting_clonotypes_count(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    return int(np.sum((values1 > 0) & (values2 > 0)))


def relative_diversity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    denominator = np.sum(values1 > 0) * np.sum(values2 > 0)
    return intersecting_clonotypes_count(values1, values2) / denominator if denominator else np.nan


def pearson_correlation(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    mask = (values1 > 0) & (values2 > 0)
    if np.sum(mask) < 2 or np.std(values1[mask]) == 0 or np.std(values2[mask]) == 0:
        return np.nan
    return float(np.corrcoef(values1[mask], values2[mask])[0, 1])


def f1_similarity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    mask = (values1 > 0) & (values2 > 0)
    return float(np.sqrt(np.sum(values1[mask]) * np.sum(values2[mask])))


def f2_similarity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    return float(np.sum(np.sqrt(values1 * values2)))


def jaccard_similarity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    union = np.sum((values1 > 0) | (values2 > 0))
    return intersecting_clonotypes_count(values1, values2) / union if union else np.nan


def jaccard_distance(values1, values2):
    return 1 - jaccard_similarity(values1, values2)


def dice_similarity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    denominator = np.sum(values1 > 0) + np.sum(values2 > 0)
    return 2 * intersecting_clonotypes_count(values1, values2) / denominator if denominator else np.nan


def dice_distance(values1, values2):
    return 1 - dice_similarity(values1, values2)


def szymkiewicz_simpson_similarity(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    denominator = min(np.sum(values1 > 0), np.sum(values2 > 0))
    return intersecting_clonotypes_count(values1, values2) / denominator if denominator else np.nan


def l1_distance(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    return float(np.sum(np.abs(values1 - values2)))


def total_variation_distance(values1, values2):
    values1, values2 = _normalized_metric_vectors(values1, values2)
    return float(0.5 * np.sum(np.abs(values1 - values2)))


def l2_distance(values1, values2):
    values1, values2 = _validate_metric_vectors(values1, values2)
    return float(np.linalg.norm(values1 - values2))


def morisita_horn_similarity(values1, values2):
    values1, values2 = _normalized_metric_vectors(values1, values2)
    denominator = np.sum(values1 ** 2) + np.sum(values2 ** 2)
    return float(2 * np.sum(values1 * values2) / denominator) if denominator else np.nan


def hellinger_distance(values1, values2):
    values1, values2 = _normalized_metric_vectors(values1, values2)
    return float(np.linalg.norm(np.sqrt(values1) - np.sqrt(values2)) / np.sqrt(2))


def decide_count_and_frac_columns(colnames, by_umi, suppress_warnings=False):
    count_column = colnames["count_column"]
    fraction_column = colnames["fraction_column"]
    if by_umi:
        if colnames["umi"]:
            count_column = colnames["umi_column"]
            fraction_column = colnames["umi_fraction_column"]
        elif not suppress_warnings:
            print("WARNING! Clonoset does not contain UMI column. Using reads for clone count instead.\nTo avoid this warning set parameter 'by_umi=False'")
    return count_column, fraction_column


def filter_by_functionality(clonoset_in, colnames=None, functional=True):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    cdr3aa_column = colnames["cdr3aa_column"]
    if functional:
        clonoset = clonoset.loc[~clonoset[cdr3aa_column].str.contains(r"\*|_", na=False)]
        clonoset = clonoset.loc[clonoset[cdr3aa_column] != ""]
    else:
        clonoset = clonoset.loc[(clonoset[cdr3aa_column].str.contains(r"\*|_", na=False)) | (clonoset[cdr3aa_column] == "")]

    return clonoset


def get_column_names_from_clonoset(clonoset, *, normalize=str.lower, strict=False):

    # all possible names for column types
    column_alias_map = OrderedDict({
        "umi_column": ["uniqueumicount", "uniquemoleculecount"],
        "umi_fraction_column": ["uniqueumifraction", "uniquemoleculefraction"],
        "count_column": ["count", "clonecount", "readcount", "read.count"],
        "fraction_column": ["freq", "clonefraction", "frequency", "readfraction"],
        "v_column": ["v", "allvhitswithscore", "bestvgene", "v_call"],
        "d_column": ["d", "alldhitswithscore", "bestdgene", "d_call"],
        "j_column": ["j", "alljhitswithscore", "bestjgene", "j_call"],
        "c_column": ["c", "allchitswithscore", "bestcgene"],
        "cdr3aa_column": ["cdr3aa", "aaseqcdr3", "cdr3.amino.acid.sequence", "junction_aa"],
        "cdr3nt_column": ["cdr3nt", "nseqcdr3", "cdr3.nucleotide.sequence", "junction"]
    })

    cols = {normalize(c): c for c in clonoset.columns}

    colnames = {}
    missing = []

    for required_name, aliases in column_alias_map.items():
        match = next((cols[normalize(a)] for a in aliases if normalize(a) in cols), None)
        if match is None:
            missing.append(required_name)
            colnames[required_name] = None
        else:
            colnames[required_name] = match

    if strict and missing:
        raise KeyError(f"Missing required columns: {missing}")

    # detect if there is a separate column for UMI counts
    colnames["umi"] = colnames["umi_column"] is not None

    return colnames


# def combine_metadata_from_folders(folders, metadata_filename="metadata.txt"):
#     if isinstance(folders, str):
#         folders = [folders]
#     list_of_metadata_dfs = []
#     for folder in folders:
#         full_paths = False
#         curr_metadata = pd.read_csv(os.path.join(folder, metadata_filename), sep="\t")
#         if "#file.name" in curr_metadata.columns:
#             if os.path.exists(curr_metadata.iloc[0]["#file.name"]):
#                 full_paths = True
#             if not full_paths:
#                 curr_metadata["#file.name"] = curr_metadata["#file.name"].apply(lambda x: os.path.join(folder, x))
#         list_of_metadata_dfs.append(curr_metadata)
#     return pd.concat(list_of_metadata_dfs)
