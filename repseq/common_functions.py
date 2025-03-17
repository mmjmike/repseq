import os
import pandas as pd
import concurrent.futures
import numpy as np
import math

def print_progress_bar(samples_done, samples_total, program_name="", object_name="sample(s)"):
    done = int(samples_done/samples_total*50)
    bar = '#'*done + '-'*(50-done)
    progress_bar = f"{program_name} |{bar}| {samples_done}/{samples_total} {object_name} processed"
    end = "\r"
    if samples_done == samples_total:
        end = "\n"
    print(progress_bar, end=end)
    
def combine_metadata_from_folders(folders, metadata_filename="metadata.txt"):
    if isinstance(folders, str):
        folders = [folders]
    list_of_metadata_dfs = []
    for folder in folders:
        filenames = []
        full_paths = False
        curr_metadata = pd.read_csv(os.path.join(folder, metadata_filename), sep="\t")
        if "#file.name" in curr_metadata.columns:
            if os.path.exists(curr_metadata.iloc[0]["#file.name"]):
                full_paths = True
            if not full_paths:
                curr_metadata["#file.name"] = curr_metadata["#file.name"].apply(lambda x: os.path.join(folder, x))
        list_of_metadata_dfs.append(curr_metadata)
    return pd.concat(list_of_metadata_dfs)

def run_parallel_calculation(function, tasks, program_name, object_name="tasks", verbose=True, cpu=None):
    result_list = []
    tasks_total = len(tasks)
    tasks_done = 0
    if verbose:
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
    list_of_numbers = list(list_of_numbers)
    total_size = sum(list_of_numbers)
    freqs = [s/total_size for s in list_of_numbers]
    diversity = len(list_of_numbers)
    sw = -sum([f*np.log(f) for f in freqs])
    sw_norm = sw/np.log(diversity)
    clonality = 1 - sw_norm
    chao1 = calc_chao1_index(list_of_numbers)

    results = {"shannon_wiener": sw,
               "norm_shannon_wiener": sw_norm,
               "diversity": diversity,
               "clonality": clonality,
               "chao1": chao1
               }

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
    return string[math.ceil(len(string)/2)-3:math.ceil(len(string)/2)+2]

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

def kl_divergence(p, q):
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def jensen_shannon_divergence(p, q):
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    
    # Normalize the distributions
    p /= np.sum(p)
    q /= np.sum(q)
    
    # Calculate the average distribution
    m = 0.5 * (p + q)
    
    # Calculate KLD for each distribution
    kld_p_m = kl_divergence(p, m)
    kld_q_m = kl_divergence(q, m)
    
    # Calculate JSD
    jsd = 0.5 * kld_p_m + 0.5 * kld_q_m
    return jsd