import pandas as pd
import numpy as np
import os
import re
from .io import read_mixcr_clonoset
from .common_functions import print_progress_bar
import random


CHAIN_VARIANTS = {"TRA": {"TRAD", "TRA"},
                  "TRB": {"TRB"},
                  "TRD": {"TRAD", "TRD"},
                  "TRG": {"TRG"},
                  "TRAD": {"TRA", "TRD", "TRAD"},
                  "TCR": {"TRA", "TRB", "TRG", "TRD", "TRAD"},
                  "IGH": {"IGH"},
                  "IGK": {"IGK"},
                  "IGL": {"IGL"},
                  "BCR": {"IGH", "IGK", "IGL"}}

def find_all_exported_clonosets(folders, chain=None):
    if isinstance(folders, str):
        folders = [folders]
    clonosets_dfs = []
    for folder in folders:
        clonosets_dfs.append(find_all_exported_clonosets_in_folder(folder, chain=chain))
    return pd.concat(clonosets_dfs)
    
def find_all_exported_clonosets_in_folder(folder, chain=None):

    if chain is not None:
        if chain.upper() not in CHAIN_VARIANTS:
            chains = ",".join(list(CHAIN_VARIANTS.keys()))
            print(f"ERROR: Unsupported chain name '{chain}'. Possible variants: {chains}")
    
    all_files = os.listdir(folder)
    files = []
    for f in all_files:
        old_match = re.match("(\S+)\.clonotypes\.([A-Z]+)\.txt", f)
        if old_match is not None:
            sample_id = old_match.group(1)
            sample_chain = old_match.group(2)
            # print(f"{f}: old_match")
        else:
            new_match = re.match("(\S+)\.clones_([A-Z]+)\.tsv", f)
            if new_match is not None:
                sample_id = new_match.group(1)
                sample_chain = new_match.group(2)
                # print(f"{f}: new_match")
            else:
                # print(f"{f}: no_match")
                continue
                
        if chain is None:
            if sample_chain in CHAIN_VARIANTS:
                files.append([sample_id, sample_chain, os.path.join(folder, f)])
            else:
                continue
        else:
            if sample_chain not in CHAIN_VARIANTS[chain.upper()]:
                continue
            files.append([sample_id, sample_chain, os.path.join(folder, f)])
        
    files_df = pd.DataFrame(files, columns=["sample_id", "chain", "filename"])
    return files_df

def get_clonoset_stats_from_df(clonosets_df):
    results = []
    columns = ["sample_id", "clones", "clones_func","clones_func_singletons", "clones_func_non_singletons", "reads", "reads_func", "umi", "umi_func"]
    
    samples_total = len(clonosets_df)
    samples_done = 0
    program_name = "Clonoset stats"
    print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    for i,r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        clonoset = read_mixcr_clonoset(filename)
        colnames = get_column_names_from_clonoset(clonoset)
        count_column = colnames["count_column"]
        umi_column = colnames["umi_column"]
                
        # all stats
        clones_num = len(clonoset)
        read_num = clonoset[count_column].sum()
        umi_count = None
        if colnames["umi"] is not None:
            umi_count = clonoset[umi_column].sum()
        
        # stats for functional clones
        clonoset = filter_nonfunctional_clones(clonoset, colnames=colnames)
        func_clones_num = len(clonoset)
        func_read_num = clonoset[count_column].sum()
        func_umi_count = None
        if colnames["umi"] is not None:
            func_umi_count = clonoset[umi_column].sum()
            func_singletons = len(clonoset.loc[clonoset[umi_column] == 1])
        else:
            func_singletons = len(clonoset.loc[clonoset[count_column] == 1])
        results.append([sample_id, clones_num, func_clones_num, func_singletons, func_clones_num-func_singletons, read_num, func_read_num, umi_count, func_umi_count])
        samples_done += 1
        print_progress_bar(samples_done, samples_total, program_name=program_name)
        
    return pd.DataFrame(results, columns=columns)
        
def filter_clonosets_by_sample_list(clonosets_df, samples_list):
    if samples_list is not None:
        clonosets_df = clonosets_df.loc[clonosets_df["sample_id"].isin(samples_list)]
    return clonosets_df

def filter_nonfunctional_clones(clonoset, colnames=None):
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    return clonoset.loc[~clonoset[colnames["cdr3aa_column"]].str.contains("\*|_")]


def get_clonoset_stats(folders, samples_list=None, chain=None):
    clonosets_df = find_all_exported_clonosets(folders, chain=chain)
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    # print(clonosets_df)
    stats = get_clonoset_stats_from_df(clonosets_df)
    return stats


def get_column_names_from_clonoset(clonoset):

    colnames = {"count_column": None,
                "fraction_column": None,
                "umi": None,
                "umi_column": None,
                "umi_fraction_column": None,
                "cdr3aa_column": None,
                "cdr3nt_column": None,
                "v_column": None,
                "d_column": None,
                "j_column": None,
                "c_column": None
                }

    if "uniqueUMICount" in clonoset.columns:
        colnames["umi_column"] = "uniqueUMICount"
        colnames["umi_fraction_column"] = "uniqueUMIFraction"
    if "uniqueMoleculeCount" in clonoset.columns:
        colnames["umi_column"] = "uniqueMoleculeCount"
        colnames["umi_fraction_column"] = "uniqueMoleculeFraction"
    if colnames["umi_column"] is not None:
        colnames["umi"] = True

    if "cloneCount" in clonoset.columns:
        colnames["count_column"] = "cloneCount"
        colnames["fraction_column"] = "cloneFraction"
    if "count" in clonoset.columns:
        colnames["count_column"] = "count"
        colnames["fraction_column"] = "freq"
    if "readCount" in clonoset.columns:
        colnames["count_column"] = "readCount"
        colnames["fraction_column"] = "readFraction"

    if "aaSeqCDR3" in clonoset.columns:
        colnames["cdr3aa_column"] = "aaSeqCDR3"
    if "cdr3aa" in clonoset.columns:
        colnames["cdr3aa_column"] = "cdr3aa"

    if "nSeqCDR3" in clonoset.columns:
        colnames["cdr3nt_column"] = "nSeqCDR3"
    if "cdr3nt" in clonoset.columns:
        colnames["cdr3nt_column"] = "cdr3nt"

    if "allVHitsWithScore" in clonoset.columns:
        colnames["v_column"] = "allVHitsWithScore"
    if "v" in clonoset.columns:
        colnames["v_column"] = "v"

    if "allDHitsWithScore" in clonoset.columns:
        colnames["d_column"] = "allDHitsWithScore"
    if "d" in clonoset.columns:
        colnames["d_column"] = "d"

    if "alljHitsWithScore" in clonoset.columns:
        colnames["j_column"] = "allJHitsWithScore"
    if "j" in clonoset.columns:
        colnames["j_column"] = "j"

    if "allCHitsWithScore" in clonoset.columns:
        colnames["c_column"] = "allCHitsWithScore"
    if "c" in clonoset.columns:
        colnames["c_column"] = "c"

    return colnames


def downsample_clonoset(clonoset, downsample_size, seed=None, by_umi=False):
    clonoset= clonoset.copy()
    if seed is not None:
        random.seed(seed)
    
    colnames = get_column_names_from_clonoset(clonoset)
    count_column = colnames["count_column"]
    if by_umi:
        if colnames["umi"] is not None:
            count_column = colnames["umi_column"]
        else:
            print("WARNING! This clonoset does not contain UMI column. Using reads for downsample instead.")
            
    total_count = int(clonoset[count_column].sum())
    
    if total_count < downsample_size:
        return f"total count {total_count} is less than downsample size {downsample_size}"
    
    sample = sorted(random.sample(range(total_count), downsample_size))
    curr_sum = 0
    i = 0
    new_counts_dict = {}
    for index,r in clonoset.iterrows():
        curr_sum+=r[count_column]
        new_count = 0
        if i == downsample_size:
            break
        while(sample[i]<curr_sum):
            new_count+=1
            i+=1
            if i == downsample_size:
                break
        if new_count > 0:
            new_counts_dict[index]=new_count
    (indices,counts) = zip(*new_counts_dict.items())
    clonoset = clonoset.loc[clonoset.index.isin(indices)]
    clonoset[count_column] = counts
    
    clonoset = recount_fractions_for_clonoset(clonoset, colnames=colnames)
    
    return clonoset.reset_index(drop=True)


def recount_fractions_for_clonoset(clonoset, colnames=None):
    clonoset = clonoset.copy()

    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    count_column = colnames["count_column"]
    fraction_column = colnames["fraction_column"]
    umi_column = colnames["umi_column"]
    umi_fraction_column = colnames["umi_fraction_column"]

    clonoset[fraction_column] = clonoset[count_column]/clonoset[count_column].sum()
    if colnames["umi"]:
        clonoset[umi_fraction_column] = clonoset[umi_column]/clonoset[umi_column].sum()
    return clonoset


##################### unchecked functions ########################

# dirty
def diversity_estimation(folders, downsample_size, chain=None, seed=None, only_functional=True, samples_list=None, iterations=3, by_umi=False):
    
    clonosets_df = find_all_exported_clonosets(folders, chain=chain)
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    
    samples_total = len(clonosets_df)
    samples_done = 0
    
    program_name = "Diversity estimate"
    
    cf.print_progress_bar(samples_done, samples_total, program_name=program_name)
    results = []
    for i, r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        clonoset = read_mixcr_clonoset(filename).rename(columns={'uniqueUMICount': 'uniqueMoleculeCount',
                                                                 'uniqueUMIFraction': 'uniqueMoleculeFraction'})
        if only_functional:
            clonoset = filter_nonfunctional_clones(clonoset)
        
        diversity_values = []
        sw_values = []
        n_sw_values = []

        for i in range(iterations):
            downsampled = downsample_clonoset(clonoset, downsample_size)
            if isinstance(clonoset, str):
                print(f"Error in clonoset '{filename}': {clonoset}")
                return
            if 'cloneCount' in downsampled.columns:
                column = 'cloneCount'
            elif by_umi and 'uniqueMoleculeCount' in downsampled.columns:
                column = "uniqueMoleculeCount"
            else:
                column = "readCount"
            sw, sw_norm, diversity = shannon_wiener(downsampled[column])
            diversity_values.append(diversity)
            sw_values.append(sw)
            n_sw_values.append(sw_norm)
        results.append([sample_id, np.mean(diversity_values), np.mean(sw_values), np.mean(n_sw_values)])
        samples_done += 1
        cf.print_progress_bar(samples_done, samples_total, program_name=program_name)
    return pd.DataFrame(results, columns = ["sample_id", "observed_diversity", "shannon_wiener", "norm_shannon_wiener"])


# dirty
def diversity_estimation_parallel(folders, downsample_size, chain=None, seed=None, only_functional=True, samples_list=None, iterations=3, by_umi=False):
    
    clonosets_df = find_all_exported_clonosets(folders, chain=chain)
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    
    program_name = "Diversity estimate"
    
    tasks=[]
    for i, r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        task=(sample_id, filename, iterations, only_functional, downsample_size, by_umi)
        tasks.append(task)
        
    result_list = run_parallel_calculation(diversity_estimation_mp, tasks, program_name)
    df = pd.DataFrame(results, columns = ["sample_id", "observed_diversity", "shannon_wiener", "norm_shannon_wiener"])
    return df

# dirty
def diversity_estimation_mp(args):
    (sample_id, filename, iterations, only_functional, downsample_size, by_umi) = args
    clonoset = read_mixcr_clonoset(filename).rename(columns={'uniqueUMICount': 'uniqueMoleculeCount',
                                                                 'uniqueUMIFraction': 'uniqueMoleculeFraction'})
    
    if only_functional:
        clonoset = filter_nonfunctional_clones(clonoset)
        
    diversity_values = []
    sw_values = []
    n_sw_values = []

    for i in range(iterations):
        downsampled = downsample_clonoset(clonoset, downsample_size)
        if isinstance(clonoset, str):
            print(f"Error in clonoset '{filename}': {clonoset}")
            return
        if 'cloneCount' in downsampled.columns:
            column = 'cloneCount'
        elif by_umi and 'uniqueMoleculeCount' in downsampled.columns:
            column = "uniqueMoleculeCount"
        else:
            column = "readCount"
        sw, sw_norm, diversity = shannon_wiener(downsampled[column])
        diversity_values.append(diversity)
        sw_values.append(sw)
        n_sw_values.append(sw_norm)
    return (sample_id, np.mean(diversity_values), np.mean(sw_values), np.mean(n_sw_values))


# dirty
def take_top_clonotypes(folders, output_folder, top=0, chain=None, samples_list=None, only_functional=True, mix_tails=True, count_by_umi=False):
    # create output_folder if it doesnot exist
    os.makedirs(output_folder, exist_ok=True)
    
    clonosets_df = find_all_exported_clonosets(folders, chain=chain)
    if samples_list is not None:
        clonosets_df = clonosets_df.loc[clonosets_df["sample_id"].isin(samples_list)].reset_index(drop=True)
    
    samples_total = len(clonosets_df)
    samples_done = 0
    program_name = "Take top clonotypes"
    cf.print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    for i, r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        new_filename = os.path.join(output_folder, os.path.basename(filename))
        clonoset = pd.read_csv(filename, sep="\t")
        clonoset = take_top_clonotypes_in_clonoset(clonoset, top=top, only_functional=only_functional, mix_tails=mix_tails, count_by_umi=count_by_umi)
        clonoset.to_csv(new_filename, index=False, sep="\t")
        samples_done += 1
        cf.print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    print(f"Saved {samples_total} sample(s) to: {output_folder}")

# dirty
def take_top_clonotypes_in_clonoset(clonoset, top=0, only_functional=True, mix_tails=True):
    
    count_column = "cloneCount"
    
    if only_functional:
        clonoset=filter_nonfunctional_clones(clonoset)
        
    if mix_tails:
        clonoset=clonoset.sample(frac=1, random_state=1) #shuffle
        
    clonoset=clonoset.sort_values(by=count_column, ascending=False)
    
    if top > 0:
        clonoset=clonoset.iloc[:top]
    clonoset = recount_fractions_for_clonoset(clonoset)
    
    return clonoset

    
# dirty
def pool_clonosets(folders, samples_list=None, only_functional=False):
    clonosets_df = find_all_exported_clonosets(folders)
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    samples = len(clonosets_df)
    clonosets = []
    for i,r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        clonoset = pd.read_csv(filename, sep="\t")
        if only_functional:
            clonoset = filter_nonfunctional_clones(clonoset)
            clonoset = recount_fractions_for_clonoset(clonoset)
        clonoset["sample_id"] = sample_id
        clonosets.append(clonoset)
    pooled_clonoset = pd.concat(clonosets).reset_index(drop=True)
    total_clones = len(pooled_clonoset)
    print(f"Pooled {total_clones} clones from {samples} samples")
    return pooled_clonoset

# dirty
def convert_mixcr_clonoset(filename, output_filename, filter_nonfunctional=False, by_umi=True):
    m_clonoset = read_mixcr_clonoset(filename)
    count_column, freq_column = None, None
    if "cloneCount" in m_clonoset.columns and "cloneFraction" in m_clonoset.columns:
        count_column = "cloneCount"
        freq_column = "cloneFraction"
    elif "readCount" in m_clonoset.columns and "readFraction" in m_clonoset.columns:
        count_column = "readCount"
        freq_column = "readFraction"
    if by_umi:
        if "uniqueMoleculeCount" in m_clonoset.columns and "uniqueMoleculeFraction" in m_clonoset.columns:
            count_column = "uniqueMoleculeCount"
            freq_column = "uniqueMoleculeFraction"
        elif "uniqueUMICount" in m_clonoset.columns and "uniqueUMIFraction" in m_clonoset.columns:
            count_column = "uniqueUMICount"
            freq_column = "uniqueUMIFraction"
    if count_column is None:
        raise KeyError(f"No count column in file {filename}")
    
    
    
    m_clonoset = m_clonoset.rename(columns={count_column:"count",
                                            freq_column:"freq",
                                            "nSeqCDR3": "cdr3nt",
                                            "aaSeqCDR3": "cdr3aa"})
    m_clonoset["v"] = m_clonoset["allVHitsWithScore"].apply(lambda x: extract_segment(x))
    m_clonoset["d"] = m_clonoset["allDHitsWithScore"].apply(lambda x: extract_segment(x))
    m_clonoset["j"] = m_clonoset["allJHitsWithScore"].apply(lambda x: extract_segment(x))
    m_clonoset["VEnd"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
    m_clonoset["DStart"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
    m_clonoset["DEnd"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
    m_clonoset["JStart"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
    v_clonoset = m_clonoset.sort_values(by="count", ascending=False).reset_index(drop=True)[["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart"]]
    v_clonoset.to_csv(output_filename, index=False, sep = "\t")
    return

# dirty
def convert_mixcr_to_vdjtools(file_list, output_folder, filter_nonfunctional=False, by_umi=True):
    metadata_list = []
    for f in file_list:
        basename = os.path.splitext(os.path.basename(f))[0]
        sample_id = basename.split(".")[0]
        new_filename = f"vdjtools.{sample_id}.txt"
        new_path = os.path.join(output_folder, new_filename)
        convert_mixcr_clonoset(f, new_path, filter_nonfunctional=filter_nonfunctional, by_umi=by_umi)
        metadata_list.append([new_filename, sample_id])
    metadata_filename = os.path.join(output_folder, "metadata.txt")
    metadata = pd.DataFrame(metadata_list, columns=["#file.name", "sample.id"])
    metadata.to_csv(metadata_filename, index=False, sep="\t")