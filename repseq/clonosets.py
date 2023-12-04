import pandas as pd
import numpy as np
import os
import re
from .io import read_mixcr_clonoset, read_clonoset
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

def find_all_exported_clonosets(folders, chain=None, remove_non_target=False, non_target_threshold=0.01):
    """
    Main method of the Filter object - application of it to a clonoset

    Args:
        folders (path or list of paths): clonoset in the form of Pandas DataFrame in
            MiXCR(3 or 4+ version), VDJtools or Bioadaptive formats.
        colnames (dict, optional): Dictionary of available specific column names.
            Defaults to None - colnames imputed automatically.

    Returns:
        clonoset (pd.DataFrame): clonoset after converting to common (VDJtools-like)
            format and applying functionality filtration and downsampling or taking top

    """
    if isinstance(folders, str):
        folders = [folders]
    clonosets_dfs = []
    for folder in folders:
        clonosets_dfs.append(find_all_exported_clonosets_in_folder(folder, chain=chain, remove_non_target=remove_non_target, non_target_threshold=non_target_threshold))
    return pd.concat(clonosets_dfs)
   
def find_all_exported_clonosets_in_folder(folder, chain=None, remove_non_target=False, non_target_threshold=0.01):
    result_columns = ["sample_id", "chain", "filename"]
    if chain is not None:
        if chain.upper() not in CHAIN_VARIANTS:
            chains = ",".join(list(CHAIN_VARIANTS.keys()))
            print(f"ERROR: Unsupported chain name '{chain}'. Possible variants: {chains}")
    
    all_files = os.listdir(folder)
    files = []
    for f in all_files:
        old_match = re.match("(\S+)\.clonotypes\.([A-Z]+)(?:\W\S+)*\.txt", f)
        if old_match is not None:
            sample_id = old_match.group(1)
            sample_chain = old_match.group(2)
            # print(f"{f}: old_match")
        else:
            new_match = re.match("(\S+)\.clones_([A-Z]+)(?:\W\S+)*\.tsv", f)
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
        
    files_df = pd.DataFrame(files, columns=result_columns)
    if remove_non_target:
        files_df = remove_non_target_clonosets_from_files_df(files_df, threshold=non_target_threshold)

    return files_df

def remove_non_target_clonosets_from_files_df(files_df, threshold=0.01):
    columns = files_df.columns
    sample_counts = files_df.sample_id.value_counts()
    files_df["read_count"] = 1
    files_df["total_read_count"] = 1
    for sample,count in sample_counts.items():
        if count > 1:
            sample_reads = 0
            for i, r in files_df.loc[files_df.sample_id == sample].iterrows():
                clonoset = read_mixcr_clonoset(r["filename"])
                colnames = get_column_names_from_clonoset(clonoset)
                clonoset_reads = clonoset[colnames["count_column"]].sum()
                files_df.loc[i,"read_count"] = clonoset_reads
                sample_reads += clonoset_reads
            files_df.loc[files_df.sample_id == sample, "total_read_count"] = sample_reads
    files_df = files_df.loc[files_df.read_count/files_df.total_read_count > threshold]
    return files_df[columns].reset_index(drop=True)

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


def filter_nonfunctional_clones(clonoset_in, colnames=None,):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    clonoset = clonoset.loc[~clonoset[colnames["cdr3aa_column"]].str.contains("\*|_")]
    # clonoset = clonoset.loc[clonoset[colnames["cdr3aa_column"]] != ""]
    return clonoset

def filter_by_functionality(clonoset_in, colnames=None, functional=True):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    cdr3aa_column = colnames["cdr3aa_column"]
    if functional:
        clonoset = clonoset.loc[~clonoset[cdr3aa_column].str.contains("\*|_")]
        clonoset = clonoset.loc[clonoset[cdr3aa_column] != ""]
    else:
        clonoset = clonoset.loc[(clonoset[cdr3aa_column].str.contains("\*|_")) | (clonoset[cdr3aa_column] == "")]

    return clonoset



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
    if "Read.count" in clonoset.columns:
        colnames["count_column"] = "Read.count"

    if "aaSeqCDR3" in clonoset.columns:
        colnames["cdr3aa_column"] = "aaSeqCDR3"
    if "cdr3aa" in clonoset.columns:
        colnames["cdr3aa_column"] = "cdr3aa"
    if "CDR3.amino.acid.sequence" in clonoset.columns:
        colnames["cdr3aa_column"] = "cdr3aCDR3.amino.acid.sequencea"

    if "nSeqCDR3" in clonoset.columns:
        colnames["cdr3nt_column"] = "nSeqCDR3"
    if "cdr3nt" in clonoset.columns:
        colnames["cdr3nt_column"] = "cdr3nt"
    if "CDR3.nucleotide.sequence" in clonoset.columns:
        colnames["cdr3nt_column"] = "CDR3.nucleotide.sequence"
    

    if "allVHitsWithScore" in clonoset.columns:
        colnames["v_column"] = "allVHitsWithScore"
    if "v" in clonoset.columns:
        colnames["v_column"] = "v"
    if "bestVGene" in clonoset.columns:
        colnames["v_column"] = "bestVGene"

    if "allDHitsWithScore" in clonoset.columns:
        colnames["d_column"] = "allDHitsWithScore"
    if "d" in clonoset.columns:
        colnames["d_column"] = "d"

    if "allJHitsWithScore" in clonoset.columns:
        colnames["j_column"] = "allJHitsWithScore"
    if "j" in clonoset.columns:
        colnames["j_column"] = "j"
    if "bestJGene" in clonoset.columns:
        colnames["j_column"] = "bestJGene"

    if "allCHitsWithScore" in clonoset.columns:
        colnames["c_column"] = "allCHitsWithScore"
    if "c" in clonoset.columns:
        colnames["c_column"] = "c"

    return colnames



def downsample_clonoset(clonoset, downsample_size, seed=None, by_umi=False, colnames=None):
    clonoset= clonoset.copy()
    if seed is not None:
        random.seed(seed)
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    count_column, fraction_column = decide_count_and_frac_columns(colnames, by_umi)
            
    total_count = int(clonoset[count_column].sum())
    
    if total_count < downsample_size:
        return f"total count {total_count} is less than downsample size {downsample_size}"
    elif total_count == downsample_size:
        clonoset = recount_fractions_for_clonoset(clonoset, colnames=colnames)
        return clonoset

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

def take_top_clonotypes_in_clonoset(clonoset, top, only_functional=True, mix_tails=True, seed=None, by_umi=False, colnames=None):
    clonoset= clonoset.copy()
    if only_functional:
        clonoset=filter_nonfunctional_clones(clonoset)
        
    if mix_tails:
        clonoset=clonoset.sample(frac=1, random_state=seed) #shuffle
    
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    
    count_column, fraction_column = decide_count_and_frac_columns(colnames, by_umi)

    clonoset=clonoset.sort_values(by=count_column, ascending=False)
    
    if top > len(clonoset):
        print(f"Warning! Clonoset size - {len(clonoset)} - is less than required top - {top}")
    if top > 0:
        clonoset=clonoset.iloc[:top]
    clonoset = recount_fractions_for_clonoset(clonoset)
    
    return clonoset

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




def pool_clonotypes_from_clonosets_df_old(clonosets_df, samples_list=None, top=None, downsample=None, only_functional=True, by_umi=False, exclude_singletons=False, seed=None, rename_columns_for_clustering=False):
    
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    
    clonotypes_dfs = []
    
    for index, row in clonosets_df.iterrows():
        sample_id = row["sample_id"]
        filename = row["filename"]
        clonoset=read_mixcr_clonoset(filename)
        colnames = get_column_names_from_clonoset(clonoset)
        if only_functional:
            clonoset = filter_nonfunctional_clones(clonoset, colnames=colnames)
        
        count_column, fraction_column = decide_count_and_frac_columns(colnames, by_umi)
        if exclude_singletons:
            clonoset = clonoset.loc[clonoset[count_column] > 1]
        
        clonoset_size = clonoset[count_column].sum()
        if downsample is not None:
            if downsample > clonoset_size:
                raise(ValueError, f"Downsample size ({downsample}) exceeds the clonoset size ({clonoset_size})")
            clonoset = downsample_clonoset(clonoset, downsample, seed=seed, by_umi=by_umi)
        clonoset = recount_fractions_for_clonoset(clonoset, colnames=colnames)

        if rename_columns_for_clustering:
            clonoset = clonoset.rename(columns={colnames["v_column"]: "v",
                                    colnames["j_column"]: "j",
                                    colnames["cdr3aa_column"]: "cdr3aa",
                                    colnames["cdr3nt_column"]: "cdr3nt",
                                    colnames["fraction_column"]:"freq",
                                    colnames["count_column"]: "count"})

        clonoset["sample_id"] = sample_id
        clonotypes_dfs.append(clonoset)
    result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)
    clonotypes_number = len(result_df)
    samples_number = len(result_df["sample_id"].unique())
    print(f"Pooled {clonotypes_number} clonotypes from {samples_number} samples")
    return result_df

def decide_count_and_frac_columns(colnames, by_umi, suppress_warnings=False):
    count_column = colnames["count_column"]
    fraction_column = colnames["fraction_column"]
    if by_umi:
        if colnames["umi"] is not None:
            count_column = colnames["umi_column"]
            fraction_column = colnames["umi_fraction_column"]
        elif not suppress_warnings:
            print("WARNING! Clonoset does not contain UMI column. Using reads for clone count instead.\nTo avoid this warning set parameter 'by_umi=False'")
    return count_column, fraction_column

##################### unchecked functions ########################



    
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
# def convert_mixcr_clonoset(filename, output_filename, filter_nonfunctional=False, by_umi=True):
#     m_clonoset = read_mixcr_clonoset(filename)
#     count_column, freq_column = None, None
#     if "cloneCount" in m_clonoset.columns and "cloneFraction" in m_clonoset.columns:
#         count_column = "cloneCount"
#         freq_column = "cloneFraction"
#     elif "readCount" in m_clonoset.columns and "readFraction" in m_clonoset.columns:
#         count_column = "readCount"
#         freq_column = "readFraction"
#     if by_umi:
#         if "uniqueMoleculeCount" in m_clonoset.columns and "uniqueMoleculeFraction" in m_clonoset.columns:
#             count_column = "uniqueMoleculeCount"
#             freq_column = "uniqueMoleculeFraction"
#         elif "uniqueUMICount" in m_clonoset.columns and "uniqueUMIFraction" in m_clonoset.columns:
#             count_column = "uniqueUMICount"
#             freq_column = "uniqueUMIFraction"
#     if count_column is None:
#         raise KeyError(f"No count column in file {filename}")
    
    
#     m_clonoset = m_clonoset.rename(columns={count_column:"count",
#                                             freq_column:"freq",
#                                             "nSeqCDR3": "cdr3nt",
#                                             "aaSeqCDR3": "cdr3aa"})
#     m_clonoset["v"] = m_clonoset["allVHitsWithScore"].apply(lambda x: extract_segment(x))
#     m_clonoset["d"] = m_clonoset["allDHitsWithScore"].apply(lambda x: extract_segment(x))
#     m_clonoset["j"] = m_clonoset["allJHitsWithScore"].apply(lambda x: extract_segment(x))
#     m_clonoset["VEnd"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
#     m_clonoset["DStart"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
#     m_clonoset["DEnd"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
#     m_clonoset["JStart"] = m_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
#     v_clonoset = m_clonoset.sort_values(by="count", ascending=False).reset_index(drop=True)[["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart"]]
#     v_clonoset.to_csv(output_filename, index=False, sep = "\t")
#     return

# dirty
# def convert_mixcr_to_vdjtools(file_list, output_folder, filter_nonfunctional=False, by_umi=True):
#     metadata_list = []
#     for f in file_list:
#         basename = os.path.splitext(os.path.basename(f))[0]
#         sample_id = basename.split(".")[0]
#         new_filename = f"vdjtools.{sample_id}.txt"
#         new_path = os.path.join(output_folder, new_filename)
#         convert_mixcr_clonoset(f, new_path, filter_nonfunctional=filter_nonfunctional, by_umi=by_umi)
#         metadata_list.append([new_filename, sample_id])
#     metadata_filename = os.path.join(output_folder, "metadata.txt")
#     metadata = pd.DataFrame(metadata_list, columns=["#file.name", "sample.id"])
#     metadata.to_csv(metadata_filename, index=False, sep="\t")