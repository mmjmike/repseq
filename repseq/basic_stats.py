import numpy as np
import pandas as pd
import math

from .io import read_mixcr_clonoset
from .common_functions import print_progress_bar, extract_refpoint_position
from .clonosets import (get_column_names_from_clonoset, filter_nonfunctional_clones,
                        recount_fractions_for_clonoset, filter_clonosets_by_sample_list,
                        get_clonoset_stats_from_df)


def calculate_basic_stats(clonosets_df, samples_list=None, by_umi=True, add_5_central=False, top=None, downsample=None, iterations=3, only_functional=True, seed=None):
    # filter by sample list
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    
    stats = get_clonoset_stats_from_df(clonosets_df)
    nt_len_insert_size_and_convergence = calc_nt_len_insert_size_and_convergence_for_df(clonosets_df, by_umi=by_umi, only_functional=only_functional)
    
    # diversity
    # aa_properties
    # aa_5_central_properties
    # merge
    
    
    ### ADD CHECK FOR MULTIPLE CHAINS FOR THE SAME SAMPLE_ID
    
    clonosets_df = clonosets_df.merge(stats)
    clonosets_df = clonosets_df.merge(nt_len_insert_size_and_convergence)

    
    return clonosets_df

def calc_nt_len_insert_size_and_convergence_for_df(clonosets_df, by_umi=True, only_functional=True):
    results = []
    columns = ["sample_id", "mean_nt_len", "mean_nt_insert_size","convergence"]
    
    samples_total = len(clonosets_df)
    samples_done = 0
    program_name = "CDR3nt_len, insert_size and Convergence stats"
    print_progress_bar(samples_done, samples_total, program_name=program_name)
    
    for i,r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        clonoset = read_mixcr_clonoset(filename)
        
        colnames = get_column_names_from_clonoset(clonoset)
        cdr3nt_column = colnames["cdr3nt_column"]
        cdr3aa_column = colnames["cdr3aa_column"]
        fraction_column = colnames["fraction_column"]
        umi_fraction_column = colnames["umi_fraction_column"]
        
        if by_umi:
            if colnames["umi"]:
                fraction_column = umi_fraction_column
            else:
                print("WARNING! This clonoset does not contain UMI column. Using reads for downsample instead.")
                
        if only_functional:
            clonoset = filter_nonfunctional_clones(clonoset, colnames=colnames)
            clonoset = recount_fractions_for_clonoset(clonoset, colnames=colnames)
        
        
        if "refPoints" in clonoset.columns:
            clonoset["VEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
            clonoset["DStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
            clonoset["DEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
            clonoset["JStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
        
        clonoset["nt_len"] = clonoset[cdr3nt_column].apply(lambda x: len(x))
        clonoset["insert_size"] = clonoset.apply(lambda x: calc_insert_size(x.VEnd, x.DStart, x.DEnd, x.JStart), axis=1)
        
        nt_len_mean = np.average(clonoset["nt_len"], weights=clonoset[fraction_column])
        insert_size_mean = np.average(clonoset["insert_size"], weights=clonoset[fraction_column])
        clones = len(clonoset)
        unique_aa_cdr3 = len(clonoset[cdr3aa_column].unique())
        convergence = round(clones/unique_aa_cdr3, 8)
        
        results.append([sample_id, nt_len_mean, insert_size_mean, convergence])
        samples_done += 1
        print_progress_bar(samples_done, samples_total, program_name=program_name)
        
    return pd.DataFrame(results, columns=columns)

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

def center_5(string):
    return string[math.ceil(len(string)/2)-3:math.ceil(len(string)/2)+2]

def center_52(string):
    return string[int(len(string)/2)-3:int(len(string)/2)+2]

