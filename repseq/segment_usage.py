from .common_functions import print_progress_bar, combine_metadata_from_folders
from .processing_stats import pool_metadata
from .clonosets import (find_all_exported_clonosets, filter_clonosets_by_sample_list,
                        filter_nonfunctional_clones, recount_fractions_for_clonoset,
                        get_column_names_from_clonoset)
from .io import read_mixcr_clonoset
import pandas as pd
import re
import os

def calculate_segment_usage_batch(folders, segment="V", only_functional=True, by_umi=False, samples_list=None):
    segment_variants=["V", "J", "C"]
    segment = segment.upper()
    if segment not in segment_variants:
        raise ValueError("ERROR: unsupported segment variant '{}'. Possible values: {}".format(segment, ", ".join(segment_variants)))
    segment_column = f"all{segment}HitsWithScore"
    
    clonosets_df = find_all_exported_clonosets(folders)
    clonosets_df = filter_clonosets_by_sample_list(clonosets_df, samples_list)
    
    print(f"\nCalculating {segment}-usage for {len(clonosets_df)} sample(s)\n"+"-"*50)
    
    samples_finished = 0
    samples_num = len(clonosets_df)
    print_progress_bar(samples_finished, samples_num, program_name="Calculate {}-usage".format(segment.upper()))
    list_of_usages = []
    
    for index, row in clonosets_df.iterrows():
        sample_id=row["sample_id"]
        file_name=row["filename"]
        
        clonoset = read_mixcr_clonoset(file_name)
        
        if only_functional:
            clonoset = filter_nonfunctional_clones(clonoset)
            
        # filter out clones with undetermined segment of interest
        clonoset = clonoset[clonoset[segment_column].notnull()]
        clonoset = recount_fractions_for_clonoset(clonoset)
        
        colnames = get_column_names_from_clonoset(clonoset)
        fraction_column = colnames["fraction_column"]
        umi_fraction_column = colnames["umi_fraction_column"]
        
        if by_umi:
            if colnames["umi"]:
                fraction_column = umi_fraction_column
            else:
                print("WARNING! This clonoset does not contain UMI column. Using reads for downsample instead.")

        #take only the first segment
        clonoset[segment_column]=clonoset[segment_column].apply(lambda x: x.split("*")[0])
        clonoset[segment_column]=clonoset[segment_column].apply(lambda x: x.split("(")[0])
        
        usage_df=clonoset.groupby([segment_column])[[fraction_column]].sum()
        usage_df[segment.lower()] = usage_df.index
        usage_df= usage_df.reset_index(drop=True)
        usage_df["sample_id"] = sample_id
        usage_df = usage_df.rename(columns={fraction_column: "usage"})[[segment.lower(), "usage", "sample_id"]]
        list_of_usages.append(usage_df)      

        samples_finished += 1
        print_progress_bar(samples_finished, samples_num, program_name="Calculate {}-usage".format(segment.upper())) 
    
    all_usage_df = pd.concat(list_of_usages)
    sample_chains = determine_chain_for_samples(all_usage_df)
    all_usage_df = all_usage_df.merge(sample_chains)
    return fill_zero_usage(all_usage_df, segment).reset_index(drop=True)

def fill_zero_usage(df, segment):
    fill_zeros_lines = []
    for chain in list(df.chain.unique()):
        segments = df.loc[df.chain == chain][segment.lower()].unique()
        samples = df.loc[df.chain == chain]["sample_id"].unique()
        pairs_set = set()
        for i, r in df.loc[df.chain == chain].iterrows():
            pairs_set.add((r["sample_id"], r[segment.lower()]))
        for seg in segments:
            for sample in samples:
                if (sample, seg) not in pairs_set:
                    fill_zeros_lines.append((seg, 0, sample, chain))
    return pd.concat([df, pd.DataFrame(fill_zeros_lines, columns=df.columns)])

def get_chain(df, segment):
    df1 = df.copy()
    df1["short_segm"] = df1[segment.lower()].apply(lambda x: x[0:3])
    return df1.groupby(["short_segm"])[["usage"]].sum().sort_values("usage", ascending=False).index[0]

def determine_chain_for_samples(df):
    segment = df.columns[0]
    sample_chains = []
    for sample in list(df.sample_id.unique()):
        chain = get_chain(df.loc[df.sample_id == sample], segment)
        sample_chains.append((sample, chain))
    return pd.DataFrame(sample_chains, columns=["sample_id", "chain"])

