from .common_functions import print_progress_bar, combine_metadata_from_folders
from .processing_stats import pool_metadata
import pandas as pd
import re
import os

def calculate_segment_usage_batch(list_of_folders_after_mixcr, segment="V", functional=True):
    if isinstance(list_of_folders_after_mixcr, str):
        list_of_folders_after_mixcr = [list_of_folders_after_mixcr]
    
    segment_variants=["V", "J"]
    
    if segment.upper() not in segment_variants:
        raise ValueError("ERROR: unsupported segment variant '{}'. Possible values: {}".format(segment, ", ".join(segment_variants)))
    list_of_usages = []
    for folder in list_of_folders_after_mixcr:
        print("\nCalculating {}-usage for samples from: {}\n--------------------------".format(segment.upper(), folder))
        usage = get_segment_usage_from_metafile(folder, segment=segment, functional=functional)
        usage_file = os.path.join(folder, "all_samples_{}_usage.tsv".format(segment.lower()))
        usage.to_csv(usage_file, sep="\t", index=False)
        list_of_usages.append(usage)
        print("Saved {}-usage to file: {}".format(segment.upper(), usage_file))
    return list_of_usages
              
def get_segment_usage_from_metafile(folder, segment="V", functional=True):
    list_of_all_samples_usages = []
    
    if isinstance(folder, str):
        folder = [folder]
    metadata_filename = "metadata.txt"
    metadata = combine_metadata_from_folders(folder, metadata_filename=metadata_filename)
            
    samples_num = len(metadata)
    samples_finished = 0
    first = True
    
    segment_column_name = "all{}HitsWithScore".format(segment.upper())
    
    for index, row in metadata.iterrows():
        sample_id=row["sample.id"]
        file_name=row["#file.name"]
        
        clonoset_df = pd.read_csv(file_name, sep="\t")
        if functional:
            clonoset_df = clonoset_df.loc[~clonoset_df["aaSeqCDR3"].str.contains("\*|_")]
        else:
            clonoset_df = clonoset_df.loc[clonoset_df["aaSeqCDR3"].str.contains("\*|_")]
        
        clonoset_df = clonoset_df[clonoset_df[segment_column_name].notnull()]
        #sort by counts
        clonoset_df.sort_values(by="cloneCount", ascending=False, inplace=True) 
        
        #adjust frequencies in subset table (make them sum to 1)
        clonoset_df.loc[:,"cloneFraction"]=clonoset_df["cloneCount"]/clonoset_df["cloneCount"].sum()
        
        #take only the first segment
        clonoset_df[segment_column_name]=clonoset_df[segment_column_name].apply(lambda x: x.split("*")[0])
        
        sample_usage=clonoset_df.groupby([segment_column_name]).sum()["cloneFraction"]
        
        usage_df = pd.DataFrame({segment.lower():sample_usage.index, 'usage':sample_usage.values})
        usage_df["sample_id"] = sample_id
        list_of_all_samples_usages.append(usage_df)      

        samples_finished += 1
        print_progress_bar(samples_finished, samples_num, program_name="Calculate {}-usage".format(segment.upper())) 
    
    all_usage_df = pd.concat(list_of_all_samples_usages)
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
    df["short_segm"] = df[segment.lower()].apply(lambda x: x[0:3])
    return df.groupby(["short_segm"]).sum().sort_values("usage", ascending=False).index[0]

def determine_chain_for_samples(df):
    segment = df.columns[0]
    sample_chains = []
    for sample in list(df.sample_id.unique()):
        chain = get_chain(df.loc[df.sample_id == sample], segment)
        sample_chains.append((sample, chain))
    return pd.DataFrame(sample_chains, columns=["sample_id", "chain"])
