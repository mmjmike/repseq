import pandas as pd
import os
from .migec import checkout_stats

def check_number_of_clonotypes(folders, samples_list=[], metadata_filename="vdjtools_metadata.txt", add_total=True):
    all_metadata = pool_metadata(folders, metadata_filename, samples_list)
    
    sample_ids=[]
    number_of_functional_clonotypes=[]
    number_of_clonotypes_total=[]
    number_of_singletons=[]
    number_of_migs_functional=[]
    
    for index, row in all_metadata.iterrows():
        clonoset_data=pd.read_csv(row["#file.name"],sep="\t")
        sample_ids.append(row["sample.id"])
        number_of_clonotypes_total.append(clonoset_data.shape[0])

        clonoset_data = clonoset_data.rename(columns={"bestVGene": "v",
                                        "bestJGene": "j",
                                        "CDR3.amino.acid.sequence": "cdr3aa",
                                        "CDR3.nucleotide.sequence": "cdr3nt",
                                        "allVHitsWithScore": "v",
                                        "allJHitsWithScore": "j",
                                        "aaSeqCDR3": "cdr3aa",
                                        "nSeqCDR3": "cdr3nt",
                                        "Sample":"sample_id",
                                        "cloneFraction":"freq",
                                        "Read.count": "count",
                                        "cloneCount": "count"})
    
        clonoset_data["v"] = clonoset_data["v"].apply(lambda x: x.split("*")[0])
        clonoset_data["j"] = clonoset_data["j"].apply(lambda x: x.split("*")[0])

        clonoset_data=clonoset_data.loc[~clonoset_data["cdr3aa"].str.contains("\*|_", na=False)]
        number_of_functional_clonotypes.append(clonoset_data.shape[0])
        number_of_singletons.append(clonoset_data.loc[clonoset_data["count"]==1].shape[0])
        number_of_migs_functional.append(clonoset_data["count"].sum())
    clonotypes_num_df = pd.DataFrame({"sample_id":  sample_ids, "functional_clonotypes": number_of_functional_clonotypes,
                                      "functional_migs": number_of_migs_functional,
                                      "total_clonotypes": number_of_clonotypes_total, "singletons":number_of_singletons})
    clonotypes_num_df.sort_values(by="functional_clonotypes", ascending=True, inplace=True)
    clonotypes_num_df.reset_index(drop=True, inplace=True)
    clonotypes_num_df["functional_not_singletons"] = clonotypes_num_df["functional_clonotypes"] - clonotypes_num_df["singletons"]
    
    old_samples = False
    if "old.sample.id" in all_metadata.columns:
        old_samples = True
        clonotypes_num_df = clonotypes_num_df.merge(all_metadata.rename(columns={"sample.id": "sample_id",
                         "old.sample.id": "old_sample_id"})[["sample_id", "old_sample_id"]])

    # add TOTAL 
    if add_total:
        total_row = ["TOTAL"] + list(clonotypes_num_df.sum(numeric_only=True, axis=0))
        if old_samples:
            total_row += [""]
        clonotypes_num_df.loc[len(clonotypes_num_df)] = total_row
    return clonotypes_num_df

def pool_metadata(folders, metadata_filename, sample_list=[]):
    # combine metadata from several folders if a list of them is specified
    if isinstance(folders, str):
        folders = [folders]
    all_metadata_dfs = []
    for folder in folders:

        vdjtools_metadata_filename = os.path.join(folder, metadata_filename)
        metadata = pd.read_csv(vdjtools_metadata_filename, sep="\t")

        full_paths = False
        if os.path.exists(metadata.iloc[0]["#file.name"]):
            full_paths = True
        if not full_paths:
            metadata["#file.name"] = metadata["#file.name"].apply(lambda x: os.path.join(folder, x))

        all_metadata_dfs.append(metadata)

    all_metadata = pd.concat(all_metadata_dfs).reset_index(drop=True)
    
    multichain = False
    if "old.sample.id" in all_metadata.columns:
        multichain = True
        all_metadata["old.sample.id"] = all_metadata["old.sample.id"].fillna(all_metadata["sample.id"])

    # filter samples according to sample_list if needed
    if sample_list != []:
        if multichain:
            all_metadata = all_metadata.loc[all_metadata["sample.id"].isin(sample_list) | all_metadata["old.sample.id"].isin(sample_list)]
        else:
            all_metadata = all_metadata.loc[all_metadata["sample.id"].isin(sample_list)]
    columns = ["#file.name", "sample.id"]
    if multichain:
        columns += ["old.sample.id"]
    return all_metadata[columns]


def processing_stats(folders):
    all_assemble_results = pool_metadata(folders, "assemble.log.txt")
    all_assemble_results = all_assemble_results[["#SAMPLE_ID", "MIG_COUNT_THRESHOLD", "READS_TOTAL","READS_GOOD_TOTAL","MIGS_TOTAL","MIGS_GOOD_TOTAL"]]
    all_assemble_results = all_assemble_results.rename(columns={"#SAMPLE_ID":"sample_id",
                                                                "MIG_COUNT_THRESHOLD": "overseq_thr",
                                                                "READS_TOTAL":"reads",
                                                                "READS_GOOD_TOTAL":"reads_good",
                                                                "MIGS_TOTAL":"migs",
                                                                "MIGS_GOOD_TOTAL":"migs_good"})

    return all_assemble_results

def show_all_stats(folder, suffix="", vdjtools_folder=None):
    checkout_folder = os.path.join(folder, "checkout")
    assemble_folder = os.path.join(folder, "assemble")
    if vdjtools_folder is None:
        clonotypes_folder = assemble_folder
    else:
        clonotypes_folder = vdjtools_folder

    if suffix != "":
        checkout_folder += "_" + suffix
        assemble_folder += "_" + suffix
        
    checkout = checkout_stats(checkout_folder)
    processing = processing_stats(assemble_folder)
    clonotypes = check_number_of_clonotypes(clonotypes_folder, add_total=False)
    all_stats = checkout.merge(processing)

    if "old_sample_id" in clonotypes.columns:
        all_stats = clonotypes.merge(all_stats.rename(columns={"sample_id": "old_sample_id"}))
    else:
        all_stats = clonotypes.merge(all_stats)
    all_stats["mean_reads_per_umi"] = all_stats["reads_with_UMI"]/all_stats["migs"]
    return all_stats

def show_important_stats(folder, suffix="", vdjtools_folder=None):
    all_stats = show_all_stats(folder, suffix=suffix, vdjtools_folder=vdjtools_folder)
    return all_stats[["sample_id", "reads_total", "reads_with_UMI", "migs", "migs_good", "total_clonotypes",
                      "functional_clonotypes", "reads_with_umi_percent", "mean_reads_per_umi"]]
