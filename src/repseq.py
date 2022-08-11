import pandas as pd
import os


def check_number_of_clonotypes(folders, samples_list=[], metadata_filename="vdjtools_metadata.txt"):
    all_metadata = pool_metadata(folders, metadata_filename, samples_list)

    sample_ids=[]
    number_of_functional_clonotypes=[]
    number_of_clonotypes_total=[]
    number_of_singletons=[]

    for index, row in all_metadata.iterrows():
        clonoset_data=pd.read_csv(row["#file.name"],sep="\t")
        sample_ids.append(row["sample.id"])
        number_of_clonotypes_total.append(clonoset_data.shape[0])
        clonoset_data=clonoset_data.loc[~clonoset_data["cdr3aa"].str.contains("\*|_", na=False)]
        number_of_functional_clonotypes.append(clonoset_data.shape[0])
        number_of_singletons.append(clonoset_data.loc[clonoset_data["count"]==1].shape[0])
    clonotypes_num_df = pd.DataFrame({"sample_id":  sample_ids, "functional_clonotypes": number_of_functional_clonotypes,
                                      "total_clonotypes": number_of_clonotypes_total, "singletons":number_of_singletons})
    clonotypes_num_df.sort_values(by="functional_clonotypes", ascending=True, inplace=True)
    return clonotypes_num_df

def pool_top_clonotypes_to_df(folders, samples_list=[], top=0, functional=True, exclude_singletons=False, cdr3aa_len_range=[], metadata_filename="vdjtools_metadata.txt"):
    all_metadata = pool_metadata(folders, metadata_filename,samples_list)

    clonotypes_dfs = []
    for index, metadata_line in all_metadata.iterrows():
        clonoset_data=pd.read_csv(metadata_line["#file.name"],sep="\t")
        if exclude_singletons:
            clonoset_data=clonoset_data.loc[~clonoset_data["count"]>1]
        if functional:
            clonoset_data=clonoset_data.loc[~clonoset_data["cdr3aa"].str.contains("\*|_")]
            clonoset_data=clonoset_data.sample(frac=1, random_state=1) #shuffle
            clonoset_data=clonoset_data.sort_values(by="count", ascending=False) #sort by counts "back"
        if top > 0:
            clonoset_data=clonoset_data.iloc[:top]
        if cdr3aa_len_range:
            clonoset_data=clonoset_data.loc[(clonoset_data["cdr3aa"].str.len() <= cdr3aa_len_range[-1])
                                            & (clonoset_data["cdr3aa"].str.len() >= cdr3aa_len_range[0])]
        clonoset_data["freq"]=clonoset_data["count"]/clonoset_data["count"].sum()
        clonoset_data["sample_id"] = sample_id
        clonotypes_dfs.append(clonoset_data)
    return pd.concat(clonotypes_dfs).reset_index(drop=True)

def pool_metadata(folders, metadata_filename, sample_list):
    # combine metadata from several folders if a list of them is specified
    if isinstance(folders, str):
        folders = [folders]
    all_metadata_dfs = []
    for folder in folders:
        vdjtools_metadata_filename = os.path.join(folder, metadata_filename)
        metadata = pd.read_csv(vdjtools_metadata_filename, sep="\t")
        all_metadata_dfs.append(metadata)
    all_metadata = pd.concat(all_metadata_dfs).reset_index(drop=True)

    # filter samples according to sample_list if needed
    if sample_list != []:
        all_metadata = all_metadata.loc[all_metadata["sample.id"].isin(sample_list)]
    return all_metadata
