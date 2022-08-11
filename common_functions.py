import os
import pandas as pd

def print_progress_bar(samples_done, samples_total, program_name="", object_name="sample(s)"):
    done = int(samples_done/samples_total*50)
    bar = '#'*done + '-'*(50-done)
    progress_bar = f"{program_name} |{bar}| {samples_done}/{samples_total} {object_name} processed"
    end = "\r"
    if samples_done == samples_total:
        end = "\n"
    print(progress_bar, end=end)
    
def combine_metadata_from_folders(folders, metadata_filename="metadata.txt"):
    list_of_metadata_dfs = []
    for folder in folders:
        filenames = []
        full_paths = False
        curr_metadata = pd.read_csv(os.path.join(folder, metadata_filename), sep="\t")
        if os.path.exists(curr_metadata.iloc[0]["#file.name"]):
            full_paths = True
        if not full_paths:
            curr_metadata["#file.name"] = curr_metadata["#file.name"].apply(lambda x: os.path.join(folder, x))
        list_of_metadata_dfs.append(curr_metadata)
    return pd.concat(list_of_metadata_dfs)