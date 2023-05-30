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

def run_parallel_calculation(function, tasks, program_name, object_name="tasks"):
    result_list = []
    tasks_total = len(tasks)
    tasks_done = 0
    print_progress_bar(tasks_done, tasks_total, program_name, object_name=object_name)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in executor.map(function, tasks):
            result_list.append(result)
            tasks_done+=1
            print_progress_bar(tasks_done, tasks_total, program_name, object_name=object_name)
    return result_list

def shannon_wiener(list_of_numbers):
    list_of_numbers = list(list_of_numbers)
    total_size = sum(list_of_numbers)
    freqs = [s/total_size for s in list_of_numbers]
    diversity = len(list_of_numbers)
    sw = -sum([f*np.log(f) for f in freqs])
    sw_norm = sw/np.log(diversity)
    return sw, sw_norm, diversity

def extract_segment(s):
    segm = str(s).split("*")[0]
    segm = str(s).split("(")[0]
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