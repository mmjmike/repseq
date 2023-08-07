import pandas as pd
import numpy as np
import random
import os

from .clonosets import (get_column_names_from_clonoset,
                        filter_by_functionality)
from .common_functions import (calc_insert_size,
                               extract_refpoint_position,
                               center_5,
                               shannon_wiener,
                               run_parallel_calculation)
from .io import read_clonoset
from repseq.clone_filter import Filter

REPSEQ_PATH = os.path.join(os.path.expanduser("~"), "soft", "repseq")
AA_PROPS_PATH = os.path.join(REPSEQ_PATH, "repseq", "resourses", "aa_property_table.txt")


def calc_clonoset_stats(clonosets_df, cl_filter=None):
    df = generic_calculation(clonosets_df, calculate_clonoset_stats_cl, clonoset_filter=cl_filter, program_name="CalcClonosetStats")
    convert_dict = {"clones": int,
                    "clones_func": int,
                    "clones_func_singletons": int,
                    "clones_func_non_singletons": int,
                    "clones_nonfunc": int,
                    "reads": int,
                    "reads_func": int,
                    "reads_nonfunc": int}
    if not df["umi"].isnull().values.any():
        convert_dict.update({"umi": "umi_count",
                        "umi_func": "func_umi_count",
                        "umi_nonfunc": "umi_nonfunc"})
 
    df = df.astype(convert_dict)
    return df


def calc_diversity_stats(clonosets_df, cl_filter=None, iterations=3, seed=None, drop_small_samples=True):
    df = generic_calculation(clonosets_df, calculate_diversity_stats_cl, clonoset_filter=cl_filter,
                             program_name="CalcDiversityStats", iterations=iterations, seed=seed, drop_small_samples=drop_small_samples)
    return df


def calc_convergence(clonosets_df, cl_filter=None, iterations=1, seed=None, drop_small_samples=True):
    df = generic_calculation(clonosets_df, calculate_convergence_cl, clonoset_filter=cl_filter,
                             program_name="CalcConvergence", iterations=iterations, seed=seed, drop_small_samples=drop_small_samples)
    return df


def calc_cdr3_properties(clonosets_df, cl_filter=None, iterations=1, seed=None, drop_small_samples=True):
    if cl_filter is None:
        cl_filter = Filter(functionality="f")
    df = generic_calculation(clonosets_df, calculate_cdr3_properties_cl, clonoset_filter=cl_filter,
                             program_name="CalcCDR3aaProperties", iterations=iterations, seed=seed, drop_small_samples=drop_small_samples)
    return df


def calculate_clonoset_stats_cl(clonoset):
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
    clonoset = filter_by_functionality(clonoset, colnames=colnames)
    func_clones_num = len(clonoset)
    func_read_num = clonoset[count_column].sum()
    func_umi_count = None
    umi_nonfunc = None
    umi_nonfunc_freq = None
    if colnames["umi"] is not None:
        func_umi_count = clonoset[umi_column].sum()
        func_singletons = len(clonoset.loc[clonoset[umi_column] == 1])
        umi_nonfunc = umi_count-func_umi_count
        umi_nonfunc_freq = (umi_count-func_umi_count)/umi_count
        
    else:
        func_singletons = len(clonoset.loc[clonoset[count_column] == 1])

    result = {"clones": clones_num,
              "clones_func": func_clones_num,
              "clones_func_singletons": func_singletons,
              "clones_func_non_singletons": func_clones_num-func_singletons,
              "clones_nonfunc": clones_num-func_clones_num,
              "clones_nonfunc_freq": (clones_num-func_clones_num)/clones_num,
              "reads": read_num,
              "reads_func": func_read_num,
              "reads_nonfunc": read_num-func_read_num,
              "reads_nonfunc_freq": (read_num-func_read_num)/read_num,
              "umi": umi_count,
              "umi_func": func_umi_count,
              "umi_nonfunc": umi_nonfunc,
              "umi_nonfunc_freq": umi_nonfunc_freq}

    return result


def calculate_cdr3_properties_cl(clonoset_in, colnames=None):
    
    # copy input clonoset
    clonoset = clonoset_in.copy()
    # read properties table for amino acids
    aa_properties_dict = pd.read_csv(AA_PROPS_PATH, sep="\t",comment='#').set_index('amino_acid').to_dict()
    list_of_properties = [p for p in aa_properties_dict.keys() if p != "count"]
    
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    
    # calc mean nt_len, mean freq and, if possible, insert size
    cdr3aa_column = colnames["cdr3aa_column"]
    cdr3nt_column = colnames["cdr3nt_column"]
    fraction_column = colnames["fraction_column"]

    clonoset["nt_len"] = clonoset[cdr3nt_column].apply(lambda x: len(x))
    nt_len_mean = np.average(clonoset["nt_len"], weights=clonoset[fraction_column])
    mean_frequency = clonoset[fraction_column].mean()
    
    insert_size_columns = ["VEnd", "DStart", "DEnd", "JStart"]
    insert_size_possible = True
    if "refPoints" in clonoset.columns:
        clonoset["VEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
        clonoset["DStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
        clonoset["DEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
        clonoset["JStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
    for col in insert_size_columns:
        if col not in clonoset.columns:
            insert_size_possible = False
            break
    if insert_size_possible:
        clonoset["insert_size"] = clonoset.apply(lambda x: calc_insert_size(x.VEnd, x.DStart, x.DEnd, x.JStart), axis=1)
        insert_size_mean = np.average(clonoset["insert_size"], weights=clonoset[fraction_column])
    else:
        insert_size_mean = None

    result = {"mean_cdr3nt_len": nt_len_mean,
              "mean_insert_size": insert_size_mean,
              "mean_frequency": mean_frequency}


    # calc cdr3_aa properties both for full cdr3aa seq and for only central 5 aa's in cdr3
    clonoset["center_aa"] = clonoset["cdr3aa"].apply(lambda x: center_5(x))

    for aa_property in list_of_properties:
        prop_5 = "cdr3_5_" + aa_property
        prop_full = "cdr3_full_" + aa_property
        
        weights = clonoset[fraction_column]
        prop_5_values = clonoset["center_aa"].apply(lambda x: sum([aa_properties_dict[aa_property][aa] for aa in x]))
        prop_full_values = clonoset[cdr3aa_column].apply(lambda x: sum([aa_properties_dict[aa_property][aa] for aa in x]))

        result[prop_5] = np.average(prop_5_values, weights=weights)
        result[prop_full] = np.average(prop_full_values, weights=weights)
        
    return result


def calculate_convergence_cl(clonoset_in, colnames=None):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)

    nt_unique = len(clonoset[colnames["cdr3nt_column"]].unique())
    aa_unique = len(clonoset[colnames["cdr3aa_column"]].unique())
    convergernce = round(nt_unique/aa_unique, 8)
    result = {"convergence": convergernce}
    
    return result

def calculate_diversity_stats_cl(clonoset_in, colnames=None):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)

    counts = clonoset[colnames["count_column"]]
    
    sw, sw_norm, diversity = shannon_wiener(counts)
    
    result = {"diversity": diversity,
              "shannon_wiener": sw,
              "normalized_shannon_wiener": sw_norm}
                    
    return result


def generic_calculation(clonosets_df_in, calc_function, clonoset_filter=None, program_name="Calculation", iterations=1, seed=None, drop_small_samples=False):
    columns_retain = ["sample_id"]
    clonosets_df = clonosets_df_in.copy()
    
    if "sample_id" not in clonosets_df.columns:
        raise ValueError("Clonoset_df does not contain required column 'sample_id'")
    if "filename" not in clonosets_df.columns:
        raise ValueError("Clonoset_df does not contain required column 'filename'")
    
    split_chain_after_calculation = False
    if "chain" not in clonosets_df.columns:
        if len(clonosets_df) != len(clonosets_df.sample_id.unique()):
            raise ValueError("Clonoset_df contains nonunique sample_ids")
    else:
        columns_retain.append("chain")
        if len(clonosets_df[["sample_id", "chain"]].drop_duplicates()) != len(clonosets_df):
            raise ValueError("Clonoset_df contains nonunique sample_id+chain combinations")
        if len(clonosets_df) != len(clonosets_df.sample_id.unique()):
            clonosets_df["sample_id"] = clonosets_df["sample_id"] + "_" + clonosets_df["chain"]
            split_chain_after_calculation = False

    random_filter = False
    need_downsample = False
    need_top = False

    count_column_by_umi_and_functionality = {
            True: {"a": "umi",
                   "f": "umi_func",
                   "n": "umi_nonfunc"},
            False: {"a": "reads",
                    "f": "reads_func",
                    "n": "reads_nonfunc"}
        }
    
    clone_column_by_functionality = {
        "a": "clones",
        "f": "clones_func",
        "n": "clones_nonfunc"
    }

    exclude_samples = set()

    if clonoset_filter is not None:
        if isinstance(clonoset_filter.downsample_size, int):
            need_downsample = True
        if isinstance(clonoset_filter.top, int):
            need_top = True
        if need_downsample or need_top:
            print("Calcultating stats for original clonosets\n" + "_"*41)
            stats = calc_clonoset_stats(clonosets_df)
            downsample_column = count_column_by_umi_and_functionality[clonoset_filter.by_umi][clonoset_filter.functionality]
            top_column = clone_column_by_functionality[clonoset_filter.functionality]
        if need_downsample:
            if stats[downsample_column].isnull().any().any():
                nan_downsample_samples = list(stats[stats[downsample_column].isna()].sample_id)
                print(f"WARNING! Following samples have NaN downsample counts ('{downsample_column}'): {', '.join(nan_downsample_samples)}")
                print("These samples will be excluded from further calculations.")
                exclude_samples.update(nan_downsample_samples)
            not_enough_count_df = stats[stats[downsample_column] < clonoset_filter.downsample_size]
            if len(not_enough_count_df) > 0:
                not_enough_count_samples = list(not_enough_count_df.sample_id)
                if not drop_small_samples:
                    print(f"WARNING! Following samples have not enough downsample counts ('{downsample_column}' < {clonoset_filter.downsample_size}): {', '.join(not_enough_count_samples)}")
                    print("These samples will be excluded from further calculations.")
                    print("To suppress the warnings set drop_small_samples=True")
                exclude_samples.update(not_enough_count_samples)
        if need_top:
            if stats[downsample_column].isnull().any().any():
                nan_downsample_samples = list(stats[stats[downsample_column].isna()].sample_id)
                print(f"WARNING! Following samples have NaN counts ('{downsample_column}'): {', '.join(nan_downsample_samples)}")
                print("These samples will be excluded from further calculations.")
                exclude_samples.update(nan_downsample_samples)
            not_enough_clones_df = stats[stats[top_column] < clonoset_filter.top]
            if len(not_enough_clones_df) > 0:
                not_enough_clones_samples = list(not_enough_clones_df.sample_id)
                if not drop_small_samples:
                    print(f"WARNING! Following samples have not enough clonotypes ('{top_column}' < {clonoset_filter.top}): {', '.join(not_enough_clones_samples)}")
                    print("These samples will be excluded from further calculations.")
                    print("To suppress the warnings set drop_small_samples=True")
                exclude_samples.update(not_enough_clones_samples)


        if isinstance(clonoset_filter.downsample_size, int) or (isinstance(clonoset_filter.top, int) and clonoset_filter.mix_tails):
            random_filter = True


    if random_filter and seed is None:
        print("WARNING! Random filter is applied, but random seed is not set. This may lead to non-reproducible results.")
        print("You may set the seed (of any hashable type) by specifying 'seed='")
    
    tasks = []
    for i,r in clonosets_df.iterrows():
        sample_id = r["sample_id"]
        filename = r["filename"]
        if sample_id in exclude_samples:
            continue
        if clonoset_filter is not None:
            task = (sample_id, filename, calc_function, clonoset_filter.spawn(), iterations, seed, program_name, random_filter)
        else:
            task = (sample_id, filename, calc_function, Filter(), iterations, seed, program_name, random_filter)
        tasks.append(task)
    
    results = run_parallel_calculation(perform_generic_calculation_mp, tasks, program_name, object_name="calcultaion(s)")
    clonosets_df = clonosets_df[columns_retain]
    df = clonosets_df.merge(pd.DataFrame(results), how="left")
    if split_chain_after_calculation:
        df["sample_id"] = df["sample_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
    return df
    
def perform_generic_calculation_mp(args):
    (sample_id, filename, calc_function, clonoset_filter, iterations, seed, program_name, random_filter) = args
    clonoset = read_clonoset(filename)
    colnames = get_column_names_from_clonoset(clonoset)
        
    if random_filter and isinstance(seed, int):
        random.seed(seed)
      
    clonoset_result = {"sample_id": sample_id}
    clonoset_results = []
    filtered_clonosets = []
    
    for i in range(iterations):
        if clonoset_filter is not None:
            filtered_clonoset = clonoset_filter.apply(clonoset, colnames=colnames)
        else:
            filtered_clonoset = clonoset
        filtered_clonosets.append(filtered_clonoset)
    
    for filtered_clonoset in filtered_clonosets:
        clonoset_results.append(calc_function(filtered_clonoset))
    clonoset_result.update(pd.DataFrame(clonoset_results).mean().to_dict())
    
    return clonoset_result


