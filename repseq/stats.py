import pandas as pd
import numpy as np
import random
import os
import math


from .common_functions import (get_column_names_from_clonoset,
                               filter_by_functionality,
                               calc_insert_size,
                               extract_refpoint_position,
                               center_5,
                               diversity_metrics,
                               run_parallel_calculation)
from .io import read_clonoset
from repseq.clone_filter import Filter

REPSEQ_PATH = os.path.join(os.path.expanduser("~"), "soft", "repseq")
AA_PROPS_PATH = os.path.join(REPSEQ_PATH, "repseq", "resourses", "aa_property_table.txt")


def calc_clonoset_stats(clonosets_df, cl_filter=None, verbose=True, cpu=None):
    """
    Calculates statistics for clonosets regarding clonotype, read and UMI counts.
    Also gives counts for functional clonotypes and non-singletons: clonotypes, 
    having only one count (UMI count if present, read count in other cases).
    Clonosets are given to the function in a form of pd.DataFrame

    Args:
        clonosets_df (pd.DataFrame): dataframe, containing two required columns: 
            `sample_id` and `filename`. Also recommended to have `chain` column in this DF.
        cl_filter (Filter): clonoset filter - object from `clone_filter.py` module.
        verbose (bool): if `True`, show progress information.
        cpu (int, optional): number of worker processes. Use `1` to run
            sequentially in a for-loop for easier debugging. `None` uses the
            default `ProcessPoolExecutor` worker count.

    Returns:
        pd.DataFrame: dataframe with clonotype statistics for each sample in clonosets_df
    """

    df = generic_calculation(
        clonosets_df,
        calculate_clonoset_stats_cl,
        clonoset_filter=cl_filter,
        program_name="CalcClonosetStats",
        verbose=verbose,
        cpu=cpu,
    )
    convert_dict = {"clones": int,
                    "clones_func": int,
                    "clones_func_singletons": int,
                    "clones_func_non_singletons": int,
                    "clones_nonfunc": int,
                    "reads": int,
                    "reads_func": int,
                    "reads_nonfunc": int}
    if not df["umi"].isnull().values.any():
        convert_dict.update({"umi": int,
                        "umi_func": int,
                        "umi_nonfunc": int})
 
    df = df.astype(convert_dict)
    df["reads_per_umi"] = (df["reads"] / df["umi"]).round(2).where(
        df["umi"].notna() & df["umi"].ne(0)
    )
    return df

def calc_segment_usage(clonosets_df, segment="v", cl_filter=None, table="long", by_count=False,
                       cpu=None, drop_small_samples=False, verbose=True):
    """
    Calculates segment (`V`, `J`, or `C`) usage for several samples. By default outputs
    a long table with segment identifiers, `usage`, `sample_id`, and `chain`.
    For `vj` and `vjlen`, identifiers are pipe-delimited strings and the long
    table also contains separate `v`, `j`, and (for `vjlen`) integer `len`
    columns.
    It also may take a clone_filter as input: `cl_filter` from `clone_filter` module.


    Args:
        clonosets_df (pd.DataFrame): dataframe, containing two required columns: 
            `sample_id` and `filename`. Also recommended to have `chain` column in this DF.
        segment (str, optional): possible values are `v`, `j` or `c`. Defaults to "v".
        cl_filter (Filter, optional): clonoset filter - object from `clone_filter.py` module.
        table (str, optional): table type - `long` or `wide`. Defaults to "long".
        cpu (int, optional): number of worker processes.
        drop_small_samples (bool): if `True`, drop samples that are smaller
            than `cl_filter.top` or `cl_filter.downsample`.
        verbose (bool): if `True`, show progress information.


    Returns:
        pd.DataFrame: 'long' or 'wide'. If 'long' it contains four columns, as stated in
            the function description. If 'wide' - then it has all possible segments in 
            column names, sample_id's - in rows and usage in each cell in the table.
    """


    if cl_filter is None:
        cl_filter = Filter()
    
    table_options = ["long", "wide"]
    if table not in table_options:
        raise ValueError(f"Unknown value for 'table' parameter. Possible options: {', '.join(table_options)}")

    possible_segments = ["v", "j", "c", "vj", "vlen", "vjlen"]
    segment = segment.lower()
    if segment not in possible_segments:
        raise ValueError(f"Wrong segment value. Possible values: {', '.join(possible_segments)}")
    df = generic_calculation(
        clonosets_df,
        calc_segment_usage_cl,
        clonoset_filter=cl_filter,
        program_name="CalcSegmentUsage",
        segment=segment,
        by_count=by_count,
        cpu=cpu,
        drop_small_samples=drop_small_samples,
        verbose=verbose,
    )
    df = df.fillna(0)
    if table == "wide":
        return df

    long_df = df.melt(id_vars=["sample_id", "chain"]).rename(
        columns={"value": "usage", "variable": segment}
    )
    if segment not in {"vj", "vjlen"}:
        return long_df

    parts = long_df[segment].astype("string").str.split("|", expand=True)
    expected_parts = 2 if segment == "vj" else 3
    if parts.shape[1] != expected_parts or parts.isna().any().any():
        raise ValueError(f"Invalid pipe-delimited {segment} identifier")
    long_df[segment] = long_df[segment].astype("string")
    long_df["v"] = parts[0].astype("string")
    long_df["j"] = parts[1].astype("string")
    if segment == "vjlen":
        long_df["len"] = pd.to_numeric(parts[2], errors="raise").astype(int)
        return long_df[["sample_id", "chain", "v", "j", "len", segment, "usage"]]
    return long_df[["sample_id", "chain", "v", "j", segment, "usage"]]


def cdr3_length_distributions(
    clonosets_df,
    cl_filter=None,
    cpu=None,
    count_by_freq=True,
    seq_type="aa",
    table="long",
    zero_fill=True,
    verbose=True,
    drop_small_samples=False,
):
    """
    Calculate CDR3 length distributions for multiple clonosets.

    Args:
        clonosets_df (pd.DataFrame): dataframe containing `sample_id` and
            `filename` columns. If `chain` is present, it is retained in the
            result.
        cl_filter (Filter, optional): clonoset filter from `clone_filter.py`.
        cpu (int, optional): number of worker processes. Use `1` to run
            sequentially in a for-loop for easier debugging. `None` uses the
            default `ProcessPoolExecutor` worker count.
        count_by_freq (bool): if `True`, sum clonotype frequencies. If
            `False`, sum clonotype counts.
        seq_type (str): `aa` for amino-acid CDR3 lengths or `nt` for
            nucleotide CDR3 lengths.
        table (str): `long` or `wide`.
        zero_fill (bool): in long tables, include zero-valued rows for CDR3
            lengths found in other samples but absent from a given sample. Wide
            tables are always zero-filled.
        verbose (bool): if `True`, show progress information.
        drop_small_samples (bool): if `True`, drop samples that are smaller
            than `cl_filter.top` or `cl_filter.downsample`.

    Returns:
        pd.DataFrame: long table with `sample_id`, optional `chain`,
            `cdr3_length`, and `freq`/`count`; or wide table with one column
            per CDR3 length.
    """
    table_options = ["long", "wide"]
    if table not in table_options:
        raise ValueError(f"Unknown value for 'table' parameter. Possible options: {', '.join(table_options)}")
    seq_type_options = ["aa", "nt"]
    seq_type = seq_type.lower()
    if seq_type not in seq_type_options:
        raise ValueError(f"Unknown value for 'seq_type' parameter. Possible options: {', '.join(seq_type_options)}")

    df = generic_calculation(
        clonosets_df,
        calc_cdr3_length_distribution_cl,
        clonoset_filter=cl_filter,
        program_name="CalcCDR3LengthDistribution",
        verbose=verbose,
        cpu=cpu,
        drop_small_samples=drop_small_samples,
        count_by_freq=count_by_freq,
        seq_type=seq_type,
    )
    id_vars = ["sample_id"]
    if "chain" in df.columns:
        id_vars.append("chain")
    length_columns = sorted([column for column in df.columns if column not in id_vars])
    df = df[id_vars + length_columns]
    value_column = "freq" if count_by_freq else "count"

    if table == "wide":
        return df.fillna(0)

    if zero_fill:
        df = df.fillna(0)
    result = df.melt(
        id_vars=id_vars,
        value_vars=length_columns,
        var_name="cdr3_length",
        value_name=value_column,
    )
    if not zero_fill:
        result = result.dropna(subset=[value_column])
    result["cdr3_length"] = result["cdr3_length"].astype(int)
    sort_columns = id_vars + ["cdr3_length"]
    return result.sort_values(by=sort_columns).reset_index(drop=True)


def calc_cdr3_length_distribution_cl(clonoset_in, colnames=None, count_by_freq=True, seq_type="aa"):
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset_in)
    if seq_type == "aa":
        sequence_column = colnames["cdr3aa_column"]
    elif seq_type == "nt":
        sequence_column = colnames["cdr3nt_column"]
    else:
        raise ValueError("seq_type must be either 'aa' or 'nt'")
    if sequence_column is None:
        raise ValueError(f"Could not find CDR3 {seq_type} sequence column")

    weight_column = colnames["fraction_column"] if count_by_freq else colnames["count_column"]
    if weight_column is None:
        raise ValueError("Could not find clonotype frequency column" if count_by_freq else "Could not find clonotype count column")

    clonoset = clonoset_in.copy()
    clonoset["cdr3_length"] = clonoset[sequence_column].fillna("").apply(len)
    return clonoset[[weight_column, "cdr3_length"]].groupby("cdr3_length").sum().to_dict()[weight_column]


def calc_segment_usage_cl(clonoset_in, segment="v", colnames=None, by_count=False):
    if segment == "vj":
        return calc_vjlen_usage_cl(clonoset_in, colnames=None, include_j=True, include_len=False, by_count=by_count)
    elif segment == "vlen":
        return calc_vjlen_usage_cl(clonoset_in, colnames=None, include_j=False, include_len=True, by_count=by_count)
    elif segment == "vjlen":
        return calc_vjlen_usage_cl(clonoset_in, colnames=None, include_j=True, include_len=True, by_count=by_count)
    colnames = get_column_names_from_clonoset(clonoset_in)
    freq_column = colnames["fraction_column"]
    if by_count:
        freq_column = colnames["count_column"]
    segment_column = colnames[f"{segment}_column"]
    result = clonoset_in[[freq_column, segment_column]].groupby(segment_column).sum().to_dict()[freq_column]
    return result

def calc_vjlen_usage_cl(clonoset_in, colnames=None, include_j=True, include_len=False, by_count=False):
    if not include_j and not include_len:
        print("WARNING! 'calc_vjlen_usage_cl' must have at least one of the flags 'include_j' or 'include_len' equal to 'True'. Calling 'calc_segment_usage_cl' with segment='v' instead")
        return calc_segment_usage_cl(clonoset_in)
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset_in)
    freq_column = colnames["fraction_column"]
    if by_count:
        freq_column = colnames["count_column"]
    clonoset = clonoset_in.copy()
    v_column = colnames["v_column"]
    columns_to_join = [v_column]
    if include_j:
        j_column = colnames["j_column"]
        columns_to_join.append(j_column)
    if include_len:
        cdr3aa_column = colnames["cdr3aa_column"]
        aalen_column = "cdr3aa_len"
        columns_to_join.append(aalen_column)
        clonoset[aalen_column] = clonoset[cdr3aa_column].apply(lambda x: len(x))
    
    result = clonoset[[freq_column] + columns_to_join].groupby(columns_to_join).sum().to_dict()[freq_column]
    if include_j:
        return {
            "|".join(map(str, key if isinstance(key, tuple) else (key,))): value
            for key, value in result.items()
        }
    if include_len:
        return {
            "|".join(map(str, key if isinstance(key, tuple) else (key,))): value
            for key, value in result.items()
        }
    return result

def _print_normalization_message(function_name, cl_filter, seed, iterations):
    downsample = None if cl_filter is None else cl_filter.downsample_size
    top = None if cl_filter is None else cl_filter.top
    if _should_print_downsampling_recommendation(cl_filter):
        print(
            f"{function_name}: for the most honest comparison, use a filter with "
            "equal downsampling for all samples. This is preferred to top=N filters "
            "and much preferred to calculations without normalization."
        )
    print(
        f"{function_name} settings: seed={seed}, downsample={downsample}, "
        f"top={top}, iterations={iterations}"
    )


def _should_print_downsampling_recommendation(cl_filter):
    if cl_filter is None:
        return True
    if cl_filter.downsample_size is None:
        return True
    return (
        len(cl_filter.white_list) > 0
        or len(cl_filter.black_list) > 0
        or cl_filter.top is not None
        or cl_filter.count_threshold is not None
    )


def calc_diversity_stats(clonosets_df, cl_filter=None, iterations=3, seed=None,
                         drop_small_samples=False, cpu=None, verbose=True):
    """
    Calculates richness, diversity, evenness, and dominance metrics for each
    clonoset in `clonosets_df`.

    The first metric columns are `diversity`, `norm_shannon_wiener`,
    `clonality`, `shannon_wiener`, and `chao1`. Additional columns include
    `richness`, `ace`, `goods_coverage`, `d50`, `simpson`,
    `inverse_simpson`, `gini_simpson`, `berger_parker`, and
    `gini_coefficient`.

    It is highly recommended to use equal downsampling for all input clonosets
    before comparing diversity metrics.
    
    Args:
        clonosets_df (pd.DataFrame): dataframe, containing two required columns: 
            `sample_id` and `filename`. Also recommended to have `chain` column in this DF.
        segment (str, optional): possible values are `v`, `j` or `c`. Defaults to "v".
        cl_filter (Filter, optional): clonoset filter - object from `clone_filter.py` module.
        table (str, optional): table type - `long` or `wide`. Defaults to "long".


    Returns:
        pd.DataFrame: 'long' or 'wide'. If 'long' it contains four columns, as stated in
            the function description. If 'wide' - then it has all possible segments in 
            column names, sample_id's - in rows and usage in each cell in the table.
    """
    if verbose:
        _print_normalization_message("CalcDiversityStats", cl_filter, seed, iterations)
    df = generic_calculation(
        clonosets_df,
        calculate_diversity_stats_cl,
        clonoset_filter=cl_filter,
        program_name="CalcDiversityStats",
        iterations=iterations,
        seed=seed,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
    )
    return df



def cluster_properties(
    clonosets_df,
    cl_filter=None,
    overlap_type="aaVJ",
    mismatches=1,
    min_cluster_size=2,
    cpu=None,
    verbose=True,
):
    """Calculate per-sample statistics of clonotype cluster sizes.

    Each filtered clonoset is clustered independently. Clusters smaller than
    ``min_cluster_size`` nodes are removed before calculating the mean cluster
    size and standard diversity metrics over retained cluster sizes.

    Args:
        clonosets_df (pd.DataFrame): Sample table containing ``sample_id`` and
            ``filename``; ``chain`` is retained when present.
        cl_filter (Filter | None): Optional clonotype filter.
        overlap_type (str): Clonotype comparison definition used for clustering.
        mismatches (int): Maximum sequence mismatches between neighbours.
        min_cluster_size (int): Minimum retained cluster node count.
        cpu (int | None): Worker count used across samples.
        verbose (bool): Show calculation progress.

    Returns:
        pd.DataFrame: One row per clonoset with ``mean_cluster_size`` and the
        standard diversity metrics calculated from retained cluster sizes.
    """
    if (
        not isinstance(min_cluster_size, (int, np.integer))
        or isinstance(min_cluster_size, bool)
        or min_cluster_size < 1
    ):
        raise ValueError("min_cluster_size must be a positive integer.")
    return generic_calculation(
        clonosets_df,
        calculate_cluster_properties_cl,
        clonoset_filter=cl_filter,
        program_name="ClusterProperties",
        overlap_type=overlap_type,
        mismatches=mismatches,
        min_cluster_size=int(min_cluster_size),
        cpu=cpu,
        verbose=verbose,
    )


def calculate_cluster_properties_cl(
    clonoset_in,
    colnames=None,
    overlap_type="aaVJ",
    mismatches=1,
    min_cluster_size=2,
):
    """Calculate cluster-size statistics for one already-filtered clonoset."""
    from .clustering import Clusters, cluster_size

    if (
        not isinstance(min_cluster_size, (int, np.integer))
        or isinstance(min_cluster_size, bool)
        or min_cluster_size < 1
    ):
        raise ValueError("min_cluster_size must be a positive integer.")
    clonoset = clonoset_in.copy()
    if clonoset.empty:
        result = {"mean_cluster_size": np.nan}
        result.update(diversity_metrics([]))
        return result
    clonoset["sample_id"] = "sample"
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(clonoset, verbosity=False)
    clusters.create_clusters(
        overlap_type=overlap_type,
        mismatches=mismatches,
        cpu=1,
        verbosity=False,
    )
    retained = clusters.filter(cluster_size >= min_cluster_size)
    cluster_sizes = retained.custom_properties([cluster_size])[
        "cluster_size"
    ].to_numpy()

    result = {
        "mean_cluster_size": (
            float(np.mean(cluster_sizes)) if len(cluster_sizes) else np.nan
        )
    }
    result.update(diversity_metrics(cluster_sizes))
    return result


def calc_rarefaction_points(sample_df, cl_filter=None, iterations=3, seed=0,
                            drop_small_samples=False, cpu=None, verbose=True):
    """Calculate observed-diversity points for rarefaction curves.

    Each clonoset is prefiltered once with ``cl_filter``. It is then
    downsampled repeatedly at depths 32, 100, 316, 1000, ... (half orders of
    magnitude) below the filtered total count. The final point always contains
    the complete filtered clonoset count and its observed clonotype diversity.

    Args:
        sample_df (pd.DataFrame): Table containing ``sample_id`` and
            ``filename``. ``chain`` is retained when present.
        cl_filter (Filter, optional): Filter applied before rarefaction.
        iterations (int): Number of deterministic downsampling iterations per
            non-final depth. Defaults to 3.
        seed (hashable): Base seed used to derive a distinct seed for every
            iteration. Defaults to 0.
        drop_small_samples (bool): Passed to :func:`generic_calculation` for
            prefilters containing ``top`` or ``downsample`` limits.
        cpu (int, optional): Number of parallel worker processes.
        verbose (bool): Show progress information.

    Returns:
        pd.DataFrame: Long table with ``sample_id``, optional ``chain``,
        ``rarefaction_depth``, and mean observed ``diversity``.
    """
    if not isinstance(iterations, int) or iterations < 1:
        raise ValueError("iterations must be a positive integer")

    wide = generic_calculation(
        sample_df,
        calculate_rarefaction_points_cl,
        clonoset_filter=cl_filter,
        program_name="CalcRarefactionPoints",
        iterations=1,
        seed=seed,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
        rarefaction_iterations=iterations,
        rarefaction_seed=seed,
    )
    id_columns = ["sample_id"] + (["chain"] if "chain" in wide.columns else [])
    depth_columns = [column for column in wide.columns if column not in id_columns]
    result = wide.melt(
        id_vars=id_columns,
        value_vars=depth_columns,
        var_name="rarefaction_depth",
        value_name="diversity",
    ).dropna(subset=["diversity"])
    result["rarefaction_depth"] = result["rarefaction_depth"].astype(int)
    return result.sort_values(id_columns + ["rarefaction_depth"]).reset_index(drop=True)


def calc_rarefaction_curve(sample_df, cl_filter=None, iterations=3, seed=0,
                           drop_small_samples=False, cpu=None, verbose=True):
    """Alias for :func:`calc_rarefaction_points`."""
    return calc_rarefaction_points(
        sample_df,
        cl_filter=cl_filter,
        iterations=iterations,
        seed=seed,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
    )

def calc_convergence(clonosets_df, cl_filter=None, iterations=3, seed=None,
                     drop_small_samples=False, cpu=None, verbose=True):
    if verbose:
        _print_normalization_message("CalcConvergence", cl_filter, seed, iterations)
    df = generic_calculation(
        clonosets_df,
        calculate_convergence_cl,
        clonoset_filter=cl_filter,
        program_name="CalcConvergence",
        iterations=iterations,
        seed=seed,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
    )
    return df


def calc_cdr3_properties(clonosets_df, cl_filter=None, iterations=1, seed=None,
                         drop_small_samples=False, cpu=None, verbose=True):
    if cl_filter is None:
        print("Clonoset Filter is not set. CDR3 stats will be calculated only for functional clonotypes")
        cl_filter = Filter(functionality="f")
    df = generic_calculation(
        clonosets_df,
        calculate_cdr3_properties_cl,
        clonoset_filter=cl_filter,
        program_name="CalcCDR3aaProperties",
        iterations=iterations,
        seed=seed,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
    )
    return df


def calculate_clonoset_stats_cl(clonoset):
    colnames = get_column_names_from_clonoset(clonoset)
    count_column = colnames["count_column"]
    umi_column = colnames["umi_column"]

    # all stats
    clones_num = len(clonoset)
    read_num = clonoset[count_column].sum()
    umi_count = None
    if colnames["umi"]:
        umi_count = clonoset[umi_column].sum()

    # stats for functional clones
    clonoset = filter_by_functionality(clonoset, colnames=colnames)
    func_clones_num = len(clonoset)
    func_read_num = clonoset[count_column].sum()
    func_umi_count = None
    umi_nonfunc = None
    umi_nonfunc_freq = None
    if colnames["umi"]:
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
        zero_insert_freq = clonoset.loc[clonoset["insert_size"] == 0][fraction_column].sum()
    else:
        insert_size_mean = None
        zero_insert_freq = None

    result = {"mean_cdr3nt_len": nt_len_mean,
              "mean_insert_size": insert_size_mean,
              "zero_insert_freq": zero_insert_freq,
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

    cdr3nt_column = colnames["cdr3nt_column"]
    cdr3aa_column = colnames["cdr3aa_column"]
    v_column = colnames["v_column"]
    j_column = colnames["j_column"]

    convergence = _unique_tuple_ratio(clonoset, [cdr3nt_column], [cdr3aa_column])
    convergence_v = _unique_tuple_ratio(clonoset, [cdr3nt_column, v_column], [cdr3aa_column, v_column])
    convergence_vj = _unique_tuple_ratio(clonoset, [cdr3nt_column, v_column, j_column], [cdr3aa_column, v_column, j_column])
    result = {
        "convergence": convergence,
        "convergence_v": convergence_v,
        "convergence_vj": convergence_vj,
    }
    
    return result


def _unique_tuple_ratio(clonoset, numerator_columns, denominator_columns):
    if any(column is None for column in numerator_columns + denominator_columns):
        return np.nan
    numerator = len(clonoset[numerator_columns].drop_duplicates())
    denominator = len(clonoset[denominator_columns].drop_duplicates())
    if denominator == 0:
        return np.nan
    return round(numerator / denominator, 8)

def calculate_diversity_stats_cl(clonoset_in, colnames=None):
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)

    counts = clonoset[colnames["count_column"]]
    
    result = diversity_metrics(counts)
                    
    return result



def _rarefaction_depths(total_count):
    depths = []
    exponent = 1.5
    while True:
        depth = int(round(10 ** exponent))
        if depth >= total_count:
            break
        if not depths or depth != depths[-1]:
            depths.append(depth)
        exponent += 0.5
    depths.append(int(total_count))
    return depths


def _rarefied_diversity(counts, depth, seed):
    total_count = int(np.sum(counts))
    if depth >= total_count:
        return int(np.sum(counts > 0))
    sampled_positions = random.Random(seed).sample(range(total_count), depth)
    clone_indices = np.searchsorted(np.cumsum(counts), sampled_positions, side="right")
    return int(len(np.unique(clone_indices)))


def calculate_rarefaction_points_cl(clonoset_in, rarefaction_iterations=3,
                                    rarefaction_seed=0, colnames=None):
    """Calculate rarefaction points for one prefiltered clonoset."""
    clonoset = clonoset_in.copy()
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset)
    counts = clonoset[colnames["count_column"]].to_numpy(dtype=float)
    if np.any(counts < 0) or not np.allclose(counts, np.round(counts)):
        raise ValueError("Rarefaction requires non-negative integer clonotype counts")
    counts = np.round(counts).astype(np.int64)
    total_count = int(counts.sum())
    if total_count < 1:
        return {}

    result = {}
    iteration_seeds = [
        rarefaction_seed + iteration if isinstance(rarefaction_seed, int)
        else f"{rarefaction_seed}_{iteration}"
        for iteration in range(rarefaction_iterations)
    ]
    for depth in _rarefaction_depths(total_count):
        if depth == total_count:
            result[depth] = int(np.sum(counts > 0))
        else:
            diversities = [
                _rarefied_diversity(counts, depth, iteration_seed)
                for iteration_seed in iteration_seeds
            ]
            result[depth] = float(np.mean(diversities))
    return result


def clonotypes_coverage(sample_df, cl_filter=None, by_counts=False,
                        drop_small_samples=False, cpu=None, verbose=True):
    """Calculate a half-order histogram of clonotype counts.

    Count bins have upper-bound labels ``1, 3, 10, 32, 100, ...``.
    By default, each clonotype contributes one to its bin. With
    ``by_counts=True``, each clonotype contributes its count instead.

    Args:
        sample_df (pd.DataFrame): Table containing ``sample_id`` and
            ``filename``. ``chain`` is retained when present.
        cl_filter (Filter, optional): Filter applied before calculation.
        by_counts (bool): Sum clonotype counts instead of clonotype numbers.
        drop_small_samples (bool): Passed to :func:`generic_calculation`.
        cpu (int, optional): Number of parallel worker processes.
        verbose (bool): Show progress information.

    Returns:
        pd.DataFrame: Long table with ``sample_id``, optional ``chain``,
        string ``bin``, and numeric ``value`` columns.
    """
    wide = generic_calculation(
        sample_df,
        clonotypes_coverage_cl,
        clonoset_filter=cl_filter,
        program_name="ClonotypesCoverage",
        by_counts=by_counts,
        drop_small_samples=drop_small_samples,
        cpu=cpu,
        verbose=verbose,
    )
    id_columns = ["sample_id"] + (["chain"] if "chain" in wide.columns else [])
    bin_columns = [column for column in wide.columns if column not in id_columns]
    if not bin_columns:
        return pd.DataFrame(columns=id_columns + ["bin", "value"])

    wide[bin_columns] = wide[bin_columns].fillna(0)
    result = wide.melt(
        id_vars=id_columns,
        value_vars=bin_columns,
        var_name="bin",
        value_name="value",
    )
    result["bin"] = result["bin"].astype(str)
    result["_bin_number"] = pd.to_numeric(result["bin"])
    return (
        result.sort_values(id_columns + ["_bin_number"])
        .drop(columns="_bin_number")
        .reset_index(drop=True)
    )


def clonotypes_coverage_cl(clonoset_in, colnames=None, by_counts=False):
    """Calculate half-order clonotype-count bins for one clonoset."""
    if colnames is None:
        colnames = get_column_names_from_clonoset(clonoset_in)
    count_column = colnames["count_column"]
    if count_column is None:
        raise ValueError("Could not find clonotype count column")

    histogram = {}
    for count in clonoset_in[count_column]:
        if pd.isna(count):
            raise ValueError("Clonotype counts must not contain missing values")
        if count < 0:
            raise ValueError("Negative clonotype counts are not supported")
        if count <= 1:
            upper_bound = "1"
        else:
            bin_index = math.floor(2 * math.log10(count))
            upper_bound = str(int(round(10 ** ((bin_index + 1) / 2))))
        increment = count if by_counts else 1
        histogram[upper_bound] = histogram.get(upper_bound, 0) + increment
    return histogram


clonotype_coverage = clonotypes_coverage
clonotype_coverage_cl = clonotypes_coverage_cl

def generic_calculation(clonosets_df_in, calc_function, clonoset_filter=None, program_name="Calculation",
                         iterations=1, seed=None, drop_small_samples=False, verbose=True,
                         skip_checks=False, cpu=None, **kwargs):
    '''
    Main function that applies batch calculations for multiple clonosets
    using a `calc_function`. It checks inputs, checks if clonotype counts are 
    coherent with downsample or top numbers in clonoset filter.
    All filters and calc_function are applied to each clonoset in parallel, if
    several cores are available.
    After all calculations are finished, this function combines the results into 
    one pd.DataFrame.

    Args:
        clonosets_df (pd.DataFrame): dataframe, containing two required columns: 
            `sample_id` and `filename`. Also recommended to have `chain` column in this DF.
            `sample_id`'s or `sample_id`+`chain` combinations must be unique in this DF.
        calc_function (function_name): function, applicable to a single clonoset (pd.DataFrame)
            in VDJtools-like format that returns a dictionary of properties+values as output.
        cl_filter (Filter, optional): clonoset filter - object from `clone_filter.py` module.
        program_name (str): the name of applied calculation, it is shown in the progress-bar
        iterations (int): number of iterations to obtain mean values for calculations when
            random processes (downsampling/mix_tails) in clonoset filter are applied.
            Recommended to use 3-5 iterations.
        seed (hashable): seed for random events (downsampling/mix_tails). It overrides the 
            values, specified in `cl_filter`.
        drop_small_samples (bool): `True` - samples, that can't be downsampled or top-cropped
            because of lack of counts/clonotypes will be dropped before the calculation.
            `False` - small samples will be taken into account, but with fewer counts/clonotypes 
            than those with enough counts/clonotypes.
        cpu (int, optional): number of worker processes. Use `1` to run
            sequentially in a for-loop for easier debugging. `None` uses the
            default `ProcessPoolExecutor` worker count.

    Returns:
        df (pd.DataFrame): resulting DataFrame, with `sample_id` and `chain` columns and properties
            columns for each clonoset.

    '''
    
    
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
            split_chain_after_calculation = True

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

    
    if clonoset_filter is not None and not skip_checks:
        if isinstance(clonoset_filter.downsample_size, int):
            need_downsample = True
        if isinstance(clonoset_filter.top, int):
            need_top = True
        if need_downsample or need_top:
            print("Calcultating stats for original clonosets\n" + "_"*41)
            stats = calc_clonoset_stats(clonosets_df, verbose=verbose, cpu=cpu)
            downsample_column = count_column_by_umi_and_functionality[clonoset_filter.by_umi][clonoset_filter.functionality]
            read_column = count_column_by_umi_and_functionality[False][clonoset_filter.functionality]
            top_column = clone_column_by_functionality[clonoset_filter.functionality]

        if need_downsample:
            count_by_reads_samples = set()
            if stats[downsample_column].isnull().any().any():
                nan_downsample_samples = list(stats[stats[downsample_column].isna()].sample_id)
                print(f"WARNING! Following samples have NaN downsample counts ('{downsample_column}'): {', '.join(nan_downsample_samples)}")
                if clonoset_filter.by_umi:
                    print("These samples will be counted by reads instead")
                    count_by_reads_samples = set(nan_downsample_samples)
                else:
                    print("These samples will be excluded from further calculations.")
                    exclude_samples.update(nan_downsample_samples)
            not_enough_count_df = stats[(stats[downsample_column] < clonoset_filter.downsample_size) & (~stats.sample_id.isin(count_by_reads_samples)) |
                                        (stats[read_column] < clonoset_filter.downsample_size) & (stats.sample_id.isin(count_by_reads_samples))]
            if len(not_enough_count_df) > 0:
                not_enough_count_samples = list(not_enough_count_df.sample_id)
                print(f"WARNING! downsample={clonoset_filter.downsample_size} exceeds available counts for samples: {', '.join(not_enough_count_samples)}")
                if drop_small_samples:
                    print("These samples will be excluded from further calculations.")
                    exclude_samples.update(not_enough_count_samples)
                else:
                    print("These samples will be kept and calculated without downsampling.")
        if need_top:
            if stats[downsample_column].isnull().any().any():
                nan_downsample_samples = list(stats[stats[downsample_column].isna()].sample_id)
                print(f"WARNING! Following samples have NaN counts ('{downsample_column}'): {', '.join(nan_downsample_samples)}")
                if clonoset_filter.by_umi:
                    print("These samples will be counted by reads instead")
                else:
                    print("These samples will be excluded from further calculations.")
                    exclude_samples.update(nan_downsample_samples)
            not_enough_clones_df = stats[stats[top_column] < clonoset_filter.top]
            if len(not_enough_clones_df) > 0:
                not_enough_clones_samples = list(not_enough_clones_df.sample_id)
                print(f"WARNING! top={clonoset_filter.top} exceeds available clonotypes for samples: {', '.join(not_enough_clones_samples)}")
                if drop_small_samples:
                    print("These samples will be excluded from further calculations.")
                    exclude_samples.update(not_enough_clones_samples)
                else:
                    print("These samples will be kept and calculated without top filtering.")


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
            task_filter = clonoset_filter.spawn()
            if not drop_small_samples:
                task_filter.ignore_small_clonosets = True
            task = (sample_id, filename, calc_function, task_filter, iterations, seed, program_name, random_filter, kwargs)
        else:
            task = (sample_id, filename, calc_function, None, iterations, seed, program_name, random_filter, kwargs)
        tasks.append(task)
    
    results = run_parallel_calculation(
        perform_generic_calculation_mp,
        tasks,
        program_name,
        object_name="calcultaion(s)",
        verbose=verbose,
        cpu=cpu,
    )
    clonosets_df = clonosets_df[columns_retain]
    df = clonosets_df.merge(pd.DataFrame(results), how="left")
    if split_chain_after_calculation:
        df["sample_id"] = df["sample_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
    return df
    
def perform_generic_calculation_mp(args):
    '''
    A single-core "worker" for `generic_calculation` function.
    It applies clonoset filter (several times if `iterations` > 1) and performes a
    calculation for a single filtered clonoset.


    Args:
        args (tuple): tuple, containing all required parameters
    
    Returns:
        clonoset_result (dict): a dictionary of parameters, calculated for given clonoset
            and averaged for iterations (if > 1)
    '''
    
    
    (sample_id, filename, calc_function, clonoset_filter, iterations, seed, program_name, random_filter, kwargs) = args
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
        clonoset_results.append(calc_function(filtered_clonoset, **kwargs))
    clonoset_result.update(pd.DataFrame(clonoset_results).mean().to_dict())
    
    return clonoset_result
