import scipy
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
from .common_functions import run_parallel_calculation


def wilcox_diff_expression(count_table, sample_metadata, min_samples=2, count_threshold=0.00005, pval_cutoff=None, cpu=None):
    
    # sample_metadata must contain column 'group' and 'sample_id'
    # sample_id's in sample_metadata must pair with count table column_names
    # sample_id's in sample_metadata must be unique - no repeats
    if "group" not in sample_metadata.columns or "sample_id" not in sample_metadata.columns:
        raise ValueError("'sample_metadata' must contain 'group' column with at least two groups and 'sample_id' column")
    
    # subset sample_metadata for only the samples, that are present in count_table columns
    sample_metadata_subset = sample_metadata[sample_metadata["sample_id"].isin(count_table.columns)].copy()
    
    # select and count unique groups among 
    group_list = sample_metadata_subset["group"].unique()
    groups_num = len(group_list)
    
    if groups_num < 2:
        print("There are less than 2 groups in the dataset. Nothing to calculate")
        return None
    
    elif groups_num == 2:
        
        group1 = group_list[0]
        group2 = group_list[1]
        group1_samples = list(sample_metadata_subset[sample_metadata_subset["group"] == group1]["sample_id"])
        group2_samples = list(sample_metadata_subset[sample_metadata_subset["group"] == group2]["sample_id"])
        pair = f"{group1}_vs_{group2}"
        
        df = calc_mann_whitney_for_group_pair(count_table, pair, group1_samples, group2_samples, count_threshold, min_samples)
        
        result_df_main = df.sort_values(by="p_adj")
    
    # more than two groups
    else:
        
        group_pairs_dict = dict()
        
        # group vs all pairs
        for group in group_list:
            pair = f"{group}_vs_all"
            group_samples = list(sample_metadata_subset[sample_metadata_subset["group"] == group]["sample_id"])
            group2_samples = list(sample_metadata_subset[sample_metadata_subset["group"] != group]["sample_id"])
            group_pairs_dict[pair] = [group_samples, group2_samples]
        
        # create all group pairs
        for i in range(len(group_list)):
            group1 = group_list[i]
            for j in range(len(group_list)-i-1):
                group2 = group_list[i+j+1]
                pair = f"{group1}_vs_{group2}"
                group1_samples = list(sample_metadata_subset[sample_metadata_subset["group"] == group1]["sample_id"])
                group2_samples = list(sample_metadata_subset[sample_metadata_subset["group"] == group2]["sample_id"])
                group_pairs_dict[pair] = [group1_samples, group2_samples]

        # create tasks
        tasks = []
        for pair, group_sample_list in group_pairs_dict.items():
            group1_samples = group_sample_list[0]
            group2_samples = group_sample_list[1]
            
            task = (count_table, pair, group1_samples, group2_samples, count_threshold, min_samples)
            tasks.append(task)
            
        result_dfs = run_parallel_calculation(calc_mann_whitney_for_group_pair_mp, tasks, "Calc Mann-Whitney For Pairs", object_name="pairs", cpu=cpu)
        
        result_df = pd.concat(result_dfs).reset_index(drop=True).merge(count_table)
        
        result_df_main = result_df[result_df.pair.str.contains("_vs_all")].copy()

        # perform FDR muliple test correction
        result_df_main["p_adj"] = multipletests(result_df_main["p"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]

        # In the Group vs. All setting we are interested in only those clonotypes that are greater
        # in the group rather than in All other groups
        result_df_main = result_df_main.query("logFC > 0")
        
        if len(result_df_main) > 0:
            pairs_compare_results = []
            pairs_compare_columns = ["i"]
            for group in group_list:
                pairs_compare_columns += [f"{group}_p", f"{group}_logFC"]
            
            for i,r in result_df_main.iterrows():
                curr_group = r["pair"].split("_vs_all")[0]
                feature_id = r["feature_id"]
                row = [i]
                for group2 in group_list:
                    row += find_p_val_and_logFC_for_pair(feature_id, curr_group, group2, result_df)
                    # print(feature_id, curr_group, group2)
                pairs_compare_results.append(row)
                
            pairs_compare_results_df = pd.DataFrame(pairs_compare_results, columns=pairs_compare_columns).set_index("i")
            result_df_main = result_df_main.merge(pairs_compare_results_df, left_index=True, right_index=True)
    
    # cut the table by pval_cutoff
    if pval_cutoff is not None:
        if pval_cutoff <= 0 or pval_cutoff > 1:
            raise ValueError("pval_cutoff must be in range (0,1]")
        result_df_main = result_df_main[result_df_main["p_adj"] <= pval_cutoff]

    return result_df_main


def find_p_val_and_logFC_for_pair(feature_id, group1, group2, full_result_df):
    if group1 == group2:
        return None, None
    df = full_result_df.copy()
    df = df[(df.pair.str.contains(group1)) & (df.pair.str.contains(group2)) & (df["feature_id"] == feature_id)]
    # print(len(df))
    if len(df) > 0:
        
        # take order of comparison into accout and multiply logFC by -1 if order is reversed
        group1_first = 1
        if df.iloc[0]["pair"].split("_vs_")[1] == group1:
            group1_first = -1
        
        return df.iloc[0]["p"], df.iloc[0]["logFC"]*group1_first
    else:
        return None, None
    
    
def calc_mann_whitney_for_group_pair_mp(args):
    return calc_mann_whitney_for_group_pair(*args)

def calc_mann_whitney_for_group_pair(count_table, pair, group1_samples, group2_samples, count_threshold, min_samples):
    # Alternative hypothesis for mann-whitney test:
    # if it is between group comparison, then perform two-sided test
    alternative = "two-sided"
    # if it is SomeGroup vs. All comparison we are expecting that the group has 'greater' clone count
    if "_vs_all" in pair:
        alternative = "greater"

    df = count_table[group1_samples + group2_samples].copy()
    df = df[df.select_dtypes(include=['float', 'int']).apply(lambda x: x > count_threshold).sum(axis=1) >= min_samples]
    df[["logFC", "U", "p"]] = df.apply(lambda x: mann_whitney_two_groups(x, group1_samples, group2_samples,
                                                                         alternative=alternative), axis=1,
                                                                                                 result_type='expand')
    df["pair"] = pair
    df["p_adj"] = multipletests(df["p"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
    df = df.merge(count_table[["feature_id"]], left_index=True, right_index=True)[["feature_id", "logFC", "U", "p", "p_adj", "pair"]]
    return df

def mann_whitney_two_groups(row, group1_samples, group2_samples, alternative="two-sided"):
    values1 = row[group1_samples]
    values2 = row[group2_samples]
    U,p = scipy.stats.mannwhitneyu(values1, values2, alternative=alternative)
    mean1 = np.mean(values1)
    mean2 = np.mean(values2)
    logFC = logFC_for_two_means(mean1, mean2, default_value=100)
    return logFC, U, p

def logFC_for_two_means(mean1, mean2, default_value=100):
    if mean1 == 0 and mean2 == 0:
        logFC = 0
    elif mean1 == 0:
        logFC = -default_value
    elif mean2 == 0:
        logFC = default_value
    else:
        logFC = np.log(mean1/mean2)
    return logFC