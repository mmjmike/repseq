import logomaker
import pandas as pd
from .clone_filter import Filter

def create_motif_dict(seq, seq_type="dna", weight=1):
    if seq_type == 'dna':
        variants = ["A", "T", "G", "C"]
    elif seq_type == 'rna':
        variants = ["A", "U", "G", "C"]
    elif seq_type == 'prot':
        variants = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    else:
        print(f"unknown seq type '{seq_type}'")
        return 
    motif_dict = {var: [weight if seq[i]==var else 0 for i in range(len(seq))] for var in variants}
    return motif_dict    
        
def sum_motif_dicts(list_of_dicts):
    variants = list(list_of_dicts[0].keys())
    seq_len = len(list_of_dicts[0][list(list_of_dicts[0].keys())[0]])
    sum_motif_dict = {var: [0 for i in range(seq_len)] for var in variants}
    for motif_dict in list_of_dicts:
        for var in variants:
            sum_motif_dict[var] = [a + b for a, b in zip(sum_motif_dict[var], motif_dict[var])]
    return sum_motif_dict

def normalize_motif_dict(motif_dict):
    #print(motif_dict)
    variants = list(motif_dict.keys())
    seq_len = len(motif_dict[list(motif_dict.keys())[0]])
    sum_list = [sum([motif_dict[var][i] for var in variants]) for i in range(seq_len)]
    norm_motif_dict = dict()
    for var in variants:
        norm_motif_dict[var] = [a/b for a, b in zip(motif_dict[var], sum_list)]
    return norm_motif_dict
            
def get_logo_for_clonoset(clonoset_df, weight_freq=False, seq_type="dna", plot=True):
    motif_dicts = []
    clonoset = Filter(by_umi=True).apply(clonoset_df)

    for i, r in clonoset_df.iterrows():
        weight = 1
        if weight_freq:
            weight = r["freq"]
        column = "cdr3nt"
        if seq_type == "prot":
            column = "cdr3aa"
        seq = r[column]
        motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
        motif_dicts.append(motif_dict)
    #print(len(motif_dicts))
    motif_dict_sum = sum_motif_dicts(motif_dicts)
    motif_dict_sum = normalize_motif_dict(motif_dict_sum)
    info_matrix = logomaker.transform_matrix(pd.DataFrame(motif_dict_sum), 
                                  from_type='probability', 
                                  to_type='information')
    if plot:
        logomaker.Logo(info_matrix)
        #plt.show()
    else:
        return pd.DataFrame(motif_dict_sum)

def get_logo_for_clonoset(clonoset_df, weight_freq=False, seq_type="dna", plot=True):
    
    clonoset = Filter(by_umi=True).apply(clonoset_df)

    column = "cdr3nt"
    if seq_type == "prot":
        column = "cdr3aa"

    list_of_clonotypes = []

    for i, r in clonoset.iterrows():
        seq = r[column]
        if weight_freq:
            weight = r["freq"]
            clone = (seq, weight)
        else:
            clone = (seq,)
        list_of_clonotypes.append(clone)

    get_logo_for_list_of_clonotypes(list_of_clonotypes, seq_type, plot=plot)

def get_logo_for_list_of_clonotypes(list_of_clonotypes, seq_type, plot=True):
    motif_dicts = []
    for clone in list_of_clonotypes:
        weight = 1
        if len(clone) > 1:
            weight = clone[1]
        seq = clone[0]
        motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
        motif_dicts.append(motif_dict)
    #print(len(motif_dicts))
    motif_dict_sum = sum_motif_dicts(motif_dicts)
    motif_dict_sum = normalize_motif_dict(motif_dict_sum)
    info_matrix = logomaker.transform_matrix(pd.DataFrame(motif_dict_sum), 
                                  from_type='probability', 
                                  to_type='information')
    if plot:
        logomaker.Logo(info_matrix)
        #plt.show()
    else:
        return pd.DataFrame(motif_dict_sum)

def get_consensus_from_motif_dict(motif_dict):
    seq_len = len(motif_dict[list(motif_dict.keys())[0]])
    seq = ""
    for i in range(seq_len):
        best_char = "-"
        max_probability = 0
        for char, probabilities in motif_dict.items():
            if probabilities[i] > max_probability:
                max_probability = probabilities[i]
                best_char = char
        seq += best_char
    return seq