import pandas as pd
import numpy as np
from itertools import combinations
from math import comb
from .common_functions import run_parallel_calculation, overlap_type_to_flags
from .io import read_clonoset
from repseq import clone_filter as clf
import olga
import olga.load_model as load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen
import os
import sys


def calculate_clonotypes_pgen(clonosets_df, cl_filter=None, overlap_type='aaVJ', mismatches=1, generation_model='human_T_beta', olga_warnings=False):
    
    '''
    calculates probability of clonotype generation using OLGA model
    
    Args:
        clonosets_df (pd.DataFrame): An output of read_clonoset() function 
        overlap_type (str): Possible values are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`. aa/nt define which CDR3 sequence
            to use (amino acid or nucleotide). V/J in the overlap_type define whether to check V or J segments
            to decide if clonotypes are equal
        mismatches (int): number of mismatches in heighbors CDR3 sequence, default=1
        generation_model (str): generation model used by OLGA. Possible values are `human_B_heavy`, `human_B_kappa`, 
            `human_B_lambda`, `human_T_alpha`, `human_T_beta`, `mouse_T_alpha`, `mouse_T_beta`, default='human_T_beta'
            

    Returns: 
        clonosets_df (pd.DataFrame): a dataframe with calculated p_gens for each clonotype
    '''
    aa, check_v, check_j = overlap_type_to_flags(overlap_type)
    pgen_model = create_olga_model(generation_model)
    
    if not cl_filter:
        cl_filter = clf.Filter(functionality="a")
    clns_olga = cl_filter.apply(clonosets_df)
    
    # creating tuples of size 1-3 each according to the overlap_type 
    overlap_type_vars = []
    if aa:
        overlap_type_vars.append('cdr3aa')
    else:
        overlap_type_vars.append('cdr3nt')
    if check_v:
        overlap_type_vars.append('v')
    if check_j:
        overlap_type_vars.append('j')
        
    # creating neighbours, a {(cdr3, v if required, j if required): [neighbours]} dict, and a {neighbour: [(cdr3, v if required, j if required)]} dict
    clonoset_seqs = list(clns_olga[overlap_type_vars].itertuples(index=False, name=None))
    clonoset_seqs_for_calculation = []
    if len(overlap_type_vars) != 3:
        for cln in clonoset_seqs:
            if len(cln) == 2:
                seq, segment = cln
                v, j = segment, segment
            elif len(cln) == 1:
                seq = cln[0]
            if not check_v:
                v = None
            if not check_j:
                j = None
            clonoset_seqs_for_calculation.append((seq, v, j))
    else:
        clonoset_seqs_for_calculation = clonoset_seqs
        
    # create a temporary column to merge p_gen values with respective clonotypes
    clonoset_total_pgen = clns_olga.copy()
    clonoset_total_pgen['temp'] = clonoset_seqs
    neighbours, cln_neighbours_dict = create_neighbours(clonoset_seqs, mismatches, check_v, check_j)
        
    clonotype_pgen_dict = {}
    neighbours_pgen_dict = {}
    
    n_chunks = 40
    chunk_size = len(clonoset_seqs_for_calculation) // n_chunks + 1
    tasks_clns = []
    for i in range(0, len(clonoset_seqs_for_calculation), chunk_size):
        if i + chunk_size < len(clonoset_seqs_for_calculation):
            task = (clonoset_seqs_for_calculation[i:i + chunk_size], aa, pgen_model, olga_warnings)
        else:
            task = (clonoset_seqs_for_calculation[i:], aa, pgen_model, olga_warnings)
        tasks_clns.append(task)
        
    clonotype_pgen = run_parallel_calculation(calculate_pgen_mp, 
                                              tasks_clns, 
                                              'Calculating p_gen using OLGA', 
                                              'chunks_clonoset')
    for el in clonotype_pgen:
        clonotype_pgen_dict.update(el)
    if len(overlap_type_vars) != 3:
        clonotype_pgen_dict = {tuple([x for x in ct if x is not None]):pgen for ct, pgen in clonotype_pgen_dict.items()}
    
    n_chunks = 40
    chunk_size = len(neighbours) // n_chunks + 1
    tasks_nb = []
    for i in range(0, len(neighbours), chunk_size):
        if i + chunk_size < len(neighbours):
            task = (neighbours[i:i + chunk_size], aa, pgen_model, olga_warnings)
        else:
            task = (neighbours[i:], aa, pgen_model, olga_warnings)
        tasks_nb.append(task)
        
    neighbours_pgen = run_parallel_calculation(calculate_pgen_mp, 
                                                tasks_nb, 
                                                'Calculating p_gen using OLGA', 
                                                'chunks_neighbors')
    for el in neighbours_pgen:
        neighbours_pgen_dict.update(el)
    
    # matching (neighbor, pgen) pairs with their respective clonotypes
    clonotype_neighbours_pgen = {ct: [(n, neighbours_pgen_dict[n]) for n in ns] 
                                 for ct, ns in cln_neighbours_dict.items()}
    
    # pgen = sum(neighbors_pgen) - (C(len(cdr3),mismatches) - 1) * pgen(clonotype)
    clonoset_total_pgen_to_add = {'temp': [], 'p_gen': []}
    for ct, ns in clonotype_neighbours_pgen.items():
        ns_pgen_sum = sum([pgen[1] for pgen in ns])
        n_overlaps = comb(len(ct[0]), mismatches) - 1
        pgen_ct = clonotype_pgen_dict[ct]
        clonoset_total_pgen_to_add['temp'].append(ct)
        clonoset_total_pgen_to_add['p_gen'].append(ns_pgen_sum - n_overlaps * pgen_ct)

    clonoset_total_pgen_to_add = pd.DataFrame(clonoset_total_pgen_to_add)
    
    result_df = pd.merge(clonoset_total_pgen_to_add, clonoset_total_pgen, how='right', on=['temp'])
    result_df.drop('temp', axis=1, inplace=True)
    pgen = result_df.pop('p_gen')
    result_df['pgen'] = pgen
    return result_df
    
    
def create_neighbours(data, mismatches=1, check_v=True, check_j=True):
    '''
    Generates clonotype neighbours with the desired number of  mismatches (default=1) and overlap_type.
    Creates a dictionary of (cdr3 + segments)-heighbours pairs and a list of unique neighbours.
    Note that OLGA does not support regular expressions or custom_alphabet in pgen calculation for nucleotide sequences!
    '''
    res = set()
    cln_neigh_dict = {key: set() for key in data}
    for cln in data:
        if len(cln) == 3:
            seq, v, j = cln
        elif len(cln) == 2:
            seq, segment = cln
            v, j = segment, segment
        elif len(cln) == 1:
            seq = cln[0]
        mismatch_pos = combinations(range(len(seq)), mismatches)
        for pos in mismatch_pos:
            seq_list = list(seq)
            for mm in pos:
                seq_list[mm] = 'X' 
            seq_neighbour = ''.join(seq_list)
            # example: if overlap_type is 'aaV', a single neighbor will be denoted as (aa sequence, V, None) tuple
            if not check_v:
                v = None
            if not check_j:
                j = None
            neighbour_x = (seq_neighbour, v, j)
            res.add(neighbour_x)
            cln_neigh_dict[cln].add(neighbour_x)
    cln_neigh_dict = {cln: list(neigh) for cln, neigh in cln_neigh_dict.items()}
    print(f'Created {len(res)} neighbors')
    return list(res), cln_neigh_dict


def create_olga_model(generation_model='human_T_beta'):
        
    olga_path = os.path.dirname(sys.modules['olga'].__file__)
    params_file_name = os.path.join(olga_path, f'default_models/{generation_model}/model_params.txt')
    marginals_file_name = os.path.join(olga_path, f'default_models/{generation_model}/model_marginals.txt')
    V_anchor_pos_file = os.path.join(olga_path, f'default_models/{generation_model}/V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(olga_path, f'default_models/{generation_model}/J_gene_CDR3_anchors.csv')

    #loading OLGA model
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)

    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file_name)
        
    pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)
    
    return pgen_model
    
    
def calculate_pgen_mp(args):
    olga_inputs, aa, pgen_model, olga_warnings = args
    result = {}
    for olga_input in olga_inputs:
        if aa:
            res = {olga_input: pgen_model.compute_aa_CDR3_pgen(*olga_input, olga_warnings)}
        else:
            res = {olga_input: pgen_model.compute_nt_CDR3_pgen(*olga_input, olga_warnings)}
        result.update(res)
    return result


