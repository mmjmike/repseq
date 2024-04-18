import pandas as pd
from .common_functions import run_parallel_calculation, combine_metadata_from_folders, print_progress_bar
from .logo import create_motif_dict, sum_motif_dicts, get_consensus_from_motif_dict, get_logo_for_list_of_clonotypes
from .clone_filter import Filter
from .io import read_clonoset

import networkx as nx
from networkx.algorithms import community
import numpy as np
from scipy.sparse import csr_matrix
import os
import json


TCRDIST_BLOSUM = {('A', 'A'): 0,  ('A', 'C'): 4,  ('A', 'D'): 4,  ('A', 'E'): 4,  ('A', 'F'): 4,  ('A', 'G'): 4,  ('A', 'H'): 4,  ('A', 'I'): 4,  ('A', 'K'): 4,  ('A', 'L'): 4,  ('A', 'M'): 4,  ('A', 'N'): 4,  ('A', 'P'): 4,  ('A', 'Q'): 4,  ('A', 'R'): 4,  ('A', 'S'): 3,  ('A', 'T'): 4,  ('A', 'V'): 4,  ('A', 'W'): 4,  ('A', 'Y'): 4,  ('C', 'A'): 4,  ('C', 'C'): 0,  ('C', 'D'): 4,  ('C', 'E'): 4,  ('C', 'F'): 4,  ('C', 'G'): 4,  ('C', 'H'): 4,  ('C', 'I'): 4,  ('C', 'K'): 4,  ('C', 'L'): 4,  ('C', 'M'): 4,  ('C', 'N'): 4,  ('C', 'P'): 4,  ('C', 'Q'): 4,  ('C', 'R'): 4,  ('C', 'S'): 4,  ('C', 'T'): 4,  ('C', 'V'): 4,  ('C', 'W'): 4,  ('C', 'Y'): 4,  ('D', 'A'): 4,  ('D', 'C'): 4,  ('D', 'D'): 0,  ('D', 'E'): 2,  ('D', 'F'): 4,  ('D', 'G'): 4,  ('D', 'H'): 4,  ('D', 'I'): 4,  ('D', 'K'): 4,  ('D', 'L'): 4,  ('D', 'M'): 4,  ('D', 'N'): 3,  ('D', 'P'): 4,  ('D', 'Q'): 4,  ('D', 'R'): 4,  ('D', 'S'): 4,  ('D', 'T'): 4,  ('D', 'V'): 4,  ('D', 'W'): 4,  ('D', 'Y'): 4,  ('E', 'A'): 4,  ('E', 'C'): 4,  ('E', 'D'): 2,  ('E', 'E'): 0,  ('E', 'F'): 4,  ('E', 'G'): 4,  ('E', 'H'): 4,  ('E', 'I'): 4,  ('E', 'K'): 3,  ('E', 'L'): 4,  ('E', 'M'): 4,  ('E', 'N'): 4,  ('E', 'P'): 4,  ('E', 'Q'): 2,  ('E', 'R'): 4,  ('E', 'S'): 4,  ('E', 'T'): 4,  ('E', 'V'): 4,  ('E', 'W'): 4,  ('E', 'Y'): 4,  ('F', 'A'): 4,  ('F', 'C'): 4,  ('F', 'D'): 4,  ('F', 'E'): 4,  ('F', 'F'): 0,  ('F', 'G'): 4,  ('F', 'H'): 4,  ('F', 'I'): 4,  ('F', 'K'): 4,  ('F', 'L'): 4,  ('F', 'M'): 4,  ('F', 'N'): 4,  ('F', 'P'): 4,  ('F', 'Q'): 4,  ('F', 'R'): 4,  ('F', 'S'): 4,  ('F', 'T'): 4,  ('F', 'V'): 4,  ('F', 'W'): 3,  ('F', 'Y'): 1,  ('G', 'A'): 4,  ('G', 'C'): 4,  ('G', 'D'): 4,  ('G', 'E'): 4,  ('G', 'F'): 4,  ('G', 'G'): 0,  ('G', 'H'): 4,  ('G', 'I'): 4,  ('G', 'K'): 4,  ('G', 'L'): 4,  ('G', 'M'): 4,  ('G', 'N'): 4,  ('G', 'P'): 4,  ('G', 'Q'): 4,  ('G', 'R'): 4,  ('G', 'S'): 4,  ('G', 'T'): 4,  ('G', 'V'): 4,  ('G', 'W'): 4,  ('G', 'Y'): 4,  ('H', 'A'): 4,  ('H', 'C'): 4,  ('H', 'D'): 4,  ('H', 'E'): 4,  ('H', 'F'): 4,  ('H', 'G'): 4,  ('H', 'H'): 0,  ('H', 'I'): 4,  ('H', 'K'): 4,  ('H', 'L'): 4,  ('H', 'M'): 4,  ('H', 'N'): 3,  ('H', 'P'): 4,  ('H', 'Q'): 4,  ('H', 'R'): 4,  ('H', 'S'): 4,  ('H', 'T'): 4,  ('H', 'V'): 4,  ('H', 'W'): 4,  ('H', 'Y'): 2,  ('I', 'A'): 4,  ('I', 'C'): 4,  ('I', 'D'): 4,  ('I', 'E'): 4,  ('I', 'F'): 4,  ('I', 'G'): 4,  ('I', 'H'): 4,  ('I', 'I'): 0,  ('I', 'K'): 4,  ('I', 'L'): 2,  ('I', 'M'): 3,  ('I', 'N'): 4,  ('I', 'P'): 4,  ('I', 'Q'): 4,  ('I', 'R'): 4,  ('I', 'S'): 4,  ('I', 'T'): 4,  ('I', 'V'): 1,  ('I', 'W'): 4,  ('I', 'Y'): 4,  ('K', 'A'): 4,  ('K', 'C'): 4,  ('K', 'D'): 4,  ('K', 'E'): 3,  ('K', 'F'): 4,  ('K', 'G'): 4,  ('K', 'H'): 4,  ('K', 'I'): 4,  ('K', 'K'): 0,  ('K', 'L'): 4,  ('K', 'M'): 4,  ('K', 'N'): 4,  ('K', 'P'): 4,  ('K', 'Q'): 3,  ('K', 'R'): 2,  ('K', 'S'): 4,  ('K', 'T'): 4,  ('K', 'V'): 4,  ('K', 'W'): 4,  ('K', 'Y'): 4,  ('L', 'A'): 4,  ('L', 'C'): 4,  ('L', 'D'): 4,  ('L', 'E'): 4,  ('L', 'F'): 4,  ('L', 'G'): 4,  ('L', 'H'): 4,  ('L', 'I'): 2,  ('L', 'K'): 4,  ('L', 'L'): 0,  ('L', 'M'): 2,  ('L', 'N'): 4,  ('L', 'P'): 4,  ('L', 'Q'): 4,  ('L', 'R'): 4,  ('L', 'S'): 4,  ('L', 'T'): 4,  ('L', 'V'): 3,  ('L', 'W'): 4,  ('L', 'Y'): 4,  ('M', 'A'): 4,  ('M', 'C'): 4,  ('M', 'D'): 4,  ('M', 'E'): 4,  ('M', 'F'): 4,  ('M', 'G'): 4,  ('M', 'H'): 4,  ('M', 'I'): 3,  ('M', 'K'): 4,  ('M', 'L'): 2,  ('M', 'M'): 0,  ('M', 'N'): 4,  ('M', 'P'): 4,  ('M', 'Q'): 4,  ('M', 'R'): 4,  ('M', 'S'): 4,  ('M', 'T'): 4,  ('M', 'V'): 3,  ('M', 'W'): 4,  ('M', 'Y'): 4,  ('N', 'A'): 4,  ('N', 'C'): 4,  ('N', 'D'): 3,  ('N', 'E'): 4,  ('N', 'F'): 4,  ('N', 'G'): 4,  ('N', 'H'): 3,  ('N', 'I'): 4,  ('N', 'K'): 4,  ('N', 'L'): 4,  ('N', 'M'): 4,  ('N', 'N'): 0,  ('N', 'P'): 4,  ('N', 'Q'): 4,  ('N', 'R'): 4,  ('N', 'S'): 3,  ('N', 'T'): 4,  ('N', 'V'): 4,  ('N', 'W'): 4,  ('N', 'Y'): 4,  ('P', 'A'): 4,  ('P', 'C'): 4,  ('P', 'D'): 4,  ('P', 'E'): 4,  ('P', 'F'): 4,  ('P', 'G'): 4,  ('P', 'H'): 4,  ('P', 'I'): 4,  ('P', 'K'): 4,  ('P', 'L'): 4,  ('P', 'M'): 4,  ('P', 'N'): 4,  ('P', 'P'): 0,  ('P', 'Q'): 4,  ('P', 'R'): 4,  ('P', 'S'): 4,  ('P', 'T'): 4,  ('P', 'V'): 4,  ('P', 'W'): 4,  ('P', 'Y'): 4,  ('Q', 'A'): 4,  ('Q', 'C'): 4,  ('Q', 'D'): 4,  ('Q', 'E'): 2,  ('Q', 'F'): 4,  ('Q', 'G'): 4,  ('Q', 'H'): 4,  ('Q', 'I'): 4,  ('Q', 'K'): 3,  ('Q', 'L'): 4,  ('Q', 'M'): 4,  ('Q', 'N'): 4,  ('Q', 'P'): 4,  ('Q', 'Q'): 0,  ('Q', 'R'): 3,  ('Q', 'S'): 4,  ('Q', 'T'): 4,  ('Q', 'V'): 4,  ('Q', 'W'): 4,  ('Q', 'Y'): 4,  ('R', 'A'): 4,  ('R', 'C'): 4,  ('R', 'D'): 4,  ('R', 'E'): 4,  ('R', 'F'): 4,  ('R', 'G'): 4,  ('R', 'H'): 4,  ('R', 'I'): 4,  ('R', 'K'): 2,  ('R', 'L'): 4,  ('R', 'M'): 4,  ('R', 'N'): 4,  ('R', 'P'): 4,  ('R', 'Q'): 3,  ('R', 'R'): 0,  ('R', 'S'): 4,  ('R', 'T'): 4,  ('R', 'V'): 4,  ('R', 'W'): 4,  ('R', 'Y'): 4,  ('S', 'A'): 3,  ('S', 'C'): 4,  ('S', 'D'): 4,  ('S', 'E'): 4,  ('S', 'F'): 4,  ('S', 'G'): 4,  ('S', 'H'): 4,  ('S', 'I'): 4,  ('S', 'K'): 4,  ('S', 'L'): 4,  ('S', 'M'): 4,  ('S', 'N'): 3,  ('S', 'P'): 4,  ('S', 'Q'): 4,  ('S', 'R'): 4,  ('S', 'S'): 0,  ('S', 'T'): 3,  ('S', 'V'): 4,  ('S', 'W'): 4,  ('S', 'Y'): 4,  ('T', 'A'): 4,  ('T', 'C'): 4,  ('T', 'D'): 4,  ('T', 'E'): 4,  ('T', 'F'): 4,  ('T', 'G'): 4,  ('T', 'H'): 4,  ('T', 'I'): 4,  ('T', 'K'): 4,  ('T', 'L'): 4,  ('T', 'M'): 4,  ('T', 'N'): 4,  ('T', 'P'): 4,  ('T', 'Q'): 4,  ('T', 'R'): 4,  ('T', 'S'): 3,  ('T', 'T'): 0,  ('T', 'V'): 4,  ('T', 'W'): 4,  ('T', 'Y'): 4,  ('V', 'A'): 4,  ('V', 'C'): 4,  ('V', 'D'): 4,  ('V', 'E'): 4,  ('V', 'F'): 4,  ('V', 'G'): 4,  ('V', 'H'): 4,  ('V', 'I'): 1,  ('V', 'K'): 4,  ('V', 'L'): 3,  ('V', 'M'): 3,  ('V', 'N'): 4,  ('V', 'P'): 4,  ('V', 'Q'): 4,  ('V', 'R'): 4,  ('V', 'S'): 4,  ('V', 'T'): 4,  ('V', 'V'): 0,  ('V', 'W'): 4,  ('V', 'Y'): 4,  ('W', 'A'): 4,  ('W', 'C'): 4,  ('W', 'D'): 4,  ('W', 'E'): 4,  ('W', 'F'): 3,  ('W', 'G'): 4,  ('W', 'H'): 4,  ('W', 'I'): 4,  ('W', 'K'): 4,  ('W', 'L'): 4,  ('W', 'M'): 4,  ('W', 'N'): 4,  ('W', 'P'): 4,  ('W', 'Q'): 4,  ('W', 'R'): 4,  ('W', 'S'): 4,  ('W', 'T'): 4,  ('W', 'V'): 4,  ('W', 'W'): 0,  ('W', 'Y'): 2,  ('Y', 'A'): 4,  ('Y', 'C'): 4,  ('Y', 'D'): 4,  ('Y', 'E'): 4,  ('Y', 'F'): 1,  ('Y', 'G'): 4,  ('Y', 'H'): 2,  ('Y', 'I'): 4,  ('Y', 'K'): 4,  ('Y', 'L'): 4,  ('Y', 'M'): 4,  ('Y', 'N'): 4,  ('Y', 'P'): 4,  ('Y', 'Q'): 4,  ('Y', 'R'): 4,  ('Y', 'S'): 4,  ('Y', 'T'): 4,  ('Y', 'V'): 4,  ('Y', 'W'): 2,  ('Y', 'Y'): 0}
TCRDIST_CDR3_N_CUT = 3
TCRDIST_CDR3_C_CUT = 2
TCRDIST_CDR3_SCORE_MULTIPLIER = 3
with open(os.path.join(os.path.dirname(__file__), 'tcrdist_ab_v_segments.json')) as f:
    for jsonObj in f:
        TCRDIST_V_DIST = json.loads(jsonObj)

class Node:
    
    def __init__(self, node_id, seq_nt, seq_aa, v, j, sample_id, size=1):
        self.id = node_id
        self.v = v
        self.j = j
        self.seq_aa = seq_aa
        self.seq_nt = seq_nt
        self.sample_id = sample_id
        self.size = size
        self.additional_properties = {}
        
    def is_neighbour_of(self, other, mismatches=1, aa=True, check_v=False, check_j=False):
        """function compare two strings and return
        True if their are equal
            or if they have one mismatch and equal length
        False in all other conditions
        """
        
        if aa:
            string1 = self.seq_aa
            string2 = other.seq_aa
        else:
            string1 = self.seq_nt
            string2 = other.seq_nt
        if len(string1) != len(string2):
             return False
        if check_v and self.v != other.v:
            return False
        if check_j and self.j != other.j:
            return False
        hamm_dist = sum([a != b for a,b in zip(string1,string2)]) 
        if hamm_dist > mismatches:
            return False     
        return True
    
    def __str__(self):
        return "{}_{}".format(self.seq_aa, self.id)

    def add_properties(self, metadata):
        if self.sample_id not in metadata:
            for property in list(metadata[list(metadata)[0]].keys()):
                self.additional_properties[property] = None
        else:
            self.additional_properties.update(metadata[self.sample_id])


def pool_clonotypes_from_clonosets_df(clonosets_df, cl_filter=None):
    
    if cl_filter is None:
        cl_filter = Filter()
    
    clonotypes_dfs = []
    
    for index, row in clonosets_df.iterrows():
        sample_id = row["sample_id"]
        filename = row["filename"]
        clonoset = read_clonoset(filename)
        clonoset = cl_filter.apply(clonoset)
        clonoset["sample_id"] = sample_id
        clonotypes_dfs.append(clonoset)

    result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)

    print(f"Pooled {len(result_df)} clonotypes from {len(clonosets_df)} samples")

    return result_df

        
def find_nodes_and_edges(clonoset_input, mismatches=1, overlap_type="aaV", igh=False, count_by_freq=False):
    if isinstance(clonoset_input, str):
        clonoset=pd.read_csv(clonoset_input,sep="\t")
    else:
        clonoset = clonoset_input
        
    possible_overlap_types = ["aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ"]
    if overlap_type not in possible_overlap_types:
        print("Incorrect overlap type. Possible values: {}".format(", ".join(possible_overlap_types)))    
        return None
    aa = False
    if overlap_type[0:2] == "aa":
        aa = True
    check_v = False
    if "V" in overlap_type:
        check_v = True
    check_j = False
    if "J" in overlap_type:
        check_j = True
    
    clonoset = Filter(by_umi=True).apply(clonoset)
    
    # clonoset["v"] = clonoset["v"].apply(lambda x: x.split("*")[0].split("(")[0])
    # clonoset["j"] = clonoset["j"].apply(lambda x: x.split("*")[0].split("(")[0])
    if igh:
        def _split_c(c_segm):
            try:
                return c_segm.split("*")[0]
            except AttributeError:
                return "None"
        clonoset["c"] = clonoset["c"].apply(lambda x: _split_c(x).split("(")[0])

    nodes_by_len={}
    list_of_all_nodes = []

    
    size_column_name = "freq"
    if not count_by_freq:
        size_column_name = "count"
    
    for index, row in clonoset.iterrows():
        v = row["v"]
        j = row["j"]
        cdr3aa = row["cdr3aa"]
        if "cdr3nt" in clonoset.columns:
            cdr3nt = row["cdr3nt"]
        else:
            cdr3nt = "-"
        size = row[size_column_name]
        sample_id = row["sample_id"]
        len_cdr3aa = len(cdr3aa)
        if len_cdr3aa not in nodes_by_len:
            nodes_by_len[len_cdr3aa] = []
        node = Node(index, cdr3nt, cdr3aa, v, j, sample_id, size=size)
        if igh:
            node.additional_properties["c"] = row["c"]
        nodes_by_len[len_cdr3aa].append(node)
        list_of_all_nodes.append(node)
    print("Nodes list created: {} nodes".format(len(list_of_all_nodes)))
    
    tasks = []
    for aa_len in nodes_by_len:
        nodes_list=nodes_by_len[aa_len]
        task = (nodes_list, mismatches, aa, check_v, check_j)
        tasks.append(task)
    
    program_name = "Find neighbour clonotypes"
    result_list = run_parallel_calculation(find_edges_in_nodes_set_mp, tasks, program_name)
    edges=[j for i in result_list for j in i]
    print("found {} edges".format(len(edges)))
    nodes_set = set(list_of_all_nodes)
    connected_nodes_set = set()
    for edge in edges:
        connected_nodes_set.add(edge[0])
        connected_nodes_set.add(edge[1])
    single_nodes_set = nodes_set.difference(connected_nodes_set)
    
    return list(single_nodes_set), edges

def find_edges_in_nodes_set_mp(args):
    (nodes_list, mismatches, aa, check_v, check_j) = args
    edges = []
    
    for i in range(len(nodes_list)-1):
        for j in range(i+1,len(nodes_list)):
            node_1=nodes_list[i]
            node_2=nodes_list[j]
            if node_1.is_neighbour_of(node_2, mismatches, aa, check_v, check_j):
                edges.append((node_1, node_2, 1)) #save unique codes
    return edges

def find_nodes_and_edges_tcrdist_no_gaps(clonoset_input, radius=16, count_by_freq=True, igh=False):
    if isinstance(clonoset_input, str):
        clonoset=pd.read_csv(clonoset_input,sep="\t")
    else:
        clonoset = clonoset_input
    
    clonoset = clonoset.rename(columns={"bestVGene": "v",
                                        "bestJGene": "j",
                                        "CDR3.amino.acid.sequence": "cdr3aa",
                                        "CDR3.nucleotide.sequence": "cdr3nt",
                                        "allVHitsWithScore": "v",
                                        "allJHitsWithScore": "j",
                                        "allCHitsWithScore": "c",
                                        "aaSeqCDR3": "cdr3aa",
                                        "nSeqCDR3": "cdr3nt",
                                        "Sample":"sample_id",
                                        "cloneFraction":"freq",
                                        "Read.count": "count",
                                        "cloneCount": "count",
                                        "uniqueUMICount":"uniqueMoleculeCount",
                                        "uniqueUMIFraction":"uniqueMoleculeFraction"})
    
    count_column = "count"
    fraction_column = "freq"
    if "uniqueMoleculeCount" in clonoset.columns:
        count_column = "uniqueMoleculeCount"
        fraction_column = "uniqueMoleculeFraction"
    elif "readCount" in clonoset.columns:
        count_column = "readCount"
        fraction_column = "readFraction"
    clonoset = clonoset.rename(columns={count_column: "count", fraction_column: "freq"})
    
    clonoset["v"] = clonoset["v"].apply(lambda x: x.split("*")[0])
    clonoset["j"] = clonoset["j"].apply(lambda x: x.split("*")[0])
    if igh:
        def _split_c(c_segm):
            try:
                return c_segm.split("*")[0]
            except AttributeError:
                return "None"
        clonoset["c"] = clonoset["c"].apply(lambda x: _split_c(x))



    nodes_by_len = {}
    list_of_all_nodes = []

    if count_by_freq and "freq" in clonoset.columns:
        size_column_name = "freq"
    else:
        size_column_name = "count"
        if size_column_name not in clonoset.columns:
            clonoset["count"] = 1
    print(size_column_name)

    for index, row in clonoset.iterrows():
        v = row["v"]
        if v not in TCRDIST_V_DIST:
            print(f"Warning! {v} gene was not recognized in reference db no cdr seq could be inferred. The clone was skipped")
            continue
        j = row["j"]
        cdr3aa = row["cdr3aa"]
        if "cdr3nt" in clonoset.columns:
            cdr3nt = row["cdr3nt"]
        else:
            cdr3nt = "-"
        size = row[size_column_name]
        sample_id = row["sample_id"]
        len_cdr3aa = len(cdr3aa)
        if len_cdr3aa not in nodes_by_len:
            nodes_by_len[len_cdr3aa] = []
        node = Node(index, cdr3nt, cdr3aa, v, j, sample_id, size=size)
        if igh:
            node.additional_properties["c"] = row["c"]
        nodes_by_len[len_cdr3aa].append(node)
        list_of_all_nodes.append(node)
    print("Nodes list created: {} nodes".format(len(list_of_all_nodes)))

    tasks = []
    for aa_len in nodes_by_len:
        nodes_list=nodes_by_len[aa_len]
        task = (nodes_list, radius)
        tasks.append(task)
    
    program_name = "Find neighbour clonotypes (TCRdist with no CDR3 gaps)"
    result_list = run_parallel_calculation(find_nodes_and_edges_tcrdist_no_gaps_mp, tasks, program_name)
    edges=[j for i in result_list for j in i]
    print("found {} edges".format(len(edges)))
    nodes_set = set(list_of_all_nodes)
    connected_nodes_set = set()
    for edge in edges:
        connected_nodes_set.add(edge[0])
        connected_nodes_set.add(edge[1])
    single_nodes_set = nodes_set.difference(connected_nodes_set)
    
    return list(single_nodes_set), edges

def find_nodes_and_edges_tcrdist_no_gaps_mp(args):
    (nodes_list, radius) = args
    edges = []
    
    for i in range(len(nodes_list)-1):
        for j in range(i+1,len(nodes_list)):
            node_1=nodes_list[i]
            node_2=nodes_list[j]
            dist = tcr_dist(node_1, node_2, radius)
            if dist >= 0:
                edges.append((node_1, node_2, dist)) #save unique codes
    return edges

def tcr_dist(node_1, node_2, radius):
    dist=TCRDIST_V_DIST[node_1.v][node_2.v]
    if dist > radius:
        return -1
    dist += TCRDIST_CDR3_SCORE_MULTIPLIER*sum([TCRDIST_BLOSUM[(a,b)] for a,b in zip(node_1.seq_aa[TCRDIST_CDR3_N_CUT:-TCRDIST_CDR3_C_CUT],
                                                     node_2.seq_aa[TCRDIST_CDR3_N_CUT:-TCRDIST_CDR3_C_CUT])])
    if dist <= radius:
        return dist
    return -1

def create_clusters_old(clonoset_input, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True):
    tcr_dist=False
    if isinstance(tcrdist_radius, int):
        tcr_dist=True
        nodes, edges = find_nodes_and_edges_tcrdist_no_gaps(clonoset_input, radius=tcrdist_radius, count_by_freq=count_by_freq, igh=igh)
    else:
        nodes, edges = find_nodes_and_edges(clonoset_input, mismatches=mismatches, overlap_type=overlap_type, igh=igh, count_by_freq=count_by_freq)
    
    main_graph = nx.Graph()
    main_graph.add_nodes_from(nodes)
    print("-----------------------------\nNexworkX graph created")

    program_name = "Adding edges"
    edges_done = 0
    edges_total = len(edges)
    print_progress_bar(edges_done, edges_total, program_name=program_name, object_name="edge(s)")
    node_id_dict = {node.id: node for node in main_graph}
    for edge in edges:
        node1 = node_id_dict[edge[0].id]
        node2 = node_id_dict[edge[1].id]
        length = edge[2]
        main_graph.add_edge(node1, node2, length=length)
        edges_done += 1
        print_progress_bar(edges_done, edges_total, program_name=program_name, object_name="edge(s)")
    clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
    total_clusters = len(clusters)
    cluster_num = len(filter_one_node_clusters(clusters))
    singletons = total_clusters - cluster_num
    print(f"Found {cluster_num} clusters (2 or more nodes) and {singletons} single nodes. Total: {total_clusters}")
    clusters.sort(key=lambda x: len(x), reverse=True)
    return clusters


def create_clusters(clonoset_df, cl_filter=None, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True):
    
    if cl_filter is None:
        cl_filter = Filter()

    clonoset_input = pool_clonotypes_from_clonosets_df(clonoset_df, cl_filter=cl_filter)

    tcr_dist=False

    if isinstance(tcrdist_radius, int):
        tcr_dist=True
        nodes, edges = find_nodes_and_edges_tcrdist_no_gaps(clonoset_input, radius=tcrdist_radius, count_by_freq=count_by_freq, igh=igh)
    else:
        nodes, edges = find_nodes_and_edges(clonoset_input, mismatches=mismatches, overlap_type=overlap_type, igh=igh, count_by_freq=count_by_freq)
    
    main_graph = nx.Graph()
    main_graph.add_nodes_from(nodes)
    print("-----------------------------\nNexworkX graph created")

    program_name = "Adding edges..."
    edges_done = 0
    edges_total = len(edges)
    # print_progress_bar(edges_done, edges_total, program_name=program_name, object_name="edge(s)")
    node_id_dict = {node.id: node for node in main_graph}
    for edge in edges:
        node1 = node_id_dict[edge[0].id]
        node2 = node_id_dict[edge[1].id]
        length = edge[2]
        main_graph.add_edge(node1, node2, length=length)
        edges_done += 1
        # print_progress_bar(edges_done, edges_total, program_name=program_name, object_name="edge(s)")
    clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
    total_clusters = len(clusters)
    cluster_num = len(filter_one_node_clusters(clusters))
    singletons = total_clusters - cluster_num
    print(f"Found {cluster_num} clusters (2 or more nodes) and {singletons} single nodes. Total: {total_clusters}")
    clusters.sort(key=lambda x: len(x), reverse=True)
    write_cluster_no_to_nodes(clusters)
    return clusters

def write_cluster_no_to_nodes(clusters):
    cluster_no = 0
    for cluster in clusters:
        for node in cluster:
            node.additional_properties["cluster_no"] = cluster_no
        cluster_no += 1

def find_cluster_communities_louvain(clusters, resolution=1, threshold=1e-07, seed=1):
    total_communities = 0
    new_clusters = []
    for cluster in clusters:
        if len(cluster) < 2:
            for node in cluster:
                node.additional_properties["community"] = total_communities
            total_communities += 1 
            new_clusters.append(cluster)
        else:
            cluster_communities = community.louvain_communities(cluster, resolution=resolution, threshold=threshold, seed=seed)
            for com in cluster_communities:
                com_nodes = []
                for node in com:
                    node.additional_properties["community"] = total_communities
                    com_nodes.append(node)
                total_communities += 1  
                new_clusters.append(cluster.subgraph(com_nodes))
    new_clusters.sort(key=lambda x: len(x), reverse=True)
    return new_clusters

def filter_one_node_clusters(clusters):
    return [c for c in clusters if len(c)>1]

def save_clusters_for_cytoscape(clusters, output_prefix, sample_metadata=None):
    sif_filename = output_prefix + ".sif"
    #properties_filename = output_prefix + ".prop.tsv"
    properties_metadata_filename = output_prefix + ".prop.metadata.tsv"
    edges = []
    nodes = []

    additional_properties=[]
    for node in clusters[0]:
        additional_properties = list(node.additional_properties.keys())
        break

    for cluster in clusters:
        attributes = nx.get_edge_attributes(cluster,'length')
        for u,v in cluster.edges():
            node1_id = str(u)
            node2_id = str(v)
            try:
                length = attributes[(u,v)]
                edges.append(f"{node1_id}\t{length}\t{node2_id}")
            except KeyError:
                try:
                    length = attributes[(v,u)]
                    edges.append(f"{node1_id}\t{length}\t{node2_id}")
                except KeyError:
                    edges.append(f"{node1_id}\ttneighbour\t{node2_id}")
        if len(cluster) == 1:
            list(cluster.nodes())[0].id
            edges.append(str(list(cluster.nodes())[0]))
        for node in cluster:
            add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
            nodes.append((str(node), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.size, *add_properties_values))

    with open(sif_filename, "w") as f:
        f.write("\n".join(edges))
    print("Saved edges to: {}".format(sif_filename))
    
    properties_names = ["code", "cdr3aa", "v", "j", "cdr3nt", "sample_id", "size"] + additional_properties
    properties_df = pd.DataFrame(nodes, columns=properties_names)
    if sample_metadata is not None:
        properties_df = properties_df.merge(sample_metadata)
    properties_df.to_csv(properties_metadata_filename, index=False, sep="\t")
    print("Saved node properties and metadata to: {}".format(properties_metadata_filename))

def add_alice_hits_to_clusters(clusters, alice_hits_df, check_samples=True):
    nt_seq_colname = "CDR3.nucleotide.sequence"
    if nt_seq_colname not in alice_hits_df.columns:
        nt_seq_colname = "cdr3nt"


    if check_samples:
        samples = list(alice_hits_df["sample_id"].unique())
        alice_hits_dict = {}
        for sample in samples:
            hits = set()
            for index, row in alice_hits_df.loc[alice_hits_df["sample_id"] == sample].iterrows():
                v = row["bestVGene"]
                j = row["bestJGene"]
                cdr3nt = row[nt_seq_colname]
                hit = (cdr3nt, v, j)
                hits.add(hit)
            alice_hits_dict[sample] = hits
        for cluster in clusters:
            for node in cluster:
                if node.sample_id in samples and (node.seq_nt, node.v, node.j) in alice_hits_dict[node.sample_id]:
                    node.additional_properties["alice_hit"] = True
                else:
                    node.additional_properties["alice_hit"] = False
    else:
        hits = set()
        for index, row in alice_hits_df.iterrows():
            v = row["bestVGene"]
            j = row["bestJGene"]
            cdr3nt = row[nt_seq_colname]
            hit = (cdr3nt, v, j)
            hits.add(hit)
            for cluster in clusters:
                for node in cluster:
                    if (node.seq_nt, node.v, node.j) in hits:
                        node.additional_properties["alice_hit"] = True
                    else:
                        node.additional_properties["alice_hit"] = False

def filter_clusters_with_alice_hits(clusters):
    filtered_clusters = []
    for cluster in clusters:
        good_node = False
        for node in cluster:
            try:
               a = node.additional_properties["alice_hit"]
            except KeyError:
                print ("ALICE hits are not specified for these clusters. Specify this by ... [not still written]")
                return None
            if node.additional_properties["alice_hit"]:
                good_node = True
                break
        if good_node:
            filtered_clusters.append(cluster)
            continue
    return filtered_clusters

def pool_clonotypes_to_df(folders, samples_list=None, top=0, functional=True, exclude_singletons=False, cdr3aa_len_range=[], metadata_filename="vdjtools_metadata.txt"):
    all_metadata = combine_metadata_from_folders(folders, metadata_filename=metadata_filename)
    #pool_metadata(folders, metadata_filename,samples_list)
    
    clonotypes_dfs = []
    for index, row in all_metadata.iterrows():
        sample_id = row["sample.id"]
        if samples_list is not None:
            if sample_id not in samples_list:
                continue
        clonoset_data=pd.read_csv(row["#file.name"],sep="\t")

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

        if exclude_singletons:
            clonoset_data=clonoset_data.loc[clonoset_data["count"]>1]
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
        sample_id = row["sample.id"]
        # clonotypes_num = clonoset_data.shape[0]
        clonoset_data["sample_id"] = sample_id
        clonotypes_dfs.append(clonoset_data)
#         print("Added {} clonotypes from {}".format(clonotypes_num, sample_id))
    result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)
    clonotypes_number = len(result_df)
    samples_number = len(result_df["sample_id"].unique())
    print("Pooled {} clonotypes from {} samples".format(clonotypes_number, samples_number))
    return result_df

def add_metadata(clusters, metadata):
    columns = metadata.columns
    if "sample_id" not in columns:
        print("Error: 'sample_id' column is compulsory, but is not present in given metadata")
        return 
    metadata = metadata.set_index("sample_id")
    metadata_dict = metadata.to_dict("index")
    for cluster in clusters:
        for node in cluster:
            node.add_properties(metadata_dict)

def pool_alice_hits_to_df(folders, samples_list=None, metadata_filename="vdjtools_metadata.txt"):
    all_metadata = combine_metadata_from_folders(folders, metadata_filename=metadata_filename)
    #pool_metadata(folders, metadata_filename,samples_list)
    
    clonotypes_dfs = []
    for index, row in all_metadata.iterrows():
        sample_id = row["sample.id"]
        if samples_list is not None:
            if sample_id not in samples_list:
                continue
        clonoset_data=pd.read_csv(row["#file.name"],sep="\t")
        clonoset_data["freq"]=clonoset_data["Read.count"]/clonoset_data["Read.count"].sum()
        sample_id = row["sample.id"]
        # clonotypes_num = clonoset_data.shape[0]
        clonoset_data["sample_id"] = sample_id
        clonotypes_dfs.append(clonoset_data)
#         print("Added {} clonotypes from {}".format(clonotypes_num, sample_id))
    result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)
    clonotypes_number = len(result_df)
    samples_number = len(result_df["sample_id"].unique())
    print("Pooled {} clonotypes from {} samples".format(clonotypes_number, samples_number))
    return result_df

def cluster_properties(clusters, weighed=False):
    properties_list = ["nodes", "edges", "diameter", "density", "eccentricity",
                       "concensus_cdr3aa", "concensus_cdr3nt", "concensus_v", "concensus_j"]
    results = []
    for cluster in clusters:
        average_eccentricity = np.mean(list(nx.eccentricity(cluster).values()))
        aa_consensus = calc_cluster_consensus(cluster, seq_type="prot", weighed=weighed)
        nt_consensus = calc_cluster_consensus(cluster, seq_type="dna", weighed=weighed)
        v_consensus = calc_cluster_consensus_segment(cluster, segment_type="v", weighed=weighed)
        j_consensus = calc_cluster_consensus_segment(cluster, segment_type="j", weighed=weighed)
        result = (len(cluster), 
                  nx.number_of_edges(cluster), 
                  nx.diameter(cluster),
                  nx.density(cluster), 
                  average_eccentricity,
                  aa_consensus,
                  nt_consensus,
                  v_consensus,
                  j_consensus)
        results.append(result)
    return pd.DataFrame(results, columns=properties_list)
        
def calc_cluster_consensus(cluster, seq_type="dna", weighed=False):
    motif_dicts = []
    for node in cluster:
        weight = 1
        if weighed:
            weight = node.size
        if seq_type == "dna":
            seq = node.seq_nt
            if seq == "-":
                return "-"
        else:
            seq = node.seq_aa
        node_motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
        motif_dicts.append(node_motif_dict)
    motif_dict = sum_motif_dicts(motif_dicts)
    consensus_seq = get_consensus_from_motif_dict(motif_dict)
    return consensus_seq

def plot_cluster_logo(cluster, seq_type="prot", weighed=False):
    list_of_clonotypes = []
    seq_types = ["prot", "dna"]
    if seq_type not in seq_types:
        raise ValueError(f"Wrong 'seq_type'! Possible values: {', '.join(seq_types)}")
    for node in cluster:
        if seq_type == "dna":
            seq = node.seq_nt
            if seq == "-":
                raise ValueError(f"Node '{node.id}' does not have specified 'cdr3nt' value. Unable to create Logo for 'dna' seq_type")
        else:
            seq = node.seq_aa
        if weighed:
            weight = node.size
            clone = (seq, weight)
        else:
            clone = (seq)
        list_of_clonotypes.append(clone)
    get_logo_for_list_of_clonotypes(list_of_clonotypes, seq_type, plot=True)


def calc_cluster_consensus_segment(cluster, segment_type="v", weighed=False):
    segments = {}
    for node in cluster:
        if segment_type == "v":
            segment = node.v
        else:
            segment = node.j
        weight = 1
        if weighed:
            weight = node.size
        if segment not in segments:
            segments[segment] = weight
        else:
            segments[segment] += weight
    best_segment = ""
    best_score = 0
    for segment in segments:
        score = segments[segment]
        if segments[segment] > best_score:
            best_score = score
            best_segment = segment
    return best_segment

def weight_function(length):
    return 1/(1+length)

def export_clusters_to_gae(clusters):
    
    # additional_properties=[]
    # for node in clusters[0]:
    #     additional_properties = list(node.additional_properties.keys())
    #     break
    
    node_id = 0
    clone_id_dict = {}
    weights = []
    rows = []
    columns = []
    clone_list = []
    for cluster in clusters:
        for node in cluster:
            # add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
            # clone_list.append((str(node), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.size, *add_properties_values))
            clone_list.append((node_id, node.v, node.j, node.seq_aa, node.seq_nt, node.sample_id, node.size))
            clone_id_dict[node.id] = node_id
            weights.append(1)
            rows.append(node_id)
            columns.append(node_id)
            node_id += 1

        attributes = nx.get_edge_attributes(cluster,'length')
        for u,v in cluster.edges():
            length = attributes[(u,v)]
            weight = weight_function(length)
            u_id = clone_id_dict[u.id]
            v_id = clone_id_dict[v.id]
            weights.append(weight)
            weights.append(weight)
            rows.append(u_id)
            rows.append(v_id)
            columns.append(v_id)
            columns.append(u_id)

    col_names = ["node_id", "v", "j", "cdr3aa", "cdr3nt", "sample_id", "size"]
    # col_names += additional_properties
    clonoset = pd.DataFrame(clone_list, columns=col_names)
    clone_count = len(clonoset)

    # create sparse adjacency matrix with floats
    adjacency_matrix = csr_matrix((np.array(weights), (np.array(rows), np.array(columns))),
                          shape = (clone_count, clone_count), 
                          dtype = float).toarray()
    return (adjacency_matrix, clonoset, weights, rows, columns)