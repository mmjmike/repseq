import pandas as pd
from .common_functions import run_parallel_calculation
import networkx as nx


class Node:
    
    def __init__(self, node_id, seq_nt, seq_aa, v, j, sample_id, size=1):
        self.id = node_id
        self.v = v
        self.j = j
        self.seq_aa = seq_aa
        self.seq_nt = seq_nt
        self.sample_id = sample_id
        self.group = None
        self.size = size
        self.additional_properties = {}
        
    def get_group_from_metadata_df(self, sample_metadata_df):
        # obtain experimental group from 
        sub_df = sample_metadata_df.loc[sample_metadata_df["sample_id"] == self.sample_id]
        group = sub_df.iloc[0]["experimental_group"]
        self.group = group
        
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


def find_nodes_and_edges(clonoset_input, sample_metadata, mismatches=1, overlap_type="aa", alice=False):
    if isinstance(clonoset_input, str):
        clonoset=pd.read_csv(clonoset_input,sep="\t")
    else:
        clonoset = clonoset_input
        
    possible_overlap_types = ["aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ"]
    if overlap_type not in possible_overlap_types:
        print("Incorrect overlap type. Possible values: {}".format(", ".join(possible_overlap_types)))    
        return None
    if overlap_type[0:2] == "aa":
        aa = True
    check_v = False
    if "V" in overlap_type:
        check_v = True
    check_j = False
    if "J" in overlap_type:
        check_j = True
    
    v_column = "v"
    j_column = "j"
    aa_column = "cdr3aa"
    nt_column = "cdr3nt"
    if alice:
        v_column = "bestVGene"
        j_column = "bestJGene"
        aa_column = "CDR3.amino.acid.sequence"
        nt_column = "CDR3.nucleotide.sequence"
    
    clonoset = clonoset.merge(sample_metadata[["sample_id", "experimental_group"]])
    nodes_by_len={}
    list_of_all_nodes = []
    for index, row in clonoset.iterrows():
        v = row[v_column]
        j = row[j_column]
        cdr3aa = row[aa_column]
        cdr3nt = row[nt_column]
        freq = row["freq"]
        sample_id = row["sample_id"]
        group = row["experimental_group"]
        len_cdr3aa = len(cdr3aa)
        if len_cdr3aa not in nodes_by_len:
            nodes_by_len[len_cdr3aa] = []
        node = Node(index, cdr3nt, cdr3aa, v, j, sample_id, size=freq)
        node.group = group
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
                edges.append((node_1, node_2)) #save unique codes
    return edges

def create_clusters(nodes, edges):
    main_graph = nx.Graph()
    main_graph.add_nodes_from(nodes)

    for edge in edges:
        node1_id = edge[0].id
        node2_id = edge[1].id
        for node in main_graph:
            if node.id == node1_id:
                node1 = node
            if node.id == node2_id:
                node2 = node
        main_graph.add_edge(node1, node2)
    clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
    return clusters

def filter_one_node_clusters(clusters):
    return [c for c in clusters if len(c)>1]

def save_clusters_for_cytoscape(clusters, output_prefix):
    sif_filename = output_prefix + ".sif"
    #properties_filename = output_prefix + ".prop.tsv"
    properties_metadata_filename = output_prefix + ".prop.metadata.tsv"
    edges = []
    nodes = []
    for cluster in clusters:
        for u,v in cluster.edges():
            edges.append("\tneighbour\t".join([str(u),str(v)]))
        if len(cluster) == 1:
            list(cluster.nodes())[0].id
            edges.append(str(list(cluster.nodes())[0]))
        for node in cluster:
            nodes.append((str(node), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.group, node.size))
    with open(sif_filename, "w") as f:
        f.write("\n".join(edges))
    print("Saved edges to: {}".format(sif_filename))
    
    properties_df = pd.DataFrame(nodes, columns=["code", "cdr3aa", "v", "j", "cdr3nt", "sample_id", "experimental_group", "size"])
    properties_df.to_csv(properties_metadata_filename, index=False, sep="\t")
    print("Saved node properties and metadata to: {}".format(properties_metadata_filename))

def add_alice_hits_to_clusters(clusters, alice_hits_df, check_samples=True):
    if check_samples:
        samples = list(alice_hits_df["sample_id"].unique())
        alice_hits_dict = {}
        for sample in samples:
            hits = set()
            for index, row in alice_hits_df.loc[alice_hits_df["sample_id"] == sample].iterrows():
                v = row["bestVGene"]
                j = row["bestJGene"]
                cdr3nt = row["CDR3.nucleotide.sequence"]
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
            cdr3nt = row["CDR3.nucleotide.sequence"]
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