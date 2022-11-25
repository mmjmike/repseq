import pandas as pd
from .common_functions import run_parallel_calculation, combine_metadata_from_folders, print_progress_bar
from .logo import create_motif_dict, sum_motif_dicts, get_consensus_from_motif_dict
import networkx as nx
import numpy as np


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

        
def find_nodes_and_edges(clonoset_input, mismatches=1, overlap_type="aaV", igh=False):
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
                                        "cloneCount": "count"})
    
    clonoset["v"] = clonoset["v"].apply(lambda x: x.split("*")[0])
    clonoset["j"] = clonoset["j"].apply(lambda x: x.split("*")[0])
    if igh:
        def _split_c(c_segm):
            try:
                return c_segm.split("*")[0]
            except AttributeError:
                return "None"
        clonoset["c"] = clonoset["c"].apply(lambda x: _split_c(x))

    

    nodes_by_len={}
    list_of_all_nodes = []

    size_column_name = "freq"
    if size_column_name not in clonoset.columns:
        size_column_name = "count"
        if size_column_name not in clonoset.columns:
            clonoset["count"] = 1

    for index, row in clonoset.iterrows():
        v = row["v"]
        j = row["j"]
        cdr3aa = row["cdr3aa"]
        cdr3nt = row["cdr3nt"]
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
                edges.append((node_1, node_2)) #save unique codes
    return edges

def create_clusters(clonoset_input, mismatches=1, overlap_type="aaV", igh=False):
    nodes, edges = find_nodes_and_edges(clonoset_input, mismatches=mismatches, overlap_type=overlap_type, igh=igh)
    
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
        main_graph.add_edge(node1, node2)
        edges_done += 1
        print_progress_bar(edges_done, edges_total, program_name=program_name, object_name="edge(s)")
    clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
    total_clusters = len(clusters)
    cluster_num = len(filter_one_node_clusters(clusters))
    singletons = total_clusters - cluster_num
    print(f"Found {cluster_num} clusters (2 or more nodes) and {singletons} single nodes. Total: {total_clusters}")
    return clusters

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
        else:
            seq = node.seq_aa
        node_motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
        motif_dicts.append(node_motif_dict)
    motif_dict = sum_motif_dicts(motif_dicts)
    consensus_seq = get_consensus_from_motif_dict(motif_dict)
    return consensus_seq

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