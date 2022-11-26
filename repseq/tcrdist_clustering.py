import os
import json
import pandas as pd
import networkx as nx

from tcrdist.rep_funcs import compute_pw_sparse_out_of_memory, compute_n_tally_out_of_memory
from tcrdist.repertoire import TCRrep
from .slurm import run_slurm_command_from_jupyter
from .clustering import Node, save_clusters_for_cytoscape
from .io import save_dill_dump, read_dill_dump


REPSEQ_PATH = os.path.join(os.path.expanduser("~"), "soft", "repseq")


def build_tcr_dist_clusters_slurm(clonoset_filename, radius, output_prefix,
                                  chain="beta", species="human", group_colname="group",
                                  cpus=1, time_estimate=4, memory=100, append_command=None):
    # print(os.path.abspath(a_module.__file__))
    
    script_path = os.path.join(REPSEQ_PATH, "repseq", "tcrdist_clusters_slurm.py")
    log_filename = f"{output_prefix}.log"
    command = f"python {script_path} {clonoset_filename} {radius} {output_prefix} --chain {chain} --species {species} --cpus {cpus} --group_colname {group_colname}"
    if isinstance(append_command, str):
        command += f"; {append_command}"
    jn = os.path.basename(output_prefix)
    jobname = f"TCRdist_clusters_{jn}"

    print(f"Running slurm command: {command}")
    run_slurm_command_from_jupyter(command, jobname, cpus, time_estimate, memory, log_filename=log_filename)


def build_tcr_dist_clusters(clonoset_filename, radius, output_prefix, chain="beta", species="human", cpus=1, group_colname="group"):
    
    # Create TCRdist repertoire object
    print("Creating tcrdist repertoire object")
    rep = create_tcr_dist_rep_from_file(clonoset_filename, chain, species)
    print("Tcrdist repertoire object created")


    clone_df_filename = f"{output_prefix}.clone_df.tsv"
    rep.clone_df.to_csv(clone_df_filename, sep="\t", index=False)
    clones = len(rep.clone_df)
    print(f"Clonoset DataFrame ({clones} clones) written to: {clone_df_filename}")
    
    # Calculate TCR distance matrix
    print(f"Computing sparse distances with radius {radius}...")
    sparse_matrix, fragments = compute_pw_sparse_out_of_memory(tr = rep,
                    row_size      = 670,
                    pm_processes  = cpus,
                    pm_pbar       = True,
                    max_distance  = radius,
                    matrix_name   = f'rw_{chain}',
                    reassemble    = True,
                    cleanup       = False)
    
    sparse_matrix_filename = f"{output_prefix}.sparse.matrix"
    save_dill_dump(sparse_matrix, sparse_matrix_filename)
    print(f"Sparse distance matrix written to: {sparse_matrix_filename}")
    
    matrix_fragments_filename = f"{output_prefix}.fragments.matrix"
    save_dill_dump(fragments, matrix_fragments_filename)
    print(f"Sparse matrix fragments written to: {matrix_fragments_filename}")
    
    
    # Calculate neighbours table
    print(f"Computing neighbors with radius {radius}...")
    nhood_df = compute_n_tally_out_of_memory(fragments,
                    matrix_name  = f"rw_{chain}",
                    pm_processes = cpus,
                    to_memory    = True,
                    x_cols       = [group_colname],
                    knn_radius   = radius)
    nhood_df_filename = f"{output_prefix}.neighbours_df.tsv"
    nhood_df.to_csv(nhood_df_filename, sep="\t", index=False)
    print(f"Neighbours DataFrame written to: {nhood_df_filename}")
    nhood_df_read = pd.read_csv(nhood_df_filename, sep="\t")
    
    # create list of NetworkX clusters
    clusters = create_tcr_dist_clusters(rep.clone_df, sparse_matrix, nhood_df_read)
    clusters_filename = f"{output_prefix}.clusters"
    save_dill_dump(clusters, clusters_filename)
    print(f"Clusters list written to: {clusters_filename}")
    
    save_clusters_for_cytoscape([c for c in clusters if len(c) > 1], output_prefix)

    return clusters


def save_pooled_mixcr_clonoset_for_tcr_dist(clonoset, output_filename, chain="beta", count_by_umi=False):
    if chain == "beta":
        cdr3_aa_colname = "cdr3_b_aa"
        cdr3_nt_colname = "cdr3_b_nucseq"
        v_colname = "v_b_gene"
        j_colname = "j_b_gene"
    elif chain == "alpha":
        cdr3_aa_colname = "cdr3_a_aa"
        cdr3_nt_colname = "cdr3_a_nucseq"
        v_colname = "v_a_gene"
        j_colname = "j_a_gene"
    clonoset = clonoset.rename(columns={"aaSeqCDR3": cdr3_aa_colname,
                                    "nSeqCDR3": cdr3_nt_colname,
                                    "allVHitsWithScore": v_colname,
                                    "allJHitsWithScore": j_colname})
    clonoset[v_colname] = clonoset[v_colname].apply(lambda x: x.split("*")[0] + "*01")
    clonoset[j_colname] = clonoset[j_colname].apply(lambda x: x.split("*")[0] + "*01")
    if count_by_umi:
        clonoset["count"] = clonoset["uniqueUMICount"]
    else:
        clonoset["count"] = clonoset["readCount"]
    
    clonoset = clonoset[['count', 'sample_id', cdr3_aa_colname, cdr3_nt_colname, v_colname, j_colname]]

    clonoset.to_csv(output_filename, sep="\t", index=False)
    print(f"Saved pooled clonoset to {output_filename}")


def create_tcr_dist_rep_from_file(clonoset_filename, chain, species):
    clonoset = pd.read_csv(clonoset_filename, sep="\t")
    rep = TCRrep(
        cell_df = clonoset, 
        chains = [chain],
        organism = species,
        db_file = 'alphabeta_gammadelta_db.tsv',
        compute_distances=False)
    return rep


def create_tcr_dist_clusters(clone_df, dist_matrix, nhood_df):
    nodes = create_cluster_nodes_from_tcrdist_df(clone_df)
    print("Nodes total: {}".format(len(nodes)))
    main_graph = nx.Graph()
    main_graph.add_nodes_from(nodes)
    # main_graph = add_edges_from_tcrdist_neighbor_df(main_graph, dist_matrix, nhood_df)
    add_edges_from_tcrdist_neighbor_df(main_graph, dist_matrix, nhood_df)
    clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
    non_single_clusters = len([c for c in clusters if len(c) > 1])
    print("Clusters: {}; single nodes: {}".format(non_single_clusters, len(clusters)-non_single_clusters))
    return clusters.sort(key=lambda x: len(x), reverse=True)


def create_cluster_nodes_from_tcrdist_df(df):
    nodes = []
    for i, r in df.iterrows():
        try:
            v = r["v_b_gene"].split("*")[0]
            j = r["j_b_gene"].split("*")[0]
        except KeyError:
            v = r["v_a_gene"].split("*")[0]
            j = r["j_a_gene"].split("*")[0]
        node = Node(str(i), r["cdr3_b_nucseq"], r["cdr3_b_aa"], v, j, r["sample_id"], r["count"])
        nodes.append(node)
    return nodes

   
def add_edges_from_tcrdist_neighbor_df(graph, matrix, df):
    node_id_dict = {node.id: node for node in graph}
#     print(node_id_dict)
    for i, r in df.iterrows():
        node1_id = str(i)
        node1 = node_id_dict[node1_id]
        neighbors = r.neighbors
        if isinstance(neighbors, str):
            neighbors = json.loads(neighbors)
        for n in neighbors:
            if n != i:
                node2 = node_id_dict[str(n)]
                length = matrix[i,n]
                if length == -1:
                    length = 0
                graph.add_edge(node1, node2, length=length)
    # return graph