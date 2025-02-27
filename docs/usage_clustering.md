# Usage: clustering

[Clustering](functions.md#clustering) find clusters in given clonosets

```py
from repseq import clustering

clusters = clustering.create_clusters(clonosets_df, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
```

Save clusters in Cytoscape format. Filter for clusters of size > 1 (does not include single nodes). Metadata is merged to the node metadata by sample_id column. To read metadata, one could use `read_yaml_metadata` function from [io module](functions.md#io).

```py
clusters_output_prefix = os.path.join(output_dir, "clusters")
# here, one-node clusters are filtered out
clusters_filtered = clustering.filter_one_node_clusters(clusters)
clustering.save_clusters_for_cytoscape(clusters_filtered, clusters_output_prefix, sample_metadata=metadata)
```

Metadata could also be added to node properties prior to saving to Cytoscape. The info will be added to `node.additional_properties` dictionary.
```py
clustering.add_metadata(clusters, metadata)
```

`cluster_properties` outputs a DataFrame with the following columns: ["cluster_id", "nodes", "edges", "diameter", "density", "eccentricity", "concensus_cdr3aa", "concensus_cdr3nt", "concensus_v", "concensus_j"]. To get cluster properties and save the output, use:

```py
cluster_properties = clustering.cluster_properties(clusters_filtered, weighed=True)
cluster_properties.to_csv('clusters.tsv', sep='\t')
```

To create a counnt table for the clusters, use `count_table_by_cluster`. 

```py
count_table_by_cluster = intersections.count_table_by_cluster(clonosets_df, clusters_list, cl_filter=downsample_filter, overlap_type="aaV", mismatches=1)
```

To visualize cluster's CDR3 consensus sequence, use:

```py
clustering.plot_cluster_logo(clusters[0])
```