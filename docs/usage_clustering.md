# Usage: clustering

[Clustering](functions.md#clustering) finds clusters in given clonosets.


## How to create clusters

`create_clusters` function args:

* `cl_filter`: see [the previous page](usage_stats.md#filtering-clonosets) and [Filter](functions.md#-clone_filter) for further explanations.
* `mismatches`: specifies the maximum number of mismatches allowed for clonotypes to qualify as neighbours (adjacent).
* `overlap_type`: Possible overlap types are [`aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`], `aa`/`nt` stands for an amino acid or nucleotide sequence, and `V`/`J`/`VJ` denote a segment type.
* If `igh`=True, the constant (C) segment is kept.
* If `tcrdist_radius` is not None, only edges between clones with a tsrdist metric less than or equal to the specified radius are built. It overrides `overlap_type` and `mismatches`. A TCRdist modification is introduced: no gaps are allowed in CDR3, weights are in 3:1 ratio for CDR3 compared to other regions. Also, V-segment distances are pre-calculated and currently are available only for <i>Homo sapiens</i>.
* If `by_freq`=True (default), clonotype frequencies are used instead of counts.


```py
from repseq import clustering

clusters = clustering.create_clusters(clonosets_df, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
```

Output is a list of NetworkX Graph() objects — separate clusters and single nodes. The list is sorted by cluster size.

```py
clusters[:3]
[<networkx.classes.graph.Graph at 0x7fa133b1bd00>,
 <networkx.classes.graph.Graph at 0x7fa15200cb80>,
 <networkx.classes.graph.Graph at 0x7fa154aa3c10>]
```

<br>

## Object state and calculation caching

`Clusters` tracks completed workflow stages in `clusters.state` and records the
parameters used for each stage in `clusters.state_parameters`. Printing the
object reports whether clonotypes were read, clusters were created, metadata was
added, node Pgen was calculated, and ALICE was calculated. It also reports the
input source, filter settings, metadata columns, clustering parameters, graph
size, and cluster/singleton counts.

```py
print(clusters)
clusters.state
clusters.state_parameters
```

The expected order is to read clonotypes, create clusters, and then optionally
run ALICE. Clonotype-reading methods and `create_clusters` accept
`verbosity=False` to suppress progress messages; the older `verbose=False`
keyword remains supported. Metadata may be supplied any time after clonotypes
are read; metadata
added before clustering is retained and applied when nodes are created. Methods
that require clusters raise an actionable error if clustering has not yet been
run.

`Clusters.properties(cpu=None)` calculates one cluster per worker with a progress
bar, caches the result after its first calculation, and returns a copy of the
cached dataframe. Use `cpu=1` for sequential calculation. Single-node clusters
reuse their node values directly; V and J consensus calculations are skipped
when the clustering `overlap_type` already guarantees a single V or J value.
Repeating `create_clusters` with unchanged clonotypes
and identical mismatch/TCRdist parameters reuses the existing graph and prints
the previously obtained cluster summary instead of recalculating it. Reading new
clonotypes or changing/filtering the cluster collection invalidates applicable
caches and downstream ALICE state.

## Plotting cluster networks

Use `Clusters.plot_cluster` to plot one cluster, a list of up to 50
clusters, or every cluster as facets. Explicit selections use persistent integer
`cluster_no` values or string IDs such as `"cluster_235"`, with the same mixed
identifier behavior as `Clusters.select`; they are not positions in the current
Python list. Pass `cluster_no=None` to select every cluster. This automatic mode
defaults to `max_clusters=50`; if the object has
more clusters, no figure is created until you either increase `max_clusters`
or pass a smaller explicit selection. By default, node size uses linear `count`
scaling with `min_size + linear_scale * count`, `min_size=50`, and
`linear_scale=1`. Set `size` to another numeric node property, use `size=None`
for uniform nodes, or set `log_scaled=True` for
`min_size + linear_scale * log2(value + 1)`. The
default spring layout can be replaced with `kamada_kawai`, `circular`, `shell`,
or `spectral`. Each panel title shows the persistent `cluster_id`. `height`
sets each facet's height in inches and `aspect` sets width divided by height; an
explicit `figsize` still overrides these automatic dimensions.

```py
fig = clusters.plot_cluster("cluster_10", label="seq_aa")

fig = clusters.plot_cluster(None)

fig = clusters.plot_cluster(
    ["cluster_10", 11, "cluster_12"],
    layout="kamada_kawai",
    color="sample_id",
    palette={"sample_1": "#4C78A8", "sample_2": "#F58518"},
    label="id",
    shape="group",
    size="count",
    min_size=50,
    log_scaled=False,
    linear_scale=1,
    height=3.5,
    aspect=1.2,
)
```

`color`, `label`, and `shape` may name built-in node attributes such as `id`,
`sample_id`, `v`, or `j`, or values previously added to
`node.additional_properties`. Numeric color properties use a continuous gradient
by default, while boolean, categorical, and string properties use discrete
colors. Set `color_mode="discrete"` to treat numeric values as levels, or
`color_mode="continuous"` to explicitly request a gradient. A palette name or
color sequence controls continuous gradient colors. Missing color values are
displayed in light grey. Shape grouping supports five levels, displayed as
circle, triangle, rhombus, hexagon, and square. Discrete color and shape legends
are placed below the faceted panel; continuous colors use a color bar.

Calculate OLGA generation probabilities before plotting `log10_pgen` as a
gradient:

```py
clusters.calc_pgen(hum_pgen_model, cpu=4)
fig = clusters.plot_cluster(
    "cluster_10",
    color="log10_pgen",
    palette="magma",
)
```

Each cluster is processed as a separate parallel job. Use `cpu=1` for a
sequential calculation or omit `cpu` to use the executor default. Each node
receives `pgen` and `log10_pgen` additional properties, where `log10_pgen` is
`-log10(pgen)`. Nodes with `TRBV21-1` or `TRBV7-5` receive
missing values because these genes are unfamiliar to the standard OLGA human
TRB model. For IGH clusters, constant-gene calls are
stored as the `isotype` node property with labels such as `IgM`, `IgG1`, and
`IgA2`. Non-IGH constant calls, including IGK, IGL, TRA, TRB, TRG, and TRD,
produce `isotype=None` rather than an error. Coloring by `isotype` displays
these values as `NA` and uses the same ordering and palette as
`rsplot.isotype_fraction`.

## Plotting cluster sequence logos

Use `Clusters.plot_logo` to create protein or DNA sequence logos for one or
more clusters. The first argument accepts persistent cluster numbers, `cluster_N`
IDs, or mixed iterables using the same rules as `select`. Logo weights may come
from `count`, `freq`, `nodes`, or any custom non-negative numeric node property.
`nodes` assigns equal weight to every node. Set `plot=False` to return the
normalized motif dataframe; multiple clusters return a dictionary keyed by
`cluster_id`.

```py
clusters.plot_logo("cluster_10")
clusters.plot_logo(10, seq_type="dna", weight="freq")
motif = clusters.plot_logo("cluster_10", weight="custom_weight", plot=False)
motifs = clusters.plot_logo([10, "cluster_11"], plot=False)
```

## Filtering clusters and calculating custom properties

Cluster expressions provide a compact interface for filtering without manually
iterating over every graph. Import the standard metrics and node predicates:

```py
from repseq.clustering import (
    all_nodes,
    any_nodes,
    cluster_size,
    proportion,
    total_count,
    top_clusters,
)
```

Select known clusters directly with `Clusters.select`. Integers are interpreted
as `cluster_no`, while strings such as `"cluster_235"` are interpreted as
`cluster_id`. Ranges, dataframe/Series columns, and mixed identifier lists are
accepted. Requested order is preserved and duplicate identifiers are returned
only once.

```py
selected = clusters.select(["cluster_24", 25])
selected = clusters.select(range(10, 20))
selected = clusters.select(
    cluster_table.loc[cluster_table["cluster_size"] >= 5, "cluster_id"]
)
```

Use `top_clusters(N)` inside `filter` to rank the collection by node count
(descending), total node `count` (descending), and amino-acid consensus
(alphabetically), then retain the first `N` clusters. The returned collection is
in ranking order. `clusters.top_clusters(N)` is an equivalent convenience
method.

```py
largest = clusters.filter(top_clusters(100))
largest = clusters.top_clusters(100)
```

Metrics can be compared directly. `cluster_size` is the number of nodes and
`total_count` is the sum of node `count` values. Filtering returns a new
`Clusters` collection and preserves the original cluster numbers.

```py
large_clusters = clusters.filter(cluster_size >= 10)
high_count_clusters = clusters.filter(total_count >= 10)
```

Use `&`, `|`, and `~` to combine conditions. Parenthesize every comparison
because Python's bitwise operators have different precedence from comparisons.

```py
selected = clusters.filter(
    (cluster_size >= 5)
    & (
        proportion(
            "group_property_name",
            ["group_1", "group_2"],
            weight="count",
        )
        > 0.25
    )
    & (
        total_count(
            "group_property_name",
            "control",
            weight="freq",
        )
        == 0
    )
)
```

`all_nodes` requires every node to have one of the selected values, while
`any_nodes` requires at least one matching node. Properties may be built-in
node attributes or values from `node.additional_properties`. The aliases
`cdr3aa` and `cdr3nt` refer to `node.seq_aa` and `node.seq_nt`; these aliases
are also accepted by `plot_cluster`.

```py
naive_clusters = clusters.filter(all_nodes("isotype", ["IgM", "IgD"]))
specific_clusters = clusters.filter(any_nodes("specificity", "specific"))
ighv1_3_clusters = clusters.filter(any_nodes("v", "IGHV1-3"))
sequence_clusters = clusters.filter(
    any_nodes("cdr3aa", "CARSRKDCSGGSCYSGGFDYW")
)
```

`weight` may be `count`, `freq`, `nodes`, or any non-negative numeric node
property. `nodes` gives every matching node a weight of one. `proportion`
divides the matching weight by the total weight in the cluster. Callable
`total_count(...)` sums the selected nodes using the requested weight.

Use `custom_properties` to evaluate metrics for every cluster. The result always
starts with `cluster_no` and `cluster_id`. Generated aggregate names are
dataframe-friendly; use `.alias(...)` when a shorter name is preferred.

```py
group_1_proportion = proportion(
    "group_property_name",
    "group_1",
    weight="count",
)

cluster_table = clusters.custom_properties([
    cluster_size,
    total_count,
    group_1_proportion.alias("group_1_count_proportion"),
])
```

## Intersecting clusters with new clonosets

`Clusters.intersect_with_clonosets` measures how much of each target clonoset is
similar to each existing cluster. Matching follows the directional rules used by
`intersections.similarity`: `overlap_type` selects sequence and V/J constraints,
and `mismatches` sets the maximum Hamming distance for sequence-based overlap
types. A target clonotype contributes only once to a cluster even if it matches
several nodes in that cluster.

```py
cluster_counts = clusters.intersect_with_clonosets(
    target_clonosets,
    cl_filter=func_filter,
    overlap_type="aaVJ",
    mismatches=1,
    cpu=4,
)
```

The result is a wide cluster count table with `cluster_id`, `consensus`,
`concensus_cdr3aa`, `concensus_v`, and `concensus_j`, followed by one column per
target `sample_id`. By default, sample columns contain the summed filtered target
counts matching each cluster. With `by_freq=True`, they contain matched target
counts divided by the filtered target sample's total count, matching
`intersections.similarity(result="freq")`.

```py
cluster_frequencies = clusters.intersect_with_clonosets(
    target_clonosets,
    overlap_type="aaV",
    mismatches=1,
    by_freq=True,
)
```

## Clusters from a pooled DataFrame

Alternatively, one can create clusters from a dataframe with clonotypes. Mandatory columns are [`freq`, `count`, `v`, `j`, `cdr3aa`, `cdr3nt`, `sample_id`].

`pooled_df` example:

|    |   count |        freq | cdr3nt                                           | cdr3aa           | v        | d     | j       | c     |   VEnd |   DStart |   DEnd |   JStart | sample_id          |
|---:|--------:|------------:|:-------------------------------------------------|:-----------------|:---------|:------|:--------|:------|-------:|---------:|-------:|---------:|:-------------------|
|  0 |     117 | 7.46674e-05 | TGTGCCAGCAGTCGCCACAGTTACAGGGATGGCTACACCTTC       | CASSRHSYRDGYTF   | TRBV12-3 | TRBD1 | TRBJ1-2 | TRBC1 |     11 |       22 |     27 |       28 | sample1_nCD4_1_TRB |
|  1 |     109 | 6.95619e-05 | TGTGCCAGCAGTTTAGCGCATCAGGGAGGCAGCTATGGCTACACCTTC | CASSLAHQGGSYGYTF | TRBV12-4 | TRBD2 | TRBJ1-2 | TRBC1 |     18 |       23 |     28 |       32 | sample1_nCD4_1_TRB |
|  2 |     105 | 6.70092e-05 | TGTGCCAGCAGCCCGGGACTGGCCTACAATGAGCAGTTCTTC       | CASSPGLAYNEQFF   | TRBV12-3 | TRBD2 | TRBJ2-1 | TRBC2 |     10 |       11 |     19 |       22 | sample1_nCD4_1_TRB |


```py
create_clusters_from_pooled_df(pooled_df, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None,
                                   count_by_freq=True, _run_from_create_clusters=False)
```

## Adding metadata to clusters

Metadata, if present, could also be added to node properties prior to saving to Cytoscape. The info will be added to `node.additional_properties` dictionary. Note that the metadata should contain the same `sample_id`s that were used in the `clonosets_df` when creating the clusters. 

Metadata example:

|    | sample_id          |   group | type   |
|---:|:-------------------|--------:|:-------|
|  0 | sample1_nCD4_1_TRB |       1 | nCD4   |
|  1 | sample1_nCD8_1_TRB |       1 | nCD8   |
|  3 | sample2_nCD4_1_TRB |       2 | nCD4   |
|  4 | sample2_nCD8_1_TRB |       2 | nCD8   |


```py
clustering.add_metadata(clusters, metadata)
```

<br>

## Save clusters in Cytoscape format

Here, we filter out single-node clusters. Clusters are exported in two forms: edges are saved in a .sif file and cluster properties are in a tab-separated .csv file.  In the case of TCRdist, the edges are also assigned a length (radius).

```py
clusters_output_prefix = os.path.join(output_dir, "clusters")
# here, one-node clusters are filtered out
clusters_filtered = clustering.filter_one_node_clusters(clusters)
clustering.save_clusters_for_cytoscape(clusters_filtered, clusters_output_prefix, sample_metadata=metadata)
```

<br>

## Cluster properties

!!! tip "node size"
    If `weighed` is set to True, the weight of a node is determined by its size. The size of the node is defined by the `by_freq` parameter in the `create_clusters` function, which indicates whether the size is calculated based on counts or frequencies.

Cluster properties include consensus CDR3, v- and j-segment sequences, as well as some properties of clusters as graphs:

* diameter: the maximum eccentricity in a graph.
* [density](https://networkx.org/documentation/stable/reference/generated/networkx.classes.function.density.html): The density is 0 for a graph without edges and 1 for a complete graph. The density of multigraphs can be higher than 1.
* eccentricity: for a node v, it is the maximum distance from v to all other nodes in a graph. Cluster-wise, it is the average eccentricity of all the nodes within the cluster.

```py
cluster_properties = clustering.cluster_properties(clusters_filtered, weighed=True)
cluster_properties.to_csv('clusters.tsv', sep='\t')
```

|    | cluster_id   |   nodes |   edges |   diameter |   density |   eccentricity | concensus_cdr3aa   | concensus_cdr3nt                        | concensus_v   | concensus_j   |
|---:|:-------------|--------:|--------:|-----------:|----------:|---------------:|:-------------------|:----------------------------------------|:--------------|:--------------|
|  0 | cluster_0    |      58 |     164 |         10 | 0.0992136 |        7.96552 | CASSLTGSYEQYF      | TGCGCCAGCAGCTTGGCAGGGTCCTACGAGCAGTACTTC | TRBV5-1       | TRBJ2-7       |
|  1 | cluster_1    |      47 |     188 |          7 | 0.173913  |        5.40426 | CASSLGGNTEAFF      | TGCGCCAGCAGCTTGGCAGGGAACACTGAAGCTTTCTTT | TRBV5-1       | TRBJ1-1       |
|  2 | cluster_2    |      38 |     154 |          7 | 0.219061  |        5.5     | CASSLDTYEQYF       | TGCGCCAGCAGCTTGGACACCTACGAGCAGTACTTC    | TRBV5-1       | TRBJ2-7       |
|  3 | cluster_3    |      35 |     167 |          6 | 0.280672  |        4.45714 | CASSLSYEQYF        | TGTGCCAGCAGTTTAGCCTACGAGCAGTACTTC       | TRBV12-3      | TRBJ2-7       |
|  4 | cluster_4    |      33 |      79 |          7 | 0.149621  |        5.42424 | CASSLGTDTQYF       | TGCGCCAGCAGCTTGGGCACAGATACGCAGTATTTT    | TRBV5-1       | TRBJ2-3       |


For creating a table with counts or frequencies by cluster, see [the previous page](usage_intersections.md#count_table_by_cluster).

<br>


## Sequence logo

To visualize cluster's CDR3 consensus sequence, use `plot_cluster_logo`. Possible `seq_type` are `prot` and `dna`.

```py
clustering.plot_cluster_logo(clusters[0])
```

![logo](images_docs/logo.png)

```py
clustering.plot_cluster_logo(clusters[0], seq_type='dna', weighed=True)
```
![logo](images_docs/logo_dna.png)

<br>


##  Custom cluster metric example

```py
top_filter = clf.Filter(functionality="f", top=4000, by_umi=True, mix_tails=True, seed=100)
cd4_clusters = clustering.create_clusters(clonosets_df, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
clustering.add_metadata(clonosets_df, metadata)
```

Metadata in this example:

|    | sample_id        | experimental_group   | subset   |
|---:|:-----------------|:---------------------|:---------|
|  0 | UCB4_nCD4_1_TRB  | late                 | nCD4     |
|  4 | UCB11_nCD4_1_TRB | preterm              | nCD4     |
|  7 | UCB10_nCD4_1_TRB | term                 | nCD4     |
|  8 | UCB2_nCD4_1_TRB  | term                 | nCD4     |

This function calculates the total frequency of all clonotypes within a cluster, as well as the percentage of clonotypes from different experimental groups present in the cluster.

```py
def calc_custom_clusters_properties(clusters):
    results = []
    properties = ["nodes", "sum_frequency", "preterm_percent", "term_percent", "late_percent"]
    
    for cluster in clusters:
        
        prop1_nodes = len(cluster)
        prop2_sum_size = sum([node.size for node in cluster])
        prop3_sum_size = calc_cluster_preterm_clones_percent(cluster)
        
        result = [prop1_nodes, prop2_sum_size, *prop3_sum_size]
    
        results.append(result)
    df = pd.DataFrame(results, columns = properties)
    return df

def calc_cluster_clones_percent(cluster):
    total_size = sum([node.size for node in cluster])
    preterm = sum([node.size for node in cluster if node.additional_properties["experimental_group"] == "preterm"])
    term = sum([node.size for node in cluster if node.additional_properties["experimental_group"] == "term"])
    late = sum([node.size for node in cluster if node.additional_properties["experimental_group"] == "late"])
    percent_preterm = round(preterm/total_size*100, 2)
    percent_term = round(term/total_size*100, 2)
    percent_late = round(late/total_size*100, 2)
    return percent_preterm, percent_term, percent_late
```

|    |   nodes |   sum_frequency |   preterm_percent |   term_percent |   late_percent |
|---:|--------:|----------------:|------------------:|---------------:|---------------:|
|  0 |      73 |      0.0210043  |             24.1  |          62.04 |          13.86 |
|  1 |      41 |      0.0115223  |             44.21 |          41.43 |          14.36 |
|  2 |      39 |      0.0114498  |             25.18 |          65.35 |           9.47 |
|  3 |      38 |      0.0119551  |             32.43 |          55.64 |          11.93 |
|  4 |      29 |      0.00817272 |             36.06 |          55.56 |           8.38 |

<br>To see another example of a custom function for stats calculation, visit [stats](usage_stats.md#custom-stats) page.

<br>


## Community detection

Split branchy clusters into independently numbered community clusters with Leiden
(the default) or Louvain:

```py
community_clusters = clusters.find_cluster_communities(
    algorithm="leiden",
    resolution=1,
    seed=1,
)
```

Clusters with fewer than four nodes are copied without community detection. The
returned clusters and nodes store `previous_community_id`, while the source
clusters and nodes store `splitted_community_id` for coloring and comparison.
Calling the method again overwrites `splitted_community_id` on the source object.

<br>

## [ALICE](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000314) (Antigen-specific Lymphocyte Identification by Clustering of Expanded sequences)

ALICE works by treating clonotypes as graph vertices, with edges connecting sequences which differ by at most 1 CDR3 amino acid. It identifies clonotypes with a higher numbers of neighbors than expected by a null model of recombination, separating clusters of antigen-responding clonotypes from clusters arising from recombination statistics.<br>Currently, it is implemented for <i>H.sapiens</i> only.

```py
alice(clusters, overlap_type='aaVJ', mismatches=1, species="hs", olga_warnings=False)
```
