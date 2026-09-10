# Usage: intersections between clonosets

## Directional clonotype similarity

`intersections.similarity` calculates the directional S-metric for exact or
mismatch-tolerant clonotype comparisons. Rows are target samples and columns
are comparison samples. A target clonotype contributes once to a matrix value
when it is similar to at least one comparison clonotype, even if several
comparison clonotypes match it.

```py
from repseq import intersections

similarity_frequency = intersections.similarity(
    clonosets_df,
    cl_filter=downsample_filter,
    overlap_type="aaV",
    mismatches=1,
    result="freq",
    cpu=4,
)
```

Similarity definitions are:

- `aa` or `nt`: compare CDR3 amino-acid or nucleotide sequences;
- `aaV` or `ntV`: additionally require the same V segment;
- `aaVJ` or `ntVJ`: additionally require the same V and J segments;
- `VJ`: compare V/J combinations and ignore CDR3 sequence;
- `VJlen`: compare V/J combinations and amino-acid CDR3 length.

For sequence-based overlap types, `mismatches` is the maximum CDR3 Hamming
distance. It is ignored for `VJ` and `VJlen`.

The `result` argument controls the output:

- `freq`: total target-sample frequency of clonotypes with at least one match;
- `count`: total target-sample count of clonotypes with at least one match;
- `number`: number of distinct target clonotypes with at least one match;
- `table`: every matching target/comparison clonotype pair, including counts,
  frequencies, Hamming distance, sample identifiers, and pair identifier.

```py
similarity_counts = intersections.similarity(
    clonosets_df,
    overlap_type="ntVJ",
    mismatches=2,
    result="count",
)

similarity_pairs = intersections.similarity(
    clonosets_df,
    overlap_type="aa",
    mismatches=1,
    result="table",
)
```

Pass `clonosets_df2` and optionally `cl_filter2` for a rectangular comparison.
The first dataframe remains the target set in matrix rows, and the second
becomes the comparison set in columns.

To see further details, check the [Intersections](functions.md#intersections) module.

## Clonoset intersection

`intersect_clones_in_samples_batch` function performs pairwise clonotype overlapping for all clonosets.
<br>Possible overlap types are [aa, aaV, aaVJ, nt, ntV, ntVJ, VJ, VJlen]. aa/nt selects an amino acid or nucleotide sequence; VJ uses a segment pair and VJlen additionally uses amino-acid CDR3 length.
An output table contains clonotype sequence and V/J segments if required, overlapping clonotypes for each pair, and clonosets they belong to.    

!!! tip "clonosets_df and clonosets_df2"
    If `clonosets_df2` is None (default), samples within `clonosets_df` are compared with each other; 
    Otherwise, the comparison is performed exclusively between samples from `clonosets_df` and `clonosets_df2`. 
    <br>If `cl_filter2` is set, it is applied to clonosets in `clonosets_df2`. If there are samples with non-unique sample_ids between the two dataframes, both filters will be applied to those samples.


```py
from repseq import intersections
from repseq import clone_filter as clf
from repseq import clustering

downsample_filter = clf.Filter(functionality="f", downsample=15000, by_umi=True, seed=100)
intersect_df = intersections.intersect_clones_in_samples_batch(clonosets_df, cl_filter=downsample_filter, overlap_type="aaV")
```

|    | cdr3aa           | v        | sample1_count | sample2_count | sample1_freq | sample2_freq | sample1 | sample2 | pair |
|---:|:-----------------|:---------|--------------:|--------------:|-------------:|-------------:|:--------|:--------|:-----|
|  0 | CASSLGQVNTEAFF   | TRBV12-3 |            10 |             0 |     0.666667 |            0 | sample1 | sample2 | sample1_vs_sample2 |
|  1 | CSARDPASGRVDTQYF | TRBV20-1 |             5 |            12 |     0.333333 |            1 | sample1 | sample2 | sample1_vs_sample2 |

<br>

## Overlap distances between clonosets

Calculate overlap distances between clonosets. F, F2, C, J, BC or JCD [metric](https://mixcr.com/mixcr/reference/mixcr-postanalysis/?h=pairwise#pairwise-distance-metrics) can be used. The mismatches option specifies the maximum number of mismatches allowed for clonotypes to be considered similar. 

* F2 - clonotype-wise sum of geometric mean frequencies
* F -  geometric mean of relative overlap frequencies
* C - total frequency of clonotypes in sample1 that are similar to clonotypes in sample2
* BC ([Bray-Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity)) - sum of differences between clonotype frequencies or counts (`by_freq`=False) in sample1 and sample2 divided by the total counts in sample1 and sample2  
* J ([Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)) - size of sample1 and sample2 intersection divided by the size of their union
* JCD ([Jensen-Shannon divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence))

```py
f2_ntVJ = intersections.overlap_distances(clonosets, cl_filter=downsample_filter, overlap_type="ntVJ", mismatches=0, metric="F2")
f_cd4_aaV = intersections.overlap_distances(clonosets_df.query("subset=='nCD4'"), cl_filter=downsample_filter, overlap_type="aaV", mismatches=0, metric="F")
```

<br>

## `count_table`

Create a table containing the number of times each clonotype appears in each clonoset in `clonosets_df`. 
<br>For `overlap_type`, possible overlap types are [aa, aaV, aaVJ, nt, ntV, ntVJ], aa/nt stands for an amino acid or nucleotide sequence, and V/J/VJ denote a segment type. 

```py
count_table = intersections.count_table(clonosets, cl_filter=downsample_filter, overlap_type="aaV", mismatches=0)
```

|                                |   sample1_nCD4_1_TRB |   sample1_nCD8_1_TRB |   sample1_nTreg_1_TRB |   sample2_nCD4_1_TRB |   sample2_nCD8_1_TRB |   sample2_nTreg_1_TRB |
|:-------------------------------|---------------------:|---------------------:|----------------------:|---------------------:|---------------------:|----------------------:|
| ('CASSLGQVNTEAFF', 'TRBV12-3') |                    1 |                    0 |                     0 |                    0 |                    0 |                     0 |
| ('CASSPKQGNPYEQYF', 'TRBV18')  |                    0 |                    1 |                     0 |                    0 |                    0 |                     0 |
| ('CASSLLAGGTDTQYF', 'TRBV7-2') |                    0 |                    1 |                     0 |                    1 |                    1 |                     0 |
| ('CASSHGEGTQYF', 'TRBV3-1')    |                    2 |                    0 |                     0 |                    0 |                    0 |                     0 |
| ('CASSDREGYTEAFF', 'TRBV6-5')  |                    0 |                    1 |                     0 |                    0 |                    0 |                     0 |


## `count_table_by_cluster`

Create a count table for clusters as opposed to single clonotypes (clusters are provided by the user). They can be created with `create_clusters` function from [clustering](functions.md#clustering) module. 

```py
clusters = clustering.create_clusters(clonosets, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
clusters_filtered = clustering.filter_one_node_clusters(clusters)
```

```py
count_table_by_cluster = intersections.count_table_by_cluster(clonosets_df, clusters_list, cl_filter=downsample_filter, overlap_type="aaV", mismatches=1)
```

|    | feature_id   |   sample1_nCD4_1_TRB |   sample1_nCD8_1_TRB |   sample2_nCD4_1_TRB |   sample2_nCD8_1_TRB |
|---:|:-------------|---------------------:|---------------------:|---------------------:|---------------------:|
|  0 | cluster_0    |          0.000133333 |          0.0008      |          0.000866667 |          0.0014      |
|  1 | cluster_1    |          0.000333333 |          0.000333333 |          0.000666667 |          0.000866667 |
|  2 | cluster_2    |          0.000133333 |          0.000333333 |          0.000666667 |          0.000666667 |
|  3 | cluster_3    |          6.66667e-05 |          0.000333333 |          0.0008      |          0.0008      |
|  4 | cluster_4    |          0.000333333 |          0.000133333 |          0.00106667  |          6.66667e-05 |

<br>

## TCRnet

TCRnet compares two datasets with their respective clonosets, typically an experimental dataset and a control one. 

```py
clonoset_df_exp = ...
clonoset_df_control = ...
tcrnet_compared_clns = intersections.tcrnet(clonosets_df_exp, clonoset_df_control, cl_filter=downsampling, overlap_type="aaVJ", mismatches=1)
``` 

|    | clone                                      |   count_exp |   count_control |   group_count_exp |   group_count_control |    fold |   p_value_b |   p_value_p |   p_value_b_adj |   p_value_p_adj |   log10_b_adj |   log10_p_adj |   log2_fc |
|---:|:-------------------------------------------|------------:|----------------:|------------------:|----------------------:|--------:|------------:|------------:|----------------:|----------------:|--------------:|--------------:|----------:|
|  0 | ('CASSPGVGFVEKLFF', 'TRBV11-2', 'TRBJ1-4') |           1 |               1 |                 7 |                    29 | 4.14286 |    0.211254 |     0.20811 |        0.250473 |        0.247205 |      0.601239 |      0.606943 |   2.05063 |
|  1 | ('CASSLMKTENEKLFF', 'TRBV11-2', 'TRBJ1-4') |           1 |               0 |                 7 |                    29 | 8.28571 |    0        |     0       |        0        |        0        |    inf        |    inf        |   3.05063 |
|  2 | ('CASSLGGHPNEKLFF', 'TRBV11-2', 'TRBJ1-4') |           1 |               0 |                 7 |                    29 | 8.28571 |    0        |     0       |        0        |        0        |    inf        |    inf        |   3.05063 |


Plot the TCRnet result as a volcano plot. Points are blue only when they pass both the adjusted p-value and fold-change thresholds; the other three pass states use different grey shades. Infinite significance values are plotted at one unit above the largest finite value.

```py
from repseq import plot as rsplot

fig = rsplot.tcrnet_volcano(
    tcrnet_compared_clns,
    y="log10_b_adj",  # or "log10_p_adj"
    p_threshold=0.05,
    log2_fc_threshold=1,
)
```


??? info "Visualization"
    Properties from proc_table can be visualized in Jupyter notebook using %%R cell magic. 
    ![intersections](images_docs/intersections_table.png)
    
    ```py
    %%R -i intersect_df -h 600 -w 700
    intersect_df %>% 
        ggplot(aes(x=sample1_count, y=sample2_count)) +
            geom_point()+
            theme_bw()+
            facet_wrap(vars(pair)) +
            scale_x_log10(limits = c(1e-5, 3.5e-03)) +
            scale_y_log10(limits = c(1e-5, 3.5e-03))
    ```


## Beta diversity

The `repseq.beta` module calculates many metrics from one count-first full
intersection table. See the dedicated [Beta diversity usage page](usage_beta.md)
for the complete metric list, formulas, normalization rules, and examples.
