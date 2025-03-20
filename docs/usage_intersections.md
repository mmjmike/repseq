# Usage: intersections between clonosets

To see the details on all possible parameters in functions, check the [Intersections](functions.md#intersections) module.
    


`intersect_clones_in_samples_batch` function performs pairwise clonotype overlapping for all clonosets.
<br>Possible overlap types are [aa, aaV, aaVJ, nt, ntV, ntVJ], aa/nt stands for an amino acid or nucleotide sequence, and V/J/VJ denote a segment type. 
An output table contains clonotype sequence and V/J segments if required, overlapping clonotypes for each pair, and clonosets they belong to.    

```py
from repseq import intersections
from repseq import clone_filter as clf
from repseq import clustering

downsample_filter = clf.Filter(functionality="f", downsample=15000, by_umi=True, seed=100)
intersect_df = intersections.intersect_clones_in_samples_batch(clonosets_df, cl_filter=downsample_filter, overlap_type="aaV", by_freq=True)
```

|    | cdr3aa           | v        |   sample1_count |   sample2_count | sample1            | sample2            | pair                                     |
|---:|:-----------------|:---------|----------------:|----------------:|:-------------------|:-------------------|:-----------------------------------------|
|  0 | CASSLGQVNTEAFF   | TRBV12-3 |     6.66667e-05 |     0           | sample1_nCD4_1_TRB | sample2_nCD4_1_TRB | sample1_nCD4_1_TRB_vs_sample2_nCD4_1_TRB |
|  1 | CSARDPASGRVDTQYF | TRBV20-1 |     0           |     6.66667e-05 | sample1_nCD4_1_TRB | sample2_nCD4_1_TRB | sample1_nCD4_1_TRB_vs_sample2_nCD4_1_TRB |
|  2 | CASSPKQGNPYEQYF  | TRBV18   |     0           |     6.66667e-05 | sample1_nCD4_1_TRB | sample2_nCD4_1_TRB | sample1_nCD4_1_TRB_vs_sample2_nCD4_1_TRB |
|  3 | CASSWNPTGGTEAFF  | TRBV5-6  |     0           |     6.66667e-05 | sample1_nCD4_1_TRB | sample2_nCD4_1_TRB | sample1_nCD4_1_TRB_vs_sample2_nCD4_1_TRB |
|  4 | CASSLLAGGTDTQYF  | TRBV7-2  |     0           |     6.66667e-05 | sample1_nCD4_1_TRB | sample2_nCD4_1_TRB | sample1_nCD4_1_TRB_vs_sample2_nCD4_1_TRB |


Calculate overlap distances between clonosets. F, F2 or C metric can be used. The mismatches option specifies the maximum number of mismatches allowed for clonotypes to be considered similar.

* F2 - sum of sqrt of product of similar clonotype frequencies in two clonosets. 
* F - sqrt of the sum of frequency products. 
* C - total frequency of clonotypes in sample1 that are similar to clonotypes in sample2

```py
f2_ntVJ = intersections.overlap_distances(clonosets, cl_filter=downsample_filter, overlap_type="ntVJ", mismatches=0, metric="F2")
f_cd4_aaV = intersections.overlap_distances(clonosets_df.query("subset=='nCD4'"), cl_filter=downsample_filter, overlap_type="aaV", mismatches=0, metric="F")
```

Create a table containing the number of times each clonotype appears in each clonoset in clonosets_df.

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

Create a counnt table for clusters (are provided by the user). Clusters can be created with `create_clusters` function from clustering module. 

```py
# cluster_list = clustering.create_clusters(clonosets, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
count_table_by_cluster = intersections.count_table_by_cluster(clonosets_df, clusters_list, cl_filter=downsample_filter, overlap_type="aaV", mismatches=1)
```

TCRnet compares two datasets with their respective clonosets, typically an experimental dataset and a control one. 

```py
clonoset_df_exp = ...
clonosets_df_control = ...
tcrnet_compared_clns = intersections.tcrnet(clonosets_df_exp, clonosets_df_control, cl_filter=downsampling, overlap_type="aaVJ", mismatches=1)
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
            scale_x_log10() +
            scale_y_log10()
    ```
