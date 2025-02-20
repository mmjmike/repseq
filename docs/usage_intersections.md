# Usage: intersections between clonosets

To see the details on all possible parameters in functions, check the [Intersections](functions.md#intersections) module.
    
```py
from repseq import intersections
from repseq import clone_filter as clf
from repseq import clustering

# This functions performs pairwise clonotype overlapping for all clonosets.
# Possible overlap types are [aa, aaV, aaVJ, nt, ntV, ntVJ], aa/nt stands for an amino acid or nucleotide sequence, and V/J/VJ denote a segment type. 
# An output table contains clonotype sequence and V/J segments if required, overlapping clonotypes for each pair, and clonosets they belong to.
downsample_filter = clf.Filter(functionality="f", downsample=15000, by_umi=True, seed=100)
intersect_df = intersections.intersect_clones_in_samples_batch(clonosets_df, cl_filter=downsample_filter, overlap_type="aaV", by_freq=True)

# Calculate overlap distances between clonosets. F, F2 or C metric can be used. The mismatches option specifies the maximum number of mismatches allowed for clonotypes to be considered similar.
# F2 - sum of sqrt of product of similar clonotype frequencies in two clonosets. 
# F - sqrt of the sum of frequency products. 
# C - total frequency of clonotypes in sample1 that are similar to clonotypes in sample2
f2_ntVJ = intersections.overlap_distances(clonosets, cl_filter=downsample_filter, overlap_type="ntVJ", mismatches=0, metric="F2")
f_cd4_aaV = intersections.overlap_distances(clonosets_df.query("subset=='nCD4'"), cl_filter=downsample_filter, overlap_type="aaV", mismatches=0, metric="F")

# Create a table containing the number of times each clonotype appears in each clonoset in clonosets_df
count_table = intersections.count_table(clonosets, cl_filter=downsample_filter, overlap_type="aaV", mismatches=0)

# Create a counnt table for clusters (are provided by the user). Clusters can be created with `create_clusters` function from clustering module. 
# cluster_list = clustering.create_clusters(clonosets, cl_filter=top_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=True)
count_table_by_cluster = intersections.count_table_by_cluster(clonosets_df, clusters_list, cl_filter=downsample_filter, overlap_type="aaV", mismatches=1)

# tcrnet allows to compare two clonosets, typically an experimental dataset and a control one
clonoset_df_exp = ...
clonosets_df_control = ...
tcrnet_compared_clns = intersections.tcrnet(clonosets_df_exp, clonosets_df_control, cl_filter=downsampling, overlap_type="aaVJ", mismatches=1)

```