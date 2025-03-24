# Usage: differential expression analysis 

!!! note "Setting up the environment"
    Before getting started, make sure that main_repseq environment is chosen. Otherwise, check the [installation guide](installation.md)

```py
import sys
REPSEQ_PATH = 'your_path_to_repseq'
sys.path.append(REPSEQ_PATH)
from repseq import io as repseqio
from repseq import common_functions as cf
from repseq import mixcr as mx
from repseq import slurm
from repseq import clonosets as cl
from repseq import stats
from repseq import clone_filter as clf
from repseq import intersections
from repseq import clustering
from repseq import diffexp
import os

import pandas as pd
import random
import json
import re
```

## Input data

To run differential expression analysis, two dataframes are required:

* `metadata` dataframe. Note that <b>sample_id</b> (must be unique and match `count_table` column names) and <b>group</b> columns are <b>nessesary</b> for this pipeline. Example:

|    | sample_id         | Subset   |   Replica |   Number_of_cells | group          |
|---:|:------------------|:---------|----------:|------------------:|:---------------|
|  0 | Luk_TRB_AdV_CD4_1 | CD4      |         1 |              1376 | PepTivator_AdV |
|  1 | Luk_TRB_AdV_CD4_2 | CD4      |         2 |              1504 | PepTivator_AdV |
|  2 | Luk_TRB_AdV_CD4_3 | CD4      |         3 |              1623 | PepTivator_AdV |
|  3 | Luk_TRB_AdV_CD8_1 | CD8      |         1 |             37698 | PepTivator_AdV |
|  4 | Luk_TRB_AdV_CD8_2 | CD8      |         2 |             21359 | PepTivator_AdV |
|  5 | Luk_TRB_AdV_CD8_3 | CD8      |         3 |             21031 | PepTivator_AdV |
|  6 | Luk_TRB_CMV_CD4_1 | CD4      |         1 |              1018 | PepTivator_CMV |
|  7 | Luk_TRB_CMV_CD4_2 | CD4      |         2 |              1032 | PepTivator_CMV |
|  8 | Luk_TRB_CMV_CD4_3 | CD4      |         3 |              2788 | PepTivator_CMV |
|  9 | Luk_TRB_CMV_CD8_1 | CD8      |         1 |             19617 | PepTivator_CMV |
| 10 | Luk_TRB_CMV_CD8_2 | CD8      |         2 |             17490 | PepTivator_CMV |
| 11 | Luk_TRB_CMV_CD8_3 | CD8      |         3 |             16662 | PepTivator_CMV |
| 12 | Luk_TRB_EBV_CD4_1 | CD4      |         1 |               453 | PepTivator_EBV |
| 13 | Luk_TRB_EBV_CD4_2 | CD4      |         2 |               840 | PepTivator_EBV |
| 14 | Luk_TRB_EBV_CD4_3 | CD4      |         3 |               770 | PepTivator_EBV |
| 15 | Luk_TRB_EBV_CD8_1 | CD8      |         1 |             41274 | PepTivator_EBV |
| 16 | Luk_TRB_EBV_CD8_2 | CD8      |         2 |             22020 | PepTivator_EBV |
| 17 | Luk_TRB_EBV_CD8_3 | CD8      |         3 |             39014 | PepTivator_EBV |


*  `count_table` which can be created with [Intersections](usage_intersections.md#count_table) `count_table` (for clonotypes of [clonotype groups](usage_diffexp.md#for-clonotype-groups)) or `count_table_by_cluster` (for [clusters](usage_diffexp.md#count-table-by-cluster)) functions. Example:

|    | feature_id   |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |
|---:|:-------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
|  0 | cluster_0    |          0.329744   |           0.323807  |           0.31638   |          0.178582   |           0.190506  |          0.183408   |          0.194283   |          0.314896   |          0.266949   |
|  1 | cluster_1    |          0.0875153  |           0.0826597 |           0.12148   |          0.0544798  |           0.0355812 |          0.0637625  |          0.0852074  |          0.0706873  |          0.075124   |
|  2 | cluster_2    |          0.0126884  |           0.0273287 |           0.0147922 |          0.386493   |           0.449984  |          0.361615   |          0.0247404  |          0.0136364  |          0.0211886  |
|  3 | cluster_3    |          0.20035    |           0.172857  |           0.192064  |          0.0404682  |           0.0345952 |          0.0217206  |          0.214885   |          0.180312   |          0.163982   |
|  4 | cluster_4    |          0.00825616 |           0.0143059 |           0.0245092 |          0.00313302 |           0.0107374 |          0.00292042 |          0.00648297 |          0.00469606 |          0.00751628 |


!!! tip "Tips"
    For the `count_table`, it’s best to use UMI counts. If deep sequencing samples (>30 000 UMIs per sample) are available, it is recommended to downsample them to the same number of UMIs for normalization and faster computation. For details on downsampling filter usage, visit [stats usage](usage_stats.md#filtering-clonosets) page.

Also, `clonosets_df` is required for creating a count table. Example (here, a subset of CD8 cells is selected):

```py
count_table_by_cluster_cd8
```

|    | sample_id         | chain   | filename         | Subset   |   Replica |   Number_of_cells | Peptid         |
|---:|:------------------|:--------|:-----------------|:---------|----------:|------------------:|:---------------|
|  0 | Luk_TRB_AdV_CD8_3 | TRB     | path_to_sample_0 | CD8      |         3 |             21031 | PepTivator_AdV |
|  2 | Luk_TRB_CMV_CD8_2 | TRB     | path_to_sample_1 | CD8      |         2 |             17490 | PepTivator_CMV |
|  3 | Luk_TRB_AdV_CD8_1 | TRB     | path_to_sample_2 | CD8      |         1 |             37698 | PepTivator_AdV |
|  5 | Luk_TRB_CMV_CD8_3 | TRB     | path_to_sample_3 | CD8      |         3 |             16662 | PepTivator_CMV |
|  6 | Luk_TRB_EBV_CD8_2 | TRB     | path_to_sample_4 | CD8      |         2 |             22020 | PepTivator_EBV |
| 12 | Luk_TRB_EBV_CD8_1 | TRB     | path_to_sample_5 | CD8      |         1 |             41274 | PepTivator_EBV |
| 14 | Luk_TRB_EBV_CD8_3 | TRB     | path_to_sample_6 | CD8      |         3 |             39014 | PepTivator_EBV |
| 16 | Luk_TRB_AdV_CD8_2 | TRB     | path_to_sample_7 | CD8      |         2 |             21359 | PepTivator_AdV |
| 17 | Luk_TRB_CMV_CD8_1 | TRB     | path_to_sample_8 | CD8      |         1 |             19617 | PepTivator_CMV |

<br>

## Count table by cluster

### Creating clusters

Here, only functional clonotypes are used (no frameshifts or premature stop codons) and one-node clusters are filtered out. 

```py
func_filter = clf.Filter(functionality="f", by_umi=True)

clusters_cd8 = clustering.create_clusters(cd8_clonosets, cl_filter=func_filter, mismatches=1, overlap_type="aaV", igh=False, tcrdist_radius=None, count_by_freq=False)

clusters_cd8_nonsingle = clustering.filter_one_node_clusters(clusters_cd8)
```

### `count_table_by_cluster`

!!! tip "Number of mismatches"
    It is recommended to set `mismatches` to 0, as the results with either 0 or 1 `mismatches` will be identical, however, using the former will result in faster computation. When the same clonoset is used for both clustering and intersections, the clone will be already present in the cluster.

```py
count_table_by_cluster_cd8 = intersections.count_table_by_cluster(cd8_clonosets, clusters_cd8_nonsingle, cl_filter=func_filter, overlap_type="aaV", mismatches=0, by_freq=True)

count_table_by_cluster_cd8.head()
```

|    | feature_id   |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |
|---:|:-------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
|  0 | cluster_0    |          0.329744   |           0.323807  |           0.31638   |          0.178582   |           0.190506  |          0.183408   |          0.194283   |          0.314896   |          0.266949   |
|  1 | cluster_1    |          0.0875153  |           0.0826597 |           0.12148   |          0.0544798  |           0.0355812 |          0.0637625  |          0.0852074  |          0.0706873  |          0.075124   |
|  2 | cluster_2    |          0.0126884  |           0.0273287 |           0.0147922 |          0.386493   |           0.449984  |          0.361615   |          0.0247404  |          0.0136364  |          0.0211886  |
|  3 | cluster_3    |          0.20035    |           0.172857  |           0.192064  |          0.0404682  |           0.0345952 |          0.0217206  |          0.214885   |          0.180312   |          0.163982   |
|  4 | cluster_4    |          0.00825616 |           0.0143059 |           0.0245092 |          0.00313302 |           0.0107374 |          0.00292042 |          0.00648297 |          0.00469606 |          0.00751628 |

<br>

## Differential expression calculation

 `wilcox_diff_expression` has two operating modes: 

  * If there are 2 group present, a pairwise comparison is done;
  * If the number of groups is higher than 2, each group is compared to all others, with additional columns for pairwise comparison between groups.

 Features (clonotypes or clusters) that appear in `min_samples` or more samples and have a `count_threshold` or higher counts or frequencies. This filters out unreliable differences, as the final results undergo FDR correction by the BH method. The default values for `min_samples` and `count_threshold` are `2` and `5e-04`, respectively. 

### For clusters

```py
diff_expression_by_cluster = diffexp.wilcox_diff_expression(count_table_by_cluster_cd8, sample_metadata, min_samples=2, count_threshold=0.00005)
```

|      | feature_id   |   logFC |   U |           p |     p_adj | pair                  |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |
|-----:|:-------------|--------:|----:|------------:|----------:|:----------------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
|  959 | cluster_71   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_EBV_vs_all |                   0 |                   0 |                   0 |           0         |          0          |          0          |         0.00337752  |         0.00498357  |         0.000248729 |
| 1007 | cluster_76   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_EBV_vs_all |                   0 |                   0 |                   0 |           0         |          0          |          0          |         0.000309606 |         0.000533954 |         0.000287593 |
| 1072 | cluster_83   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_CMV_vs_all |                   0 |                   0 |                   0 |           0.0018276 |          0.00372521 |          0.00638842 |         0           |         0           |         0           |



### For clonotype groups

Here, a clonotype group is defined based on the `overlap_type` used during `count_table` calculation. <br>Possible overlap types are [`aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`], `aa`/`nt` stands for an amino acid or nucleotide sequence, and `V`/`J`/`VJ` denote a segment type. 

```py
count_table_by_clonotypes_aaV_cd8 = intersections.count_table(cd8_clonosets, cl_filter=func_filter, overlap_type="aaV", mismatches=0, by_freq=False))
```

|    | feature_id                       |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |
|---:|:---------------------------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
|  0 | ('CASSPGLAAPYSEQFF', 'TRBV12-3') |                   1 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |
|  1 | ('CASRQGAGEQYF', 'TRBV19')       |                   1 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |
|  2 | ('CSAGTTYGTDIISQHF', 'TRBV20-1') |                   0 |                   4 |                   0 |                   0 |                   0 |                   0 |                  10 |                   0 |                   0 |


```py
diff_expression_by_clonotypes_aaV = diffexp.wilcox_diff_expression(count_table_by_clonotypes_aaV_cd8, sample_metadata, min_samples=2, count_threshold=2)
```
|     | feature_id                     |     logFC |   U |           p |     p_adj | pair                  |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |
|----:|:-------------------------------|----------:|----:|------------:|----------:|:----------------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|
|   2 | ('CASSSGLLGEQFF', 'TRBV6-5')   |   9.62006 |  36 | 0.00118822  | 0.0481036 | PepTivator_EBV_vs_all |                   0 |                   1 |                   0 |                   0 |                   0 |                   0 |                 421 |                 204 |                3141 |
| 111 | ('CASSSGLLDTQYF', 'TRBV6-5')   |   7.93852 |  36 | 0.00118822  | 0.0481036 | PepTivator_EBV_vs_all |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |                1428 |                1462 |                 614 |
| 118 | ('CASSRTSGSFLFEQYF', 'TRBV14') | 100       |  36 | 0.000304934 | 0.0201417 | PepTivator_AdV_vs_all |                 223 |                   3 |                   7 |                   0 |                   0 |                   0 |                   0 |                   0 |                   0 |


Output has the following columns:

*  <b>LogFC</b>: Log fold change of counts / frequencies between two groups. If all the samples in one group have zero counts or frequencies, logFC is set to 100 or -100.

* <b>U</b>: Mann-Whitney U statistic

* <b>p</b>: p-value

* <b>p_adj</b>: p-value after multiple testing correction using Benjamini-Hochberg (BH) procedure

* Pairwise group comparisons and their raw p-values (_p in a column name) and log fold changes (_logFC in a column name).

<br>

## Merge cluster properties by `feature_id`

* Create cluster properties

```py
cd8_cluster_props = clustering.cluster_properties(clusters_cd8_nonsingle, weighed=True).rename(columns = {"cluster_id": "feature_id"})

diff_expression_by_cluster.merge(cd8_cluster_props)
```

|    | feature_id   |   logFC |   U |           p |     p_adj | pair                  |   Luk_TRB_AdV_CD8_1 |   Luk_TRB_AdV_CD8_2 |   Luk_TRB_AdV_CD8_3 |   Luk_TRB_CMV_CD8_1 |   Luk_TRB_CMV_CD8_2 |   Luk_TRB_CMV_CD8_3 |   Luk_TRB_EBV_CD8_1 |   Luk_TRB_EBV_CD8_2 |   Luk_TRB_EBV_CD8_3 |   nodes |   edges |   diameter |   density |   eccentricity | concensus_cdr3aa   | concensus_cdr3nt                                    | concensus_v   | concensus_j   |
|---:|:-------------|--------:|----:|------------:|----------:|:----------------------|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------:|--------:|-----------:|----------:|---------------:|:-------------------|:----------------------------------------------------|:--------------|:--------------|
|  0 | cluster_71   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_EBV_vs_all |                   0 |                   0 |                   0 |           0         |          0          |          0          |         0.00337752  |         0.00498357  |         0.000248729 |       7 |      21 |          1 |  1        |        1       | CASSGASGSFNEQFF    | TGTGCCAGTAGTGGGGCTAGCGGGAGTTTTAATGAGCAGTTCTTC       | TRBV19        | TRBJ2-1       |
|  1 | cluster_76   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_EBV_vs_all |                   0 |                   0 |                   0 |           0         |          0          |          0          |         0.000309606 |         0.000533954 |         0.000287593 |       7 |      21 |          1 |  1        |        1       | CASTSHGTSKDNEQFF   | TGTGCCAGCACGTCTCACGGGACTAGCAAAGACAATGAGCAGTTCTTC    | TRBV7-6       | TRBJ2-1       |
|  2 | cluster_83   |     100 |  36 | 0.000304934 | 0.0227611 | PepTivator_CMV_vs_all |                   0 |                   0 |                   0 |           0.0018276 |          0.00372521 |          0.00638842 |         0           |         0           |         0           |       6 |      14 |          2 |  0.933333 |        1.33333 | CASSFLLGQGADYEQYF  | TGTGCCAGCAGCTTCCTCCTGGGACAGGGGGCCGACTACGAGCAGTACTTC | TRBV11-2      | TRBJ2-7       |

??? info "Visualization"