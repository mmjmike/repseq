# Usage: differential enrichment analysis

The `diff_enrichment` module identifies clonotypes or other repertoire features
that are enriched in one experimental group. It is designed to work directly
with the wide tables returned by `intersections.count_table`.

```py
import pandas as pd

from repseq import diff_enrichment as rsde
from repseq import intersections
from repseq import plot as rsplot
```

## Stateful analysis with `Analyzer`

`Analyzer` manages the complete count-table, prefilter, statistics, postfilter,
and optional chain-pairing workflow while caching each completed step.

```python
analyzer = rsde.Analyzer(
    samples_df=samples_df,
    method="mann_whitney",
    min_samples=3,
)
analyzer.run()

statistics = analyzer.statistics_df
filtered = analyzer.postfiltered
```

Parameters can also be supplied as a dictionary or updated later. Updating a
parameter clears only results downstream of the affected step.

```python
analyzer.update_parameters({"min_count": 3, "max_p_adj": 0.05})
analyzer.run()
```

For paired chains, select a branch explicitly or define chain-specific
overrides with a suffix. A callable result table retrieves another branch
without changing the active chain.

```python
analyzer.update_parameters({"overlap_type_TRA": "aaVJ"})
analyzer.select_chain("TRB")
trb_statistics = analyzer.statistics_df
tra_statistics = analyzer.statistics_df("TRA")
```

Supported paired-chain combinations are `TRA`–`TRB`, `TRG`–`TRD`, and
`IGH`–`IGKL`. The analyzer normalizes supported `TRAD`, `IGK`, and `IGL`
aliases when their paired chain identifies the intended branch.

Pairing scores can be shown as a heatmap. The analyzer uses `-log10` colors
for JSD scores and untransformed colors for other pairing methods.

```python
pairing_figure = analyzer.plot_pairing()
```

## Input data

`rsde.calc_statistics` requires two dataframes:

1. A wide count table with one row per feature and numeric sample columns.
2. A sample metadata table containing `sample_id` and `group` columns.

Numeric count-table sample columns are matched to metadata by `sample_id`.
`samples_metadata` may contain additional samples that are absent from the
count table; they are reported and ignored. Numeric table columns that are not
listed in the metadata are also reported and ignored. Group detection and the
minimum group-size checks use only samples represented in both inputs. Each
represented group must contain at least two samples, and at least two groups
are required.

Example metadata:

| sample_id | group |
|:----------|:------|
| Pep1_1 | Pep1 |
| Pep1_2 | Pep1 |
| Pep1_3 | Pep1 |
| Control_1 | Control |
| Control_2 | Control |
| Control_3 | Control |

```py
samples_metadata = pd.DataFrame(
    {
        "sample_id": [
            "Pep1_1",
            "Pep1_2",
            "Pep1_3",
            "Control_1",
            "Control_2",
            "Control_3",
        ],
        "group": [
            "Pep1",
            "Pep1",
            "Pep1",
            "Control",
            "Control",
            "Control",
        ],
    }
)
```

## Creating a count table

Use counts rather than frequencies for count-based methods:

```py
count_table = intersections.count_table(
    clonosets_df,
    overlap_type="aaVJ",
    by_freq=False,
)
```

A typical table contains feature-description columns followed by one numeric
column per sample:

| clonotype | cdr3aa | v | j | Pep1_1 | Pep1_2 | Pep1_3 | Control_1 | Control_2 | Control_3 |
|:----------|:-------|:--|:--|-------:|-------:|-------:|----------:|----------:|----------:|
| CASSAAA\|TRBV1\|TRBJ1 | CASSAAA | TRBV1 | TRBJ1 | 5 | 4 | 7 | 0 | 0 | 1 |
| CASSBBB\|TRBV2\|TRBJ2 | CASSBBB | TRBV2 | TRBJ2 | 0 | 1 | 0 | 3 | 4 | 2 |

## Prefiltering features

`rsde.prefilter` marks features that have sufficient total abundance and
recurrence across samples. It does not remove rows.

```py
count_table = rsde.prefilter(
    count_table,
    min_samples=3,
    min_count=2,
    min_total_count=10,
)
```

A feature passes when both conditions are true:

- At least `min_samples` numeric columns have a value greater than or equal to
  `min_count`.
- The sum across all numeric columns is greater than or equal to
  `min_total_count`.

The boolean `prefilter_pass` column is inserted before the numeric columns. If
this column is present, `calc_statistics` tests only rows where it is `True`.
Rows that fail remain in the result with missing statistical values.

```py
passed_count_table = count_table[count_table["prefilter_pass"]]
```

!!! warning
    `prefilter` uses every numeric column. Numeric annotation columns should be
    converted to a non-numeric dtype or removed before prefiltering if they
    must not contribute to the thresholds.

## Running the default analysis

The default method is a Mann–Whitney U test:

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="mann_whitney",
    cpu=4,
)
```

When `feature_column=None`, the first count-table column is used as the feature
identifier. Its values must be unique. An explicit column can be selected when
the first column is not an appropriate identifier:

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    feature_column="clonotype",
    method="mann_whitney",
)
```

At startup, the function reports:

- the selected method and feature column;
- the number of groups;
- each group and its sample IDs;
- numeric sample columns ignored because they are absent from the metadata;
- whether prefiltering is active and how many features will be tested.

Set `verbose=False` to suppress this output and parallel progress reporting.

## Sorting results

Both `rsde.calc_statistics` and `rsde.postfilter` sort their output by
default using:

```py
[
    "prefilter_pass",
    "postfilter_pass",
    "enriched_in",
    "mean_group_count",
    "log2FC",
    "p_adj",
]
```

Filter columns place `True` first, `mean_group_count` and `log2FC` sort
descending, and all other columns sort ascending. Categorical
`enriched_in` values follow their category order; other values sort
alphabetically. Columns absent from the output are silently ignored.

Pass another column sequence to `sort`, or disable sorting with
`sort=False`, `sort=None`, `sort=[]`, or a list containing only columns
that are absent from the output. A single column name is also accepted:

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    sort=["p_adj", "feature"],
)

unsorted_result = rsde.postfilter(result, sort=False)
```

## Two-group and multi-group analyses

With exactly two groups, one two-sided comparison is calculated. The
`enriched_in` value is the group with the greater mean count, and `log2FC` is
reported in that group's direction.

With more than two groups, each group is compared with all remaining samples:

```text
Pep1 vs all
Pep2 vs all
Control vs all
```

These tests use the one-sided alternative that the target group is greater.
Group-vs-all calculations are dispatched through
`run_parallel_calculation`. Pairwise group comparisons are not calculated.

## Statistical methods

`presence_threshold` is a common analysis parameter for every method. Before
statistics are calculated, represented sample counts lower than this threshold
are set to zero in an internal analysis matrix. The original count-table values
are not modified and are returned unchanged. `log2FC` and statistical tests use
the thresholded values, while `mean_group_count` reports the mean of the real,
unmodified counts for the selected `enriched_in` group.

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="mann_whitney",
    presence_threshold=2,
)
```

| Method | Test | Method-specific parameters | Suggested role |
|:-------|:-----|:---------------------------|:---------------|
| `mann_whitney` | Mann–Whitney U test on replicate counts | None beyond common parameters | General replicate-level abundance comparison; default |
| `fisher` | Fisher exact test on detected/not-detected replicates | `presence_threshold` | Robust recurrence evidence for sparse clonotypes |
| `fisher_count` | Fisher exact test on aggregated feature and non-feature counts | `sample_totals` | Secondary aggregated UMI/count evidence |
| `hurdle` | Recurrence Fisher test plus positive-abundance Mann–Whitney test | `presence_threshold`, `sample_totals`, `hurdle_combine_method`, `cpm_scale`, `pseudocount` | Combined recurrence and abundance evidence |
| `quasi_binomial` | Quasi-binomial GLM using sample totals | `sample_totals` | Depth-aware replicate-level proportion model |
| `negative_binomial` | Negative-binomial GLM with a library-size offset | `sample_totals`, `negative_binomial_alpha` | Advanced depth-aware count model |
| `permutation` | Group-label permutation test on log-CPM differences | `sample_totals`, `n_permutations`, `random_state`, `cpm_scale`, `pseudocount` | Low-assumption enrichment test |

Parameters that belong to unselected methods are ignored.

### Recurrence with Fisher's exact test

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="fisher",
    presence_threshold=2,
)
```

A sample is considered positive when its count is greater than or equal to
`presence_threshold`. This method tests whether a feature is detected in more
replicates of the target group.

### Aggregated count enrichment

```py
sample_totals = {
    "Pep1_1": 25_000,
    "Pep1_2": 27_500,
    "Pep1_3": 24_800,
    "Control_1": 31_000,
    "Control_2": 29_500,
    "Control_3": 30_200,
}

result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="fisher_count",
    sample_totals=sample_totals,
)
```

`fisher_count` requires non-negative integer counts and integer sample totals.
If `sample_totals` is omitted for a count-aware method, totals are calculated
as the sum of each represented sample column across the complete input table.

!!! important
    Supply external library totals when the count table contains only a subset
    of the repertoire. Aggregated Fisher testing treats counts as independent
    observations and can produce very small p-values, so it is best used as
    secondary evidence rather than the only enrichment test.

### Hurdle analysis

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="hurdle",
    presence_threshold=2,
    sample_totals=sample_totals,
    hurdle_combine_method="fisher",
)
```

The hurdle method combines a recurrence p-value with a Mann–Whitney p-value
calculated from log-CPM values among detected samples. Set
`hurdle_combine_method="max"` for a more conservative combination.

### Permutation analysis

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="permutation",
    sample_totals=sample_totals,
    n_permutations=10_000,
    random_state=1,
    cpu=4,
)
```

`random_state` makes the result reproducible. More permutations improve
p-value resolution but increase runtime.

## Multiple-testing correction

P-values are adjusted together after all group comparisons are concatenated.
Benjamini–Hochberg FDR correction is the default:

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    p_adjust_method="fdr_bh",
)
```

Other methods accepted by `statsmodels.stats.multitest.multipletests` can be
selected, such as `bonferroni`, `holm`, or `fdr_by`.

## Result columns

The following columns are inserted immediately before the original numeric
columns:

| Column | Description |
|:-------|:------------|
| `enriched_in` | Group associated with the selected comparison |
| `method` | Statistical method used |
| `mean_group_count` | Mean feature count across samples belonging to `enriched_in` |
| `log2FC` | Log2 fold change for `enriched_in` against the comparison samples |
| `p_val` | Raw p-value |
| `p_adj` | Multiple-testing-adjusted p-value |

When exactly one comparison mean is zero, `log2FC` uses the finite value from
`log2fc_zero_value`, which defaults to `100`. When both means are zero, it is
`0`.

Candidate features can be selected using both statistical significance and an
effect-size threshold:

```py
candidates = result[
    result["prefilter_pass"]
    & (result["p_adj"] < 0.05)
    & (result["log2FC"] >= 2)
].copy()
```

## Postfiltering statistical results

`rsde.postfilter` adds a boolean `postfilter_pass` column without removing any
rows. The column is placed after `p_adj`, immediately before the original
count-table columns. P-value maxima are optional and use strict comparisons;
minimum effect-size thresholds are inclusive:

```py
result = rsde.postfilter(
    result,
    max_p_adj=None,
    max_p_val=None,
    min_group_mean=2,
    min_logfc=1,
    verbose=True,
)
```

A row passes when all of the following are true:

- `p_adj < max_p_adj`, when `max_p_adj` is not `None`;
- `p_val < max_p_val`, when `max_p_val` is not `None`;
- `mean_group_count >= min_group_mean`;
- `log2FC >= min_logfc`.

Values exactly equal to a p-value maximum fail, while values exactly equal to a
minimum mean or log2FC threshold pass. `None` means that any value is accepted
for that p-value column. Rows missing a value required by an active filter
receive `postfilter_pass=False`.

By default, `enriched_in` is not filtered. Select one group with a string or
several groups with a list:

```py
pep1_result = rsde.postfilter(
    result,
    max_p_adj=0.05,
    min_group_mean=5,
    min_logfc=2,
    groups="Pep1",
)

selected_groups = rsde.postfilter(
    result,
    groups=["Pep1", "Pep2"],
)

without_controls = rsde.postfilter(
    result,
    groups_exclude=["Control"],
)
```

`groups_exclude` accepts the same string or sequence forms as `groups`, but
removes matching `enriched_in` groups instead. Values absent from `enriched_in`
are silently ignored. If both options are provided, `groups` takes precedence
and verbose output reports that `groups_exclude` was skipped.

Retrieve the passing rows with:

```py
postfiltered_result = result[result["postfilter_pass"]].copy()
```

With `verbose=True`, the function prints every active comparison explicitly,
using `<` for p-values and `>=` for minimum thresholds. Disabled p-value and
group filters are printed as `any`. If `prefilter_pass` is present, only its
`True` rows can pass and the output reports both the prefilter count out of the
whole table and the postfilter count out of prefilter-passing features. Without
`prefilter_pass`, the postfilter count is reported out of the full table.

## Simplified and expanded results

`simplify=True` is the default. For analyses with more than two groups, it
keeps the comparison with the lowest raw p-value for each tested feature:

```py
result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    simplify=True,
)
```

Use `simplify=False` to retain every group-vs-all result. Each passing feature
then occurs once per group, while a feature that failed prefiltering remains a
single row with missing statistical values:

```py
long_result = rsde.calc_statistics(
    count_table,
    samples_metadata,
    method="fisher",
    presence_threshold=2,
    simplify=False,
)
```

## Differential-enrichment heatmap

After filtering the statistics table to the features of interest, plot their
original sample counts with `rsplot.de_heatmap`:

```py
filtered_result = result[
    result["prefilter_pass"]
    & (result["p_adj"] < 0.05)
    & (result["log2FC"] >= 2)
].copy()

fig = rsplot.de_heatmap(
    filtered_result,
    samples_metadata,
    feature_column="clonotype",
    log_values=True,
)
```

If `postfilter_pass` is present, only its `True` rows are plotted. The same rule
applies to `prefilter_pass`. When both columns are present, a feature must be
`True` in both columns to appear in the heatmap. Filtering is applied before
checking `enriched_in`, so excluded rows may contain missing statistical
values.

The function identifies count columns by matching numeric table columns to
`samples_metadata["sample_id"]`. Statistical columns such as `log2FC`,
`p_val`, and `p_adj` are therefore not plotted.

Samples are displayed together by `group`. If the count-table columns already
contain contiguous group blocks, their existing order is preserved. If groups
are interleaved, samples are stably regrouped without changing their order
within each group. White vertical gaps separate adjacent sample groups.

The top annotation strip shows the sample group and the left annotation strip
shows each feature's `enriched_in` group. Both strips use the same group-to-color
mapping. Every non-missing `enriched_in` value must occur in
`samples_metadata["group"]`.

`log_values=True` is the default and follows the same rule as
`rsplot.beta_metric`: zeros are
replaced by one tenth of the smallest positive plotted count before applying
log10. Cell labels, when `show_values=True`, always display the original count
values. Integer-valued counts are written in full without a trailing `.0` or
scientific notation, so a count such as `456789` remains `456789`. Genuine
decimal values retain compact decimal formatting. The default heatmap colors
are the R `pheatmap` blue-to-red palette.

```py
fig = rsplot.de_heatmap(
    filtered_result,
    samples_metadata,
    show_values=False,
    group_palette={
        "Pep1": "#1f77b4",
        "Control": "#ff7f0e",
    },
)
```

## Pairing features between chains

`rsde.pair_chains` compares feature abundance patterns between two
chain-specific `calc_statistics` tables. For example, it can rank candidate
TRA/TRB feature pairs across matched biological samples.

The pairing metadata requires:

- `sample_id`: the chain-specific count-table column name;
- `sample`: the shared biological sample identifier used to pair the chains.

Each represented biological sample must have exactly one sample ID in each
count table.

| sample_id | sample |
|:----------|:-------|
| donor1_TRA | donor1 |
| donor1_TRB | donor1 |
| donor2_TRA | donor2 |
| donor2_TRB | donor2 |
| donor3_TRA | donor3 |
| donor3_TRB | donor3 |

```py
paired_metadata = pd.DataFrame(
    {
        "sample_id": [
            "donor1_TRA",
            "donor1_TRB",
            "donor2_TRA",
            "donor2_TRB",
            "donor3_TRA",
            "donor3_TRB",
        ],
        "sample": [
            "donor1",
            "donor1",
            "donor2",
            "donor2",
            "donor3",
            "donor3",
        ],
    }
)
```

Calculate Jensen–Shannon divergence with:

```py
pairing_matrix = rsde.pair_chains(
    tra_statistics,
    trb_statistics,
    paired_metadata,
    method="jsd",
)
```

Rows from each table are normalized to sum to one across paired samples before
scoring. For JSD, smaller values indicate more similar abundance patterns and
`0` indicates identical normalized patterns.

Pearson correlation is also available:

```py
pairing_matrix = rsde.pair_chains(
    tra_statistics,
    trb_statistics,
    paired_metadata,
    method="pearson",
)
```

For Pearson scores, larger values indicate better agreement. Features with
constant normalized profiles receive `NaN` correlations.

The matrix columns are features from the first count table and rows are
features from the second count table:

```text
                         count_table1 features
                       TRA_1    TRA_2    TRA_3
count_table2  TRB_1     0.00     0.54     0.31
features      TRB_2     0.48     0.02     0.60
              TRB_3     0.25     0.41     0.01
```

When `prefilter_pass` is present, only `True` rows are eligible. Without filter
lists, all eligible features are returned. Filter lists retain requested IDs
and add their best-scoring partners from the opposite chain:

```py
selected_pairs = rsde.pair_chains(
    tra_statistics,
    trb_statistics,
    paired_metadata,
    method="jsd",
    filter_ids1=["TRA_feature_1", "TRA_feature_2"],
    filter_ids2=["TRB_feature_5"],
)
```

In this example, the result includes both requested TRA features, the requested
TRB feature, the best TRB partner for each requested TRA feature, and the best
TRA partner for the requested TRB feature. JSD selects the minimum score;
Pearson selects the maximum score.

By default, the first column of each table is used as its feature ID. If both
tables use the same feature-column name, pass one value. Different names can be
provided as a two-item tuple:

```py
pairing_matrix = rsde.pair_chains(
    tra_statistics,
    trb_statistics,
    paired_metadata,
    feature_column=("tra_feature", "trb_feature"),
)
```

## Differential-enrichment volcano plot

`rsplot.de_volcano` plots effect size against statistical significance. The
default vertical axis uses adjusted p-values, point color represents
`enriched_in`, and point area is proportional to
`log1p(mean_group_count)` by default. The `log1p` transformation keeps zero
counts valid while reducing domination by very abundant features.

```py
fig = rsplot.de_volcano(
    result,
    p_column="p_adj",
    size_range=(20, 300),
    log_sizes=True,
)
```

Rows with missing `log2FC`, selected p-value, `mean_group_count`, or
`enriched_in` values are silently omitted. Use raw p-values with:

```py
fig = rsplot.de_volcano(
    result,
    p_column="p_val",
)
```

The y-axis is `-log10(p_adj)` or `-log10(p_val)`. Zero p-values are replaced by
one tenth of the smallest positive plotted p-value before transformation, so
they remain finite on the plot.

Repeated identical `log2FC` values greater than `10` are adjusted for plotting.
For a repeated value `X`, the largest plotted value below `X` is found and all
copies of `X` are replaced by `min(max_non_X + 2, X)`. This changes only the
displayed x-coordinate; the input dataframe is not modified.

Group colors and point-size limits can be customized:

```py
fig = rsplot.de_volcano(
    result,
    group_palette={
        "Pep1": "#1f77b4",
        "Control": "#ff7f0e",
    },
    size_range=(30, 500),
    alpha=0.8,
)
```

Use raw `mean_group_count` values for point sizing when needed:

```py
fig = rsplot.de_volcano(
    result,
    log_sizes=False,
)
```

## API reference

### `prefilter`

::: repseq.diff_enrichment.prefilter
    options:
      show_root_heading: true
      show_source: false

### `calc_statistics`

::: repseq.diff_enrichment.calc_statistics
    options:
      show_root_heading: true
      show_source: false

### `postfilter`

::: repseq.diff_enrichment.postfilter
    options:
      show_root_heading: true
      show_source: false
