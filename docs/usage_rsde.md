# Differential enrichment with `rsde`

The `repseq.diff_enrichment` module, conventionally imported as `rsde`, finds
clonotypes, clusters, or other repertoire features enriched in one sample
group. It supports two equivalent analysis styles:

- **Object-oriented workflow:** `rsde.Analyzer` manages count-table creation,
  filtering, statistics, paired chains, cached results, and plots.
- **Procedural workflow:** call `rsde.prefilter`, `rsde.calc_statistics`,
  `rsde.postfilter`, and `rsde.pair_chains` directly for maximum control.

```py
import pandas as pd

from repseq import diff_enrichment as rsde
from repseq import intersections
from repseq import plot as rsplot
```

## Choosing a workflow

Use `rsde.Analyzer` when you want:

- one object that stores parameters and intermediate tables;
- automatic count-table or cluster-table construction;
- selective cache invalidation after parameter changes;
- separate analysis branches for paired receptor chains;
- a single `run()` command and plotting wrappers.

Use the procedural API when you want:

- to start from an existing count table;
- complete control over when each function runs;
- custom transformations between analysis steps;
- to integrate the functions into another pipeline framework.

Both styles produce the same core table formats and use the same statistical
functions.

## Shared concepts

### Analysis stages

The standard differential-enrichment workflow is:

1. Read sample metadata.
2. Create a feature-by-sample count table.
3. Mark sufficiently abundant features with `prefilter_pass`.
4. Calculate group enrichment statistics.
5. Mark statistically and biologically relevant rows with `postfilter_pass`.
6. Optionally pair features from two receptor chains.

`prefilter` and `postfilter` add boolean columns but do not remove rows. This
keeps original counts, failed features, and statistical results in one table.

### Sample metadata

For a single-chain analysis, metadata must contain:

| Column | Required | Meaning |
|:--|:--:|:--|
| `sample_id` | yes | Unique identifier matching the count-table sample column |
| `filename` | for count-table creation | Path to the clonotype table |
| `group` | yes | Experimental group used by differential enrichment |
| `sample` | no | Biological sample name; required for paired-chain analysis |
| `chain` | no | Receptor chain; omitted metadata is treated as chain `XCR` |

Example single-chain metadata:

```py
samples_df = pd.DataFrame(
    {
        "sample_id": [
            "stim_1",
            "stim_2",
            "stim_3",
            "control_1",
            "control_2",
            "control_3",
        ],
        "filename": [
            "data/stim_1.tsv",
            "data/stim_2.tsv",
            "data/stim_3.tsv",
            "data/control_1.tsv",
            "data/control_2.tsv",
            "data/control_3.tsv",
        ],
        "group": [
            "stimulated",
            "stimulated",
            "stimulated",
            "control",
            "control",
            "control",
        ],
    }
)
```

For two chains, `sample` identifies the biological pair and `sample_id` is
unique within each chain:

```py
paired_samples_df = pd.DataFrame(
    {
        "sample": [
            "case_1",
            "case_2",
            "control_1",
            "control_2",
            "case_1",
            "case_2",
            "control_1",
            "control_2",
        ],
        "sample_id": [
            "case1_TRA",
            "case2_TRA",
            "control1_TRA",
            "control2_TRA",
            "case1_TRB",
            "case2_TRB",
            "control1_TRB",
            "control2_TRB",
        ],
        "filename": [
            "data/case1_TRA.tsv",
            "data/case2_TRA.tsv",
            "data/control1_TRA.tsv",
            "data/control2_TRA.tsv",
            "data/case1_TRB.tsv",
            "data/case2_TRB.tsv",
            "data/control1_TRB.tsv",
            "data/control2_TRB.tsv",
        ],
        "group": [
            "case",
            "case",
            "control",
            "control",
            "case",
            "case",
            "control",
            "control",
        ],
        "chain": ["TRA", "TRA", "TRA", "TRA", "TRB", "TRB", "TRB", "TRB"],
    }
)
```

### Count-table format

A count table contains feature-description columns followed by one numeric
column per represented `sample_id`:

| clonotype | cdr3aa | v | j | stim_1 | stim_2 | control_1 | control_2 |
|:--|:--|:--|:--|--:|--:|--:|--:|
| CASSAAA\|TRBV1\|TRBJ1 | CASSAAA | TRBV1 | TRBJ1 | 8 | 5 | 0 | 1 |
| CASSBBB\|TRBV2\|TRBJ2 | CASSBBB | TRBV2 | TRBJ2 | 0 | 1 | 6 | 4 |

The first column is used as the feature identifier by default. It must contain
unique values. Pass `feature_column` to `calc_statistics` when another column
should identify features.

!!! warning
    `prefilter` treats every numeric column as a sample column. Remove numeric
    annotation columns or convert them to a non-numeric dtype before
    prefiltering.

## Object-oriented workflow with `Analyzer`

### Minimal complete analysis

Create an analyzer with metadata and any non-default parameters, then call
`run()`:

```py
analyzer = rsde.Analyzer(
    samples_df=samples_df,
    overlap_type="aaVJ",
    min_samples=3,
    min_count=2,
    min_total_count=10,
    method="mann_whitney",
    max_p_adj=0.05,
    min_logfc=1,
    cpu=4,
)

analyzer.run()
```

`run()` executes count-table creation, prefiltering, statistics, and
postfiltering. For two chains it runs one branch per chain and optionally
merges them with feature pairing.

Printing or evaluating the analyzer shows its current state:

```py
print(analyzer)
```

The state includes groups, chains, active chain, completed stages, and the
number of features passing prefilter and postfilter.

### Constructing with a parameter dictionary

Parameters may be supplied directly, through a dictionary, or with both forms.
Direct keyword arguments override duplicate dictionary entries.

```py
parameters = {
    "overlap_type": "aaVJ",
    "min_samples": 3,
    "method": "hurdle",
    "max_p_adj": 0.05,
}

analyzer = rsde.Analyzer(
    parameters=parameters,
    samples_df=samples_df,
    cpu=4,
)
```

Inspect all effective global and chain-specific settings with:

```py
parameters = analyzer.get_parameters()
```

### Reading or replacing metadata

Metadata can be supplied later:

```py
analyzer = rsde.Analyzer(verbose=True)
analyzer.read_samples(samples_df)
```

Reading different metadata clears every cached result. Reading an equal
dataframe again leaves cached results intact.

The `samples_df` property returns the effective metadata used by the analyzer.
It includes the default `XCR` chain when no `chain` column was supplied and may
include normalized chain names.

### Running one stage at a time

The workflow can be executed interactively:

```py
analyzer.run_count_table()
analyzer.run_prefilter()
analyzer.run_statistics()
analyzer.run_postfilter()
```

Each stage accepts its relevant parameters:

```py
analyzer.run_count_table(
    overlap_type="aaVJ",
    count_by_freq=False,
)

analyzer.run_prefilter(
    min_samples=2,
    min_count=3,
    min_total_count=15,
)

analyzer.run_statistics(
    method="permutation",
    n_permutations=20_000,
    p_adjust_method="fdr_bh",
)

analyzer.run_postfilter(
    max_p_adj=0.05,
    min_group_mean=2,
    min_logfc=1,
)
```

Calling a stage again with the same effective parameters reuses the cached
table and reports that no calculation is needed.

### Count table versus clustering

By default, `run_count_table()` calls `intersections.count_table` with
`mismatches=0`:

```py
analyzer.update_parameters(
    {
        "clustering": False,
        "overlap_type": "aaV",
        "count_by_freq": False,
    }
)
analyzer.run_count_table()
```

Enable clustering to pool clonotypes into connected components before making
the table:

```py
analyzer.run_count_table(
    clustering=True,
    overlap_type="aaVJ",
    mismatches=1,
)
```

`mismatches` affects only clustering. The direct count-table path always uses
exact matches.

### Updating parameters and cache invalidation

Use `update_parameters` with a dictionary or keyword arguments:

```py
analyzer.update_parameters(
    {
        "min_count": 3,
        "max_p_adj": 0.01,
    }
)
```

Only affected stages and their descendants are cleared:

| Changed parameter category | Cleared results |
|:--|:--|
| Count-table parameters | count table, prefilter, statistics, postfilter, pairing |
| Prefilter parameters | prefilter, statistics, postfilter, pairing |
| Statistics parameters | statistics, postfilter, pairing |
| Postfilter parameters | postfilter, pairing |
| Pairing parameters | pairing only |
| `cpu`, `verbose`, `sort` | no immediate cache deletion |

After an update, call `run()` to calculate only missing or changed stages:

```py
analyzer.update_parameters(max_p_adj=0.01)
analyzer.run()
```

### Accessing result tables

The active chain's tables are available as properties:

```py
count_table = analyzer.count_table
prefiltered = analyzer.prefiltered
statistics = analyzer.statistics_df
postfiltered = analyzer.postfiltered
pairing_matrix = analyzer.pairing_matrix
```

If a table has not been calculated, the analyzer prints a message and returns
`None`.

With two chains, switch the active branch:

```py
analyzer.select_chain("TRA")
tra_statistics = analyzer.statistics_df
```

Calculated tables are callable, so another chain can be retrieved without
changing the active branch:

```py
tra_statistics = analyzer.statistics_df("TRA")
trb_statistics = analyzer.statistics_df("TRB")
```

The same syntax works for `count_table`, `prefiltered`, and `postfiltered`.

### Supported paired chains and aliases

`Analyzer` supports at most two chains and recognizes these pairs:

- `TRA`–`TRB`
- `TRG`–`TRD`
- `IGH`–`IGKL`

Context-dependent aliases are normalized when metadata is read:

- `TRAD` with `TRB` becomes `TRA`.
- `TRAD` with `TRG` becomes `TRD`.
- `IGK` with `IGH` becomes `IGKL`.
- `IGL` with `IGH` becomes `IGKL`.

Metadata containing `IGH`, `IGK`, and `IGL` simultaneously is rejected. Merge
the light-chain repertoires manually and label them `IGKL`.

When chain pairing is enabled, two-chain metadata must contain `sample` and
each biological sample must have one represented sample ID per chain.

Inspect available chains and the current branch with:

```py
chains = analyzer.chains
analyzer.select_chain("TRB")
```

### Chain-specific parameters

Most stage parameters can override the global value for one chain by appending
the chain name:

```py
analyzer.update_parameters(
    {
        "overlap_type": "aaV",
        "overlap_type_TRA": "aaVJ",
        "min_count_TRA": 3,
        "min_count_TRB": 2,
    }
)
```

If a chain-specific value is absent, that branch uses the global value.

The following parameters support chain suffixes:

```text
cl_filter, overlap_type, mismatches, clustering, count_by_freq,
min_samples, min_count, min_total_count,
method, p_adjust_method, presence_threshold, sample_totals,
hurdle_combine_method, cpm_scale, pseudocount, n_permutations,
negative_binomial_alpha,
max_p_adj, max_p_val, min_group_mean, min_logfc, sort
```

When a stage method is called with new parameters in a two-chain analysis, a
chain-specific value is created for the active chain:

```py
analyzer.select_chain("TRA")
analyzer.run_count_table(overlap_type="aaVJ")

# TRB still uses the global overlap_type.
analyzer.select_chain("TRB")
analyzer.run_count_table()
```

### Pairing two chains

By default, `run()` pairs two postfiltered branches using Jensen–Shannon
divergence:

```py
analyzer = rsde.Analyzer(
    samples_df=paired_samples_df,
    pairing_method="jsd",
    pair_chains=True,
)
analyzer.run()

pairing_matrix = analyzer.pairing_matrix
```

Pairing can be run separately or recalculated with Pearson correlation:

```py
analyzer.pair_chains(pairing_method="pearson")
```

Set `pair_chains=False` to run both differential-enrichment branches without
creating a pairing matrix.

### Analyzer plotting methods

#### Volcano plot

`plot_volcano` uses the active chain's postfiltered table when available and
falls back to `statistics_df` otherwise:

```py
fig = analyzer.plot_volcano()
```

Force the unpostfiltered statistics table with `postfiltered=False`:

```py
fig = analyzer.plot_volcano(
    postfiltered=False,
    p_column="p_val",
)
```

Plot mean group count on the vertical axis and significance as point size.
With the default `log_sizes=True`, the vertical values use
`log1p(mean_group_count)`; set `log_sizes=False` to show raw counts:

```py
fig = analyzer.plot_volcano(by_mean_count=True)

raw_count_fig = analyzer.plot_volcano(
    by_mean_count=True,
    log_sizes=False,
)
```

Select a chain without changing the active branch:

```py
fig = analyzer.plot_volcano(chain="TRA")
```

All additional keyword arguments are forwarded to `rsplot.de_volcano`.

#### Differential-enrichment heatmap

`plot_heatmap` uses the postfiltered table and chain-specific metadata:

```py
fig = analyzer.plot_heatmap(
    chain="TRB",
    show_values=False,
    height=7,
)
```

Additional keyword arguments are forwarded to `rsplot.de_heatmap`.

#### Pairing heatmap

```py
fig = analyzer.plot_pairing(
    show_values=True,
    hclust=False,
)
```

The wrapper automatically uses `-log10` colors for JSD and untransformed
colors for other pairing methods.

### Analyzer parameter reference

The principal global parameters are:

| Parameter | Default | Used by |
|:--|:--|:--|
| `samples_df` | `None` | input metadata |
| `cl_filter` | `Filter(functionality="f", by_umi=True)` | count table or clustering |
| `overlap_type` | `"aaV"` | count table or clustering |
| `mismatches` | `1` | clustering only |
| `clustering` | `False` | count-table strategy |
| `count_by_freq` | `False` | count table or cluster table |
| `min_samples` | `3` | prefilter |
| `min_count` | `2` | prefilter |
| `min_total_count` | `10` | prefilter |
| `method` | `"mann_whitney"` | statistics |
| `simplify` | `True` | statistics |
| `p_adjust_method` | `"fdr_bh"` | statistics |
| `log2fc_zero_value` | `100` | statistics |
| `presence_threshold` | `2` | statistics |
| `sample_totals` | `None` | count-aware statistics |
| `hurdle_combine_method` | `"fisher"` | hurdle statistics |
| `cpm_scale` | `1_000_000` | normalized count methods |
| `pseudocount` | `0.5` | normalized count methods |
| `n_permutations` | `10_000` | permutation statistics |
| `negative_binomial_alpha` | `1.0` | negative-binomial statistics |
| `max_p_adj` | `None` | postfilter |
| `max_p_val` | `None` | postfilter |
| `min_group_mean` | `2` | postfilter |
| `min_logfc` | `1` | postfilter |
| `groups` | `None` | postfilter |
| `groups_exclude` | `None` | postfilter |
| `sort` | `DEFAULT_SORT_COLUMNS` | statistics and postfilter |
| `cpu` | `None` | parallel calculations |
| `verbose` | `True` | all verbose procedures |
| `pair_chains` | `True` | workflow orchestration |
| `pairing_method` | `"jsd"` | chain pairing |

!!! note
    `Analyzer` defaults `presence_threshold` to `2`, while a direct
    `rsde.calc_statistics` call defaults it to `1`.

## Procedural workflow

### Complete procedural example

```py
count_table = intersections.count_table(
    samples_df,
    overlap_type="aaVJ",
    mismatches=0,
    by_freq=False,
    cpu=4,
)

prefiltered = rsde.prefilter(
    count_table,
    min_samples=3,
    min_count=2,
    min_total_count=10,
)

statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="mann_whitney",
    presence_threshold=2,
    p_adjust_method="fdr_bh",
    cpu=4,
)

postfiltered = rsde.postfilter(
    statistics,
    max_p_adj=0.05,
    min_group_mean=2,
    min_logfc=1,
)

significant = postfiltered[postfiltered["postfilter_pass"]]
```

### Creating a count table

Use counts for count-aware statistical methods:

```py
count_table = intersections.count_table(
    samples_df,
    cl_filter=None,
    overlap_type="aaVJ",
    mismatches=0,
    by_freq=False,
    cpu=4,
    verbose=True,
)
```

Relevant overlap types include `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, `ntVJ`,
`VJ`, and `VJlen`.

### Prefiltering

```py
prefiltered = rsde.prefilter(
    count_table,
    min_samples=3,
    min_count=2,
    min_total_count=10,
    verbose=True,
)
```

A feature passes when both conditions hold:

- at least `min_samples` numeric columns are greater than or equal to
  `min_count`;
- the sum across numeric columns is greater than or equal to
  `min_total_count`.

The result contains every input row plus `prefilter_pass`.

### Calculating statistics

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    feature_column="clonotype",
    method="mann_whitney",
    simplify=True,
    p_adjust_method="fdr_bh",
    log2fc_zero_value=100,
    presence_threshold=1,
    cpu=4,
    verbose=True,
)
```

The metadata must contain `sample_id` and `group`. Numeric count-table columns
are matched to `sample_id`; unmatched metadata rows and unrelated numeric
columns are reported and ignored by the statistical setup.

Values below `presence_threshold` are treated as zero for calculations, but
the original values remain in the returned table.

### Statistical methods

| Method | Use case | Important parameters |
|:--|:--|:--|
| `mann_whitney` | Rank-based abundance comparison | `presence_threshold` |
| `fisher` | Replicate-level presence/absence recurrence | `presence_threshold` |
| `fisher_count` | Aggregated feature counts against library totals | `sample_totals` |
| `hurdle` | Combined recurrence and positive-abundance evidence | `hurdle_combine_method`, `sample_totals`, `cpm_scale`, `pseudocount` |
| `quasi_binomial` | Overdispersed count proportion model | `sample_totals` |
| `negative_binomial` | Negative-binomial count model | `sample_totals`, `negative_binomial_alpha` |
| `permutation` | Label-permutation enrichment test | `sample_totals`, `n_permutations`, `random_state`, `cpm_scale`, `pseudocount` |

#### Fisher recurrence test

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="fisher",
    presence_threshold=2,
)
```

#### Aggregated Fisher count test

```py
library_totals = {
    "stim_1": 100_000,
    "stim_2": 95_000,
    "stim_3": 102_000,
    "control_1": 110_000,
    "control_2": 105_000,
    "control_3": 108_000,
}

statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="fisher_count",
    sample_totals=library_totals,
)
```

#### Hurdle test

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="hurdle",
    hurdle_combine_method="fisher",
    cpm_scale=1_000_000,
    pseudocount=0.5,
)
```

Set `hurdle_combine_method="max"` for the alternative conservative
combination rule.

#### Negative-binomial test

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="negative_binomial",
    negative_binomial_alpha=1.0,
    sample_totals=library_totals,
)
```

#### Permutation test

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    method="permutation",
    n_permutations=20_000,
    random_state=42,
    cpu=4,
)
```

### Two groups and multiple groups

With exactly two groups, one two-sided comparison is calculated per tested
feature. With more than two groups, each group is compared with all remaining
samples using one-sided enrichment tests.

`simplify=True` retains the lowest-p-value group result per feature.
`simplify=False` retains every group-versus-rest result:

```py
expanded_statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    simplify=False,
)
```

### Multiple-testing correction

`p_adjust_method` is forwarded to `statsmodels.stats.multitest.multipletests`.
The default is Benjamini–Hochberg false-discovery-rate correction:

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    p_adjust_method="fdr_bh",
)
```

### Statistical result columns

`calc_statistics` inserts:

| Column | Meaning |
|:--|:--|
| `enriched_in` | Group with higher abundance or recurrence |
| `method` | Statistical method used |
| `mean_group_count` | Mean original count in the enriched group |
| `log2FC` | Log2 fold change for the enriched group |
| `p_val` | Raw p-value |
| `p_adj` | Multiple-testing-adjusted p-value |

Rows failing prefilter remain present with missing statistical values.

When the denominator of the fold change is zero, `log2fc_zero_value` is used
as the finite display value.

### Sorting

Statistics and postfiltering default to `DEFAULT_SORT_COLUMNS`:

```py
statistics = rsde.calc_statistics(
    prefiltered,
    samples_df,
    sort=["p_adj", "log2FC"],
)
```

Use `sort=False` or `sort=None` to preserve row order.

### Postfiltering

```py
postfiltered = rsde.postfilter(
    statistics,
    max_p_adj=0.05,
    max_p_val=None,
    min_group_mean=2,
    min_logfc=1,
    groups=None,
    groups_exclude=None,
    verbose=True,
)
```

P-value maxima use strict comparisons; count and fold-change minima are
inclusive. If `prefilter_pass` exists, only prefilter-passing rows can pass
postfiltering.

Filter to passing rows when a reduced table is needed:

```py
significant = postfiltered[postfiltered["postfilter_pass"]].copy()
```

Restrict or exclude enrichment groups:

```py
stimulated_only = rsde.postfilter(
    statistics,
    groups="stimulated",
)

exclude_control = rsde.postfilter(
    statistics,
    groups_exclude="control",
)
```

`groups` takes precedence over `groups_exclude` when both are supplied.

### Procedural chain pairing

Calculate each chain separately, then pass the two statistics or postfiltered
tables to `pair_chains`:

```py
pairing_matrix = rsde.pair_chains(
    tra_postfiltered,
    trb_postfiltered,
    paired_samples_df,
    method="jsd",
)
```

Supported methods are:

- `jsd`: Jensen–Shannon divergence; smaller values indicate more similar
  across-sample profiles.
- `pearson`: Pearson correlation; larger values indicate more similar
  profiles.

The matrix columns are features from the first table and rows are features
from the second table.

The procedural function requires globally unique, chain-specific `sample_id`
values and a `sample` column linking biological pairs.

Select requested features while retaining their best-scoring partners:

```py
pairing_matrix = rsde.pair_chains(
    tra_postfiltered,
    trb_postfiltered,
    paired_samples_df,
    method="jsd",
    filter_ids1=["TRA_feature_1", "TRA_feature_2"],
    filter_ids2=["TRB_feature_5"],
)
```

Use `feature_column` when the feature identifier is not the first column:

```py
pairing_matrix = rsde.pair_chains(
    tra_postfiltered,
    trb_postfiltered,
    paired_samples_df,
    feature_column=("tra_feature", "trb_feature"),
)
```

## Plotting procedural results

### Volcano plot

```py
fig = rsplot.de_volcano(
    statistics,
    p_column="p_adj",
    group_palette=None,
    size_range=(20, 300),
    log_sizes=True,
)
```

The default plot uses:

- `log2FC` on the horizontal axis;
- `-log10(p_adj)` or `-log10(p_val)` on the vertical axis;
- `enriched_in` for point color;
- mean group count for point size.

Use mean count on the vertical axis and significance for point size:

```py
fig = rsplot.de_volcano(
    statistics,
    p_column="p_val",
    by_mean_count=True,
    log_sizes=True,
)
```

In this mode `log_sizes=True` applies `log1p` to the mean-count vertical axis.
Use `log_sizes=False` for raw mean counts. Point sizes remain proportional to
the selected p-value column's `-log10` values.

### Differential-enrichment heatmap

```py
fig = rsplot.de_heatmap(
    postfiltered,
    samples_df,
    feature_column="clonotype",
    log_values=True,
    show_values=True,
)
```

If `postfilter_pass` or `prefilter_pass` is present, only passing rows are
displayed.

### Pairing heatmap

```py
fig = rsplot.de_pairing(
    pairing_matrix,
    log_minus=True,
    hclust=False,
    show_values=True,
)
```

For JSD matrices, `log_minus=True` colors cells by `-log10(score)` while
displaying original scores. Zeros are placed one decade below the smallest
positive value before transformation. Pearson matrices normally use
`log_minus=False` because correlations can be negative.

## Reproducible workflow checklist

1. Verify that sample IDs match count-table columns.
2. Use counts, not frequencies, for count-aware methods.
3. Record the clone filter, overlap type, and clustering settings.
4. Set and report prefilter thresholds.
5. Select a statistical method appropriate for the study design.
6. Record `presence_threshold`, sample totals, and method-specific parameters.
7. Apply multiple-testing correction.
8. Report postfilter thresholds separately from statistical testing.
9. Preserve the full marked table as well as the passing subset.
10. Set `random_state` for reproducible permutation tests.

## Common errors

### No features pass prefilter

Inspect feature abundance and reduce `min_samples`, `min_count`, or
`min_total_count` only when scientifically justified.

### Sample columns are not detected

Check that count-table sample columns are numeric and exactly match metadata
`sample_id` values.

### A count-aware method rejects totals

Provide positive `sample_totals` for every represented sample. For
`fisher_count`, feature counts and totals must be integer-compatible and a
feature count cannot exceed its sample total.

### Pairing fails

Confirm that:

- both chain tables contain numeric sample columns;
- metadata has `sample_id` and `sample`;
- each biological sample has exactly one represented sample ID per chain;
- both chains contain the same biological samples;
- feature identifiers are unique among eligible rows.

### A plot contains no rows

Check missing values and the `prefilter_pass` or `postfilter_pass` columns.
`de_heatmap` honors available pass columns, while `Analyzer.plot_volcano`
prefers the postfiltered table when it exists.

## API reference

### Stateful analyzer

::: repseq.diff_enrichment.Analyzer
    options:
      show_root_heading: true
      members:
        - read_samples
        - update_parameters
        - get_parameters
        - select_chain
        - run_count_table
        - run_prefilter
        - run_statistics
        - run_postfilter
        - pair_chains
        - run
        - plot_volcano
        - plot_heatmap
        - plot_pairing

### Procedural functions

::: repseq.diff_enrichment.prefilter
    options:
      show_root_heading: true

::: repseq.diff_enrichment.calc_statistics
    options:
      show_root_heading: true

::: repseq.diff_enrichment.postfilter
    options:
      show_root_heading: true

::: repseq.diff_enrichment.pair_chains
    options:
      show_root_heading: true

### Plotting functions

::: repseq.plot.de_volcano
    options:
      show_root_heading: true

::: repseq.plot.de_heatmap
    options:
      show_root_heading: true

::: repseq.plot.de_pairing
    options:
      show_root_heading: true
