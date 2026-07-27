# Usage: beta diversity

The `repseq.beta` module calculates pairwise repertoire similarity and distance
matrices. It creates the computationally expensive clonotype intersection table
once and reuses it for all requested metrics.

## Count-first full table

`beta.metrics` always calls `intersect_clones_in_samples_batch` in count mode.
For every sample pair, the full table contains the union of clonotypes and these
columns:

- `sample1_count`, `sample2_count`: raw counts after filtering and pooling;
- `sample1_freq`, `sample2_freq`: count divided by the total count of that
  sample in the pair; each frequency column sums to 1 within a pair;
- `sample1`, `sample2`, `pair`: pair identifiers.

The legacy `by_freq` argument is accepted for compatibility but ignored.

```python
from repseq import beta

all_results = beta.metrics(
    clonosets_df,
    overlap_type="aaV",
    metrics=None,
    cpu=4,
)
full_table = all_results["full_table"]
```

A single metric name returns one `DataFrame`. A list, or `metrics=None`, returns
a dictionary of metric matrices. `full_table` is also a valid metric name.
Metric names accept spaces or hyphens as aliases for underscores.

```python
jaccard = beta.metrics(clonosets_df, metrics="jaccard")
selected = beta.metrics(
    clonosets_df,
    metrics=["f2", "bray-curtis", "kl-divergence", "full_table"],
)
```

Use `metrics_from_table` to avoid reading and intersecting clonosets again:

```python
reused = beta.metrics_from_table(
    full_table,
    metrics=["pearson", "morisita_horn", "hellinger"],
)
```

## Plotting beta-diversity metrics

`rsplot.beta_metric` accepts either the dictionary returned by `beta.metrics`
or a metric matrix directly. Select a dictionary entry with `metric`. When no
metric is selected, the function warns, prints the available keys, and returns
without plotting. A direct matrix is always plotted as-is and the `metric`
argument is ignored.

```python
from repseq import plot as rsplot

rsplot.beta_metric(all_results, metric="f2")
rsplot.beta_metric(jaccard)
```

The heatmap uses the default R `pheatmap` blue-to-red `RdYlBu` palette and
clusters rows and columns by default. Set `hclust=False` to preserve the matrix
order. Up to three metadata columns can annotate samples.

```python
rsplot.beta_metric(
    all_results,
    metric="bray_curtis",
    metadata=metadata,
    group=["experimental_group", "tissue"],
    hclust=True,
)
```

Set `ignore_diagonal=True` to replace cells comparing a sample with itself by
`NA`. With `log_values=True`, zero values are replaced by one tenth of the
smallest positive matrix value before applying log10. `show_values=True` writes
compact original values in the cells, including `NA`; color still represents
the transformed value.

## Plotting full beta tables

`rsplot.beta_table` accepts a beta-results dictionary and automatically selects
`full_table`, or accepts the full table directly.

With `plot_type="dots"`, each pair is a scatterplot of clonotype frequencies.
Dots have black borders and 0.5 opacity, and a grey dashed identity line is
drawn behind them. Set `log_scale=True` to use logarithmic axes; zero
frequencies are placed below the smallest positive frequency using `log_base`.
For comparisons within one sample set, pair plots occupy the lower triangle and
F2 values occupy the upper triangle. Comparisons between two sample sets use a
complete row-by-column tile matrix.

```python
rsplot.beta_table(all_results, plot_type="dots", log_scale=True)
```

With `plot_type="diff"`, each pair is a cumulative-frequency matching plot.
The `top=20` clonotypes with the largest mean frequency receive the existing
20-color palette, while all remaining clonotypes are grey. Within-set pairs use
a wrapped facet layout; comparisons between two sample sets use a tile matrix.

```python
rsplot.beta_table(all_results, plot_type="diff", top=20)
```

## Notation and normalization

For raw count vectors $x$ and $y$ over the pairwise clonotype union, define:

- $N_x=\sum_i x_i$, $N_y=\sum_i y_i$;
- $p_i=x_i/N_x$, $q_i=y_i/N_y$;
- $S_x=\sum_i I(x_i>0)$, $S_y=\sum_i I(y_i>0)$;
- $S_{xy}=\sum_i I(x_i>0)I(y_i>0)$.

Presence metrics use raw counts only to determine whether a clonotype is
present. Frequency-based metrics receive counts and normalize them internally
to $p$ and $q$.

## Metrics and formulas

| Metric | Description and formula | Input semantics |
|---|---|---|
| `number_of_intersecting_clonotypes` | Number of shared clonotypes: $S_{xy}$. | Presence from counts |
| `relative_diversity` | Shared diversity relative to the product of repertoire diversities: $S_{xy}/(S_xS_y)$. | Presence from counts |
| `pearson` | Pearson correlation $r(p_i,q_i)$, calculated over clonotypes present in both samples. Returns `NaN` for fewer than two shared clonotypes or a constant vector. | Internally normalized frequencies |
| `f1` | Geometric mean of total shared-clonotype frequencies: $\sqrt{(\sum_{i:x_i y_i>0}p_i)(\sum_{i:x_i y_i>0}q_i)}$. | Internally normalized frequencies |
| `f2` | Sum of clonotype-wise geometric mean frequencies: $\sum_i\sqrt{p_iq_i}$. | Internally normalized frequencies |
| `jaccard` | Shared presence divided by union presence: $S_{xy}/(S_x+S_y-S_{xy})$. | Presence from counts |
| `jaccard_distance` | $1-\mathrm{Jaccard}$. | Presence from counts |
| `dice` | Sørensen–Dice coefficient: $2S_{xy}/(S_x+S_y)$. | Presence from counts |
| `dice_distance` | $1-\mathrm{Dice}$. | Presence from counts |
| `szymkiewicz_simpson` | Overlap coefficient: $S_{xy}/\min(S_x,S_y)$. | Presence from counts |
| `bray_curtis` | Bray–Curtis dissimilarity: $\sum_i|p_i-q_i|/\sum_i(p_i+q_i)$. For normalized vectors this equals total variation distance. | Internally normalized frequencies |
| `l1` | Manhattan distance: $\sum_i|p_i-q_i|$. | Internally normalized frequencies |
| `total_variation` | $\frac12\sum_i|p_i-q_i|$. | Internally normalized frequencies |
| `l2` | Euclidean distance: $\sqrt{\sum_i(p_i-q_i)^2}$. | Internally normalized frequencies |
| `morisita_horn` | Morisita–Horn similarity: $2\sum_i p_iq_i/(\sum_i p_i^2+\sum_i q_i^2)$. | Internally normalized frequencies |
| `jensen_shannon` | Jensen–Shannon divergence: $\frac12 KL(p\|m)+\frac12 KL(q\|m)$, $m=(p+q)/2$. | Internally normalized frequencies |
| `kl_divergence` | Directional Kullback–Leibler divergence: $KL(p\|q)=\sum_i p_i\log(p_i/q_i)$. A small epsilon protects zero denominators. | Internally normalized frequencies |
| `hellinger` | Hellinger distance: $\frac1{\sqrt2}\sqrt{\sum_i(\sqrt{p_i}-\sqrt{q_i})^2}$. | Internally normalized frequencies |
| `full_table` | Original pairwise union table with counts and frequencies. | Not a scalar metric |

`kl_divergence` is asymmetric, so `matrix.loc[a, b]` can differ from
`matrix.loc[b, a]`. Other metric matrices are symmetric for within-table
comparisons.

## Overlap types

Sequence-based overlap types are `aa`, `aaV`, `aaVJ`, `nt`, `ntV`, and
`ntVJ`. `VJ` pools by the V/J pair without sequence identity, while `VJlen`
pools by V/J and amino-acid CDR3 length.
