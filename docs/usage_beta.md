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
