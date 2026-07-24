r"""Pairwise beta-diversity metrics for immune repertoires.

All metric functions receive two raw clonotype-count vectors, ``x`` and ``y``.
Let ``N_x = sum(x)``, ``N_y = sum(y)``, ``p_i = x_i / N_x``, and
``q_i = y_i / N_y``. Presence is represented by ``I(x_i > 0)``. Define
``S_x`` and ``S_y`` as the numbers of present clonotypes and ``S_xy`` as the
number present in both samples.

Supported metrics:

- ``number_of_intersecting_clonotypes``: ``S_xy``.
- ``relative_diversity``: ``S_xy / (S_x * S_y)``.
- ``pearson``: Pearson correlation of ``p`` and ``q`` over shared clonotypes.
- ``f1``: ``sqrt(sum_shared(p) * sum_shared(q))``.
- ``f2``: ``sum_i sqrt(p_i * q_i)``.
- ``jaccard``: ``S_xy / (S_x + S_y - S_xy)``.
- ``jaccard_distance``: ``1 - jaccard``.
- ``dice``: ``2 * S_xy / (S_x + S_y)``.
- ``dice_distance``: ``1 - dice``.
- ``szymkiewicz_simpson``: ``S_xy / min(S_x, S_y)``.
- ``bray_curtis``: ``sum_i |p_i - q_i| / sum_i (p_i + q_i)``.
- ``l1``: ``sum_i |p_i - q_i|``.
- ``total_variation``: ``0.5 * sum_i |p_i - q_i|``.
- ``l2``: ``sqrt(sum_i (p_i - q_i)^2)``.
- ``morisita_horn``: ``2 * sum_i p_i*q_i / (sum_i p_i^2 + sum_i q_i^2)``.
- ``jensen_shannon``: ``0.5*KL(p || m) + 0.5*KL(q || m)``, where
  ``m = (p + q) / 2``.
- ``kl_divergence``: ``sum_i p_i * log(p_i / q_i)``. This metric is
  directional, so opposite matrix cells may differ.
- ``hellinger``: ``sqrt(sum_i (sqrt(p_i) - sqrt(q_i))^2) / sqrt(2)``.
- ``full_table``: the count-first pairwise union table used for calculation.

Presence metrics use the raw count vectors. Frequency-based metrics normalize
those vectors internally; callers do not need to pre-normalize them.
"""

from collections import OrderedDict

import numpy as np
import pandas as pd

from .common_functions import (
    bray_curtis_dissimilarity,
    dice_distance,
    dice_similarity,
    f1_similarity,
    f2_similarity,
    hellinger_distance,
    intersecting_clonotypes_count,
    jaccard_distance,
    jaccard_similarity,
    jensen_shannon_divergence,
    kl_divergence,
    l1_distance,
    l2_distance,
    morisita_horn_similarity,
    pearson_correlation,
    relative_diversity,
    szymkiewicz_simpson_similarity,
    total_variation_distance,
)
from .intersections import intersect_clones_in_samples_batch


METRICS = OrderedDict([
    ("number_of_intersecting_clonotypes", intersecting_clonotypes_count),
    ("relative_diversity", relative_diversity),
    ("pearson", pearson_correlation),
    ("f1", f1_similarity),
    ("f2", f2_similarity),
    ("jaccard", jaccard_similarity),
    ("jaccard_distance", jaccard_distance),
    ("dice", dice_similarity),
    ("dice_distance", dice_distance),
    ("szymkiewicz_simpson", szymkiewicz_simpson_similarity),
    ("bray_curtis", bray_curtis_dissimilarity),
    ("l1", l1_distance),
    ("total_variation", total_variation_distance),
    ("l2", l2_distance),
    ("morisita_horn", morisita_horn_similarity),
    ("jensen_shannon", jensen_shannon_divergence),
    ("kl_divergence", kl_divergence),
    ("hellinger", hellinger_distance),
])

_ALIASES = {
    "intersection": "number_of_intersecting_clonotypes",
    "intersecting_clonotypes": "number_of_intersecting_clonotypes",
    "number_intersecting_clonotypes": "number_of_intersecting_clonotypes",
    "pearson_correlation": "pearson",
    "szymkiewicz_simpson_coefficient": "szymkiewicz_simpson",
    "simpson_overlap": "szymkiewicz_simpson",
    "braycurtis": "bray_curtis",
    "manhattan": "l1",
    "l1_distance": "l1",
    "total_variation_distance": "total_variation",
    "euclidean": "l2",
    "l2_distance": "l2",
    "morista_horn": "morisita_horn",
    "morisita_horn_coefficient": "morisita_horn",
    "jensen_shannon_divergence": "jensen_shannon",
    "jsd": "jensen_shannon",
    "kl": "kl_divergence",
    "kullback_leibler": "kl_divergence",
    "hellinger_distance": "hellinger",
}


def metrics(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=None,
            clonosets_df2=None, cl_filter2=None, metrics=None, cpu=None):
    """Calculate beta-diversity matrices from clonoset files.

    The function first builds a full pairwise union table containing raw counts
    and derived frequencies, then calculates the requested matrices from its
    count columns.

    Args:
        clonosets_df (pd.DataFrame): First table of sample IDs and filenames.
        cl_filter (Filter, optional): Filter for the first sample table.
        overlap_type (str): Clonotype identity definition.
        by_freq (bool, optional): Deprecated compatibility argument. Ignored;
            intersections always use counts.
        clonosets_df2 (pd.DataFrame, optional): Optional second sample table.
        cl_filter2 (Filter, optional): Filter for the second sample table.
        metrics (str or list[str], optional): Metric name or names. ``None``
            calculates every metric and returns ``full_table`` as well.
        cpu (int, optional): Number of pair-calculation worker processes.

    Returns:
        pd.DataFrame or dict[str, pd.DataFrame]: One matrix for a single metric,
        otherwise a dictionary keyed by canonical metric name.
    """
    full_table = intersect_clones_in_samples_batch(
        clonosets_df,
        cl_filter=cl_filter,
        overlap_type=overlap_type,
        by_freq=by_freq,
        clonosets_df2=clonosets_df2,
        cl_filter2=cl_filter2,
        cpu=cpu,
    )
    full_table.attrs["sample_list"] = list(clonosets_df.sort_values("sample_id").sample_id)
    full_table.attrs["sample_list2"] = (
        list(clonosets_df2.sort_values("sample_id").sample_id)
        if isinstance(clonosets_df2, pd.DataFrame) else None
    )
    return metrics_from_table(full_table, metrics=metrics)


def metrics_from_table(full_table, metrics=None):
    """Calculate beta-diversity matrices from an existing full table.

    ``full_table`` must contain ``sample1_count``, ``sample2_count``,
    ``sample1``, and ``sample2``. Frequency columns are retained in returned
    ``full_table`` output but metrics always start from raw count vectors and
    normalize internally where their definitions require frequencies.

    Args:
        full_table (pd.DataFrame): Output of
            :func:`intersect_clones_in_samples_batch`.
        metrics (str or list[str], optional): Metric name or names. ``None``
            calculates all metrics and includes the original table.

    Returns:
        pd.DataFrame or dict[str, pd.DataFrame]: One matrix for a single metric,
        otherwise a dictionary keyed by canonical metric name.
    """
    required = {"sample1", "sample2", "sample1_count", "sample2_count"}
    missing = required.difference(full_table.columns)
    if missing:
        raise ValueError(f"full_table is missing required columns: {', '.join(sorted(missing))}")

    metric_names, single = _normalize_metrics(metrics)
    result = OrderedDict()
    for metric_name in metric_names:
        if metric_name == "full_table":
            result[metric_name] = full_table
        else:
            result[metric_name] = _metric_matrix(full_table, metric_name)
    return next(iter(result.values())) if single else dict(result)


def _normalize_metrics(metrics):
    if metrics is None:
        return [*METRICS, "full_table"], False
    single = isinstance(metrics, str)
    values = [metrics] if single else list(metrics)
    if not values:
        raise ValueError("metrics must contain at least one metric name")
    normalized = []
    for value in values:
        if not isinstance(value, str):
            raise TypeError("Metric names must be strings")
        name = value.strip().lower().replace("–", "_").replace("-", "_").replace(" ", "_")
        while "__" in name:
            name = name.replace("__", "_")
        name = _ALIASES.get(name, name)
        if name not in METRICS and name != "full_table":
            possible = ", ".join([*METRICS, "full_table"])
            raise ValueError(f"Unknown beta-diversity metric '{value}'. Possible values: {possible}")
        if name not in normalized:
            normalized.append(name)
    return normalized, single or len(normalized) == 1


def _metric_matrix(full_table, metric_name):
    row_samples, column_samples = _matrix_samples(full_table)
    matrix = pd.DataFrame(np.nan, index=row_samples, columns=column_samples, dtype=float)
    function = METRICS[metric_name]

    for (_, _), pair in full_table.groupby(["sample1", "sample2"], sort=False):
        sample1 = pair["sample1"].iloc[0]
        sample2 = pair["sample2"].iloc[0]
        values1 = pair["sample1_count"].to_numpy(dtype=float)
        values2 = pair["sample2_count"].to_numpy(dtype=float)
        if sample1 in matrix.index and sample2 in matrix.columns:
            matrix.loc[sample1, sample2] = function(values1, values2)
        if sample2 in matrix.index and sample1 in matrix.columns:
            reverse = function(values2, values1) if metric_name == "kl_divergence" else matrix.loc[sample1, sample2]
            matrix.loc[sample2, sample1] = reverse

    if row_samples == column_samples:
        sample_vectors = _sample_vectors(full_table)
        for sample in row_samples:
            if sample in sample_vectors:
                values = sample_vectors[sample]
                matrix.loc[sample, sample] = function(values, values)
    matrix.index.name = "sample1"
    matrix.columns.name = "sample2"
    return matrix


def _matrix_samples(full_table):
    sample_list = full_table.attrs.get("sample_list")
    sample_list2 = full_table.attrs.get("sample_list2")
    if sample_list is not None:
        return list(sample_list), list(sample_list2) if sample_list2 is not None else list(sample_list)
    samples = list(dict.fromkeys([*full_table["sample1"], *full_table["sample2"]]))
    return samples, samples


def _sample_vectors(full_table):
    vectors = {}
    for (_, _), pair in full_table.groupby(["sample1", "sample2"], sort=False):
        sample1 = pair["sample1"].iloc[0]
        sample2 = pair["sample2"].iloc[0]
        vectors.setdefault(sample1, pair["sample1_count"].to_numpy(dtype=float))
        vectors.setdefault(sample2, pair["sample2_count"].to_numpy(dtype=float))
    return vectors


__all__ = ["METRICS", "metrics", "metrics_from_table"]
