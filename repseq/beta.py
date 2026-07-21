"""Pairwise beta-diversity metrics for immune repertoires."""

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


def metrics(clonosets_df, cl_filter=None, overlap_type="aaV", by_freq=True,
            clonosets_df2=None, cl_filter2=None, metrics=None, cpu=None):
    """Calculate beta-diversity matrices after building the full pair table."""
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
    """Calculate beta-diversity matrices from an existing full intersection table."""
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
