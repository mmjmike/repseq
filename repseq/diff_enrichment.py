"""Utilities for differential enrichment analysis."""

from numbers import Integral, Real

import pandas as pd


PREFILTER_COLUMN = "prefilter_pass"


def prefilter(
    count_table,
    min_samples=3,
    min_count=2,
    min_total_count=8,
    verbose=True,
):
    """Mark features that pass count-based abundance thresholds.

    Filtering is calculated from all numeric columns. A row passes when its
    numeric values sum to at least ``min_total_count`` and at least
    ``min_samples`` numeric columns contain a value of at least ``min_count``.
    Non-numeric columns are preserved but ignored during filtering.

    Parameters
    ----------
    count_table : pandas.DataFrame
        A wide feature-by-sample table, such as the output of
        :func:`repseq.intersections.count_table`.
    min_samples : int, default 3
        Minimum number of numeric columns with values greater than or equal to
        ``min_count``.
    min_count : real number, default 2
        Minimum value for a numeric column to count toward ``min_samples``.
    min_total_count : real number, default 8
        Minimum sum across all numeric columns.
    verbose : bool, default True
        Print the thresholds and the number of passing features.

    Returns
    -------
    pandas.DataFrame
        A copy of ``count_table`` with a boolean ``prefilter_pass`` column
        inserted immediately before the first numeric column.
    """
    _validate_prefilter_arguments(
        count_table,
        min_samples=min_samples,
        min_count=min_count,
        min_total_count=min_total_count,
    )

    numeric_columns = list(count_table.select_dtypes(include="number").columns)
    numeric_values = count_table[numeric_columns]
    sample_threshold_pass = numeric_values.ge(min_count).sum(axis=1) >= min_samples
    total_threshold_pass = numeric_values.sum(axis=1) >= min_total_count
    prefilter_pass = sample_threshold_pass & total_threshold_pass

    result = count_table.copy()
    insert_position = (
        result.columns.get_loc(numeric_columns[0])
        if numeric_columns
        else len(result.columns)
    )
    result.insert(insert_position, PREFILTER_COLUMN, prefilter_pass.astype(bool))

    if verbose:
        passed_features = int(prefilter_pass.sum())
        print("Differential enrichment prefilter\n" + "-" * 50)
        print(f"Minimum samples: {min_samples}")
        print(f"Minimum count per sample: {min_count}")
        print(f"Minimum total count: {min_total_count}")
        print(f"Features passed: {passed_features} of {len(count_table)}")

    return result


def _validate_prefilter_arguments(
    count_table,
    min_samples,
    min_count,
    min_total_count,
):
    if not isinstance(count_table, pd.DataFrame):
        raise TypeError("count_table must be a pandas DataFrame")
    if PREFILTER_COLUMN in count_table.columns:
        raise ValueError(
            f"count_table already contains the reserved column {PREFILTER_COLUMN!r}"
        )
    if (
        not isinstance(min_samples, Integral)
        or isinstance(min_samples, bool)
        or min_samples < 0
    ):
        raise ValueError("min_samples must be a non-negative integer")
    for name, value in (
        ("min_count", min_count),
        ("min_total_count", min_total_count),
    ):
        if not isinstance(value, Real) or isinstance(value, bool) or value < 0:
            raise ValueError(f"{name} must be a non-negative number")
