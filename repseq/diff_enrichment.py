"""Utilities for differential enrichment analysis."""

from numbers import Integral, Real

import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from statsmodels.stats.multitest import multipletests

from .common_functions import run_parallel_calculation


PREFILTER_COLUMN = "prefilter_pass"
POSTFILTER_COLUMN = "postfilter_pass"
STATISTICS_COLUMNS = [
    "enriched_in",
    "method",
    "mean_group_count",
    "log2FC",
    "p_val",
    "p_adj",
]
SUPPORTED_METHODS = {
    "mann_whitney",
    "fisher",
    "fisher_count",
    "hurdle",
    "quasi_binomial",
    "negative_binomial",
    "permutation",
}
PAIR_CHAIN_METHOD_ALIASES = {
    "pearson": "pearson",
    "pearson_correlation": "pearson",
    "correlation": "pearson",
    "jsd": "jsd",
    "jensen_shannon": "jsd",
    "jensen-shannon": "jsd",
    "jensen_shannon_divergence": "jsd",
    "jenson_shannon": "jsd",
    "jenson-shannon": "jsd",
}


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


def postfilter(
    statistics_table,
    max_p_adj=1,
    max_p_val=1,
    min_group_mean=2,
    min_logfc=1,
    groups=None,
):
    """Mark differential-enrichment rows that pass result thresholds.

    All thresholds are inclusive. Rows pass when ``p_adj <= max_p_adj``,
    ``p_val <= max_p_val``, ``mean_group_count >= min_group_mean``, and
    ``log2FC >= min_logfc``. When ``groups`` is provided, ``enriched_in`` must
    also match one of the requested groups.

    Parameters
    ----------
    statistics_table : pandas.DataFrame
        A table returned by :func:`calc_statistics`.
    max_p_adj, max_p_val : real number, default 1
        Inclusive maximum adjusted and raw p-values.
    min_group_mean : real number, default 2
        Inclusive minimum ``mean_group_count``.
    min_logfc : real number, default 1
        Inclusive minimum ``log2FC``.
    groups : str or sequence of str, optional
        Allowed ``enriched_in`` groups. ``None`` disables group filtering.

    Returns
    -------
    pandas.DataFrame
        A copy of ``statistics_table`` with boolean ``postfilter_pass`` added
        immediately before the first numeric column.
    """
    selected_groups = _validate_postfilter_arguments(
        statistics_table,
        max_p_adj=max_p_adj,
        max_p_val=max_p_val,
        min_group_mean=min_group_mean,
        min_logfc=min_logfc,
        groups=groups,
    )
    postfilter_pass = (
        statistics_table["p_adj"].le(max_p_adj)
        & statistics_table["p_val"].le(max_p_val)
        & statistics_table["mean_group_count"].ge(min_group_mean)
        & statistics_table["log2FC"].ge(min_logfc)
    )
    if selected_groups is not None:
        postfilter_pass &= statistics_table["enriched_in"].isin(selected_groups)
    postfilter_pass = postfilter_pass.fillna(False).astype(bool)

    result = statistics_table.copy()
    numeric_columns = list(result.select_dtypes(include="number").columns)
    insert_position = (
        result.columns.get_loc(numeric_columns[0])
        if numeric_columns
        else len(result.columns)
    )
    result.insert(insert_position, POSTFILTER_COLUMN, postfilter_pass)
    result.attrs = statistics_table.attrs.copy()
    return result


def _validate_postfilter_arguments(
    statistics_table,
    max_p_adj,
    max_p_val,
    min_group_mean,
    min_logfc,
    groups,
):
    if not isinstance(statistics_table, pd.DataFrame):
        raise TypeError("statistics_table must be a pandas DataFrame")
    if POSTFILTER_COLUMN in statistics_table.columns:
        raise ValueError(
            "statistics_table already contains the reserved column "
            f"{POSTFILTER_COLUMN!r}"
        )
    required_columns = {
        "p_adj",
        "p_val",
        "mean_group_count",
        "log2FC",
        "enriched_in",
    }
    missing_columns = required_columns.difference(statistics_table.columns)
    if missing_columns:
        raise ValueError(
            "statistics_table is missing required postfilter columns: "
            f"{sorted(missing_columns)}"
        )
    for name, value in (("max_p_adj", max_p_adj), ("max_p_val", max_p_val)):
        if (
            not isinstance(value, Real)
            or isinstance(value, bool)
            or value < 0
            or value > 1
        ):
            raise ValueError(f"{name} must be a number in [0, 1]")
    if (
        not isinstance(min_group_mean, Real)
        or isinstance(min_group_mean, bool)
        or min_group_mean < 0
    ):
        raise ValueError("min_group_mean must be a non-negative number")
    if not isinstance(min_logfc, Real) or isinstance(min_logfc, bool):
        raise ValueError("min_logfc must be a number")
    if groups is None:
        return None
    if isinstance(groups, str):
        return [groups]
    if not isinstance(groups, (list, tuple, set, pd.Index, np.ndarray)):
        raise TypeError("groups must be None, a string, or a sequence of strings")
    selected_groups = list(groups)
    if not selected_groups:
        raise ValueError("groups must not be empty")
    if any(not isinstance(group, str) for group in selected_groups):
        raise TypeError("groups must contain only strings")
    return list(dict.fromkeys(selected_groups))


def calc_statistics(
    count_table,
    samples_metadata,
    feature_column=None,
    method="mann_whitney",
    simplify=True,
    p_adjust_method="fdr_bh",
    log2fc_zero_value=100,
    presence_threshold=1,
    sample_totals=None,
    hurdle_combine_method="fisher",
    cpm_scale=1_000_000,
    pseudocount=0.5,
    n_permutations=10_000,
    random_state=None,
    negative_binomial_alpha=1.0,
    cpu=None,
    verbose=True,
):
    """Calculate differential enrichment statistics for count-table features.

    With two groups, one two-sided comparison is calculated. With more than
    two groups, each group is compared with all remaining samples using a
    one-sided enrichment test. If ``prefilter_pass`` is present, only rows
    marked ``True`` are tested; all rows remain in the returned table.

    Supported methods are ``mann_whitney``, ``fisher`` (replicate recurrence),
    ``fisher_count`` (aggregated counts), ``hurdle``, ``quasi_binomial``,
    ``negative_binomial``, and ``permutation``. Parameters specific to methods
    other than the selected method are ignored.

    Parameters
    ----------
    count_table : pandas.DataFrame
        Wide feature-by-sample count table.
    samples_metadata : pandas.DataFrame
        Table containing unique ``sample_id`` values and their ``group``.
    feature_column : hashable, optional
        Unique feature identifier column. Defaults to the first column.
    method : str, default "mann_whitney"
        Statistical method to use.
    simplify : bool, default True
        Keep only the lowest-p-value group result for each feature. If false,
        retain one row per tested group comparison.
    p_adjust_method : str, default "fdr_bh"
        Multiple-testing method accepted by
        :func:`statsmodels.stats.multitest.multipletests`.
    log2fc_zero_value : real number, default 100
        Absolute finite log2-fold-change used when exactly one mean is zero.
    presence_threshold : real number, default 1
        Counts below this threshold are treated as zero during statistical
        calculations. Original count values are preserved in the output.
    sample_totals : mapping or pandas.Series, optional
        Per-sample library totals used by count-aware methods. By default,
        totals are calculated from all rows of ``count_table``.
    hurdle_combine_method : {"fisher", "max"}, default "fisher"
        How ``hurdle`` combines presence and positive-abundance p-values.
    cpm_scale : positive real number, default 1_000_000
        Scale used for normalized abundance in ``hurdle`` and ``permutation``.
    pseudocount : positive real number, default 0.5
        Pseudocount used before log2 transformation.
    n_permutations : positive int, default 10_000
        Label permutations used by ``permutation``.
    random_state : int, optional
        Random seed used by ``permutation``.
    negative_binomial_alpha : positive real number, default 1.0
        Fixed dispersion parameter used by ``negative_binomial``.
    cpu : int, optional
        Worker count passed to ``run_parallel_calculation``.
    verbose : bool, default True
        Print analysis setup and parallel progress.

    Returns
    -------
    pandas.DataFrame
        The input table with ``enriched_in``, ``method``,
        ``mean_group_count``, ``log2FC``, ``p_val``, and ``p_adj`` inserted
        before its numeric columns.
    """
    setup = _prepare_statistics_setup(
        count_table,
        samples_metadata,
        feature_column=feature_column,
        method=method,
    )
    method = setup["method"]
    _validate_statistics_parameters(
        method=method,
        log2fc_zero_value=log2fc_zero_value,
        presence_threshold=presence_threshold,
        hurdle_combine_method=hurdle_combine_method,
        cpm_scale=cpm_scale,
        pseudocount=pseudocount,
        n_permutations=n_permutations,
        random_state=random_state,
        negative_binomial_alpha=negative_binomial_alpha,
    )
    _validate_p_adjust_method(p_adjust_method)
    if verbose:
        _print_statistics_setup(setup)

    sample_columns = setup["sample_columns"]
    values = count_table[sample_columns].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Sample count columns must contain only finite values")
    if (values < 0).any():
        raise ValueError("Sample count columns must contain non-negative values")
    effective_values = values.copy()
    thresholded_values = (effective_values < presence_threshold) & (
        effective_values != 0
    )
    effective_values[effective_values < presence_threshold] = 0
    if verbose:
        print(
            f"Analysis count threshold: {presence_threshold}. Values below "
            "this threshold are treated as zero for statistical calculations."
        )
        print(
            f"Non-zero count values replaced with zero: "
            f"{int(thresholded_values.sum())}. Original counts will be preserved "
            "in the output table."
        )

    count_aware_methods = {
        "fisher_count",
        "hurdle",
        "quasi_binomial",
        "negative_binomial",
        "permutation",
    }
    if method in count_aware_methods:
        library_totals = _resolve_sample_totals(
            count_table,
            sample_columns,
            sample_totals=sample_totals,
        )
        if (library_totals <= 0).any():
            invalid_samples = [
                sample_columns[index]
                for index in np.flatnonzero(library_totals <= 0)
            ]
            raise ValueError(
                "Positive sample totals are required for method "
                f"{method!r}; invalid samples: {invalid_samples}"
            )
    else:
        library_totals = np.ones(len(sample_columns), dtype=float)
    if method == "fisher_count":
        if not np.allclose(values, np.round(values)):
            raise ValueError("method='fisher_count' requires integer counts")
        if not np.allclose(library_totals, np.round(library_totals)):
            raise ValueError("method='fisher_count' requires integer sample totals")
        if (values > library_totals).any():
            raise ValueError("Feature counts cannot exceed their sample totals")

    pass_mask = setup["pass_mask"]
    passed_positions = np.flatnonzero(pass_mask)
    statistics = pd.DataFrame(
        columns=["_row_position", *STATISTICS_COLUMNS]
    )
    if len(passed_positions):
        tested_values = effective_values[passed_positions]
        original_tested_values = values[passed_positions]
        tasks = _build_statistics_tasks(
            tested_values=tested_values,
            original_tested_values=original_tested_values,
            passed_positions=passed_positions,
            groups=setup["groups"],
            group_sample_indices=setup["group_sample_indices"],
            method=method,
            library_totals=library_totals,
            log2fc_zero_value=log2fc_zero_value,
            presence_threshold=presence_threshold,
            hurdle_combine_method=hurdle_combine_method,
            cpm_scale=cpm_scale,
            pseudocount=pseudocount,
            n_permutations=n_permutations,
            random_state=random_state,
            negative_binomial_alpha=negative_binomial_alpha,
        )
        result_tables = run_parallel_calculation(
            _calculate_statistics_task,
            tasks,
            "Calculating differential enrichment",
            object_name="group comparisons",
            verbose=verbose,
            cpu=cpu,
        )
        if verbose:
            print("Combining group-comparison result tables.")
        statistics = pd.concat(result_tables, ignore_index=True)
        if verbose:
            print(
                "Adjusting p-values for multiple testing "
                f"using {p_adjust_method!r}."
            )
        statistics["p_adj"] = _adjust_p_values(
            statistics["p_val"],
            method=p_adjust_method,
        )
    elif verbose:
        print("No features passed prefiltering; statistical tests were skipped.")

    if simplify:
        if verbose:
            print(
                "Simplifying statistics by selecting the lowest-p-value "
                "comparison for each feature."
            )
        if len(statistics):
            statistics = simplify_statistics(statistics)
    elif verbose:
        print(
            "Statistics simplification is disabled; retaining all group "
            "comparisons."
        )

    if verbose:
        print("Assembling the final differential enrichment output table.")
    result = _assemble_statistics_output(
        count_table,
        statistics,
        pass_mask=pass_mask,
        simplify=simplify,
    )
    if verbose:
        print("Differential enrichment analysis finished successfully!")
    return result


def simplify_statistics(statistics):
    """Keep the lowest-p-value group result for each tested feature row."""
    ranked = statistics.assign(
        _p_sort=statistics["p_val"].fillna(np.inf)
    ).sort_values(["_row_position", "_p_sort"], kind="stable")
    return ranked.drop_duplicates("_row_position", keep="first").drop(
        columns="_p_sort"
    )


def simplify(statistics):
    """Return the lowest-p-value group result for each tested feature row."""
    return simplify_statistics(statistics)


def pair_chains(
    count_table1,
    count_table2,
    samples_metadata,
    method="jsd",
    feature_column=None,
    filter_ids1=None,
    filter_ids2=None,
):
    """Score feature pairs from two chains across paired biological samples.

    Parameters
    ----------
    count_table1, count_table2 : pandas.DataFrame
        Chain-specific tables returned by :func:`calc_statistics`.
    samples_metadata : pandas.DataFrame
        Metadata containing unique ``sample_id`` values and a ``sample`` value
        that pairs one sample ID from each chain.
    method : {"jsd", "pearson"}, default "jsd"
        Jensen–Shannon divergence or Pearson correlation across paired samples.
    feature_column : hashable or pair of hashable, optional
        Feature identifier column. A single value is used for both tables; a
        two-item tuple/list specifies table-specific columns. By default, each
        table's first column is used.
    filter_ids1, filter_ids2 : sequence, optional
        Requested feature IDs for either chain. Requested IDs are retained and
        the best-scoring partners required by IDs on the opposite axis are
        added. Without filters, all eligible features are returned.

    Returns
    -------
    pandas.DataFrame
        Score matrix with table-2 features on rows and table-1 features on
        columns.
    """
    method = _normalize_pair_chain_method(method)
    feature_column1, feature_column2 = _pair_chain_feature_columns(
        count_table1,
        count_table2,
        feature_column,
    )
    paired_samples = _prepare_paired_chain_samples(
        count_table1,
        count_table2,
        samples_metadata,
        feature_column1,
        feature_column2,
    )
    features1, values1 = _prepare_pair_chain_features(
        count_table1,
        feature_column1,
        paired_samples["sample_id1"],
        "count_table1",
    )
    features2, values2 = _prepare_pair_chain_features(
        count_table2,
        feature_column2,
        paired_samples["sample_id2"],
        "count_table2",
    )

    normalized1 = _normalize_pair_chain_rows(values1)
    normalized2 = _normalize_pair_chain_rows(values2)
    if method == "pearson":
        scores = _pairwise_pearson_scores(normalized1, normalized2)
    else:
        scores = _pairwise_jsd_scores(normalized1, normalized2)
    score_matrix = pd.DataFrame(
        scores,
        index=pd.Index(features2, name=feature_column2),
        columns=pd.Index(features1, name=feature_column1),
    )
    score_matrix = _filter_pair_chain_matrix(
        score_matrix,
        filter_ids1=filter_ids1,
        filter_ids2=filter_ids2,
        method=method,
    )
    score_matrix.attrs["method"] = method
    score_matrix.attrs["paired_samples"] = paired_samples.copy()
    return score_matrix


def _normalize_pair_chain_method(method):
    normalized = str(method).strip().casefold().replace(" ", "_")
    selected = PAIR_CHAIN_METHOD_ALIASES.get(normalized)
    if selected is None:
        raise ValueError("method must be one of: jsd, pearson")
    return selected


def _pair_chain_feature_columns(count_table1, count_table2, feature_column):
    for name, table in (
        ("count_table1", count_table1),
        ("count_table2", count_table2),
    ):
        if not isinstance(table, pd.DataFrame):
            raise TypeError(f"{name} must be a pandas DataFrame")
        if table.shape[1] == 0:
            raise ValueError(f"{name} must contain at least one column")
    if feature_column is None:
        columns = (count_table1.columns[0], count_table2.columns[0])
    elif isinstance(feature_column, (list, tuple)):
        if len(feature_column) != 2:
            raise ValueError("feature_column must contain exactly two column names")
        columns = tuple(feature_column)
    else:
        columns = (feature_column, feature_column)
    if columns[0] not in count_table1.columns:
        raise ValueError(f"feature_column {columns[0]!r} is not in count_table1")
    if columns[1] not in count_table2.columns:
        raise ValueError(f"feature_column {columns[1]!r} is not in count_table2")
    return columns


def _prepare_paired_chain_samples(
    count_table1,
    count_table2,
    samples_metadata,
    feature_column1,
    feature_column2,
):
    if not isinstance(samples_metadata, pd.DataFrame):
        raise TypeError("samples_metadata must be a pandas DataFrame")
    missing_columns = {"sample_id", "sample"}.difference(samples_metadata.columns)
    if missing_columns:
        raise ValueError(
            "samples_metadata must contain columns 'sample_id' and 'sample'; "
            f"missing: {sorted(missing_columns)}"
        )
    if samples_metadata[["sample_id", "sample"]].isna().any().any():
        raise ValueError(
            "samples_metadata columns 'sample_id' and 'sample' must not contain "
            "missing values"
        )
    duplicated_ids = samples_metadata["sample_id"].duplicated(keep=False)
    if duplicated_ids.any():
        duplicates = samples_metadata.loc[duplicated_ids, "sample_id"].tolist()
        raise ValueError(
            "samples_metadata['sample_id'] values must be unique; duplicated "
            f"values: {duplicates}"
        )

    metadata_ids = set(samples_metadata["sample_id"])
    numeric1 = set(count_table1.select_dtypes(include="number").columns)
    numeric2 = set(count_table2.select_dtypes(include="number").columns)
    sample_ids1 = [
        column
        for column in count_table1.columns
        if column in metadata_ids and column in numeric1 and column != feature_column1
    ]
    sample_ids2 = [
        column
        for column in count_table2.columns
        if column in metadata_ids and column in numeric2 and column != feature_column2
    ]
    non_numeric1 = [
        column
        for column in count_table1.columns
        if column in metadata_ids and column not in numeric1
    ]
    non_numeric2 = [
        column
        for column in count_table2.columns
        if column in metadata_ids and column not in numeric2
    ]
    if non_numeric1 or non_numeric2:
        raise ValueError(
            "Paired sample columns must be numeric; non-numeric columns: "
            f"count_table1={non_numeric1}, count_table2={non_numeric2}"
        )
    if not sample_ids1 or not sample_ids2:
        raise ValueError("Both count tables must contain metadata-matched sample columns")
    overlapping_ids = sorted(set(sample_ids1).intersection(sample_ids2), key=str)
    if overlapping_ids:
        raise ValueError(
            "The two count tables must use different chain-specific sample_id "
            f"values; shared IDs: {overlapping_ids}"
        )

    metadata = samples_metadata.set_index("sample_id")
    metadata1 = metadata.loc[sample_ids1, ["sample"]].reset_index()
    metadata2 = metadata.loc[sample_ids2, ["sample"]].reset_index()
    duplicated_samples1 = metadata1["sample"].duplicated(keep=False)
    duplicated_samples2 = metadata2["sample"].duplicated(keep=False)
    if duplicated_samples1.any() or duplicated_samples2.any():
        raise ValueError(
            "Each biological sample must have exactly one represented sample_id "
            "per count table"
        )

    samples1 = set(metadata1["sample"])
    samples2 = set(metadata2["sample"])
    if samples1 != samples2:
        raise ValueError(
            "Every represented sample in one count table must have a paired "
            "sample_id in the other table. Missing from count_table2: "
            f"{sorted(samples1 - samples2, key=str)}; missing from count_table1: "
            f"{sorted(samples2 - samples1, key=str)}"
        )
    sample_id2_by_sample = metadata2.set_index("sample")["sample_id"]
    return pd.DataFrame(
        {
            "sample": metadata1["sample"].tolist(),
            "sample_id1": metadata1["sample_id"].tolist(),
            "sample_id2": [
                sample_id2_by_sample[sample]
                for sample in metadata1["sample"]
            ],
        }
    )


def _prepare_pair_chain_features(
    count_table,
    feature_column,
    sample_columns,
    table_name,
):
    if PREFILTER_COLUMN in count_table.columns:
        valid_prefilter = count_table[PREFILTER_COLUMN].isin([True, False])
        if not valid_prefilter.all():
            raise ValueError(
                f"{table_name}[{PREFILTER_COLUMN!r}] must contain only True or False"
            )
        eligible = count_table.loc[count_table[PREFILTER_COLUMN].eq(True)].copy()
    else:
        eligible = count_table.copy()
    if eligible.empty:
        raise ValueError(f"{table_name} contains no eligible features")
    duplicated = eligible[feature_column].duplicated(keep=False)
    if duplicated.any():
        duplicate_value = eligible.loc[duplicated, feature_column].iloc[0]
        raise ValueError(
            f"{table_name} feature column {feature_column!r} contains duplicate "
            f"eligible value {duplicate_value!r}"
        )
    values = eligible[list(sample_columns)].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{table_name} paired sample counts must be finite")
    if (values < 0).any():
        raise ValueError(f"{table_name} paired sample counts must be non-negative")
    return eligible[feature_column].tolist(), values


def _normalize_pair_chain_rows(values):
    values = np.asarray(values, dtype=float)
    totals = values.sum(axis=1, keepdims=True)
    return np.divide(
        values,
        totals,
        out=np.zeros_like(values, dtype=float),
        where=totals != 0,
    )


def _pairwise_pearson_scores(values1, values2):
    centered1 = values1 - values1.mean(axis=1, keepdims=True)
    centered2 = values2 - values2.mean(axis=1, keepdims=True)
    denominator = np.outer(
        np.linalg.norm(centered2, axis=1),
        np.linalg.norm(centered1, axis=1),
    )
    return np.divide(
        centered2 @ centered1.T,
        denominator,
        out=np.full((len(values2), len(values1)), np.nan),
        where=denominator != 0,
    )


def _pairwise_jsd_scores(values1, values2):
    scores = np.empty((len(values2), len(values1)), dtype=float)
    for row_index, values in enumerate(values2):
        midpoint = 0.5 * (values1 + values)
        first_terms = np.zeros_like(values1)
        second_terms = np.zeros_like(values1)
        np.log(
            np.divide(
                values1,
                midpoint,
                out=np.ones_like(values1),
                where=midpoint != 0,
            ),
            out=first_terms,
            where=values1 != 0,
        )
        np.log(
            np.divide(
                values,
                midpoint,
                out=np.ones_like(values1),
                where=midpoint != 0,
            ),
            out=second_terms,
            where=values != 0,
        )
        scores[row_index] = 0.5 * (
            np.sum(values1 * first_terms, axis=1)
            + np.sum(values * second_terms, axis=1)
        )
    return scores


def _normalize_filter_ids(filter_ids, eligible_ids, name):
    if filter_ids is None:
        return None
    if isinstance(filter_ids, (str, bytes)):
        requested = [filter_ids]
    else:
        requested = list(filter_ids)
    if not requested:
        raise ValueError(f"{name} must not be empty")
    requested = list(dict.fromkeys(requested))
    eligible = set(eligible_ids)
    missing = [feature_id for feature_id in requested if feature_id not in eligible]
    if missing:
        raise ValueError(
            f"{name} contains IDs absent after prefiltering: {missing}"
        )
    return requested


def _best_pair_id(values, candidates, method):
    finite = np.isfinite(values)
    if not finite.any():
        return candidates[0]
    valid_positions = np.flatnonzero(finite)
    valid_values = values[finite]
    selected = (
        valid_positions[np.argmin(valid_values)]
        if method == "jsd"
        else valid_positions[np.argmax(valid_values)]
    )
    return candidates[int(selected)]


def _filter_pair_chain_matrix(
    score_matrix,
    filter_ids1,
    filter_ids2,
    method,
):
    selected1 = _normalize_filter_ids(
        filter_ids1,
        list(score_matrix.columns),
        "filter_ids1",
    )
    selected2 = _normalize_filter_ids(
        filter_ids2,
        list(score_matrix.index),
        "filter_ids2",
    )
    if selected1 is None and selected2 is None:
        return score_matrix
    selected1 = [] if selected1 is None else selected1
    selected2 = [] if selected2 is None else selected2
    requested1 = list(selected1)
    requested2 = list(selected2)

    for feature_id1 in requested1:
        best_id2 = _best_pair_id(
            score_matrix[feature_id1].to_numpy(dtype=float),
            list(score_matrix.index),
            method,
        )
        if best_id2 not in selected2:
            selected2.append(best_id2)
    for feature_id2 in requested2:
        best_id1 = _best_pair_id(
            score_matrix.loc[feature_id2].to_numpy(dtype=float),
            list(score_matrix.columns),
            method,
        )
        if best_id1 not in selected1:
            selected1.append(best_id1)
    return score_matrix.loc[selected2, selected1]


def _prepare_statistics_setup(
    count_table,
    samples_metadata,
    feature_column,
    method,
):
    if not isinstance(count_table, pd.DataFrame):
        raise TypeError("count_table must be a pandas DataFrame")
    if not isinstance(samples_metadata, pd.DataFrame):
        raise TypeError("samples_metadata must be a pandas DataFrame")
    missing_metadata_columns = {
        "sample_id",
        "group",
    }.difference(samples_metadata.columns)
    if missing_metadata_columns:
        raise ValueError(
            "samples_metadata must contain columns 'sample_id' and 'group'; "
            f"missing: {sorted(missing_metadata_columns)}"
        )
    if count_table.shape[1] == 0:
        raise ValueError("count_table must contain at least one column")

    method = str(method).casefold()
    if method not in SUPPORTED_METHODS:
        raise ValueError(
            "method must be one of: " + ", ".join(sorted(SUPPORTED_METHODS))
        )
    feature_column = count_table.columns[0] if feature_column is None else feature_column
    if feature_column not in count_table.columns:
        raise ValueError(f"feature_column {feature_column!r} is not in count_table")
    duplicated_features = count_table[feature_column].duplicated(keep=False)
    if duplicated_features.any():
        duplicate_value = count_table.loc[duplicated_features, feature_column].iloc[0]
        if pd.isna(duplicate_value):
            duplicate_indices = count_table.index[
                count_table[feature_column].isna()
            ].tolist()
        else:
            duplicate_indices = count_table.index[
                count_table[feature_column].eq(duplicate_value)
            ].tolist()
        raise ValueError(
            f"Column {feature_column!r} was used as feature_column, but its values "
            f"are non-unique. Example duplicated value {duplicate_value!r} occurs "
            f"at row indices {duplicate_indices}."
        )

    if samples_metadata["sample_id"].isna().any():
        raise ValueError("samples_metadata['sample_id'] must not contain missing values")
    if samples_metadata["group"].isna().any():
        raise ValueError("samples_metadata['group'] must not contain missing values")
    duplicated_samples = samples_metadata["sample_id"].duplicated(keep=False)
    if duplicated_samples.any():
        duplicates = samples_metadata.loc[duplicated_samples, "sample_id"].tolist()
        raise ValueError(
            "samples_metadata['sample_id'] values must be unique; duplicated values: "
            f"{duplicates}"
        )

    numeric_columns = list(count_table.select_dtypes(include="number").columns)
    numeric_sample_columns = [
        column for column in numeric_columns if column != feature_column
    ]
    metadata_sample_ids = samples_metadata["sample_id"].tolist()
    metadata_only_samples = [
        sample_id
        for sample_id in metadata_sample_ids
        if sample_id not in count_table.columns
    ]
    non_numeric_samples = [
        sample_id
        for sample_id in metadata_sample_ids
        if sample_id in count_table.columns
        and sample_id not in numeric_sample_columns
    ]
    if non_numeric_samples:
        raise ValueError(
            "Sample columns must be numeric; non-numeric samples: "
            f"{non_numeric_samples}"
        )
    sample_columns = [
        column
        for column in numeric_sample_columns
        if column in set(metadata_sample_ids)
    ]
    ignored_sample_columns = [
        column
        for column in numeric_sample_columns
        if column not in set(metadata_sample_ids)
    ]
    metadata = samples_metadata.set_index("sample_id").loc[sample_columns].reset_index()
    groups = metadata["group"].drop_duplicates().tolist()
    if len(groups) < 2:
        raise ValueError(
            "Differential enrichment requires at least 2 groups represented by "
            f"count_table samples; detected {len(groups)} group(s): {groups}"
        )
    group_samples = {
        group: metadata.loc[metadata["group"].eq(group), "sample_id"].tolist()
        for group in groups
    }
    undersized_groups = {
        group: samples
        for group, samples in group_samples.items()
        if len(samples) < 2
    }
    if undersized_groups:
        details = "; ".join(
            f"{group!r}: {len(samples)} sample(s) {samples}"
            for group, samples in undersized_groups.items()
        )
        raise ValueError(
            "Each group must contain at least 2 samples for statistical "
            f"comparison. Undersized groups: {details}"
        )

    sample_positions = {sample: index for index, sample in enumerate(sample_columns)}
    group_sample_indices = {
        group: np.array([sample_positions[sample] for sample in samples], dtype=int)
        for group, samples in group_samples.items()
    }
    if PREFILTER_COLUMN in count_table.columns:
        valid_prefilter_values = count_table[PREFILTER_COLUMN].isin([True, False])
        if not valid_prefilter_values.all():
            raise ValueError(
                f"{PREFILTER_COLUMN!r} must contain only True or False values"
            )
        pass_mask = count_table[PREFILTER_COLUMN].eq(True).to_numpy()
    else:
        pass_mask = np.ones(len(count_table), dtype=bool)

    return {
        "feature_column": feature_column,
        "method": method,
        "numeric_columns": numeric_columns,
        "sample_columns": sample_columns,
        "ignored_sample_columns": ignored_sample_columns,
        "metadata_only_samples": metadata_only_samples,
        "groups": groups,
        "group_samples": group_samples,
        "group_sample_indices": group_sample_indices,
        "pass_mask": pass_mask,
        "prefilter_used": PREFILTER_COLUMN in count_table.columns,
        "total_features": len(count_table),
    }


def _print_statistics_setup(setup):
    print("Differential enrichment analysis started\n" + "-" * 50)
    print(f"Method: {setup['method']}")
    print(f"Feature column: {setup['feature_column']}")
    print(f"Groups detected: {len(setup['groups'])}")
    for group in setup["groups"]:
        samples = setup["group_samples"][group]
        print(f"Group {group!r} ({len(samples)} samples): {samples}")
    if setup["ignored_sample_columns"]:
        print(
            "Numeric sample columns absent from samples_metadata will be ignored "
            f"({len(setup['ignored_sample_columns'])}): "
            f"{setup['ignored_sample_columns']}"
        )
    else:
        print("All numeric sample columns are represented in samples_metadata.")
    if setup["metadata_only_samples"]:
        print(
            "samples_metadata entries absent from count_table will be ignored "
            f"({len(setup['metadata_only_samples'])}): "
            f"{setup['metadata_only_samples']}"
        )
    else:
        print("All samples_metadata sample IDs are represented in count_table.")
    if setup["prefilter_used"]:
        print(
            "Prefilter applied: "
            f"{int(setup['pass_mask'].sum())} of {setup['total_features']} "
            "features will be tested."
        )
    else:
        print(
            "No prefilter_pass column detected: all "
            f"{setup['total_features']} features will be tested."
        )


def _validate_statistics_parameters(
    method,
    log2fc_zero_value,
    presence_threshold,
    hurdle_combine_method,
    cpm_scale,
    pseudocount,
    n_permutations,
    random_state,
    negative_binomial_alpha,
):
    parameters = [
        ("log2fc_zero_value", log2fc_zero_value, False),
        ("presence_threshold", presence_threshold, True),
    ]
    if method in {"hurdle", "permutation"}:
        parameters.extend(
            [
                ("cpm_scale", cpm_scale, False),
                ("pseudocount", pseudocount, False),
            ]
        )
    if method == "negative_binomial":
        parameters.append(
            ("negative_binomial_alpha", negative_binomial_alpha, False)
        )
    for name, value, allow_zero in parameters:
        valid = isinstance(value, Real) and not isinstance(value, bool)
        valid = valid and (value >= 0 if allow_zero else value > 0)
        if not valid:
            comparison = "non-negative" if allow_zero else "positive"
            raise ValueError(f"{name} must be a {comparison} number")
    if method == "hurdle" and hurdle_combine_method not in {"fisher", "max"}:
        raise ValueError("hurdle_combine_method must be 'fisher' or 'max'")
    if method == "permutation" and (
        not isinstance(n_permutations, Integral)
        or isinstance(n_permutations, bool)
        or n_permutations < 1
    ):
        raise ValueError("n_permutations must be a positive integer")
    if method == "permutation" and random_state is not None and (
        not isinstance(random_state, Integral) or isinstance(random_state, bool)
    ):
        raise ValueError("random_state must be an integer or None")


def _resolve_sample_totals(count_table, sample_columns, sample_totals):
    if sample_totals is None:
        return count_table[sample_columns].sum(axis=0).to_numpy(dtype=float)
    try:
        totals = pd.Series(sample_totals, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("sample_totals must be a mapping or pandas Series") from error
    missing = [sample for sample in sample_columns if sample not in totals.index]
    if missing:
        raise ValueError(f"sample_totals is missing samples: {missing}")
    values = totals.loc[sample_columns].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("sample_totals must contain only finite values")
    return values


def _build_statistics_tasks(
    tested_values,
    original_tested_values,
    passed_positions,
    groups,
    group_sample_indices,
    method,
    library_totals,
    log2fc_zero_value,
    presence_threshold,
    hurdle_combine_method,
    cpm_scale,
    pseudocount,
    n_permutations,
    random_state,
    negative_binomial_alpha,
):
    comparisons = []
    if len(groups) == 2:
        comparisons.append((groups[0], groups[1], True))
    else:
        comparisons.extend((group, None, False) for group in groups)

    all_indices = np.arange(tested_values.shape[1])
    tasks = []
    for comparison_index, (target_group, background_group, two_sided) in enumerate(comparisons):
        target_indices = group_sample_indices[target_group]
        if background_group is None:
            background_indices = np.setdiff1d(all_indices, target_indices)
        else:
            background_indices = group_sample_indices[background_group]
        seed = None
        if method == "permutation" and random_state is not None:
            seed = int(random_state) + comparison_index
        tasks.append(
            {
                "values": tested_values,
                "original_values": original_tested_values,
                "row_positions": passed_positions,
                "target_group": target_group,
                "background_group": background_group,
                "target_indices": target_indices,
                "background_indices": background_indices,
                "two_sided": two_sided,
                "method": method,
                "library_totals": library_totals,
                "log2fc_zero_value": log2fc_zero_value,
                "presence_threshold": presence_threshold,
                "hurdle_combine_method": hurdle_combine_method,
                "cpm_scale": cpm_scale,
                "pseudocount": pseudocount,
                "n_permutations": n_permutations,
                "random_state": seed,
                "negative_binomial_alpha": negative_binomial_alpha,
            }
        )
    return tasks


def _calculate_statistics_task(task):
    values = task["values"]
    original_values = task["original_values"]
    target_indices = task["target_indices"]
    background_indices = task["background_indices"]
    target_values = values[:, target_indices]
    background_values = values[:, background_indices]
    original_target_values = original_values[:, target_indices]
    original_background_values = original_values[:, background_indices]
    target_means = target_values.mean(axis=1)
    background_means = background_values.mean(axis=1)
    log2fc = np.array(
        [
            _log2fc_for_two_means(target_mean, background_mean, task["log2fc_zero_value"])
            for target_mean, background_mean in zip(target_means, background_means)
        ]
    )
    p_values = _calculate_method_p_values(
        task,
        target_values=target_values,
        background_values=background_values,
    )
    enriched_in = np.full(len(values), task["target_group"], dtype=object)
    mean_group_count = original_target_values.mean(axis=1)
    if task["two_sided"]:
        background_enriched = background_means > target_means
        enriched_in[background_enriched] = task["background_group"]
        original_background_means = original_background_values.mean(axis=1)
        mean_group_count[background_enriched] = original_background_means[
            background_enriched
        ]
        log2fc = np.abs(log2fc)
    return pd.DataFrame(
        {
            "_row_position": task["row_positions"],
            "enriched_in": enriched_in,
            "method": task["method"],
            "mean_group_count": mean_group_count,
            "log2FC": log2fc,
            "p_val": p_values,
        }
    )


def _calculate_method_p_values(task, target_values, background_values):
    alternative = "two-sided" if task["two_sided"] else "greater"
    method = task["method"]
    if method == "mann_whitney":
        return np.array(
            [
                _mann_whitney_p_value(target, background, alternative)
                for target, background in zip(target_values, background_values)
            ]
        )
    if method == "fisher":
        return _fisher_presence_p_values(
            target_values,
            background_values,
            threshold=task["presence_threshold"],
            alternative=alternative,
        )
    if method == "fisher_count":
        return _fisher_count_p_values(task, target_values, background_values, alternative)
    if method == "hurdle":
        return _hurdle_p_values(task, target_values, background_values, alternative)
    if method == "quasi_binomial":
        return _glm_p_values(task, family="quasi_binomial")
    if method == "negative_binomial":
        return _glm_p_values(task, family="negative_binomial")
    if method == "permutation":
        return _permutation_p_values(task)
    raise RuntimeError(f"Unsupported method reached calculation: {method}")


def _mann_whitney_p_value(target, background, alternative):
    p_value = scipy.stats.mannwhitneyu(
        target,
        background,
        alternative=alternative,
    ).pvalue
    return 1.0 if np.isnan(p_value) else float(p_value)


def _fisher_presence_p_values(target_values, background_values, threshold, alternative):
    p_values = []
    for target, background in zip(target_values, background_values):
        target_present = int(np.sum(target >= threshold))
        background_present = int(np.sum(background >= threshold))
        table = [
            [target_present, len(target) - target_present],
            [background_present, len(background) - background_present],
        ]
        p_values.append(scipy.stats.fisher_exact(table, alternative=alternative).pvalue)
    return np.asarray(p_values)


def _fisher_count_p_values(task, target_values, background_values, alternative):
    target_totals = task["library_totals"][task["target_indices"]].sum()
    background_totals = task["library_totals"][task["background_indices"]].sum()
    p_values = []
    for target, background in zip(target_values, background_values):
        target_count = int(round(target.sum()))
        background_count = int(round(background.sum()))
        table = [
            [target_count, int(round(target_totals - target_count))],
            [background_count, int(round(background_totals - background_count))],
        ]
        p_values.append(scipy.stats.fisher_exact(table, alternative=alternative).pvalue)
    return np.asarray(p_values)


def _hurdle_p_values(task, target_values, background_values, alternative):
    presence_p = _fisher_presence_p_values(
        target_values,
        background_values,
        threshold=task["presence_threshold"],
        alternative=alternative,
    )
    totals = task["library_totals"]
    target_totals = totals[task["target_indices"]]
    background_totals = totals[task["background_indices"]]
    abundance_p = []
    for target, background in zip(target_values, background_values):
        target_present = target >= task["presence_threshold"]
        background_present = background >= task["presence_threshold"]
        target_log_cpm = np.log2(
            target[target_present] / target_totals[target_present] * task["cpm_scale"]
            + task["pseudocount"]
        )
        background_log_cpm = np.log2(
            background[background_present]
            / background_totals[background_present]
            * task["cpm_scale"]
            + task["pseudocount"]
        )
        if len(target_log_cpm) == 0 or len(background_log_cpm) == 0:
            abundance_p.append(1.0)
        else:
            abundance_p.append(
                scipy.stats.mannwhitneyu(
                    target_log_cpm,
                    background_log_cpm,
                    alternative=alternative,
                ).pvalue
            )
    abundance_p = np.asarray(abundance_p)
    if task["hurdle_combine_method"] == "max":
        return np.maximum(presence_p, abundance_p)
    return np.array(
        [
            scipy.stats.combine_pvalues([presence, abundance], method="fisher").pvalue
            for presence, abundance in zip(presence_p, abundance_p)
        ]
    )


def _glm_p_values(task, family):
    target_mask = np.zeros(task["values"].shape[1], dtype=float)
    target_mask[task["target_indices"]] = 1.0
    design = sm.add_constant(target_mask, has_constant="add")
    totals = task["library_totals"]
    p_values = []
    for counts in task["values"]:
        try:
            if family == "quasi_binomial":
                model = sm.GLM(
                    counts / totals,
                    design,
                    family=sm.families.Binomial(),
                    var_weights=totals,
                )
                fit = model.fit(scale="X2")
            else:
                model = sm.GLM(
                    counts,
                    design,
                    family=sm.families.NegativeBinomial(
                        alpha=task["negative_binomial_alpha"]
                    ),
                    offset=np.log(totals),
                )
                fit = model.fit()
            coefficient = float(fit.params[1])
            two_sided_p = float(fit.pvalues[1])
            if task["two_sided"]:
                p_value = two_sided_p
            elif coefficient > 0:
                p_value = two_sided_p / 2
            else:
                p_value = 1 - two_sided_p / 2
        except (ValueError, np.linalg.LinAlgError, PerfectSeparationError):
            p_value = 1.0
        p_values.append(min(max(p_value, 0.0), 1.0))
    return np.asarray(p_values)


def _permutation_p_values(task):
    values = task["values"]
    totals = task["library_totals"]
    log_cpm = np.log2(
        values / totals[np.newaxis, :] * task["cpm_scale"] + task["pseudocount"]
    )
    target_indices = task["target_indices"]
    background_indices = task["background_indices"]
    observed = (
        log_cpm[:, target_indices].mean(axis=1)
        - log_cpm[:, background_indices].mean(axis=1)
    )
    exceedances = np.zeros(len(values), dtype=int)
    rng = np.random.default_rng(task["random_state"])
    sample_count = values.shape[1]
    target_count = len(target_indices)
    for _ in range(task["n_permutations"]):
        permuted_target = rng.choice(sample_count, size=target_count, replace=False)
        permuted_background = np.setdiff1d(np.arange(sample_count), permuted_target)
        score = (
            log_cpm[:, permuted_target].mean(axis=1)
            - log_cpm[:, permuted_background].mean(axis=1)
        )
        if task["two_sided"]:
            exceedances += np.abs(score) >= np.abs(observed)
        else:
            exceedances += score >= observed
    return (exceedances + 1) / (task["n_permutations"] + 1)


def _log2fc_for_two_means(mean1, mean2, zero_value):
    if mean1 == 0 and mean2 == 0:
        return 0.0
    if mean1 == 0:
        return -float(zero_value)
    if mean2 == 0:
        return float(zero_value)
    return float(np.log2(mean1 / mean2))


def _adjust_p_values(p_values, method):
    adjusted = np.full(len(p_values), np.nan)
    finite = np.isfinite(p_values.to_numpy(dtype=float))
    if finite.any():
        try:
            adjusted[finite] = multipletests(
                p_values.to_numpy(dtype=float)[finite],
                method=method,
            )[1]
        except ValueError as error:
            raise ValueError(
                f"Invalid p_adjust_method {method!r}: {error}"
            ) from error
    return adjusted


def _validate_p_adjust_method(method):
    try:
        multipletests([1.0], method=method)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Invalid p_adjust_method {method!r}: {error}"
        ) from error


def _assemble_statistics_output(
    count_table,
    statistics,
    pass_mask,
    simplify,
):
    base_columns = [
        column
        for column in count_table.columns
        if column not in STATISTICS_COLUMNS
    ]
    base_table = count_table[base_columns].copy()
    row_position_column = "__repseq_row_position__"
    while row_position_column in base_table.columns:
        row_position_column = "_" + row_position_column
    comparison_order_column = "__repseq_comparison_order__"
    while comparison_order_column in base_table.columns:
        comparison_order_column = "_" + comparison_order_column

    base_table[row_position_column] = np.arange(len(base_table), dtype=int)
    pass_mask = np.asarray(pass_mask, dtype=bool)
    passing_table = base_table.loc[pass_mask]
    failing_table = base_table.loc[~pass_mask].copy()

    statistics_for_merge = statistics[
        ["_row_position", *STATISTICS_COLUMNS]
    ].copy()
    statistics_for_merge[comparison_order_column] = np.arange(
        len(statistics_for_merge),
        dtype=int,
    )
    statistics_for_merge = statistics_for_merge.rename(
        columns={"_row_position": row_position_column}
    )
    passing_result = passing_table.merge(
        statistics_for_merge,
        how="left",
        on=row_position_column,
        sort=False,
        validate="one_to_one" if simplify else "one_to_many",
    )

    for column in STATISTICS_COLUMNS:
        failing_table[column] = None
    failing_table[comparison_order_column] = -1

    result = pd.concat(
        [passing_result, failing_table],
        axis=0,
        ignore_index=True,
        sort=False,
    )
    result = result.sort_values(
        [row_position_column, comparison_order_column],
        kind="stable",
    )
    row_positions = result[row_position_column].to_numpy(dtype=int)
    result = result.drop(
        columns=[row_position_column, comparison_order_column]
    )
    result.index = count_table.index.take(row_positions)

    original_numeric_columns = [
        column
        for column in count_table.select_dtypes(include="number").columns
        if column in base_columns
    ]
    insert_position = (
        base_columns.index(original_numeric_columns[0])
        if original_numeric_columns
        else len(base_columns)
    )
    ordered_columns = [
        *base_columns[:insert_position],
        *STATISTICS_COLUMNS,
        *base_columns[insert_position:],
    ]
    result = result[ordered_columns]
    result.attrs = count_table.attrs.copy()
    return result
