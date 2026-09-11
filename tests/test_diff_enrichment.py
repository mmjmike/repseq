import numpy as np
import pandas as pd
import pytest

from repseq import diff_enrichment as rsde


def _count_table():
    return pd.DataFrame(
        {
            "clonotype": ["AAA|V1", "BBB|V2", "CCC|V3", "DDD|V4"],
            "cdr3aa": ["AAA", "BBB", "CCC", "DDD"],
            "v": ["V1", "V2", "V3", "V4"],
            "sample1": [2, 4, 2, 1],
            "sample2": [2, 4, 2, 1],
            "sample3": [4, 0, 3, 6],
        }
    )


def test_prefilter_marks_rows_and_inserts_status_before_numeric_columns():
    count_table = _count_table()

    result = rsde.prefilter(count_table, verbose=False)

    assert list(result.columns) == [
        "clonotype",
        "cdr3aa",
        "v",
        "prefilter_pass",
        "sample1",
        "sample2",
        "sample3",
    ]
    assert result["prefilter_pass"].tolist() == [False, False, False, False]
    pd.testing.assert_frame_equal(
        result.drop(columns="prefilter_pass"),
        count_table,
    )
    assert result is not count_table


def test_prefilter_uses_inclusive_thresholds_and_ignores_non_numeric_columns():
    count_table = pd.DataFrame(
        {
            "feature": ["passes exactly", "fails sample threshold"],
            "annotation": ["100", "100"],
            "sample1": [2, 7],
            "sample2": [2, 1],
            "sample3": [6, 0],
        }
    )

    result = rsde.prefilter(count_table, verbose=False)

    assert result["prefilter_pass"].tolist() == [True, False]


def test_prefilter_prints_all_settings_and_passed_feature_count(capsys):
    rsde.prefilter(
        _count_table(),
        min_samples=2,
        min_count=3,
        min_total_count=7,
    )

    output = capsys.readouterr().out
    assert "Minimum samples: 2" in output
    assert "Minimum count per sample: 3" in output
    assert "Minimum total count: 7" in output
    assert "Features passed: 1 of 4" in output


def test_prefilter_places_status_at_end_when_there_are_no_numeric_columns():
    count_table = pd.DataFrame({"feature": ["a", "b"]})

    result = rsde.prefilter(
        count_table,
        min_samples=0,
        min_total_count=0,
        verbose=False,
    )

    assert list(result.columns) == ["feature", "prefilter_pass"]
    assert result["prefilter_pass"].tolist() == [True, True]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"min_samples": 1.5}, "min_samples"),
        ({"min_samples": -1}, "min_samples"),
        ({"min_count": -1}, "min_count"),
        ({"min_total_count": -1}, "min_total_count"),
    ],
)
def test_prefilter_rejects_invalid_thresholds(kwargs, message):
    with pytest.raises(ValueError, match=message):
        rsde.prefilter(_count_table(), verbose=False, **kwargs)


def test_prefilter_rejects_existing_status_column():
    count_table = _count_table().assign(prefilter_pass=True)

    with pytest.raises(ValueError, match="already contains"):
        rsde.prefilter(count_table, verbose=False)


def _postfilter_table():
    table = pd.DataFrame(
        {
            "feature": ["exact", "low_mean", "low_logfc", "high_p", "missing"],
            "prefilter_pass": [True, True, True, True, False],
            "enriched_in": ["A", "A", "B", "B", "A"],
            "method": ["mann_whitney"] * 5,
            "mean_group_count": [2, 1.99, 3, 3, 10],
            "log2FC": [1, 2, 0.99, 2, 10],
            "p_val": [1, 0.01, 0.01, 0.999, 0.001],
            "p_adj": [1, 0.01, 0.01, 0.999, 0.001],
            "sample1": [2, 3, 4, 5, 0],
        }
    )
    table.attrs["source"] = "statistics"
    return table


def test_postfilter_uses_strict_p_values_and_inclusive_minimums():
    statistics_table = _postfilter_table()

    result = rsde.postfilter(
        statistics_table,
        max_p_adj=1,
        max_p_val=1,
        min_group_mean=2,
        min_logfc=1,
        sort=False,
        verbose=False,
    )

    assert result["postfilter_pass"].tolist() == [False, False, False, True, False]
    assert list(result.columns) == [
        "feature",
        "prefilter_pass",
        "enriched_in",
        "method",
        "mean_group_count",
        "log2FC",
        "p_val",
        "p_adj",
        "postfilter_pass",
        "sample1",
    ]
    pd.testing.assert_frame_equal(
        result.drop(columns="postfilter_pass"),
        statistics_table,
    )
    assert result.attrs == {"source": "statistics"}


@pytest.mark.parametrize(
    ("groups", "expected"),
    [
        (None, [True, False, False, True, False]),
        ("A", [True, False, False, False, False]),
        (["B"], [False, False, False, True, False]),
        (["A", "B"], [True, False, False, True, False]),
    ],
)
def test_postfilter_optionally_filters_enriched_groups(groups, expected):
    result = rsde.postfilter(
        _postfilter_table(), groups=groups, sort=False, verbose=False
    )

    assert result["postfilter_pass"].tolist() == expected


@pytest.mark.parametrize(
    ("kwargs", "error_type", "message"),
    [
        ({"max_p_adj": 1.1}, ValueError, "max_p_adj"),
        ({"max_p_val": -0.1}, ValueError, "max_p_val"),
        ({"min_group_mean": -1}, ValueError, "min_group_mean"),
        ({"min_logfc": "one"}, ValueError, "min_logfc"),
        ({"groups": []}, ValueError, "groups must not be empty"),
        ({"groups": ["A", 2]}, TypeError, "only strings"),
        ({"groups_exclude": []}, ValueError, "groups_exclude must not be empty"),
        ({"groups_exclude": ["A", 2]}, TypeError, "only strings"),
    ],
)
def test_postfilter_rejects_invalid_arguments(kwargs, error_type, message):
    with pytest.raises(error_type, match=message):
        rsde.postfilter(_postfilter_table(), verbose=False, **kwargs)


@pytest.mark.parametrize(
    ("groups_exclude", "expected"),
    [
        (None, [True, False, False, True, False]),
        ("A", [False, False, False, True, False]),
        (["B"], [True, False, False, False, False]),
        (["A", "B"], [False, False, False, False, False]),
    ],
)
def test_postfilter_optionally_excludes_enriched_groups(groups_exclude, expected):
    result = rsde.postfilter(
        _postfilter_table(),
        groups_exclude=groups_exclude,
        sort=False,
        verbose=False,
    )

    assert result["postfilter_pass"].tolist() == expected


def test_postfilter_groups_override_groups_exclude(capsys):
    result = rsde.postfilter(
        _postfilter_table(),
        groups="A",
        groups_exclude=["A", "not_present"],
        sort=False,
    )

    assert result["postfilter_pass"].tolist() == [True, False, False, False, False]
    output = capsys.readouterr().out
    assert "groups_exclude: skipped because groups overrides groups_exclude" in output
    assert "not_present" not in output


def test_postfilter_silently_ignores_unmatched_groups_exclude(capsys):
    result = rsde.postfilter(
        _postfilter_table(),
        groups_exclude=["not_present"],
        sort=False,
    )

    assert result["postfilter_pass"].tolist() == [True, False, False, True, False]
    output = capsys.readouterr().out
    assert "groups_exclude" not in output
    assert "not_present" not in output


def test_postfilter_rejects_existing_status_column():
    statistics_table = _postfilter_table().assign(postfilter_pass=True)

    with pytest.raises(ValueError, match="already contains"):
        rsde.postfilter(statistics_table, verbose=False)


def test_postfilter_verbose_reports_comparisons_and_prefilter_denominator(capsys):
    rsde.postfilter(
        _postfilter_table(),
        max_p_adj=0.05,
        max_p_val=None,
        min_group_mean=2,
        min_logfc=1,
        groups="A",
    )

    output = capsys.readouterr().out
    assert "Adjusted p-value: p_adj < 0.05" in output
    assert "Raw p-value: any" in output
    assert "Group mean count: mean_group_count >= 2" in output
    assert "Log2 fold change: log2FC >= 1" in output
    assert "Enriched-in groups: enriched_in in ['A']" in output
    assert "Features passed prefilter: 4 of 5" in output
    assert "Features passed postfilter: 0 of 4 prefilter-passing features" in output


def test_postfilter_verbose_reports_full_denominator_without_prefilter(capsys):
    statistics_table = _postfilter_table().drop(columns="prefilter_pass")

    result = rsde.postfilter(statistics_table, sort=False)

    assert result["postfilter_pass"].tolist() == [True, False, False, True, True]
    output = capsys.readouterr().out
    assert "Adjusted p-value: any" in output
    assert "Raw p-value: any" in output
    assert "Features passed postfilter: 3 of 5" in output


def test_postfilter_sorts_by_default_priorities():
    result = rsde.postfilter(_postfilter_table(), verbose=False)

    assert result["feature"].tolist() == [
        "exact",
        "high_p",
        "low_mean",
        "low_logfc",
        "missing",
    ]


def test_postfilter_sorts_custom_and_categorical_columns():
    statistics_table = _postfilter_table()
    statistics_table.loc[statistics_table["feature"].eq("low_logfc"), "log2FC"] = 2
    statistics_table["enriched_in"] = pd.Categorical(
        statistics_table["enriched_in"],
        categories=["B", "A"],
        ordered=True,
    )

    categorical = rsde.postfilter(
        statistics_table,
        sort=["enriched_in"],
        verbose=False,
    )
    arbitrary = rsde.postfilter(
        statistics_table,
        sort=["sample1"],
        verbose=False,
    )
    directional = rsde.postfilter(
        statistics_table,
        sort=["mean_group_count", "log2FC", "p_adj"],
        verbose=False,
    )

    assert categorical["feature"].tolist() == [
        "low_logfc",
        "high_p",
        "exact",
        "low_mean",
        "missing",
    ]
    assert arbitrary["feature"].tolist() == [
        "missing",
        "exact",
        "low_mean",
        "low_logfc",
        "high_p",
    ]
    assert directional["feature"].tolist() == [
        "missing",
        "low_logfc",
        "high_p",
        "exact",
        "low_mean",
    ]


@pytest.mark.parametrize("sort", [False, None, [], ["not_a_column"]])
def test_postfilter_skips_sorting_without_usable_columns(sort):
    result = rsde.postfilter(_postfilter_table(), sort=sort, verbose=False)

    assert result["feature"].tolist() == _postfilter_table()["feature"].tolist()
    assert result.attrs == {"source": "statistics"}


def _two_group_metadata():
    return pd.DataFrame(
        {
            "sample_id": ["a1", "a2", "b1", "b2"],
            "group": ["A", "A", "B", "B"],
        }
    )


def _statistics_count_table(prefilter=True):
    data = {
        "feature": ["f1", "f2", "f3"],
        "annotation": ["first", "second", "filtered"],
        "a1": [8, 0, 100],
        "a2": [8, 0, 100],
        "b1": [2, 4, 0],
        "b2": [2, 4, 0],
        "ignored_sample": [20, 20, 20],
    }
    if prefilter:
        data = {
            "feature": data.pop("feature"),
            "annotation": data.pop("annotation"),
            "prefilter_pass": [True, True, False],
            **data,
        }
    return pd.DataFrame(data)


def test_calc_statistics_two_groups_preserves_rows_and_prefilter_status(capsys):
    count_table = _statistics_count_table()

    result = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        cpu=1,
    )

    assert list(result.columns) == [
        "feature",
        "annotation",
        "prefilter_pass",
        "enriched_in",
        "method",
        "mean_group_count",
        "log2FC",
        "p_val",
        "p_adj",
        "a1",
        "a2",
        "b1",
        "b2",
        "ignored_sample",
    ]
    assert result["feature"].tolist() == ["f1", "f2", "f3"]
    assert result.loc[0, "enriched_in"] == "A"
    assert result.loc[0, "mean_group_count"] == 8
    assert result.loc[0, "log2FC"] == pytest.approx(2)
    assert result.loc[1, "enriched_in"] == "B"
    assert result.loc[1, "mean_group_count"] == 4
    assert result.loc[1, "log2FC"] == 100
    assert result.loc[:1, "method"].tolist() == ["mann_whitney"] * 2
    assert result.loc[2, rsde.STATISTICS_COLUMNS].isna().all()
    pd.testing.assert_frame_equal(
        result.drop(columns=rsde.STATISTICS_COLUMNS),
        count_table,
    )

    output = capsys.readouterr().out
    assert "Differential enrichment analysis started" in output
    assert "Groups detected: 2" in output
    assert "Group 'A' (2 samples): ['a1', 'a2']" in output
    assert "Group 'B' (2 samples): ['b1', 'b2']" in output
    assert "ignored_sample" in output
    assert "2 of 3 features will be tested" in output
    assert "Analysis count threshold: 1" in output
    assert "Original counts will be preserved" in output
    assert "Combining group-comparison result tables" in output
    assert "Adjusting p-values for multiple testing using 'fdr_bh'" in output
    assert "Simplifying statistics" in output
    assert "Assembling the final differential enrichment output table" in output
    assert "Differential enrichment analysis finished successfully" in output


def test_calc_statistics_reports_when_simplification_is_disabled(capsys):
    rsde.calc_statistics(
        _statistics_count_table(),
        _two_group_metadata(),
        simplify=False,
        cpu=1,
    )

    output = capsys.readouterr().out
    assert "Statistics simplification is disabled" in output
    assert "Assembling the final differential enrichment output table" in output


def test_calc_statistics_without_prefilter_analyzes_all_features():
    result = rsde.calc_statistics(
        _statistics_count_table(prefilter=False),
        _two_group_metadata(),
        cpu=1,
        verbose=False,
    )

    assert result["method"].notna().all()


def test_calc_statistics_applies_default_and_custom_sorting():
    count_table = _statistics_count_table().iloc[[1, 0, 2]]

    default_sorted = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        cpu=1,
        verbose=False,
    )
    custom_sorted = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        sort=["feature"],
        cpu=1,
        verbose=False,
    )
    unsorted = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        sort=["not_a_column"],
        cpu=1,
        verbose=False,
    )

    assert default_sorted["feature"].tolist() == ["f1", "f2", "f3"]
    assert custom_sorted["feature"].tolist() == ["f1", "f2", "f3"]
    assert unsorted["feature"].tolist() == ["f2", "f1", "f3"]


@pytest.mark.parametrize("ordered", [False, True])
def test_calc_statistics_preserves_categorical_group_dtype(ordered):
    samples_metadata = _two_group_metadata()
    samples_metadata["group"] = pd.Categorical(
        samples_metadata["group"],
        categories=["B", "A", "unused"],
        ordered=ordered,
    )

    result = rsde.calc_statistics(
        _statistics_count_table(),
        samples_metadata,
        cpu=1,
        verbose=False,
    )

    assert isinstance(result["enriched_in"].dtype, pd.CategoricalDtype)
    assert result["enriched_in"].dtype.categories.tolist() == ["B", "A", "unused"]
    assert result["enriched_in"].dtype.ordered is ordered
    assert result["feature"].tolist() == ["f2", "f1", "f3"]


def test_presence_threshold_changes_analysis_but_preserves_real_counts():
    count_table = pd.DataFrame(
        {
            "feature": ["low_count_feature"],
            "a1": [1.0],
            "a2": [1.0],
            "b1": [0.0],
            "b2": [0.0],
        }
    )

    unthresholded = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        method="mann_whitney",
        presence_threshold=0,
        cpu=1,
        verbose=False,
    )
    thresholded = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        method="mann_whitney",
        presence_threshold=2,
        cpu=1,
        verbose=False,
    )

    assert unthresholded.loc[0, "log2FC"] == 100
    assert thresholded.loc[0, "log2FC"] == 0
    assert thresholded.loc[0, "p_val"] == 1
    assert thresholded.loc[0, "mean_group_count"] == 1
    pd.testing.assert_frame_equal(
        thresholded[["feature", "a1", "a2", "b1", "b2"]],
        count_table,
    )


def test_calc_statistics_preserves_dataframe_attributes():
    count_table = _statistics_count_table()
    count_table.attrs["source"] = "count_table"

    result = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        cpu=1,
        verbose=False,
    )

    assert result.attrs == {"source": "count_table"}


def _three_group_inputs():
    count_table = pd.DataFrame(
        {
            "feature": ["enriched_a", "enriched_b", "filtered"],
            "prefilter_pass": [True, True, False],
            "a1": [10, 0, 3],
            "a2": [10, 0, 3],
            "b1": [0, 10, 3],
            "b2": [0, 10, 3],
            "c1": [0, 0, 3],
            "c2": [0, 0, 3],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["a1", "a2", "b1", "b2", "c1", "c2"],
            "group": ["A", "A", "B", "B", "C", "C"],
        }
    )
    return count_table, metadata


def test_calc_statistics_multigroup_simplifies_to_best_group():
    count_table, metadata = _three_group_inputs()

    result = rsde.calc_statistics(
        count_table,
        metadata,
        method="fisher",
        presence_threshold=2,
        cpu=1,
        verbose=False,
    )

    assert result["feature"].tolist() == ["enriched_a", "enriched_b", "filtered"]
    assert result.loc[:1, "enriched_in"].tolist() == ["A", "B"]
    assert result.loc[:1, "mean_group_count"].tolist() == [10, 10]
    assert result.loc[:1, "log2FC"].tolist() == [100, 100]
    assert result.loc[2, rsde.STATISTICS_COLUMNS].isna().all()


def test_calc_statistics_multigroup_expands_passed_features_and_parallelizes(monkeypatch):
    count_table, metadata = _three_group_inputs()
    seen_tasks = []

    def run_sequential(function, tasks, *args, **kwargs):
        seen_tasks.extend(tasks)
        return [function(task) for task in tasks]

    monkeypatch.setattr(rsde, "run_parallel_calculation", run_sequential)

    result = rsde.calc_statistics(
        count_table,
        metadata,
        method="fisher",
        simplify=False,
        verbose=False,
    )

    assert len(seen_tasks) == 3
    assert len(result) == 7
    assert result.loc[result["feature"].eq("enriched_a"), "enriched_in"].tolist() == [
        "A",
        "B",
        "C",
    ]
    assert result.loc[
        result["feature"].eq("enriched_a"), "mean_group_count"
    ].tolist() == [10, 0, 0]
    assert result.loc[result["feature"].eq("filtered"), "method"].isna().all()


def test_calc_statistics_expanded_output_preserves_original_index_order():
    count_table, metadata = _three_group_inputs()
    count_table.index = pd.Index([10, 10, 30], name="source_row")
    count_table.attrs["source"] = "custom-index-table"

    result = rsde.calc_statistics(
        count_table,
        metadata,
        method="fisher",
        simplify=False,
        cpu=1,
        sort=False,
        verbose=False,
    )

    assert result.index.tolist() == [10, 10, 10, 10, 10, 10, 30]
    assert result.index.name == "source_row"
    assert result["feature"].tolist() == [
        "enriched_a",
        "enriched_a",
        "enriched_a",
        "enriched_b",
        "enriched_b",
        "enriched_b",
        "filtered",
    ]
    assert result.attrs == {"source": "custom-index-table"}


def test_calc_statistics_reports_default_feature_column_duplicates():
    count_table = _statistics_count_table().copy()
    count_table.index = [10, 11, 12]
    count_table["feature"] = ["duplicate", "duplicate", "unique"]

    with pytest.raises(
        ValueError,
        match=r"'feature'.*non-unique.*'duplicate'.*\[10, 11\]",
    ):
        rsde.calc_statistics(
            count_table,
            _two_group_metadata(),
            verbose=False,
        )


def test_calc_statistics_reports_duplicate_missing_feature_indices():
    count_table = _statistics_count_table().copy()
    count_table.index = [4, 5, 6]
    count_table["feature"] = [None, None, "unique"]

    with pytest.raises(ValueError, match=r"non-unique.*\[4, 5\]"):
        rsde.calc_statistics(
            count_table,
            _two_group_metadata(),
            verbose=False,
        )


def test_calc_statistics_accepts_explicit_unique_feature_column():
    count_table = _statistics_count_table().copy()
    count_table.insert(1, "unique_id", [1, 2, 3])
    count_table["feature"] = "same annotation"

    result = rsde.calc_statistics(
        count_table,
        _two_group_metadata(),
        feature_column="unique_id",
        cpu=1,
        verbose=False,
    )

    assert result["unique_id"].tolist() == [1, 2, 3]


@pytest.mark.parametrize(
    ("metadata", "message"),
    [
        (pd.DataFrame({"sample_id": ["a1"]}), "must contain columns"),
        (
            pd.DataFrame(
                {
                    "sample_id": ["a1", "a2", "b1", "b1"],
                    "group": ["A", "A", "B", "B"],
                }
            ),
            "must be unique",
        ),
        (
            pd.DataFrame(
                {
                    "sample_id": ["a1", "a2", "b1", "b2", "annotation"],
                    "group": ["A", "A", "B", "B", "annotation_group"],
                }
            ),
            "non-numeric samples",
        ),
        (
            pd.DataFrame(
                {
                    "sample_id": ["a1", "a2", "b1", "b2"],
                    "group": ["A", "A", "A", "A"],
                }
            ),
            "at least 2 groups",
        ),
        (
            pd.DataFrame(
                {
                    "sample_id": ["a1", "a2", "b1"],
                    "group": ["A", "A", "B"],
                }
            ),
            "Each group must contain at least 2 samples",
        ),
    ],
)
def test_calc_statistics_validates_metadata(metadata, message):
    with pytest.raises(ValueError, match=message):
        rsde.calc_statistics(
            _statistics_count_table(),
            metadata,
            verbose=False,
        )


def test_calc_statistics_ignores_metadata_samples_absent_from_count_table(capsys):
    metadata = pd.concat(
        [
            _two_group_metadata(),
            pd.DataFrame(
                {
                    "sample_id": ["extra1", "extra2"],
                    "group": ["ExtraGroup", "ExtraGroup"],
                }
            ),
        ],
        ignore_index=True,
    )

    result = rsde.calc_statistics(
        _statistics_count_table(),
        metadata,
        cpu=1,
    )

    assert result.loc[:1, "method"].tolist() == ["mann_whitney"] * 2
    output = capsys.readouterr().out
    assert "samples_metadata entries absent from count_table will be ignored" in output
    assert "['extra1', 'extra2']" in output
    assert "Groups detected: 2" in output


@pytest.mark.parametrize(
    "method",
    [
        "mann_whitney",
        "fisher",
        "fisher_count",
        "hurdle",
        "quasi_binomial",
        "negative_binomial",
        "permutation",
    ],
)
def test_calc_statistics_supported_methods_return_standard_columns(method):
    result = rsde.calc_statistics(
        _statistics_count_table(),
        _two_group_metadata(),
        method=method,
        n_permutations=20,
        random_state=7,
        cpu=1,
        verbose=False,
    )

    assert list(result.loc[:1, "method"]) == [method, method]
    assert result.loc[:1, "p_val"].between(0, 1).all()
    assert result.loc[:1, "p_adj"].between(0, 1).all()


@pytest.mark.parametrize("method", ["fisher", "fisher_count", "hurdle"])
def test_calc_statistics_accepts_tuple_scipy_results(monkeypatch, method):
    monkeypatch.setattr(
        rsde.scipy.stats,
        "fisher_exact",
        lambda *args, **kwargs: (1.0, 0.25),
    )
    monkeypatch.setattr(
        rsde.scipy.stats,
        "mannwhitneyu",
        lambda *args, **kwargs: (1.0, 0.5),
    )
    monkeypatch.setattr(
        rsde.scipy.stats,
        "combine_pvalues",
        lambda *args, **kwargs: (1.0, 0.4),
    )

    result = rsde.calc_statistics(
        _statistics_count_table(),
        _two_group_metadata(),
        method=method,
        cpu=1,
        verbose=False,
    )

    assert result.loc[:1, "p_val"].between(0, 1).all()


def test_calc_statistics_ignores_parameters_for_other_methods():
    result = rsde.calc_statistics(
        _statistics_count_table(),
        _two_group_metadata(),
        method="mann_whitney",
        sample_totals={"not": "usable"},
        hurdle_combine_method="not-a-method",
        cpm_scale=-1,
        pseudocount=-1,
        n_permutations=0,
        random_state="not-an-integer",
        negative_binomial_alpha=-1,
        cpu=1,
        verbose=False,
    )

    assert result.loc[:1, "method"].tolist() == ["mann_whitney"] * 2


def test_calc_statistics_validates_presence_threshold_for_mann_whitney():
    with pytest.raises(ValueError, match="presence_threshold"):
        rsde.calc_statistics(
            _statistics_count_table(),
            _two_group_metadata(),
            method="mann_whitney",
            presence_threshold=-1,
            verbose=False,
        )


def test_calc_statistics_validates_p_adjust_method_when_no_features_pass():
    count_table = _statistics_count_table().assign(prefilter_pass=False)

    with pytest.raises(ValueError, match="Invalid p_adjust_method"):
        rsde.calc_statistics(
            count_table,
            _two_group_metadata(),
            p_adjust_method="not-a-method",
            verbose=False,
        )


def test_simplify_keeps_lowest_p_value_per_feature_position():
    statistics = pd.DataFrame(
        {
            "_row_position": [0, 0, 1, 1],
            "enriched_in": ["A", "B", "A", "B"],
            "method": ["fisher"] * 4,
            "mean_group_count": [5, 4, 8, 7],
            "log2FC": [2, 1, 3, 4],
            "p_val": [0.2, 0.1, 0.01, 0.02],
            "p_adj": [0.2, 0.2, 0.04, 0.04],
        }
    )

    result = rsde.simplify(statistics)

    assert result["enriched_in"].tolist() == ["B", "A"]


def _paired_chain_inputs():
    count_table1 = pd.DataFrame(
        {
            "feature1": ["a1", "a2", "a_filtered"],
            "prefilter_pass": [True, True, False],
            "enriched_in": ["group1", "group2", "group1"],
            "mean_group_count": [5.0, 4.0, 100.0],
            "tra_s1": [10, 0, 100],
            "tra_s2": [0, 5, 100],
            "tra_s3": [0, 5, 100],
        }
    )
    count_table2 = pd.DataFrame(
        {
            "feature2": ["b1", "b2", "b_filtered"],
            "prefilter_pass": [True, True, False],
            "enriched_in": ["group1", "group2", "group1"],
            "mean_group_count": [4.0, 3.0, 100.0],
            "trb_s1": [8, 0, 100],
            "trb_s2": [0, 4, 100],
            "trb_s3": [0, 4, 100],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": [
                "tra_s1",
                "trb_s1",
                "tra_s2",
                "trb_s2",
                "tra_s3",
                "trb_s3",
                "unused_sample",
            ],
            "sample": ["s1", "s1", "s2", "s2", "s3", "s3", "unused"],
        }
    )
    return count_table1, count_table2, metadata


def test_pair_chains_jsd_normalizes_rows_and_excludes_prefiltered_features():
    count_table1, count_table2, metadata = _paired_chain_inputs()

    result = rsde.pair_chains(count_table1, count_table2, metadata)

    assert result.index.tolist() == ["b1", "b2"]
    assert result.columns.tolist() == ["a1", "a2"]
    assert result.index.name == "feature2"
    assert result.columns.name == "feature1"
    assert result.loc["b1", "a1"] == pytest.approx(0)
    assert result.loc["b2", "a2"] == pytest.approx(0)
    assert result.loc["b1", "a2"] == pytest.approx(np.log(2))
    assert result.loc["b2", "a1"] == pytest.approx(np.log(2))
    assert result.attrs["method"] == "jsd"
    assert result.attrs["paired_samples"].to_dict("records") == [
        {"sample": "s1", "sample_id1": "tra_s1", "sample_id2": "trb_s1"},
        {"sample": "s2", "sample_id1": "tra_s2", "sample_id2": "trb_s2"},
        {"sample": "s3", "sample_id1": "tra_s3", "sample_id2": "trb_s3"},
    ]


def test_pair_chains_pearson_returns_correlation_matrix():
    count_table1, count_table2, metadata = _paired_chain_inputs()

    result = rsde.pair_chains(
        count_table1,
        count_table2,
        metadata,
        method="pearson",
    )

    assert result.loc["b1", "a1"] == pytest.approx(1)
    assert result.loc["b2", "a2"] == pytest.approx(1)
    assert result.loc["b1", "a2"] == pytest.approx(-1)
    assert result.loc["b2", "a1"] == pytest.approx(-1)


def test_pair_chains_filters_requested_ids_and_adds_best_opposite_partners():
    count_table1, count_table2, metadata = _paired_chain_inputs()

    result = rsde.pair_chains(
        count_table1,
        count_table2,
        metadata,
        filter_ids1=["a1"],
        filter_ids2=["b2"],
    )

    assert result.columns.tolist() == ["a1", "a2"]
    assert result.index.tolist() == ["b2", "b1"]


def test_pair_chains_single_filter_returns_requested_ids_and_best_partners():
    count_table1, count_table2, metadata = _paired_chain_inputs()

    result1 = rsde.pair_chains(
        count_table1,
        count_table2,
        metadata,
        filter_ids1=["a1"],
    )
    result2 = rsde.pair_chains(
        count_table1,
        count_table2,
        metadata,
        filter_ids2=["b2"],
    )

    assert result1.columns.tolist() == ["a1"]
    assert result1.index.tolist() == ["b1"]
    assert result2.columns.tolist() == ["a2"]
    assert result2.index.tolist() == ["b2"]


def test_pair_chains_rejects_unpaired_count_table_samples():
    count_table1, count_table2, metadata = _paired_chain_inputs()
    count_table2 = count_table2.drop(columns="trb_s3")

    with pytest.raises(ValueError, match="must have a paired sample_id"):
        rsde.pair_chains(count_table1, count_table2, metadata)


def test_pair_chains_rejects_filter_ids_removed_by_prefilter():
    count_table1, count_table2, metadata = _paired_chain_inputs()

    with pytest.raises(ValueError, match="absent after prefiltering"):
        rsde.pair_chains(
            count_table1,
            count_table2,
            metadata,
            filter_ids1=["a_filtered"],
        )
