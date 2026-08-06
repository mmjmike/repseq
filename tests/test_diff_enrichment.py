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
    assert result["prefilter_pass"].tolist() == [True, False, False, False]
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
            "sample3": [4, 0],
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
