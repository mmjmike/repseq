import pandas as pd
import pytest

from repseq import intersections


def _run_sequential(function, tasks, *args, **kwargs):
    return [function(task) for task in tasks]


def _sequence_clonoset_dicts():
    return {
        "s1": {
            (3, "V1"): [("AAA", 6), ("AAT", 4)],
        },
        "s2": {
            (3, "V1"): [("AAA", 5), ("AAC", 3), ("GGG", 2)],
        },
    }


def _mock_preparation(monkeypatch, clonoset_dicts, sample_list, sample_list2=None):
    def fake_prepare(*args, **kwargs):
        two_dataframes = sample_list2 is not None
        return (
            clonoset_dicts,
            len(sample_list) + (len(sample_list2) if sample_list2 else 0),
            two_dataframes,
            sample_list,
            sample_list2,
        )

    monkeypatch.setattr(
        intersections,
        "prepare_clonotypes_dfs_for_intersections",
        fake_prepare,
    )
    monkeypatch.setattr(
        intersections,
        "run_parallel_calculation",
        _run_sequential,
    )


def _dummy_samples(sample_ids):
    return pd.DataFrame(
        [
            {"sample_id": sample_id, "filename": f"{sample_id}.tsv"}
            for sample_id in sample_ids
        ]
    )


def test_similarity_frequency_is_directional_and_avoids_double_counting(monkeypatch):
    _mock_preparation(
        monkeypatch,
        _sequence_clonoset_dicts(),
        ["s1", "s2"],
    )

    matrix = intersections.similarity(
        _dummy_samples(["s1", "s2"]),
        overlap_type="aaV",
        mismatches=1,
        result="freq",
        cpu=3,
    )

    assert matrix.loc["s1", "s1"] == 1
    assert matrix.loc["s1", "s2"] == 1
    assert matrix.loc["s2", "s1"] == pytest.approx(0.8)
    assert matrix.loc["s2", "s2"] == 1
    assert matrix.index.name == "sample1"
    assert matrix.columns.name == "sample2"


def test_similarity_count_and_number_use_unique_target_clonotypes(monkeypatch):
    _mock_preparation(
        monkeypatch,
        _sequence_clonoset_dicts(),
        ["s1", "s2"],
    )
    samples = _dummy_samples(["s1", "s2"])

    counts = intersections.similarity(
        samples,
        overlap_type="aaV",
        mismatches=1,
        result="count",
    )
    numbers = intersections.similarity(
        samples,
        overlap_type="aaV",
        mismatches=1,
        result="number",
    )

    assert counts.loc["s1", "s2"] == 10
    assert counts.loc["s2", "s1"] == 8
    assert numbers.loc["s1", "s2"] == 2
    assert numbers.loc["s2", "s1"] == 2
    assert numbers.loc["s2", "s2"] == 3
    assert pd.api.types.is_integer_dtype(counts.dtypes.iloc[0])
    assert pd.api.types.is_integer_dtype(numbers.dtypes.iloc[0])


def test_similarity_table_keeps_all_matching_pairs_with_counts_and_frequencies(
    monkeypatch,
):
    _mock_preparation(
        monkeypatch,
        _sequence_clonoset_dicts(),
        ["s1", "s2"],
    )

    table = intersections.similarity(
        _dummy_samples(["s1", "s2"]),
        overlap_type="aaV",
        mismatches=1,
        result="table",
    )

    assert list(table.columns) == intersections._similarity_table_columns()
    assert len(table) == 8
    assert set(table["pair"]) == {"s1_vs_s2", "s2_vs_s1"}
    assert set(table["mismatches"]) == {0, 1}
    s1_aaa = table.loc[
        (table["sample1"] == "s1") & (table["clone1"] == ("AAA", "V1"))
    ]
    assert len(s1_aaa) == 2
    assert set(s1_aaa["sample1_count"]) == {6}
    assert all(value == pytest.approx(0.6) for value in s1_aaa["sample1_freq"])


def test_similarity_rectangular_matrices_keep_target_samples_in_rows(monkeypatch):
    clonoset_dicts = {
        "target": {(3,): [("AAA", 7), ("GGG", 3)]},
        "comparison": {(3,): [("AAA", 10)]},
    }
    _mock_preparation(
        monkeypatch,
        clonoset_dicts,
        ["target"],
        ["comparison"],
    )

    matrix = intersections.similarity(
        _dummy_samples(["target"]),
        clonosets_df2=_dummy_samples(["comparison"]),
        overlap_type="aa",
        mismatches=0,
        result="freq",
    )

    assert list(matrix.index) == ["target"]
    assert list(matrix.columns) == ["comparison"]
    assert matrix.loc["target", "comparison"] == pytest.approx(0.7)


def test_similarity_vj_ignores_mismatch_argument(monkeypatch):
    clonoset_dicts = {
        "s1": {("V1", "J1"): 10},
        "s2": {("V1", "J1"): 5, ("V2", "J2"): 5},
    }
    _mock_preparation(monkeypatch, clonoset_dicts, ["s1", "s2"])

    matrix = intersections.similarity(
        _dummy_samples(["s1", "s2"]),
        overlap_type="VJ",
        result="freq",
    )

    assert matrix.loc["s1", "s2"] == 1
    assert matrix.loc["s2", "s1"] == pytest.approx(0.5)
    assert matrix.attrs["mismatches"] == 0


def test_similarity_validates_arguments_and_deprecates_by_freq(monkeypatch):
    _mock_preparation(monkeypatch, _sequence_clonoset_dicts(), ["s1", "s2"])
    samples = _dummy_samples(["s1", "s2"])

    with pytest.raises(ValueError, match="result must be one of"):
        intersections.similarity(samples, result="bad")
    with pytest.raises(ValueError, match="non-negative integer"):
        intersections.similarity(samples, mismatches=-1)
    with pytest.raises(TypeError, match="non-negative integer"):
        intersections.similarity(samples, mismatches=1.5)
    with pytest.warns(DeprecationWarning, match="deprecated and ignored"):
        intersections.similarity(samples, by_freq=True)
