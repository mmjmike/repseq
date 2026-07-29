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


@pytest.mark.parametrize(
    ("overlap_type", "clonoset_dicts", "expected_columns", "expected_row"),
    [
        (
            "aaVJ",
            {
                "s1": {("AAA", "V1", "J1"): 5},
                "s2": {("AAA", "V1", "J1"): 3},
            },
            ["clonotype", "cdr3aa", "v", "j", "s1", "s2"],
            {
                "clonotype": "AAA|V1|J1",
                "cdr3aa": "AAA",
                "v": "V1",
                "j": "J1",
                "s1": 5,
                "s2": 3,
            },
        ),
        (
            "VJ",
            {
                "s1": {("V1", "J1"): 5},
                "s2": {("V1", "J1"): 3},
            },
            ["clonotype", "v", "j", "s1", "s2"],
            {
                "clonotype": "V1|J1",
                "v": "V1",
                "j": "J1",
                "s1": 5,
                "s2": 3,
            },
        ),
        (
            "VJlen",
            {
                "s1": {("V1", "J1", 15): 5},
                "s2": {("V1", "J1", 15): 3},
            },
            ["clonotype", "v", "j", "len", "s1", "s2"],
            {
                "clonotype": "V1|J1|15",
                "v": "V1",
                "j": "J1",
                "len": 15,
                "s1": 5,
                "s2": 3,
            },
        ),
    ],
)
def test_count_table_returns_string_clonotype_and_component_columns(
    monkeypatch,
    overlap_type,
    clonoset_dicts,
    expected_columns,
    expected_row,
):
    conversion_args = {}

    def fake_convert(*args, **kwargs):
        conversion_args.update(kwargs)
        return clonoset_dicts

    monkeypatch.setattr(intersections, "convert_clonosets_to_compact_dicts", fake_convert)
    monkeypatch.setattr(intersections, "run_parallel_calculation", _run_sequential)

    table = intersections.count_table(
        _dummy_samples(["s1", "s2"]),
        overlap_type=overlap_type,
        mismatches=2 if overlap_type in {"VJ", "VJlen"} else 0,
    )

    assert list(table.columns) == expected_columns
    assert isinstance(table.index, pd.RangeIndex)
    assert table.iloc[0].to_dict() == expected_row
    if overlap_type in {"VJ", "VJlen"}:
        assert conversion_args["strict"] is True


def test_count_table_uses_only_custom_clonotypes_in_dataframe_order(monkeypatch):
    clonoset_dicts = {
        "s1": {
            ("AAA", "V1", "J1"): 5,
            ("IGNORED", "V2", "J2"): 9,
        },
        "s2": {
            ("AAA", "V1", "J1"): 3,
        },
    }
    monkeypatch.setattr(
        intersections,
        "convert_clonosets_to_compact_dicts",
        lambda *args, **kwargs: clonoset_dicts,
    )
    monkeypatch.setattr(intersections, "run_parallel_calculation", _run_sequential)
    custom_clonotypes = pd.DataFrame(
        [
            {"cdr3aa": "ABSENT", "v": "V3", "j": "J3"},
            {"cdr3aa": "AAA", "v": "V1", "j": "J1"},
            {"cdr3aa": "AAA", "v": "V1", "j": "J1"},
        ]
    )

    table = intersections.count_table(
        _dummy_samples(["s1", "s2"]),
        overlap_type="aaVJ",
        custom_clonotypes_df=custom_clonotypes,
    )

    assert table.to_dict("records") == [
        {
            "clonotype": "ABSENT|V3|J3",
            "cdr3aa": "ABSENT",
            "v": "V3",
            "j": "J3",
            "s1": 0,
            "s2": 0,
        },
        {
            "clonotype": "AAA|V1|J1",
            "cdr3aa": "AAA",
            "v": "V1",
            "j": "J1",
            "s1": 5,
            "s2": 3,
        },
    ]


@pytest.mark.parametrize(
    ("overlap_type", "custom_row", "expected_clonotype"),
    [
        ("aa", {"cdr3aa": "AAA"}, ("AAA",)),
        ("aaV", {"cdr3aa": "AAA", "v": "V1"}, ("AAA", "V1")),
        ("aaVJ", {"cdr3aa": "AAA", "v": "V1", "j": "J1"}, ("AAA", "V1", "J1")),
        ("nt", {"cdr3nt": "GCT"}, ("GCT",)),
        ("ntV", {"cdr3nt": "GCT", "v": "V1"}, ("GCT", "V1")),
        ("ntVJ", {"cdr3nt": "GCT", "v": "V1", "j": "J1"}, ("GCT", "V1", "J1")),
        ("VJ", {"v": "V1", "j": "J1"}, ("V1", "J1")),
        ("VJlen", {"cdr3aa": "AAA", "v": "V1", "j": "J1"}, ("V1", "J1", 3)),
    ],
)
def test_custom_clonotype_columns_follow_overlap_type(
    overlap_type,
    custom_row,
    expected_clonotype,
):
    custom_clonotypes = pd.DataFrame([custom_row])

    assert intersections._clonotypes_from_custom_dataframe(
        custom_clonotypes,
        overlap_type,
    ) == [expected_clonotype]


def test_count_table_custom_clonotypes_reports_schema_errors(monkeypatch):
    def fail_if_clonosets_are_read(*args, **kwargs):
        raise AssertionError("clonosets should not be read after invalid custom input")

    monkeypatch.setattr(
        intersections,
        "convert_clonosets_to_compact_dicts",
        fail_if_clonosets_are_read,
    )

    with pytest.raises(ValueError) as error:
        intersections.count_table(
            _dummy_samples(["s1"]),
            overlap_type="aaVJ",
            custom_clonotypes_df=pd.DataFrame({"cdr3aa": ["AAA"], "v": ["V1"]}),
        )

    message = str(error.value)
    assert "overlap_type='aaVJ'" in message
    assert "Missing required columns: j" in message
    assert "Required columns: cdr3aa, v, j" in message
    assert "Available columns: cdr3aa, v" in message
    assert "Supported schemas:" in message


@pytest.mark.parametrize(
    ("custom_clonotypes", "error_type", "message"),
    [
        ([{"cdr3aa": "AAA"}], TypeError, "must be a pandas DataFrame"),
        (pd.DataFrame({"cdr3aa": [None]}), ValueError, "contains missing values"),
        (pd.DataFrame({"cdr3aa": [123]}), TypeError, "must be strings"),
    ],
)
def test_custom_clonotypes_validates_dataframe_values(
    custom_clonotypes,
    error_type,
    message,
):
    with pytest.raises(error_type, match=message):
        intersections._clonotypes_from_custom_dataframe(custom_clonotypes, "aa")


def test_tcrnet_returns_string_clonotype_and_component_columns(monkeypatch):
    monkeypatch.setattr(
        intersections,
        "pool_clonotypes_from_clonosets_df",
        lambda clonosets_df, cl_filter=None: clonosets_df.attrs["pool_name"],
    )

    def fake_prepare(clonoset, **kwargs):
        if clonoset == "experimental":
            return {(3, "V1", "J1"): [("AAA", 2), ("AAT", 1)]}
        return {(3, "V1", "J1"): [("AAA", 1)]}

    monkeypatch.setattr(intersections, "prepare_clonoset_for_intersection", fake_prepare)
    monkeypatch.setattr(intersections, "run_parallel_calculation", _run_sequential)

    experimental = _dummy_samples(["experimental"])
    experimental.attrs["pool_name"] = "experimental"
    control = _dummy_samples(["control"])
    control.attrs["pool_name"] = "control"

    table = intersections.tcrnet(
        experimental,
        control,
        overlap_type="aaVJ",
        mismatches=1,
    )

    assert list(table.columns[:4]) == ["clonotype", "cdr3aa", "v", "j"]
    assert "clone" not in table.columns
    assert set(table["clonotype"]) == {"AAA|V1|J1", "AAT|V1|J1"}
    assert all(isinstance(clonotype, str) for clonotype in table["clonotype"])
