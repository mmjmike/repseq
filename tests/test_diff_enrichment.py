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
