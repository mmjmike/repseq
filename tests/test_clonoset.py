import pandas as pd

from repseq import Clonoset
from repseq import io
from repseq.clonoset import standardize_to_vdjtools_columns


def _mixcr_like_table():
    return pd.DataFrame(
        {
            "cloneCount": [10, 5],
            "cloneFraction": [0.667, 0.333],
            "nSeqCDR3": ["TGTGCC", "TGTGCT"],
            "aaSeqCDR3": ["CASSLG", "CASSQG"],
            "allVHitsWithScore": ["TRBV7-9*01(100)", "TRBV5-1*01(90)"],
            "allJHitsWithScore": ["TRBJ2-7*01(80)", "TRBJ1-1*01(70)"],
        }
    )


def test_clonoset_standardizes_mixcr_like_columns():
    clonoset = Clonoset(_mixcr_like_table(), sample_id="sample1", chain="TRB")

    assert clonoset.sample_id == "sample1"
    assert clonoset.chain == "TRB"
    assert {"count", "freq", "cdr3nt", "cdr3aa", "v", "j"}.issubset(
        clonoset.columns
    )
    assert clonoset["count"].tolist() == [10, 5]
    assert clonoset["cdr3aa"].tolist() == ["CASSLG", "CASSQG"]


def test_standardize_to_vdjtools_columns_can_extract_segment_names():
    table = standardize_to_vdjtools_columns(
        _mixcr_like_table(),
        extract_segments=True,
        copy=True,
    )

    assert table["v"].tolist() == ["TRBV7-9", "TRBV5-1"]
    assert table["j"].tolist() == ["TRBJ2-7", "TRBJ1-1"]


def test_read_clonoset_keeps_dataframe_default(tmp_path):
    path = tmp_path / "sample.tsv"
    _mixcr_like_table().to_csv(path, sep="\t", index=False)

    result = io.read_clonoset(path)

    assert isinstance(result, pd.DataFrame)
    assert not isinstance(result, Clonoset)
    assert "cloneCount" in result.columns


def test_read_clonoset_can_return_clonoset(tmp_path):
    path = tmp_path / "sample.tsv"
    _mixcr_like_table().to_csv(path, sep="\t", index=False)

    result = io.read_clonoset(
        path,
        as_clonoset=True,
        sample_id="sample1",
        chain="TRB",
        metadata={"subject": "subject1"},
    )

    assert isinstance(result, Clonoset)
    assert result.sample_id == "sample1"
    assert result.metadata == {"subject": "subject1"}
    assert {"count", "freq", "cdr3nt", "cdr3aa", "v", "j"}.issubset(
        result.columns
    )


def test_clonoset_normalize_freq_returns_new_object():
    clonoset = Clonoset(_mixcr_like_table())

    normalized = clonoset.normalize_freq()

    assert normalized is not clonoset
    assert normalized["freq"].tolist() == [10 / 15, 5 / 15]
