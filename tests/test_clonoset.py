import pandas as pd
import zipfile

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


def test_read_clonoset_detects_comma_and_mixcr4_dtypes(tmp_path):
    path = tmp_path / "sample.csv"
    pd.DataFrame(
        {
            "cloneId": [1, 2],
            "readCount": [10, 5],
            "readFraction": [2 / 3, 1 / 3],
            "nSeqCDR3": ["00123", "00456"],
            "aaSeqCDR3": ["CASSLG", "CASSQG"],
        }
    ).to_csv(path, index=False)

    result = io.read_clonoset(path)

    assert str(result["cloneId"].dtype) == "Int64"
    assert result["readFraction"].dtype == "float64"
    assert str(result["nSeqCDR3"].dtype) == "string"
    assert result["nSeqCDR3"].tolist() == ["00123", "00456"]


def test_read_clonoset_accepts_mixcr_region_not_covered_values(tmp_path):
    path = tmp_path / "sample.tsv"
    pd.DataFrame(
        {
            "cloneId": [1],
            "readCount": [10],
            "readFraction": [1.0],
            "nSeqFR1": ["region_not_covered"],
            "minQualFR1": ["region_not_covered"],
            "nSeqCDR1": ["region_not_covered"],
            "minQualCDR1": ["region_not_covered"],
            "nSeqCDR3": ["region_not_covered"],
            "minQualCDR3": ["region_not_covered"],
            "aaSeqCDR3": ["region_not_covered"],
        }
    ).to_csv(path, sep="\t", index=False)

    result = io.read_clonoset(path)

    for column in (
        "nSeqFR1",
        "minQualFR1",
        "nSeqCDR1",
        "minQualCDR1",
        "nSeqCDR3",
        "minQualCDR3",
        "aaSeqCDR3",
    ):
        assert str(result[column].dtype) == "string"
        assert result.loc[0, column] == "region_not_covered"


def test_read_clonoset_detects_airr_types_inside_zip(tmp_path):
    airr_path = tmp_path / "sample.tsv"
    pd.DataFrame(
        {
            "sequence_id": ["clone1", "clone2"],
            "productive": ["T", "F"],
            "v_call": ["TRBV1*01", "TRBV2*01"],
            "j_call": ["TRBJ1*01", "TRBJ2*01"],
            "junction": ["TGTGCC", "TGTGCT"],
            "duplicate_count": [10, None],
        }
    ).to_csv(airr_path, sep="\t", index=False)
    zip_path = tmp_path / "sample.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.write(airr_path, arcname="inside.tsv")

    result = io.read_clonoset(zip_path)

    assert str(result["sequence_id"].dtype) == "string"
    assert str(result["productive"].dtype) == "boolean"
    assert result["productive"].tolist() == [True, False]
    assert str(result["duplicate_count"].dtype) == "Int64"
    assert pd.isna(result.loc[1, "duplicate_count"])


def test_detect_clonoset_format_signatures():
    assert io._detect_clonoset_format(["cloneId", "cloneCount", "cloneFraction"]) == "MiXCR3"
    assert io._detect_clonoset_format(["cloneId", "readCount", "readFraction"]) == "MiXCR4"
    assert io._detect_clonoset_format(["count", "freq", "cdr3nt", "cdr3aa"]) == "VDJtools"
    assert io._detect_clonoset_format(["sequence_id", "v_call", "j_call"]) == "AIRR"
    assert io._detect_clonoset_format(["SEQUENCE_ID", "SEQUENCE_INPUT", "V_CALL"]) == "IgBLAST"
    assert io._detect_clonoset_format(["#count", "CDR3nt", "CDR3aa", "V", "J"]) == "TRUST4"


def test_clonoset_normalize_freq_returns_new_object():
    clonoset = Clonoset(_mixcr_like_table())

    normalized = clonoset.normalize_freq()

    assert normalized is not clonoset
    assert normalized["freq"].tolist() == [10 / 15, 5 / 15]
