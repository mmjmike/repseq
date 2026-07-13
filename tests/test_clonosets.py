import warnings

from repseq import clonosets


def test_find_all_mixcr_clonosets_finds_mixcr_exported_files(tmp_path):
    filename = tmp_path / "sample1.clones_TRB.tsv"
    filename.write_text("cloneId\tcloneCount\n0\t1\n")

    result = clonosets.find_all_mixcr_clonosets(str(tmp_path))

    assert result.loc[0, "sample_id"] == "sample1"
    assert result.loc[0, "chain"] == "TRB"
    assert result.loc[0, "filename"] == str(filename)


def test_find_all_exported_clonosets_warns_and_delegates(tmp_path):
    filename = tmp_path / "sample1.clones_TRB.tsv"
    filename.write_text("cloneId\tcloneCount\n0\t1\n")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = clonosets.find_all_exported_clonosets(str(tmp_path))

    assert result.loc[0, "sample_id"] == "sample1"
    assert any(item.category is DeprecationWarning for item in caught)
    assert any("find_all_mixcr_clonosets" in str(item.message) for item in caught)
