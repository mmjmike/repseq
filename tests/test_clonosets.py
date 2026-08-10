import warnings

from repseq import clonosets


def test_find_all_mixcr_clonosets_finds_mixcr_exported_files(tmp_path):
    filename = tmp_path / "sample1.clones_TRB.tsv"
    filename.write_text("cloneId\tcloneCount\n0\t1\n")

    result = clonosets.find_all_mixcr_clonosets(str(tmp_path))

    assert result.loc[0, "sample_id"] == "sample1"
    assert result.loc[0, "chain"] == "TRB"
    assert result.loc[0, "filename"] == str(filename)


def test_find_all_mixcr_clonosets_finds_gzipped_exported_files(tmp_path):
    filenames = {
        tmp_path / "sample1.clonotypes.TRB.txt.gz",
        tmp_path / "sample2.clones_TRA.tsv.gz",
    }
    for filename in filenames:
        filename.write_bytes(b"")

    result = clonosets.find_all_mixcr_clonosets(str(tmp_path))

    assert set(result["sample_id"]) == {"sample1", "sample2"}
    assert set(result["chain"]) == {"TRB", "TRA"}
    assert set(result["filename"]) == {str(filename) for filename in filenames}


def test_find_all_mixcr_clonosets_extracts_mix_id(tmp_path):
    mixed_filename = tmp_path / "mix1.sample1.clones_TRA.tsv"
    plain_filename = tmp_path / "sample2.clones_TRB.tsv"
    mixed_filename.write_text("cloneId\tcloneCount\n0\t1\n")
    plain_filename.write_text("cloneId\tcloneCount\n0\t1\n")

    result = clonosets.find_all_mixcr_clonosets(str(tmp_path)).set_index("sample_id")

    assert result.loc["sample1", "mix_id"] == "mix1"
    assert result.loc["sample1", "chain"] == "TRA"
    assert result.loc["sample1", "filename"] == str(mixed_filename)
    assert result.loc["sample2", "mix_id"] is None


def test_find_all_mixcr_clonosets_keeps_mix_id_first_across_folders(tmp_path):
    plain_folder = tmp_path / "plain"
    mixed_folder = tmp_path / "mixed"
    plain_folder.mkdir()
    mixed_folder.mkdir()
    (plain_folder / "sample1.clones_TRB.tsv").write_text("")
    (mixed_folder / "mix1.sample2.clones_TRA.tsv").write_text("")

    result = clonosets.find_all_mixcr_clonosets(
        [str(plain_folder), str(mixed_folder)]
    )

    assert list(result.columns) == ["mix_id", "sample_id", "chain", "filename"]


def test_find_all_exported_clonosets_warns_and_delegates(tmp_path):
    filename = tmp_path / "sample1.clones_TRB.tsv"
    filename.write_text("cloneId\tcloneCount\n0\t1\n")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = clonosets.find_all_exported_clonosets(str(tmp_path))

    assert result.loc[0, "sample_id"] == "sample1"
    assert any(item.category is DeprecationWarning for item in caught)
    assert any("find_all_mixcr_clonosets" in str(item.message) for item in caught)
