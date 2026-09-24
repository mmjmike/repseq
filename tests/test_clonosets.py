import warnings

import pandas as pd
import pytest

from repseq import clonosets
from repseq.clone_filter import Filter


@pytest.mark.parametrize("cpu", [1, 2])
def test_pool_clonotypes_reads_samples_in_order_with_filter_options(tmp_path, cpu):
    samples = []
    for sample_id, count in [("first", 2), ("second", 3)]:
        filename = tmp_path / f"{sample_id}.tsv"
        pd.DataFrame({
            "count": [count], "freq": [1.0], "cdr3nt": ["TGT"],
            "cdr3aa": ["C"], "v": ["TRBV1"], "j": ["TRBJ1"],
            "extra": [sample_id],
        }).to_csv(filename, sep="\t", index=False)
        samples.append((sample_id, str(filename)))

    result = clonosets.pool_clonotypes_from_clonosets_df(
        pd.DataFrame(samples, columns=["sample_id", "filename"]),
        cl_filter=Filter(convert=False), cpu=cpu,
    )

    assert result["sample_id"].tolist() == ["first", "second"]
    assert result["extra"].tolist() == ["first", "second"]
    assert result["count"].tolist() == [2, 3]


def test_pool_clonotypes_passes_cpu_to_parallel_runner(monkeypatch):
    calls = []

    def fake_parallel(function, tasks, program_name, **kwargs):
        calls.append((program_name, kwargs, len(tasks)))
        return [function(task) for task in tasks]

    monkeypatch.setattr(clonosets, "run_parallel_calculation", fake_parallel)
    monkeypatch.setattr(
        clonosets, "read_clonoset", lambda filename: pd.DataFrame({"count": [1]})
    )
    clonosets.pool_clonotypes_from_clonosets_df(
        pd.DataFrame({"sample_id": ["sample"], "filename": ["sample.tsv"]}),
        cl_filter=Filter(convert=False), cpu=3,
    )

    assert calls == [("Pooling clonotypes", {
        "object_name": "samples", "verbose": False, "cpu": 3,
    }, 1)]


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


def test_find_all_mixcr_clonosets_keeps_sample_id_first_across_folders(tmp_path):
    plain_folder = tmp_path / "plain"
    mixed_folder = tmp_path / "mixed"
    plain_folder.mkdir()
    mixed_folder.mkdir()
    (plain_folder / "sample1.clones_TRB.tsv").write_text("")
    (mixed_folder / "mix1.sample2.clones_TRA.tsv").write_text("")

    result = clonosets.find_all_mixcr_clonosets(
        [str(plain_folder), str(mixed_folder)]
    )

    assert list(result.columns) == ["sample_id", "mix_id", "chain", "filename"]


def test_find_all_exported_clonosets_warns_and_delegates(tmp_path):
    filename = tmp_path / "sample1.clones_TRB.tsv"
    filename.write_text("cloneId\tcloneCount\n0\t1\n")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = clonosets.find_all_exported_clonosets(str(tmp_path))

    assert result.loc[0, "sample_id"] == "sample1"
    assert any(item.category is DeprecationWarning for item in caught)
    assert any("find_all_mixcr_clonosets" in str(item.message) for item in caught)
