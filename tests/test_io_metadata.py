import warnings
import json

import pandas as pd

from repseq import io
from repseq import vdjtools


def test_read_ngsik_metadata_returns_empty_dataframe_for_missing_file(tmp_path):
    result = io.read_ngsik_metadata(str(tmp_path), verbose=False)

    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_read_yaml_metadata_warns_and_delegates(tmp_path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = io.read_yaml_metadata(str(tmp_path), verbose=False)

    assert isinstance(result, pd.DataFrame)
    assert result.empty
    assert any(item.category is DeprecationWarning for item in caught)
    assert any("read_ngsik_metadata" in str(item.message) for item in caught)


def test_open_json_report_reads_last_valid_json_record(tmp_path):
    report_filename = tmp_path / "sample.assemble.report.json"
    reports = [
        {"type": "assemblerReport", "clones": 1},
        {"type": "assemblerReport", "clones": 2},
        {"type": "assemblerReport", "clones": 3},
    ]
    report_filename.write_text(
        "\n".join(
            [
                json.dumps(reports[0]),
                "",
                "not a json record",
                json.dumps(reports[1], indent=2),
                json.dumps(reports[2]),
            ]
        )
    )

    assert io.open_json_report(report_filename)["clones"] == 3


def _write_mixcr_clonoset(path, count):
    pd.DataFrame(
        [
            {
                "cloneCount": count,
                "cloneFraction": 1.0,
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV1*01(100)",
                "allDHitsWithScore": "TRBD1*01(10)",
                "allJHitsWithScore": "TRBJ1*01(80)",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def test_save_to_vdjtools_uses_chain_in_output_filename(tmp_path):
    input_file = tmp_path / "sample.clones_TRB.tsv"
    output_folder = tmp_path / "vdjtools"
    _write_mixcr_clonoset(input_file, 5)
    samples = pd.DataFrame(
        [{"sample_id": "sample1", "chain": "TRB", "filename": str(input_file)}]
    )

    metadata = io.save_to_vdjtools(samples, str(output_folder))

    output_file = output_folder / "vdjtools.sample1.TRB.txt"
    assert output_file.exists()
    assert (output_folder / "metadata.txt").exists()
    assert metadata.loc[0, "#file.name"] == "vdjtools.sample1.TRB.txt"
    assert metadata.loc[0, "sample.id"] == "sample1"


def test_save_to_vdjtools_keeps_old_filename_without_chain_column(tmp_path):
    input_file = tmp_path / "sample.clones_TRB.tsv"
    output_folder = tmp_path / "vdjtools"
    _write_mixcr_clonoset(input_file, 5)
    samples = pd.DataFrame([{"sample_id": "sample1", "filename": str(input_file)}])

    io.save_to_vdjtools(samples, str(output_folder))

    assert (output_folder / "vdjtools.sample1.txt").exists()


def test_save_to_vdjtools_does_not_write_any_files_when_conflict_exists(tmp_path, capsys):
    input_file_1 = tmp_path / "sample1.clones_TRB.tsv"
    input_file_2 = tmp_path / "sample2.clones_TRB.tsv"
    output_folder = tmp_path / "vdjtools"
    output_folder.mkdir()
    _write_mixcr_clonoset(input_file_1, 5)
    _write_mixcr_clonoset(input_file_2, 7)
    conflict = output_folder / "vdjtools.sample1.TRB.txt"
    conflict.write_text("existing\n")
    samples = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "filename": str(input_file_1)},
            {"sample_id": "sample2", "chain": "TRB", "filename": str(input_file_2)},
        ]
    )

    result = io.save_to_vdjtools(samples, str(output_folder))

    captured = capsys.readouterr().out
    assert result is None
    assert "force_overwrite=True" in captured
    assert str(conflict) in captured
    assert conflict.read_text() == "existing\n"
    assert not (output_folder / "vdjtools.sample2.TRB.txt").exists()
    assert not (output_folder / "metadata.txt").exists()


def test_save_to_vdjtools_force_overwrite_writes_conflicting_files(tmp_path):
    input_file = tmp_path / "sample.clones_TRB.tsv"
    output_folder = tmp_path / "vdjtools"
    output_folder.mkdir()
    _write_mixcr_clonoset(input_file, 5)
    output_file = output_folder / "vdjtools.sample1.TRB.txt"
    output_file.write_text("existing\n")
    samples = pd.DataFrame(
        [{"sample_id": "sample1", "chain": "TRB", "filename": str(input_file)}]
    )

    io.save_to_vdjtools(samples, str(output_folder), force_overwrite=True)

    assert output_file.read_text() != "existing\n"


def test_vdjtools_save_to_vdjtools_warns_and_delegates(tmp_path):
    input_file = tmp_path / "sample.clones_TRB.tsv"
    output_folder = tmp_path / "vdjtools"
    _write_mixcr_clonoset(input_file, 5)
    samples = pd.DataFrame([{"sample_id": "sample1", "filename": str(input_file)}])

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        vdjtools.save_to_vdjtools(samples, str(output_folder))

    assert any(item.category is DeprecationWarning for item in caught)
    assert any("repseq.io.save_to_vdjtools" in str(item.message) for item in caught)
