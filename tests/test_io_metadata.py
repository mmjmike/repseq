import warnings
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

from repseq import io
from repseq import vdjtools


def test_mixcr_import_does_not_require_requests(tmp_path):
    (tmp_path / "requests.py").write_text(
        'raise AttributeError("requests dependency is broken")\n'
    )
    project_root = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path), str(project_root), environment.get("PYTHONPATH", "")]
    )

    result = subprocess.run(
        [sys.executable, "-c", "from repseq import mixcr"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


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
        [
            {
                "sample_id": "sample1",
                "chain": "TRB",
                "filename": str(input_file),
                "group": "case",
                "subject": "subject1",
            }
        ]
    )

    metadata = io.save_to_vdjtools(samples, str(output_folder))
    metadata_from_file = pd.read_csv(output_folder / "metadata.txt", sep="\t")

    output_file = output_folder / "vdjtools.sample1.TRB.txt"
    assert output_file.exists()
    assert (output_folder / "metadata.txt").exists()
    assert metadata.loc[0, "#file.name"] == "vdjtools.sample1.TRB.txt"
    assert metadata.loc[0, "sample.id"] == "sample1"
    assert "original_filename" in metadata.columns
    assert "filename" not in metadata.columns
    assert metadata.loc[0, "original_filename"] == str(input_file)
    assert metadata.loc[0, "group"] == "case"
    assert metadata.loc[0, "subject"] == "subject1"
    assert metadata_from_file.loc[0, "original_filename"] == str(input_file)
    assert metadata_from_file.loc[0, "group"] == "case"
    assert metadata_from_file.loc[0, "subject"] == "subject1"


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


def test_save_to_vdjtools_force_overwrite_merges_existing_metadata_for_existing_files(tmp_path):
    old_input_file = tmp_path / "old.clones_TRB.tsv"
    new_input_file = tmp_path / "new.clones_TRA.tsv"
    missing_input_file = tmp_path / "missing.clones_TRG.tsv"
    output_folder = tmp_path / "vdjtools"
    output_folder.mkdir()
    _write_mixcr_clonoset(old_input_file, 5)
    _write_mixcr_clonoset(new_input_file, 7)
    _write_mixcr_clonoset(missing_input_file, 9)
    old_output = output_folder / "vdjtools.old.TRB.txt"
    old_output.write_text("old clonoset\n")
    stale_output_name = "vdjtools.missing.TRG.txt"
    existing_metadata = pd.DataFrame(
        [
            {
                "#file.name": old_output.name,
                "sample.id": "old",
                "sample_id": "old",
                "chain": "TRB",
                "filename": str(old_input_file),
                "group": "old_group",
            },
            {
                "#file.name": stale_output_name,
                "sample.id": "missing",
                "sample_id": "missing",
                "chain": "TRG",
                "filename": str(missing_input_file),
                "group": "stale_group",
            },
        ]
    )
    existing_metadata.to_csv(output_folder / "metadata.txt", sep="\t", index=False)
    samples = pd.DataFrame(
        [
            {
                "sample_id": "new",
                "chain": "TRA",
                "filename": str(new_input_file),
                "group": "new_group",
            }
        ]
    )

    metadata = io.save_to_vdjtools(samples, str(output_folder), force_overwrite=True)

    assert set(metadata["#file.name"]) == {
        "vdjtools.old.TRB.txt",
        "vdjtools.new.TRA.txt",
    }
    old_row = metadata.loc[metadata["sample_id"] == "old"].iloc[0]
    new_row = metadata.loc[metadata["sample_id"] == "new"].iloc[0]
    assert old_row["original_filename"] == str(old_input_file)
    assert old_row["group"] == "old_group"
    assert new_row["original_filename"] == str(new_input_file)
    assert new_row["group"] == "new_group"


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


def test_save_to_airr_writes_populated_columns_and_original_metadata_names(tmp_path):
    input_folder = tmp_path / "inputs"
    input_folder.mkdir()
    input_file = input_folder / "sample.clones_TRB.tsv"
    pd.DataFrame(
        [
            {
                "cloneId": 7,
                "readCount": 12,
                "uniqueUMICount": 3,
                "targetSequences": "TGTGCCAGC",
                "targetQualities": "IIIIIIIII",
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV1*01(100)",
                "allJHitsWithScore": "TRBJ1*01(80)",
            }
        ]
    ).to_csv(input_file, sep="\t", index=False)
    output_folder = tmp_path / "airr"
    samples = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "TRB",
                "filename": str(input_file),
                "group": "case",
            }
        ]
    )

    metadata = io.save_to_airr(samples, output_folder)

    output_name = "sample.clones_TRB.tsv.airr.tsv"
    output_file = output_folder / output_name
    airr = pd.read_csv(output_file, sep="\t")
    metadata_from_file = pd.read_csv(output_folder / "metadata.csv")

    assert list(metadata.columns) == list(samples.columns)
    assert metadata.loc[0, "filename"] == output_name
    assert metadata.loc[0, "sample_id"] == "sample1"
    assert metadata.loc[0, "group"] == "case"
    pd.testing.assert_frame_equal(metadata, metadata_from_file)
    assert airr.loc[0, "sequence_id"] == 7
    assert airr.loc[0, "sequence"] == "TGTGCCAGC"
    assert airr.loc[0, "quality"] == "IIIIIIIII"
    assert airr.loc[0, "v_call"] == "TRBV1*01"
    assert airr.loc[0, "j_call"] == "TRBJ1*01"
    assert airr.loc[0, "junction"] == "TGTGCC"
    assert airr.loc[0, "junction_aa"] == "CASSLG"
    assert airr.loc[0, "duplicate_count"] == 12
    assert airr.loc[0, "umi_count"] == 3
    assert bool(airr.loc[0, "productive"])
    assert airr.loc[0, "junction_length"] == 6
    assert airr.loc[0, "junction_aa_length"] == 6
    assert airr.loc[0, "locus"] == "TRB"
    assert "d_call" not in airr.columns
    assert "germline_alignment" not in airr.columns


def test_save_to_airr_stops_before_writing_when_any_target_exists(tmp_path, capsys):
    input_file = tmp_path / "sample.tsv"
    _write_mixcr_clonoset(input_file, 5)
    output_folder = tmp_path / "airr"
    output_folder.mkdir()
    conflict = output_folder / "sample.tsv.airr.tsv"
    conflict.write_text("existing\n")
    samples = pd.DataFrame([{"sample_id": "sample1", "filename": str(input_file)}])

    result = io.save_to_airr(samples, output_folder)

    assert result is None
    assert conflict.read_text() == "existing\n"
    assert not (output_folder / "metadata.csv").exists()
    assert "force_overwrite=True" in capsys.readouterr().out
