import warnings
import json

import pandas as pd

from repseq import io


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
