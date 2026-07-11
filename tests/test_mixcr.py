import inspect
import json
import shlex
import sys
import types

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from repseq import mixcr


def test_mixcr4_analyze_batch_preserves_custom_tag_pattern(tmp_path):
    tag_pattern = r"^N{0:2}tggtatcaacgcagagt(SMPL:N{5})(UMI:N{14})N{1}gctN{16}(R1:*)\^N{20}(R2:*)"
    sample_df = pd.DataFrame(
        [
            {
                "sample_id": "sample_1",
                "R1": "R1.fastq.gz",
                "R2": "R2.fastq.gz",
                "tag_pattern": tag_pattern,
            }
        ]
    )

    jobs = mixcr.mixcr4_analyze_batch(
        sample_df,
        str(tmp_path),
        command_template="mixcr analyze test-preset -f r1 r2 output_prefix",
        mixcr_path="echo",
        memory=16,
        custom_tag_pattern_column="tag_pattern",
        backend="local",
    )

    command_parts = shlex.split(jobs.loc[0, "command"])
    tag_pattern_index = command_parts.index("--tag-pattern") + 1
    assert command_parts[tag_pattern_index] == tag_pattern
    assert jobs.loc[0, "status"] == "finished"
    assert (tmp_path / "mixcr_analyze_batch.log").exists()
    assert (tmp_path / "logs" / "mixcr_analyze_batch_jobs.csv").exists()


def test_mixcr_public_batch_functions_do_not_expose_max_workers():
    assert "max_workers" not in inspect.signature(mixcr.mixcr4_analyze_batch).parameters
    assert "max_workers" not in inspect.signature(mixcr.mixcr_7genes_run_batch).parameters
    assert "max_workers" not in inspect.signature(mixcr.mixcr4_reports).parameters


def test_mixcr4_reports_uses_batch_log_without_csv_table(tmp_path):
    jobs = mixcr.mixcr4_reports(str(tmp_path), mixcr_path="echo", backend="local", memory=16)

    assert set(jobs["status"]) == {"finished"}
    batch_log = tmp_path / "mixcr_reports_slurm_batch.log"
    assert batch_log.exists()
    assert batch_log.read_text().splitlines()[0] == "# MiXCR 4 Reports"
    assert not (tmp_path / "logs" / "mixcr_reports_slurm_batch_jobs.csv").exists()


def test_check_batch_progress_prints_failed_jobs(tmp_path, capsys):
    batch_filename = tmp_path / "mixcr_analyze_batch.log"
    table = pd.DataFrame(
        [
            {
                "jobname": "mixcr_analyze_ok",
                "sample_id": "ok",
                "backend": "local",
                "job_id": "",
                "status": "finished",
                "returncode": "0",
                "cwd": str(tmp_path),
                "log_filename": str(tmp_path / "ok.log"),
                "command": "true",
            },
            {
                "jobname": "mixcr_analyze_bad",
                "sample_id": "bad",
                "backend": "slurm",
                "job_id": "12345",
                "status": "failed",
                "returncode": "1",
                "cwd": str(tmp_path),
                "log_filename": str(tmp_path / "bad.log"),
                "command": "false",
            },
        ]
    )
    mixcr._write_batch_table(str(batch_filename), table, program_name="MIXCR4 Analyze Batch")

    mixcr.check_batch_progress(str(tmp_path))

    captured = capsys.readouterr().out
    assert "Failed jobs:" in captured
    assert "mixcr_analyze_bad" in captured
    assert "12345" in captured


def test_check_batch_progress_updates_slurm_status_from_query(tmp_path, monkeypatch):
    batch_filename = tmp_path / "mixcr_analyze_batch.log"
    table = pd.DataFrame(
        [
            {
                "jobname": "mixcr_analyze_sample",
                "sample_id": "sample",
                "backend": "slurm",
                "job_id": "42",
                "status": "submitted",
                "returncode": "",
                "cwd": str(tmp_path),
                "log_filename": str(tmp_path / "sample.log"),
                "command": "mixcr",
            }
        ]
    )
    mixcr._write_batch_table(str(batch_filename), table, program_name="MIXCR4 Analyze Batch")

    monkeypatch.setattr(mixcr, "_query_squeue", lambda job_ids: {})
    monkeypatch.setattr(mixcr, "_query_sacct", lambda job_ids: {"42": {"status": "finished", "returncode": "0"}})
    monkeypatch.setattr(mixcr, "_query_scontrol", lambda job_ids: {})

    mixcr.check_batch_progress(str(tmp_path))

    _, updated_table = mixcr._read_batch_table(str(batch_filename))
    assert updated_table.loc[0, "status"] == "finished"
    assert str(updated_table.loc[0, "returncode"]) == "0"


def test_show_report_images_reports_missing_images(tmp_path, monkeypatch, capsys):
    display_module = types.ModuleType("IPython.display")
    display_module.Image = lambda filename: ("Image", filename)
    display_module.SVG = lambda filename: ("SVG", filename)
    display_module.display = lambda obj: None
    ipython_module = types.ModuleType("IPython")
    ipython_module.display = display_module

    monkeypatch.setitem(sys.modules, "IPython", ipython_module)
    monkeypatch.setitem(sys.modules, "IPython.display", display_module)

    mixcr.show_report_images(str(tmp_path))

    captured = capsys.readouterr().out
    assert "No alignQc image found" in captured
    assert "No chainQc image found" in captured


def _write_align_report(folder, sample_id, chains):
    report = {
        "totalReadsProcessed": 100,
        "aligned": 80,
        "overlappedAligned": 70,
        "notAlignedReasons": {
            "NoHits": 10,
            "NoCDR3Parts": 2,
            "NoVHits": 1,
            "NoJHits": 1,
            "VAndJOnDifferentTargets": 1,
            "LowTotalScore": 1,
            "NoBarcode": 4,
        },
        "chainUsage": {"chains": chains},
    }
    (folder / f"{sample_id}.align.report.json").write_text(json.dumps(report) + "\n")


def _write_assemble_report(folder, sample_id, chains):
    report = {
        "type": "assemblerReport",
        "totalReadsProcessed": 100,
        "clonalChainUsage": {"type": "chainUsage", "chains": chains},
    }
    (folder / f"{sample_id}.assemble.report.json").write_text(json.dumps(report) + "\n")


def test_show_qc_plot_align_uses_two_legend_columns(tmp_path, monkeypatch):
    _write_align_report(
        tmp_path,
        "sample1",
        {"TRB": {"total": 80, "nonFunctional": 5, "hasStops": 2, "isOOF": 3}},
    )
    legend_calls = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    def fake_legend(*args, **kwargs):
        legend_calls.append(kwargs)

    monkeypatch.setattr(Figure, "legend", fake_legend)

    mixcr.show_qc_plot(str(tmp_path), chart_type="align")

    assert legend_calls[-1]["ncol"] == 2


def test_show_qc_plot_chains_uses_chain_count_legend_columns(tmp_path, monkeypatch):
    _write_assemble_report(
        tmp_path,
        "sample1",
        {
            "TRA": {"total": 40, "nonFunctional": 5, "hasStops": 2, "isOOF": 3},
            "TRB": {"total": 60, "nonFunctional": 6, "hasStops": 4, "isOOF": 2},
        },
    )
    legend_calls = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    def fake_legend(*args, **kwargs):
        legend_calls.append(kwargs)

    monkeypatch.setattr(Figure, "legend", fake_legend)

    mixcr.show_qc_plot(str(tmp_path), chart_type="chains")

    assert legend_calls[-1]["ncol"] == 2


def test_show_qc_plot_chains_reads_clonal_chain_usage_from_assemble_report(tmp_path, monkeypatch):
    _write_align_report(
        tmp_path,
        "sample1",
        {"IGH": {"total": 100, "nonFunctional": 0, "hasStops": 0, "isOOF": 0}},
    )
    _write_assemble_report(
        tmp_path,
        "sample1",
        {"TRB": {"total": 60, "nonFunctional": 6, "hasStops": 4, "isOOF": 2}},
    )
    plotted_columns = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    def fake_barh(*args, **kwargs):
        plotted_columns.append(kwargs["label"])

    monkeypatch.setattr(Axes, "barh", fake_barh)

    mixcr.show_qc_plot(str(tmp_path), chart_type="chains")

    assert plotted_columns == ["TRB", "TRB (OOF)", "TRB (stops)"]


def test_show_qc_plot_chains_requires_assemble_report(tmp_path, monkeypatch, capsys):
    _write_align_report(
        tmp_path,
        "sample1",
        {"TRB": {"total": 60, "nonFunctional": 6, "hasStops": 4, "isOOF": 2}},
    )
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    mixcr.show_qc_plot(str(tmp_path), chart_type="chains")

    assert "No chains report data found" in capsys.readouterr().out


def test_show_qc_plot_coverage_accepts_processing_table(monkeypatch):
    processing_table = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "extracted_chain": "TRB",
                "reads_per_umi": 4.5,
                "overseq_threshold": 7,
            }
        ]
    )
    legend_calls = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    def fake_legend(*args, **kwargs):
        legend_calls.append(kwargs)

    monkeypatch.setattr(Figure, "legend", fake_legend)

    mixcr.show_qc_plot(".", chart_type="coverage", processing_table=processing_table)

    assert legend_calls[-1]["title"] == "Coverage"


def test_show_qc_plot_coverage_reads_refine_reports_directly(tmp_path, monkeypatch):
    report = {
        "correctionReport": {
            "outputRecords": 45,
            "steps": [{"outputDiversity": 12}],
            "filterReport": {
                "numberOfGroupsAccepted": 10,
                "operatorReports": [
                    {"operatorReport": {"threshold": 7}},
                ],
            },
        }
    }
    (tmp_path / "sample1.refine.report.json").write_text(json.dumps(report) + "\n")
    legend_calls = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)
    monkeypatch.setattr(
        mixcr,
        "get_processing_table",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("coverage plot should not build full processing table")
        ),
    )

    def fake_legend(*args, **kwargs):
        legend_calls.append(kwargs)

    monkeypatch.setattr(Figure, "legend", fake_legend)

    mixcr.show_qc_plot(str(tmp_path), chart_type="coverage")

    assert legend_calls[-1]["title"] == "Coverage"
