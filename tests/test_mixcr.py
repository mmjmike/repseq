import inspect
import json
import shlex
import sys
import types

import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from repseq import mixcr


SMALL_PROCESSING_TABLE_COLUMNS = [
    "sample_id",
    "extracted_chain",
    "reads_total",
    "reads_with_umi_pc",
    "reads_aligned_pc",
    "reads_overlapped_aln_pc",
    "reads_per_umi",
    "overseq_threshold",
    "clones_func",
    "umi_in_func_clones",
]


def test_get_processing_table_small_returns_selected_columns(tmp_path):
    table = mixcr.get_processing_table(str(tmp_path), small=True)

    assert table.columns.tolist() == SMALL_PROCESSING_TABLE_COLUMNS


def test_get_processing_table_small_applies_to_folder_lists(tmp_path):
    first_folder = tmp_path / "first"
    second_folder = tmp_path / "second"
    first_folder.mkdir()
    second_folder.mkdir()

    table = mixcr.get_processing_table(
        [str(first_folder), str(second_folder)],
        small=True,
    )

    assert table.columns.tolist() == SMALL_PROCESSING_TABLE_COLUMNS


def test_get_processing_table_uses_mix_and_sample_ids_for_reports(
    tmp_path, monkeypatch
):
    clonosets = pd.DataFrame(
        [{
            "sample_id": "sample1",
            "mix_id": "mix1",
            "chain": "TRA",
            "filename": str(tmp_path / "mix1.sample1.clones_TRA.tsv"),
        }]
    )
    report_calls = []

    monkeypatch.setattr(
        mixcr,
        "find_all_exported_clonosets_in_folder",
        lambda folder, chain=None: clonosets,
    )

    def read_report(sample_id, folder, report_type):
        report_calls.append((sample_id, report_type))
        if report_type == "refine":
            raise FileNotFoundError
        if report_type == "align":
            return {
                "totalReadsProcessed": 100,
                "notAlignedReasons": {"NoBarcode": 0},
                "aligned": 80,
                "overlappedAligned": 40,
            }
        return {"clones": 2, "readsInClones": 10}

    monkeypatch.setattr(mixcr, "read_json_report", read_report)
    monkeypatch.setattr(
        mixcr,
        "read_clonoset",
        lambda filename: pd.DataFrame({"readCount": [6, 4]}),
    )
    monkeypatch.setattr(mixcr, "filter_by_functionality", lambda clonoset: clonoset)

    table = mixcr.get_processing_table(
        str(tmp_path),
        show_offtarget=True,
        small=True,
    )

    assert table.columns.tolist() == [
        "sample_id", "mix_id", *SMALL_PROCESSING_TABLE_COLUMNS[1:]
    ]
    assert table.loc[0, "sample_id"] == "sample1"
    assert table.loc[0, "mix_id"] == "mix1"
    assert report_calls == [
        ("mix1", "align"),
        ("mix1.sample1", "refine"),
        ("mix1.sample1", "assemble"),
    ]


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


def test_mixcr4_analyze_batch_groups_sample_barcoded_files(tmp_path):
    sample_df = pd.DataFrame(
        [
            {
                "sample_id": "sample_1",
                "R1": "mix_R1.fastq.gz",
                "R2": "mix_R2.fastq.gz",
                "mix_id": "mix_1",
                "SMPL": "AAAAA",
                "miNNNPattern": "",
            },
            {
                "sample_id": "sample_2",
                "R1": "mix_R1.fastq.gz",
                "R2": "mix_R2.fastq.gz",
                "mix_id": "mix_1",
                "SMPL": "CCCCC",
                "miNNNPattern": None,
            },
            {
                "sample_id": "sample_3",
                "R1": "sample_3_R1.fastq.gz",
                "R2": "sample_3_R2.fastq.gz",
                "mix_id": "",
                "SMPL": "",
                "miNNNPattern": "",
            },
        ]
    )

    jobs = mixcr.mixcr4_analyze_batch(
        sample_df,
        str(tmp_path),
        command_template="mixcr analyze test-preset -f r1 r2 output_prefix",
        mixcr_path="echo",
        memory=16,
        backend="local",
    )

    assert jobs["jobname"].tolist() == ["mixcr_analyze_mix_1", "mixcr_analyze_sample_3"]
    mix_command = shlex.split(jobs.loc[jobs["jobname"] == "mixcr_analyze_mix_1", "command"].iloc[0])
    tag_pattern_index = mix_command.index("--tag-pattern") + 1
    sample_table_index = mix_command.index("--sample-table") + 1
    assert mix_command[tag_pattern_index] == mixcr.DEFAULT_SAMPLE_BARCODE_TAG_PATTERN
    assert "--split-by-sample" in mix_command
    assert mix_command[sample_table_index] == str(tmp_path / "mix_1_samples.tsv")
    assert mix_command[-3:] == ["mix_R1.fastq.gz", "mix_R2.fastq.gz", "mix_1"]

    barcode_table = pd.read_csv(tmp_path / "mix_1_samples.tsv", sep="	", keep_default_na=False)
    assert barcode_table.to_dict("records") == [
        {"Sample": "sample_1", "TagPattern": "", "SMPL": "AAAAA"},
        {"Sample": "sample_2", "TagPattern": "", "SMPL": "CCCCC"},
    ]


def test_mixcr4_analyze_batch_uses_minnn_pattern_for_sample_barcoded_files(tmp_path):
    tag_pattern = r"^(SMPL:N{6})(UMI:N{12})(R1:*)\^(R2:*)"
    sample_df = pd.DataFrame(
        [
            {
                "sample_id": "sample_1",
                "R1": "mix_R1.fastq.gz",
                "R2": "mix_R2.fastq.gz",
                "mix_id": "mix_1",
                "SMPL": "AAAAAA",
                "miNNNPattern": tag_pattern,
            },
            {
                "sample_id": "sample_2",
                "R1": "mix_R1.fastq.gz",
                "R2": "mix_R2.fastq.gz",
                "mix_id": "mix_1",
                "SMPL": "CCCCCC",
                "miNNNPattern": tag_pattern,
            },
        ]
    )

    jobs = mixcr.mixcr4_analyze_batch(
        sample_df,
        str(tmp_path),
        mixcr_path="echo",
        memory=16,
        backend="local",
    )

    command_parts = shlex.split(jobs.loc[0, "command"])
    assert command_parts[command_parts.index("--tag-pattern") + 1] == tag_pattern


@pytest.mark.parametrize(
    ("column", "values", "message"),
    [
        ("mix_id", None, "mix_id"),
        ("mix_id", ["mix_1", "mix_2"], "same non-empty mix_id"),
        ("SMPL", None, "SMPL"),
        ("SMPL", ["AAAAA", "AAAAA"], "distinct"),
        ("SMPL", ["AAAAA", "CCCC"], "same length"),
        ("SMPL", ["AAAAA", "CCCNC"], "only ATGC"),
        ("SMPL", [123, "CCCCC"], "must be strings"),
        ("miNNNPattern", ["pattern_1", "pattern_2"], "must be the same"),
    ],
)
def test_mixcr4_analyze_batch_validates_sample_barcoded_rows(tmp_path, column, values, message):
    sample_df = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "R1": ["mix_R1.fastq.gz", "mix_R1.fastq.gz"],
            "R2": ["mix_R2.fastq.gz", "mix_R2.fastq.gz"],
            "mix_id": ["mix_1", "mix_1"],
            "SMPL": ["AAAAA", "CCCCC"],
            "miNNNPattern": ["", ""],
        }
    )
    if values is None:
        sample_df = sample_df.drop(columns=column)
    else:
        sample_df[column] = values

    with pytest.raises(ValueError, match=message):
        mixcr.mixcr4_analyze_batch(sample_df, str(tmp_path), mixcr_path="echo", memory=16)


def test_mixcr4_analyze_batch_passes_constraint_to_slurm(tmp_path, monkeypatch):
    captured = {}

    def fake_submit(*args, **kwargs):
        captured["constraint"] = kwargs["constraint"]
        return b"Submitted batch job 42\n", b""

    monkeypatch.setattr(mixcr, "run_slurm_command_from_jupyter", fake_submit)
    sample_df = pd.DataFrame(
        [{"sample_id": "sample_1", "R1": "R1.fastq.gz", "R2": "R2.fastq.gz"}]
    )

    jobs = mixcr.mixcr4_analyze_batch(
        sample_df,
        str(tmp_path),
        mixcr_path="mixcr",
        memory=16,
        backend="slurm",
        constraint="hpc",
    )

    assert captured["constraint"] == "hpc"
    assert str(jobs.loc[0, "job_id"]) == "42"
    assert jobs.loc[0, "status"] == "submitted"


def test_find_alleles_groups_clns_files_by_donor_and_exports_results(tmp_path, capsys):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    for filename in ["sample_1.clns", "mix.sample_2.clns", "sample_10.clns", "unlisted.clns"]:
        (input_dir / filename).touch()
    sample_df = pd.DataFrame([
        {"sample_id": "sample_1", "donor_id": "donor_a"},
        {"sample_id": "sample_2", "donor_id": "donor_a"},
        {"sample_id": "sample_10", "donor_id": "donor_b"},
        {"sample_id": "missing", "donor_id": "donor_c"},
    ])

    jobs = mixcr.find_alleles(
        str(input_dir),
        str(output_dir),
        "echo",
        sample_df,
        backend="local",
    )

    assert jobs["jobname"].tolist() == [
        "mixcr_find_alleles_donor_a",
        "mixcr_find_alleles_donor_b",
    ]
    donor_a_command = jobs.loc[
        jobs["jobname"] == "mixcr_find_alleles_donor_a", "command"
    ].iloc[0]
    command_parts = [shlex.split(command) for command in donor_a_command.split(" && ")]
    assert command_parts[0][:4] == ["echo", "-Xmx32g", "findAlleles", "-f"]
    assert command_parts[0][-2:] == [
        str(input_dir / "mix.sample_2.clns"),
        str(input_dir / "sample_1.clns"),
    ]
    assert command_parts[1][-2:] == [
        str(output_dir / "mix.sample_2.clns"),
        str(output_dir / "mix.sample_2.clones.tsv"),
    ]
    assert command_parts[2][-2:] == [
        str(output_dir / "sample_1.clns"),
        str(output_dir / "sample_1.clones.tsv"),
    ]
    assert (output_dir / "find_alleles_batch.log").exists()
    assert (output_dir / "logs" / "find_alleles_batch_jobs.csv").exists()
    assert (output_dir / "logs" / "mixcr_find_alleles_donor_a.log").exists()

    output = capsys.readouterr().out
    assert "Donor donor_a: 2 .clns file(s)" in output
    assert str(input_dir / "sample_1.clns") in output
    assert str(input_dir / "mix.sample_2.clns") in output
    assert "Donor donor_c: 0 .clns file(s)" in output


def test_find_alleles_uses_default_slurm_resources(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    (input_dir / "sample_1.clns").touch()
    captured = {}

    def fake_submit(command, jobname, cpus, time_estimate, memory, **kwargs):
        captured.update({
            "command": command,
            "jobname": jobname,
            "cpus": cpus,
            "time_estimate": time_estimate,
            "memory": memory,
            "log_filename": kwargs["log_filename"],
        })
        return b"Submitted batch job 42\n", b""

    monkeypatch.setattr(mixcr, "run_slurm_command_from_jupyter", fake_submit)
    jobs = mixcr.find_alleles(
        str(input_dir),
        str(output_dir),
        sample_df=pd.DataFrame([{"sample_id": "sample_1", "donor_id": "donor_a"}]),
        backend="slurm",
    )

    assert captured["jobname"] == "mixcr_find_alleles_donor_a"
    assert captured["cpus"] == 4
    assert captured["time_estimate"] == 0.5
    assert captured["memory"] == 32
    assert captured["log_filename"] == str(
        output_dir / "logs" / "mixcr_find_alleles_donor_a.log"
    )
    assert "findAlleles" in captured["command"]
    assert "exportClones" in captured["command"]
    assert jobs.loc[0, "status"] == "submitted"
    assert str(jobs.loc[0, "job_id"]) == "42"


@pytest.mark.parametrize("missing_column", ["sample_id", "donor_id"])
def test_find_alleles_requires_sample_and_donor_columns(tmp_path, missing_column):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    sample_df = pd.DataFrame([{"sample_id": "sample_1", "donor_id": "donor_a"}])

    with pytest.raises(ValueError, match=missing_column):
        mixcr.find_alleles(
            str(input_dir),
            str(tmp_path / "output"),
            sample_df=sample_df.drop(columns=missing_column),
        )


def test_mixcr_public_batch_functions_do_not_expose_max_workers():
    assert "max_workers" not in inspect.signature(mixcr.mixcr4_analyze_batch).parameters
    assert "max_workers" not in inspect.signature(mixcr.find_alleles).parameters
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

    assert plotted_columns == ["TRB", "TRB (stops)", "TRB (OOF)"]


def test_show_qc_plot_chains_requires_assemble_report(tmp_path, monkeypatch, capsys):
    _write_align_report(
        tmp_path,
        "sample1",
        {"TRB": {"total": 60, "nonFunctional": 6, "hasStops": 4, "isOOF": 2}},
    )
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    mixcr.show_qc_plot(str(tmp_path), chart_type="chains")

    assert "No chains report data found" in capsys.readouterr().out


def test_show_qc_plot_chains_plots_samples_missing_some_chain_columns(tmp_path, monkeypatch):
    _write_assemble_report(
        tmp_path,
        "sample_TRA_TRB",
        {
            "TRA": {"total": 40, "nonFunctional": 5, "hasStops": 2, "isOOF": 3},
            "TRB": {"total": 60, "nonFunctional": 6, "hasStops": 4, "isOOF": 2},
        },
    )
    _write_assemble_report(
        tmp_path,
        "sample_TRB_only",
        {"TRB": {"total": 15562, "nonFunctional": 569, "hasStops": 42, "isOOF": 527}},
    )
    plotted = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    def fake_barh(*args, **kwargs):
        plotted.append((kwargs["label"], list(kwargs["width"])))

    monkeypatch.setattr(Axes, "barh", fake_barh)

    mixcr.show_qc_plot(str(tmp_path), chart_type="chains", count_type="abs")

    trb_widths = dict(plotted)["TRB"]
    assert trb_widths == [14993.0, 54.0]


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


def test_show_qc_plot_coverage_labels_actual_values(monkeypatch):
    processing_table = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "reads_per_umi": 4.5,
                "overseq_threshold": 7,
            }
        ]
    )
    annotations = []
    marker_positions = []
    monkeypatch.setattr(mixcr.plt, "show", lambda: None)

    original_annotate = Axes.annotate
    original_vlines = Axes.vlines

    def capture_annotation(self, text, *args, **kwargs):
        annotations.append((text, kwargs))
        return original_annotate(self, text, *args, **kwargs)

    def capture_vlines(self, x, *args, **kwargs):
        marker_positions.append(x)
        return original_vlines(self, x, *args, **kwargs)

    monkeypatch.setattr(Axes, "annotate", capture_annotation)
    monkeypatch.setattr(Axes, "vlines", capture_vlines)

    mixcr.show_qc_plot(
        ".",
        chart_type="coverage",
        processing_table=processing_table,
    )

    annotations_by_text = {text: kwargs for text, kwargs in annotations}
    assert marker_positions == [6.0]
    assert annotations_by_text["4.5"]["xy"][0] == 4.5
    assert annotations_by_text["7"]["xy"][0] == 6.0
    assert annotations_by_text["7"]["color"] == "#d62728"


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
