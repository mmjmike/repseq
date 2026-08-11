import math

import pytest

from repseq import slurm


SINFO_OUTPUT = """short*|02:00:00|node[01-02]|hpc
medium|16:00:00|node[03-04]|hpc,gpu
long|6-00:00:00|node05|bigmem
infinite|infinite|node06|(null)
"""


def test_slurm_time_to_hours_parses_sinfo_limits():
    assert slurm._slurm_time_to_hours("02:30:00") == 2.5
    assert slurm._slurm_time_to_hours("1-02:30:00") == 26.5
    assert slurm._slurm_time_to_hours("30:00") == 0.5
    assert slurm._slurm_time_to_hours("90") == 1.5
    assert math.isinf(slurm._slurm_time_to_hours("infinite"))


def test_update_partition_limits_uses_sinfo(monkeypatch):
    monkeypatch.setattr(slurm, "_run_sinfo", lambda: SINFO_OUTPUT)
    monkeypatch.setattr(slurm, "PARTITION_LIMIT_HOURS", {})

    limits = slurm.update_partition_limits()

    assert limits == {
        "short": 2,
        "medium": 16,
        "long": 144,
        "infinite": math.inf,
    }
    assert slurm.PARTITION_LIMIT_HOURS == limits


@pytest.mark.parametrize(
    ("hours", "expected_partition"),
    [
        (1.5, "short"),
        (2, "short"),
        (2.1, "medium"),
        (16, "medium"),
        (16.1, "long"),
        (144, "long"),
        (145, "infinite"),
    ],
)
def test_partition_by_time_uses_discovered_limits(hours, expected_partition):
    limits = {
        "short": "02:00:00",
        "medium": "16:00:00",
        "long": "6-00:00:00",
        "infinite": "infinite",
    }

    assert slurm.partition_by_time(hours, limits) == expected_partition


def test_partition_by_time_rejects_requests_above_all_limits():
    limits = {"short": 2, "medium": 16, "long": 144}

    with pytest.raises(ValueError, match="exceeds available"):
        slurm.partition_by_time(145, limits)


def test_run_slurm_command_writes_selected_partition_and_constraint(tmp_path, monkeypatch):
    monkeypatch.setattr(slurm, "ALDAN_TEMP_SLURM_DIR", str(tmp_path))
    monkeypatch.setattr(
        slurm,
        "update_partition_limits",
        lambda: {"short": 2, "medium": 16, "long": 144, "infinite": math.inf},
    )

    class FakeProcess:
        returncode = 0

        def communicate(self):
            return b"Submitted batch job 42\n", b""

    monkeypatch.setattr(slurm, "Popen", lambda *args, **kwargs: FakeProcess())

    stdout, stderr = slurm.run_slurm_command_from_jupyter(
        "echo hello",
        jobname="constraint_test",
        time_estimate=4.25,
        constraint="hpc",
        verbose=False,
    )

    script = next(tmp_path.glob("*_constraint_test.sh")).read_text()
    assert stdout == b"Submitted batch job 42\n"
    assert stderr == b""
    assert "#SBATCH --partition=medium" in script
    assert "#SBATCH --time=4:15:00" in script
    assert "#SBATCH --constraint=hpc" in script


def test_run_slurm_command_omits_constraint_by_default(tmp_path, monkeypatch):
    monkeypatch.setattr(slurm, "ALDAN_TEMP_SLURM_DIR", str(tmp_path))
    monkeypatch.setattr(slurm, "update_partition_limits", lambda: {"short": 2})

    class FakeProcess:
        returncode = 0

        def communicate(self):
            return b"Submitted batch job 42\n", b""

    monkeypatch.setattr(slurm, "Popen", lambda *args, **kwargs: FakeProcess())

    slurm.run_slurm_command_from_jupyter(
        "echo hello",
        jobname="no_constraint_test",
        verbose=False,
    )

    script = next(tmp_path.glob("*_no_constraint_test.sh")).read_text()
    assert "#SBATCH --constraint" not in script


def test_info_prints_partitions_nodes_and_constraints(monkeypatch, capsys):
    monkeypatch.setattr(slurm, "_run_sinfo", lambda: SINFO_OUTPUT)

    slurm.info()

    output = capsys.readouterr().out
    assert "SLURM partitions and nodes" in output
    assert "short" in output
    assert "02:00:00" in output
    assert "node[01-02]" in output
    assert "Possible constraints: bigmem, gpu, hpc" in output


def test_info_reports_when_slurm_is_not_installed(monkeypatch, capsys):
    monkeypatch.setattr(slurm.shutil, "which", lambda command: None)

    slurm.info()

    assert slurm.SLURM_NOT_AVAILABLE_MESSAGE in capsys.readouterr().out


def test_run_slurm_command_reports_when_slurm_is_not_installed(monkeypatch):
    monkeypatch.setattr(slurm.shutil, "which", lambda command: None)

    with pytest.raises(RuntimeError, match="SLURM is not installed"):
        slurm.run_slurm_command_from_jupyter("echo hello", jobname="missing_slurm")
