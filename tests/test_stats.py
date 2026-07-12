import pandas as pd

from repseq import stats


def _write_clonoset(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_calc_clonoset_stats_handles_clonoset_without_umi_columns(tmp_path):
    filename = tmp_path / "sample.tsv"
    _write_clonoset(
        filename,
        [
            {
                "count": 2,
                "freq": 2 / 3,
                "cdr3nt": "TGTGCC",
                "cdr3aa": "CASSLG",
                "v": "TRBV1",
                "d": ".",
                "j": "TRBJ1",
            },
            {
                "count": 1,
                "freq": 1 / 3,
                "cdr3nt": "TGTGCT",
                "cdr3aa": "CAS*LG",
                "v": "TRBV2",
                "d": ".",
                "j": "TRBJ2",
            },
        ],
    )
    clonosets = pd.DataFrame(
        [{"sample_id": "sample1", "chain": "TRB", "filename": str(filename)}]
    )

    result = stats.calc_clonoset_stats(clonosets, verbose=False, cpu=1)

    assert result.loc[0, "clones"] == 2
    assert result.loc[0, "reads"] == 3
    assert result.loc[0, "clones_func"] == 1
    assert result.loc[0, "reads_func"] == 2
    assert pd.isna(result.loc[0, "umi"])


def test_generic_calculation_passes_cpu_to_parallel_runner(monkeypatch, tmp_path):
    filename = tmp_path / "sample.tsv"
    _write_clonoset(
        filename,
        [
            {
                "count": 1,
                "freq": 1.0,
                "cdr3nt": "TGTGCC",
                "cdr3aa": "CASSLG",
                "v": "TRBV1",
                "d": ".",
                "j": "TRBJ1",
            }
        ],
    )
    clonosets = pd.DataFrame(
        [{"sample_id": "sample1", "chain": "TRB", "filename": str(filename)}]
    )
    seen_cpu = []

    def fake_run_parallel(function, tasks, program_name, object_name="tasks", verbose=True, cpu=None):
        seen_cpu.append(cpu)
        return [function(task) for task in tasks]

    monkeypatch.setattr(stats, "run_parallel_calculation", fake_run_parallel)

    result = stats.generic_calculation(
        clonosets,
        lambda clonoset: {"clones": len(clonoset)},
        verbose=False,
        cpu=1,
    )

    assert seen_cpu == [1]
    assert result.loc[0, "clones"] == 1
