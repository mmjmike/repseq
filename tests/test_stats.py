import pandas as pd
import pytest

from repseq import stats
from repseq.clone_filter import Filter


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


def test_cdr3_length_distributions_long_zero_fill_with_chain(tmp_path):
    file1 = tmp_path / "sample1.tsv"
    file2 = tmp_path / "sample2.tsv"
    _write_clonoset(
        file1,
        [
            {"count": 2, "freq": 0.7, "cdr3nt": "TGTGCC", "cdr3aa": "CASS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.3, "cdr3nt": "TGTGCT", "cdr3aa": "CASSLG", "v": "TRBV2", "d": ".", "j": "TRBJ2"},
        ],
    )
    _write_clonoset(
        file2,
        [
            {"count": 5, "freq": 1.0, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
        ],
    )
    clonosets = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "filename": str(file1)},
            {"sample_id": "sample2", "chain": "TRB", "filename": str(file2)},
        ]
    )

    result = stats.cdr3_length_distributions(
        clonosets,
        cpu=1,
        verbose=False,
        table="long",
        zero_fill=True,
    )

    assert list(result.columns) == ["sample_id", "chain", "cdr3_length", "freq"]
    values = {
        (row.sample_id, row.cdr3_length): row.freq
        for row in result.itertuples(index=False)
    }
    assert values[("sample1", 3)] == 0
    assert values[("sample1", 4)] == 0.7
    assert values[("sample1", 6)] == 0.3
    assert values[("sample2", 3)] == 1.0
    assert values[("sample2", 4)] == 0
    assert values[("sample2", 6)] == 0


def test_cdr3_length_distributions_long_without_zero_fill(tmp_path):
    file1 = tmp_path / "sample1.tsv"
    file2 = tmp_path / "sample2.tsv"
    _write_clonoset(file1, [{"count": 2, "freq": 1.0, "cdr3nt": "TGTGCC", "cdr3aa": "CASS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}])
    _write_clonoset(file2, [{"count": 5, "freq": 1.0, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}])
    clonosets = pd.DataFrame(
        [
            {"sample_id": "sample1", "filename": str(file1)},
            {"sample_id": "sample2", "filename": str(file2)},
        ]
    )

    result = stats.cdr3_length_distributions(
        clonosets,
        cpu=1,
        verbose=False,
        zero_fill=False,
    )

    assert "chain" not in result.columns
    assert set(zip(result["sample_id"], result["cdr3_length"])) == {
        ("sample1", 4),
        ("sample2", 3),
    }


def test_cdr3_length_distributions_wide_is_zero_filled_and_counted_by_count(tmp_path):
    file1 = tmp_path / "sample1.tsv"
    file2 = tmp_path / "sample2.tsv"
    _write_clonoset(file1, [{"count": 2, "freq": 1.0, "cdr3nt": "TGTGCC", "cdr3aa": "CASS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}])
    _write_clonoset(file2, [{"count": 5, "freq": 1.0, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}])
    clonosets = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "filename": str(file1)},
            {"sample_id": "sample2", "chain": "TRB", "filename": str(file2)},
        ]
    )

    result = stats.cdr3_length_distributions(
        clonosets,
        cpu=1,
        verbose=False,
        count_by_freq=False,
        table="wide",
    )

    assert result.loc[result["sample_id"] == "sample1", 3].iloc[0] == 0
    assert result.loc[result["sample_id"] == "sample1", 4].iloc[0] == 2
    assert result.loc[result["sample_id"] == "sample2", 3].iloc[0] == 5
    assert result.loc[result["sample_id"] == "sample2", 4].iloc[0] == 0


def test_cdr3_length_distributions_nt_lengths(tmp_path):
    filename = tmp_path / "sample.tsv"
    _write_clonoset(
        filename,
        [
            {"count": 2, "freq": 1.0, "cdr3nt": "TGTGCC", "cdr3aa": "C", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
        ],
    )
    clonosets = pd.DataFrame([{"sample_id": "sample1", "filename": str(filename)}])

    result = stats.cdr3_length_distributions(
        clonosets,
        cpu=1,
        verbose=False,
        seq_type="nt",
    )

    assert result.loc[0, "cdr3_length"] == 6
    assert result.loc[0, "freq"] == 1.0


def test_calc_diversity_stats_passes_cpu_to_generic_calculation(monkeypatch):
    seen = {}

    def fake_generic_calculation(*args, **kwargs):
        seen["cpu"] = kwargs["cpu"]
        return pd.DataFrame([{"sample_id": "sample1", "diversity": 1}])

    monkeypatch.setattr(stats, "generic_calculation", fake_generic_calculation)

    result = stats.calc_diversity_stats(
        pd.DataFrame([{"sample_id": "sample1", "filename": "sample.tsv"}]),
        cpu=1,
        verbose=False,
    )

    assert seen["cpu"] == 1
    assert result.loc[0, "diversity"] == 1


def test_calc_cdr3_properties_passes_cpu_to_generic_calculation(monkeypatch):
    seen = {}

    def fake_generic_calculation(*args, **kwargs):
        seen["cpu"] = kwargs["cpu"]
        return pd.DataFrame([{"sample_id": "sample1", "mean_cdr3nt_len": 6}])

    monkeypatch.setattr(stats, "generic_calculation", fake_generic_calculation)

    result = stats.calc_cdr3_properties(
        pd.DataFrame([{"sample_id": "sample1", "filename": "sample.tsv"}]),
        cpu=1,
        verbose=False,
    )

    assert seen["cpu"] == 1
    assert result.loc[0, "mean_cdr3nt_len"] == 6


def test_calc_convergence_returns_cdr3_v_and_vj_convergence(tmp_path):
    filename = tmp_path / "sample.tsv"
    _write_clonoset(
        filename,
        [
            {"count": 1, "freq": 0.25, "cdr3nt": "AAA", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.25, "cdr3nt": "AAT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.25, "cdr3nt": "AAC", "cdr3aa": "CAS", "v": "TRBV2", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.25, "cdr3nt": "AAG", "cdr3aa": "CAT", "v": "TRBV2", "d": ".", "j": "TRBJ2"},
        ],
    )
    clonosets = pd.DataFrame([{"sample_id": "sample1", "filename": str(filename)}])

    result = stats.calc_convergence(clonosets, cpu=1, verbose=False)

    assert result.loc[0, "convergence"] == 2
    assert result.loc[0, "convergence_v"] == pytest.approx(4 / 3)
    assert result.loc[0, "convergence_vj"] == pytest.approx(4 / 3)


def test_calc_convergence_defaults_to_three_iterations_and_prints_settings(monkeypatch, capsys):
    seen = {}

    def fake_generic_calculation(*args, **kwargs):
        seen["iterations"] = kwargs["iterations"]
        return pd.DataFrame([{"sample_id": "sample1", "convergence": 1}])

    monkeypatch.setattr(stats, "generic_calculation", fake_generic_calculation)

    stats.calc_convergence(
        pd.DataFrame([{"sample_id": "sample1", "filename": "sample.tsv"}]),
        cl_filter=Filter(downsample=100),
        seed=123,
    )

    captured = capsys.readouterr().out
    assert seen["iterations"] == 3
    assert "equal downsampling for all samples" not in captured
    assert "seed=123, downsample=100, top=None, iterations=3" in captured


def test_calc_diversity_stats_prints_normalization_settings(monkeypatch, capsys):
    def fake_generic_calculation(*args, **kwargs):
        return pd.DataFrame([{"sample_id": "sample1", "diversity": 1}])

    monkeypatch.setattr(stats, "generic_calculation", fake_generic_calculation)

    stats.calc_diversity_stats(
        pd.DataFrame([{"sample_id": "sample1", "filename": "sample.tsv"}]),
        cl_filter=Filter(top=50),
        seed=7,
        iterations=5,
    )

    captured = capsys.readouterr().out
    assert "equal downsampling for all samples" in captured
    assert "seed=7, downsample=None, top=50, iterations=5" in captured


def test_calc_diversity_stats_prints_recommendation_when_downsample_is_combined_with_extra_filters(monkeypatch, capsys):
    def fake_generic_calculation(*args, **kwargs):
        return pd.DataFrame([{"sample_id": "sample1", "diversity": 1}])

    monkeypatch.setattr(stats, "generic_calculation", fake_generic_calculation)

    stats.calc_diversity_stats(
        pd.DataFrame([{"sample_id": "sample1", "filename": "sample.tsv"}]),
        cl_filter=Filter(downsample=100, count_threshold=2),
        seed=7,
    )

    captured = capsys.readouterr().out
    assert "equal downsampling for all samples" in captured
    assert "seed=7, downsample=100, top=None, iterations=3" in captured


def test_generic_calculation_warns_and_keeps_small_samples_by_default(tmp_path, capsys):
    small_file = tmp_path / "small.tsv"
    large_file = tmp_path / "large.tsv"
    _write_clonoset(
        small_file,
        [{"count": 1, "freq": 1.0, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}],
    )
    _write_clonoset(
        large_file,
        [
            {"count": 1, "freq": 0.5, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.5, "cdr3nt": "TGTGCC", "cdr3aa": "CASS", "v": "TRBV2", "d": ".", "j": "TRBJ2"},
        ],
    )
    clonosets = pd.DataFrame(
        [
            {"sample_id": "small", "filename": str(small_file)},
            {"sample_id": "large", "filename": str(large_file)},
        ]
    )

    result = stats.calc_diversity_stats(
        clonosets,
        cl_filter=Filter(top=2),
        cpu=1,
        verbose=False,
    )

    captured = capsys.readouterr().out
    assert "WARNING! top=2 exceeds available clonotypes for samples: small" in captured
    assert "kept and calculated without top filtering" in captured
    assert result.loc[result["sample_id"] == "small", "diversity"].iloc[0] == 1


def test_generic_calculation_warns_and_drops_small_samples_when_requested(tmp_path, capsys):
    small_file = tmp_path / "small.tsv"
    large_file = tmp_path / "large.tsv"
    _write_clonoset(
        small_file,
        [{"count": 1, "freq": 1.0, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"}],
    )
    _write_clonoset(
        large_file,
        [
            {"count": 1, "freq": 0.5, "cdr3nt": "TGT", "cdr3aa": "CAS", "v": "TRBV1", "d": ".", "j": "TRBJ1"},
            {"count": 1, "freq": 0.5, "cdr3nt": "TGTGCC", "cdr3aa": "CASS", "v": "TRBV2", "d": ".", "j": "TRBJ2"},
        ],
    )
    clonosets = pd.DataFrame(
        [
            {"sample_id": "small", "filename": str(small_file)},
            {"sample_id": "large", "filename": str(large_file)},
        ]
    )

    result = stats.calc_diversity_stats(
        clonosets,
        cl_filter=Filter(top=2),
        drop_small_samples=True,
        cpu=1,
        verbose=False,
    )

    captured = capsys.readouterr().out
    assert "WARNING! top=2 exceeds available clonotypes for samples: small" in captured
    assert "excluded from further calculations" in captured
    assert pd.isna(result.loc[result["sample_id"] == "small", "diversity"].iloc[0])
