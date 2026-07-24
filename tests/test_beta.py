import numpy as np
import pandas as pd
import pytest

from repseq import beta
from repseq.intersections import intersect_clones_in_samples_batch


def _full_table():
    return pd.DataFrame({
        "cdr3aa": ["A", "B", "C"],
        "sample1_count": [5, 3, 2],
        "sample2_count": [8, 0, 12],
        "sample1_freq": [0.5, 0.3, 0.2],
        "sample2_freq": [0.4, 0.0, 0.6],
        "sample1": ["s1"] * 3,
        "sample2": ["s2"] * 3,
        "pair": ["s1_vs_s2"] * 3,
    })


def _write_clonoset(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_metrics_from_table_returns_single_matrix_and_aliases():
    result = beta.metrics_from_table(_full_table(), metrics="jaccard distance")

    assert list(result.index) == ["s1", "s2"]
    assert result.loc["s1", "s2"] == pytest.approx(1 / 3)
    assert result.loc["s2", "s1"] == pytest.approx(1 / 3)
    assert result.loc["s1", "s1"] == 0


def test_metrics_from_table_calculates_requested_metrics():
    result = beta.metrics_from_table(
        _full_table(),
        metrics=[
            "number of intersecting clonotypes", "relative_diversity", "f1", "f2",
            "dice", "Szymkiewicz–Simpson coefficient", "bray-curtis", "L1",
            "total variation", "L2 distance", "morista horn", "hellinger",
            "kl-divergence", "full_table",
        ],
    )

    assert result["number_of_intersecting_clonotypes"].loc["s1", "s2"] == 2
    assert result["relative_diversity"].loc["s1", "s2"] == pytest.approx(1 / 3)
    assert result["f1"].loc["s1", "s2"] == pytest.approx(np.sqrt(0.7))
    assert result["f2"].loc["s1", "s2"] == pytest.approx(np.sqrt(0.2) + np.sqrt(0.12))
    assert result["dice"].loc["s1", "s2"] == pytest.approx(0.8)
    assert result["szymkiewicz_simpson"].loc["s1", "s2"] == 1
    assert result["bray_curtis"].loc["s1", "s2"] == pytest.approx(0.4)
    assert result["l1"].loc["s1", "s2"] == pytest.approx(0.8)
    assert result["total_variation"].loc["s1", "s2"] == pytest.approx(0.4)
    assert result["l2"].loc["s1", "s2"] == pytest.approx(np.sqrt(0.26))
    assert result["morisita_horn"].loc["s1", "s1"] == 1
    assert result["hellinger"].loc["s1", "s1"] == 0
    assert result["kl_divergence"].loc["s1", "s2"] != result["kl_divergence"].loc["s2", "s1"]
    assert result["full_table"] is not None


def test_metrics_none_returns_all_metrics_and_full_table():
    result = beta.metrics_from_table(_full_table())

    assert list(result) == [*beta.METRICS, "full_table"]


def test_intersection_supports_vj_and_vjlen(tmp_path):
    rows1 = [
        {"count": 6, "freq": 0.6, "cdr3nt": "AAA", "cdr3aa": "CASS", "v": "V1", "j": "J1"},
        {"count": 4, "freq": 0.4, "cdr3nt": "BBB", "cdr3aa": "CASST", "v": "V1", "j": "J1"},
    ]
    rows2 = [
        {"count": 5, "freq": 0.5, "cdr3nt": "CCC", "cdr3aa": "XXXX", "v": "V1", "j": "J1"},
        {"count": 5, "freq": 0.5, "cdr3nt": "DDD", "cdr3aa": "YYYYYY", "v": "V2", "j": "J2"},
    ]
    file1, file2 = tmp_path / "s1.tsv", tmp_path / "s2.tsv"
    _write_clonoset(file1, rows1)
    _write_clonoset(file2, rows2)
    clonosets = pd.DataFrame([
        {"sample_id": "s1", "filename": str(file1)},
        {"sample_id": "s2", "filename": str(file2)},
    ])

    with pytest.warns(DeprecationWarning, match="always returns counts"):
        vj = intersect_clones_in_samples_batch(
            clonosets, overlap_type="VJ", by_freq=True, cpu=1
        )
    vjlen = intersect_clones_in_samples_batch(clonosets, overlap_type="VJlen", cpu=1)

    assert list(vj.columns[:2]) == ["v", "j"]
    assert list(vj.columns[2:6]) == [
        "sample1_count", "sample2_count", "sample1_freq", "sample2_freq"
    ]
    shared_vj = vj[(vj["v"] == "V1") & (vj["j"] == "J1")].iloc[0]
    assert shared_vj["sample1_count"] == 10
    assert shared_vj["sample2_count"] == 5
    assert shared_vj["sample1_freq"] == 1
    assert shared_vj["sample2_freq"] == 0.5
    assert vj.groupby("pair")["sample1_freq"].sum().iloc[0] == pytest.approx(1)
    assert vj.groupby("pair")["sample2_freq"].sum().iloc[0] == pytest.approx(1)
    assert list(vjlen.columns[:3]) == ["v", "j", "len"]
    assert set(vjlen["len"]) == {4, 5, 6}


def test_metrics_passes_cpu_to_intersection(monkeypatch):
    captured = {}

    def fake_intersection(*args, **kwargs):
        captured.update(kwargs)
        return _full_table()

    monkeypatch.setattr(beta, "intersect_clones_in_samples_batch", fake_intersection)
    clonosets = pd.DataFrame([
        {"sample_id": "s1", "filename": "s1.tsv"},
        {"sample_id": "s2", "filename": "s2.tsv"},
    ])

    beta.metrics(clonosets, metrics="jaccard", cpu=3)

    assert captured["cpu"] == 3
