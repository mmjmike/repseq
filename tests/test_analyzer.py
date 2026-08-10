import pandas as pd
import pytest

import repseq.diff_enrichment as rsde
from repseq import intersections


def _single_chain_samples():
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "filename": ["a.tsv", "b.tsv"],
            "group": ["A", "B"],
        }
    )


def _paired_samples(chains=("TRA", "TRB")):
    return pd.DataFrame(
        {
            "sample": ["p1", "p2", "p1", "p2"],
            "sample_id": ["a1", "a2", "b1", "b2"],
            "filename": ["a1.tsv", "a2.tsv", "b1.tsv", "b2.tsv"],
            "group": ["A", "B", "A", "B"],
            "chain": [chains[0], chains[0], chains[1], chains[1]],
        }
    )


def _count_table(sample_ids):
    return pd.DataFrame(
        {
            "clonotype": ["CASS|TRBV1"],
            "cdr3aa": ["CASS"],
            "v": ["TRBV1"],
            **{sample_id: [3] for sample_id in sample_ids},
        }
    )


def test_analyzer_defaults_and_parameter_updates():
    analyzer = rsde.Analyzer(parameters={"min_samples": 5}, verbose=False)

    parameters = analyzer.get_parameters()

    assert parameters["min_samples"] == 5
    assert parameters["overlap_type"] == "aaV"
    assert parameters["presence_threshold"] == 2
    assert parameters["pairing_method"] == "jsd"
    assert parameters["cl_filter"].functionality == "f"
    assert parameters["cl_filter"].by_umi is True


def test_read_samples_normalizes_supported_chain_aliases(capsys):
    analyzer = rsde.Analyzer(verbose=True)
    analyzer.read_samples(_paired_samples(("TRAD", "TRB")))

    assert analyzer.chains == ["TRA", "TRB"]
    assert analyzer.samples_df["chain"].tolist() == ["TRA", "TRA", "TRB", "TRB"]
    assert "Renaming chain TRAD to TRA" in capsys.readouterr().out


def test_read_samples_rejects_unsupported_chain_combinations():
    analyzer = rsde.Analyzer(verbose=False)

    with pytest.raises(ValueError, match="Unsupported chain combination"):
        analyzer.read_samples(_paired_samples(("TRA", "TRD")))


def test_read_samples_rejects_igh_with_separate_light_chains():
    samples = pd.DataFrame(
        {
            "sample": ["p1", "p1", "p1"],
            "sample_id": ["h1", "k1", "l1"],
            "filename": ["h1.tsv", "k1.tsv", "l1.tsv"],
            "group": ["A", "A", "A"],
            "chain": ["IGH", "IGK", "IGL"],
        }
    )

    with pytest.raises(ValueError, match="merge IGK and IGL manually"):
        rsde.Analyzer(verbose=False).read_samples(samples)


def test_count_table_forwards_runtime_arguments_and_ignores_mismatches(monkeypatch):
    calls = []

    def fake_count_table(samples_df, **kwargs):
        calls.append((samples_df.copy(), kwargs))
        return _count_table(samples_df["sample_id"])

    monkeypatch.setattr(intersections, "count_table", fake_count_table)
    analyzer = rsde.Analyzer(
        samples_df=_single_chain_samples(), cpu=1, verbose=False, mismatches=7
    )

    result = analyzer.run_count_table()
    repeated = analyzer.run_count_table()

    assert result is repeated
    assert len(calls) == 1
    assert calls[0][1]["mismatches"] == 0
    assert calls[0][1]["cpu"] == 1
    assert calls[0][1]["verbose"] is False


def test_chain_specific_parameters_and_callable_results(monkeypatch):
    calls = []

    def fake_count_table(samples_df, **kwargs):
        calls.append((samples_df["chain"].iloc[0], kwargs["overlap_type"]))
        return _count_table(samples_df["sample_id"])

    monkeypatch.setattr(intersections, "count_table", fake_count_table)
    analyzer = rsde.Analyzer(samples_df=_paired_samples(), verbose=False)
    analyzer.run_count_table(overlap_type="aaVJ")
    analyzer.select_chain("TRB")
    analyzer.run_count_table()

    assert calls == [("TRA", "aaVJ"), ("TRB", "aaV")]
    assert analyzer.get_parameters()["overlap_type_TRA"] == "aaVJ"
    assert analyzer.count_table("TRA")["a1"].tolist() == [3]


def test_parameter_updates_invalidate_only_downstream_results(monkeypatch):
    monkeypatch.setattr(
        intersections,
        "count_table",
        lambda samples_df, **kwargs: _count_table(samples_df["sample_id"]),
    )
    analyzer = rsde.Analyzer(samples_df=_single_chain_samples(), verbose=False)
    analyzer.run_count_table()
    analyzer.run_prefilter(min_samples=1, min_total_count=1)
    count_table = analyzer.count_table

    analyzer.update_parameters(min_count=3)

    assert analyzer.count_table is count_table
    assert analyzer.prefiltered is None


def test_run_executes_both_branches_and_pairs(monkeypatch):
    monkeypatch.setattr(
        intersections,
        "count_table",
        lambda samples_df, **kwargs: _count_table(samples_df["sample_id"]),
    )

    def fake_statistics(count_table, samples_metadata, **kwargs):
        result = count_table.copy()
        insert_at = 3
        for name, value in reversed(
            [
                ("enriched_in", samples_metadata["group"].iloc[0]),
                ("method", kwargs["method"]),
                ("mean_group_count", 3.0),
                ("log2FC", 2.0),
                ("p_val", 0.01),
                ("p_adj", 0.02),
            ]
        ):
            result.insert(insert_at, name, value)
        return result

    monkeypatch.setattr(rsde, "calc_statistics", fake_statistics)
    analyzer = rsde.Analyzer(
        samples_df=_paired_samples(),
        verbose=False,
        min_samples=1,
        min_total_count=1,
        max_p_adj=0.05,
    )

    analyzer.run()

    assert analyzer.postfiltered("TRA")["postfilter_pass"].all()
    assert analyzer.postfiltered("TRB")["postfilter_pass"].all()
    assert analyzer.pairing_matrix.shape == (1, 1)
    assert "Pairing matrix: 1 x 1" in repr(analyzer)


def test_pairing_supports_sample_ids_reused_between_chains(monkeypatch):
    samples = _paired_samples()
    samples.loc[samples["chain"].eq("TRB"), "sample_id"] = ["a1", "a2"]
    monkeypatch.setattr(
        intersections,
        "count_table",
        lambda samples_df, **kwargs: _count_table(samples_df["sample_id"]),
    )

    def fake_statistics(count_table, samples_metadata, **kwargs):
        result = count_table.copy()
        for name, value in (
            ("enriched_in", "A"),
            ("method", "mann_whitney"),
            ("mean_group_count", 3.0),
            ("log2FC", 2.0),
            ("p_val", 0.01),
            ("p_adj", 0.02),
        ):
            result.insert(3, name, value)
        return result

    monkeypatch.setattr(rsde, "calc_statistics", fake_statistics)
    analyzer = rsde.Analyzer(
        samples_df=samples,
        verbose=False,
        min_samples=1,
        min_total_count=1,
        max_p_adj=0.05,
    )

    analyzer.run()

    assert analyzer.pairing_matrix.shape == (1, 1)
    assert analyzer.samples_df["sample_id"].tolist() == ["a1", "a2", "a1", "a2"]


def test_analyzer_has_only_correctly_spelled_prefiltered_property():
    assert hasattr(rsde.Analyzer, "prefiltered")
    assert not hasattr(rsde.Analyzer, "prefilered")


def test_plot_volcano_prefers_postfiltered_and_can_force_statistics(monkeypatch):
    from repseq import plot as rsplot

    analyzer = rsde.Analyzer(samples_df=_single_chain_samples(), verbose=False)
    statistics = pd.DataFrame({"feature": ["statistics"]})
    postfiltered = pd.DataFrame({"feature": ["postfiltered"]})
    analyzer._store("statistics_df", statistics, {}, "XCR")
    analyzer._store("postfiltered", postfiltered, {}, "XCR")
    plotted = []
    monkeypatch.setattr(
        rsplot,
        "de_volcano",
        lambda table, **kwargs: plotted.append((table, kwargs)) or table,
    )

    assert analyzer.plot_volcano(alpha=0.5).iloc[0, 0] == "postfiltered"
    assert analyzer.plot_volcano(postfiltered=False).iloc[0, 0] == "statistics"
    assert plotted[0][1] == {"alpha": 0.5}


def test_plot_volcano_falls_back_to_statistics(monkeypatch):
    from repseq import plot as rsplot

    analyzer = rsde.Analyzer(samples_df=_single_chain_samples(), verbose=False)
    statistics = pd.DataFrame({"feature": ["statistics"]})
    analyzer._store("statistics_df", statistics, {}, "XCR")
    monkeypatch.setattr(rsplot, "de_volcano", lambda table, **kwargs: table)

    assert analyzer.plot_volcano().iloc[0, 0] == "statistics"
