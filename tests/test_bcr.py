import matplotlib
import pandas as pd
import pytest
from pathlib import Path

matplotlib.use("Agg")

from repseq import bcr
from repseq import plot as rsplot


REGIONS = ["aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"]


def _write_tree_inputs(tmp_path):
    rows = [
        {
            "treeId": 1, "nodeId": 1, "isObserved": True,
            "fileName": "mix.sample.one.clns", "cloneId": 1,
            "readCount": 5, "isotype": "IgM", "bestVHit": "IGHV1",
            "bestJHit": "IGHJ1", "mutationRate": 0.1,
            "distanceFromGermline": 2,
        },
        {
            "treeId": 1, "nodeId": 2, "isObserved": True,
            "fileName": "mix.sample.one.clns", "cloneId": 2,
            "readCount": 2, "isotype": "IgG", "bestVHit": "IGHV2",
            "bestJHit": "IGHJ2", "mutationRate": 0.2,
            "distanceFromGermline": 4,
        },
        {
            "treeId": 1, "nodeId": 3, "isObserved": False,
            "fileName": pd.NA, "cloneId": 999, "readCount": 0,
            "isotype": pd.NA, "bestVHit": "IGHV1", "bestJHit": "IGHJ1",
            "mutationRate": 0.3, "distanceFromGermline": 5,
        },
        {
            "treeId": 2, "nodeId": 4, "isObserved": True,
            "fileName": "sample.two.clns", "cloneId": 4,
            "readCount": 20, "uniqueMoleculeCount": 3, "isotype": "IgD",
            "bestVHit": "IGHV3", "bestJHit": "IGHJ3", "mutationRate": 0.05,
            "distanceFromGermline": 1,
        },
    ]
    sequences = {
        1: ["CAR", "WAA", "GG", "TTT", "AAA", "WG"],
        2: ["CAS", "WAA", "GA", "TTA", "ABA", "WG"],
        3: ["CAT", "WAA", "GA", "TTA", "ACA", "WG"],
        4: ["CCC", "FFF", "GG", "HHH", "ZZZ", "WW"],
    }
    for row in rows:
        for column, sequence in zip(REGIONS, sequences[row["nodeId"]]):
            row[column] = sequence
    trees_dir = tmp_path / "trees"
    trees_dir.mkdir()
    trees_filename = trees_dir / "donor_trees.tsv"
    pd.DataFrame(rows).to_csv(trees_filename, sep="\t", index=False)
    repertoires_dir = tmp_path / "repertoires"
    repertoires_dir.mkdir()
    (repertoires_dir / "mix.sample.one.clns").touch()
    pd.DataFrame(
        {"cloneId": [1, 2], "readCount": [5, 2], "readFraction": [0.7, 0.3],
         "uniqueMoleculeCount": [10, 2]}
    ).to_csv(repertoires_dir / "mix.sample.one.clones_IGH.tsv", sep="\t", index=False)
    newick_dir = tmp_path / "newick"
    newick_dir.mkdir()
    (newick_dir / "1.tree").write_text("(1:0.1,2:0.2)3;")
    (newick_dir / "2.tree").write_text("4;")
    return trees_filename, newick_dir


def test_tree_analyzer_properties_are_enriched_sorted_and_cached(tmp_path, capsys):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    assert analyzer.read_trees_table(trees_filename) is None
    assert "Read 2 trees from 2 samples: 4 nodes, 3 observed nodes" in capsys.readouterr().out
    assert analyzer.trees_df.loc[
        analyzer.trees_df["nodeId"] == 1, "uniqueMoleculeCount"
    ].iloc[0] == 10
    assert analyzer.read_trees_newick(newick_dir) is None
    assert "Read 2 Newick tree filenames" in capsys.readouterr().out
    analyzer.read_metadata(pd.DataFrame({"sample_id": ["sample.one", "two"], "timepoint": [1, 2]}))

    properties = analyzer.trees_properties

    assert analyzer.trees_properties is properties
    assert properties["treeId"].tolist() == [1, 2]
    first = properties.iloc[0]
    assert first["nodes"] == 3
    assert first["nodes_obs"] == 2
    assert first["reads"] == 7
    assert first["umi"] == 12
    assert first["isotypes"] == ["IgG", "IgM"]
    assert first["v"] == "IGHV1"
    assert first["j"] == "IGHJ1"
    assert first["consensus_CDR1"] == "CAR"
    assert first["consensus_CDR3"] == "AAA"
    assert first["mean_mutation_rate"] == pytest.approx(0.15)
    assert bool(first["isotype_switched"])
    assert first["max_distance_from_germline"] == 5
    assert analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 1, "sample_id"].iloc[0] == "sample.one"
    assert analyzer.trees_df.loc[
        analyzer.trees_df["treeId"] == 1, "newick_filename"
    ].unique().tolist() == [str((newick_dir / "1.tree").resolve())]


def test_read_trees_newick_attaches_paths_without_reading(tmp_path, monkeypatch, capsys):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)

    def fail_if_read(_path, *args, **kwargs):
        raise AssertionError("read_trees_newick must not read tree contents")

    with monkeypatch.context() as patcher:
        patcher.setattr(Path, "read_text", fail_if_read)
        result = analyzer.read_trees_newick(newick_dir)

    assert result is None
    assert analyzer.newick_trees["1"] == str((newick_dir / "1.tree").resolve())
    assert "Read 2 Newick tree filenames" in capsys.readouterr().out
    assert analyzer.trees_df.loc[
        analyzer.trees_df["treeId"] == 1, "newick_filename"
    ].nunique() == 1


def test_umi_enrichment_discovers_and_pools_clonosets_once(tmp_path, monkeypatch):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    find_calls = []
    pool_calls = []
    original_find = bcr.clonosets.find_all_mixcr_clonosets
    original_pool = bcr.clonosets.pool_clonotypes_from_clonosets_df

    def tracked_find(folders, *args, **kwargs):
        find_calls.append(folders)
        return original_find(folders, *args, **kwargs)

    def tracked_pool(clonosets_df, cl_filter=None):
        pool_calls.append((clonosets_df.copy(), cl_filter))
        return original_pool(clonosets_df, cl_filter=cl_filter)

    monkeypatch.setattr(bcr.clonosets, "find_all_mixcr_clonosets", tracked_find)
    monkeypatch.setattr(bcr.clonosets, "pool_clonotypes_from_clonosets_df", tracked_pool)

    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.trees_properties

    assert len(find_calls) == 1
    assert len(pool_calls) == 1
    assert pool_calls[0][0]["sample_id"].tolist() == ["sample.one"]
    assert pool_calls[0][1].convert is False
    assert analyzer.trees_df.loc[
        analyzer.trees_df["nodeId"] == 2, "uniqueMoleculeCount"
    ].iloc[0] == 2


def test_tree_analyzer_draw_tree_wrapper(tmp_path, monkeypatch):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    metadata = pd.DataFrame({"sample_id": ["one"]})
    analyzer.read_metadata(metadata)
    captured = {}

    def fake_draw_tree(trees_df, tree_id, **kwargs):
        captured.update(trees_df=trees_df, tree_id=tree_id, **kwargs)
        return "axis"

    monkeypatch.setattr(bcr.rsplot, "draw_tree", fake_draw_tree)
    assert analyzer.draw_tree(1) == "axis"
    assert captured == {
        "trees_df": analyzer.trees_df,
        "tree_id": 1,
        "metadata": analyzer.metadata,
        "group": "isotype",
        "label": "timepoint",
    }


def test_draw_tree_returns_axis(tmp_path):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.read_trees_newick(newick_dir)
    analyzer.read_metadata(pd.DataFrame({"sample_id": ["sample.one"], "timepoint": ["day 1"]}))
    analyzer.trees_properties

    axis = rsplot.draw_tree(
        analyzer.trees_df, 1, metadata=analyzer.metadata, group="isotype", label="timepoint"
    )

    assert axis.get_title() == "Tree 1"
    assert len(axis.collections) == 3
