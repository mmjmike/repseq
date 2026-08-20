import matplotlib
import numpy as np
import pandas as pd
import pytest
import sys
import types
from pathlib import Path
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgba
from io import StringIO
from Bio import Phylo

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
            "DistanceFromGermline": 2,
        },
        {
            "treeId": 1, "nodeId": 2, "isObserved": True,
            "fileName": "mix.sample.one.clns", "cloneId": 2,
            "readCount": 2, "isotype": "IgG", "bestVHit": "IGHV2",
            "bestJHit": "IGHJ2", "mutationRate": 0.2,
            "DistanceFromGermline": 4,
        },
        {
            "treeId": 1, "nodeId": 3, "isObserved": False,
            "fileName": pd.NA, "cloneId": 999, "readCount": 0,
            "isotype": pd.NA, "bestVHit": "IGHV1", "bestJHit": "IGHJ1",
            "mutationRate": 0.3, "DistanceFromGermline": 5,
        },
        {
            "treeId": 2, "nodeId": 4, "isObserved": True,
            "fileName": "sample.two.clns", "cloneId": 4,
            "readCount": 20, "uniqueMoleculeCount": 3, "isotype": "IgD",
            "uniqueMoleculeFraction": 0.3,
            "bestVHit": "IGHV3", "bestJHit": "IGHJ3", "mutationRate": 0.05,
            "DistanceFromGermline": 1,
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
         "uniqueMoleculeCount": [10, 2], "uniqueMoleculeFraction": [0.1, 0.02],
         "customClonosetColumn": ["first", "second"]}
    ).to_csv(repertoires_dir / "mix.sample.one.clones_IGH.tsv", sep="\t", index=False)
    newick_dir = tmp_path / "newick"
    newick_dir.mkdir()
    (newick_dir / "1.tree").write_text("(1:0.1,2:0.2)3;")
    (newick_dir / "2.tree").write_text("4;")
    return trees_filename, newick_dir


def _mutation_ref_points():
    points = [""] * 22
    positions = {
        5: 2,
        6: 4,
        7: 6,
        8: 8,
        9: 10,
        18: 13,
        19: 15,
    }
    for index, position in positions.items():
        points[index] = str(position)
    return ":".join(points)


def test_tree_analyzer_properties_are_enriched_sorted_and_cached(tmp_path, capsys):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    assert analyzer.read_trees_table(trees_filename) is None
    assert "Read 2 trees from 2 samples: 4 nodes, 3 observed nodes" in capsys.readouterr().out
    assert analyzer.trees_df.loc[
        analyzer.trees_df["nodeId"] == 1, "uniqueMoleculeCount"
    ].iloc[0] == 10
    assert analyzer.trees_df.loc[
        analyzer.trees_df["nodeId"] == 1, "uniqueMoleculeFraction"
    ].iloc[0] == pytest.approx(0.1)
    assert analyzer.read_trees_newick(newick_dir) is None
    assert "Read 2 Newick tree filenames" in capsys.readouterr().out
    assert analyzer.read_metadata(
        pd.DataFrame({"sample_id": ["sample.one", "two"], "timepoint": [1, 2]})
    ) is None
    assert "Metadata updated for 2 samples and 3 observed nodes" in capsys.readouterr().out
    assert analyzer.trees_df.loc[
        analyzer.trees_df["nodeId"] == 1, "timepoint"
    ].iloc[0] == 1
    assert pd.isna(
        analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 3, "timepoint"].iloc[0]
    )

    properties = analyzer.trees_properties
    output = capsys.readouterr().out

    assert analyzer.trees_properties is properties
    assert "Calculating properties for 2 trees (CDR3 consensus only)" in output
    assert "Finished calculating tree properties" in output
    assert capsys.readouterr().out == ""
    assert properties["treeId"].tolist() == [1, 2]
    first = properties.iloc[0]
    assert first["nodes"] == 3
    assert first["nodes_obs"] == 2
    assert first["reads"] == 7
    assert first["umi"] == 12
    assert first["isotypes"] == ["IgG", "IgM"]
    assert first["v"] == "IGHV1"
    assert first["j"] == "IGHJ1"
    assert "consensus_CDR1" not in properties.columns
    assert first["consensus_CDR3"] == "AAA"
    assert first["mean_mutation_rate"] == pytest.approx(0.15)
    assert bool(first["isotype_switched"])
    assert first["max_distance_from_germline"] == 5
    assert analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 1, "sample_id"].iloc[0] == "sample.one"
    assert analyzer.trees_df.loc[
        analyzer.trees_df["treeId"] == 1, "newick_filename"
    ].unique().tolist() == [str((newick_dir / "1.tree").resolve())]

    full_properties = properties(all_consensuses=True)
    output = capsys.readouterr().out
    assert full_properties is properties
    assert "Calculating additional consensus sequences for 2 trees" in output
    assert "Finished calculating additional consensus sequences" in output
    assert full_properties.loc[0, "consensus_CDR1"] == "CAR"
    assert full_properties.loc[0, "consensus_FR2"] == "WAA"
    assert full_properties.loc[0, "consensus_CDR2"] == "GG"
    assert full_properties.loc[0, "consensus_FR3"] == "TTT"
    assert full_properties.loc[0, "consensus_FR4"] == "WG"
    assert analyzer.trees_properties is full_properties
    assert capsys.readouterr().out == ""


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
        "metadata": None,
        "group": "isotype",
        "label": "timepoint",
    }


def test_read_metadata_before_trees_table_does_nothing(capsys):
    analyzer = bcr.TreeAnalyzer()

    assert analyzer.read_metadata(pd.DataFrame({"sample_id": ["sample.one"]})) is None

    assert analyzer.metadata is None
    assert analyzer.trees_df is None
    assert "Trees table has not been read. Metadata was not loaded" in capsys.readouterr().out


def test_to_count_table_uses_metadata_order_and_umi_counts(tmp_path):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.read_metadata(
        pd.DataFrame(
            {
                "sample_id": ["two", "sample.one"],
                "timepoint": [2, 1],
            }
        )
    )

    count_table = analyzer.to_count_table()

    assert count_table.columns.tolist() == [
        "treeId",
        "v",
        "j",
        "consensus_cdr3",
        "two",
        "sample.one",
    ]
    assert count_table.loc[count_table["treeId"] == 1, "sample.one"].iloc[0] == 12
    assert count_table.loc[count_table["treeId"] == 1, "two"].iloc[0] == 0
    assert count_table.loc[count_table["treeId"] == 2, "two"].iloc[0] == 3


def test_to_count_table_sorts_samples_and_falls_back_to_reads(tmp_path):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.trees_df["uniqueMoleculeCount"] = pd.NA

    count_table = analyzer.to_count_table()

    assert count_table.columns.tolist() == [
        "treeId",
        "v",
        "j",
        "consensus_cdr3",
        "sample.one",
        "two",
    ]
    assert count_table.loc[count_table["treeId"] == 1, "sample.one"].iloc[0] == 7
    assert count_table.loc[count_table["treeId"] == 2, "two"].iloc[0] == 20


def test_to_count_table_after_full_properties_cache_uses_numeric_numpy_values(
    tmp_path, monkeypatch
):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.trees_properties(all_consensuses=True)
    original_isclose = np.isclose

    def checked_isclose(left, right, *args, **kwargs):
        assert np.asarray(left).dtype.kind in "fiu"
        assert np.asarray(right).dtype.kind in "fiu"
        return original_isclose(left, right, *args, **kwargs)

    monkeypatch.setattr(bcr.np, "isclose", checked_isclose)

    count_table = analyzer.to_count_table()

    assert count_table.loc[count_table["treeId"] == 1, "sample.one"].iloc[0] == 12
    assert count_table.loc[count_table["treeId"] == 2, "two"].iloc[0] == 3


@pytest.mark.parametrize("requested_tree_id", [1, "1", 1.0])
def test_get_logo_for_tree_passes_tree_cdr3_sequences(
    tmp_path, monkeypatch, requested_tree_id
):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    captured = {}
    logo_module = types.ModuleType("repseq.logo")

    def fake_logo(list_of_clonotypes, sequence_type):
        captured["clonotypes"] = list_of_clonotypes
        captured["sequence_type"] = sequence_type
        return "logo"

    logo_module.get_logo_for_list_of_clonotypes = fake_logo
    monkeypatch.setitem(sys.modules, "repseq.logo", logo_module)
    monkeypatch.setattr(bcr, "logo", logo_module, raising=False)

    result = analyzer.get_logo_for_tree(requested_tree_id)

    assert result == "logo"
    assert captured == {
        "clonotypes": [("AAA",), ("ABA",), ("ACA",)],
        "sequence_type": "prot",
    }


def test_get_tree_clonotypes_returns_full_source_rows_and_selected_tree_data(tmp_path):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.read_metadata(
        pd.DataFrame(
            {
                "sample_id": ["sample.one", "two"],
                "timepoint": [1, 2],
                "donor_id": ["donor", "donor"],
            }
        )
    )

    clonotypes = analyzer.get_tree_clonotypes("1")

    assert len(clonotypes) == 2
    assert clonotypes.columns[:7].tolist() == [
        "treeId",
        "sample_id",
        "cloneId",
        "nMutationsRate",
        "timepoint",
        "donor_id",
        "isotype",
    ]
    assert clonotypes["cloneId"].tolist() == [1, 2]
    assert clonotypes["nMutationsRate"].tolist() == pytest.approx([0.1, 0.2])
    assert clonotypes["timepoint"].tolist() == [1, 1]
    assert clonotypes["donor_id"].tolist() == ["donor", "donor"]
    assert clonotypes["isotype"].tolist() == ["IgM", "IgG"]
    assert clonotypes["uniqueMoleculeCount"].tolist() == [10, 2]
    assert clonotypes["customClonosetColumn"].tolist() == ["first", "second"]
    assert "nodeId" not in clonotypes.columns
    assert "aaSeqCDR3" not in clonotypes.columns


def test_get_mutation_positions_parses_substitutions_deletions_and_insertions():
    positions, mutation_types = bcr.get_mutation_positions("SA1GDC2I3T")

    assert positions == [1, 2, 3]
    assert mutation_types == [[1], [2], [3]]


def test_plot_mutations_rate_uses_observed_nodes_and_region_boundaries(tmp_path):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    empty_alignment = "0|0|0|0|0||0"
    analyzer.trees_df["refPoints"] = _mutation_ref_points()
    analyzer.trees_df["allVAlignments"] = [
        "0|15|15|0|15|SC2TSA3G|100;0|15|15|0|15|SA9G|90",
        "0|15|15|0|15|SA3G|100",
        "0|15|15|0|15|SC2T|100",
        empty_alignment,
    ]
    analyzer.trees_df["allDAlignments"] = empty_alignment
    analyzer.trees_df["allJAlignments"] = empty_alignment

    mutation_rates = analyzer._get_mutation_rate_df("1")
    axis = analyzer.plot_mutations_rate(1)

    assert len(mutation_rates) == 13
    assert mutation_rates.loc[mutation_rates["position"] == 0, "rate"].iloc[0] == 0.5
    assert mutation_rates.loc[mutation_rates["position"] == 1, "rate"].iloc[0] == 1.0
    assert mutation_rates.loc[mutation_rates["position"] == 7, "rate"].iloc[0] == 0
    assert mutation_rates.loc[mutation_rates["position"] == 0, "region"].iloc[0] == "CDR1"
    assert mutation_rates.loc[mutation_rates["position"] == 2, "region"].iloc[0] == "FR2"
    assert mutation_rates.loc[mutation_rates["position"] == 8, "region"].iloc[0] == "CDR3"
    assert mutation_rates.loc[mutation_rates["position"] == 11, "region"].iloc[0] == "FR4"
    assert len(axis.patches) == 13
    assert len(axis.lines) == 5
    assert all(line.get_linestyle() == "--" for line in axis.lines)
    assert axis.get_title() == "Tree 1 mutation frequencies"
    assert axis.get_ylabel() == "Mutation frequency"


def _trajectory_analyzer(tmp_path, timepoint_column="timepoint"):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 4, "treeId"] = 1

    early_replicate = analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 1].copy()
    early_replicate["fileName"] = "mix.early.rep.clns"
    early_replicate["sample_id"] = "early.rep"
    early_replicate["cloneId"] = 101
    early_replicate["uniqueMoleculeCount"] = 4
    early_replicate["uniqueMoleculeFraction"] = 0.2

    late_replicate = analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 4].copy()
    late_replicate["fileName"] = "mix.late.rep.clns"
    late_replicate["sample_id"] = "late.rep"
    late_replicate["cloneId"] = 102
    late_replicate["uniqueMoleculeCount"] = 1
    late_replicate["uniqueMoleculeFraction"] = 0.1
    global_only = analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 4].copy()
    global_only["treeId"] = 2
    global_only["fileName"] = "mix.global.only.clns"
    global_only["sample_id"] = "global.only"
    global_only["cloneId"] = 103
    global_only["uniqueMoleculeCount"] = 5
    global_only["uniqueMoleculeFraction"] = 0.4
    analyzer.trees_df = pd.concat(
        [analyzer.trees_df, early_replicate, late_replicate, global_only],
        ignore_index=True,
    )

    metadata = pd.DataFrame(
        {
            "sample_id": ["sample.one", "early.rep", "two", "late.rep", "global.only"],
            timepoint_column: [1, 1, 2, 2, 3],
        }
    )
    analyzer.read_metadata(metadata)
    return analyzer


def test_timepoint_trajectory_summarizes_fraction_and_dispersion(tmp_path):
    analyzer = _trajectory_analyzer(tmp_path)

    trajectory = analyzer._get_timepoint_trajectory_df(1)
    axis = analyzer.timepoint_trajectory(1)

    assert trajectory["timepoint"].tolist() == [1, 2, 3]
    assert trajectory["mean"].tolist() == pytest.approx([0.16, 0.2, 0])
    assert trajectory["minimum"].tolist() == pytest.approx([0.12, 0.1, 0])
    assert trajectory["maximum"].tolist() == pytest.approx([0.2, 0.3, 0])
    assert trajectory["samples"].tolist() == [2, 2, 0]
    assert axis.lines[0].get_ydata().tolist() == pytest.approx([0.16, 0.2, 0])
    assert [tick.get_text() for tick in axis.get_xticklabels()] == ["1", "2", "3"]
    assert axis.get_ylabel() == "Lineage fraction in repertoire by UMI count"


def test_timepoint_trajectory_supports_counts_and_custom_feature(tmp_path):
    analyzer = _trajectory_analyzer(tmp_path, timepoint_column="Timepoints")

    trajectory = analyzer._get_timepoint_trajectory_df(
        1,
        timepoint_feature="Timepoints",
        by_freq=False,
    )
    axis = analyzer.timepoint_trajectory(
        1,
        timepoint_feature="Timepoints",
        by_freq=False,
    )

    assert trajectory["Timepoints"].tolist() == [1, 2, 3]
    assert trajectory["mean"].tolist() == pytest.approx([8, 2, 0])
    assert trajectory["minimum"].tolist() == pytest.approx([4, 1, 0])
    assert trajectory["maximum"].tolist() == pytest.approx([12, 3, 0])
    assert axis.get_xlabel() == "Timepoints"
    assert axis.get_ylabel() == "Lineage UMI count"


def test_timepoint_isotypes_averages_sample_values_and_zero_fills(tmp_path):
    analyzer = _trajectory_analyzer(tmp_path)

    trajectory = analyzer._get_timepoint_isotypes_df(1)
    axis = analyzer.timepoint_isotypes(1)
    values = trajectory.pivot(index="timepoint", columns="isotype", values="mean")

    assert list(trajectory["isotype"].cat.categories) == ["IgM", "IgD", "IgG"]
    assert values.loc[1, "IgM"] == pytest.approx(0.15)
    assert values.loc[1, "IgD"] == pytest.approx(0)
    assert values.loc[1, "IgG"] == pytest.approx(0.01)
    assert values.loc[2, "IgM"] == pytest.approx(0)
    assert values.loc[2, "IgD"] == pytest.approx(0.2)
    assert values.loc[2, "IgG"] == pytest.approx(0)
    assert values.loc[3].tolist() == pytest.approx([0, 0, 0])
    assert [line.get_label() for line in axis.lines] == ["IgM", "IgD", "IgG"]
    assert axis.lines[0].get_ydata().tolist() == pytest.approx([0.15, 0, 0])
    assert axis.lines[1].get_ydata().tolist() == pytest.approx([0, 0.2, 0])
    assert axis.lines[2].get_ydata().tolist() == pytest.approx([0.01, 0, 0])
    assert np.allclose(to_rgba(axis.lines[0].get_color()), to_rgba("#E41A1C"))
    assert np.allclose(to_rgba(axis.lines[1].get_color()), to_rgba("#FF7F00"))
    assert np.allclose(to_rgba(axis.lines[2].get_color()), to_rgba("#4DAF4A"))
    assert axis.get_ylabel() == "Mean isotype fraction by UMI count"
    assert not axis.collections


def test_timepoint_isotypes_supports_counts(tmp_path):
    analyzer = _trajectory_analyzer(tmp_path, timepoint_column="Timepoints")

    trajectory = analyzer._get_timepoint_isotypes_df(
        1,
        timepoint_feature="Timepoints",
        by_freq=False,
    )
    axis = analyzer.timepoint_isotypes(
        1,
        timepoint_feature="Timepoints",
        by_freq=False,
    )
    values = trajectory.pivot(index="Timepoints", columns="isotype", values="mean")

    assert values.loc[1, "IgM"] == pytest.approx(7)
    assert values.loc[1, "IgG"] == pytest.approx(1)
    assert values.loc[2, "IgD"] == pytest.approx(2)
    assert values.loc[3].tolist() == pytest.approx([0, 0, 0])
    assert axis.get_xlabel() == "Timepoints"
    assert axis.get_ylabel() == "Mean isotype UMI count"


def test_timepoint_trajectory_without_metadata_prints_instructions(tmp_path, capsys):
    trees_filename, _ = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    capsys.readouterr()

    result = analyzer.timepoint_trajectory(1)

    assert result is None
    assert "run ta.read_metadata(metadata)" in capsys.readouterr().out


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
    point_collections = [
        collection
        for collection in axis.collections
        if isinstance(collection, PathCollection)
    ]
    assert len(point_collections) == 3
    observed_points = [
        collection for collection in point_collections
        if collection.get_sizes()[0] > 22
    ]
    reconstructed_points = [
        collection for collection in point_collections
        if collection.get_sizes()[0] == 22
    ]
    observed_x = [collection.get_offsets()[0, 0] for collection in observed_points]
    observed_xy = [tuple(collection.get_offsets()[0]) for collection in observed_points]
    reconstructed_x = reconstructed_points[0].get_offsets()[0, 0]
    reconstructed_y = reconstructed_points[0].get_offsets()[0, 1]
    assert sorted(observed_x) == pytest.approx([0.1, 0.2])
    assert sorted(observed_xy) == pytest.approx([(0.1, 1.0), (0.2, 2.0)])
    assert reconstructed_x == pytest.approx(0)
    assert reconstructed_y == pytest.approx(1.5)
    assert axis.get_xlabel() == "Branch length"


def test_draw_tree_uses_isotype_fraction_palette_and_order(tmp_path):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.read_trees_newick(newick_dir)

    axis = analyzer.draw_tree(1, group="isotype", label=None)
    observed = [
        collection
        for collection in axis.collections
        if isinstance(collection, PathCollection) and collection.get_sizes()[0] > 22
    ]
    colors_by_x = {
        float(collection.get_offsets()[0, 0]): collection.get_facecolors()[0]
        for collection in observed
    }

    assert np.allclose(colors_by_x[0.1], to_rgba("#E41A1C"))
    assert np.allclose(colors_by_x[0.2], to_rgba("#4DAF4A"))
    assert [text.get_text() for text in axis.get_legend().get_texts()] == [
        "IgM",
        "IgG",
    ]


@pytest.mark.parametrize("requested_tree_id", [6388, "6388"])
def test_draw_tree_normalizes_integer_string_and_float_tree_ids(
    tmp_path, requested_tree_id
):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    analyzer.trees_df["treeId"] = analyzer.trees_df["treeId"].astype(float)
    analyzer.trees_df.loc[analyzer.trees_df["treeId"] == 1.0, "treeId"] = 6388.0
    analyzer.trees_df["nodeId"] = analyzer.trees_df["nodeId"].astype(float)
    (newick_dir / "1.tree").rename(newick_dir / "6388.tree")
    analyzer.read_trees_newick(newick_dir)

    axis = analyzer.draw_tree(requested_tree_id)

    assert axis.get_title() == f"Tree {requested_tree_id}"
    assert sum(isinstance(collection, PathCollection) for collection in axis.collections) == 3


def test_phylo_positions_preserve_actual_branch_lengths():
    tree = Phylo.read(StringIO("((2:1)1:1)3;"), "newick")
    rsplot._restore_numeric_internal_node_names(tree, {"1", "2", "3"})
    x_positions, _ = rsplot._phylo_positions(tree)

    observed_internal = next(clade for clade in tree.find_clades() if clade.name == "1")
    reconstructed_root = next(clade for clade in tree.find_clades() if clade.name == "3")
    observed_terminal = next(clade for clade in tree.find_clades() if clade.name == "2")
    assert x_positions[reconstructed_root] == pytest.approx(0)
    assert x_positions[observed_internal] == pytest.approx(1)
    assert x_positions[observed_terminal] == pytest.approx(2)
    assert not observed_internal.is_terminal()
    assert not reconstructed_root.is_terminal()


def test_draw_tree_supports_multiple_observations_for_one_node(tmp_path):
    trees_filename, newick_dir = _write_tree_inputs(tmp_path)
    analyzer = bcr.TreeAnalyzer()
    analyzer.read_trees_table(trees_filename)
    duplicate = analyzer.trees_df.loc[analyzer.trees_df["nodeId"] == 1].copy()
    duplicate["sample_id"] = "sample.duplicate"
    duplicate["fileName"] = "mix.sample.duplicate.clns"
    duplicate["isotype"] = "IgA"
    duplicate["uniqueMoleculeCount"] = 4
    analyzer.trees_df = pd.concat([analyzer.trees_df, duplicate], ignore_index=True)
    analyzer.read_trees_newick(newick_dir)

    axis = analyzer.draw_tree(1)

    assert axis.get_title() == "Tree 1"
    assert sum(isinstance(collection, PathCollection) for collection in axis.collections) == 4
