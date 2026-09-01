import matplotlib

matplotlib.use("Agg")

import copy
import matplotlib.pyplot as plt
import numpy as np
import sys
import types
import warnings
import pandas as pd
import pytest
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgba

import repseq.clustering as clustering_module
import repseq.plot as rsplot
from repseq.clustering import (
    Cluster,
    Clusters,
    Node,
    all_nodes,
    any_nodes,
    cluster_size,
    proportion,
    total_count,
    top_clusters,
)
from repseq.clone_filter import Filter


def _clusters_with_two_samples():
    clusters = Clusters()

    cluster_0 = Cluster()
    nodes_0 = [
        Node(0, "TGTGCT", "CASS", "TRBV1", "TRBJ1", "sample_1", 0.1, 3),
        Node(1, "TGTGCT", "CASS", "TRBV1", "TRBJ1", "sample_1", 0.2, 7),
        Node(2, "TGTGCT", "CASS", "TRBV1", "TRBJ1", "sample_2", 0.4, 5),
    ]
    cluster_0.add_nodes_from(nodes_0)
    cluster_0.add_edges_from([(nodes_0[0], nodes_0[1]), (nodes_0[1], nodes_0[2])])

    cluster_1 = Cluster()
    node_1 = Node(3, "TGTGCC", "CATS", "TRBV2", "TRBJ2", "sample_2", 0.3, 11)
    cluster_1.add_node(node_1)

    clusters.clusters = [cluster_0, cluster_1]
    clusters.clonosets_df = pd.DataFrame(
        {"sample_id": ["sample_1", "sample_2", "sample_3"]}
    )
    for cluster_no, cluster in enumerate(clusters):
        for node in cluster:
            node.additional_properties["cluster_no"] = cluster_no

    return clusters


def test_to_count_table_sums_counts_by_cluster_and_sample():
    count_table = _clusters_with_two_samples().to_count_table()

    expected = pd.DataFrame(
        {
            "cluster_id": ["cluster_0", "cluster_1"],
            "consensus": ["CASS|TRBV1|TRBJ1", "CATS|TRBV2|TRBJ2"],
            "concensus_cdr3aa": ["CASS", "CATS"],
            "concensus_v": ["TRBV1", "TRBV2"],
            "concensus_j": ["TRBJ1", "TRBJ2"],
            "sample_1": [10, 0],
            "sample_2": [5, 11],
            "sample_3": [0, 0],
        }
    )

    pd.testing.assert_frame_equal(count_table, expected)


def test_to_count_table_sums_frequencies_when_requested():
    count_table = _clusters_with_two_samples().to_count_table(by_freq=True)

    assert count_table.columns.tolist() == [
        "cluster_id",
        "consensus",
        "concensus_cdr3aa",
        "concensus_v",
        "concensus_j",
        "sample_1",
        "sample_2",
        "sample_3",
    ]
    np.testing.assert_allclose(count_table["sample_1"], [0.3, 0])
    np.testing.assert_allclose(count_table["sample_2"], [0.4, 0.3])
    np.testing.assert_allclose(count_table["sample_3"], [0, 0])


@pytest.mark.parametrize(
    ("overlap_type", "expected_groups", "expected_edges"),
    [
        (
            "aaV",
            [
                {("TRBV1", "TRBJ1"), ("TRBV1", "TRBJ2")},
                {("TRBV2", "TRBJ1")},
            ],
            3,
        ),
        (
            "aaVJ",
            [
                {("TRBV1", "TRBJ1")},
                {("TRBV1", "TRBJ2")},
                {("TRBV2", "TRBJ1")},
            ],
            1,
        ),
    ],
)
def test_find_nodes_and_edges_groups_by_required_segments(
    monkeypatch, overlap_type, expected_groups, expected_edges
):
    clusters = Clusters()
    clusters.clonotypes = pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CATS", "CARS", "CAGS"],
            "cdr3nt": ["TGTGCT", "TGTGCC", "TGTGCA", "TGTGCG"],
            "v": ["TRBV1", "TRBV1", "TRBV1", "TRBV2"],
            "j": ["TRBJ1", "TRBJ1", "TRBJ2", "TRBJ1"],
            "sample_id": ["sample_1"] * 4,
            "freq": [0.25] * 4,
            "count": [1] * 4,
        }
    )
    compared_groups = []

    def capture_tasks(function, tasks, *args, **kwargs):
        compared_groups.extend(
            [{(node.v, node.j) for node in task[0]} for task in tasks]
        )
        return [function(task) for task in tasks]

    monkeypatch.setattr(clustering_module, "run_parallel_calculation", capture_tasks)

    _, edges = clusters.find_nodes_and_edges(
        mismatches=1, overlap_type=overlap_type, cpu=1, verbose=False
    )

    assert compared_groups == expected_groups
    assert len(edges) == expected_edges


def test_clusters_str_reports_state_graph_and_cluster_counts():
    clusters = _clusters_with_two_samples()
    clusters.clonotypes = clusters.clonosets_df

    summary = str(clusters)

    assert "Completed: clonotypes read, clusters created" in summary
    assert "Graph: 4 nodes and 2 edges" in summary
    assert (
        "Clusters: 1 cluster (2 or more nodes) and 1 single node. Total: 2"
        in summary
    )
    assert "Metadata: not added" in summary
    assert "Node Pgen: not calculated" in summary
    assert "ALICE: not calculated" in summary


def test_plot_cluster_facets_style_nodes_and_add_legends():
    clusters = _clusters_with_two_samples()
    color_groups = ["control", "control", "case", "case"]
    shape_groups = ["alpha", "beta", "alpha", "beta"]
    for node, color_group, shape_group in zip(
        [node for cluster in clusters for node in cluster],
        color_groups,
        shape_groups,
    ):
        node.additional_properties["color_group"] = color_group
        node.additional_properties["shape_group"] = shape_group

    figure = clusters.plot_cluster(
        [0, 1],
        layout="circular",
        color="color_group",
        palette={"control": "#112233", "case": "#AABBCC"},
        label="id",
        shape="shape_group",
        ncols=2,
        size="count",
        min_size=100,
        log_scaled=True,
        linear_scale=10,
    )

    assert [axis.get_title() for axis in figure.axes] == ["cluster_0", "cluster_1"]
    assert [legend.get_title().get_text() for legend in figure.legends] == [
        "color_group",
        "shape_group",
    ]
    assert [text.get_text() for text in figure.axes[0].texts] == ["0", "1", "2"]
    assert [text.get_text() for text in figure.axes[1].texts] == ["3"]

    plotted_sizes = []
    for axis in figure.axes:
        for collection in axis.collections:
            if isinstance(collection, PathCollection):
                plotted_sizes.extend(collection.get_sizes())
    expected_sizes = [
        100 + 10 * np.log2(count + 1) for count in [3, 7, 5, 11]
    ]
    np.testing.assert_allclose(sorted(plotted_sizes), sorted(expected_sizes))
    assert figure.number not in plt.get_fignums()


def test_plot_cluster_orders_color_and_shape_levels_by_source_columns():
    clusters = Clusters()
    cluster = Cluster()
    nodes = [
        Node(index, "TGTGCT", "CASS", "TRBV1", "TRBJ1", "sample_1", 0.1, 1)
        for index in range(3)
    ]
    color_groups = ["first", "second", "third"]
    shape_groups = ["circle", "triangle", "diamond"]
    for node, color_group, shape_group in zip(nodes, color_groups, shape_groups):
        node.additional_properties.update(
            {"color_group": color_group, "shape_group": shape_group}
        )
    cluster.add_nodes_from([nodes[2], nodes[0], nodes[1]])
    clusters.clusters = [cluster]
    clusters.clonotypes = pd.DataFrame(
        {"color_group": color_groups, "shape_group": shape_groups}
    )

    figure = clusters.plot_cluster(
        0, color="color_group", shape="shape_group", size=None
    )

    assert [text.get_text() for text in figure.legends[0].get_texts()] == color_groups
    assert [text.get_text() for text in figure.legends[1].get_texts()] == shape_groups
    assert [
        handle.get_marker() for handle in figure.legends[1].legend_handles
    ] == ["o", "^", "D"]
    plt.close(figure)


@pytest.mark.parametrize(
    "layout", ["spring", "kamada_kawai", "circular", "shell", "spectral"]
)
def test_plot_cluster_supports_common_layouts(layout):
    figure = _clusters_with_two_samples().plot_cluster(0, layout=layout)

    assert figure.axes[0].get_title() == "cluster_0"
    plt.close(figure)


def test_plot_cluster_none_plots_all_clusters():
    figure = _clusters_with_two_samples().plot_cluster(None)

    assert [axis.get_title() for axis in figure.axes] == [
        "cluster_0",
        "cluster_1",
    ]
    plt.close(figure)


def test_plot_cluster_none_respects_max_clusters_without_plotting():
    clusters = _clusters_with_two_samples()
    open_figures = plt.get_fignums()

    with pytest.raises(
        ValueError,
        match=(
            r"would plot all 2 clusters.*max_clusters=1.*No plot was created.*"
            r"Increase max_clusters.*smaller selection"
        ),
    ):
        clusters.plot_cluster(None, max_clusters=1)

    assert plt.get_fignums() == open_figures


def test_plot_cluster_ignores_max_clusters_for_explicit_selection():
    figure = _clusters_with_two_samples().plot_cluster([0, 1], max_clusters=1)

    assert [axis.get_title() for axis in figure.axes] == [
        "cluster_0",
        "cluster_1",
    ]
    plt.close(figure)


def test_plot_cluster_none_can_raise_max_clusters_above_default():
    clusters = Clusters()
    for cluster_no in range(52):
        cluster = Cluster()
        cluster.id = cluster_no
        node = Node(
            cluster_no,
            "TGT",
            "CASS",
            "TRBV1",
            "TRBJ1",
            "sample_1",
            1.0,
            1,
        )
        node.additional_properties["cluster_no"] = cluster_no
        cluster.add_node(node)
        clusters.clusters.append(cluster)

    figure = clusters.plot_cluster(
        None,
        max_clusters=52,
        layout="circular",
        ncols=10,
        figsize=(10, 6),
        size=None,
    )

    assert len(figure.axes) == 60
    assert figure.axes[51].get_title() == "cluster_51"
    assert all(not axis.get_visible() for axis in figure.axes[52:])
    plt.close(figure)


@pytest.mark.parametrize("max_clusters", [0, -1, 1.5, True, None])
def test_plot_cluster_none_requires_positive_integer_max_clusters(max_clusters):
    with pytest.raises(
        ValueError, match="max_clusters must be a positive integer"
    ):
        _clusters_with_two_samples().plot_cluster(None, max_clusters=max_clusters)


def test_plot_cluster_rejects_more_than_fifty_facets():
    with pytest.raises(ValueError, match="At most 50 clusters"):
        _clusters_with_two_samples().plot_cluster(list(range(51)))


def test_plot_cluster_reads_style_only_from_known_node_properties():
    with pytest.raises(ValueError, match="Node property 'unknown'"):
        _clusters_with_two_samples().plot_cluster(0, color="unknown")


def test_plot_cluster_uses_five_default_shape_levels():
    clusters = Clusters()
    cluster = Cluster()
    shape_levels = ["circle", "triangle", "rhombus", "hexagon", "square"]
    for index, shape_level in enumerate(shape_levels):
        node = Node(
            index,
            "TGTGCT",
            "CASS",
            "TRBV1",
            "TRBJ1",
            "sample_1",
            0.2,
            index + 1,
        )
        node.additional_properties["shape_group"] = shape_level
        cluster.add_node(node)
    clusters.clusters = [cluster]

    figure = clusters.plot_cluster(0, shape="shape_group")

    markers = [handle.get_marker() for handle in figure.legends[0].legend_handles]
    assert markers == ["o", "^", "D", "h", "s"]
    plt.close(figure)


def test_plot_cluster_rejects_more_than_five_shape_levels():
    clusters = Clusters()
    cluster = Cluster()
    for index in range(6):
        node = Node(
            index,
            "TGTGCT",
            "CASS",
            "TRBV1",
            "TRBJ1",
            "sample_1",
            1 / 6,
            1,
        )
        node.additional_properties["shape_group"] = index
        cluster.add_node(node)
    clusters.clusters = [cluster]

    with pytest.raises(ValueError, match="shape supports at most five levels"):
        clusters.plot_cluster(0, shape="shape_group")


def _plotted_node_sizes(figure):
    sizes = []
    for axis in figure.axes:
        for collection in axis.collections:
            if isinstance(collection, PathCollection):
                sizes.extend(collection.get_sizes())
    return sizes


def test_plot_cluster_supports_linear_and_uniform_node_sizes():
    clusters = _clusters_with_two_samples()

    linear_figure = clusters.plot_cluster(
        0,
        size="freq",
        min_size=20,
        log_scaled=False,
        linear_scale=100,
    )
    uniform_figure = clusters.plot_cluster(0, size=None, min_size=45)

    np.testing.assert_allclose(
        sorted(_plotted_node_sizes(linear_figure)),
        sorted([30, 40, 60]),
    )
    np.testing.assert_allclose(_plotted_node_sizes(uniform_figure), [45, 45, 45])
    assert linear_figure.number not in plt.get_fignums()
    assert uniform_figure.number not in plt.get_fignums()


def test_find_nodes_and_edges_recodes_igh_c_as_isotype():
    clusters = Clusters()
    clusters.clonotypes = pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CATS"],
            "cdr3nt": ["TGTGCT", "TGTGCC"],
            "v": ["IGHV1", "IGHV1"],
            "j": ["IGHJ1", "IGHJ1"],
            "c": ["IGHG1*01", "IGHM*02"],
            "sample_id": ["sample_1", "sample_2"],
            "freq": [0.4, 0.6],
            "count": [4, 6],
        }
    )

    single_nodes, edges = clusters.find_nodes_and_edges(
        mismatches=1, overlap_type="aaVJ", cpu=1, verbose=False
    )
    nodes = set(single_nodes)
    for node_1, node_2, _ in edges:
        nodes.update([node_1, node_2])

    assert sorted(node.additional_properties["isotype"] for node in nodes) == [
        "IgG1",
        "IgM",
    ]
    assert all("c" not in node.additional_properties for node in nodes)


def test_plot_cluster_uses_isotype_fraction_order_and_colors():
    clusters = _clusters_with_two_samples()
    isotypes = ["IgG2", "IgM", "IgA1", "IgG1"]
    for node, isotype in zip(
        [node for cluster in clusters for node in cluster], isotypes
    ):
        node.additional_properties["isotype"] = isotype

    figure = clusters.plot_cluster([0, 1], color="isotype", size=None)

    legend = figure.legends[0]
    labels = [text.get_text() for text in legend.get_texts()]
    expected_order = rsplot._isotype_order(isotypes)
    expected_colors = rsplot._isotype_colors(expected_order)
    plotted_colors = [
        to_rgba(handle.get_markerfacecolor()) for handle in legend.legend_handles
    ]
    assert labels == expected_order
    for isotype, plotted_color in zip(expected_order, plotted_colors):
        assert np.allclose(plotted_color, to_rgba(expected_colors[isotype]))


def _clusters_with_filter_properties():
    clusters = _clusters_with_two_samples()
    nodes = [node for cluster in clusters for node in cluster]
    groups = ["group_1", "group_2", "other", "control"]
    isotypes = ["IgM", "IgD", "IgM", "IgG1"]
    specificities = ["other", "other", "other", "specific"]
    custom_weights = [2.0, 3.0, 5.0, 7.0]
    for node, group, isotype, specificity, custom_weight in zip(
        nodes, groups, isotypes, specificities, custom_weights
    ):
        node.additional_properties.update(
            {
                "group_property_name": group,
                "isotype": isotype,
                "specificity": specificity,
                "custom_weight": custom_weight,
            }
        )
    return clusters


def test_clusters_filter_supports_metrics_and_weighted_proportions():
    clusters = _clusters_with_filter_properties()

    filtered = clusters.filter(
        (cluster_size >= 3)
        & (
            proportion(
                "group_property_name",
                ["group_1", "group_2"],
                weight="count",
            )
            > 0.25
        )
        & (
            total_count(
                "group_property_name", "control", weight="freq"
            )
            == 0
        )
    )

    assert list(filtered) == [clusters[0]]
    assert filtered[0].id is clusters[0].id
    assert len(clusters) == 2


def test_clusters_filter_supports_all_any_or_and_not():
    clusters = _clusters_with_filter_properties()

    all_igm_or_igd = clusters.filter(all_nodes("isotype", ["IgM", "IgD"]))
    specific_or_large = clusters.filter(
        any_nodes("specificity", "specific") | (cluster_size >= 3)
    )
    without_trbv2 = clusters.filter(~any_nodes("v", "TRBV2"))

    assert list(all_igm_or_igd) == [clusters[0]]
    assert list(specific_or_large) == [clusters[0], clusters[1]]
    assert list(without_trbv2) == [clusters[0]]


def test_cluster_filter_uses_cdr3_sequence_aliases():
    clusters = _clusters_with_filter_properties()

    aa_match = clusters.filter(any_nodes("cdr3aa", "CATS"))
    nt_match = clusters.filter(any_nodes("cdr3nt", "TGTGCC"))

    assert list(aa_match) == [clusters[1]]
    assert clusters[1] in list(nt_match)


def test_cluster_aggregates_support_nodes_and_custom_numeric_weights():
    clusters = _clusters_with_filter_properties()

    by_nodes = total_count(
        "group_property_name", ["group_1", "group_2"], weight="nodes"
    )
    by_custom_weight = proportion(
        "group_property_name", "group_1", weight="custom_weight"
    )
    properties = clusters.custom_properties(
        [by_nodes.alias("selected_nodes"), by_custom_weight]
    )

    assert properties["selected_nodes"].tolist() == [2, 0]
    np.testing.assert_allclose(
        properties[by_custom_weight.name],
        [2 / 10, 0],
    )


def test_custom_properties_returns_identifiers_and_requested_metrics():
    clusters = _clusters_with_filter_properties()
    group_1_proportion = proportion(
        "group_property_name", "group_1", weight="count"
    )

    properties = clusters.custom_properties(
        [cluster_size, total_count, group_1_proportion]
    )

    assert properties.columns.tolist() == [
        "cluster_no",
        "cluster_id",
        "cluster_size",
        "total_count",
        group_1_proportion.name,
    ]
    assert properties["cluster_no"].tolist() == [0, 1]
    assert properties["cluster_id"].tolist() == ["cluster_0", "cluster_1"]
    assert properties["cluster_size"].tolist() == [3, 1]
    assert properties["total_count"].tolist() == [15, 11]
    np.testing.assert_allclose(properties[group_1_proportion.name], [3 / 15, 0])


def test_plot_cluster_understands_cdr3aa_and_cdr3nt_aliases():
    clusters = _clusters_with_two_samples()

    aa_figure = clusters.plot_cluster(1, label="cdr3aa", size=None)
    nt_figure = clusters.plot_cluster(1, label="cdr3nt", size=None)

    assert [text.get_text() for text in aa_figure.axes[0].texts] == ["CATS"]
    assert [text.get_text() for text in nt_figure.axes[0].texts] == ["TGTGCC"]


def test_cluster_filter_requires_boolean_expression():
    with pytest.raises(TypeError, match="boolean cluster expression"):
        _clusters_with_filter_properties().filter(cluster_size)


def test_filtered_custom_properties_preserve_original_cluster_number():
    clusters = _clusters_with_filter_properties()
    filtered = clusters.filter(any_nodes("cdr3aa", "CATS"))

    properties = filtered.custom_properties([cluster_size])

    assert properties["cluster_no"].tolist() == [1]
    assert properties["cluster_id"].tolist() == ["cluster_1"]


def test_cluster_boolean_expressions_short_circuit():
    clusters = _clusters_with_filter_properties()

    selected = clusters.filter(
        (cluster_size >= 1) | any_nodes("unknown_property", "value")
    )

    assert list(selected) == list(clusters)


def test_properties_includes_total_count_between_edges_and_diameter():
    properties = _clusters_with_two_samples().properties(cpu=1)

    assert properties.columns.tolist()[:7] == [
        "cluster_no",
        "cluster_id",
        "nodes",
        "edges",
        "total_count",
        "diameter",
        "density",
    ]
    assert properties["total_count"].tolist() == [15, 11]


def test_properties_is_method_without_legacy_property_behavior():
    clusters = _clusters_with_two_samples()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        properties_method = clusters.properties

    assert isinstance(Clusters.__dict__["properties"], types.FunctionType)
    assert isinstance(properties_method, types.MethodType)
    assert caught == []
    with pytest.raises(AttributeError):
        _ = properties_method.columns


@pytest.mark.parametrize(
    ("weight", "expected_weights"),
    [
        ("count", [3, 7, 5]),
        ("freq", [0.1, 0.2, 0.4]),
        ("nodes", [1, 1, 1]),
        ("custom_weight", [2.0, 3.0, 5.0]),
    ],
)
def test_plot_logo_supports_builtin_and_custom_weights(
    monkeypatch, weight, expected_weights
):
    clusters = _clusters_with_filter_properties()
    calls = []

    def fake_logo(clonotypes, seq_type, plot=True):
        calls.append((clonotypes, seq_type, plot))
        return "logo_result"

    monkeypatch.setattr(
        clustering_module, "get_logo_for_list_of_clonotypes", fake_logo
    )

    result = clusters.plot_logo(0, weight=weight, plot=False)

    assert result == "logo_result"
    assert calls[0][1:] == ("prot", False)
    assert [sequence for sequence, _ in calls[0][0]] == ["CASS"] * 3
    np.testing.assert_allclose(
        [node_weight for _, node_weight in calls[0][0]], expected_weights
    )


def test_plot_logo_supports_dna_sequences(monkeypatch):
    clusters = _clusters_with_two_samples()
    captured = {}

    def fake_logo(clonotypes, seq_type, plot=True):
        captured["clonotypes"] = clonotypes
        captured["seq_type"] = seq_type
        captured["plot"] = plot
        return "dna_logo"

    monkeypatch.setattr(
        clustering_module, "get_logo_for_list_of_clonotypes", fake_logo
    )

    result = clusters.plot_logo(1, seq_type="dna", weight="nodes")

    assert result == "dna_logo"
    assert captured == {
        "clonotypes": [("TGTGCC", 1)],
        "seq_type": "dna",
        "plot": True,
    }


def _pooled_clonotypes_for_state_tests():
    return pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CATS"],
            "cdr3nt": ["TGTGCT", "TGTGCC"],
            "v": ["TRBV1", "TRBV1"],
            "j": ["TRBJ1", "TRBJ1"],
            "sample_id": ["sample_1", "sample_2"],
            "freq": [0.4, 0.6],
            "count": [4, 6],
        }
    )


def test_empty_state_and_prerequisite_messages():
    clusters = Clusters()

    summary = str(clusters)

    assert "State: empty" in summary
    assert "Clonotypes: not read" in summary
    assert "Clusters: not created" in summary
    with pytest.raises(RuntimeError, match="clonotypes have not been read"):
        clusters.create_clusters(verbose=False)
    with pytest.raises(RuntimeError, match="clonotypes have not been read"):
        clusters.add_metadata(pd.DataFrame({"sample_id": ["sample_1"]}))
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        clusters.plot_cluster(0)
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        clusters.as_dataframe()
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        clusters.filter(cluster_size >= 1)
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        clusters.save_to_cytoscape("unused")
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        _ = clusters.properties(cpu=1)
    with pytest.raises(RuntimeError, match="clusters have not been created"):
        clusters.alice(overlap_type="aaVJ", mismatches=1)


def test_pooled_read_tracks_source_counts_and_state():
    clusters = Clusters()
    result = clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests()
    )
    assert result is None

    assert clusters.state["clonotypes_read"]
    assert not clusters.state["clusters_created"]
    assert clusters.state_parameters["clonotypes_read"] == {
        "source": "pooled_clonoset",
        "clonotypes": 2,
        "samples": 2,
        "converted": False,
    }
    summary = str(clusters)
    assert "source: pooled_clonoset" in summary
    assert "clonotypes: 2" in summary
    assert "samples: 2" in summary


def test_clonosets_df_read_tracks_filter_parameters(monkeypatch):
    source = _pooled_clonotypes_for_state_tests().drop(columns="sample_id")
    monkeypatch.setattr(
        clustering_module, "read_clonoset", lambda filename: source.copy()
    )
    clonosets_df = pd.DataFrame(
        {"sample_id": ["sample_1"], "filename": ["sample.tsv"]}
    )
    cl_filter = Filter(count_threshold=2)
    clusters = Clusters()
    result = clusters.read_from_clonosets_df(
        clonosets_df, cl_filter=cl_filter, verbose=False
    )
    assert result is None

    parameters = clusters.state_parameters["clonotypes_read"]
    assert parameters["source"] == "clonosets_df"
    assert parameters["input_samples"] == 1
    assert parameters["filter"]["count_threshold"] == 2
    assert "count_threshold: 2" in str(clusters)


def test_metadata_can_be_added_before_clustering_and_is_applied_later():
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests()
    )
    clusters.add_metadata(
        pd.DataFrame(
            {
                "sample_id": ["sample_1", "sample_2"],
                "group": ["control", "case"],
                "timepoint": [0, 1],
            }
        )
    )

    assert clusters.state["metadata_added"]
    assert clusters.state_parameters["metadata_added"]["columns"] == [
        "group",
        "timepoint",
    ]
    clusters.create_clusters(
        overlap_type="aaVJ", mismatches=1, cpu=1, verbose=False
    )

    groups = {
        node.sample_id: node.additional_properties["group"]
        for cluster in clusters
        for node in cluster
    }
    assert groups == {"sample_1": "control", "sample_2": "case"}
    assert "Metadata columns: group, timepoint" in str(clusters)


def test_create_clusters_tracks_parameters_and_reuses_identical_result(
    monkeypatch, capsys
):
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests()
    )
    clusters.create_clusters(
        overlap_type="aaVJ", mismatches=1, cpu=1, verbose=False
    )
    original_clusters = clusters.clusters

    def fail_if_recalculated(*args, **kwargs):
        raise AssertionError("clustering was recalculated")

    monkeypatch.setattr(clusters, "find_nodes_and_edges", fail_if_recalculated)
    result = clusters.create_clusters(
        overlap_type="aaVJ", mismatches=1, cpu=1, verbosity=True
    )

    assert result is None
    assert clusters.clusters is original_clusters
    assert clusters.state_parameters["clusters_created"] == {
        "method": "mismatches",
        "overlap_type": "aaVJ",
        "mismatches": 1,
    }
    output = capsys.readouterr().out
    assert "already created from the same clonotypes" in output


def test_properties_are_cached_and_return_defensive_copies(monkeypatch):
    clusters = _clusters_with_two_samples()
    original_eccentricity = clustering_module.nx.eccentricity
    calls = []

    def counted_eccentricity(cluster):
        calls.append(cluster)
        return original_eccentricity(cluster)

    monkeypatch.setattr(
        clustering_module.nx, "eccentricity", counted_eccentricity
    )
    first = clusters.properties(cpu=1)
    second = clusters.properties(cpu=1)
    first.loc[0, "total_count"] = -1
    third = clusters.properties(cpu=1)

    assert len(calls) == 1
    assert second["total_count"].tolist() == [15, 11]
    assert third["total_count"].tolist() == [15, 11]


def test_alice_tracks_pgen_and_parameters(monkeypatch):
    clusters = _clusters_with_two_samples()
    fake_module = types.ModuleType("repseq.pgen_calculation")

    def fake_calculate_clonotypes_pgen(clonosets_df, **kwargs):
        result = clonosets_df.copy()
        result["pgen"] = 0.001
        result["alice_neighbour_count"] = 0
        result["n_neighbours"] = 0
        return result

    fake_module.calculate_clonotypes_pgen = fake_calculate_clonotypes_pgen
    monkeypatch.setitem(
        sys.modules, "repseq.pgen_calculation", fake_module
    )

    clusters.alice(
        overlap_type="aaVJ",
        mismatches=1,
        generation_model="test_model",
        Q=1,
        alpha=0.05,
    )

    assert clusters.state["node_pgen_calculated"]
    assert clusters.state["alice_calculated"]
    assert (
        clusters.state_parameters["alice_calculated"]["generation_model"]
        == "test_model"
    )
    assert "Node Pgen: calculated" in str(clusters)
    assert "ALICE: calculated" in str(clusters)


def test_plot_cluster_defaults_to_linear_count_sizes():
    figure = _clusters_with_two_samples().plot_cluster(0)

    np.testing.assert_allclose(
        sorted(_plotted_node_sizes(figure)),
        sorted([53, 57, 55]),
    )


def test_plot_cluster_log_scaling_supports_zero_values():
    clusters = _clusters_with_two_samples()
    next(iter(clusters[0])).count = 0

    figure = clusters.plot_cluster(
        0, log_scaled=True, min_size=50, linear_scale=10
    )

    assert min(_plotted_node_sizes(figure)) == 50


def test_select_accepts_mixed_numbers_ids_ranges_and_dataframe_columns():
    clusters = _clusters_with_two_samples()

    mixed = clusters.select(["cluster_1", 0, "cluster_1"])
    by_range = clusters.select(range(2))
    properties = clusters.custom_properties([cluster_size])
    by_dataframe_ids = clusters.select(
        properties.loc[properties["cluster_no"] == 1, "cluster_id"]
    )

    assert list(mixed) == [clusters[1], clusters[0]]
    assert list(by_range) == [clusters[0], clusters[1]]
    assert list(by_dataframe_ids) == [clusters[1]]
    assert mixed.custom_properties([cluster_size])["cluster_no"].tolist() == [1, 0]


def test_select_rejects_unknown_or_malformed_identifiers():
    clusters = _clusters_with_two_samples()

    with pytest.raises(KeyError, match="cluster_25"):
        clusters.select(["cluster_0", 25])
    with pytest.raises(KeyError, match="bad_id"):
        clusters.select(["bad_id"])


def _clusters_for_top_ranking():
    clusters = Clusters()
    specifications = [
        (0, 2, 2, "BBB"),
        (1, 1, 100, "ZZZ"),
        (2, 2, 5, "CCC"),
        (3, 2, 5, "AAA"),
    ]
    for cluster_no, node_count, total_node_count, sequence in specifications:
        cluster = Cluster()
        counts = [1] * node_count
        counts[0] += total_node_count - node_count
        for node_index, count in enumerate(counts):
            node = Node(
                f"{cluster_no}_{node_index}",
                "TGT",
                sequence,
                "TRBV1",
                "TRBJ1",
                "sample_1",
                1 / node_count,
                count,
            )
            node.additional_properties["cluster_no"] = cluster_no
            cluster.add_node(node)
        cluster.id = cluster_no
        clusters.clusters.append(cluster)
    return clusters


def test_filter_top_clusters_uses_size_count_and_consensus_ranking():
    clusters = _clusters_for_top_ranking()

    selected = clusters.filter(top_clusters(3))
    convenient = clusters.top_clusters(2)

    assert [cluster.id for cluster in selected] == [3, 2, 0]
    assert [cluster.id for cluster in convenient] == [3, 2]
    assert clusters.custom_properties([cluster_size])["cluster_no"].tolist() == [
        0,
        1,
        2,
        3,
    ]


def test_top_clusters_validates_requested_number():
    with pytest.raises(ValueError, match="non-negative integer"):
        top_clusters(-1)


def test_create_clusters_keeps_one_original_node_per_clonotype_after_parallel_copy(
    monkeypatch,
):
    pooled = pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CASS", "CATS", "CARS"],
            "cdr3nt": ["TGTGCT", "TGTGCT", "TGTGCC", "TGTGCA"],
            "v": ["TRBV1", "TRBV1", "TRBV1", "TRBV2"],
            "j": ["TRBJ1", "TRBJ1", "TRBJ1", "TRBJ1"],
            "sample_id": ["sample_1"] * 4,
            "freq": [0.25] * 4,
            "count": [1, 2, 3, 4],
        },
        index=[9, 9, 9, 9],
    )
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(pooled)

    def copied_parallel_results(function, tasks, *args, **kwargs):
        return copy.deepcopy([function(task) for task in tasks])

    monkeypatch.setattr(
        clustering_module,
        "run_parallel_calculation",
        copied_parallel_results,
    )

    clusters.create_clusters(
        overlap_type="aaVJ", mismatches=1, cpu=2, verbose=False
    )

    plotted_nodes = [node for cluster in clusters for node in cluster]
    assert len(plotted_nodes) == len(pooled)
    assert sorted(node.id for node in plotted_nodes) == [0, 1, 2, 3]
    assert len({node.id for node in plotted_nodes}) == len(plotted_nodes)
    assert len({id(node) for node in plotted_nodes}) == len(plotted_nodes)
    assert sorted(len(cluster) for cluster in clusters) == [1, 3]

    trbv1_cluster_ids = {
        cluster.id
        for cluster in clusters
        for node in cluster
        if node.v == "TRBV1"
    }
    identical_cass_cluster_ids = {
        cluster.id
        for cluster in clusters
        for node in cluster
        if node.seq_aa == "CASS" and node.v == "TRBV1" and node.j == "TRBJ1"
    }
    assert len(trbv1_cluster_ids) == 1
    assert len(identical_cass_cluster_ids) == 1


def test_mutating_cluster_workflow_methods_do_not_return_self():
    clusters = Clusters()
    assert clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests()
    ) is None
    assert clusters.add_metadata(
        pd.DataFrame(
            {
                "sample_id": ["sample_1", "sample_2"],
                "group": ["control", "case"],
            }
        )
    ) is None
    assert clusters.create_clusters(
        overlap_type="aaVJ", mismatches=1, cpu=1, verbose=False
    ) is None
    assert clusters.filter(cluster_size >= 1, inplace=True) is None


def test_verbosity_false_suppresses_read_and_cached_clustering_output(capsys):
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests(), verbosity=False
    )
    clusters.create_clusters(
        overlap_type="aaVJ",
        mismatches=1,
        cpu=1,
        verbosity=False,
    )
    capsys.readouterr()

    clusters.create_clusters(
        overlap_type="aaVJ",
        mismatches=1,
        cpu=1,
        verbosity=False,
    )

    assert capsys.readouterr().out == ""


def test_verbosity_overrides_legacy_verbose_keyword(capsys):
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(
        _pooled_clonotypes_for_state_tests(),
        verbose=True,
        verbosity=False,
    )
    clusters.create_clusters(
        overlap_type="aaVJ",
        mismatches=1,
        cpu=1,
        verbose=True,
        verbosity=False,
    )

    assert capsys.readouterr().out == ""


@pytest.mark.parametrize(
    "constant_call",
    ["TRBC2", "TRAC", "TRGC1", "TRDC", "IGKC", "IGLC1"],
)
def test_non_igh_constant_segments_have_none_isotype(constant_call):
    pooled = pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CATS"],
            "cdr3nt": ["TGTGCT", "TGTGCC"],
            "v": ["TRBV1", "TRBV1"],
            "j": ["TRBJ1", "TRBJ1"],
            "c": [constant_call, constant_call],
            "sample_id": ["sample_1", "sample_1"],
            "freq": [0.5, 0.5],
            "count": [2, 1],
        }
    )
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(pooled, verbosity=False)

    clusters.create_clusters(
        overlap_type="aaV", mismatches=1, cpu=1, verbosity=False
    )

    nodes = [node for cluster in clusters for node in cluster]
    assert all(node.additional_properties["isotype"] is None for node in nodes)
    figure = clusters.plot_cluster(0, color="isotype", size=None)
    assert [text.get_text() for text in figure.legends[0].get_texts()] == ["NA"]


def test_mixed_igh_and_tcr_constants_recode_independently():
    pooled = pd.DataFrame(
        {
            "cdr3aa": ["CASS", "CATS"],
            "cdr3nt": ["TGTGCT", "TGTGCC"],
            "v": ["TRBV1", "TRBV2"],
            "j": ["TRBJ1", "TRBJ1"],
            "c": ["IGHG1*01", "TRBC2"],
            "sample_id": ["sample_1", "sample_1"],
            "freq": [0.5, 0.5],
            "count": [2, 1],
        }
    )
    clusters = Clusters()
    clusters.read_from_pooled_clonoset(pooled, verbosity=False)
    clusters.create_clusters(
        overlap_type="aaV", mismatches=1, cpu=1, verbosity=False
    )

    isotypes = sorted(
        (node.additional_properties["isotype"] for cluster in clusters for node in cluster),
        key=lambda value: "" if value is None else value,
    )
    assert isotypes == [None, "IgG1"]


def _clusters_for_clonoset_intersection():
    clusters = Clusters()
    specifications = [
        (10, [("CASS", "AAA"), ("CATS", "AAT")], "V1", "J1"),
        (20, [("GGGG", "GGG")], "V2", "J2"),
    ]
    for cluster_no, sequences, v_gene, j_gene in specifications:
        cluster = Cluster()
        cluster.id = cluster_no
        for node_index, (cdr3aa, cdr3nt) in enumerate(sequences):
            node = Node(
                f"{cluster_no}_{node_index}",
                cdr3nt,
                cdr3aa,
                v_gene,
                j_gene,
                "comparison",
                1 / len(sequences),
                1,
            )
            node.additional_properties["cluster_no"] = cluster_no
            cluster.add_node(node)
        clusters.clusters.append(cluster)
    return clusters


def _target_clonosets_for_cluster_intersection():
    return {
        "sample_1.tsv": pd.DataFrame(
            {
                "cdr3aa": ["CARS", "CASS", "GGGG"],
                "cdr3nt": ["AAC", "AAA", "GGG"],
                "v": ["V1", "V1", "V2"],
                "j": ["J1", "J1", "J2"],
                "count": [5, 3, 2],
                "freq": [0.5, 0.3, 0.2],
            }
        ),
        "sample_2.tsv": pd.DataFrame(
            {
                "cdr3aa": ["CASS", "AAAA"],
                "cdr3nt": ["AAA", "CCC"],
                "v": ["V1", "V3"],
                "j": ["J1", "J3"],
                "count": [4, 6],
                "freq": [0.4, 0.6],
            }
        ),
    }


def test_intersect_with_clonosets_counts_unique_target_clonotypes(
    monkeypatch,
):
    clusters = _clusters_for_clonoset_intersection()
    targets = _target_clonosets_for_cluster_intersection()
    monkeypatch.setattr(
        clustering_module,
        "read_clonoset",
        lambda filename: targets[filename].copy(),
    )
    samples = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "filename": ["sample_1.tsv", "sample_2.tsv"],
        }
    )

    table = clusters.intersect_with_clonosets(
        samples, overlap_type="aaVJ", mismatches=1, cpu=1
    )

    assert table.columns.tolist() == [
        "cluster_id",
        "consensus",
        "concensus_cdr3aa",
        "concensus_v",
        "concensus_j",
        "sample_1",
        "sample_2",
    ]
    assert table["cluster_id"].tolist() == ["cluster_10", "cluster_20"]
    np.testing.assert_allclose(table["sample_1"], [8, 2])
    np.testing.assert_allclose(table["sample_2"], [4, 0])


def test_intersect_with_clonosets_calculates_directional_frequencies(
    monkeypatch,
):
    clusters = _clusters_for_clonoset_intersection()
    targets = _target_clonosets_for_cluster_intersection()
    monkeypatch.setattr(
        clustering_module,
        "read_clonoset",
        lambda filename: targets[filename].copy(),
    )
    samples = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "filename": ["sample_1.tsv", "sample_2.tsv"],
        }
    )

    table = clusters.intersect_with_clonosets(
        samples,
        overlap_type="aaVJ",
        mismatches=1,
        by_freq=True,
        cpu=1,
    )

    np.testing.assert_allclose(table["sample_1"], [0.8, 0.2])
    np.testing.assert_allclose(table["sample_2"], [0.4, 0])


def test_intersect_with_clonosets_respects_mismatch_threshold(monkeypatch):
    clusters = _clusters_for_clonoset_intersection()
    targets = _target_clonosets_for_cluster_intersection()
    monkeypatch.setattr(
        clustering_module,
        "read_clonoset",
        lambda filename: targets[filename].copy(),
    )
    samples = pd.DataFrame(
        {"sample_id": ["sample_1"], "filename": ["sample_1.tsv"]}
    )

    table = clusters.intersect_with_clonosets(
        samples, overlap_type="aaVJ", mismatches=0, cpu=1
    )

    np.testing.assert_allclose(table["sample_1"], [3, 2])


def test_intersect_with_clonosets_forwards_cpu(monkeypatch):
    clusters = _clusters_for_clonoset_intersection()
    targets = _target_clonosets_for_cluster_intersection()
    monkeypatch.setattr(
        clustering_module,
        "read_clonoset",
        lambda filename: targets[filename].copy(),
    )
    captured = {}

    def fake_parallel(function, tasks, *args, **kwargs):
        captured["cpu"] = kwargs["cpu"]
        return [function(task) for task in tasks]

    monkeypatch.setattr(
        clustering_module, "run_parallel_calculation", fake_parallel
    )
    samples = pd.DataFrame(
        {"sample_id": ["sample_1"], "filename": ["sample_1.tsv"]}
    )

    clusters.intersect_with_clonosets(samples, cpu=3)

    assert captured["cpu"] == 3


def _renumber_clusters(clusters, cluster_numbers):
    for cluster, cluster_no in zip(clusters, cluster_numbers):
        cluster.id = cluster_no
        for node in cluster:
            node.additional_properties["cluster_no"] = cluster_no
    return clusters


def test_plot_cluster_uses_persistent_ids_after_filtering_and_scales_facets():
    clusters = _renumber_clusters(_clusters_with_two_samples(), [10, 11])
    filtered = clusters.select([10, 11])

    figure = filtered.plot_cluster(
        ["cluster_10", 11],
        ncols=2,
        height=3,
        aspect=1.5,
        size=None,
    )

    assert [axis.get_title() for axis in figure.axes] == [
        "cluster_10",
        "cluster_11",
    ]
    np.testing.assert_allclose(figure.get_size_inches(), [9, 3])
    with pytest.raises(KeyError, match="cluster_0"):
        filtered.plot_cluster(0)


def test_plot_logo_accepts_cluster_ids_and_mixed_lists(monkeypatch):
    clusters = _renumber_clusters(_clusters_with_two_samples(), [10, 11])
    calls = []

    def fake_logo(clonotypes, seq_type, plot=True):
        calls.append((clonotypes, seq_type, plot))
        return f"logo_{len(calls)}"

    monkeypatch.setattr(
        clustering_module, "get_logo_for_list_of_clonotypes", fake_logo
    )

    single = clusters.plot_logo("cluster_10", plot=False)
    multiple = clusters.plot_logo([10, "cluster_11"], plot=False)

    assert single == "logo_1"
    assert multiple == {"cluster_10": "logo_2", "cluster_11": "logo_3"}
    assert len(calls[0][0]) == 3
    assert len(calls[1][0]) == 3
    assert len(calls[2][0]) == 1


@pytest.mark.parametrize(("height", "aspect"), [(0, 1), (3, 0), (-1, 1), (3, -1)])
def test_plot_cluster_validates_height_and_aspect(height, aspect):
    with pytest.raises(ValueError, match="positive number"):
        _clusters_with_two_samples().plot_cluster(
            0, height=height, aspect=aspect
        )


def test_properties_forwards_cpu_and_uses_parallel_runner(monkeypatch):
    clusters = _clusters_with_two_samples()
    captured = {}

    def fake_parallel(function, tasks, *args, **kwargs):
        captured["cpu"] = kwargs["cpu"]
        captured["task_count"] = len(tasks)
        return [function(task) for task in tasks]

    monkeypatch.setattr(
        clustering_module, "run_parallel_calculation", fake_parallel
    )

    properties = clusters.properties(cpu=3)

    assert captured == {"cpu": 3, "task_count": 2}
    assert properties["cluster_id"].tolist() == ["cluster_0", "cluster_1"]


def test_properties_singleton_shortcut_skips_all_consensus_calculations(
    monkeypatch,
):
    clusters = _clusters_with_two_samples()
    original_consensus = Cluster.calc_cluster_consensus
    original_segment = Cluster.calc_cluster_consensus_segment

    def guarded_consensus(cluster, *args, **kwargs):
        if len(cluster) == 1:
            raise AssertionError("singleton sequence consensus was calculated")
        return original_consensus(cluster, *args, **kwargs)

    def guarded_segment(cluster, *args, **kwargs):
        if len(cluster) == 1:
            raise AssertionError("singleton segment consensus was calculated")
        return original_segment(cluster, *args, **kwargs)

    monkeypatch.setattr(Cluster, "calc_cluster_consensus", guarded_consensus)
    monkeypatch.setattr(Cluster, "calc_cluster_consensus_segment", guarded_segment)

    properties = clusters.properties(cpu=1)
    singleton = properties.loc[properties["nodes"] == 1].iloc[0]

    assert singleton["concensus_cdr3aa"] == "CATS"
    assert singleton["concensus_cdr3nt"] == "TGTGCC"
    assert singleton["concensus_v"] == "TRBV2"
    assert singleton["concensus_j"] == "TRBJ2"


def test_properties_uses_first_v_and_j_when_overlap_requires_them(monkeypatch):
    clusters = _clusters_with_two_samples()
    clusters.overlap_type = "aaVJ"

    def fail_segment_consensus(cluster, *args, **kwargs):
        raise AssertionError("V/J modal consensus should have been skipped")

    monkeypatch.setattr(
        Cluster, "calc_cluster_consensus_segment", fail_segment_consensus
    )

    properties = clusters.properties(cpu=1)

    assert properties["concensus_v"].tolist() == ["TRBV1", "TRBV2"]
    assert properties["concensus_j"].tolist() == ["TRBJ1", "TRBJ2"]


def test_properties_prints_progress_on_first_calculation(capsys):
    clusters = _clusters_with_two_samples()

    clusters.properties(cpu=1)

    output = capsys.readouterr().out
    assert "Calculating cluster properties" in output
    assert "2/2 clusters processed" in output
