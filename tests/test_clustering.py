import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgba

import repseq.clustering as clustering_module
import repseq.plot as rsplot
from repseq.clustering import Cluster, Clusters, Node


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


def test_clusters_str_reports_cluster_node_and_singleton_counts():
    clusters = _clusters_with_two_samples()
    clusters.clonotypes = clusters.clonosets_df

    summary = str(clusters)

    assert (
        "Clusters from 3 samples with 2 clusters and 4 nodes, "
        "of which 1 is a single node."
    ) in summary


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
        log_power=2,
    )

    assert [axis.get_title() for axis in figure.axes] == ["Cluster 0", "Cluster 1"]
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
        100 + np.log2(count + 1) ** 2 for count in [3, 7, 5, 11]
    ]
    np.testing.assert_allclose(sorted(plotted_sizes), sorted(expected_sizes))
    assert figure.number not in plt.get_fignums()


@pytest.mark.parametrize(
    "layout", ["spring", "kamada_kawai", "circular", "shell", "spectral"]
)
def test_plot_cluster_supports_common_layouts(layout):
    figure = _clusters_with_two_samples().plot_cluster(0, layout=layout)

    assert figure.axes[0].get_title() == "Cluster 0"
    plt.close(figure)


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
