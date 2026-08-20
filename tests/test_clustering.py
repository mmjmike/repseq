import numpy as np
import pandas as pd
import pytest

import repseq.clustering as clustering_module
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
