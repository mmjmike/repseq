import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from repseq import plot as rsplot


def _stats_df():
    return pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "TRA",
                "diversity": 10,
                "norm_shannon_wiener": 0.8,
                "chao1": 12,
                "berger_parker": 0.2,
                "d50": 0.1,
                "convergence": 1.1,
                "convergence_v": 1.2,
                "convergence_vj": 1.3,
            },
            {
                "sample_id": "sample2",
                "chain": "TRB",
                "diversity": 20,
                "norm_shannon_wiener": 0.7,
                "chao1": 22,
                "berger_parker": 0.3,
                "d50": 0.2,
                "convergence": 1.4,
                "convergence_v": 1.5,
                "convergence_vj": 1.6,
            },
        ]
    )


def test_plot_module_import_and_diversity_subset_warning():
    metadata = pd.DataFrame(
        [{"sample_id": "sample1", "chain": "TRA", "condition": "treated"}]
    )

    with pytest.warns(UserWarning, match="1 of 2 samples"):
        grid = rsplot.diversity_stats(
            _stats_df(),
            metadata=metadata,
            properties=["diversity"],
            group="condition",
        )

    assert grid.fig is not None


def test_plot_stats_rejects_three_group_columns():
    metadata = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "TRA",
                "condition": "treated",
                "batch": "b1",
                "donor": "d1",
            }
        ]
    )

    with pytest.raises(ValueError, match="group can contain at most 2"):
        rsplot.plot_stats(
            _stats_df().iloc[:1],
            metadata=metadata,
            properties=["diversity"],
            group=["condition", "batch", "donor"],
        )


def test_plot_stats_checks_group_columns_are_in_metadata():
    metadata = pd.DataFrame([{"sample_id": "sample1", "chain": "TRA"}])

    with pytest.raises(ValueError, match="group column"):
        rsplot.plot_stats(
            _stats_df().iloc[:1],
            metadata=metadata,
            properties=["diversity"],
            group="condition",
        )


def test_plot_stats_preserves_ordered_group_categories():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "condition": "treated"},
            {"sample_id": "sample2", "chain": "TRB", "condition": "control"},
        ]
    )
    metadata["condition"] = pd.Categorical(
        metadata["condition"],
        categories=["control", "treated"],
        ordered=True,
    )

    grid = rsplot.plot_stats(
        _stats_df(),
        metadata=metadata,
        properties=["diversity"],
        group="condition",
    )

    tick_labels = [tick.get_text() for tick in grid.axes.flat[0].get_xticklabels()]
    assert tick_labels == ["control", "treated"]


def test_convergence_defaults_create_property_facets():
    grid = rsplot.convergence(_stats_df())

    titles = [ax.get_title() for ax in grid.axes.flat]
    assert titles == ["Convergence", "Convergence V", "Convergence Vj"]


def test_two_split_columns_create_interaction_columns():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "condition": "treated", "batch": "b1"},
            {"sample_id": "sample2", "chain": "TRB", "condition": "control", "batch": "b2"},
        ]
    )

    grid = rsplot.diversity_stats(
        _stats_df(),
        metadata=metadata,
        properties=["diversity", "chao1"],
        split=["condition", "batch"],
    )

    titles = [ax.get_title() for ax in grid.axes.flat]
    assert any("treated | b1" in title for title in titles)
    assert any("control | b2" in title for title in titles)
