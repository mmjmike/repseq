import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba

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


def test_single_group_does_not_pass_new_seaborn_legend_argument(monkeypatch):
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "condition": "treated"},
            {"sample_id": "sample2", "chain": "TRB", "condition": "control"},
        ]
    )
    original_boxplot = rsplot.sns.boxplot
    original_stripplot = rsplot.sns.stripplot

    def compatible_boxplot(*args, **kwargs):
        assert "legend" not in kwargs
        return original_boxplot(*args, **kwargs)

    def compatible_stripplot(*args, **kwargs):
        assert "legend" not in kwargs
        return original_stripplot(*args, **kwargs)

    monkeypatch.setattr(rsplot.sns, "boxplot", compatible_boxplot)
    monkeypatch.setattr(rsplot.sns, "stripplot", compatible_stripplot)

    grid = rsplot.plot_stats(
        _stats_df(),
        metadata=metadata,
        properties=["diversity"],
        group="condition",
    )

    assert grid.axes.flat[0].legend_ is None


def test_split_panels_only_show_samples_from_their_subset():
    stats_df = pd.concat(
        [
            _stats_df(),
            pd.DataFrame(
                [
                    {
                        "sample_id": "sample3",
                        "chain": "TRA",
                        "diversity": 30,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "condition": "treated"},
            {"sample_id": "sample2", "chain": "TRB", "condition": "control"},
            {"sample_id": "sample3", "chain": "TRA", "condition": "treated"},
        ]
    )

    grid = rsplot.plot_stats(
        stats_df,
        metadata=metadata,
        properties=["diversity"],
        split="condition",
    )

    labels_by_panel = {
        ax.get_title(): [tick.get_text() for tick in ax.get_xticklabels()]
        for ax in grid.axes.flat
    }
    assert labels_by_panel["treated"] == ["sample1", "sample3"]
    assert labels_by_panel["control"] == ["sample2"]


def _clonoset_stats_df():
    return pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "TRB",
                "reads": 150,
                "reads_func": 120,
                "reads_per_umi": 3.0,
                "clones": 15,
                "clones_func": 12,
                "umi": 50,
                "umi_func": 40,
            },
            {
                "sample_id": "sample2",
                "chain": "TRB",
                "reads": 100,
                "reads_func": 90,
                "reads_per_umi": np.nan,
                "clones": 10,
                "clones_func": 8,
                "umi": np.nan,
                "umi_func": np.nan,
            },
        ]
    )


def test_clonoset_stats_default_overlays_total_and_functional_bars():
    grid = rsplot.clonoset_stats(_clonoset_stats_df())
    axes = {ax.get_title(): ax for ax in grid.axes.flat}

    assert set(axes) == {"Reads", "Reads Per Umi", "Clones", "Umi"}
    read_patches = axes["Reads"].patches
    np.testing.assert_allclose(
        [patch.get_height() for patch in read_patches],
        [150, 100, 120, 90],
    )
    np.testing.assert_allclose([patch.get_width() for patch in read_patches], 0.8)
    total_color = read_patches[0].get_facecolor()[:3]
    functional_color = read_patches[2].get_facecolor()[:3]
    assert sum(total_color) > sum(functional_color)
    assert total_color[2] > total_color[0]
    assert functional_color[2] > functional_color[0]
    assert [text.get_text() for text in axes["Reads"].texts] == [
        "150(120)",
        "100(90)",
    ]
    assert [tick.get_text() for tick in axes["Umi"].get_xticklabels()] == [
        "sample1"
    ]
    assert [text.get_text() for text in axes["Umi"].texts] == ["50(40)"]


def test_clonoset_stats_skips_umi_panel_when_all_values_are_missing():
    stats_df = _clonoset_stats_df()
    stats_df[["umi", "umi_func", "reads_per_umi"]] = np.nan

    grid = rsplot.clonoset_stats(stats_df)

    assert {ax.get_title() for ax in grid.axes.flat} == {
        "Reads",
        "Reads Per Umi",
        "Clones",
    }


def test_clonoset_stats_grouped_defaults_follow_available_columns(monkeypatch):
    captured = []

    def capture_plot_stats(stats_df, **kwargs):
        captured.append(kwargs)
        return "grid"

    monkeypatch.setattr(rsplot, "plot_stats", capture_plot_stats)
    stats_df = _clonoset_stats_df()

    assert rsplot.clonoset_stats(stats_df, group="condition") == "grid"
    assert captured[-1]["properties"] == [
        "reads",
        "reads_per_umi",
        "clones_func",
        "umi_func",
    ]

    stats_df["umi_func"] = np.nan
    assert rsplot.clonoset_stats(stats_df, group="condition") == "grid"
    assert captured[-1]["properties"] == [
        "reads",
        "reads_per_umi",
        "clones_func",
    ]


def test_clonoset_stats_custom_properties_use_standard_plot(monkeypatch):
    captured = {}

    def capture_plot_stats(stats_df, **kwargs):
        captured.update(kwargs)
        return "grid"

    monkeypatch.setattr(rsplot, "plot_stats", capture_plot_stats)

    assert (
        rsplot.clonoset_stats(_clonoset_stats_df(), properties=["reads_func"])
        == "grid"
    )
    assert captured["properties"] == ["reads_func"]
    assert captured["group"] is None


def test_clonoset_stats_default_supports_two_ordered_splits():
    metadata = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "TRB",
                "tissue": "blood",
                "sex": "M",
            },
            {
                "sample_id": "sample2",
                "chain": "TRB",
                "tissue": "tumor",
                "sex": "F",
            },
        ]
    )
    metadata["tissue"] = pd.Categorical(
        metadata["tissue"], categories=["tumor", "blood"], ordered=True
    )
    metadata["sex"] = pd.Categorical(
        metadata["sex"], categories=["F", "M"], ordered=True
    )

    grid = rsplot.clonoset_stats(
        _clonoset_stats_df(),
        metadata=metadata,
        split=["tissue", "sex"],
    )

    column_titles = [ax.get_title() for ax in grid.axes[0]]
    assert column_titles == ["Reads | tumor | F", "Reads | blood | M"]


def test_clonoset_stats_grouped_boxplots_start_at_zero():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "condition": "control"},
            {"sample_id": "sample2", "chain": "TRB", "condition": "control"},
        ]
    )

    grid = rsplot.clonoset_stats(
        _clonoset_stats_df(),
        metadata=metadata,
        group="condition",
    )

    assert all(ax.get_ylim()[0] == 0 for ax in grid.axes.flat)


def test_diversity_stats_always_start_at_zero():
    grid = rsplot.diversity_stats(_stats_df(), properties=["diversity"])

    assert grid.axes.flat[0].get_ylim()[0] == 0


def test_other_stats_boxplots_keep_automatic_y_limits():
    stats_df = pd.DataFrame(
        [
            {"sample_id": "sample1", "metric": 10},
            {"sample_id": "sample2", "metric": 20},
        ]
    )
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "condition": "control"},
            {"sample_id": "sample2", "condition": "control"},
        ]
    )

    grid = rsplot.cdr3aa_stats(
        stats_df,
        metadata=metadata,
        properties=["metric"],
        group="condition",
    )

    assert grid.axes.flat[0].get_ylim()[0] > 0

def _v_usage_long(sample_count=4):
    rows = []
    for sample_number in range(sample_count):
        for segment, usage in [("TRBV10", 0.6), ("TRBV2", 0.4)]:
            rows.append(
                {
                    "sample_id": f"sample{sample_number}",
                    "chain": "TRB",
                    "v": segment,
                    "usage": usage + sample_number * 0.01,
                }
            )
    return pd.DataFrame(rows)


def test_parse_gene_name_handles_isotypes_dual_genes_and_natural_sorting():
    assert rsplot.parse_gene_name("IGHD")["gene_type"] == "C"
    assert rsplot.parse_gene_name("IGHD3-10")["gene_type"] == "D"
    assert rsplot.parse_gene_name("TRAV8-2DV6")["dual_designation"] == "DV6"
    assert rsplot.parse_gene_name("TRAV8-2/DV6")["dual_designation"] == "DV6"
    assert rsplot._segment_family("TRBJ2-1", "j") == "TRBJ2"
    assert rsplot._segment_family("IGHG2", "c") == "IGHG"
    assert rsplot._sort_gene_names(
        ["TRBV12-3-2", "TRBV7-8", "TRBV2", "TRBV12-3-1", "TRBV7-3"]
    ) == ["TRBV2", "TRBV7-3", "TRBV7-8", "TRBV12-3-1", "TRBV12-3-2"]


def test_segment_usage_detects_wide_table_and_keeps_chain_axes_independent():
    usage = pd.DataFrame(
        [
            {
                "sample_id": "sample_TRA",
                "chain": "TRA",
                "TRAV10": 0.7,
                "TRAV2": 0.3,
                "TRBV2": 0,
            },
            {
                "sample_id": "sample_TRB",
                "chain": "TRB",
                "TRAV10": 0,
                "TRAV2": 0,
                "TRBV2": 1.0,
            },
        ]
    )

    fig = rsplot.segment_usage(usage)

    heatmaps = {ax.get_title(): ax for ax in fig.axes if ax.get_title()}
    assert set(heatmaps) == {"TRA", "TRB"}
    assert heatmaps["TRA"].collections[0].cmap.name == "pheatmap_default"
    assert fig.number not in plt.get_fignums()
    assert [tick.get_text() for tick in heatmaps["TRA"].get_xticklabels()] == [
        "TRAV2",
        "TRAV10",
    ]
    assert [tick.get_text() for tick in heatmaps["TRB"].get_xticklabels()] == [
        "TRBV2"
    ]


def test_segment_usage_combines_v_families_within_each_sample():
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "v": "TRBV7-2", "usage": 0.2},
            {"sample_id": "sample1", "chain": "TRB", "v": "TRBV7-3", "usage": 0.3},
            {"sample_id": "sample1", "chain": "TRB", "v": "TRBV10-1", "usage": 0.5},
            {"sample_id": "sample2", "chain": "TRB", "v": "TRBV7-2", "usage": 0.1},
            {"sample_id": "sample2", "chain": "TRB", "v": "TRBV7-3", "usage": 0.6},
            {"sample_id": "sample2", "chain": "TRB", "v": "TRBV10-1", "usage": 0.3},
        ]
    )

    fig = rsplot.segment_usage(usage, combine_families=True)
    heatmap = next(ax for ax in fig.axes if ax.get_title() == "TRB")

    assert [tick.get_text() for tick in heatmap.get_xticklabels()] == [
        "TRBV7",
        "TRBV10",
    ]
    assert heatmap.get_xlabel() == "Segment Family"
    np.testing.assert_allclose(
        np.asarray(heatmap.collections[0].get_array()).reshape(2, 2),
        [[0.5, 0.5], [0.7, 0.3]],
    )


def test_segment_usage_combines_c_segments_by_isotype():
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "IGH", "c": "IGHG1", "usage": 0.2},
            {"sample_id": "sample1", "chain": "IGH", "c": "IGHG2", "usage": 0.3},
            {"sample_id": "sample1", "chain": "IGH", "c": "IGHA1", "usage": 0.4},
            {"sample_id": "sample1", "chain": "IGH", "c": ".", "usage": 0.1},
        ]
    )

    fig = rsplot.segment_usage(
        usage,
        plot_type="barplot",
        combine_families=True,
    )
    ax = fig.axes[0]

    assert [tick.get_text() for tick in ax.get_xticklabels()] == [
        "IGHA",
        "IGHG",
        "NA",
    ]
    assert ax.get_xlabel() == "Segment Family"
    np.testing.assert_allclose(
        [patch.get_height() for patch in ax.patches],
        [0.4, 0.5, 0.1],
    )


def test_segment_usage_renames_dot_segment_in_wide_tables():
    usage = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "chain": "IGH",
                "IGHG1": 0.8,
                ".": 0.2,
            }
        ]
    )

    fig = rsplot.segment_usage(usage)
    heatmap = next(ax for ax in fig.axes if ax.get_title() == "IGH")

    assert [tick.get_text() for tick in heatmap.get_xticklabels()] == [
        "IGHG1",
        "NA",
    ]
    np.testing.assert_allclose(
        np.asarray(heatmap.collections[0].get_array()),
        [0.8, 0.2],
    )


def test_segment_usage_combines_families_in_boxplots():
    usage = _v_usage_long()
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "control" if sample_number < 2 else "treated",
            }
            for sample_number in range(4)
        ]
    )

    fig = rsplot.segment_usage(
        usage,
        metadata=metadata,
        plot_type="boxplot",
        group="condition",
        combine_families=True,
    )
    ax = fig.axes[0]

    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["TRBV2", "TRBV10"]
    assert ax.get_xlabel() == "Segment Family"


def test_segment_usage_grouped_barplot_uses_standard_deviation(monkeypatch):
    usage = _v_usage_long()
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "control" if sample_number < 2 else "treated",
            }
            for sample_number in range(4)
        ]
    )
    error_bars = []
    original_bar = Axes.bar

    def capture_bar(self, *args, **kwargs):
        error_bars.append(kwargs.get("yerr"))
        return original_bar(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)

    fig = rsplot.segment_usage(
        usage,
        metadata=metadata,
        plot_type="barplot",
        group="condition",
    )

    assert fig is not None
    assert fig.number not in plt.get_fignums()
    assert len(error_bars) == 2
    assert all(np.all(np.asarray(values) > 0) for values in error_bars)
    assert [tick.get_text() for tick in fig.axes[0].get_xticklabels()] == [
        "TRBV2",
        "TRBV10",
    ]


def test_segment_usage_ungrouped_boxplot_stops_above_ten_samples():
    with pytest.warns(UserWarning, match="at most 10 groups"):
        fig = rsplot.segment_usage(
            _v_usage_long(sample_count=11),
            plot_type="boxplot",
        )

    assert fig is None


def test_segment_usage_heatmap_supports_three_ordered_annotations():
    usage = _v_usage_long()
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "treated" if sample_number >= 2 else "control",
                "batch": f"b{sample_number % 2 + 1}",
                "sex": "F" if sample_number % 2 else "M",
            }
            for sample_number in range(4)
        ]
    )
    metadata["condition"] = pd.Categorical(
        metadata["condition"],
        categories=["control", "treated"],
        ordered=True,
    )

    fig = rsplot.segment_usage(
        usage,
        metadata=metadata,
        group=["condition", "batch", "sex"],
    )

    annotation_axes = [ax for ax in fig.axes if len(ax.images) == 1]
    assert len(annotation_axes) == 1
    assert [tick.get_text() for tick in annotation_axes[0].get_xticklabels()] == [
        "Condition",
        "Batch",
        "Sex",
    ]
    legend_labels = [text.get_text() for text in fig.legends[0].get_texts()]
    assert "control" in legend_labels
    assert "treated" in legend_labels
    assert all(":" not in label for label in legend_labels)


def test_segment_usage_split_order_is_preserved_in_rows():
    usage = _v_usage_long()
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "tissue": "blood" if sample_number < 2 else "tumor",
            }
            for sample_number in range(4)
        ]
    )
    metadata["tissue"] = pd.Categorical(
        metadata["tissue"], categories=["tumor", "blood"], ordered=True
    )

    fig = rsplot.segment_usage(
        usage,
        metadata=metadata,
        plot_type="barplot",
        split="tissue",
    )

    assert [ax.get_title() for ax in fig.axes] == ["TRB | tumor", "TRB | blood"]


def test_segment_usage_rejects_too_many_boxplot_groups():
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "control",
                "batch": "b1",
            }
            for sample_number in range(4)
        ]
    )

    with pytest.raises(ValueError, match="group can contain at most 1"):
        rsplot.segment_usage(
            _v_usage_long(),
            metadata=metadata,
            plot_type="boxplot",
            group=["condition", "batch"],
        )


def _cdr3_lengths_long(sample_count=2):
    rows = []
    for sample_number in range(sample_count):
        rows.extend(
            [
                {
                    "sample_id": f"sample{sample_number}",
                    "chain": "TRB",
                    "cdr3_length": 3,
                    "freq": 0.2 + sample_number * 0.2,
                },
                {
                    "sample_id": f"sample{sample_number}",
                    "chain": "TRB",
                    "cdr3_length": 4,
                    "freq": 0.8 - sample_number * 0.2,
                },
            ]
        )
    return pd.DataFrame(rows)


def test_cdr3_length_distributions_plots_long_sample_frequencies():
    fig = rsplot.cdr3_length_distributions(_cdr3_lengths_long())
    ax = fig.axes[0]

    assert ax.get_title() == "TRB"
    assert ax.get_xlabel() == "CDR3 Length"
    assert ax.get_ylabel() == "Frequency"
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["3", "4"]
    np.testing.assert_allclose(
        [patch.get_height() for patch in ax.patches],
        [0.2, 0.8, 0.4, 0.6],
    )
    assert fig.number not in plt.get_fignums()


def test_cdr3_length_distributions_group_means_have_no_error_bars(monkeypatch):
    lengths = _cdr3_lengths_long(sample_count=4)
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "control" if sample_number < 2 else "treated",
                "tissue": "blood",
                "sex": "F",
            }
            for sample_number in range(4)
        ]
    )
    bar_kwargs = []
    original_bar = Axes.bar

    def capture_bar(self, *args, **kwargs):
        bar_kwargs.append(kwargs)
        return original_bar(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)
    fig = rsplot.cdr3_length_distributions(
        lengths,
        metadata=metadata,
        group="condition",
        split=["tissue", "sex"],
    )

    assert [ax.get_title() for ax in fig.axes] == ["TRB | blood | F"]
    assert all("yerr" not in kwargs for kwargs in bar_kwargs)
    np.testing.assert_allclose(
        [patch.get_height() for patch in fig.axes[0].patches],
        [0.3, 0.7, 0.7, 0.3],
    )


def test_cdr3_length_distributions_accepts_wide_count_table():
    lengths = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", 3: 2, 4: 8},
            {"sample_id": "sample2", "chain": "TRB", 3: 4, 4: 6},
        ]
    )

    fig = rsplot.cdr3_length_distributions(lengths)

    assert fig.axes[0].get_ylabel() == "Count"
    np.testing.assert_allclose(
        [patch.get_height() for patch in fig.axes[0].patches],
        [2, 8, 4, 6],
    )


def test_cdr3_length_distributions_preserves_ordered_two_split_panels():
    lengths = _cdr3_lengths_long()
    metadata = pd.DataFrame(
        [
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "tissue": "blood",
                "sex": "M",
            },
            {
                "sample_id": "sample1",
                "chain": "TRB",
                "tissue": "tumor",
                "sex": "F",
            },
        ]
    )
    metadata["tissue"] = pd.Categorical(
        metadata["tissue"], categories=["tumor", "blood"], ordered=True
    )
    metadata["sex"] = pd.Categorical(
        metadata["sex"], categories=["F", "M"], ordered=True
    )

    fig = rsplot.cdr3_length_distributions(
        lengths,
        metadata=metadata,
        split=["tissue", "sex"],
    )

    assert [ax.get_title() for ax in fig.axes] == [
        "TRB | tumor | F",
        "TRB | blood | M",
    ]


def test_cdr3_length_distributions_enforces_group_and_series_limits():
    metadata = pd.DataFrame(
        [
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "condition": "control",
                "batch": "b1",
            }
        ]
    )
    with pytest.raises(ValueError, match="group can contain at most 1"):
        rsplot.cdr3_length_distributions(
            _cdr3_lengths_long(sample_count=1),
            metadata=metadata,
            group=["condition", "batch"],
        )

    with pytest.warns(UserWarning, match="at most 10 groups"):
        fig = rsplot.cdr3_length_distributions(_cdr3_lengths_long(sample_count=11))
    assert fig is None

def _vj_usage_long(sample_count=2):
    rows = []
    for sample_number in range(sample_count):
        rows.extend(
            [
                {
                    "sample_id": f"sample{sample_number}",
                    "chain": "TRB",
                    "v": "TRBV10",
                    "j": "TRBJ2-1",
                    "vj": "TRBV10|TRBJ2-1",
                    "usage": 0.2 + sample_number * 0.05,
                },
                {
                    "sample_id": f"sample{sample_number}",
                    "chain": "TRB",
                    "v": "TRBV2",
                    "j": "TRBJ1-1",
                    "vj": "TRBV2|TRBJ1-1",
                    "usage": 0.6 - sample_number * 0.05,
                },
            ]
        )
    return pd.DataFrame(rows)


def _vjlen_usage_long(sample_ids=("sample0", "sample1")):
    rows = []
    for sample_number, sample_id in enumerate(sample_ids):
        rows.extend(
            [
                {
                    "sample_id": sample_id,
                    "chain": "TRB",
                    "v": "TRBV2",
                    "j": "TRBJ1-1",
                    "len": 15,
                    "vjlen": "TRBV2|TRBJ1-1|15",
                    "usage": 0.2 + sample_number * 0.2,
                },
                {
                    "sample_id": sample_id,
                    "chain": "TRB",
                    "v": "TRBV10",
                    "j": "TRBJ2-1",
                    "len": 14,
                    "vjlen": "TRBV10|TRBJ2-1|14",
                    "usage": 0.1 if sample_number == 0 else 0,
                },
            ]
        )
    return pd.DataFrame(rows)


def test_vj_usage_draws_largest_black_edged_dots_first_and_repels_series():
    fig = rsplot.vj_usage(_vj_usage_long())

    collection = fig.axes[0].collections[0]
    sizes = collection.get_sizes()
    offsets = collection.get_offsets()
    assert list(sizes) == sorted(sizes, reverse=True)
    assert np.allclose(collection.get_edgecolors()[0], to_rgba("black"))
    assert any(not np.isclose(value, round(value)) for value in offsets.ravel())
    assert fig.number not in plt.get_fignums()


def test_vj_usage_keeps_single_dot_at_its_vj_center():
    usage = _vj_usage_long()
    usage.loc[
        (usage["sample_id"] == "sample1")
        & (usage["vj"] == "TRBV10|TRBJ2-1"),
        "usage",
    ] = 0

    fig = rsplot.vj_usage(usage)

    offsets = np.asarray(fig.axes[0].collections[0].get_offsets())
    assert any(np.allclose(point, [1, 1]) for point in offsets)


def test_vj_usage_grouped_values_are_means_and_limit_is_eight():
    usage = _vj_usage_long(sample_count=4)
    metadata = pd.DataFrame(
        [
            {
                "sample_id": f"sample{sample_number}",
                "chain": "TRB",
                "condition": "control" if sample_number < 2 else "treated",
            }
            for sample_number in range(4)
        ]
    )
    data, groups, _ = rsplot._prepare_combination_plot_data(
        usage, metadata, "condition", None, "vj"
    )
    means = rsplot._aggregate_vj_panel(data, groups[0], grouped=True)
    control_v2 = means.loc[
        (means["_v"] == "TRBV2") & (means["condition"] == "control"),
        "_value",
    ].iloc[0]
    assert control_v2 == pytest.approx(0.575)

    with pytest.raises(ValueError, match="at most 8 samples"):
        rsplot.vj_usage(_vj_usage_long(sample_count=9))


def test_vj_usage_accepts_pipe_delimited_wide_columns():
    usage = pd.DataFrame(
        [
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "TRBV2|TRBJ1-1": 0.7,
                "TRBV10|TRBJ2-1": 0.3,
            }
        ]
    )

    fig = rsplot.vj_usage(usage)

    assert [tick.get_text() for tick in fig.axes[0].get_xticklabels()] == [
        "TRBV2",
        "TRBV10",
    ]


def test_vj_usage_combines_v_and_j_families():
    usage = pd.DataFrame(
        [
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "v": "TRBV7-2",
                "j": "TRBJ2-1",
                "vj": "TRBV7-2|TRBJ2-1",
                "usage": 0.2,
            },
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "v": "TRBV7-3",
                "j": "TRBJ2-2",
                "vj": "TRBV7-3|TRBJ2-2",
                "usage": 0.3,
            },
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "v": "TRBV10-1",
                "j": "TRBJ1-1",
                "vj": "TRBV10-1|TRBJ1-1",
                "usage": 0.5,
            },
        ]
    )

    legacy_usage = usage.drop(columns=["v", "j"])
    data, _, _ = rsplot._prepare_combination_plot_data(
        legacy_usage,
        None,
        None,
        None,
        "vj",
        combine_families=True,
    )
    combined = data.loc[(data["_v"] == "TRBV7") & (data["_j"] == "TRBJ2")]
    assert combined["_value"].iloc[0] == pytest.approx(0.5)

    fig = rsplot.vj_usage(legacy_usage, combine_families=True)
    ax = fig.axes[0]
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["TRBV7", "TRBV10"]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["TRBJ1", "TRBJ2"]
    assert ax.get_xlabel() == "V segment family"
    assert ax.get_ylabel() == "J segment family"
    assert len(ax.collections[0].get_offsets()) == 2


def test_vj_usage_accepts_component_only_long_tables():
    usage = _vj_usage_long().drop(columns="vj")

    fig = rsplot.vj_usage(usage)

    assert [tick.get_text() for tick in fig.axes[0].get_xticklabels()] == [
        "TRBV2",
        "TRBV10",
    ]


def test_vj_usage_validates_optional_component_columns():
    usage = _vj_usage_long()
    usage.loc[0, "v"] = "TRBV3"

    with pytest.raises(ValueError, match="must match the separate v, j"):
        rsplot.vj_usage(usage)


def test_vj_usage_rejects_tuple_columns_and_tuple_identifiers():
    tuple_wide = pd.DataFrame(
        [{"sample_id": "sample0", "chain": "TRB", ("TRBV2", "TRBJ1-1"): 1.0}]
    )
    with pytest.raises(ValueError, match="pipe-delimited string column names"):
        rsplot.vj_usage(tuple_wide)

    tuple_long = pd.DataFrame(
        [
            {
                "sample_id": "sample0",
                "chain": "TRB",
                "v": "TRBV2",
                "j": "TRBJ1-1",
                "vj": ("TRBV2", "TRBJ1-1"),
                "usage": 1.0,
            }
        ]
    )
    with pytest.raises(ValueError, match="pipe-delimited strings"):
        rsplot.vj_usage(tuple_long)


def test_vjlen_usage_plots_two_sample_frequencies_and_default_style():
    fig = rsplot.vjlen_usage(_vjlen_usage_long())

    collection = fig.axes[0].collections[0]
    offsets = np.asarray(collection.get_offsets())
    assert any(np.allclose(point, [0.2, 0.4]) for point in offsets)
    assert np.allclose(collection.get_edgecolors()[0], to_rgba("black", 0.6))
    assert np.allclose(collection.get_facecolors()[0], to_rgba("#d62728", 0.6))
    assert fig.axes[0].get_xlabel() == "sample0"
    assert fig.axes[0].get_ylabel() == "sample1"


def test_vjlen_usage_combines_families_within_each_length():
    usage = pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "chain": "TRB",
                "v": v_gene,
                "j": j_gene,
                "len": 15,
                "vjlen": f"{v_gene}|{j_gene}|15",
                "usage": usage_value,
            }
            for sample_id, values in {
                "sample0": [
                    ("TRBV7-2", "TRBJ2-1", 0.2),
                    ("TRBV7-3", "TRBJ2-2", 0.3),
                ],
                "sample1": [
                    ("TRBV7-2", "TRBJ2-1", 0.4),
                    ("TRBV7-3", "TRBJ2-2", 0.1),
                ],
            }.items()
            for v_gene, j_gene, usage_value in values
        ]
    )

    legacy_usage = usage.drop(columns=["v", "j", "len"])
    fig = rsplot.vjlen_usage(
        legacy_usage,
        combine_families=True,
        labels=True,
    )
    ax = fig.axes[0]
    offsets = np.asarray(ax.collections[0].get_offsets())

    assert offsets.shape == (1, 2)
    assert np.allclose(offsets[0], [0.5, 0.5])
    assert {text.get_text() for text in ax.texts} == {"V7|J2|15"}


def test_vjlen_usage_grouped_axes_are_group_means():
    usage = _vjlen_usage_long(("s1", "s2", "s3", "s4"))
    metadata = pd.DataFrame(
        [
            {"sample_id": "s1", "chain": "TRB", "condition": "control"},
            {"sample_id": "s2", "chain": "TRB", "condition": "control"},
            {"sample_id": "s3", "chain": "TRB", "condition": "treated"},
            {"sample_id": "s4", "chain": "TRB", "condition": "treated"},
        ]
    )

    fig = rsplot.vjlen_usage(usage, metadata=metadata, group="condition")

    offsets = np.asarray(fig.axes[0].collections[0].get_offsets())
    assert any(np.allclose(point, [0.3, 0.7]) for point in offsets)
    assert fig.axes[0].get_xlabel() == "control"
    assert fig.axes[0].get_ylabel() == "treated"


def test_vjlen_usage_supports_two_split_facet_grid():
    sample_ids = [f"s{index}" for index in range(8)]
    usage = _vjlen_usage_long(tuple(sample_ids))
    metadata_rows = []
    for row_index, row_value in enumerate(["r1", "r2"]):
        for column_index, column_value in enumerate(["c1", "c2"]):
            start = (row_index * 2 + column_index) * 2
            for sample_id in sample_ids[start:start + 2]:
                metadata_rows.append(
                    {
                        "sample_id": sample_id,
                        "chain": "TRB",
                        "row_group": row_value,
                        "column_group": column_value,
                    }
                )
    metadata = pd.DataFrame(metadata_rows)

    fig = rsplot.vjlen_usage(
        usage,
        metadata=metadata,
        split=["row_group", "column_group"],
    )

    assert len(fig.axes) == 4
    assert [ax.get_title() for ax in fig.axes] == [
        "r1 | c1",
        "r1 | c2",
        "r2 | c1",
        "r2 | c2",
    ]


def test_vjlen_usage_log_scale_adds_compact_isolated_labels():
    fig = rsplot.vjlen_usage(
        _vjlen_usage_long(),
        log_scale=True,
        labels=True,
    )

    assert fig.axes[0].get_xscale() == "log"
    assert fig.axes[0].get_yscale() == "log"
    labels = {text.get_text() for text in fig.axes[0].texts}
    assert "V2|J1-1|15" in labels
    assert "V10|J2-1|14" in labels


def test_vjlen_usage_rejects_multiple_chains_and_wrong_panel_size():
    multiple_chains = pd.concat(
        [
            _vjlen_usage_long(),
            pd.DataFrame(
                [
                    {
                        "sample_id": "sample0",
                        "chain": "TRA",
                        "v": "TRAV2",
                        "j": "TRAJ1",
                        "len": 15,
                        "vjlen": "TRAV2|TRAJ1|15",
                        "usage": 0.1,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="exactly one chain"):
        rsplot.vjlen_usage(multiple_chains)

    with pytest.raises(ValueError, match="exactly 2 samples"):
        rsplot.vjlen_usage(_vjlen_usage_long(("s1", "s2", "s3")))


def _beta_metric_matrix():
    matrix = pd.DataFrame(
        [
            [1.0, 0.0004, 2.5e-6],
            [0.0004, 1.0, 0.0],
            [2.5e-6, 0.0, 1.0],
        ],
        index=["s1", "s2", "s3"],
        columns=["s1", "s2", "s3"],
    )
    matrix.index.name = "sample1"
    matrix.columns.name = "sample2"
    return matrix


def test_beta_metric_formats_one_without_decimal():
    assert rsplot._format_beta_value(1.0) == "1"
    assert rsplot._format_beta_value(2.0) == "2.0"


def _beta_full_table_same_set():
    rows = []
    pair_values = {
        ("s1", "s2"): [(0.4, 0.2), (0.2, 0.4), (0.4, 0.0), (0.0, 0.4)],
        ("s1", "s3"): [(0.4, 0.2), (0.2, 0.4), (0.4, 0.0), (0.0, 0.4)],
        ("s2", "s3"): [(0.4, 0.2), (0.2, 0.4), (0.4, 0.0), (0.0, 0.4)],
    }
    for (sample1, sample2), values in pair_values.items():
        for clone_number, (frequency1, frequency2) in enumerate(values):
            rows.append(
                {
                    "cdr3aa": f"CASS{clone_number}",
                    "sample1": sample1,
                    "sample2": sample2,
                    "sample1_freq": frequency1,
                    "sample2_freq": frequency2,
                }
            )
    table = pd.DataFrame(rows)
    table.attrs["sample_list"] = ["s1", "s2", "s3"]
    table.attrs["sample_list2"] = None
    return table


def _beta_full_table_two_sets():
    rows = []
    for row_sample in ["r1", "r2"]:
        for column_sample in ["c1", "c2"]:
            rows.extend(
                [
                    {
                        "cdr3aa": "CASS1",
                        "sample1": row_sample,
                        "sample2": column_sample,
                        "sample1_freq": 0.75,
                        "sample2_freq": 0.25,
                    },
                    {
                        "cdr3aa": "CASS2",
                        "sample1": row_sample,
                        "sample2": column_sample,
                        "sample1_freq": 0.25,
                        "sample2_freq": 0.75,
                    },
                ]
            )
    table = pd.DataFrame(rows)
    table.attrs["sample_list"] = ["r1", "r2"]
    table.attrs["sample_list2"] = ["c1", "c2"]
    return table


def test_beta_metric_dictionary_requires_metric_selection(capsys):
    beta_results = {"f2": _beta_metric_matrix(), "jaccard": _beta_metric_matrix()}

    with pytest.warns(UserWarning, match="No beta-diversity metric was selected"):
        fig = rsplot.beta_metric(beta_results)

    assert fig is None
    output = capsys.readouterr().out
    assert "f2" in output
    assert "jaccard" in output
    assert "metric='<key>'" in output


def test_beta_metric_logs_colors_and_displays_original_values_with_metadata():
    matrix = _beta_metric_matrix()
    metadata = pd.DataFrame(
        [
            {"sample_id": "s1", "condition": "control"},
            {"sample_id": "s2", "condition": "treated"},
            {"sample_id": "s3", "condition": "treated"},
        ]
    )

    fig = rsplot.beta_metric(
        {"f2": matrix},
        metric="f2",
        metadata=metadata,
        group="condition",
        hclust=False,
        ignore_diagonal=True,
        log_values=True,
        show_values=True,
    )
    heatmap = next(ax for ax in fig.axes if ax.get_title() == "f2")
    labels = {text.get_text() for text in heatmap.texts}

    assert "NA" in labels
    assert "0.0004" in labels
    assert "2.5e-6" in labels
    assert "0" in labels
    assert heatmap.collections[0].cmap.name == "pheatmap_default"
    assert any(ax.get_ylabel() == "log10(f2)" for ax in fig.axes)
    assert {text.get_text() for text in fig.legends[0].get_texts()} == {
        "Condition: control",
        "Condition: treated",
    }


def test_beta_metric_rejects_more_than_three_annotation_groups():
    metadata = pd.DataFrame(
        [{"sample_id": "s1", "a": "a", "b": "b", "c": "c", "d": "d"}]
    )

    with pytest.raises(ValueError, match="group can contain at most 3"):
        rsplot.beta_metric(
            _beta_metric_matrix().iloc[:1, :1],
            metadata=metadata,
            group=["a", "b", "c", "d"],
        )


def test_beta_table_dots_uses_lower_triangle_and_upper_f2_values():
    fig = rsplot.beta_table(
        {"full_table": _beta_full_table_same_set()},
        plot_type="dots",
        log_scale=True,
    )

    scatter_axes = [ax for ax in fig.axes if ax.collections]
    assert len(scatter_axes) == 3
    assert all(
        collection.get_alpha() == 0.5
        for ax in scatter_axes
        for collection in ax.collections
    )
    assert all(
        ax.get_xscale() == "log" and ax.get_yscale() == "log"
        for ax in scatter_axes
    )
    assert sum(
        text.get_text().startswith("F2")
        for ax in fig.axes
        for text in ax.texts
    ) == 3
    assert all(
        any(line.get_linestyle() == "--" for line in ax.lines)
        for ax in scatter_axes
    )


def test_beta_table_dots_tiles_two_sample_sets():
    fig = rsplot.beta_table(_beta_full_table_two_sets(), plot_type="dots")

    assert len(fig.axes) == 4
    assert all(len(ax.collections) == 1 for ax in fig.axes)


def test_beta_table_diff_uses_facets_and_matrix_tiles():
    same_set = rsplot.beta_table(
        _beta_full_table_same_set(),
        plot_type="diff",
        top=1,
    )
    visible_same_set = [ax for ax in same_set.axes if ax.get_visible()]
    assert len(visible_same_set) == 3
    assert all(ax.patches for ax in visible_same_set)
    for ax in visible_same_set:
        labels = {text.get_text() for text in ax.texts}
        assert "NonOverlapping" in labels
        assert "NotShown" in labels
        assert "CASS0" in labels
        assert ax.get_ylabel() == "Cumulative abundance"

    two_sets = rsplot.beta_table(
        _beta_full_table_two_sets(),
        plot_type="diff",
        top=2,
    )
    assert len(two_sets.axes) == 4
    assert all(ax.patches for ax in two_sets.axes)


def test_rarefaction_curve_plots_chain_aware_sample_labels():
    rarefaction = pd.DataFrame([
        {"sample_id": "ucb_ntreg", "chain": "TRA", "rarefaction_depth": 32, "diversity": 10},
        {"sample_id": "ucb_ntreg", "chain": "TRA", "rarefaction_depth": 100, "diversity": 20},
        {"sample_id": "ucb_ntreg", "chain": "TRB", "rarefaction_depth": 32, "diversity": 12},
        {"sample_id": "ucb_ntreg", "chain": "TRB", "rarefaction_depth": 100, "diversity": 25},
        {"sample_id": "unique", "chain": "TRB", "rarefaction_depth": 32, "diversity": 8},
        {"sample_id": "unique", "chain": "TRB", "rarefaction_depth": 100, "diversity": 15},
    ])

    fig = rsplot.rarefaction_curve(rarefaction)
    ax = fig.axes[0]
    labels = {text.get_text() for text in ax.get_legend().get_texts()}

    assert labels == {"ucb_ntreg(TRA)", "ucb_ntreg(TRB)", "unique"}
    legend = ax.get_legend()
    assert legend._loc == 2
    assert legend.get_bbox_to_anchor()._bbox.x0 > 1
    assert to_rgba(ax.lines[0].get_color()) == to_rgba("#4e79a7")
    assert to_rgba(ax.lines[1].get_color()) == to_rgba("#a0cde8")
    assert ax.get_xscale() == "log"
    assert ax.get_xlabel() == "Rarefaction depth"
    assert ax.get_ylabel() == "Observed diversity"


def test_rarefaction_curve_validates_columns():
    with pytest.raises(ValueError, match="must contain columns"):
        rsplot.rarefaction_curve(pd.DataFrame({"sample_id": ["s1"]}))
