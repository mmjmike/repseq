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
    assert heatmaps["TRA"].collections[0].cmap.name == "RdBu_r"
    assert fig.number not in plt.get_fignums()
    assert [tick.get_text() for tick in heatmaps["TRA"].get_xticklabels()] == [
        "TRAV2",
        "TRAV10",
    ]
    assert [tick.get_text() for tick in heatmaps["TRB"].get_xticklabels()] == [
        "TRBV2"
    ]


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
