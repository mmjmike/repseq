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


def _beta_mds_table():
    return pd.DataFrame(
        [
            {"sample_id": "s1", "MDS1": -2.0, "MDS2": -1.0},
            {"sample_id": "s2", "MDS1": -1.0, "MDS2": 0.0},
            {"sample_id": "s3", "MDS1": 0.0, "MDS2": -1.0},
            {"sample_id": "s4", "MDS1": 1.0, "MDS2": 1.0},
            {"sample_id": "s5", "MDS1": 2.0, "MDS2": 0.0},
            {"sample_id": "s6", "MDS1": 3.0, "MDS2": 1.0},
        ]
    )


def test_beta_mds_draws_groups_centroids_spokes_and_dispersion():
    metadata = pd.DataFrame(
        {
            "sample_id": [f"s{number}" for number in range(1, 7)],
            "condition": ["control"] * 3 + ["treated"] * 3,
        }
    )
    fig = rsplot.beta_mds(
        _beta_mds_table(),
        metadata=metadata,
        group="condition",
        centroids=True,
        dispersion=True,
    )
    axis = fig.axes[0]

    assert fig.number not in plt.get_fignums()
    assert axis.get_xlabel() == "MDS1"
    assert axis.get_ylabel() == "MDS2"
    assert len(axis.lines) == 6
    assert len(axis.patches) == 2
    assert len(axis.collections) == 4
    assert all(patch.get_facecolor()[-1] == pytest.approx(0.2) for patch in axis.patches)
    assert {text.get_text() for text in fig.legends[0].get_texts()} == {
        "control",
        "treated",
    }


def test_beta_mds_supports_two_ordered_splits_and_one_group_maximum():
    metadata = pd.DataFrame(
        {
            "sample_id": [f"s{number}" for number in range(1, 7)],
            "condition": ["control"] * 3 + ["treated"] * 3,
            "tissue": ["blood", "tumor", "blood", "tumor", "blood", "tumor"],
            "sex": ["F", "M", "F", "M", "F", "M"],
        }
    )
    metadata["tissue"] = pd.Categorical(
        metadata["tissue"], categories=["tumor", "blood"], ordered=True
    )
    metadata["sex"] = pd.Categorical(
        metadata["sex"], categories=["M", "F"], ordered=True
    )
    fig = rsplot.beta_mds(
        _beta_mds_table(), metadata=metadata, split=["tissue", "sex"]
    )

    assert [axis.get_title() for axis in fig.axes if axis.get_visible()] == [
        "tumor | M",
        "blood | F",
    ]
    with pytest.raises(ValueError, match="group can contain at most 1"):
        rsplot.beta_mds(
            _beta_mds_table(),
            metadata=metadata,
            group=["condition", "sex"],
        )


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


def test_plot_stats_uses_group_and_split_columns_from_stats_table():
    stats_df = _stats_df().assign(
        condition=["treated", "control"],
        batch=["b1", "b2"],
    )

    grouped = rsplot.diversity_stats(
        stats_df,
        properties=["diversity"],
        group="condition",
    )
    tick_labels = [
        tick.get_text() for tick in grouped.axes.flat[0].get_xticklabels()
    ]
    assert tick_labels == ["treated", "control"]

    split = rsplot.diversity_stats(
        stats_df,
        properties=["diversity"],
        split=["condition", "batch"],
    )
    assert [ax.get_title() for ax in split.axes.flat] == [
        "treated | b1",
        "control | b2",
    ]


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


def test_plot_stats_ignores_metadata_groups_without_plotted_samples():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "experimental_group": "term"},
            {"sample_id": "sample2", "chain": "TRB", "experimental_group": "term"},
            {"sample_id": "absent", "chain": "TRB", "experimental_group": "preterm"},
        ]
    )
    metadata["experimental_group"] = pd.Categorical(
        metadata["experimental_group"],
        categories=["term", "preterm"],
        ordered=True,
    )

    grouped = rsplot.plot_stats(
        _stats_df(),
        metadata=metadata,
        properties=["diversity"],
        group="experimental_group",
    )

    assert [tick.get_text() for tick in grouped.axes.flat[0].get_xticklabels()] == [
        "term"
    ]


def test_plot_stats_ignores_metadata_splits_without_plotted_samples():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRA", "experimental_group": "term"},
            {"sample_id": "sample2", "chain": "TRB", "experimental_group": "term"},
            {"sample_id": "absent", "chain": "TRB", "experimental_group": "preterm"},
        ]
    )
    metadata["experimental_group"] = pd.Categorical(
        metadata["experimental_group"],
        categories=["term", "preterm"],
        ordered=True,
    )

    split = rsplot.plot_stats(
        _stats_df(),
        metadata=metadata,
        properties=["diversity"],
        split="experimental_group",
    )

    assert [ax.get_title() for ax in split.axes.flat] == ["term"]


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


def _processing_table():
    return pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "extracted_chain": "TRA",
                "reads_aligned_pc": 91.0,
                "reads_per_umi": 3.0,
                "clones_func": 12,
                "umi_in_func_clones": 40,
            },
            {
                "sample_id": "sample1",
                "extracted_chain": "TRB",
                "reads_aligned_pc": 92.0,
                "reads_per_umi": 4.0,
                "clones_func": 18,
                "umi_in_func_clones": 55,
            },
            {
                "sample_id": "sample2",
                "extracted_chain": "TRB",
                "reads_aligned_pc": 88.0,
                "reads_per_umi": 2.5,
                "clones_func": 9,
                "umi_in_func_clones": 30,
            },
        ]
    )


def test_processing_uses_default_properties_and_chain_sample_labels():
    grid = rsplot.processing(_processing_table())

    assert [ax.get_title() for ax in grid.axes.flat] == [
        "Reads Aligned Pc",
        "Reads Per Umi",
        "Clones Func",
        "Umi In Func Clones",
    ]
    expected_labels = ["sample1 (TRA)", "sample1 (TRB)", "sample2 (TRB)"]
    visible_labels_by_axis = [
        [tick.get_text() for tick in ax.get_xticklabels() if tick.get_text()]
        for ax in grid.axes.flat
        if ax.get_xticklabels()
    ]
    assert visible_labels_by_axis
    assert all(labels == expected_labels for labels in visible_labels_by_axis)
    assert all(ax.get_ylim()[0] == 0 for ax in grid.axes.flat)


def test_processing_keeps_extracted_chain_separate_from_metadata_chain(monkeypatch):
    captured = {}

    def capture_plot_stats(stats_df, **kwargs):
        captured["stats_df"] = stats_df
        captured.update(kwargs)
        return "grid"

    monkeypatch.setattr(rsplot, "plot_stats", capture_plot_stats)
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TCR", "condition": "control"},
            {"sample_id": "sample2", "chain": "BCR", "condition": "treated"},
        ]
    )
    palette = {"control": "#336699", "treated": "#cc5500"}

    assert rsplot.processing(
        _processing_table(),
        metadata=metadata,
        group="condition",
        palette=palette,
    ) == "grid"
    assert "extracted_chain" in captured["stats_df"].columns
    assert "chain" not in captured["stats_df"].columns
    assert captured["properties"] == rsplot.PROCESSING_PROPERTIES
    assert captured["metadata"] is metadata
    assert captured["group"] == "condition"
    assert captured["palette"] is palette
    assert captured["zero_bottom"] is True


def test_processing_merges_unique_sample_metadata_by_sample_id():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TCR", "condition": "control"},
            {"sample_id": "sample2", "chain": "BCR", "condition": "treated"},
        ]
    )

    grid = rsplot.processing(
        _processing_table(),
        metadata=metadata,
        properties=["clones_func"],
        group="condition",
    )

    tick_labels = [tick.get_text() for tick in grid.axes.flat[0].get_xticklabels()]
    assert tick_labels == ["control", "treated"]


def test_processing_combines_metadata_group_with_table_split():
    metadata = pd.DataFrame(
        [
            {"sample_id": "sample1", "Primers": "primer_a"},
            {"sample_id": "sample2", "Primers": "primer_b"},
        ]
    )

    grid = rsplot.processing(
        _processing_table(),
        metadata=metadata,
        properties=["clones_func"],
        group="Primers",
        split="extracted_chain",
        aspect=1.5,
    )

    assert [ax.get_title() for ax in grid.axes.flat] == ["TRA", "TRB"]


def test_stats_wrappers_accept_custom_property_columns():
    custom_stats = _stats_df().assign(custom_score=[2.5, 4.5])
    custom_processing = _processing_table().assign(custom_score=[2.5, 3.5, 4.5])
    custom_clonosets = _clonoset_stats_df().assign(custom_score=[2.5, 4.5])

    plot_calls = [
        (rsplot.processing, custom_processing),
        (rsplot.cdr3aa_stats, custom_stats),
        (rsplot.diversity_stats, custom_stats),
        (rsplot.convergence, custom_stats),
        (rsplot.clonoset_stats, custom_clonosets),
    ]

    for plot_function, stats_table in plot_calls:
        grid = plot_function(stats_table, properties=["custom_score"])
        assert grid.axes.flat[0].get_title() == "Custom Score"
        plt.close(grid.fig)


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


def test_clonoset_stats_default_uses_split_from_stats_table():
    stats_df = _clonoset_stats_df().assign(tissue=["blood", "tumor"])

    grid = rsplot.clonoset_stats(stats_df, split="tissue")

    assert [ax.get_title() for ax in grid.axes[0]] == [
        "Reads | blood",
        "Reads | tumor",
    ]


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


def test_segment_usage_parses_and_sorts_alleles_alphabetically():
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "TRB", "v": allele, "usage": value}
            for allele, value in [
                ("TRBV6-2*A1", 0.2),
                ("TRBV6-2*02", 0.3),
                ("TRBV6-2*010", 0.5),
            ]
        ]
    )

    assert rsplot.parse_gene_name("TRBV6-2*LONG")["allele"] == "LONG"

    fig = rsplot.segment_usage(usage)
    heatmap = next(ax for ax in fig.axes if ax.get_title() == "TRB")

    assert [tick.get_text() for tick in heatmap.get_xticklabels()] == [
        "TRBV6-2*010",
        "TRBV6-2*02",
        "TRBV6-2*A1",
    ]


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


def test_segment_usage_keeps_missing_calls_as_na_family():
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "IGH", "c": "IGHG1", "usage": 0.5},
            {"sample_id": "sample1", "chain": "IGH", "c": pd.NA, "usage": 0.2},
            {"sample_id": "sample1", "chain": "IGH", "c": "<NA>", "usage": 0.3},
        ]
    )

    fig = rsplot.segment_usage(
        usage,
        plot_type="barplot",
        combine_families=True,
    )
    ax = fig.axes[0]

    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["IGHG", "NA"]
    np.testing.assert_allclose(
        [patch.get_height() for patch in ax.patches],
        [0.5, 0.5],
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


@pytest.mark.parametrize(
    ("gene", "expected"),
    [
        ("IGHA", "IgA"),
        ("IGHA1", "IgA1"),
        ("IGHA2", "IgA2"),
        ("IGHD", "IgD"),
        ("IGHE", "IgE"),
        ("IGHEP1", "NA"),
        ("IGHG1", "IgG1"),
        ("IGHG1A", "IgG1A"),
        ("IGHG1B", "IgG1B"),
        ("IGHG2", "IgG2"),
        ("IGHG2B", "IgG2B"),
        ("IGHG2B_hinge", "IgG2B"),
        ("IGHG2C", "IgG2C"),
        ("IGHG2C_hinge", "IgG2C"),
        ("IGHG3", "IgG3"),
        ("IGHG4", "IgG4"),
        ("IGHGP", "NA"),
        ("IGHM", "IgM"),
        ("IGHM1", "IgM1"),
        ("IGHM2", "IgM2"),
        (".", "NA"),
    ],
)
def test_recode_isotype_supports_all_mixcr_igh_constant_genes(gene, expected):
    assert rsplot._recode_isotype(gene) == expected


def test_isotype_fraction_keeps_missing_calls_as_na_isotype():
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "IGH", "c": "IGHM", "usage": 0.5},
            {"sample_id": "sample1", "chain": "IGH", "c": pd.NA, "usage": 0.2},
            {"sample_id": "sample1", "chain": "IGH", "c": "<NA>", "usage": 0.3},
        ]
    )

    fig = rsplot.isotype_fraction(usage)
    containers = {
        container.get_label(): container for container in fig.axes[0].containers
    }

    assert [text.get_text() for text in fig.axes[0].get_legend().get_texts()] == [
        "IgM",
        "NA",
    ]
    assert containers["IgM"].patches[0].get_width() == pytest.approx(0.5)
    assert containers["NA"].patches[0].get_width() == pytest.approx(0.5)


def test_isotype_fraction_uses_requested_order_palette_and_right_to_left_stack():
    genes = [
        "IGHM",
        "IGHD",
        "IGHG1",
        "IGHG2",
        "IGHG3",
        "IGHG4",
        "IGHA1",
        "IGHA2",
        "IGHE",
        ".",
    ]
    usage = pd.DataFrame(
        [
            {"sample_id": "sample1", "chain": "IGH", "c": gene, "usage": 1}
            for gene in genes
        ]
    )

    fig = rsplot.isotype_fraction(usage)
    ax = fig.axes[0]
    expected_order = [
        "IgM",
        "IgD",
        "IgG1",
        "IgG2",
        "IgG3",
        "IgG4",
        "IgA1",
        "IgA2",
        "IgE",
        "NA",
    ]
    expected_colors = [
        "#E41A1C",
        "#FF7F00",
        "#4DAF4A",
        "#74C476",
        "#A1D99B",
        "#D9F0D3",
        "#377EB8",
        "#6BAED6",
        "#984EA3",
        "#999999",
    ]

    assert [text.get_text() for text in ax.get_legend().get_texts()] == expected_order
    container_colors = {
        container.get_label(): container.patches[0].get_facecolor()
        for container in ax.containers
    }
    for isotype, color in zip(expected_order, expected_colors):
        assert np.allclose(container_colors[isotype], to_rgba(color))
    assert ax.containers[-1].get_label() == "IgM"
    assert ax.containers[-1].patches[0].get_x() == pytest.approx(0.9)
    assert ax.containers[-1].patches[0].get_width() == pytest.approx(0.1)
    assert ax.get_xlabel() == "Fraction"
    assert ax.get_ylabel() == "Sample ID"


def test_isotype_fraction_keeps_canonical_color_when_variants_are_absent():
    usage = pd.DataFrame(
        [
            {"sample_id": "s1", "chain": "IGH", "c": "IGHG4", "usage": 1},
            {"sample_id": "s1", "chain": "IGH", "c": "IGHA2", "usage": 1},
        ]
    )

    fig = rsplot.isotype_fraction(usage)
    containers = {
        container.get_label(): container for container in fig.axes[0].containers
    }

    assert np.allclose(
        containers["IgG4"].patches[0].get_facecolor(), to_rgba("#D9F0D3")
    )
    assert np.allclose(
        containers["IgA2"].patches[0].get_facecolor(), to_rgba("#6BAED6")
    )


def test_isotype_fraction_uses_unique_metadata_labels():
    usage = pd.DataFrame(
        [
            {"sample_id": "s1", "chain": "IGH", "c": "IGHM", "usage": 1},
            {"sample_id": "s2", "chain": "IGH", "c": "IGHA1", "usage": 1},
        ]
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "display_name": ["Donor 1", "Donor 2"]}
    )

    fig = rsplot.isotype_fraction(
        usage, metadata=metadata, label="display_name"
    )

    assert [tick.get_text() for tick in fig.axes[0].get_yticklabels()] == [
        "Donor 1",
        "Donor 2",
    ]
    assert fig.axes[0].get_ylabel() == "display_name"

    metadata["display_name"] = "duplicate"
    with pytest.raises(ValueError, match="label values must be unique"):
        rsplot.isotype_fraction(
            usage, metadata=metadata, label="display_name"
        )


def test_isotype_count_preserves_values_and_combines_families():
    usage = pd.DataFrame(
        [
            {"sample_id": "s1", "chain": "IGH", "c": "IGHM1", "usage": 2},
            {"sample_id": "s1", "chain": "IGH", "c": "IGHM2", "usage": 3},
            {"sample_id": "s1", "chain": "IGH", "c": "IGHG1A", "usage": 5},
            {"sample_id": "s1", "chain": "IGH", "c": "IGHG2B_hinge", "usage": 7},
            {"sample_id": "s1", "chain": "IGH", "c": "IGHEP1", "usage": 11},
        ]
    )

    fig = rsplot.isotype_count(usage, combine_families=True)
    ax = fig.axes[0]
    containers = {container.get_label(): container for container in ax.containers}

    assert [text.get_text() for text in ax.get_legend().get_texts()] == [
        "IgM",
        "IgG",
        "NA",
    ]
    assert containers["IgM"].patches[0].get_width() == pytest.approx(5)
    assert containers["IgG"].patches[0].get_width() == pytest.approx(12)
    assert containers["NA"].patches[0].get_width() == pytest.approx(11)
    assert np.allclose(
        containers["IgM"].patches[0].get_facecolor(), to_rgba("#E41A1C")
    )
    assert np.allclose(
        containers["IgG"].patches[0].get_facecolor(), to_rgba("#4DAF4A")
    )
    assert ax.get_xlabel() == "Count"


def test_isotype_fraction_accepts_wide_hinge_columns():
    usage = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "chain": "IGH",
                "IGHG2B_hinge": 0.75,
                "IGHEP1": 0.25,
            }
        ]
    )

    fig = rsplot.isotype_fraction(usage)

    assert [text.get_text() for text in fig.axes[0].get_legend().get_texts()] == [
        "IgG2B",
        "NA",
    ]


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


def _pairing_matrix():
    matrix = pd.DataFrame(
        [[0.0, 0.01], [0.1, np.nan]],
        index=pd.Index(["chain2_a", "chain2_b"], name="chain2_feature"),
        columns=pd.Index(["chain1_a", "chain1_b"], name="chain1_feature"),
    )
    matrix.attrs["method"] = "jsd"
    return matrix


def test_de_pairing_uses_negative_log10_colors_and_original_cell_values():
    fig = rsplot.de_pairing(_pairing_matrix(), log_minus=True)
    heatmap = next(ax for ax in fig.axes if ax.get_title() == "Chain pairing")
    plotted = np.ma.filled(
        np.ma.asarray(heatmap.collections[0].get_array(), dtype=float),
        np.nan,
    )

    np.testing.assert_allclose(
        plotted,
        [[3.0, 2.0], [1.0, np.nan]],
        equal_nan=True,
    )
    assert {text.get_text() for text in heatmap.texts} == {
        "0",
        "0.01",
        "0.1",
        "NA",
    }
    assert heatmap.collections[0].cmap.name == "pheatmap_default"
    assert heatmap.get_xlabel() == "chain1_feature"
    assert heatmap.get_ylabel() == "chain2_feature"
    assert any(ax.get_ylabel() == "-log10(JSD)" for ax in fig.axes)


def test_de_pairing_rejects_negative_log_transformed_values():
    matrix = _pairing_matrix()
    matrix.iloc[0, 0] = -0.1

    with pytest.raises(ValueError, match="cannot contain negatives"):
        rsplot.de_pairing(matrix, log_minus=True)


def _de_heatmap_inputs(sample_order=None):
    sample_order = sample_order or ["a1", "b1", "a2", "b2"]
    sample_values = {
        "a1": [456789, 0],
        "a2": [4, 0],
        "b1": [0, 6],
        "b2": [0, 7],
    }
    statistics_table = pd.DataFrame(
        {
            "feature": ["feature_a", "feature_b"],
            "enriched_in": ["A", "B"],
            "method": ["mann_whitney", "mann_whitney"],
            "mean_group_count": [4.5, 6.5],
            "log2FC": [100, 100],
            "p_val": [0.01, 0.02],
            "p_adj": [0.02, 0.02],
            **{sample: sample_values[sample] for sample in sample_order},
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["a1", "a2", "b1", "b2"],
            "group": ["A", "A", "B", "B"],
        }
    )
    return statistics_table, metadata


def test_de_heatmap_groups_interleaved_samples_and_matches_annotation_colors():
    statistics_table, metadata = _de_heatmap_inputs()

    fig = rsplot.de_heatmap(
        statistics_table,
        metadata,
        log_values=True,
        show_values=True,
    )

    heatmap = next(
        ax for ax in fig.axes if ax.get_title() == "Differential enrichment"
    )
    row_annotation = next(
        ax for ax in fig.axes if [tick.get_text() for tick in ax.get_xticklabels()] == ["Enriched in"]
    )
    column_annotation = next(
        ax for ax in fig.axes if [tick.get_text() for tick in ax.get_yticklabels()] == ["Group"]
    )

    assert [tick.get_text() for tick in heatmap.get_xticklabels()] == [
        "a1",
        "a2",
        "b1",
        "b2",
    ]
    assert [tick.get_text() for tick in heatmap.get_yticklabels()] == [
        "feature_a",
        "feature_b",
    ]
    assert heatmap.collections[0].cmap.name == "pheatmap_default"
    assert any(np.array_equal(line.get_xdata(), [2, 2]) for line in heatmap.lines)
    assert any(
        np.array_equal(line.get_xdata(), [2, 2])
        for line in column_annotation.lines
    )
    np.testing.assert_allclose(
        row_annotation.collections[0].get_facecolors(),
        column_annotation.collections[0].get_facecolors()[[0, 2]],
    )
    assert {text.get_text() for text in heatmap.texts} >= {
        "0",
        "4",
        "6",
        "7",
        "456789",
    }
    assert any(ax.get_ylabel() == "log10(Count)" for ax in fig.axes)
    assert [text.get_text() for text in fig.legends[0].get_texts()] == ["A", "B"]


def test_de_heatmap_preserves_samples_when_groups_are_already_contiguous():
    statistics_table, metadata = _de_heatmap_inputs(
        sample_order=["b1", "b2", "a1", "a2"]
    )

    fig = rsplot.de_heatmap(
        statistics_table,
        metadata,
        show_values=False,
    )
    heatmap = next(
        ax for ax in fig.axes if ax.get_title() == "Differential enrichment"
    )

    assert [tick.get_text() for tick in heatmap.get_xticklabels()] == [
        "b1",
        "b2",
        "a1",
        "a2",
    ]


def test_de_heatmap_integer_formatter_preserves_real_decimals():
    assert rsplot._format_de_heatmap_value(456789.0) == "456789"
    assert rsplot._format_de_heatmap_value(2.0) == "2"
    assert rsplot._format_de_heatmap_value(0.25) == "0.25"


def test_de_heatmap_requires_enriched_groups_in_metadata():
    statistics_table, metadata = _de_heatmap_inputs()
    statistics_table.loc[0, "enriched_in"] = "MissingGroup"

    with pytest.raises(ValueError, match="enriched_in groups are absent"):
        rsplot.de_heatmap(statistics_table, metadata)


@pytest.mark.parametrize(
    ("pass_column", "pass_values", "expected_feature"),
    [
        ("postfilter_pass", [True, False], "feature_a"),
        ("prefilter_pass", [False, True], "feature_b"),
    ],
)
def test_de_heatmap_uses_only_features_passing_available_filter_column(
    pass_column,
    pass_values,
    expected_feature,
):
    statistics_table, metadata = _de_heatmap_inputs()
    statistics_table[pass_column] = pass_values
    statistics_table.loc[
        statistics_table[pass_column].eq(False), "enriched_in"
    ] = None

    fig = rsplot.de_heatmap(statistics_table, metadata, show_values=False)
    heatmap = next(
        ax for ax in fig.axes if ax.get_title() == "Differential enrichment"
    )

    assert [tick.get_text() for tick in heatmap.get_yticklabels()] == [
        expected_feature
    ]


def test_de_heatmap_combines_prefilter_and_postfilter_pass_columns():
    statistics_table, metadata = _de_heatmap_inputs()
    statistics_table["prefilter_pass"] = [True, True]
    statistics_table["postfilter_pass"] = [True, False]

    fig = rsplot.de_heatmap(statistics_table, metadata, show_values=False)
    heatmap = next(
        ax for ax in fig.axes if ax.get_title() == "Differential enrichment"
    )

    assert [tick.get_text() for tick in heatmap.get_yticklabels()] == [
        "feature_a"
    ]


def _de_volcano_table():
    return pd.DataFrame(
        {
            "log2FC": [100, 100, 11, 5, np.nan],
            "p_adj": [0.01, 0.02, 0.03, 0.04, 0.05],
            "p_val": [0.1, 0.2, 0.3, 0.4, 0.5],
            "mean_group_count": [10, 20, 5, 40, 50],
            "enriched_in": ["A", "A", "B", "B", "A"],
        }
    )


def test_de_volcano_uses_adjusted_p_values_sizes_colors_and_drops_na():
    fig = rsplot.de_volcano(_de_volcano_table())
    ax = fig.axes[0]
    points = ax.collections[0]
    offsets = np.asarray(points.get_offsets())

    np.testing.assert_allclose(offsets[:, 0], [13, 13, 11, 5])
    np.testing.assert_allclose(
        offsets[:, 1],
        -np.log10([0.01, 0.02, 0.03, 0.04]),
    )
    assert len(offsets) == 4
    np.testing.assert_allclose(
        points.get_sizes(),
        rsplot._scaled_dot_sizes(np.log1p([10, 20, 5, 40]), (20, 300)),
    )
    assert points.get_sizes()[3] == points.get_sizes().max()
    np.testing.assert_allclose(
        points.get_facecolors()[0],
        points.get_facecolors()[1],
    )
    np.testing.assert_allclose(
        points.get_facecolors()[2],
        points.get_facecolors()[3],
    )
    assert not np.allclose(
        points.get_facecolors()[0],
        points.get_facecolors()[2],
    )
    assert ax.get_xlabel() == "log2FC"
    assert ax.get_ylabel() == "-log10(p_adj)"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["A", "B"]


def test_de_volcano_can_plot_raw_p_values():
    fig = rsplot.de_volcano(_de_volcano_table(), p_column="p_val")
    ax = fig.axes[0]
    offsets = np.asarray(ax.collections[0].get_offsets())

    np.testing.assert_allclose(
        offsets[:, 1],
        -np.log10([0.1, 0.2, 0.3, 0.4]),
    )
    assert ax.get_ylabel() == "-log10(p_val)"


def test_de_volcano_can_size_points_by_raw_mean_group_count():
    fig = rsplot.de_volcano(_de_volcano_table(), log_sizes=False)
    points = fig.axes[0].collections[0]

    np.testing.assert_allclose(
        points.get_sizes(),
        rsplot._scaled_dot_sizes([10, 20, 5, 40], (20, 300)),
    )


@pytest.mark.parametrize(
    ("p_column", "p_values"),
    [
        ("p_adj", [0.01, 0.02, 0.03, 0.04]),
        ("p_val", [0.1, 0.2, 0.3, 0.4]),
    ],
)
def test_de_volcano_by_mean_count_uses_counts_on_y_and_p_values_for_size(
    p_column, p_values
):
    fig = rsplot.de_volcano(
        _de_volcano_table(),
        p_column=p_column,
        by_mean_count=True,
        log_sizes=False,
    )
    ax = fig.axes[0]
    points = ax.collections[0]

    np.testing.assert_allclose(
        points.get_offsets(),
        [[13, 10], [13, 20], [11, 5], [5, 40]],
    )
    np.testing.assert_allclose(
        points.get_sizes(),
        rsplot._scaled_dot_sizes(-np.log10(p_values), (20, 300)),
    )
    assert ax.get_xlabel() == "log2FC"
    assert ax.get_ylabel() == "Mean group count"


def test_de_volcano_by_mean_count_log_sizes_transforms_vertical_axis():
    fig = rsplot.de_volcano(
        _de_volcano_table(),
        by_mean_count=True,
        log_sizes=True,
    )
    ax = fig.axes[0]

    np.testing.assert_allclose(
        ax.collections[0].get_offsets()[:, 1],
        np.log1p([10, 20, 5, 40]),
    )
    assert ax.get_ylabel() == "log1p(Mean group count)"


def test_de_volcano_draws_postfiltered_points_in_background():
    statistics_table = _de_volcano_table().iloc[:4].assign(
        enriched_in=["A", "C", "B", "B"],
        prefilter_pass=[True, True, False, True],
        postfilter_pass=[True, False, False, False],
    )

    fig = rsplot.de_volcano(statistics_table)
    ax = fig.axes[0]
    filtered_points, retained_points = ax.collections

    np.testing.assert_allclose(
        np.asarray(filtered_points.get_offsets()),
        [[13, -np.log10(0.02)], [5, -np.log10(0.04)]],
    )
    assert filtered_points.get_sizes().tolist() == [
        rsplot._scaled_dot_sizes(np.log1p([10, 20, 5, 40]), (20, 300)).min()
    ]
    np.testing.assert_allclose(filtered_points.get_facecolors()[0, :3], [0.6] * 3)
    assert filtered_points.get_alpha() == 0.3
    assert filtered_points.get_edgecolors().size == 0
    assert filtered_points.get_zorder() < retained_points.get_zorder()
    np.testing.assert_allclose(
        np.asarray(retained_points.get_offsets()),
        [[13, -np.log10(0.01)], [11, -np.log10(0.03)]],
    )
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [
        "A",
        "B",
        "filtered_out",
    ]


def test_de_volcano_repeated_high_log2fc_adjustment_uses_maximum_below_value():
    adjusted = rsplot._adjust_de_volcano_log2fc([100, 100, 11, 4])

    assert adjusted.tolist() == [13, 13, 11, 4]


def test_beta_table_dots_uses_lower_triangle_and_upper_f2_values():
    fig = rsplot.beta_table(
        {"full_table": _beta_full_table_same_set()},
        plot_type="dots",
        log_scale=True,
        matrix_layout=True,
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


def test_beta_table_dots_defaults_to_wrapped_pair_facets_without_f2_values():
    fig = rsplot.beta_table(_beta_full_table_same_set(), plot_type="dots")

    visible_axes = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible_axes) == 3
    assert all(len(ax.collections) == 1 for ax in visible_axes)
    assert {ax.get_title() for ax in visible_axes} == {
        "s1 vs s2",
        "s1 vs s3",
        "s2 vs s3",
    }
    assert not any(
        text.get_text().startswith("F2")
        for ax in fig.axes
        for text in ax.texts
    )


def test_beta_table_dots_pads_zero_values_inside_axes():
    linear = rsplot.beta_table(_beta_full_table_same_set(), plot_type="dots")
    linear_axes = [ax for ax in linear.axes if ax.collections]
    assert all(ax.get_xlim()[0] < 0 and ax.get_ylim()[0] < 0 for ax in linear_axes)

    logarithmic = rsplot.beta_table(
        _beta_full_table_same_set(),
        plot_type="dots",
        log_scale=True,
    )
    for ax in [axis for axis in logarithmic.axes if axis.collections]:
        offsets = np.asarray(ax.collections[0].get_offsets())
        assert ax.get_xlim()[0] < offsets[:, 0].min()
        assert ax.get_ylim()[0] < offsets[:, 1].min()


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


def _clonotypes_coverage_df():
    rows = []
    for sample_id, chain, values in [
        ("sample1", "TRA", [1, 3, 2, 1, 0, 0]),
        ("sample2", "TRB", [0, 1, 2, 3, 4, 0]),
    ]:
        for bin_label, value in zip(["1", "3", "10", "32", "100", "316"], values):
            rows.append({
                "sample_id": sample_id,
                "chain": chain,
                "bin": bin_label,
                "value": value,
            })
    return pd.DataFrame(rows)


def test_clonotypes_coverage_plots_numeric_bins_at_equal_spacing():
    grid = rsplot.clonotypes_coverage(_clonotypes_coverage_df())
    ax = grid.axes.flat[0]

    assert [tick.get_text() for tick in ax.get_xticklabels()] == [
        "1", "3", "10", "32", "100", "316"
    ]
    np.testing.assert_allclose(ax.get_xticks(), [0, 1, 2, 3, 4, 5])
    assert ax.get_xlabel() == "Clonotype size"
    assert {text.get_text() for text in grid.legend.texts} == {"sample1", "sample2"}


def test_clonotypes_coverage_supports_two_metadata_splits_and_rejects_group():
    metadata = pd.DataFrame([
        {"sample_id": "sample1", "chain": "TRA", "tissue": "blood", "sex": "F"},
        {"sample_id": "sample2", "chain": "TRB", "tissue": "tumor", "sex": "M"},
    ])

    grid = rsplot.clonotypes_coverage(
        _clonotypes_coverage_df(), metadata=metadata, split=["tissue", "sex"]
    )
    assert [ax.get_title() for ax in grid.axes.flat] == ["blood | F", "tumor | M"]

    with pytest.raises(ValueError, match="group is not supported"):
        rsplot.clonotypes_coverage(
            _clonotypes_coverage_df(), metadata=metadata, group="tissue"
        )


def test_clonotypes_coverage_separate_trims_panel_specific_high_zero_bins():
    coverage = _clonotypes_coverage_df()
    grid = rsplot.clonotypes_coverage(coverage, separate=True)
    axes = {ax.get_title(): ax for ax in grid.axes.flat}

    assert [tick.get_text() for tick in axes["sample1"].get_xticklabels()] == [
        "1", "3", "10", "32", "100"
    ]
    assert [tick.get_text() for tick in axes["sample2"].get_xticklabels()] == [
        "1", "3", "10", "32", "100"
    ]
    assert all(ax.get_xlabel() == "Clonotype size" for ax in axes.values())
    assert all(
        to_rgba(patch.get_facecolor()) == to_rgba("#CCCCCC")
        for ax in axes.values()
        for patch in ax.patches
    )

    untrimmed = rsplot.clonotypes_coverage(
        coverage, separate=True, trim_high_zero_bins=False
    )
    assert [tick.get_text() for tick in untrimmed.axes.flat[0].get_xticklabels()] == [
        "1", "3", "10", "32", "100", "316"
    ]


def test_clonotype_coverage_aliases_and_keeps_figure_open():
    grid = rsplot.clonotype_coverage(_clonotypes_coverage_df())

    assert plt.fignum_exists(grid.fig.number)


def test_cluster_properties_plot_uses_default_metrics(monkeypatch):
    captured = {}
    stats_df = pd.DataFrame(
        [
            {
                "sample_id": "sample1",
                "mean_cluster_size": 2.5,
                "diversity": 10,
                "norm_shannon_wiener": 0.8,
            }
        ]
    )

    def fake_plot_stats(data, **kwargs):
        captured["data"] = data
        captured["kwargs"] = kwargs
        return "cluster_plot"

    monkeypatch.setattr(rsplot, "plot_stats", fake_plot_stats)

    result = rsplot.cluster_properties(stats_df)

    assert result == "cluster_plot"
    assert captured["data"] is stats_df
    assert captured["kwargs"]["properties"] == [
        "mean_cluster_size",
        "diversity",
        "norm_shannon_wiener",
    ]
    assert captured["kwargs"]["zero_bottom"] is True
