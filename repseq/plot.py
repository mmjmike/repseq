"""Plotting helpers for repseq statistics tables.

The functions in this module take wide statistics tables produced by
``repseq.stats`` and draw matplotlib/seaborn categorical summaries.  The
plotting API is intentionally dataframe-oriented so that sample metadata can be
joined immediately before plotting.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage


CDR3AA_STATS_PROPERTIES = [
    "mean_cdr3nt_len",
    "mean_insert_size",
    "cdr3_5_charge",
    "cdr3_5_kf4",
    "cdr3_5_volume",
    "cdr3_5_strength",
]

DIVERSITY_STATS_PROPERTIES = [
    "diversity",
    "norm_shannon_wiener",
    "chao1",
    "berger_parker",
    "d50",
]

CONVERGENCE_PROPERTIES = [
    "convergence",
    "convergence_v",
    "convergence_vj",
]

PROCESSING_PROPERTIES = [
    "reads_aligned_pc",
    "reads_per_umi",
    "clones_func",
    "umi_in_func_clones",
]



RAREFACTION_COLORS_20 = [
    "#4e79a7", "#a0cde8", "#f28e2b", "#ffbe7d", "#59a14f",
    "#8cd17d", "#b6992d", "#f1ce63", "#499894", "#86bcb6",
    "#e15759", "#ff9d9a", "#79706e", "#bab0ac", "#d37295",
    "#fabfd2", "#b07aa1", "#d4a6c8", "#9d7660", "#d7b5a6",
]

PHEATMAP_COLORS_7 = [
    "#4575B4",
    "#91BFDB",
    "#E0F3F8",
    "#FFFFBF",
    "#FEE090",
    "#FC8D59",
    "#D73027",
]
PHEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "pheatmap_default",
    PHEATMAP_COLORS_7,
    N=100,
)

GENE_RE = re.compile(
    r"""
    ^
    (?P<system>TR|IG)
    (?P<chain>[ABGDHKL])
    (?P<head>[A-Z])
    (?P<body>[A-Z0-9/-]*)
    (?:\*(?P<allele>[^(),;|\s]+))?
    $
    """,
    re.VERBOSE | re.IGNORECASE,
)


BODY_RE = re.compile(
    r"""
    ^
    (?P<family>\d+)
    (?P<subfamily>[A-Z]+)?

    (?:
        -?
        (?P<segment>\d+)
        (?P<attached_dual>[A-Z]+\d+)?
        (?:-(?P<subsegment>[A-Z0-9]+))?
    )?

    (?:/(?P<slash_dual>[A-Z0-9-]+))?
    $
    """,
    re.VERBOSE | re.IGNORECASE,
)


def parse_gene_name(gene):
    """Parse a V, D, J, or C gene name into sortable components."""
    if not isinstance(gene, str):
        return None
    gene_original = gene
    gene = gene.strip().upper()

    # AIRR/MiXCR multi-calls are ordered by preference; plot the first call.
    gene = gene.split(",", 1)[0]
    gene = gene.split(";", 1)[0]
    gene = gene.split("|", 1)[0]
    gene = gene.split("(", 1)[0]

    match = GENE_RE.match(gene)
    if not match:
        return None

    parsed = match.groupdict()
    system = parsed["system"].upper()
    chain = parsed["chain"].upper()
    head = parsed["head"].upper()
    body = parsed["body"].upper()

    # IGHD is the delta constant gene, while IGHD3-10 is a D segment.
    if system == "IG" and chain == "H" and (
        head in {"M", "G", "A", "E"} or (head == "D" and body == "")
    ):
        gene_type = "C"
        isotype = head
    else:
        gene_type = head
        isotype = None

    family = None
    subfamily = None
    segment = None
    subsegment = None
    dual_designation = None
    if body:
        body_match = BODY_RE.match(body)
        if body_match:
            body_parts = body_match.groupdict()
            family = body_parts["family"]
            subfamily = body_parts["subfamily"]
            segment = body_parts["segment"]
            subsegment = body_parts["subsegment"]
            dual_designation = (
                body_parts["attached_dual"] or body_parts["slash_dual"]
            )

    return {
        "original": gene_original,
        "system": system,
        "chain": chain,
        "gene_type": gene_type,
        "isotype": isotype,
        "family": family,
        "subfamily": subfamily,
        "segment": segment,
        "subsegment": subsegment,
        "dual_designation": dual_designation,
        "allele": parsed["allele"],
    }


def _natural_sort_key(value):
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in re.split(r"(\d+)", str(value))
        if part
    )


def _optional_number(value):
    return (1, 0) if value is None else (0, int(value))


def _optional_text(value):
    return (1, "") if value is None else (0, value.casefold())


def _gene_sort_key(gene):
    parsed = parse_gene_name(gene)
    if parsed is None:
        return (1, _natural_sort_key(gene))
    return (
        0,
        parsed["system"],
        parsed["chain"],
        parsed["gene_type"],
        parsed["isotype"] or "",
        _optional_number(parsed["family"]),
        parsed["subfamily"] or "",
        _optional_number(parsed["segment"]),
        _natural_sort_key(parsed["dual_designation"] or ""),
        _natural_sort_key(parsed["subsegment"] or ""),
        _optional_text(parsed["allele"]),
        _natural_sort_key(gene),
    )


def _sort_gene_names(genes):
    return sorted(pd.unique(pd.Series(genes).dropna().astype(str)), key=_gene_sort_key)


def _as_list(value, name, max_len=None):
    if value is None:
        values = []
    elif isinstance(value, str):
        values = [value]
    elif isinstance(value, Iterable):
        values = list(value)
    else:
        raise TypeError(f"{name} must be None, a string, or a list/tuple of strings")

    if max_len is not None and len(values) > max_len:
        raise ValueError(f"{name} can contain at most {max_len} column(s)")

    if any(not isinstance(item, str) for item in values):
        raise TypeError(f"{name} must contain only column names as strings")

    return values


def _caption(value):
    return " ".join(part.capitalize() for part in str(value).split("_"))


def _category_order(series):
    if isinstance(series.dtype, pd.CategoricalDtype):
        return [item for item in series.cat.categories if item in set(series.dropna())]
    return list(pd.unique(series.dropna()))


def _interaction_order(data, split_columns):
    if all(isinstance(data[column].dtype, pd.CategoricalDtype) for column in split_columns):
        levels = [_category_order(data[column]) for column in split_columns]
        observed = set(zip(*(data[column] for column in split_columns)))
        return [
            " | ".join(map(str, combination))
            for combination in _cartesian_product(levels)
            if combination in observed
        ]
    return list(pd.unique(data["_split_panel"].dropna()))


def _cartesian_product(levels):
    if not levels:
        return [()]
    result = [()]
    for level in levels:
        result = [prefix + (item,) for prefix in result for item in level]
    return result


def _merge_stats_metadata(stats_df, metadata):
    if "sample_id" not in stats_df.columns:
        raise ValueError("stats_df must contain a 'sample_id' column")

    if metadata is None:
        merge_keys = ["sample_id"]
        if "chain" in stats_df.columns:
            merge_keys.append("chain")
        data = stats_df.copy()
        for column in data.select_dtypes(include="category").columns:
            data[column] = data[column].cat.remove_unused_categories()
        return data, set(data.columns), merge_keys

    if "sample_id" not in metadata.columns:
        raise ValueError("metadata must contain a 'sample_id' column")

    merge_keys = ["sample_id"]
    if "chain" in stats_df.columns and "chain" in metadata.columns:
        merge_keys.append("chain")

    duplicated = metadata.duplicated(merge_keys, keep=False)
    if duplicated.any():
        duplicated_keys = metadata.loc[duplicated, merge_keys].drop_duplicates()
        raise ValueError(
            "metadata must contain one row per merge key. Duplicated keys: "
            f"{duplicated_keys.to_dict(orient='records')}"
        )

    stats_keys = stats_df[merge_keys].drop_duplicates()
    merged = stats_df.merge(metadata, on=merge_keys, how="left", indicator=True)
    matched_keys = merged.loc[merged["_merge"] == "both", merge_keys].drop_duplicates()

    if len(matched_keys) < len(stats_keys):
        warnings.warn(
            "Metadata contains "
            f"{len(matched_keys)} of {len(stats_keys)} samples from stats_df. "
            "Plotting a subset of samples.",
            UserWarning,
            stacklevel=2,
        )
        merged = merged.loc[merged["_merge"] == "both"].copy()
        if merged.empty:
            raise ValueError("metadata does not match any samples in stats_df")

    merged = merged.drop(columns="_merge")
    for column in merged.select_dtypes(include="category").columns:
        merged[column] = merged[column].cat.remove_unused_categories()
    return merged, set(merged.columns), merge_keys


def _validate_metadata_columns(columns, available_columns, label):
    missing = [column for column in columns if column not in available_columns]
    if missing:
        raise ValueError(
            f"{label} column(s) must be present in metadata or stats_df: {missing}"
        )


def _make_sample_labels(data):
    chain_column = next(
        (column for column in ["extracted_chain", "chain"] if column in data.columns),
        None,
    )
    use_chain = (
        chain_column is not None
        and data["sample_id"].duplicated(keep=False).any()
    )
    if use_chain:
        labels = (
            data["sample_id"].astype(str)
            + " ("
            + data[chain_column].astype(str)
            + ")"
        )
    else:
        labels = data["sample_id"].astype(str)
    order = list(pd.unique(labels))
    return pd.Categorical(labels, categories=order, ordered=True)


def _prepare_plot_data(stats_df, metadata, properties, group, split):
    group_columns = _as_list(group, "group", max_len=2)
    split_columns = _as_list(split, "split", max_len=2)
    if len(group_columns) == 0:
        group_columns = []

    data, metadata_columns, merge_keys = _merge_stats_metadata(stats_df, metadata)
    _validate_metadata_columns(group_columns, metadata_columns, "group")
    _validate_metadata_columns(split_columns, metadata_columns, "split")

    missing_properties = [column for column in properties if column not in data.columns]
    if missing_properties:
        raise ValueError(f"stats_df does not contain property column(s): {missing_properties}")

    label_columns = [
        column for column in ["extracted_chain", "chain"] if column in data.columns
    ]
    id_columns = list(
        dict.fromkeys(merge_keys + label_columns + group_columns + split_columns)
    )
    plot_data = data[id_columns + list(properties)].copy()
    plot_data["_sample_label"] = _make_sample_labels(plot_data)

    long_data = plot_data.melt(
        id_vars=id_columns + ["_sample_label"],
        value_vars=list(properties),
        var_name="property",
        value_name="value",
    )
    long_data["property_label"] = pd.Categorical(
        long_data["property"].map(_caption),
        categories=[_caption(property_name) for property_name in properties],
        ordered=True,
    )

    panel_column = None
    if len(split_columns) == 1:
        panel_column = split_columns[0]
    elif len(split_columns) == 2:
        panel_column = "_split_panel"
        long_data[panel_column] = (
            long_data[split_columns[0]].astype(str)
            + " | "
            + long_data[split_columns[1]].astype(str)
        )
        long_data[panel_column] = pd.Categorical(
            long_data[panel_column],
            categories=_interaction_order(long_data, split_columns),
            ordered=True,
        )

    return long_data, group_columns, panel_column


def _facet_layout(properties, panel_column):
    many_properties = len(properties) >= 2
    if many_properties and panel_column is None:
        return None, "property_label", 3
    if many_properties:
        return "property_label", panel_column, None
    if panel_column is not None:
        return None, panel_column, 3
    return None, None, None


def _draw_category_panel(
    data,
    x_column,
    hue_column,
    palette,
    x_order=None,
    hue_order=None,
    **kwargs,
):
    ax = kwargs.get("ax", plt.gca())
    order = x_order or _category_order(data[x_column])
    hue_order = hue_order or (_category_order(data[hue_column]) if hue_column else None)

    if hue_column is None and x_column == "_sample_label":
        sns.barplot(
            data=data,
            x=x_column,
            y="value",
            order=order,
            errorbar=None,
            ax=ax,
        )
        ax.tick_params(axis="x", rotation=90)
        return

    if hue_column == x_column:
        sns.boxplot(
            data=data,
            x=x_column,
            y="value",
            hue=hue_column,
            order=order,
            hue_order=order,
            palette=palette,
            dodge=False,
            showfliers=False,
            ax=ax,
        )
        sns.stripplot(
            data=data,
            x=x_column,
            y="value",
            hue=hue_column,
            order=order,
            hue_order=order,
            palette=palette,
            dodge=False,
            jitter=0.15,
            linewidth=0.3,
            edgecolor="black",
            alpha=0.75,
            ax=ax,
        )
        if ax.legend_ is not None:
            ax.legend_.remove()
        return

    sns.boxplot(
        data=data,
        x=x_column,
        y="value",
        hue=hue_column,
        order=order,
        hue_order=hue_order,
        palette=palette,
        dodge=True,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=data,
        x=x_column,
        y="value",
        hue=hue_column,
        order=order,
        hue_order=hue_order,
        palette=palette,
        dodge=True,
        jitter=0.15,
        linewidth=0.3,
        edgecolor="black",
        alpha=0.75,
        ax=ax,
    )


def _deduplicate_legend(grid, hue_column):
    if hue_column is None:
        return
    handles = []
    labels = []
    for ax in grid.axes.flat:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        if ax.legend_ is not None:
            ax.legend_.remove()
        for handle, label in zip(ax_handles, ax_labels):
            if label and label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        grid.fig.legend(
            handles,
            labels,
            title=_caption(hue_column),
            loc="upper center",
            ncol=min(max(1, len(labels)), 5),
            frameon=False,
        )
        grid.fig.subplots_adjust(top=0.9)


def plot_stats(
    stats_df,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
    zero_bottom=False,
):
    """Plot a wide statistics table with optional sample metadata.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        Statistics table. It must contain ``sample_id`` and the selected
        property columns. If both ``stats_df`` and ``metadata`` contain
        ``chain``, the merge uses both ``sample_id`` and ``chain``. Grouping
        and splitting columns may also be included directly in this table.
    metadata : pandas.DataFrame, optional
        Sample metadata with one row per ``sample_id`` or ``sample_id`` +
        ``chain`` key. If metadata covers only part of ``stats_df``, a warning
        is emitted and only matched samples are plotted.
    properties : str or sequence of str
        Numeric columns to plot. Multiple properties are shown as separate
        panels.
    group : str or sequence of str, optional
        Column(s) used to group samples. They may come from ``metadata`` or,
        when metadata is omitted, directly from ``stats_df``. No group draws
        one bar per sample. One group draws boxplots with jittered points by
        that group. Two groups use the first column on the x-axis and the
        second column as color.
    split : str or sequence of str, optional
        One or two columns used to split plots into panels. They may come from
        ``metadata`` or, when metadata is omitted, directly from ``stats_df``.
        Two split columns are combined into an interaction panel.
    palette : optional
        Any seaborn/matplotlib-compatible palette specification.
    height, aspect : float
        Facet size arguments passed to :class:`seaborn.FacetGrid`.
    zero_bottom : bool
        If ``True``, set the lower y-axis limit of every panel to zero.

    Returns
    -------
    seaborn.FacetGrid
        The created grid.
    """
    if properties is None:
        raise ValueError("properties must be specified for plot_stats")

    properties = _as_list(properties, "properties")
    if not properties:
        raise ValueError("properties must contain at least one column")

    plot_data, group_columns, panel_column = _prepare_plot_data(
        stats_df, metadata, properties, group, split
    )

    x_column = "_sample_label" if not group_columns else group_columns[0]
    hue_column = None
    x_order = None
    hue_order = None
    if len(group_columns) == 1:
        hue_column = group_columns[0]
        x_order = _category_order(plot_data[x_column])
        hue_order = x_order
    elif len(group_columns) == 2:
        hue_column = group_columns[1]
        x_order = _category_order(plot_data[x_column])
        hue_order = _category_order(plot_data[hue_column])

    row, col, col_wrap = _facet_layout(properties, panel_column)
    grid = sns.FacetGrid(
        plot_data,
        row=row,
        col=col,
        col_wrap=col_wrap,
        sharex=not (panel_column is not None and not group_columns),
        sharey=False,
        height=height,
        aspect=aspect,
        despine=True,
    )
    grid.map_dataframe(
        _draw_category_panel,
        x_column=x_column,
        hue_column=hue_column,
        palette=palette,
        x_order=x_order,
        hue_order=hue_order,
    )

    x_label = "Sample" if x_column == "_sample_label" else _caption(x_column)
    y_label = _caption(properties[0]) if len(properties) == 1 else "Value"
    grid.set_axis_labels(x_label, y_label)

    if row is not None and col is not None:
        grid.set_titles(row_template="{row_name}", col_template="{col_name}")
    elif col is not None:
        grid.set_titles("{col_name}")
    elif row is not None:
        grid.set_titles("{row_name}")
    else:
        grid.axes.flat[0].set_title(_caption(properties[0]))

    _deduplicate_legend(grid, hue_column if len(group_columns) == 2 else None)
    if zero_bottom:
        for ax in grid.axes.flat:
            ax.set_ylim(bottom=0)
    grid.tight_layout()
    return grid


def _default_grouped_clonoset_properties(stats_df):
    properties = ["reads"]
    if "reads_per_umi" in stats_df.columns:
        properties.append("reads_per_umi")
    properties.append("clones_func")
    if "umi_func" in stats_df.columns and stats_df["umi_func"].notna().any():
        properties.append("umi_func")
    return properties


def _format_clonoset_bar_value(value):
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def _prepare_default_clonoset_plot_data(stats_df, metadata, split):
    split_columns = _as_list(split, "split", max_len=2)
    data, metadata_columns, merge_keys = _merge_stats_metadata(stats_df, metadata)
    _validate_metadata_columns(split_columns, metadata_columns, "split")

    descriptors = [
        ("reads", "reads_func", "Reads"),
    ]
    if "reads_per_umi" in data.columns:
        descriptors.append(("reads_per_umi", None, "Reads Per Umi"))
    descriptors.append(("clones", "clones_func", "Clones"))
    if "umi" in data.columns and data["umi"].notna().any():
        descriptors.append(("umi", "umi_func", "Umi"))

    required_columns = {
        column
        for total_column, functional_column, _ in descriptors
        for column in (total_column, functional_column)
        if column is not None
    }
    missing_columns = sorted(required_columns.difference(data.columns))
    if missing_columns:
        raise ValueError(
            f"stats_df does not contain clonoset-stat column(s): {missing_columns}"
        )

    id_columns = list(dict.fromkeys(merge_keys + split_columns))
    sample_labels = _make_sample_labels(data)
    sample_order = list(sample_labels.categories)
    plot_parts = []
    for total_column, functional_column, property_label in descriptors:
        columns = id_columns + [total_column]
        if functional_column is not None:
            columns.append(functional_column)
        part = data[columns].copy()
        part["_sample_label"] = sample_labels
        part["property_label"] = property_label
        part["total"] = pd.to_numeric(part[total_column], errors="coerce")
        part["functional"] = (
            pd.to_numeric(part[functional_column], errors="coerce")
            if functional_column is not None
            else np.nan
        )
        part["_paired"] = functional_column is not None
        if functional_column is not None:
            part = part.loc[part["total"].notna() & part["functional"].notna()]
        plot_parts.append(
            part[
                id_columns
                + [
                    "_sample_label",
                    "property_label",
                    "total",
                    "functional",
                    "_paired",
                ]
            ]
        )

    plot_data = pd.concat(plot_parts, ignore_index=True)
    property_order = [descriptor[2] for descriptor in descriptors]
    plot_data["_sample_label"] = pd.Categorical(
        plot_data["_sample_label"], categories=sample_order, ordered=True
    )
    plot_data["property_label"] = pd.Categorical(
        plot_data["property_label"], categories=property_order, ordered=True
    )

    panel_column = None
    if len(split_columns) == 1:
        panel_column = split_columns[0]
    elif len(split_columns) == 2:
        panel_column = "_split_panel"
        plot_data[panel_column] = (
            plot_data[split_columns[0]].astype(str)
            + " | "
            + plot_data[split_columns[1]].astype(str)
        )
        plot_data[panel_column] = pd.Categorical(
            plot_data[panel_column],
            categories=_interaction_order(plot_data, split_columns),
            ordered=True,
        )

    return plot_data, property_order, panel_column


def _draw_default_clonoset_panel(data, colors, **kwargs):
    ax = kwargs.get("ax", plt.gca())
    sample_order = _category_order(data["_sample_label"])
    paired = bool(data["_paired"].iloc[0])

    if not paired:
        sns.barplot(
            data=data,
            x="_sample_label",
            y="total",
            order=sample_order,
            errorbar=None,
            color=colors["Total"],
            ax=ax,
        )
        ax.tick_params(axis="x", rotation=90)
        return

    values = (
        data.groupby("_sample_label", observed=True)[["total", "functional"]]
        .first()
        .reindex(sample_order)
    )
    x = np.arange(len(sample_order), dtype=float)
    ax.bar(
        x,
        values["total"].to_numpy(),
        width=0.8,
        color=colors["Total"],
        label="Total",
    )
    ax.bar(
        x,
        values["functional"].to_numpy(),
        width=0.8,
        color=colors["Functional"],
        label="Functional",
    )
    for position, total, functional in zip(
        x, values["total"], values["functional"]
    ):
        label = (
            f"{_format_clonoset_bar_value(total)}"
            f"({_format_clonoset_bar_value(functional)})"
        )
        ax.annotate(
            label,
            xy=(position, total),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=90)
    maximum = values["total"].max()
    if pd.notna(maximum) and maximum >= 0:
        ax.set_ylim(top=max(1, maximum * 1.15))


def _plot_default_clonoset_stats(
    stats_df,
    metadata,
    split,
    palette,
    height,
    aspect,
):
    plot_data, property_order, panel_column = _prepare_default_clonoset_plot_data(
        stats_df, metadata, split
    )
    if palette is None:
        shades = sns.color_palette("Blues", n_colors=4)
        color_map = {"Total": shades[1], "Functional": shades[3]}
    else:
        color_map = _palette_mapping(["Total", "Functional"], palette)
    row, col, col_wrap = _facet_layout(property_order, panel_column)
    grid = sns.FacetGrid(
        plot_data,
        row=row,
        col=col,
        col_wrap=col_wrap,
        sharex=False,
        sharey=False,
        height=height,
        aspect=aspect,
        despine=True,
    )
    grid.map_dataframe(_draw_default_clonoset_panel, colors=color_map)
    grid.set_axis_labels("Sample", "Value")

    if row is not None and col is not None:
        grid.set_titles(row_template="{row_name}", col_template="{col_name}")
    elif col is not None:
        grid.set_titles("{col_name}")
    elif row is not None:
        grid.set_titles("{row_name}")
    else:
        grid.axes.flat[0].set_title(property_order[0])

    handles = [
        Patch(facecolor=color_map[level], label=level)
        for level in ("Total", "Functional")
    ]
    grid.fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
    )
    grid.fig.tight_layout(rect=(0, 0, 1, 0.92))
    return grid


def _column_by_name(columns, names):
    names = {name.casefold() for name in names}
    matches = [
        column
        for column in columns
        if isinstance(column, str) and column.casefold() in names
    ]
    return matches[0] if matches else None


def _infer_segment_type(values):
    gene_types = {
        parsed["gene_type"].lower()
        for parsed in (parse_gene_name(value) for value in values)
        if parsed is not None
    }
    if len(gene_types) != 1 or not gene_types.issubset({"v", "j", "c"}):
        raise ValueError(
            "Could not detect one V, J, or C segment type from the gene names"
        )
    return gene_types.pop()


def _segment_chain(gene):
    parsed = parse_gene_name(gene)
    if parsed is None:
        return None
    return parsed["system"] + parsed["chain"]


def _segment_family(gene, segment_type):
    if str(gene).strip() in {".", "NA"}:
        return "NA"
    parsed = parse_gene_name(gene)
    if parsed is None or parsed["gene_type"].casefold() != segment_type:
        raise ValueError(f"Could not determine the segment family for {gene!r}")
    if segment_type == "c":
        return parsed["system"] + parsed["chain"] + (
            parsed["isotype"] or parsed["gene_type"]
        )
    if parsed["family"] is None:
        raise ValueError(f"Could not determine the segment family for {gene!r}")
    return (
        parsed["system"]
        + parsed["chain"]
        + parsed["gene_type"]
        + parsed["family"]
    )


def _combine_segment_family_usage(data, segment_type):
    data = data.copy()
    data["_segment"] = data["_segment"].map(
        lambda gene: _segment_family(gene, segment_type)
    )
    return (
        data.groupby(["sample_id", "chain", "_segment"], as_index=False, sort=False)[
            "_value"
        ]
        .sum()
    )


def _normalize_segment_usage_table(segment_usage_df):
    if not isinstance(segment_usage_df, pd.DataFrame):
        raise TypeError("segment_usage_df must be a pandas DataFrame")
    if "sample_id" not in segment_usage_df.columns:
        raise ValueError("segment_usage_df must contain a 'sample_id' column")

    data = segment_usage_df.copy()
    segment_columns = [
        column
        for column in data.columns
        if isinstance(column, str) and column.casefold() in {"v", "j", "c"}
    ]
    generic_segment_column = _column_by_name(data.columns, {"segment", "gene"})
    value_column = _column_by_name(
        data.columns, {"usage", "value", "freq", "frequency", "count"}
    )

    if segment_columns:
        if len(segment_columns) != 1:
            raise ValueError(
                "Long segment-usage tables must contain exactly one of: v, j, c"
            )
        if value_column is None:
            raise ValueError(
                "Long segment-usage tables need a usage, value, freq, frequency, "
                "or count column"
            )
        segment_column = segment_columns[0]
        segment_type = segment_column.casefold()
        keep_columns = ["sample_id", segment_column, value_column]
        if "chain" in data.columns:
            keep_columns.insert(1, "chain")
        data = data[keep_columns].rename(
            columns={segment_column: "_segment", value_column: "_value"}
        )
    elif generic_segment_column is not None:
        if value_column is None:
            raise ValueError(
                "Long segment-usage tables need a usage, value, freq, frequency, "
                "or count column"
            )
        segment_type = _infer_segment_type(data[generic_segment_column].dropna())
        keep_columns = ["sample_id", generic_segment_column, value_column]
        if "chain" in data.columns:
            keep_columns.insert(1, "chain")
        data = data[keep_columns].rename(
            columns={generic_segment_column: "_segment", value_column: "_value"}
        )
    else:
        id_columns = ["sample_id"] + (["chain"] if "chain" in data.columns else [])
        wide_segment_columns = [
            column
            for column in data.columns
            if column not in id_columns
            and (
                str(column).strip() == "."
                or parse_gene_name(str(column)) is not None
            )
        ]
        if not wide_segment_columns:
            raise ValueError(
                "Could not detect a long or wide V, J, or C segment-usage table"
            )
        segment_type = _infer_segment_type(map(str, wide_segment_columns))
        data = data.melt(
            id_vars=id_columns,
            value_vars=wide_segment_columns,
            var_name="_segment",
            value_name="_value",
        )

    data["_value"] = pd.to_numeric(data["_value"], errors="coerce")
    invalid_values = data["_value"].isna() | data["_segment"].isna()
    if invalid_values.any():
        warnings.warn(
            f"Dropped {int(invalid_values.sum())} segment-usage row(s) with "
            "missing or non-numeric values.",
            UserWarning,
            stacklevel=2,
        )
        data = data.loc[~invalid_values].copy()
    if data.empty:
        raise ValueError("No valid segment-usage values were found")

    data["_segment"] = data["_segment"].astype(str).replace({".": "NA"})
    inferred_chains = data["_segment"].map(_segment_chain)
    if "chain" not in data.columns:
        data["chain"] = inferred_chains.fillna("Unknown")
    else:
        data["chain"] = data["chain"].where(data["chain"].notna(), inferred_chains)
        data["chain"] = data["chain"].fillna("Unknown")
        # Wide batch tables contain the union of genes from every chain. Keep
        # each chain panel independent by removing those cross-chain zero-fill
        # columns after melting.
        cross_chain = inferred_chains.notna() & (
            data["chain"].astype(str) != inferred_chains.astype(str)
        )
        data = data.loc[~cross_chain].copy()
    return data, segment_type


def _prepare_segment_usage_data(
    segment_usage_df,
    metadata,
    group,
    split,
    plot_type,
    combine_families,
):
    data, segment_type = _normalize_segment_usage_table(segment_usage_df)
    if combine_families:
        data = _combine_segment_family_usage(data, segment_type)
    max_groups = 3 if plot_type == "heatmap" else 1
    group_columns = _as_list(group, "group", max_len=max_groups)
    split_columns = _as_list(split, "split", max_len=1)

    data, metadata_columns, _ = _merge_stats_metadata(data, metadata)
    if metadata is None and (group_columns or split_columns):
        raise ValueError("metadata is required when group or split columns are used")
    _validate_metadata_columns(group_columns, metadata_columns, "group")
    _validate_metadata_columns(split_columns, metadata_columns, "split")
    data["_sample_label"] = data["sample_id"].astype(str)
    return data, segment_type, group_columns, split_columns


def _panel_subsets(data, split_columns):
    chain_order = _category_order(data["chain"])
    if not split_columns:
        return [
            (str(chain), data.loc[data["chain"] == chain].copy())
            for chain in chain_order
        ]

    split_column = split_columns[0]
    split_order = _category_order(data[split_column])
    panels = []
    for chain in chain_order:
        chain_data = data.loc[data["chain"] == chain]
        for split_value in split_order:
            panel_data = chain_data.loc[chain_data[split_column] == split_value].copy()
            if not panel_data.empty:
                panels.append((f"{chain} | {split_value}", panel_data))
    return panels


def _palette_mapping(levels, palette=None):
    if isinstance(palette, dict):
        missing = [level for level in levels if level not in palette]
        if missing:
            raise ValueError(f"palette does not define colors for: {missing}")
        return {level: palette[level] for level in levels}
    colors = sns.color_palette(palette, n_colors=len(levels))
    return dict(zip(levels, colors))


def _add_figure_legend(fig, axes, title):
    handles = []
    labels = []
    for ax in axes:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        if ax.legend_ is not None:
            ax.legend_.remove()
        for handle, label in zip(ax_handles, ax_labels):
            if label and label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        fig.legend(
            handles,
            labels,
            title=title,
            loc="upper center",
            ncol=min(5, len(labels)),
            frameon=False,
        )


def _draw_segment_barplot(ax, data, group_column, palette):
    segment_order = _sort_gene_names(data["_segment"])
    x = np.arange(len(segment_order), dtype=float)

    if group_column is None:
        series_column = "_sample_label"
        series_order = _category_order(data[series_column])
        values = data.pivot_table(
            index="_segment",
            columns=series_column,
            values="_value",
            aggfunc="sum",
            fill_value=0,
            sort=False,
        ).reindex(segment_order, fill_value=0)
        color_map = _palette_mapping(series_order, palette)
        width = 0.8 / max(1, len(series_order))
        for index, series in enumerate(series_order):
            offset = (index - (len(series_order) - 1) / 2) * width
            ax.bar(
                x + offset,
                values.reindex(columns=series_order)[series].to_numpy(),
                width=width,
                color=color_map[series],
                label=str(series),
            )
    else:
        series_order = _category_order(data[group_column])
        color_map = _palette_mapping(series_order, palette)
        summary = (
            data.groupby(["_segment", group_column], observed=True)["_value"]
            .agg(["mean", "std"])
        )
        width = 0.8 / max(1, len(series_order))
        for index, series in enumerate(series_order):
            group_summary = summary.xs(series, level=group_column).reindex(segment_order)
            offset = (index - (len(series_order) - 1) / 2) * width
            ax.bar(
                x + offset,
                group_summary["mean"].fillna(0).to_numpy(),
                yerr=group_summary["std"].fillna(0).to_numpy(),
                width=width,
                capsize=2,
                color=color_map[series],
                label=str(series),
            )

    ax.set_xticks(x)
    ax.set_xticklabels(segment_order, rotation=90)


def _draw_segment_boxplot(ax, data, group_column, palette, seed):
    segment_order = _sort_gene_names(data["_segment"])
    group_order = _category_order(data[group_column])
    color_map = _palette_mapping(group_order, palette)
    sns.boxplot(
        data=data,
        x="_segment",
        y="_value",
        hue=group_column,
        order=segment_order,
        hue_order=group_order,
        palette=color_map,
        showfliers=False,
        ax=ax,
    )

    random_state = np.random.get_state()
    try:
        np.random.seed(seed)
        sns.stripplot(
            data=data,
            x="_segment",
            y="_value",
            hue=group_column,
            order=segment_order,
            hue_order=group_order,
            palette=color_map,
            dodge=True,
            jitter=0.15,
            linewidth=0.3,
            edgecolor="black",
            alpha=0.75,
            ax=ax,
        )
    finally:
        np.random.set_state(random_state)
    ax.tick_params(axis="x", rotation=90)


def _categorical_segment_usage_plot(
    data,
    panels,
    plot_type,
    group_columns,
    palette,
    height,
    aspect,
    seed,
    segment_label,
):
    group_column = group_columns[0] if group_columns else None
    series_column = group_column or "_sample_label"
    series_count = data[series_column].nunique(dropna=True)
    if series_count > 10:
        warnings.warn(
            "Categorical segment-usage plots support at most 10 groups; "
            f"found {series_count}. Nothing was plotted.",
            UserWarning,
            stacklevel=2,
        )
        return None

    max_segments = max(panel["_segment"].nunique() for _, panel in panels)
    figure_width = max(8, min(24, max_segments * 0.38 * aspect))
    fig, axes_array = plt.subplots(
        len(panels),
        1,
        figsize=(figure_width, height * len(panels)),
        squeeze=False,
    )
    axes = list(axes_array[:, 0])
    for panel_index, ((title, panel_data), ax) in enumerate(zip(panels, axes)):
        if plot_type == "barplot":
            _draw_segment_barplot(ax, panel_data, group_column, palette)
        else:
            _draw_segment_boxplot(
                ax, panel_data, group_column, palette, seed + panel_index
            )
        ax.set_title(title)
        ax.set_xlabel(segment_label)
        ax.set_ylabel("Value")
        sns.despine(ax=ax)

    _add_figure_legend(
        fig,
        axes,
        _caption(group_column) if group_column else "Sample",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def _cluster_rows(matrix):
    if len(matrix) < 2:
        return matrix, None
    linkage_matrix = linkage(matrix.to_numpy(), method="average", metric="euclidean")
    order = leaves_list(linkage_matrix)
    return matrix.iloc[order], linkage_matrix


def _heatmap_annotation_colors(data, sample_order, group_columns, palette):
    sample_metadata = (
        data[["_sample_label"] + group_columns]
        .drop_duplicates("_sample_label")
        .set_index("_sample_label")
        .reindex(sample_order)
    )
    rgb = np.zeros((len(sample_order), len(group_columns), 3))
    legend_handles = []
    for group_index, group_column in enumerate(group_columns):
        levels = _category_order(data[group_column])
        color_map = _palette_mapping(levels, palette)
        for sample_index, value in enumerate(sample_metadata[group_column]):
            rgb[sample_index, group_index] = to_rgb(
                color_map.get(value, "#d0d0d0")
            )
        legend_handles.extend(
            Patch(
                facecolor=color_map[level],
                label=str(level),
            )
            for level in levels
        )
    return rgb, legend_handles


def _heatmap_segment_usage_plot(
    panels,
    group_columns,
    palette,
    cmap,
    height,
    aspect,
    segment_label,
):
    max_segments = max(panel["_segment"].nunique() for _, panel in panels)
    max_samples = max(panel["_sample_label"].nunique() for _, panel in panels)
    figure_width = max(9, min(28, max_segments * 0.35 * aspect + 3))
    panel_height = max(height, min(12, max_samples * 0.28 + 1.5))
    fig = plt.figure(figsize=(figure_width, panel_height * len(panels)))
    outer_grid = fig.add_gridspec(len(panels), 1, hspace=0.55)
    all_legend_handles = []

    for panel_index, (title, panel_data) in enumerate(panels):
        segment_order = _sort_gene_names(panel_data["_segment"])
        sample_order = list(pd.unique(panel_data["_sample_label"]))
        matrix = panel_data.pivot_table(
            index="_sample_label",
            columns="_segment",
            values="_value",
            aggfunc="sum",
            fill_value=0,
            sort=False,
        ).reindex(index=sample_order, columns=segment_order, fill_value=0)

        linkage_matrix = None
        if group_columns:
            matrix, linkage_matrix = _cluster_rows(matrix)
            inner_grid = outer_grid[panel_index].subgridspec(
                1,
                4,
                width_ratios=[0.8, max(0.35, 0.25 * len(group_columns)), 6, 0.18],
                wspace=0.08,
            )
            dendrogram_ax = fig.add_subplot(inner_grid[0, 0])
            annotation_ax = fig.add_subplot(inner_grid[0, 1])
            heatmap_ax = fig.add_subplot(inner_grid[0, 2])
            colorbar_ax = fig.add_subplot(inner_grid[0, 3])

            if linkage_matrix is not None:
                dendrogram(
                    linkage_matrix,
                    orientation="left",
                    no_labels=True,
                    color_threshold=0,
                    above_threshold_color="#555555",
                    ax=dendrogram_ax,
                )
                dendrogram_ax.invert_yaxis()
            dendrogram_ax.axis("off")

            annotation_rgb, legend_handles = _heatmap_annotation_colors(
                panel_data,
                list(matrix.index),
                group_columns,
                palette,
            )
            annotation_ax.imshow(
                annotation_rgb,
                aspect="auto",
                interpolation="nearest",
            )
            annotation_ax.set_xticks(np.arange(len(group_columns)))
            annotation_ax.set_xticklabels(
                [_caption(column) for column in group_columns],
                rotation=90,
            )
            annotation_ax.xaxis.tick_top()
            annotation_ax.set_yticks([])
            annotation_ax.tick_params(length=0)
            for spine in annotation_ax.spines.values():
                spine.set_visible(False)
            all_legend_handles.extend(legend_handles)
        else:
            inner_grid = outer_grid[panel_index].subgridspec(
                1, 2, width_ratios=[6, 0.18], wspace=0.08
            )
            heatmap_ax = fig.add_subplot(inner_grid[0, 0])
            colorbar_ax = fig.add_subplot(inner_grid[0, 1])

        sns.heatmap(
            matrix,
            cmap=cmap,
            ax=heatmap_ax,
            cbar=True,
            cbar_ax=colorbar_ax,
            cbar_kws={"label": "Value"},
            xticklabels=True,
            yticklabels=True,
        )
        heatmap_ax.set_title(title)
        heatmap_ax.set_xlabel(segment_label)
        heatmap_ax.set_ylabel("Sample")
        heatmap_ax.tick_params(axis="x", labelrotation=90)

    if all_legend_handles:
        unique_handles = {}
        for handle in all_legend_handles:
            unique_handles.setdefault(handle.get_label(), handle)
        fig.legend(
            unique_handles.values(),
            unique_handles.keys(),
            loc="upper center",
            ncol=min(5, len(unique_handles)),
            frameon=False,
        )
        fig.subplots_adjust(top=0.93)
    return fig


def segment_usage(
    segment_usage_df,
    metadata=None,
    plot_type="heatmap",
    group=None,
    split=None,
    palette=None,
    cmap=PHEATMAP_CMAP,
    height=3.2,
    aspect=1.2,
    seed=0,
    combine_families=False,
):
    """Plot V, J, or C segment usage from a long or wide statistics table.

    Parameters
    ----------
    segment_usage_df : pandas.DataFrame
        Long output from ``stats.calc_segment_usage`` with a ``v``, ``j``, or
        ``c`` column and a value column, or its wide output with genes in
        columns. Generic long tables with ``segment`` or ``gene`` and one of
        ``usage``, ``value``, ``freq``, ``frequency``, or ``count`` are also
        accepted. The segment type and chain are detected from gene names.
    metadata : pandas.DataFrame, optional
        Metadata merged by ``sample_id`` and, when present in both tables,
        ``chain``. Group and split columns must come from this table.
    plot_type : {"heatmap", "barplot", "boxplot"}
        Heatmaps show segments in columns and samples in rows. Barplots show
        samples as separate series when ungrouped; with a group, bars show the
        mean and sample standard deviation. Boxplots require a group and add
        deterministic, horizontally jittered sample points. An ungrouped
        boxplot falls back to a barplot.
    group : str or sequence of str, optional
        One metadata column for barplots and boxplots, or up to three columns
        for heatmap sample annotations. Grouped heatmaps hierarchically cluster
        samples by their segment-usage profiles.
    split : str or one-item sequence of str, optional
        Metadata column placed in plot rows. Chains always form independent
        rows, with their own segment and sample axes.
    palette : seaborn palette or dict, optional
        Colors for samples, groups, and heatmap annotation categories.
    cmap : matplotlib colormap, optional
        Colormap for usage values in heatmaps. Defaults to the R ``pheatmap``
        palette.
    height, aspect : float
        Base panel height and width multiplier.
    seed : int, default 0
        Random seed used for horizontal jitter in boxplots.
    combine_families : bool, default False
        Sum segment usages within each sample and segment family before
        plotting. V and J genes are grouped by numeric family. C genes are
        grouped by isotype, so variants such as ``IGHG1`` and ``IGHG2`` are
        plotted together as ``IGHG``.

    Returns
    -------
    matplotlib.figure.Figure or None
        The figure, or ``None`` when a categorical plot exceeds ten series.
    """
    plot_type = str(plot_type).casefold()
    if plot_type not in {"heatmap", "barplot", "boxplot"}:
        raise ValueError("plot_type must be one of: heatmap, barplot, boxplot")

    data, _, group_columns, split_columns = _prepare_segment_usage_data(
        segment_usage_df,
        metadata,
        group,
        split,
        plot_type,
        combine_families,
    )
    segment_label = "Segment Family" if combine_families else "Segment"
    if plot_type == "boxplot" and not group_columns:
        warnings.warn(
            "boxplot requires a group column; using barplot instead.",
            UserWarning,
            stacklevel=2,
        )
        plot_type = "barplot"

    panels = _panel_subsets(data, split_columns)
    if plot_type == "heatmap":
        fig = _heatmap_segment_usage_plot(
            panels,
            group_columns,
            palette,
            cmap,
            height,
            aspect,
            segment_label,
        )
    else:
        fig = _categorical_segment_usage_plot(
            data,
            panels,
            plot_type,
            group_columns,
            palette,
            height,
            aspect,
            seed,
            segment_label,
        )

    # Inline backends display open pyplot figures and the returned object. A
    # closed Figure remains renderable and editable but appears only once.
    if fig is not None:
        plt.close(fig)
    return fig


def _normalize_cdr3_length_table(cdr3_length_df):
    if not isinstance(cdr3_length_df, pd.DataFrame):
        raise TypeError("cdr3_length_df must be a pandas DataFrame")
    if "sample_id" not in cdr3_length_df.columns:
        raise ValueError("cdr3_length_df must contain a 'sample_id' column")

    data = cdr3_length_df.copy()
    length_column = _column_by_name(data.columns, {"cdr3_length"})
    value_columns = [
        column
        for column in data.columns
        if isinstance(column, str) and column.casefold() in {"freq", "count"}
    ]
    id_columns = ["sample_id"]
    if "chain" in data.columns:
        id_columns.append("chain")

    if length_column is not None:
        if len(value_columns) != 1:
            raise ValueError(
                "Long CDR3-length tables must contain exactly one 'freq' or "
                "'count' column"
            )
        value_column = value_columns[0]
        data = data[id_columns + [length_column, value_column]].rename(
            columns={length_column: "_length", value_column: "_value"}
        )
        value_label = "Frequency" if value_column.casefold() == "freq" else "Count"
    else:
        wide_columns = [column for column in data.columns if column not in id_columns]
        parsed_lengths = pd.to_numeric(pd.Index(wide_columns), errors="coerce")
        valid_lengths = np.isfinite(parsed_lengths) & (parsed_lengths % 1 == 0)
        if not wide_columns or not valid_lengths.all():
            raise ValueError(
                "Wide CDR3-length tables must contain only integer length columns "
                "after 'sample_id' and optional 'chain'"
            )
        length_mapping = dict(zip(wide_columns, parsed_lengths.astype(int)))
        data = data.melt(
            id_vars=id_columns,
            value_vars=wide_columns,
            var_name="_length",
            value_name="_value",
        )
        data["_length"] = data["_length"].map(length_mapping)
        numeric_values = pd.to_numeric(data["_value"], errors="coerce")
        sample_totals = numeric_values.groupby(
            [data[column] for column in id_columns], observed=True
        ).sum()
        value_label = (
            "Frequency"
            if sample_totals.notna().all() and sample_totals.le(1 + 1e-8).all()
            else "Count"
        )

    data["_length"] = pd.to_numeric(data["_length"], errors="coerce")
    data["_value"] = pd.to_numeric(data["_value"], errors="coerce")
    valid_rows = (
        data["_length"].notna()
        & np.isfinite(data["_length"])
        & (data["_length"] % 1 == 0)
        & data["_value"].notna()
        & np.isfinite(data["_value"])
    )
    if not valid_rows.all():
        warnings.warn(
            f"Dropped {int((~valid_rows).sum())} CDR3-length row(s) with "
            "missing or non-numeric values.",
            UserWarning,
            stacklevel=2,
        )
        data = data.loc[valid_rows].copy()
    if data.empty:
        raise ValueError("No valid CDR3-length values were found")

    data["_length"] = data["_length"].astype(int)
    return data, value_label


def _cdr3_length_panels(data, split_columns):
    panel_columns = (["chain"] if "chain" in data.columns else []) + split_columns
    if not panel_columns:
        return [("CDR3 Length Distribution", data.copy())]

    levels = [_category_order(data[column]) for column in panel_columns]
    panels = []
    for combination in _cartesian_product(levels):
        panel_data = data
        for column, value in zip(panel_columns, combination):
            panel_data = panel_data.loc[panel_data[column] == value]
        if panel_data.empty:
            continue
        panels.append((" | ".join(map(str, combination)), panel_data.copy()))
    return panels


def _draw_cdr3_length_bars(ax, data, group_column, palette):
    length_order = sorted(pd.unique(data["_length"]))
    x = np.arange(len(length_order), dtype=float)
    series_column = group_column or "_sample_label"
    series_order = _category_order(data[series_column])
    color_map = _palette_mapping(series_order, palette)
    width = 0.8 / max(1, len(series_order))

    for index, series in enumerate(series_order):
        series_data = data.loc[data[series_column] == series]
        matrix = series_data.pivot_table(
            index="_sample_label",
            columns="_length",
            values="_value",
            aggfunc="sum",
            fill_value=0,
            sort=False,
        ).reindex(columns=length_order, fill_value=0)
        values = matrix.mean(axis=0) if group_column else matrix.sum(axis=0)
        offset = (index - (len(series_order) - 1) / 2) * width
        ax.bar(
            x + offset,
            values.to_numpy(),
            width=width,
            color=color_map[series],
            label=str(series),
        )

    ax.set_xticks(x)
    ax.set_xticklabels(length_order)


def cdr3_length_distributions(
    cdr3_length_df,
    metadata=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot CDR3-length distributions from a long or wide statistics table.

    Parameters
    ----------
    cdr3_length_df : pandas.DataFrame
        Long or wide output from :func:`stats.cdr3_length_distributions`.
    metadata : pandas.DataFrame, optional
        Sample metadata with one row per ``sample_id`` or ``sample_id`` plus
        ``chain``. Grouping and splitting columns must come from this table.
    group : str, optional
        One metadata column used to group samples. Bars show the sample mean
        for each group without error bars. Without a group, each sample is a
        separate series.
    split : str or sequence of str, optional
        Up to two metadata columns used to split plots into panels.
    palette : optional
        Any seaborn/matplotlib-compatible palette specification.
    height, aspect : float
        Panel height and width scaling.

    Returns
    -------
    matplotlib.figure.Figure or None
        The created figure, or ``None`` when more than ten series are present.
    """
    data, value_label = _normalize_cdr3_length_table(cdr3_length_df)
    group_columns = _as_list(group, "group", max_len=1)
    split_columns = _as_list(split, "split", max_len=2)
    data, metadata_columns, _ = _merge_stats_metadata(data, metadata)
    if metadata is None and (group_columns or split_columns):
        raise ValueError("metadata is required when group or split columns are used")
    _validate_metadata_columns(group_columns, metadata_columns, "group")
    _validate_metadata_columns(split_columns, metadata_columns, "split")
    data["_sample_label"] = data["sample_id"].astype(str)

    group_column = group_columns[0] if group_columns else None
    series_column = group_column or "_sample_label"
    series_count = data[series_column].nunique(dropna=True)
    if series_count > 10:
        warnings.warn(
            "CDR3-length plots support at most 10 groups; "
            f"found {series_count}. Nothing was plotted.",
            UserWarning,
            stacklevel=2,
        )
        return None

    panels = _cdr3_length_panels(data, split_columns)
    max_lengths = max(panel["_length"].nunique() for _, panel in panels)
    figure_width = max(8, min(24, max_lengths * 0.38 * aspect))
    fig, axes_array = plt.subplots(
        len(panels),
        1,
        figsize=(figure_width, height * len(panels)),
        squeeze=False,
    )
    axes = list(axes_array[:, 0])
    for (title, panel_data), ax in zip(panels, axes):
        _draw_cdr3_length_bars(ax, panel_data, group_column, palette)
        ax.set_title(title)
        ax.set_xlabel("CDR3 Length")
        ax.set_ylabel(value_label)
        sns.despine(ax=ax)

    _add_figure_legend(
        fig,
        axes,
        _caption(group_column) if group_column else "Sample",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    plt.close(fig)
    return fig


cdr3_length_distribution = cdr3_length_distributions


def _parse_pipe_usage_combination(value, expected_length):
    if not isinstance(value, str):
        return None
    combination = value.split("|")
    if len(combination) != expected_length:
        return None
    v_gene, j_gene = combination[:2]
    if parse_gene_name(str(v_gene)) is None or parse_gene_name(str(j_gene)) is None:
        return None
    if expected_length == 3:
        try:
            cdr3_length = int(combination[2])
        except (TypeError, ValueError):
            return None
        return str(v_gene), str(j_gene), cdr3_length
    return str(v_gene), str(j_gene)


def _normalize_combination_usage_table(usage_df, combination_type):
    if not isinstance(usage_df, pd.DataFrame):
        raise TypeError("usage_df must be a pandas DataFrame")
    if "sample_id" not in usage_df.columns:
        raise ValueError("usage_df must contain a 'sample_id' column")

    expected_length = 2 if combination_type == "vj" else 3
    data = usage_df.copy()
    combination_column = _column_by_name(data.columns, {combination_type})
    value_column = _column_by_name(
        data.columns, {"usage", "value", "freq", "frequency", "count"}
    )
    v_column = _column_by_name(data.columns, {"v"})
    j_column = _column_by_name(data.columns, {"j"})
    length_column = _column_by_name(data.columns, {"len"})

    id_columns = ["sample_id"] + (["chain"] if "chain" in data.columns else [])
    required_components = [v_column, j_column]
    if expected_length == 3:
        required_components.append(length_column)

    if combination_column is not None:
        if value_column is None:
            raise ValueError(
                f"Long {combination_type} tables require a value column"
            )
        component_columns = [
            column for column in required_components if column is not None
        ]
        keep_columns = id_columns + [combination_column] + component_columns
        keep_columns.append(value_column)
        data = data[keep_columns].rename(columns={value_column: "_value"})
        if not data[combination_column].map(lambda value: isinstance(value, str)).all():
            raise ValueError(
                f"{combination_type} values must use pipe-delimited strings, "
                "not tuples"
            )
        combinations = data[combination_column].map(
            lambda value: _parse_pipe_usage_combination(value, expected_length)
        )
        valid_combinations = combinations.notna()
        component_mismatch = pd.Series(False, index=data.index)
        if v_column is not None:
            parsed_v = combinations.map(
                lambda combination: combination[0] if combination is not None else None
            )
            component_mismatch |= valid_combinations & (
                data[v_column].astype(str) != parsed_v
            )
        if j_column is not None:
            parsed_j = combinations.map(
                lambda combination: combination[1] if combination is not None else None
            )
            component_mismatch |= valid_combinations & (
                data[j_column].astype(str) != parsed_j
            )
        if expected_length == 3 and length_column is not None:
            parsed_length = combinations.map(
                lambda combination: combination[2] if combination is not None else np.nan
            )
            provided_length = pd.to_numeric(data[length_column], errors="coerce")
            component_mismatch |= valid_combinations & (
                provided_length != parsed_length
            )
        if component_mismatch.any():
            raise ValueError(
                f"{combination_type} values must match the separate v, j"
                + (", and len" if expected_length == 3 else "")
                + " columns when provided"
            )
    elif value_column is not None and all(
        column is not None for column in required_components
    ):
        keep_columns = id_columns + required_components + [value_column]
        data = data[keep_columns].rename(columns={value_column: "_value"})
        identifiers = data[v_column].astype(str) + "|" + data[j_column].astype(str)
        if expected_length == 3:
            identifiers += "|" + data[length_column].astype(str)
        combinations = identifiers.map(
            lambda value: _parse_pipe_usage_combination(value, expected_length)
        )
    elif any(
        column is not None
        for column in [value_column, v_column, j_column, length_column]
    ):
        raise ValueError(
            f"Long {combination_type} tables require either a pipe-delimited "
            f"{combination_type} column or separate v, j"
            + (", and len" if expected_length == 3 else "")
            + " columns, plus a value column"
        )
    else:
        tuple_columns = [column for column in data.columns if isinstance(column, tuple)]
        if tuple_columns:
            raise ValueError(
                f"{combination_type} wide tables must use pipe-delimited string "
                "column names, not tuples"
            )
        wide_columns = [
            column
            for column in data.columns
            if column not in id_columns
            and _parse_pipe_usage_combination(column, expected_length) is not None
        ]
        if not wide_columns:
            raise ValueError(
                f"Could not detect a long or wide {combination_type} usage table"
            )
        data = data.melt(
            id_vars=id_columns,
            value_vars=wide_columns,
            var_name="_combination",
            value_name="_value",
        )
        combinations = data["_combination"].map(
            lambda value: _parse_pipe_usage_combination(value, expected_length)
        )

    data["_value"] = pd.to_numeric(data["_value"], errors="coerce")
    invalid = combinations.isna() | data["_value"].isna()
    if invalid.any():
        warnings.warn(
            f"Dropped {int(invalid.sum())} invalid {combination_type} usage row(s).",
            UserWarning,
            stacklevel=2,
        )
        data = data.loc[~invalid].copy()
        combinations = combinations.loc[~invalid]
    if data.empty:
        raise ValueError(f"No valid {combination_type} usage values were found")

    data["_v"] = [combination[0] for combination in combinations]
    data["_j"] = [combination[1] for combination in combinations]
    if expected_length == 3:
        data["_length"] = [combination[2] for combination in combinations]

    inferred_chains = data["_v"].map(_segment_chain)
    j_chains = data["_j"].map(_segment_chain)
    mismatched_gene_chains = (
        inferred_chains.notna() & j_chains.notna() & (inferred_chains != j_chains)
    )
    if mismatched_gene_chains.any():
        raise ValueError("V and J genes must belong to the same chain")

    if "chain" not in data.columns:
        data["chain"] = inferred_chains.fillna("Unknown")
    else:
        data["chain"] = data["chain"].where(data["chain"].notna(), inferred_chains)
        data["chain"] = data["chain"].fillna("Unknown")
        cross_chain = inferred_chains.notna() & (
            data["chain"].astype(str) != inferred_chains.astype(str)
        )
        data = data.loc[~cross_chain].copy()

    output_columns = ["sample_id", "chain", "_v", "_j"]
    if expected_length == 3:
        output_columns.append("_length")
    output_columns.append("_value")
    return data[output_columns]


def _combine_combination_family_usage(data, combination_type):
    data = data.copy()
    data["_v"] = data["_v"].map(lambda gene: _segment_family(gene, "v"))
    data["_j"] = data["_j"].map(lambda gene: _segment_family(gene, "j"))
    group_columns = ["sample_id", "chain", "_v", "_j"]
    if combination_type == "vjlen":
        group_columns.append("_length")
    return (
        data.groupby(group_columns, as_index=False, sort=False)["_value"]
        .sum()
    )


def _prepare_combination_plot_data(
    usage_df,
    metadata,
    group,
    split,
    combination_type,
    combine_families=False,
):
    data = _normalize_combination_usage_table(usage_df, combination_type)
    if combine_families:
        data = _combine_combination_family_usage(data, combination_type)
    group_columns = _as_list(group, "group", max_len=1)
    split_columns = _as_list(split, "split", max_len=2)
    data, metadata_columns, _ = _merge_stats_metadata(data, metadata)
    if metadata is None and (group_columns or split_columns):
        raise ValueError("metadata is required when group or split columns are used")
    _validate_metadata_columns(group_columns, metadata_columns, "group")
    _validate_metadata_columns(split_columns, metadata_columns, "split")
    data["_sample_label"] = data["sample_id"].astype(str)
    return data, group_columns, split_columns


def _ordered_observed_values(data, column):
    return _category_order(data[column])


def _combination_facet_layout(data, split_columns, include_chain):
    row_columns = (["chain"] if include_chain else []) + split_columns[:1]
    column_columns = split_columns[1:2]
    row_levels = (
        list(
            _cartesian_product(
                [_ordered_observed_values(data, column) for column in row_columns]
            )
        )
        if row_columns
        else [()]
    )
    column_levels = (
        list(
            _cartesian_product(
                [_ordered_observed_values(data, column) for column in column_columns]
            )
        )
        if column_columns
        else [()]
    )
    return row_columns, column_columns, row_levels, column_levels


def _facet_subset(data, columns, values):
    subset = data
    for column, value in zip(columns, values):
        subset = subset.loc[subset[column] == value]
    return subset.copy()


def _facet_title(row_columns, row_values, column_columns, column_values):
    values = list(row_values) + list(column_values)
    return " | ".join(map(str, values)) if values else ""


def _close_and_return(fig):
    plt.close(fig)
    return fig


def _scaled_dot_sizes(values, size_range):
    minimum_size, maximum_size = map(float, size_range)
    if minimum_size <= 0 or maximum_size < minimum_size:
        raise ValueError("size_range must contain two positive increasing values")
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return values
    maximum_value = np.nanmax(values)
    if maximum_value <= 0:
        return np.full(len(values), minimum_size)
    return minimum_size + (values / maximum_value) * (maximum_size - minimum_size)


def _aggregate_vj_panel(panel_data, series_column, grouped):
    aggregation = "mean" if grouped else "sum"
    return (
        panel_data.groupby(["_v", "_j", series_column], observed=True)["_value"]
        .agg(aggregation)
        .reset_index()
    )


def vj_usage(
    usage_df,
    metadata=None,
    group=None,
    split=None,
    palette=None,
    size_range=(20, 800),
    repel=0.18,
    height=5.0,
    aspect=1.2,
    combine_families=False,
):
    """Plot V-J usage as categorical bubbles.

    Point area represents usage. With no group, colors and radial offsets
    identify samples; with one metadata group, each point is the group mean.
    At most eight samples or group levels can be displayed. One or two split
    columns create facet rows and columns, while chains always occupy separate
    rows.

    Parameters
    ----------
    usage_df : pandas.DataFrame
        Long or wide output of ``calc_segment_usage(segment="vj")``. Explicit
        long ``v`` and ``j`` columns plus a value column are also accepted.
    metadata : pandas.DataFrame, optional
        Metadata merged by sample and chain.
    group : str, optional
        Metadata column whose sample means are plotted using distinct colors.
    split : str or sequence of up to two str, optional
        First column creates rows and second creates columns.
    palette : seaborn palette or dict, optional
        Sample or group fill colors.
    size_range : pair of float, default (20, 800)
        Minimum and maximum marker areas.
    repel : float, default 0.18
        Radial displacement around each categorical V-J center.
    height, aspect : float
        Facet dimensions.
    combine_families : bool, default False
        Sum usage within each sample by V and J segment family before plotting.

    Returns
    -------
    matplotlib.figure.Figure
        A closed figure that renders once in Jupyter.
    """
    data, group_columns, split_columns = _prepare_combination_plot_data(
        usage_df,
        metadata,
        group,
        split,
        "vj",
        combine_families,
    )
    grouped = bool(group_columns)
    series_column = group_columns[0] if grouped else "_sample_label"
    series_order = _category_order(data[series_column])
    if len(series_order) > 8:
        label = "groups" if grouped else "samples"
        raise ValueError(f"vj_usage supports at most 8 {label}; found {len(series_order)}")
    color_map = _palette_mapping(series_order, palette)

    row_columns, column_columns, row_levels, column_levels = (
        _combination_facet_layout(data, split_columns, include_chain=True)
    )
    fig, axes = plt.subplots(
        len(row_levels),
        len(column_levels),
        figsize=(height * aspect * len(column_levels), height * len(row_levels)),
        squeeze=False,
    )
    series_rank = {series: index for index, series in enumerate(series_order)}

    for row_index, row_values in enumerate(row_levels):
        row_data = _facet_subset(data, row_columns, row_values)
        for column_index, column_values in enumerate(column_levels):
            ax = axes[row_index, column_index]
            panel_data = _facet_subset(row_data, column_columns, column_values)
            if panel_data.empty:
                ax.set_visible(False)
                continue

            v_order = _sort_gene_names(panel_data["_v"])
            j_order = _sort_gene_names(panel_data["_j"])
            v_positions = {gene: index for index, gene in enumerate(v_order)}
            j_positions = {gene: index for index, gene in enumerate(j_order)}
            plotted = _aggregate_vj_panel(panel_data, series_column, grouped)
            plotted = plotted.loc[plotted["_value"] > 0].copy()
            plotted = plotted.sort_values("_value", ascending=False)

            angles = np.zeros(len(plotted), dtype=float)
            displacement = np.zeros(len(plotted), dtype=float)
            for indices in plotted.groupby(["_v", "_j"], sort=False).indices.values():
                ordered_indices = sorted(
                    indices,
                    key=lambda index: series_rank[plotted.iloc[index][series_column]],
                )
                if len(ordered_indices) > 1:
                    for position, index in enumerate(ordered_indices):
                        angles[index] = 2 * np.pi * position / len(ordered_indices)
                        displacement[index] = float(repel)
            x = plotted["_v"].map(v_positions).to_numpy(dtype=float)
            y = plotted["_j"].map(j_positions).to_numpy(dtype=float)
            x += displacement * np.cos(angles)
            y += displacement * np.sin(angles)
            ax.scatter(
                x,
                y,
                s=_scaled_dot_sizes(plotted["_value"], size_range),
                c=[color_map[value] for value in plotted[series_column]],
                edgecolors="black",
                linewidths=0.6,
            )
            ax.set_xticks(np.arange(len(v_order)))
            ax.set_xticklabels(v_order, rotation=90)
            ax.set_yticks(np.arange(len(j_order)))
            ax.set_yticklabels(j_order)
            ax.set_xlim(-0.6, len(v_order) - 0.4)
            ax.set_ylim(-0.6, len(j_order) - 0.4)
            ax.set_xlabel(
                "V segment family" if combine_families else "V segment"
            )
            ax.set_ylabel(
                "J segment family" if combine_families else "J segment"
            )
            ax.grid(color="#e5e5e5", linewidth=0.6, zorder=0)
            ax.set_axisbelow(True)
            ax.set_title(
                _facet_title(
                    row_columns, row_values, column_columns, column_values
                )
            )

    legend_handles = [
        Patch(facecolor=color_map[level], edgecolor="black", label=str(level))
        for level in series_order
    ]
    fig.legend(
        handles=legend_handles,
        title=_caption(series_column) if grouped else "Sample",
        loc="upper center",
        ncol=min(8, len(legend_handles)),
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _close_and_return(fig)


def _short_gene_label(gene):
    gene = str(gene).strip().split(",", 1)[0].split(";", 1)[0].split("|", 1)[0]
    parsed = parse_gene_name(gene)
    if parsed is None:
        return gene
    prefix = parsed["system"] + parsed["chain"]
    return gene[len(prefix):] if gene.upper().startswith(prefix) else gene


def _vjlen_label(v_gene, j_gene, cdr3_length):
    return f"{_short_gene_label(v_gene)}|{_short_gene_label(j_gene)}|{int(cdr3_length)}"


def _aggregate_vjlen_panel(panel_data, series_column, grouped, series_order):
    aggregation = "mean" if grouped else "sum"
    aggregated = (
        panel_data.groupby(
            ["_v", "_j", "_length", series_column], observed=True
        )["_value"]
        .agg(aggregation)
        .reset_index()
    )
    return aggregated.pivot_table(
        index=["_v", "_j", "_length"],
        columns=series_column,
        values="_value",
        aggfunc="sum",
        fill_value=0,
        sort=False,
    ).reindex(columns=series_order, fill_value=0)


def _log_floor(x, y):
    positive = np.concatenate([x[x > 0], y[y > 0]])
    if len(positive) == 0:
        raise ValueError("Log-scale VJ-length plots require at least one positive value")
    return float(np.min(positive)) / 10


def _add_isolated_vjlen_labels(
    ax,
    x,
    y,
    labels,
    max_labels,
    minimum_distance=35,
):
    if not labels or max_labels <= 0:
        return
    ax.figure.canvas.draw()
    display_coordinates = ax.transData.transform(np.column_stack([x, y]))
    if len(display_coordinates) == 1:
        nearest_distances = np.array([np.inf])
    else:
        differences = display_coordinates[:, None, :] - display_coordinates[None, :, :]
        distances = np.sqrt(np.sum(differences ** 2, axis=2))
        np.fill_diagonal(distances, np.inf)
        nearest_distances = distances.min(axis=1)

    candidates = [
        index
        for index in np.argsort(-nearest_distances)
        if nearest_distances[index] >= minimum_distance
    ][:max_labels]
    renderer = ax.figure.canvas.get_renderer()
    occupied = []
    offsets = [(5, 5), (5, -9), (-5, 5), (-5, -9), (10, 0), (-10, 0)]
    for index in candidates:
        for x_offset, y_offset in offsets:
            annotation = ax.annotate(
                labels[index],
                xy=(x[index], y[index]),
                xytext=(x_offset, y_offset),
                textcoords="offset points",
                ha="left" if x_offset >= 0 else "right",
                va="bottom" if y_offset >= 0 else "top",
                fontsize=7,
                arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.5},
            )
            ax.figure.canvas.draw()
            bounds = annotation.get_window_extent(renderer=renderer).expanded(1.05, 1.1)
            if not any(bounds.overlaps(existing) for existing in occupied):
                occupied.append(bounds)
                break
            annotation.remove()


def vjlen_usage(
    usage_df,
    metadata=None,
    group=None,
    split=None,
    log_scale=False,
    labels=False,
    max_labels=20,
    color="#d62728",
    alpha=0.6,
    dot_size=45,
    height=4.5,
    aspect=1.0,
    combine_families=False,
):
    """Compare V-J-CDR3-length usage between exactly two series per facet.

    The x and y axes represent usage in the two samples, or mean usage in two
    metadata groups. The input must contain one chain. Up to two metadata split
    columns form facet rows and columns. Optional labels are added only for
    isolated points and use compact forms such as ``V12-1|J1-2|15``.

    Parameters
    ----------
    usage_df : pandas.DataFrame
        Long or wide output of ``calc_segment_usage(segment="vjlen")``.
        Explicit ``v``, ``j``, length, and value columns are also accepted.
    metadata : pandas.DataFrame, optional
        Metadata merged by sample and chain.
    group : str, optional
        Metadata column containing exactly two observed levels in each facet.
        Values are averaged across samples within each level.
    split : str or sequence of up to two str, optional
        First metadata column creates rows; second creates columns.
    log_scale : bool, default False
        Use logarithmic frequency axes. Zeros are placed one decade below the
        smallest positive value in their panel.
    labels : bool, default False
        Label isolated V-J-length points with compact gene names.
    max_labels : int, default 20
        Maximum labels in each panel.
    color : matplotlib color, default "#d62728"
        Marker fill color.
    alpha : float, default 0.6
        Marker opacity.
    dot_size : float, default 45
        Marker area.
    height, aspect : float
        Facet dimensions.
    combine_families : bool, default False
        Sum usage within each sample by V family, J family, and CDR3 length
        before plotting.

    Returns
    -------
    matplotlib.figure.Figure
        A closed figure that renders once in Jupyter.
    """
    data, group_columns, split_columns = _prepare_combination_plot_data(
        usage_df,
        metadata,
        group,
        split,
        "vjlen",
        combine_families,
    )
    chains = _category_order(data["chain"])
    if len(chains) != 1:
        raise ValueError(
            "vjlen_usage requires exactly one chain in the input table; "
            f"found {len(chains)}: {', '.join(map(str, chains))}"
        )

    grouped = bool(group_columns)
    series_column = group_columns[0] if grouped else "_sample_label"
    row_columns, column_columns, row_levels, column_levels = (
        _combination_facet_layout(data, split_columns, include_chain=False)
    )
    fig, axes = plt.subplots(
        len(row_levels),
        len(column_levels),
        figsize=(height * aspect * len(column_levels), height * len(row_levels)),
        squeeze=False,
    )

    for row_index, row_values in enumerate(row_levels):
        row_data = _facet_subset(data, row_columns, row_values)
        for column_index, column_values in enumerate(column_levels):
            ax = axes[row_index, column_index]
            panel_data = _facet_subset(row_data, column_columns, column_values)
            if panel_data.empty:
                ax.set_visible(False)
                continue

            series_order = _category_order(panel_data[series_column])
            if len(series_order) != 2:
                facet_name = _facet_title(
                    row_columns, row_values, column_columns, column_values
                ) or "unsplit panel"
                label = "groups" if grouped else "samples"
                plt.close(fig)
                raise ValueError(
                    "vjlen_usage requires exactly 2 "
                    f"{label} in each panel; {facet_name} has {len(series_order)}"
                )

            comparison = _aggregate_vjlen_panel(
                panel_data, series_column, grouped, series_order
            )
            comparison = comparison.loc[(comparison > 0).any(axis=1)].copy()
            x = comparison[series_order[0]].to_numpy(dtype=float)
            y = comparison[series_order[1]].to_numpy(dtype=float)
            if log_scale:
                floor = _log_floor(x, y)
                x = np.maximum(x, floor)
                y = np.maximum(y, floor)

            ax.scatter(
                x,
                y,
                s=float(dot_size),
                facecolors=color,
                edgecolors="black",
                linewidths=0.6,
                alpha=float(alpha),
            )
            if log_scale:
                ax.set_xscale("log")
                ax.set_yscale("log")
            ax.set_xlabel(str(series_order[0]))
            ax.set_ylabel(str(series_order[1]))
            ax.set_title(
                _facet_title(
                    row_columns, row_values, column_columns, column_values
                )
                or str(chains[0])
            )
            ax.grid(color="#e5e5e5", linewidth=0.6)
            ax.set_axisbelow(True)

            if labels:
                point_labels = [
                    _vjlen_label(v_gene, j_gene, cdr3_length)
                    for v_gene, j_gene, cdr3_length in comparison.index
                ]
                _add_isolated_vjlen_labels(
                    ax,
                    x,
                    y,
                    point_labels,
                    int(max_labels),
                )

    fig.tight_layout()
    return _close_and_return(fig)


def _normalize_beta_metric_key(value):
    return str(value).strip().casefold().replace("-", "_").replace(" ", "_")


def _select_beta_metric_matrix(beta_data, metric):
    if isinstance(beta_data, dict):
        keys = list(beta_data)
        if metric is None:
            message = (
                "No beta-diversity metric was selected. Available keys: "
                + ", ".join(map(str, keys))
                + ". Pass one with metric='<key>'."
            )
            warnings.warn(message, UserWarning, stacklevel=2)
            print(message)
            return None, None
        normalized_keys = {
            _normalize_beta_metric_key(key): key
            for key in keys
        }
        selected_key = normalized_keys.get(_normalize_beta_metric_key(metric))
        if selected_key is None:
            raise ValueError(
                f"Unknown beta-diversity metric {metric!r}. Available keys: "
                + ", ".join(map(str, keys))
            )
        matrix = beta_data[selected_key]
        title = str(selected_key)
    elif isinstance(beta_data, pd.DataFrame):
        matrix = beta_data
        title = str(beta_data.attrs.get("metric", "Beta metric"))
    else:
        raise TypeError("beta_data must be a metric dictionary or pandas DataFrame")

    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("The selected beta-diversity metric must be a pandas DataFrame")
    if matrix.empty:
        raise ValueError("The beta-diversity metric matrix is empty")
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    invalid = matrix.notna() & numeric.isna()
    if invalid.any().any():
        raise ValueError("The beta-diversity metric matrix must contain numeric values")
    numeric.index = matrix.index
    numeric.columns = matrix.columns
    return numeric.astype(float), title


def _beta_log_values(values, base=10):
    values = values.copy()
    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if np.any(finite < 0):
        raise ValueError("Log-transformed beta values cannot contain negatives")
    positive = finite[finite > 0]
    floor = (
        float(np.min(positive)) / float(base)
        if len(positive)
        else np.finfo(float).tiny
    )
    values = values.mask(values == 0, floor)
    return np.log(values) / np.log(float(base))


def _format_beta_value(value):
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value == 0:
        return "0"
    if value == 1:
        return "1"
    if value.is_integer() and abs(value) < 10000:
        candidate = f"{value:.1f}"
        if len(candidate) <= 6:
            return candidate
    for precision in range(4, 0, -1):
        candidate = f"{value:.{precision}g}"
        candidate = re.sub(r"e([+-])0+(\d+)", r"e\1\2", candidate)
        candidate = candidate.replace("e+", "e")
        if len(candidate) <= 6:
            return candidate
    return candidate[:6]


def _format_de_heatmap_value(value, all_values_are_integers=False):
    if pd.isna(value):
        return "NA"
    numeric_value = float(value)
    if numeric_value.is_integer() and (
        all_values_are_integers or numeric_value > 1
    ):
        return str(int(numeric_value))
    return _format_beta_value(numeric_value)


def _beta_linkage(matrix, axis):
    values = matrix.to_numpy(dtype=float)
    if axis == 1:
        values = values.T
    if len(values) < 2:
        return None
    feature_means = np.nanmean(values, axis=0)
    overall_mean = np.nanmean(values)
    if not np.isfinite(overall_mean):
        overall_mean = 0.0
    feature_means = np.where(np.isfinite(feature_means), feature_means, overall_mean)
    missing_rows, missing_columns = np.where(~np.isfinite(values))
    values = values.copy()
    values[missing_rows, missing_columns] = feature_means[missing_columns]
    return linkage(values, method="average", metric="euclidean")


def _prepare_beta_annotation_colors(metadata, group_columns, samples):
    if not group_columns:
        return None, []
    if metadata is None:
        raise ValueError("metadata is required when group columns are used")
    if "sample_id" not in metadata.columns:
        raise ValueError("metadata must contain a 'sample_id' column")
    missing = [column for column in group_columns if column not in metadata.columns]
    if missing:
        raise ValueError("group columns not found in metadata: " + ", ".join(missing))

    sample_metadata = metadata[["sample_id"] + group_columns].drop_duplicates()
    if sample_metadata["sample_id"].duplicated().any():
        raise ValueError("metadata must have one set of group values per sample_id")
    sample_metadata = sample_metadata.set_index("sample_id")
    colors = pd.DataFrame(index=pd.Index(samples), columns=group_columns, dtype=object)
    handles = []
    for group_column in group_columns:
        levels = _category_order(metadata[group_column])
        color_map = _palette_mapping(levels)
        values = sample_metadata[group_column].reindex(samples)
        colors[group_column] = [color_map.get(value, "#d0d0d0") for value in values]
        handles.extend(
            Patch(
                facecolor=color_map[level],
                label=f"{_caption(group_column)}: {level}",
            )
            for level in levels
        )
    return colors, handles


def _groups_are_contiguous(groups):
    seen = set()
    previous = object()
    for group in groups:
        if group != previous:
            if group in seen:
                return False
            seen.add(group)
            previous = group
    return True


def _de_sample_order(sample_columns, sample_groups, group_order):
    groups = [sample_groups[sample] for sample in sample_columns]
    if _groups_are_contiguous(groups):
        return list(sample_columns)
    group_positions = {group: index for index, group in enumerate(group_order)}
    return sorted(
        sample_columns,
        key=lambda sample: group_positions[sample_groups[sample]],
    )


def _de_group_boundaries(sample_order, sample_groups):
    return [
        index
        for index in range(1, len(sample_order))
        if sample_groups[sample_order[index - 1]] != sample_groups[sample_order[index]]
    ]


def de_heatmap(
    statistics_table,
    samples_metadata,
    feature_column=None,
    log_values=True,
    show_values=True,
    cmap=PHEATMAP_CMAP,
    group_palette=None,
    height=8,
    aspect=1.0,
):
    """Plot differential-enrichment counts with matched group annotations.

    Parameters
    ----------
    statistics_table : pandas.DataFrame
        A user-filtered table returned by ``diff_enrichment.calc_statistics``.
        It must contain ``enriched_in`` and numeric sample count columns. If
        ``postfilter_pass`` or ``prefilter_pass`` is present, only ``True``
        rows are plotted.
    samples_metadata : pandas.DataFrame
        Sample metadata containing unique ``sample_id`` and ``group`` columns.
    feature_column : hashable, optional
        Column used for heatmap row labels. Defaults to the first column.
    log_values : bool, default True
        Color cells by log10 counts after replacing zeros with one tenth of the
        smallest positive value.
    show_values : bool, default True
        Display compact original count values in heatmap cells.
    cmap : matplotlib colormap, optional
        Heatmap colormap. Defaults to the R ``pheatmap`` palette.
    group_palette : palette name, sequence, or dict, optional
        Colors shared by sample-group and ``enriched_in`` annotations.
    height, aspect : float
        Figure height and width multiplier.

    Returns
    -------
    matplotlib.figure.Figure
        Closed annotated heatmap figure.
    """
    if not isinstance(statistics_table, pd.DataFrame):
        raise TypeError("statistics_table must be a pandas DataFrame")
    if statistics_table.empty:
        raise ValueError("statistics_table is empty")
    if not isinstance(samples_metadata, pd.DataFrame):
        raise TypeError("samples_metadata must be a pandas DataFrame")
    required_metadata_columns = {"sample_id", "group"}
    missing_metadata_columns = required_metadata_columns.difference(
        samples_metadata.columns
    )
    if missing_metadata_columns:
        raise ValueError(
            "samples_metadata must contain columns 'sample_id' and 'group'; "
            f"missing: {sorted(missing_metadata_columns)}"
        )
    if "enriched_in" not in statistics_table.columns:
        raise ValueError("statistics_table must contain an 'enriched_in' column")
    plotted_table = statistics_table
    for pass_column in ["postfilter_pass", "prefilter_pass"]:
        if pass_column not in plotted_table.columns:
            continue
        valid_pass_values = plotted_table[pass_column].isin([True, False])
        if not valid_pass_values.all():
            raise ValueError(
                f"statistics_table[{pass_column!r}] must contain only True or False"
            )
        plotted_table = plotted_table.loc[plotted_table[pass_column].eq(True)]
    if plotted_table.empty:
        raise ValueError(
            "No features remain after applying prefilter_pass and "
            "postfilter_pass"
        )
    if plotted_table["enriched_in"].isna().any():
        raise ValueError(
            "statistics_table['enriched_in'] contains missing values; filter "
            "out rows without differential-enrichment statistics before plotting"
        )
    if samples_metadata["sample_id"].isna().any():
        raise ValueError("samples_metadata['sample_id'] must not contain missing values")
    if samples_metadata["group"].isna().any():
        raise ValueError("samples_metadata['group'] must not contain missing values")
    duplicated_samples = samples_metadata["sample_id"].duplicated(keep=False)
    if duplicated_samples.any():
        duplicates = samples_metadata.loc[duplicated_samples, "sample_id"].tolist()
        raise ValueError(
            "samples_metadata['sample_id'] values must be unique; duplicated values: "
            f"{duplicates}"
        )

    feature_column = (
        plotted_table.columns[0]
        if feature_column is None
        else feature_column
    )
    if feature_column not in plotted_table.columns:
        raise ValueError(f"feature_column {feature_column!r} is not in statistics_table")

    metadata = samples_metadata.set_index("sample_id")
    metadata_groups = set(metadata["group"])
    enriched_groups = list(pd.unique(plotted_table["enriched_in"]))
    missing_enriched_groups = [
        group for group in enriched_groups if group not in metadata_groups
    ]
    if missing_enriched_groups:
        raise ValueError(
            "The following enriched_in groups are absent from "
            f"samples_metadata['group']: {missing_enriched_groups}"
        )

    numeric_columns = list(
        plotted_table.select_dtypes(include="number").columns
    )
    metadata_sample_ids = set(samples_metadata["sample_id"])
    sample_columns = [
        column
        for column in plotted_table.columns
        if column in metadata_sample_ids and column in numeric_columns
    ]
    if not sample_columns:
        raise ValueError(
            "No numeric count-table columns match samples_metadata['sample_id']"
        )
    non_numeric_samples = [
        column
        for column in plotted_table.columns
        if column in metadata_sample_ids and column not in numeric_columns
    ]
    if non_numeric_samples:
        raise ValueError(
            "Sample count columns must be numeric; non-numeric samples: "
            f"{non_numeric_samples}"
        )

    sample_groups = metadata.loc[sample_columns, "group"].to_dict()
    represented_metadata = metadata.loc[sample_columns].reset_index()
    group_order = _category_order(represented_metadata["group"])
    sample_order = _de_sample_order(
        sample_columns,
        sample_groups,
        group_order,
    )
    matrix = plotted_table[sample_order].apply(pd.to_numeric, errors="coerce")
    invalid = plotted_table[sample_order].notna() & matrix.isna()
    if invalid.any().any():
        raise ValueError("Sample count columns must contain numeric values")
    matrix.index = pd.Index(
        plotted_table[feature_column].astype(str),
        name=str(feature_column),
    )
    original_values = matrix.astype(float)
    color_values = original_values.copy()
    if log_values:
        color_values = _beta_log_values(color_values, base=10)

    color_levels = _category_order(samples_metadata["group"])
    color_map = _palette_mapping(color_levels, palette=group_palette)
    row_colors = pd.DataFrame(
        {
            "Enriched in": [
                color_map[group]
                for group in plotted_table["enriched_in"]
            ]
        },
        index=matrix.index,
    )
    column_colors = pd.DataFrame(
        {
            "Group": [color_map[sample_groups[sample]] for sample in sample_order]
        },
        index=pd.Index(sample_order),
    )

    grid = sns.clustermap(
        color_values,
        cmap=cmap,
        mask=color_values.isna(),
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        col_colors=column_colors,
        figsize=(height * aspect, height),
        cbar_kws={"label": "log10(Count)" if log_values else "Count"},
        xticklabels=True,
        yticklabels=True,
    )
    grid.ax_heatmap.set_title("Differential enrichment")
    grid.ax_heatmap.set_xlabel("Sample")
    grid.ax_heatmap.set_ylabel(str(feature_column))
    grid.ax_heatmap.tick_params(axis="x", labelrotation=90)

    for boundary in _de_group_boundaries(sample_order, sample_groups):
        grid.ax_heatmap.axvline(
            boundary,
            color="white",
            linewidth=4,
            zorder=10,
            clip_on=False,
        )
        if grid.ax_col_colors is not None:
            grid.ax_col_colors.axvline(
                boundary,
                color="white",
                linewidth=4,
                zorder=10,
                clip_on=False,
            )

    if show_values:
        finite_values = original_values.to_numpy(dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        all_values_are_integers = bool(len(finite_values)) and np.all(
            finite_values == np.floor(finite_values)
        )
        for row_index, row in enumerate(original_values.to_numpy(dtype=float)):
            for column_index, value in enumerate(row):
                grid.ax_heatmap.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    _format_de_heatmap_value(
                        value,
                        all_values_are_integers=all_values_are_integers,
                    ),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                )

    handles = [
        Patch(facecolor=color_map[group], label=str(group))
        for group in color_levels
    ]
    if handles:
        grid.fig.legend(
            handles=handles,
            title="Group",
            loc="upper center",
            ncol=min(4, len(handles)),
            frameon=False,
        )
        grid.fig.subplots_adjust(top=0.9)
    return _close_and_return(grid.fig)


def beta_mds(
    mds_table,
    metadata=None,
    group=None,
    split=None,
    centroids=False,
    dispersion=False,
    palette=None,
    height=4,
    aspect=1.2,
):
    """Plot the two-dimensional table returned by ``beta.mds``.

    One metadata column may color points and up to two metadata columns may
    split them into facets. Optional centroids connect to their group members,
    and dispersion draws one-standard-deviation covariance ellipses.
    """
    required = {"sample_id", "MDS1", "MDS2"}
    if not isinstance(mds_table, pd.DataFrame):
        raise TypeError("mds_table must be a pandas DataFrame")
    missing = required.difference(mds_table.columns)
    if missing:
        raise ValueError(
            f"mds_table must contain columns: {', '.join(sorted(required))}"
        )
    if mds_table["sample_id"].duplicated().any():
        raise ValueError("mds_table must contain one row per sample_id")

    group_columns = _as_list(group, "group", max_len=1)
    split_columns = _as_list(split, "split", max_len=2)
    data, available_columns, _ = _merge_stats_metadata(mds_table, metadata)
    _validate_metadata_columns(group_columns, available_columns, "group")
    _validate_metadata_columns(split_columns, available_columns, "split")
    for coordinate in ["MDS1", "MDS2"]:
        data[coordinate] = pd.to_numeric(data[coordinate], errors="coerce")
    if data[["MDS1", "MDS2"]].isna().any().any():
        raise ValueError("MDS1 and MDS2 must contain only finite numeric values")
    if not np.isfinite(data[["MDS1", "MDS2"]].to_numpy()).all():
        raise ValueError("MDS1 and MDS2 must contain only finite numeric values")

    group_column = group_columns[0] if group_columns else "_beta_mds_group"
    if not group_columns:
        data[group_column] = "Samples"
    group_order = _category_order(data[group_column])
    if not group_order:
        raise ValueError("group column must contain at least one non-missing value")
    colors = sns.color_palette(palette, n_colors=max(1, len(group_order)))
    color_map = dict(zip(group_order, colors))

    panel_column = None
    if len(split_columns) == 1:
        panel_column = split_columns[0]
    elif len(split_columns) == 2:
        panel_column = "_split_panel"
        data[panel_column] = (
            data[split_columns[0]].astype(str)
            + " | "
            + data[split_columns[1]].astype(str)
        )
        data[panel_column] = pd.Categorical(
            data[panel_column],
            categories=_interaction_order(data, split_columns),
            ordered=True,
        )
    panel_order = _category_order(data[panel_column]) if panel_column else [None]
    if not panel_order:
        raise ValueError("split columns must contain at least one observed panel")

    column_count = min(3, len(panel_order))
    row_count = int(np.ceil(len(panel_order) / column_count))
    fig, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(height * aspect * column_count, height * row_count),
        squeeze=False,
    )
    axes = list(axes.flat)
    x_limits = _padded_limits(data["MDS1"])
    y_limits = _padded_limits(data["MDS2"])

    for axis, panel_value in zip(axes, panel_order):
        panel_data = data if panel_column is None else data[data[panel_column] == panel_value]
        for group_value in group_order:
            group_data = panel_data[panel_data[group_column] == group_value]
            if group_data.empty:
                continue
            color = color_map[group_value]
            centroid = group_data[["MDS1", "MDS2"]].mean().to_numpy(dtype=float)
            if dispersion:
                _draw_beta_mds_ellipse(axis, group_data, centroid, color)
            if centroids:
                for point in group_data[["MDS1", "MDS2"]].to_numpy(dtype=float):
                    axis.plot(
                        [centroid[0], point[0]],
                        [centroid[1], point[1]],
                        color=color,
                        linewidth=0.6,
                        alpha=0.6,
                        zorder=1,
                    )
            axis.scatter(
                group_data["MDS1"],
                group_data["MDS2"],
                color=color,
                edgecolor="black",
                linewidth=0.5,
                s=38,
                zorder=3,
            )
            if centroids:
                axis.scatter(
                    [centroid[0]],
                    [centroid[1]],
                    marker="X",
                    facecolor="white",
                    edgecolor=color,
                    linewidth=1.5,
                    s=100,
                    zorder=4,
                )
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
        axis.set_xlabel("MDS1")
        axis.set_ylabel("MDS2")
        if panel_column is not None:
            axis.set_title(str(panel_value))

    for axis in axes[len(panel_order):]:
        axis.set_visible(False)
    if group_columns:
        handles = [
            Line2D(
                [], [], marker="o", linestyle="", markerfacecolor=color_map[value],
                markeredgecolor="black", label=str(value), markersize=6,
            )
            for value in group_order
        ]
        fig.legend(
            handles=handles,
            title=_caption(group_column),
            loc="upper center",
            ncol=len(handles),
        )
        fig.tight_layout(rect=(0, 0, 1, 0.9))
    else:
        fig.tight_layout()
    return _close_and_return(fig)


def _padded_limits(values):
    minimum = float(values.min())
    maximum = float(values.max())
    span = maximum - minimum
    padding = span * 0.08 if span else max(abs(minimum) * 0.08, 0.1)
    return minimum - padding, maximum + padding


def _draw_beta_mds_ellipse(axis, group_data, centroid, color):
    if len(group_data) < 2:
        return
    covariance = np.cov(group_data[["MDS1", "MDS2"]].to_numpy(dtype=float), rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = np.maximum(eigenvalues, 0)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    ellipse = Ellipse(
        centroid,
        width=2 * np.sqrt(eigenvalues[0]),
        height=2 * np.sqrt(eigenvalues[1]),
        angle=angle,
        facecolor=(*to_rgb(color), 0.2),
        edgecolor=color,
        linewidth=1.2,
        zorder=0,
    )
    axis.add_patch(ellipse)


def _pairing_negative_log10_values(matrix):
    values = matrix.astype(float)
    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if (finite < 0).any():
        raise ValueError("-log10 pairing values cannot contain negatives")
    positive = finite[finite > 0]
    floor = (
        float(np.min(positive)) / 10
        if len(positive)
        else np.finfo(float).tiny
    )
    return -np.log10(values.mask(values == 0, floor))


def de_pairing(
    pairing_matrix,
    log_minus=False,
    hclust=False,
    show_values=True,
    cmap=PHEATMAP_CMAP,
    height=8,
    aspect=1.0,
):
    """Plot a chain-pairing score matrix as a heatmap.

    Parameters
    ----------
    pairing_matrix : pandas.DataFrame
        Numeric score matrix returned by ``diff_enrichment.pair_chains``.
    log_minus : bool, default False
        Color cells by ``-log10(value)``. Zeros are replaced by one tenth of
        the smallest positive score before transformation. Original scores
        remain displayed in cells.
    hclust : bool, default False
        Hierarchically cluster rows and columns.
    show_values : bool, default True
        Display compact original values, including ``NA``, in heatmap cells.
    cmap : matplotlib colormap, optional
        Heatmap colormap. Defaults to the R ``pheatmap`` palette.
    height, aspect : float
        Figure height and width multiplier.

    Returns
    -------
    matplotlib.figure.Figure
        A closed pairing heatmap figure.
    """
    if not isinstance(pairing_matrix, pd.DataFrame):
        raise TypeError("pairing_matrix must be a pandas DataFrame")
    if pairing_matrix.empty:
        raise ValueError("pairing_matrix must not be empty")
    original_values = pairing_matrix.apply(pd.to_numeric, errors="coerce")
    invalid = pairing_matrix.notna() & original_values.isna()
    if invalid.any().any():
        raise ValueError("pairing_matrix must contain only numeric values")
    color_values = original_values.copy()
    if log_minus:
        color_values = _pairing_negative_log10_values(color_values)

    row_linkage = _beta_linkage(color_values, axis=0) if hclust else None
    column_linkage = _beta_linkage(color_values, axis=1) if hclust else None
    method = pairing_matrix.attrs.get("method")
    method_label = str(method).upper() if method else "Pairing score"
    colorbar_label = f"-log10({method_label})" if log_minus else method_label
    grid = sns.clustermap(
        color_values,
        cmap=cmap,
        mask=color_values.isna(),
        row_cluster=row_linkage is not None,
        col_cluster=column_linkage is not None,
        row_linkage=row_linkage,
        col_linkage=column_linkage,
        figsize=(height * aspect, height),
        cbar_kws={"label": colorbar_label},
        xticklabels=True,
        yticklabels=True,
    )
    grid.ax_heatmap.set_title("Chain pairing")
    grid.ax_heatmap.set_xlabel(pairing_matrix.columns.name or "Chain 1 feature")
    grid.ax_heatmap.set_ylabel(pairing_matrix.index.name or "Chain 2 feature")
    grid.ax_heatmap.tick_params(axis="x", labelrotation=90)

    if show_values:
        row_order = (
            grid.dendrogram_row.reordered_ind
            if grid.dendrogram_row is not None
            else list(range(len(pairing_matrix.index)))
        )
        column_order = (
            grid.dendrogram_col.reordered_ind
            if grid.dendrogram_col is not None
            else list(range(len(pairing_matrix.columns)))
        )
        displayed = original_values.iloc[row_order, column_order]
        for row_index, row in enumerate(displayed.to_numpy(dtype=float)):
            for column_index, value in enumerate(row):
                grid.ax_heatmap.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    _format_beta_value(value),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                )

    return _close_and_return(grid.fig)


def _adjust_de_volcano_log2fc(values):
    adjusted = pd.Series(values, copy=True, dtype=float)
    finite = adjusted[np.isfinite(adjusted)]
    repeated_high_values = finite.value_counts()
    repeated_high_values = repeated_high_values[
        (repeated_high_values.index > 10) & (repeated_high_values > 1)
    ]
    for repeated_value in repeated_high_values.index:
        values_below = finite[finite < repeated_value]
        if values_below.empty:
            continue
        max_non_repeated_value = float(values_below.max())
        replacement = min(max_non_repeated_value + 2, float(repeated_value))
        adjusted.loc[adjusted.eq(repeated_value)] = replacement
    return adjusted


def _prepare_de_volcano_data(statistics_table, p_column):
    if not isinstance(statistics_table, pd.DataFrame):
        raise TypeError("statistics_table must be a pandas DataFrame")
    if p_column not in {"p_adj", "p_val"}:
        raise ValueError("p_column must be either 'p_adj' or 'p_val'")
    required_columns = {
        "log2FC",
        "mean_group_count",
        "enriched_in",
        p_column,
    }
    missing_columns = required_columns.difference(statistics_table.columns)
    if missing_columns:
        raise ValueError(
            "statistics_table is missing required volcano columns: "
            f"{sorted(missing_columns)}"
        )

    plot_columns = ["log2FC", p_column, "mean_group_count", "enriched_in"]
    for pass_column in ["prefilter_pass", "postfilter_pass"]:
        if pass_column in statistics_table:
            plot_columns.append(pass_column)
    plotted = statistics_table[plot_columns].dropna(
        subset=["log2FC", p_column, "mean_group_count", "enriched_in"]
    ).copy()
    if plotted.empty:
        raise ValueError("No rows with complete volcano-plot values are available")
    for column in ["log2FC", p_column, "mean_group_count"]:
        converted = pd.to_numeric(plotted[column], errors="coerce")
        if converted.isna().any():
            raise ValueError(f"statistics_table[{column!r}] must be numeric")
        plotted[column] = converted.astype(float)
    if ((plotted[p_column] < 0) | (plotted[p_column] > 1)).any():
        raise ValueError(f"statistics_table[{p_column!r}] values must be in [0, 1]")
    if (plotted["mean_group_count"] < 0).any():
        raise ValueError("statistics_table['mean_group_count'] must be non-negative")

    plotted["_plot_log2FC"] = _adjust_de_volcano_log2fc(
        plotted["log2FC"]
    ).to_numpy()
    p_values = plotted[p_column].to_numpy(dtype=float)
    positive_p_values = p_values[p_values > 0]
    p_floor = (
        float(np.min(positive_p_values)) / 10
        if len(positive_p_values)
        else np.finfo(float).tiny
    )
    plotted["_negative_log10_p"] = -np.log10(
        np.where(p_values == 0, p_floor, p_values)
    )
    return plotted


def de_volcano(
    statistics_table,
    p_column="p_adj",
    group_palette=None,
    size_range=(20, 300),
    log_sizes=True,
    alpha=0.7,
    height=6,
    aspect=1.3,
    by_mean_count=False,
):
    """Plot differential-enrichment effect sizes against p-values.

    Parameters
    ----------
    statistics_table : pandas.DataFrame
        A table returned by ``diff_enrichment.calc_statistics``.
    p_column : {"p_adj", "p_val"}, default "p_adj"
        P-value column used for the vertical axis.
    group_palette : palette name, sequence, or dict, optional
        Colors assigned to ``enriched_in`` groups.
    size_range : pair of float, default (20, 300)
        Minimum and maximum scatter-point areas.
    log_sizes : bool, default True
        Scale point sizes by ``log1p(mean_group_count)``. If false, use raw
        ``mean_group_count`` values. When ``by_mean_count=True``, apply the
        same transformation to the mean-count vertical axis instead.
    alpha : float, default 0.7
        Point opacity.
    height, aspect : float
        Figure height and width multiplier.
    by_mean_count : bool, default False
        Plot ``mean_group_count`` on the vertical axis instead of
        ``-log10(p)``. In this mode point sizes are scaled by the selected
        p-value column's ``-log10`` values.

    Returns
    -------
    matplotlib.figure.Figure
        Closed volcano-plot figure.
    """
    plotted = _prepare_de_volcano_data(statistics_table, p_column)
    filtered_out = pd.Series(False, index=plotted.index)
    if {"prefilter_pass", "postfilter_pass"}.issubset(plotted.columns):
        filtered_out = plotted["prefilter_pass"].eq(True) & plotted[
            "postfilter_pass"
        ].eq(False)
    retained = plotted.loc[~filtered_out]
    group_order = _category_order(retained["enriched_in"])
    color_map = _palette_mapping(group_order, palette=group_palette)
    if by_mean_count:
        size_values = plotted["_negative_log10_p"].to_numpy(dtype=float)
    else:
        size_values = plotted["mean_group_count"].to_numpy(dtype=float)
    if log_sizes and not by_mean_count:
        size_values = np.log1p(size_values)
    point_sizes = _scaled_dot_sizes(size_values, size_range)

    fig, ax = plt.subplots(figsize=(height * aspect, height))
    if by_mean_count:
        y_values = plotted["mean_group_count"].copy()
        if log_sizes:
            y_values = np.log1p(y_values)
    else:
        y_values = plotted["_negative_log10_p"]
    if filtered_out.any():
        ax.scatter(
            plotted.loc[filtered_out, "_plot_log2FC"],
            y_values.loc[filtered_out],
            s=float(np.min(point_sizes)),
            c="#999999",
            alpha=0.3,
            edgecolors="none",
            linewidths=0,
            zorder=1,
        )
    if not retained.empty:
        ax.scatter(
            retained["_plot_log2FC"],
            y_values.loc[retained.index],
            s=point_sizes[~filtered_out.to_numpy()],
            c=[color_map[group] for group in retained["enriched_in"]],
            alpha=float(alpha),
            edgecolors="black",
            linewidths=0.5,
            zorder=2,
        )
    ax.axvline(0, color="#888888", linestyle="--", linewidth=0.8, zorder=0)
    ax.set_title("Differential enrichment volcano plot")
    ax.set_xlabel("log2FC")
    ax.set_ylabel(
        (
            "log1p(Mean group count)"
            if by_mean_count and log_sizes
            else "Mean group count"
        )
        if by_mean_count
        else f"-log10({p_column})"
    )
    ax.grid(color="#eeeeee", linewidth=0.6)
    ax.set_axisbelow(True)

    handles = [
        Patch(facecolor=color_map[group], edgecolor="black", label=str(group))
        for group in group_order
    ]
    if filtered_out.any():
        handles.append(
            Patch(
                facecolor="#999999",
                edgecolor="none",
                alpha=0.3,
                label="filtered_out",
            )
        )
    if handles:
        ax.legend(
            handles=handles,
            title="Enriched in",
            frameon=False,
            loc="best",
        )
    fig.tight_layout()
    return _close_and_return(fig)


def beta_metric(
    beta_data,
    metric=None,
    metadata=None,
    group=None,
    hclust=True,
    ignore_diagonal=False,
    log_values=False,
    show_values=True,
    cmap=PHEATMAP_CMAP,
    height=8,
    aspect=1.0,
):
    """Plot a beta-diversity metric matrix as an annotated heatmap.

    Parameters
    ----------
    beta_data : dict or pandas.DataFrame
        Dictionary returned by ``beta.metrics`` or a numeric metric matrix.
    metric : str, optional
        Dictionary key to plot. Ignored when ``beta_data`` is a matrix.
    metadata : pandas.DataFrame, optional
        Sample metadata containing ``sample_id``.
    group : str or sequence of up to three str, optional
        Metadata columns displayed as row and column annotation strips.
    hclust : bool, default True
        Hierarchically cluster rows and columns.
    ignore_diagonal : bool, default False
        Replace cells whose row and column sample names match with ``NA``.
    log_values : bool, default False
        Color cells by log10 values after replacing zeros with a small floor.
    show_values : bool, default True
        Display compact original values, including ``NA``, in heatmap cells.
    cmap : matplotlib colormap, optional
        Heatmap colormap. Defaults to the R ``pheatmap`` palette.
    height, aspect : float
        Figure height and width multiplier.

    Returns
    -------
    matplotlib.figure.Figure or None
        Closed heatmap figure, or ``None`` when a dictionary metric is omitted.
    """
    matrix, title = _select_beta_metric_matrix(beta_data, metric)
    if matrix is None:
        return None

    group_columns = _as_list(group, "group", max_len=3)
    original_values = matrix.copy()
    if ignore_diagonal:
        shared_samples = matrix.index.intersection(matrix.columns)
        for sample in shared_samples:
            original_values.loc[sample, sample] = np.nan
    color_values = original_values.copy()
    if log_values:
        color_values = _beta_log_values(color_values, base=10)

    row_colors, row_handles = _prepare_beta_annotation_colors(
        metadata,
        group_columns,
        list(matrix.index),
    )
    column_colors, column_handles = _prepare_beta_annotation_colors(
        metadata,
        group_columns,
        list(matrix.columns),
    )
    row_linkage = _beta_linkage(color_values, axis=0) if hclust else None
    column_linkage = _beta_linkage(color_values, axis=1) if hclust else None
    grid = sns.clustermap(
        color_values,
        cmap=cmap,
        mask=color_values.isna(),
        row_cluster=row_linkage is not None,
        col_cluster=column_linkage is not None,
        row_linkage=row_linkage,
        col_linkage=column_linkage,
        row_colors=row_colors,
        col_colors=column_colors,
        figsize=(height * aspect, height),
        cbar_kws={"label": f"log10({title})" if log_values else title},
        xticklabels=True,
        yticklabels=True,
    )
    grid.ax_heatmap.set_title(title)
    grid.ax_heatmap.set_xlabel(matrix.columns.name or "Sample")
    grid.ax_heatmap.set_ylabel(matrix.index.name or "Sample")
    grid.ax_heatmap.tick_params(axis="x", labelrotation=90)

    if show_values:
        row_order = (
            grid.dendrogram_row.reordered_ind
            if grid.dendrogram_row is not None
            else list(range(len(matrix.index)))
        )
        column_order = (
            grid.dendrogram_col.reordered_ind
            if grid.dendrogram_col is not None
            else list(range(len(matrix.columns)))
        )
        displayed = original_values.iloc[row_order, column_order]
        for row_index, row in enumerate(displayed.to_numpy(dtype=float)):
            for column_index, value in enumerate(row):
                grid.ax_heatmap.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    _format_beta_value(value),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                )

    handles = {}
    for handle in [*row_handles, *column_handles]:
        handles.setdefault(handle.get_label(), handle)
    if handles:
        grid.fig.legend(
            handles.values(),
            handles.keys(),
            loc="upper center",
            ncol=min(4, len(handles)),
            frameon=False,
        )
        grid.fig.subplots_adjust(top=0.9)
    return _close_and_return(grid.fig)


def _normalize_beta_full_table(beta_data):
    if isinstance(beta_data, dict):
        if "full_table" not in beta_data:
            raise ValueError("beta-diversity dictionary does not contain 'full_table'")
        table = beta_data["full_table"]
    elif isinstance(beta_data, pd.DataFrame):
        table = beta_data
    else:
        raise TypeError("beta_data must be a beta.metrics dictionary or DataFrame")
    if not isinstance(table, pd.DataFrame):
        raise TypeError("full_table must be a pandas DataFrame")
    required = {"sample1", "sample2"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError("full_table is missing columns: " + ", ".join(sorted(missing)))
    has_frequencies = {"sample1_freq", "sample2_freq"}.issubset(table.columns)
    has_counts = {"sample1_count", "sample2_count"}.issubset(table.columns)
    if not has_frequencies and not has_counts:
        raise ValueError(
            "full_table needs sample1_freq/sample2_freq or sample1_count/sample2_count"
        )

    attrs = dict(table.attrs)
    table = table.copy()
    table.attrs = attrs
    if has_frequencies:
        table["_freq1"] = pd.to_numeric(table["sample1_freq"], errors="coerce")
        table["_freq2"] = pd.to_numeric(table["sample2_freq"], errors="coerce")
    else:
        count1 = pd.to_numeric(table["sample1_count"], errors="coerce")
        count2 = pd.to_numeric(table["sample2_count"], errors="coerce")
        pair_columns = ["sample1", "sample2"]
        total1 = count1.groupby([table[column] for column in pair_columns]).transform("sum")
        total2 = count2.groupby([table[column] for column in pair_columns]).transform("sum")
        table["_freq1"] = count1.div(total1.where(total1 != 0, np.nan)).fillna(0)
        table["_freq2"] = count2.div(total2.where(total2 != 0, np.nan)).fillna(0)
    if table[["_freq1", "_freq2"]].isna().any().any():
        raise ValueError("full_table contains missing or non-numeric frequencies")
    if (table[["_freq1", "_freq2"]] < 0).any().any():
        raise ValueError("full_table frequencies cannot be negative")
    return table


def _beta_table_layout(table):
    sample_list = table.attrs.get("sample_list")
    sample_list2 = table.attrs.get("sample_list2")
    if sample_list is not None:
        rows = list(sample_list)
        if sample_list2 is None:
            return rows, rows, True
        return rows, list(sample_list2), False

    sample1_order = list(pd.unique(table["sample1"]))
    sample2_order = list(pd.unique(table["sample2"]))
    all_samples = list(dict.fromkeys([*sample1_order, *sample2_order]))
    unordered_pairs = {
        frozenset((sample1, sample2))
        for sample1, sample2 in table[["sample1", "sample2"]].itertuples(index=False)
        if sample1 != sample2
    }
    expected_pairs = len(all_samples) * (len(all_samples) - 1) // 2
    same_set = bool(expected_pairs and len(unordered_pairs) == expected_pairs)
    if same_set:
        return all_samples, all_samples, True
    return sample1_order, sample2_order, False


def _beta_pair_frequencies(table, first_sample, second_sample):
    direct = table.loc[
        (table["sample1"] == first_sample) & (table["sample2"] == second_sample)
    ]
    if not direct.empty:
        return direct, direct["_freq1"].to_numpy(float), direct["_freq2"].to_numpy(float)
    reverse = table.loc[
        (table["sample1"] == second_sample) & (table["sample2"] == first_sample)
    ]
    if not reverse.empty:
        return reverse, reverse["_freq2"].to_numpy(float), reverse["_freq1"].to_numpy(float)
    return None, np.array([], dtype=float), np.array([], dtype=float)


def _beta_pair_f2(first_values, second_values):
    if len(first_values) == 0:
        return np.nan
    return float(np.sum(np.sqrt(first_values * second_values)))


def _draw_beta_dots_panel(
    ax,
    table,
    x_sample,
    y_sample,
    log_scale,
    log_base,
):
    pair, x, y = _beta_pair_frequencies(table, x_sample, y_sample)
    if pair is None:
        ax.set_visible(False)
        return
    if log_scale:
        positive = np.concatenate([x[x > 0], y[y > 0]])
        if len(positive) == 0:
            floor = 1 / float(log_base)
        else:
            floor = float(np.min(positive)) / float(log_base)
        x = np.where(x == 0, floor, x)
        y = np.where(y == 0, floor, y)
        lower = floor
        upper = max(float(np.max(x)), float(np.max(y)), floor * float(log_base))
        ax.set_xscale("log", base=log_base)
        ax.set_yscale("log", base=log_base)
    else:
        lower = 0.0
        upper = max(float(np.max(x)), float(np.max(y)), 1e-12)
    if log_scale:
        padding_factor = float(log_base) ** 0.05
        axis_lower = lower / padding_factor
        axis_upper = upper * padding_factor
    else:
        padding = upper * 0.05
        axis_lower = lower - padding
        axis_upper = upper + padding
    ax.plot(
        [lower, upper],
        [lower, upper],
        color="#888888",
        linestyle="--",
        linewidth=0.8,
        zorder=0,
    )
    ax.scatter(
        x,
        y,
        alpha=0.5,
        facecolors="#4e79a7",
        edgecolors="black",
        linewidths=0.5,
        s=24,
        zorder=1,
    )
    ax.set_xlim(axis_lower, axis_upper)
    ax.set_ylim(axis_lower, axis_upper)
    ax.set_xlabel(str(x_sample))
    ax.set_ylabel(str(y_sample))
    ax.grid(color="#eeeeee", linewidth=0.5)
    ax.set_axisbelow(True)


def _plot_beta_dots(
    table,
    row_samples,
    column_samples,
    same_set,
    log_scale,
    log_base,
    height,
    aspect,
    matrix_layout,
):
    if log_base <= 1:
        raise ValueError("log_base must be greater than 1")
    if not matrix_layout:
        if same_set:
            pairs = [
                (row_samples[first], row_samples[second])
                for first in range(len(row_samples))
                for second in range(first + 1, len(row_samples))
            ]
        else:
            pairs = [
                (row_sample, column_sample)
                for row_sample in row_samples
                for column_sample in column_samples
            ]
        if not pairs:
            raise ValueError("beta_table dots plots require at least one sample pair")
        columns = min(3, int(np.ceil(np.sqrt(len(pairs)))))
        rows = int(np.ceil(len(pairs) / columns))
        fig, axes = plt.subplots(
            rows,
            columns,
            figsize=(height * aspect * columns, height * rows),
            squeeze=False,
        )
        for ax, (row_sample, column_sample) in zip(axes.flat, pairs):
            _draw_beta_dots_panel(
                ax,
                table,
                column_sample,
                row_sample,
                log_scale,
                log_base,
            )
            ax.set_title(f"{row_sample} vs {column_sample}")
        for ax in list(axes.flat)[len(pairs):]:
            ax.set_visible(False)
        fig.tight_layout()
        return fig

    fig, axes = plt.subplots(
        len(row_samples),
        len(column_samples),
        figsize=(height * aspect * len(column_samples), height * len(row_samples)),
        squeeze=False,
    )
    for row_index, row_sample in enumerate(row_samples):
        for column_index, column_sample in enumerate(column_samples):
            ax = axes[row_index, column_index]
            if same_set and row_index == column_index:
                ax.text(0.5, 0.5, str(row_sample), ha="center", va="center", fontsize=11)
                ax.axis("off")
            elif same_set and row_index < column_index:
                _, first, second = _beta_pair_frequencies(table, row_sample, column_sample)
                ax.text(
                    0.5,
                    0.5,
                    "F2\n" + _format_beta_value(_beta_pair_f2(first, second)),
                    ha="center",
                    va="center",
                    fontsize=10,
                )
                ax.axis("off")
            else:
                _draw_beta_dots_panel(
                    ax,
                    table,
                    column_sample,
                    row_sample,
                    log_scale,
                    log_base,
                )
    fig.tight_layout()
    return fig


def _beta_clonotype_labels(pair):
    for column_name in ["cdr3aa", "cdr3nt", "clonotype", "clone"]:
        column = _column_by_name(pair.columns, {column_name})
        if column is None:
            continue
        return pair[column].map(
            lambda value: str(value[0])
            if isinstance(value, tuple) and value
            else str(value)
        ).to_numpy(dtype=object)
    return np.asarray(
        [f"Clonotype {index + 1}" for index in range(len(pair))],
        dtype=object,
    )


def _draw_beta_abundance_label(
    ax,
    label,
    lower,
    upper,
    first_value,
    second_value,
    fontsize,
):
    side = 0 if first_value >= second_value else 1
    if upper[side] <= lower[side]:
        return
    ax.text(
        0.015 if side == 0 else 0.985,
        (lower[side] + upper[side]) / 2,
        str(label),
        ha="left" if side == 0 else "right",
        va="center",
        fontsize=fontsize,
        color="black",
        clip_on=True,
    )


def _draw_beta_diff_panel(ax, table, first_sample, second_sample, top):
    pair, first, second = _beta_pair_frequencies(table, first_sample, second_sample)
    if pair is None:
        ax.set_visible(False)
        return

    shared = (first > 0) & (second > 0)
    shared_indices = np.flatnonzero(shared)
    shared_scores = np.sqrt(first[shared_indices] * second[shared_indices])
    ranked_shared = shared_indices[
        np.argsort(-shared_scores, kind="stable")
    ]
    shown_indices = list(ranked_shared[:top])
    shown_set = set(shown_indices)
    hidden_shared = np.asarray(
        [index for index in shared_indices if index not in shown_set],
        dtype=int,
    )
    non_overlapping = ~shared

    bands = [
        (
            "NonOverlapping",
            float(first[non_overlapping].sum()),
            float(second[non_overlapping].sum()),
            "#bdbdbd",
            8,
        ),
        (
            "NotShown",
            float(first[hidden_shared].sum()) if len(hidden_shared) else 0.0,
            float(second[hidden_shared].sum()) if len(hidden_shared) else 0.0,
            "#777777",
            8,
        ),
    ]
    clonotype_labels = _beta_clonotype_labels(pair)
    color_map = {
        index: RAREFACTION_COLORS_20[color_index]
        for color_index, index in enumerate(shown_indices)
    }
    for index in reversed(shown_indices):
        bands.append(
            (
                clonotype_labels[index],
                float(first[index]),
                float(second[index]),
                color_map[index],
                6,
            )
        )

    lower = np.zeros(2, dtype=float)
    for label, first_value, second_value, color, fontsize in bands:
        values = np.asarray([first_value, second_value], dtype=float)
        if not values.any():
            continue
        upper = lower + values
        ax.fill(
            [0, 0, 1, 1],
            [lower[0], upper[0], upper[1], lower[1]],
            facecolor=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=0.9,
        )
        _draw_beta_abundance_label(
            ax,
            label,
            lower,
            upper,
            first_value,
            second_value,
            fontsize,
        )
        lower = upper

    upper_limit = max(float(lower.max()), 1.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, upper_limit)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([str(first_sample), str(second_sample)], rotation=25, ha="right")
    ax.set_ylabel("Cumulative abundance")
    ax.grid(False)


def _plot_beta_diff(table, row_samples, column_samples, same_set, top, height, aspect):
    if not 1 <= int(top) <= len(RAREFACTION_COLORS_20):
        raise ValueError(
            f"top must be between 1 and {len(RAREFACTION_COLORS_20)}"
        )
    top = int(top)
    if same_set:
        pairs = [
            (row_samples[first], row_samples[second])
            for first in range(len(row_samples))
            for second in range(first + 1, len(row_samples))
        ]
        if not pairs:
            raise ValueError("beta_table diff plots require at least two samples")
        columns = min(3, int(np.ceil(np.sqrt(len(pairs)))))
        rows = int(np.ceil(len(pairs) / columns))
        fig, axes = plt.subplots(
            rows,
            columns,
            figsize=(height * aspect * columns, height * rows),
            squeeze=False,
        )
        for ax, pair_samples in zip(axes.flat, pairs):
            _draw_beta_diff_panel(ax, table, pair_samples[0], pair_samples[1], top)
            ax.set_title(f"{pair_samples[0]} vs {pair_samples[1]}")
        for ax in list(axes.flat)[len(pairs):]:
            ax.set_visible(False)
    else:
        fig, axes = plt.subplots(
            len(row_samples),
            len(column_samples),
            figsize=(height * aspect * len(column_samples), height * len(row_samples)),
            squeeze=False,
        )
        for row_index, row_sample in enumerate(row_samples):
            for column_index, column_sample in enumerate(column_samples):
                ax = axes[row_index, column_index]
                _draw_beta_diff_panel(ax, table, row_sample, column_sample, top)
                ax.set_title(f"{row_sample} vs {column_sample}")
    fig.tight_layout()
    return fig


def beta_table(
    beta_data,
    plot_type="dots",
    log_scale=False,
    log_base=10,
    top=20,
    height=3.2,
    aspect=1.0,
    matrix_layout=False,
):
    """Plot pairwise clonotype frequencies from beta-diversity full tables.

    Parameters
    ----------
    beta_data : dict or pandas.DataFrame
        Dictionary returned by ``beta.metrics`` or its ``full_table`` directly.
    plot_type : {"dots", "diff"}, default "dots"
        Scatterplot matrix or VDJtools-style shared-abundance plots.
    log_scale : bool, default False
        Use logarithmic axes for ``dots`` plots.
    log_base : float, default 10
        Logarithm base and zero-frequency floor divisor.
    top : int, default 20
        Number of shared clonotypes ranked by geometric-mean frequency and
        displayed individually in ``diff`` plots.
    height, aspect : float
        Per-panel height and width multiplier.
    matrix_layout : bool, default False
        For ``dots`` plots, use the legacy sample-by-sample matrix with F2
        values when True. By default, show pair panels in a wrapped facet
        layout without F2 panels. This option does not affect ``diff`` plots.

    Returns
    -------
    matplotlib.figure.Figure
        Closed figure that renders once in Jupyter.
    """
    table = _normalize_beta_full_table(beta_data)
    row_samples, column_samples, same_set = _beta_table_layout(table)
    plot_type = str(plot_type).casefold()
    if plot_type == "dots":
        fig = _plot_beta_dots(
            table,
            row_samples,
            column_samples,
            same_set,
            bool(log_scale),
            float(log_base),
            float(height),
            float(aspect),
            bool(matrix_layout),
        )
    elif plot_type == "diff":
        fig = _plot_beta_diff(
            table,
            row_samples,
            column_samples,
            same_set,
            top,
            float(height),
            float(aspect),
        )
    else:
        raise ValueError("plot_type must be one of: dots, diff")
    return _close_and_return(fig)


def rarefaction_curve(
    rarefaction_df,
    palette=None,
    height=4,
    aspect=1.4,
    marker="o",
    log_x=True,
):
    """Plot rarefaction curves from ``stats.calc_rarefaction_points`` output.

    Args:
        rarefaction_df (pd.DataFrame): Table with ``sample_id``,
            ``rarefaction_depth``, and ``diversity`` columns. If duplicated
            ``sample_id`` values have distinct ``chain`` values, labels are
            shown as ``sample_id(chain)``.
        palette: Seaborn palette or dict for sample curves. By default, a
            contrasting 20-color palette is used for up to 20 samples; larger
            plots use Seaborn's existing default palette behavior.
        height (float): Figure height.
        aspect (float): Width/height ratio.
        marker (str): Matplotlib marker for observed points.
        log_x (bool): Use logarithmic x-axis.

    Returns:
        matplotlib.figure.Figure: Rarefaction curve figure.
    """
    required = {"sample_id", "rarefaction_depth", "diversity"}
    missing = required.difference(rarefaction_df.columns)
    if missing:
        raise ValueError(
            "rarefaction_df must contain columns: " + ", ".join(sorted(required))
        )
    data = rarefaction_df.copy()
    if data.empty:
        raise ValueError("rarefaction_df is empty")
    data["_sample_label"] = data["sample_id"].astype(str)
    if "chain" in data.columns:
        sample_chains = data[["sample_id", "chain"]].drop_duplicates()
        chained_sample_ids = sample_chains.groupby("sample_id").size()
        chained_sample_ids = set(chained_sample_ids[chained_sample_ids > 1].index)
        use_chain = data["sample_id"].isin(chained_sample_ids)
        data.loc[use_chain, "_sample_label"] = (
            data.loc[use_chain, "sample_id"].astype(str)
            + "("
            + data.loc[use_chain, "chain"].astype(str)
            + ")"
        )
    data = data.sort_values(["_sample_label", "rarefaction_depth"])
    sample_count = data["_sample_label"].nunique()
    if palette is None and sample_count <= len(RAREFACTION_COLORS_20):
        palette = RAREFACTION_COLORS_20[:sample_count]

    fig, ax = plt.subplots(figsize=(height * aspect, height))
    sns.lineplot(
        data=data,
        x="rarefaction_depth",
        y="diversity",
        hue="_sample_label",
        marker=marker,
        palette=palette,
        estimator=None,
        sort=True,
        ax=ax,
    )
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("Rarefaction depth")
    ax.set_ylabel("Observed diversity")
    legend = ax.get_legend()
    if legend is not None:
        ax.legend(
            title="Sample",
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            borderaxespad=0,
        )
    fig.tight_layout()
    return _close_and_return(fig)


def _coverage_bin_positions(data):
    bin_numbers = sorted(pd.unique(data["_bin_number"]))
    positions = {bin_number: position for position, bin_number in enumerate(bin_numbers)}
    labels = {
        bin_number: str(data.loc[data["_bin_number"] == bin_number, "bin"].iloc[0])
        for bin_number in bin_numbers
    }
    return positions, [labels[bin_number] for bin_number in bin_numbers]


def _coverage_sample_labels(data):
    labels = data["sample_id"].astype(str)
    if "chain" not in data.columns:
        return labels
    chain_counts = data[["sample_id", "chain"]].drop_duplicates().groupby(
        "sample_id"
    ).size()
    use_chain = data["sample_id"].map(chain_counts).gt(1)
    labels = labels.copy()
    labels.loc[use_chain] = (
        data.loc[use_chain, "sample_id"].astype(str)
        + " ("
        + data.loc[use_chain, "chain"].astype(str)
        + ")"
    )
    return labels


def _draw_clonotypes_coverage_lines(data, marker="o", **kwargs):
    ax = kwargs.get("ax", plt.gca())
    positions, labels = _coverage_bin_positions(data)
    plotted = data.assign(
        _bin_position=data["_bin_number"].map(positions)
    ).sort_values("_bin_position")
    sns.lineplot(
        data=plotted,
        x="_bin_position",
        y="value",
        estimator=None,
        marker=marker,
        sort=True,
        color=kwargs.get("color"),
        ax=ax,
    )
    ax.set_xticks(range(len(labels)), labels)


def _draw_clonotypes_coverage_bars(data, trim_high_zero_bins=True, **kwargs):
    ax = kwargs.get("ax", plt.gca())
    plotted = data.copy()
    if trim_high_zero_bins:
        nonzero_bins = plotted.loc[plotted["value"] != 0, "_bin_number"]
        highest_nonzero = nonzero_bins.max() if not nonzero_bins.empty else -np.inf
        plotted = plotted.loc[
            ~(
                (plotted["_bin_number"] > 100)
                & (plotted["_bin_number"] > highest_nonzero)
            )
        ]
    positions, labels = _coverage_bin_positions(plotted)
    plotted = plotted.assign(_bin_position=plotted["_bin_number"].map(positions))
    ax.bar(plotted["_bin_position"], plotted["value"], color="#CCCCCC")
    ax.set_xticks(range(len(labels)), labels)


def clonotypes_coverage(
    coverage_df,
    metadata=None,
    group=None,
    split=None,
    separate=False,
    trim_high_zero_bins=True,
    palette=None,
    height=3.2,
    aspect=1.3,
    marker="o",
):
    """Plot half-order clonotype-count histograms.

    With ``separate=False``, samples are drawn as lines on the same panel and
    may be split into facets by at most two metadata columns. Grouping is not
    supported for this plot type. With ``separate=True``, metadata, ``group``,
    and ``split`` are ignored and each sample is drawn as a grey barplot in a
    separate facet.

    Args:
        coverage_df (pd.DataFrame): Output of
            :func:`repseq.stats.clonotypes_coverage` with ``sample_id``,
            ``bin``, and ``value`` columns.
        metadata (pd.DataFrame, optional): Sample metadata used only for
            splitting combined line plots.
        group (optional): Unsupported for combined plots and ignored for
            separate plots.
        split (str or sequence of str, optional): One or two facet columns.
        separate (bool): Draw one barplot facet per sample.
        trim_high_zero_bins (bool): In separate plots, remove zero-valued bins
            above both 100 and the highest non-zero bin in that panel.
        palette: Seaborn palette for sample lines.
        height, aspect (float): Facet dimensions.
        marker (str): Marker used for combined line plots.

    Returns:
        seaborn.FacetGrid: The created plot grid.
    """
    required = {"sample_id", "bin", "value"}
    missing = required.difference(coverage_df.columns)
    if missing:
        raise ValueError(
            "coverage_df must contain columns: " + ", ".join(sorted(required))
        )
    if coverage_df.empty:
        raise ValueError("coverage_df is empty")

    data = coverage_df.copy()
    try:
        data["_bin_number"] = pd.to_numeric(data["bin"])
    except (TypeError, ValueError) as error:
        raise ValueError("coverage_df 'bin' values must be numeric") from error
    data["value"] = pd.to_numeric(data["value"], errors="raise")
    data["_sample_label"] = _coverage_sample_labels(data)

    if separate:
        grid = sns.FacetGrid(
            data,
            col="_sample_label",
            col_wrap=3,
            sharex=False,
            sharey=False,
            height=height,
            aspect=aspect,
            despine=True,
        )
        grid.map_dataframe(
            _draw_clonotypes_coverage_bars,
            trim_high_zero_bins=trim_high_zero_bins,
        )
        grid.set_titles("{col_name}")
    else:
        if group is not None:
            raise ValueError("group is not supported for clonotypes_coverage plots")
        split_columns = _as_list(split, "split", max_len=2)
        data, metadata_columns, _ = _merge_stats_metadata(data, metadata)
        _validate_metadata_columns(split_columns, metadata_columns, "split")
        data["_sample_label"] = _coverage_sample_labels(data)

        panel_column = None
        if len(split_columns) == 1:
            panel_column = split_columns[0]
        elif len(split_columns) == 2:
            panel_column = "_split_panel"
            data[panel_column] = (
                data[split_columns[0]].astype(str)
                + " | "
                + data[split_columns[1]].astype(str)
            )
            data[panel_column] = pd.Categorical(
                data[panel_column],
                categories=_interaction_order(data, split_columns),
                ordered=True,
            )

        grid = sns.FacetGrid(
            data,
            col=panel_column,
            col_wrap=3 if panel_column is not None else None,
            hue="_sample_label",
            palette=palette,
            sharex=True,
            sharey=False,
            height=height,
            aspect=aspect,
            despine=True,
        )
        grid.map_dataframe(_draw_clonotypes_coverage_lines, marker=marker)
        if panel_column is not None:
            grid.set_titles("{col_name}")
        grid.add_legend(title="Sample")

    grid.set_axis_labels("Clonotype size", "Value")
    grid.tight_layout()
    return grid


clonotype_coverage = clonotypes_coverage

def clonoset_stats(
    stats_df,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot clonoset size statistics with optional sample metadata.

    Grouped plots use the standard statistics renderer. Their default
    properties are ``reads``, ``reads_per_umi`` when present, ``clones_func``,
    and ``umi_func`` when it contains at least one value. Custom properties
    always use the standard renderer as well.

    Without a group or custom properties, reads, clones, and available UMI
    counts are drawn as equal-width overlaid total and functional bars. Labels
    above the
    bars show ``total(functional)``. ``reads_per_umi`` remains an ordinary bar
    plot, and samples without UMI counts are omitted from the UMI panel.
    """
    if group is not None or properties is not None:
        selected_properties = (
            _default_grouped_clonoset_properties(stats_df)
            if properties is None
            else properties
        )
        return plot_stats(
            stats_df,
            metadata=metadata,
            properties=selected_properties,
            group=group,
            split=split,
            palette=palette,
            height=height,
            aspect=aspect,
            zero_bottom=group is not None,
        )

    return _plot_default_clonoset_stats(
        stats_df,
        metadata=metadata,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
    )


def processing(
    processing_table,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot MiXCR processing statistics with optional sample metadata.

    The input is the table returned by :func:`repseq.mixcr.get_processing_table`.
    Its ``extracted_chain`` column is kept separate from a metadata ``chain``
    column and is used to distinguish extracted chains in sample labels.
    """
    return plot_stats(
        processing_table,
        metadata=metadata,
        properties=PROCESSING_PROPERTIES if properties is None else properties,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
        zero_bottom=True,
    )


def cdr3aa_stats(
    stats_df,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot CDR3 amino-acid property statistics."""
    return plot_stats(
        stats_df,
        metadata=metadata,
        properties=CDR3AA_STATS_PROPERTIES if properties is None else properties,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
    )


def diversity_stats(
    stats_df,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot diversity statistics."""
    return plot_stats(
        stats_df,
        metadata=metadata,
        properties=DIVERSITY_STATS_PROPERTIES if properties is None else properties,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
        zero_bottom=True,
    )


def convergence(
    stats_df,
    metadata=None,
    properties=None,
    group=None,
    split=None,
    palette=None,
    height=3.2,
    aspect=1.2,
):
    """Plot convergence statistics."""
    return plot_stats(
        stats_df,
        metadata=metadata,
        properties=CONVERGENCE_PROPERTIES if properties is None else properties,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
    )


def _identifier_key(value):
    if pd.isna(value):
        return None
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and np.isfinite(value) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    if re.fullmatch(r"[+-]?\d+\.0+", text):
        return text.split(".", 1)[0]
    return text


def _is_observed_value(value):
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def _restore_numeric_internal_node_names(tree, valid_node_ids=None):
    valid_node_ids = None if valid_node_ids is None else set(valid_node_ids)
    for clade in tree.get_nonterminals(order="preorder"):
        if clade.name is None and clade.confidence is not None:
            node_id = _identifier_key(clade.confidence)
            if valid_node_ids is None or node_id in valid_node_ids:
                clade.name = node_id
                clade.confidence = None
    return tree


def _phylo_positions(tree):
    x_positions = tree.depths()
    if not x_positions or max(x_positions.values()) == 0:
        x_positions = tree.depths(unit_branch_lengths=True)

    terminals = tree.get_terminals()
    y_positions = {
        terminal: float(len(terminals) - index)
        for index, terminal in enumerate(reversed(terminals))
    }

    def set_internal_y(clade):
        for child in clade.clades:
            if child not in y_positions:
                set_internal_y(child)
        if clade.clades:
            y_positions[clade] = (
                y_positions[clade.clades[0]] + y_positions[clade.clades[-1]]
            ) / 2

    set_internal_y(tree.root)
    return x_positions, y_positions


def draw_tree(trees_df, treeId, metadata=None, group=None, label=None, ax=None):
    """Draw a MiXCR SHM tree from a node table loaded by ``TreeAnalyzer``."""
    if "treeId" not in trees_df.columns:
        raise ValueError("trees_df must contain a 'treeId' column")
    tree_id_key = _identifier_key(treeId)
    tree_rows = trees_df.loc[
        trees_df["treeId"].map(_identifier_key) == tree_id_key
    ].copy()
    if tree_rows.empty:
        raise ValueError(f"treeId {treeId!r} is not present in trees_df")
    if "nodeId" not in tree_rows.columns:
        raise ValueError("trees_df must contain a 'nodeId' column")
    if "newick_filename" not in tree_rows.columns:
        raise ValueError("Newick tree filenames have not been attached to trees_df")
    filenames = tree_rows["newick_filename"].dropna().astype(str).unique()
    if len(filenames) == 0:
        raise ValueError(f"Newick tree {treeId!r} has not been attached")
    if len(filenames) > 1:
        raise ValueError(f"treeId {treeId!r} has multiple Newick filenames")
    newick_filename = Path(filenames[0])
    if not newick_filename.is_file():
        raise ValueError(f"Newick tree file does not exist: {newick_filename}")
    try:
        from Bio import Phylo
    except ImportError as error:
        raise ImportError(
            "draw_tree requires Biopython. Install repseq with the biopython dependency."
        ) from error
    valid_node_ids = set(tree_rows["nodeId"].map(_identifier_key).dropna())
    tree = Phylo.read(str(newick_filename), "newick")
    tree = _restore_numeric_internal_node_names(tree, valid_node_ids)
    node_data = tree_rows.assign(
        _node_key=tree_rows["nodeId"].map(_identifier_key)
    )
    if "sample_id" not in node_data.columns and "fileName" in node_data.columns:
        node_data["sample_id"] = node_data["fileName"].apply(
            lambda filename: str(filename).rsplit("/", 1)[-1][:-5].rsplit(".", 1)[-1]
            if str(filename).endswith(".clns") else str(filename)
        )
    if metadata is not None:
        if "sample_id" not in metadata.columns:
            raise ValueError("metadata must contain a 'sample_id' column")
        missing_metadata_columns = [
            column
            for column in metadata.columns
            if column == "sample_id" or column not in node_data.columns
        ]
        if len(missing_metadata_columns) > 1:
            node_data = node_data.merge(
                metadata.loc[:, missing_metadata_columns],
                on="sample_id",
                how="left",
            )

    if "isObserved" in node_data.columns:
        observed_mask = node_data["isObserved"].map(_is_observed_value)
    else:
        observed_mask = pd.Series(True, index=node_data.index, dtype=bool)
    read_counts = (
        pd.to_numeric(node_data["readCount"], errors="coerce")
        if "readCount" in node_data.columns
        else pd.Series(1, index=node_data.index, dtype=float)
    )
    node_data["_plot_count"] = (
        pd.to_numeric(node_data["uniqueMoleculeCount"], errors="coerce").fillna(read_counts)
        if "uniqueMoleculeCount" in node_data.columns else read_counts
    ).fillna(1)
    observed_counts = node_data.loc[observed_mask, "_plot_count"]
    max_count = max(float(observed_counts.max()) if not observed_counts.empty else 1, 1)
    x_positions, y_positions = _phylo_positions(tree)
    node_groups = {
        node_id: rows
        for node_id, rows in node_data.groupby("_node_key", sort=False)
    }

    color_values = (
        node_data.loc[observed_mask, group]
        if group is not None and group in node_data.columns else None
    )
    categories = [] if color_values is None else sorted(color_values.dropna().astype(str).unique())
    palette = dict(zip(categories, sns.color_palette(n_colors=len(categories))))

    if ax is None:
        _, ax = plt.subplots(figsize=(10, max(4, len(tree.get_terminals()) * 0.35)))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda _clade: None,
        show_confidence=False,
    )
    for clade in tree.find_clades(order="preorder"):
        x = x_positions[clade]
        y = y_positions[clade]
        node = _identifier_key(clade.name)
        rows = node_groups.get(node)
        if rows is None:
            continue
        if "isObserved" in rows.columns:
            row_observed = rows["isObserved"].map(_is_observed_value)
        else:
            row_observed = pd.Series(True, index=rows.index, dtype=bool)
        observed_rows = rows.loc[row_observed]
        non_observed_rows = rows.loc[~row_observed]

        if clade.is_terminal() and not observed_rows.empty:
            observed_rows = observed_rows.sort_values("_plot_count", ascending=False)
            label_count = len(observed_rows)
            for label_index, (_, row) in enumerate(observed_rows.iterrows()):
                category = str(row[group]) if group is not None and group in row and pd.notna(row[group]) else None
                color = palette.get(category, "0.45")
                size = 30 + 270 * np.sqrt(float(row["_plot_count"]) / max_count)
                ax.scatter(x, y, s=size, color=color, edgecolor="white", linewidth=0.6, zorder=3)
                annotations = [str(row["sample_id"])] if "sample_id" in row and pd.notna(row["sample_id"]) else []
                if label is not None and label in row and pd.notna(row[label]):
                    annotations.append(str(row[label]))
                if annotations:
                    label_y_offset = (label_index - (label_count - 1) / 2) * 11
                    ax.annotate(
                        " | ".join(annotations),
                        (x, y),
                        xytext=(6, label_y_offset),
                        textcoords="offset points",
                        fontsize=8,
                    )
        if not non_observed_rows.empty:
            ax.scatter(
                x,
                y,
                s=22,
                color="white",
                edgecolor="0.35",
                linewidth=0.8,
                zorder=3,
            )
    if categories:
        handles = [Line2D([], [], marker="o", linestyle="", color=palette[value], label=value) for value in categories]
        ax.legend(
            handles=handles,
            title=group,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.08),
            ncol=min(len(handles), 5),
        )
    ax.set_title(f"Tree {treeId}")
    ax.set_xlabel("Branch length")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    return ax


__all__ = [
    "CDR3AA_STATS_PROPERTIES",
    "DIVERSITY_STATS_PROPERTIES",
    "CONVERGENCE_PROPERTIES",
    "PROCESSING_PROPERTIES",
    "PHEATMAP_CMAP",
    "parse_gene_name",
    "plot_stats",
    "segment_usage",
    "cdr3_length_distribution",
    "cdr3_length_distributions",
    "vj_usage",
    "vjlen_usage",
    "beta_metric",
    "de_pairing",
    "beta_table",
    "rarefaction_curve",
    "clonotypes_coverage",
    "clonotype_coverage",
    "clonoset_stats",
    "processing",
    "cdr3aa_stats",
    "diversity_stats",
    "convergence",
    "draw_tree",
]
