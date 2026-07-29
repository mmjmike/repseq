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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Patch
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
    (?:\*(?P<allele>\d{2,3}))?
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
        _optional_number(parsed["allele"]),
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
        return stats_df.copy(), set(stats_df.columns), ["sample_id"]

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

    return merged.drop(columns="_merge"), set(metadata.columns), merge_keys


def _validate_metadata_columns(columns, metadata_columns, label):
    missing = [column for column in columns if column not in metadata_columns]
    if missing:
        raise ValueError(f"{label} column(s) must be present in metadata: {missing}")


def _make_sample_labels(data):
    use_chain = "chain" in data.columns and data["sample_id"].duplicated(keep=False).any()
    if use_chain:
        labels = data["sample_id"].astype(str) + " (" + data["chain"].astype(str) + ")"
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
    if metadata is None and (group_columns or split_columns):
        raise ValueError("metadata is required when group or split columns are used")
    _validate_metadata_columns(group_columns, metadata_columns, "group")
    _validate_metadata_columns(split_columns, metadata_columns, "split")

    missing_properties = [column for column in properties if column not in data.columns]
    if missing_properties:
        raise ValueError(f"stats_df does not contain property column(s): {missing_properties}")

    id_columns = list(dict.fromkeys(merge_keys + group_columns + split_columns))
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
):
    """Plot a wide statistics table with optional sample metadata.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        Statistics table. It must contain ``sample_id`` and the selected
        property columns. If both ``stats_df`` and ``metadata`` contain
        ``chain``, the merge uses both ``sample_id`` and ``chain``.
    metadata : pandas.DataFrame, optional
        Sample metadata with one row per ``sample_id`` or ``sample_id`` +
        ``chain`` key. Grouping and splitting columns must come from this
        table. If metadata covers only part of ``stats_df``, a warning is
        emitted and only matched samples are plotted.
    properties : str or sequence of str
        Numeric columns to plot. Multiple properties are shown as separate
        panels.
    group : str or sequence of str, optional
        Metadata column(s) used to group samples. No group draws one bar per
        sample. One group draws boxplots with jittered points by that group.
        Two groups use the first column on the x-axis and the second column as
        color.
    split : str or sequence of str, optional
        One or two metadata columns used to split plots into panels. Two split
        columns are combined into an interaction panel.
    palette : optional
        Any seaborn/matplotlib-compatible palette specification.
    height, aspect : float
        Facet size arguments passed to :class:`seaborn.FacetGrid`.

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
    grid.tight_layout()
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
    ax.set_xlim(lower, upper * (1.05 if not log_scale else float(log_base) ** 0.05))
    ax.set_ylim(lower, upper * (1.05 if not log_scale else float(log_base) ** 0.05))
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
):
    if log_base <= 1:
        raise ValueError("log_base must be greater than 1")
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
        properties=properties or CDR3AA_STATS_PROPERTIES,
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
        properties=properties or DIVERSITY_STATS_PROPERTIES,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
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
        properties=properties or CONVERGENCE_PROPERTIES,
        group=group,
        split=split,
        palette=palette,
        height=height,
        aspect=aspect,
    )


__all__ = [
    "CDR3AA_STATS_PROPERTIES",
    "DIVERSITY_STATS_PROPERTIES",
    "CONVERGENCE_PROPERTIES",
    "PHEATMAP_CMAP",
    "parse_gene_name",
    "plot_stats",
    "segment_usage",
    "cdr3_length_distribution",
    "cdr3_length_distributions",
    "vj_usage",
    "vjlen_usage",
    "beta_metric",
    "beta_table",
    "rarefaction_curve",
    "cdr3aa_stats",
    "diversity_stats",
    "convergence",
]
