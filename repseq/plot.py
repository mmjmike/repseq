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
from matplotlib.colors import to_rgb
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
            if column not in id_columns and parse_gene_name(str(column)) is not None
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

    data["_segment"] = data["_segment"].astype(str)
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


def _prepare_segment_usage_data(segment_usage_df, metadata, group, split, plot_type):
    data, segment_type = _normalize_segment_usage_table(segment_usage_df)
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
        ax.set_xlabel("Segment")
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
        heatmap_ax.set_xlabel("Segment")
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
    cmap="RdBu_r",
    height=3.2,
    aspect=1.2,
    seed=0,
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
    cmap : matplotlib colormap, default "RdBu_r"
        Colormap for usage values in heatmaps.
    height, aspect : float
        Base panel height and width multiplier.
    seed : int, default 0
        Random seed used for horizontal jitter in boxplots.

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
    )
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
        )

    # Inline backends display open pyplot figures and the returned object. A
    # closed Figure remains renderable and editable but appears only once.
    if fig is not None:
        plt.close(fig)
    return fig


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
    if combination_column is not None:
        required_columns = [v_column, j_column, value_column]
        if expected_length == 3:
            required_columns.append(length_column)
        if any(column is None for column in required_columns):
            raise ValueError(
                f"Long {combination_type} tables require pipe-delimited "
                f"{combination_type}, v, j"
                + (", len" if expected_length == 3 else "")
                + ", and a value column"
            )
        keep_columns = id_columns + [combination_column, v_column, j_column]
        if expected_length == 3:
            keep_columns.append(length_column)
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
        expected_identifiers = (
            data[v_column].astype(str) + "|" + data[j_column].astype(str)
        )
        if expected_length == 3:
            expected_identifiers += "|" + data[length_column].astype(str)
        if (data[combination_column].astype(str) != expected_identifiers).any():
            raise ValueError(
                f"{combination_type} values must match the separate v, j"
                + (", and len" if expected_length == 3 else "")
                + " columns"
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


def _prepare_combination_plot_data(usage_df, metadata, group, split, combination_type):
    data = _normalize_combination_usage_table(usage_df, combination_type)
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

    Returns
    -------
    matplotlib.figure.Figure
        A closed figure that renders once in Jupyter.
    """
    data, group_columns, split_columns = _prepare_combination_plot_data(
        usage_df, metadata, group, split, "vj"
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
            ax.set_xlabel("V segment")
            ax.set_ylabel("J segment")
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

    Returns
    -------
    matplotlib.figure.Figure
        A closed figure that renders once in Jupyter.
    """
    data, group_columns, split_columns = _prepare_combination_plot_data(
        usage_df, metadata, group, split, "vjlen"
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
    "parse_gene_name",
    "plot_stats",
    "segment_usage",
    "vj_usage",
    "vjlen_usage",
    "cdr3aa_stats",
    "diversity_stats",
    "convergence",
]
