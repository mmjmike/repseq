"""Plotting helpers for repseq statistics tables.

The functions in this module take wide statistics tables produced by
``repseq.stats`` and draw matplotlib/seaborn categorical summaries.  The
plotting API is intentionally dataframe-oriented so that sample metadata can be
joined immediately before plotting.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


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
    "plot_stats",
    "cdr3aa_stats",
    "diversity_stats",
    "convergence",
]
