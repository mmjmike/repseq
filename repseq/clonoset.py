from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from .common_functions import (
    decide_count_and_frac_columns,
    extract_segment,
    get_column_names_from_clonoset,
)


VDJTOOLS_COLUMNS = {
    "freq": "float64",
    "count": "int64",
    "cdr3nt": "string",
    "cdr3aa": "string",
    "v": "category",
    "d": "category",
    "j": "category",
}


def standardize_to_vdjtools_columns(
    df,
    *,
    by_umi=False,
    extract_segments=False,
    copy=False,
):
    """
    Return a clonoset table with canonical VDJtools-like column names.

    The canonical names are `freq`, `count`, `cdr3nt`, `cdr3aa`, `v`, `d`,
    and `j`. Existing canonical columns are preferred. If a canonical column
    is missing, the function renames the corresponding detected source column.

    Args:
        df (pd.DataFrame): Clonoset table in VDJtools, MiXCR, Adaptive, or
            AIRR-like format.
        by_umi (bool): If `True`, use UMI or molecule count/fraction columns
            for `count` and `freq` when they are present.
        extract_segments (bool): If `True`, remove allele and score parts from
            V/D/J values using the same parser as other repseq functions.
        copy (bool): If `True`, copy the input table before changing columns.

    Returns:
        pd.DataFrame: Table with canonical VDJtools-like column names.
    """
    clonoset = df.copy() if copy else df
    colnames = get_column_names_from_clonoset(clonoset)
    count_column, fraction_column = decide_count_and_frac_columns(
        colnames,
        by_umi,
        suppress_warnings=True,
    )

    source_by_target = {
        "count": count_column,
        "freq": fraction_column,
        "cdr3nt": colnames["cdr3nt_column"],
        "cdr3aa": colnames["cdr3aa_column"],
        "v": colnames["v_column"],
        "d": colnames["d_column"],
        "j": colnames["j_column"],
    }
    rename_dict = {
        source: target
        for target, source in source_by_target.items()
        if source is not None and source != target and target not in clonoset.columns
    }
    if rename_dict:
        clonoset = clonoset.rename(columns=rename_dict, copy=False)

    for column in ("cdr3aa", "cdr3nt"):
        if column in clonoset.columns:
            clonoset[column] = clonoset[column].fillna("")

    if extract_segments:
        for column in ("v", "d", "j"):
            if column in clonoset.columns:
                clonoset[column] = clonoset[column].apply(extract_segment)

    return clonoset


@dataclass
class Clonoset:
    """
    Lightweight wrapper around a pandas DataFrame containing clonotypes.

    `Clonoset` keeps the DataFrame as the main storage layer, so existing
    vectorized pandas operations remain fast. The wrapper carries sample-level
    attributes and provides small convenience methods while keeping access to
    the underlying table through `data` or `to_dataframe()`.

    Args:
        data (pd.DataFrame): Clonotype table.
        chain (str, optional): Receptor chain name, for example `TRA` or `TRB`.
        sample_id (str, optional): Sample identifier.
        metadata (dict, optional): Additional sample-level metadata.
        standardize (bool): Rename detected columns to canonical VDJtools-like
            names on construction.
        by_umi (bool): Use UMI or molecule columns as canonical `count`/`freq`
            when `standardize=True`.
        validate (bool): Check that `cdr3aa`, `v`, and `j` are available after
            optional standardization.

    Attributes:
        data (pd.DataFrame): Canonical or original clonoset table.
    """

    data: pd.DataFrame
    chain: Optional[str] = None
    sample_id: Optional[str] = None
    metadata: dict = field(default_factory=dict)
    standardize: bool = True
    by_umi: bool = False
    validate: bool = True

    def __post_init__(self):
        if self.standardize:
            self.data = standardize_to_vdjtools_columns(
                self.data,
                by_umi=self.by_umi,
                copy=False,
            )
        if self.validate:
            self._validate_minimal_schema()

    def _validate_minimal_schema(self):
        required = {"cdr3aa", "v", "j"}
        missing = required - set(self.data.columns)
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(f"Missing required clonoset columns: {missing_text}")

    @property
    def columns(self):
        """Columns of the underlying clonoset DataFrame."""
        return self.data.columns

    @property
    def shape(self):
        """Shape of the underlying clonoset DataFrame."""
        return self.data.shape

    def to_dataframe(self, copy=False):
        """
        Return the underlying pandas DataFrame.

        Args:
            copy (bool): If `True`, return a copy instead of the stored table.

        Returns:
            pd.DataFrame: Clonoset data.
        """
        return self.data.copy() if copy else self.data

    def copy(self, deep=True):
        """Return a copy of the clonoset object and its metadata."""
        return Clonoset(
            self.data.copy(deep=deep),
            chain=self.chain,
            sample_id=self.sample_id,
            metadata=self.metadata.copy(),
            standardize=False,
            by_umi=self.by_umi,
            validate=self.validate,
        )

    def filter(self, filter_obj):
        """Apply a repseq `Filter` object and return a new `Clonoset`."""
        return Clonoset(
            filter_obj.apply(self.data),
            chain=self.chain,
            sample_id=self.sample_id,
            metadata=self.metadata.copy(),
            standardize=False,
            by_umi=self.by_umi,
            validate=self.validate,
        )

    def top(self, n):
        """Return a new `Clonoset` with the `n` largest clonotypes by count."""
        if "count" not in self.data.columns:
            raise ValueError("Column 'count' is required for top()")
        return Clonoset(
            self.data.nlargest(n, "count"),
            chain=self.chain,
            sample_id=self.sample_id,
            metadata=self.metadata.copy(),
            standardize=False,
            by_umi=self.by_umi,
            validate=self.validate,
        )

    def normalize_freq(self):
        """Return a new `Clonoset` with `freq` recalculated from `count`."""
        if "count" not in self.data.columns:
            raise ValueError("Column 'count' is required for normalize_freq()")
        total = self.data["count"].sum()
        if total == 0:
            raise ValueError("Cannot normalize frequencies because total count is 0")
        df = self.data.copy()
        df["freq"] = df["count"] / total
        return Clonoset(
            df,
            chain=self.chain,
            sample_id=self.sample_id,
            metadata=self.metadata.copy(),
            standardize=False,
            by_umi=self.by_umi,
            validate=self.validate,
        )

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        return self.data[key]
