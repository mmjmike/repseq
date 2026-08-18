"""B-cell receptor sequence and somatic-hypermutation tree analysis."""

from __future__ import annotations

import os
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from . import plot as rsplot
from .io import read_clonoset


_SEQUENCE_COLUMNS = {
    "CDR1": "aaSeqCDR1",
    "FR2": "aaSeqFR2",
    "CDR2": "aaSeqCDR2",
    "FR3": "aaSeqFR3",
    "CDR3": "aaSeqCDR3",
    "FR4": "aaSeqFR4",
}


def _first_existing(columns, candidates):
    return next((column for column in candidates if column in columns), None)


def _id_key(value):
    if pd.isna(value):
        return None
    try:
        numeric = float(value)
        if numeric.is_integer():
            return str(int(numeric))
    except (TypeError, ValueError):
        pass
    return str(value)


def _observed_mask(data):
    if "isObserved" not in data.columns:
        raise ValueError("trees table must contain an 'isObserved' column")
    values = data["isObserved"]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype("string").str.lower().isin({"true", "t", "1", "yes"})


def _weighted_mode(values, weights):
    table = pd.DataFrame({"value": values, "weight": weights}).dropna(subset=["value"])
    if table.empty:
        return pd.NA
    table["value"] = table["value"].astype(str)
    table["weight"] = pd.to_numeric(table["weight"], errors="coerce").fillna(1)
    totals = table.groupby("value", sort=True)["weight"].sum()
    return totals[totals == totals.max()].index[0]


def _weighted_consensus(values, weights):
    sequences = pd.DataFrame({"sequence": values, "weight": weights}).dropna(subset=["sequence"])
    sequences["sequence"] = sequences["sequence"].astype(str)
    sequences = sequences[sequences["sequence"] != ""]
    if sequences.empty:
        return pd.NA
    sequences["weight"] = pd.to_numeric(sequences["weight"], errors="coerce").fillna(1)
    consensus = []
    for position in range(sequences["sequence"].str.len().max()):
        symbols = {}
        for sequence, weight in sequences.itertuples(index=False):
            if position < len(sequence):
                symbol = sequence[position]
                symbols[symbol] = symbols.get(symbol, 0) + weight
        if symbols:
            consensus.append(sorted(symbols, key=lambda symbol: (-symbols[symbol], symbol))[0])
    return "".join(consensus) if consensus else pd.NA


def _is_unswitched_isotype(isotype):
    normalized = "".join(character for character in str(isotype).upper() if character.isalnum())
    return normalized in {"M", "D", "IGM", "IGD", "IGHM", "IGHD"}


class TreeAnalyzer:
    """Analyze MiXCR SHM-tree node tables and their Newick trees."""

    def __init__(self):
        self.trees_df = None
        self.metadata = None
        self.newick_trees = {}
        self._trees_table_dir = None
        self._clonoset_cache = {}

    def _invalidate_properties(self):
        self.__dict__.pop("trees_properties", None)

    def read_trees_table(self, filename):
        """Read a MiXCR ``exportShmTreesWithNodes`` TSV table."""
        filename = Path(filename)
        self.trees_df = pd.read_csv(filename, sep="\t")
        self._trees_table_dir = filename.resolve().parent
        self._clonoset_cache.clear()
        self._add_sample_ids()
        self._attach_newick_trees()
        self._invalidate_properties()
        return self.trees_df

    def read_trees_newick(self, folder):
        """Read all ``<treeId>.tree`` Newick files from a folder."""
        folder = Path(folder)
        if not folder.is_dir():
            raise ValueError(f"Newick tree folder does not exist: {folder}")
        self.newick_trees = {
            path.stem: path.read_text().strip()
            for path in sorted(folder.iterdir())
            if path.is_file() and path.suffix == ".tree"
        }
        self._attach_newick_trees()
        return self.newick_trees

    def read_metadata(self, metadata):
        """Store sample metadata from a DataFrame or delimited text file."""
        if isinstance(metadata, pd.DataFrame):
            self.metadata = metadata.copy()
        else:
            self.metadata = pd.read_csv(metadata, sep=None, engine="python")
        self._add_sample_ids()
        return self.metadata

    def _attach_newick_trees(self):
        if self.trees_df is not None:
            self.trees_df.attrs["newick_trees"] = self.newick_trees

    def _add_sample_ids(self):
        if self.trees_df is not None and "fileName" in self.trees_df.columns:
            observed = _observed_mask(self.trees_df)
            self.trees_df["sample_id"] = pd.NA
            self.trees_df.loc[observed, "sample_id"] = self.trees_df.loc[
                observed, "fileName"
            ].map(self._sample_id)

    def _sample_id(self, filename):
        stem = Path(str(filename)).name
        if stem.endswith(".clns"):
            stem = stem[:-5]
        if self.metadata is not None and "sample_id" in self.metadata.columns:
            sample_ids = self.metadata["sample_id"].dropna().astype(str)
            matches = [sample_id for sample_id in sample_ids if stem == sample_id or stem.endswith(f".{sample_id}")]
            if matches:
                return max(matches, key=len)
        return stem.rsplit(".", 1)[-1]

    def _clonoset_filename(self, clns_filename):
        path = Path(str(clns_filename))
        if not path.is_absolute():
            path = (self._trees_table_dir or Path.cwd()) / path
        stem = path.name[:-5] if path.name.endswith(".clns") else path.stem
        candidates = [
            path.with_name(f"{stem}.tsv"),
            path.with_name(f"{stem}.clones.tsv"),
        ]
        return next((candidate for candidate in candidates if candidate.is_file()), None)

    def _read_clonoset_umi(self, clns_filename, clone_id):
        filename = self._clonoset_filename(clns_filename)
        if filename is None:
            return np.nan
        cache_key = os.fspath(filename)
        if cache_key not in self._clonoset_cache:
            clonoset = read_clonoset(filename)
            if "cloneId" not in clonoset.columns or "uniqueMoleculeCount" not in clonoset.columns:
                self._clonoset_cache[cache_key] = {}
            else:
                clone_ids = clonoset["cloneId"].map(_id_key)
                umi = pd.to_numeric(clonoset["uniqueMoleculeCount"], errors="coerce")
                self._clonoset_cache[cache_key] = dict(zip(clone_ids, umi))
        return self._clonoset_cache[cache_key].get(_id_key(clone_id), np.nan)

    def _add_umi_values(self, data):
        if "uniqueMoleculeCount" not in data.columns:
            data["uniqueMoleculeCount"] = np.nan
        else:
            data["uniqueMoleculeCount"] = pd.to_numeric(
                data["uniqueMoleculeCount"], errors="coerce"
            )
        required = {"fileName", "cloneId"}
        if required.issubset(data.columns):
            missing = _observed_mask(data) & data["uniqueMoleculeCount"].isna()
            for index, row in data.loc[missing, ["fileName", "cloneId"]].iterrows():
                data.at[index, "uniqueMoleculeCount"] = self._read_clonoset_umi(
                    row["fileName"], row["cloneId"]
                )
        return data

    def _ensure_umi_values(self):
        enriched = self._add_umi_values(self.trees_df.copy())
        self.trees_df["uniqueMoleculeCount"] = enriched["uniqueMoleculeCount"]
        return enriched

    @cached_property
    def trees_properties(self):
        """Return one cached summary row per tree."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before calculating tree properties")
        if "treeId" not in self.trees_df.columns:
            raise ValueError("trees table must contain a 'treeId' column")

        data = self._ensure_umi_values()
        mutation_column = _first_existing(
            data.columns,
            ["mutationRate", "nMutationRate", "nMutationsRate", "mutation_rate"],
        )
        distance_column = _first_existing(
            data.columns,
            ["distanceFromGermline", "distance_from_germline", "distanceFromRoot"],
        )
        rows = []
        for tree_id, tree in data.groupby("treeId", sort=False, dropna=False):
            observed = tree.loc[_observed_mask(tree)].copy()
            read_values = (
                pd.to_numeric(observed["readCount"], errors="coerce").fillna(1)
                if "readCount" in observed.columns
                else pd.Series(1, index=observed.index, dtype=float)
            )
            weights = observed["uniqueMoleculeCount"].where(
                observed["uniqueMoleculeCount"].notna(),
                read_values,
            )
            isotypes = []
            if "isotype" in observed.columns:
                isotypes = sorted(observed["isotype"].dropna().astype(str).unique().tolist())
            umi_values = pd.to_numeric(observed["uniqueMoleculeCount"], errors="coerce")
            umi = umi_values.sum(min_count=1)
            reads = (
                pd.to_numeric(observed["readCount"], errors="coerce").sum()
                if "readCount" in observed.columns else 0
            )
            row = {
                "treeId": tree_id,
                "nodes": int(len(tree)),
                "nodes_obs": int(len(observed)),
                "reads": int(reads),
                "umi": pd.NA if pd.isna(umi) else int(umi),
                "isotypes": isotypes,
                "v": _weighted_mode(observed.get("bestVHit", pd.Series(index=observed.index)), weights),
                "j": _weighted_mode(observed.get("bestJHit", pd.Series(index=observed.index)), weights),
            }
            for name, column in _SEQUENCE_COLUMNS.items():
                values = observed.get(column, pd.Series(index=observed.index, dtype="object"))
                row[f"consensus_{name}"] = _weighted_consensus(values, weights)
            if mutation_column is None:
                row["mean_mutation_rate"] = np.nan
            else:
                row["mean_mutation_rate"] = pd.to_numeric(
                    observed[mutation_column], errors="coerce"
                ).mean()
            row["isotype_switched"] = any(
                not _is_unswitched_isotype(isotype) for isotype in isotypes
            )
            if distance_column is None:
                row["max_distance_from_germline"] = np.nan
            else:
                row["max_distance_from_germline"] = pd.to_numeric(
                    tree[distance_column], errors="coerce"
                ).max()
            rows.append(row)

        properties = pd.DataFrame(rows)
        properties["umi"] = properties["umi"].astype("Int64")
        sort_count = properties["umi"].astype("Float64").fillna(properties["reads"])
        properties = properties.assign(_sort_count=sort_count).sort_values(
            ["nodes_obs", "_sort_count", "consensus_CDR3"],
            ascending=[False, False, True],
            na_position="last",
            kind="stable",
        )
        return properties.drop(columns="_sort_count").reset_index(drop=True)

    def draw_tree(self, treeId, group="isotype", label="timepoint"):
        """Draw one loaded tree, optionally colored and labeled by metadata."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before drawing a tree")
        self._ensure_umi_values()
        return rsplot.draw_tree(
            self.trees_df,
            treeId,
            metadata=self.metadata,
            group=group,
            label=label,
        )


__all__ = ["TreeAnalyzer"]
