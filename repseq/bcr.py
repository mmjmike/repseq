"""B-cell receptor sequence and somatic-hypermutation tree analysis."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import plot as rsplot
from . import clonosets
from .clone_filter import Filter
from .common_functions import print_progress_bar


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


class _TreePropertiesTable(pd.DataFrame):
    _metadata = ["_analyzer"]

    @property
    def _constructor(self):
        return pd.DataFrame

    def __call__(self, all_consensuses=False, verbose=True):
        return self._analyzer._get_trees_properties(
            all_consensuses=all_consensuses,
            verbose=verbose,
        )


class TreeAnalyzer:
    """Analyze MiXCR SHM-tree node tables and their Newick trees."""

    def __init__(self):
        self.trees_df = None
        self.metadata = None
        self.newick_trees = {}
        self._trees_table_dir = None
        self._umi_enriched = False
        self._trees_properties_cache = None
        self._all_consensuses_cached = False

    def _invalidate_properties(self):
        self._trees_properties_cache = None
        self._all_consensuses_cached = False

    def read_trees_table(self, filename):
        """Read a MiXCR ``exportShmTreesWithNodes`` TSV table."""
        filename = Path(filename)
        self.trees_df = pd.read_csv(filename, sep="\t")
        self._trees_table_dir = filename.resolve().parent
        self._umi_enriched = False
        self._add_sample_ids()
        self._attach_newick_trees()
        self._ensure_umi_values()
        self._invalidate_properties()
        observed = _observed_mask(self.trees_df)
        tree_count = self.trees_df["treeId"].nunique(dropna=True)
        sample_count = (
            self.trees_df.loc[observed, "sample_id"].nunique(dropna=True)
            if "sample_id" in self.trees_df.columns else 0
        )
        print(
            f"Read {tree_count} trees from {sample_count} samples: "
            f"{len(self.trees_df)} nodes, {int(observed.sum())} observed nodes"
        )

    def read_trees_newick(self, folder):
        """Attach ``<treeId>.tree`` filenames without reading their contents."""
        folder = Path(folder).resolve()
        if not folder.is_dir():
            raise ValueError(f"Newick tree folder does not exist: {folder}")
        self.newick_trees = {
            _id_key(path.stem): os.fspath(path)
            for path in sorted(folder.iterdir())
            if path.is_file() and path.suffix == ".tree"
        }
        self._attach_newick_trees()
        print(f"Read {len(self.newick_trees)} Newick tree filenames")

    def read_metadata(self, metadata):
        """Left-join sample metadata into the loaded trees table."""
        if self.trees_df is None:
            print("Trees table has not been read. Metadata was not loaded")
            return
        if isinstance(metadata, pd.DataFrame):
            metadata_df = metadata.copy()
        else:
            metadata_df = pd.read_csv(metadata, sep=None, engine="python")
        if "sample_id" not in metadata_df.columns:
            raise ValueError("metadata must contain a 'sample_id' column")
        if metadata_df["sample_id"].duplicated().any():
            raise ValueError("metadata must contain one row per sample_id")

        self.metadata = metadata_df
        self._add_sample_ids()
        self.trees_df["sample_id"] = self.trees_df["sample_id"].astype("string")
        metadata_df = metadata_df.copy()
        metadata_df["sample_id"] = metadata_df["sample_id"].astype("string")
        observed = _observed_mask(self.trees_df)
        tree_samples = set(self.trees_df.loc[observed, "sample_id"].dropna())
        metadata_samples = set(metadata_df["sample_id"].dropna())
        matched_samples = tree_samples & metadata_samples
        updated_nodes = int(
            (observed & self.trees_df["sample_id"].isin(matched_samples)).sum()
        )

        metadata_columns = [column for column in metadata_df.columns if column != "sample_id"]
        overlapping_columns = [
            column for column in metadata_columns if column in self.trees_df.columns
        ]
        self.trees_df = self.trees_df.merge(
            metadata_df,
            on="sample_id",
            how="left",
            suffixes=("", "_metadata"),
            sort=False,
        )
        for column in overlapping_columns:
            metadata_column = f"{column}_metadata"
            self.trees_df[column] = self.trees_df[metadata_column].combine_first(
                self.trees_df[column]
            )
            self.trees_df = self.trees_df.drop(columns=metadata_column)
        self._invalidate_properties()
        print(
            f"Metadata updated for {len(matched_samples)} samples and "
            f"{updated_nodes} observed nodes"
        )

    def _attach_newick_trees(self):
        if self.trees_df is not None:
            self.trees_df["newick_filename"] = self.trees_df["treeId"].map(
                lambda tree_id: self.newick_trees.get(_id_key(tree_id), pd.NA)
            )

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

    def _resolve_clns_filename(self, clns_filename):
        raw_path = Path(str(clns_filename))
        if raw_path.is_absolute():
            return raw_path
        candidates = [Path.cwd() / raw_path]
        if self._trees_table_dir is not None:
            candidates.append(self._trees_table_dir / raw_path)
        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()
        if self._trees_table_dir is not None:
            matches = sorted(
                {
                    match.resolve()
                    for root in [self._trees_table_dir, self._trees_table_dir.parent]
                    for match in root.rglob(raw_path.name)
                    if match.is_file()
                },
                key=lambda match: (len(match.parts), os.fspath(match)),
            )
            if matches:
                return matches[0]
        return None

    def _source_folders_and_filenames(self, data):
        filenames = data.loc[_observed_mask(data), "fileName"].dropna().astype(str).unique()
        resolved = {}
        for filename in filenames:
            resolved_filename = self._resolve_clns_filename(filename)
            if resolved_filename is not None:
                resolved[filename] = resolved_filename
        folders = sorted({os.fspath(filename.parent) for filename in resolved.values()})
        return folders, resolved

    @staticmethod
    def _matching_sample_id(clns_filename, sample_ids):
        stem = Path(str(clns_filename)).name
        if stem.endswith(".clns"):
            stem = stem[:-5]
        matches = [
            sample_id for sample_id in sample_ids
            if stem == sample_id or stem.endswith(f".{sample_id}")
        ]
        return max(matches, key=len) if matches else None

    def _add_umi_values(self, data):
        if "uniqueMoleculeCount" not in data.columns:
            data["uniqueMoleculeCount"] = np.nan
        else:
            data["uniqueMoleculeCount"] = pd.to_numeric(
                data["uniqueMoleculeCount"], errors="coerce"
            )
        required = {"fileName", "cloneId"}
        missing = _observed_mask(data) & data["uniqueMoleculeCount"].isna()
        if not required.issubset(data.columns) or not missing.any():
            return data

        folders, resolved_filenames = self._source_folders_and_filenames(data)
        if not folders:
            return data
        clonosets_df = clonosets.find_all_mixcr_clonosets(folders)
        if clonosets_df.empty:
            return data

        discovered_sample_ids = clonosets_df["sample_id"].dropna().astype(str).unique()
        filename_to_sample = {
            filename: self._matching_sample_id(resolved, discovered_sample_ids)
            for filename, resolved in resolved_filenames.items()
        }
        filename_to_sample = {
            filename: sample_id
            for filename, sample_id in filename_to_sample.items()
            if sample_id is not None
        }
        if not filename_to_sample:
            return data

        observed = _observed_mask(data)
        discovered_ids = data.loc[observed, "fileName"].astype(str).map(filename_to_sample)
        discovered_ids = discovered_ids.dropna()
        data.loc[discovered_ids.index, "sample_id"] = discovered_ids
        source_sample_ids = set(filename_to_sample.values())
        clonosets_df = clonosets_df.loc[
            clonosets_df["sample_id"].astype(str).isin(source_sample_ids)
        ].reset_index(drop=True)
        if clonosets_df.empty:
            return data

        pooled = clonosets.pool_clonotypes_from_clonosets_df(
            clonosets_df,
            cl_filter=Filter(convert=False),
        )
        umi_columns = {"sample_id", "cloneId", "uniqueMoleculeCount"}
        if not umi_columns.issubset(pooled.columns):
            return data
        umi_df = pooled.loc[:, ["sample_id", "cloneId", "uniqueMoleculeCount"]].copy()
        umi_df["sample_id"] = umi_df["sample_id"].astype(str)
        umi_df["_clone_id"] = umi_df["cloneId"].map(_id_key)
        umi_df["uniqueMoleculeCount"] = pd.to_numeric(
            umi_df["uniqueMoleculeCount"], errors="coerce"
        )
        umi_df = umi_df.drop(columns="cloneId").drop_duplicates(
            ["sample_id", "_clone_id"], keep="first"
        )

        data["_clone_id"] = data["cloneId"].map(_id_key)
        data = data.merge(
            umi_df,
            on=["sample_id", "_clone_id"],
            how="left",
            suffixes=("", "_pooled"),
            sort=False,
        )
        data["uniqueMoleculeCount"] = data["uniqueMoleculeCount"].fillna(
            data.pop("uniqueMoleculeCount_pooled")
        )
        data = data.drop(columns="_clone_id")
        return data

    def _ensure_umi_values(self):
        if self._umi_enriched:
            return self.trees_df.copy()
        enriched = self._add_umi_values(self.trees_df.copy())
        self.trees_df = enriched
        self._attach_newick_trees()
        self._umi_enriched = True
        return enriched

    @property
    def trees_properties(self):
        """Return the cached tree summary; call it to request all consensuses."""
        return self._get_trees_properties(all_consensuses=False, verbose=True)

    def _get_trees_properties(self, all_consensuses=False, verbose=True):
        if self.trees_df is None:
            raise ValueError("Read a trees table before calculating tree properties")
        if "treeId" not in self.trees_df.columns:
            raise ValueError("trees table must contain a 'treeId' column")

        if self._trees_properties_cache is None:
            self._calculate_base_tree_properties(verbose=verbose)
        if all_consensuses and not self._all_consensuses_cached:
            self._add_consensus_columns(verbose=verbose)
        return self._trees_properties_cache

    def _calculate_base_tree_properties(self, verbose=True):
        data = self._ensure_umi_values()
        grouped_trees = list(data.groupby("treeId", sort=False, dropna=False))
        tree_total = len(grouped_trees)
        if verbose:
            print(f"Calculating properties for {tree_total} trees (CDR3 consensus only)")
            print_progress_bar(0, tree_total, "Tree properties", object_name="tree(s)")

        mutation_column = _first_existing(
            data.columns,
            ["mutationRate", "nMutationRate", "nMutationsRate", "mutation_rate"],
        )
        distance_column = _first_existing(
            data.columns,
            [
                "DistanceFromGermline",
                "distanceFromGermline",
                "distance_from_germline",
                "distanceFromRoot",
            ],
        )
        rows = []
        for tree_index, (tree_id, tree) in enumerate(grouped_trees, start=1):
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
                "v": _weighted_mode(
                    observed.get(
                        "bestVHit",
                        pd.Series(index=observed.index, dtype="object"),
                    ),
                    weights,
                ),
                "j": _weighted_mode(
                    observed.get(
                        "bestJHit",
                        pd.Series(index=observed.index, dtype="object"),
                    ),
                    weights,
                ),
            }
            cdr3_values = observed.get(
                _SEQUENCE_COLUMNS["CDR3"],
                pd.Series(index=observed.index, dtype="object"),
            )
            row["consensus_CDR3"] = _weighted_consensus(cdr3_values, weights)
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
            if verbose:
                print_progress_bar(
                    tree_index,
                    tree_total,
                    "Tree properties",
                    object_name="tree(s)",
                )

        properties = pd.DataFrame(rows)
        properties["umi"] = properties["umi"].astype("Int64")
        sort_count = properties["umi"].astype("Float64").fillna(properties["reads"])
        properties = properties.assign(_sort_count=sort_count).sort_values(
            ["nodes_obs", "_sort_count", "consensus_CDR3"],
            ascending=[False, False, True],
            na_position="last",
            kind="stable",
        )
        properties = properties.drop(columns="_sort_count").reset_index(drop=True)
        cached = _TreePropertiesTable(properties)
        cached._analyzer = self
        self._trees_properties_cache = cached
        if verbose:
            print("Finished calculating tree properties")

    def _add_consensus_columns(self, verbose=True):
        data = self._ensure_umi_values()
        grouped_trees = list(data.groupby("treeId", sort=False, dropna=False))
        tree_total = len(grouped_trees)
        extra_sequences = {
            name: column
            for name, column in _SEQUENCE_COLUMNS.items()
            if name != "CDR3"
        }
        consensus_by_name = {name: {} for name in extra_sequences}
        if verbose:
            print(f"Calculating additional consensus sequences for {tree_total} trees")
            print_progress_bar(0, tree_total, "Tree consensuses", object_name="tree(s)")

        for tree_index, (tree_id, tree) in enumerate(grouped_trees, start=1):
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
            for name, column in extra_sequences.items():
                values = observed.get(
                    column,
                    pd.Series(index=observed.index, dtype="object"),
                )
                consensus_by_name[name][_id_key(tree_id)] = _weighted_consensus(
                    values,
                    weights,
                )
            if verbose:
                print_progress_bar(
                    tree_index,
                    tree_total,
                    "Tree consensuses",
                    object_name="tree(s)",
                )

        cdr3_position = self._trees_properties_cache.columns.get_loc("consensus_CDR3")
        insert_position = cdr3_position
        for name in ["CDR1", "FR2", "CDR2", "FR3"]:
            values = self._trees_properties_cache["treeId"].map(
                lambda tree_id: consensus_by_name[name].get(_id_key(tree_id), pd.NA)
            )
            self._trees_properties_cache.insert(
                insert_position,
                f"consensus_{name}",
                values,
            )
            insert_position += 1
        fr4_values = self._trees_properties_cache["treeId"].map(
            lambda tree_id: consensus_by_name["FR4"].get(_id_key(tree_id), pd.NA)
        )
        cdr3_position = self._trees_properties_cache.columns.get_loc("consensus_CDR3")
        self._trees_properties_cache.insert(
            cdr3_position + 1,
            "consensus_FR4",
            fr4_values,
        )
        self._all_consensuses_cached = True
        if verbose:
            print("Finished calculating additional consensus sequences")

    def draw_tree(self, treeId, group="isotype", label="timepoint"):
        """Draw one loaded tree, optionally colored and labeled by metadata."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before drawing a tree")
        self._ensure_umi_values()
        return rsplot.draw_tree(
            self.trees_df,
            treeId,
            metadata=None,
            group=group,
            label=label,
        )

    def get_logo_for_tree(self, treeId):
        """Plot an amino-acid CDR3 sequence logo for one tree."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before plotting a tree logo")
        if "treeId" not in self.trees_df.columns:
            raise ValueError("trees table must contain a 'treeId' column")
        if "aaSeqCDR3" not in self.trees_df.columns:
            raise ValueError("trees table must contain an 'aaSeqCDR3' column")

        tree_id_key = _id_key(treeId)
        tree_rows = self.trees_df.loc[
            self.trees_df["treeId"].map(_id_key) == tree_id_key
        ]
        if tree_rows.empty:
            raise ValueError(f"treeId {treeId!r} is not present in trees_df")
        sequences = tree_rows["aaSeqCDR3"].dropna().astype(str)
        sequences = sequences.loc[sequences != ""]
        if sequences.empty:
            raise ValueError(f"treeId {treeId!r} does not contain aaSeqCDR3 sequences")

        try:
            from . import logo
        except ModuleNotFoundError as error:
            if error.name == "logomaker":
                raise ImportError(
                    "get_logo_for_tree requires logomaker. "
                    "Install repseq with the clustering optional dependencies."
                ) from error
            raise
        list_of_clonotypes = [(sequence,) for sequence in sequences]
        return logo.get_logo_for_list_of_clonotypes(list_of_clonotypes, "prot")

    def to_count_table(self):
        """Create a wide table of tree abundance by sample."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before creating a count table")
        data = self._ensure_umi_values()
        required_columns = {"treeId", "sample_id", "isObserved"}
        missing_columns = sorted(required_columns - set(data.columns))
        if missing_columns:
            raise ValueError(
                f"trees table is missing required columns: {', '.join(missing_columns)}"
            )

        observed = data.loc[_observed_mask(data)].copy()
        observed = observed.loc[observed["sample_id"].notna()]
        read_counts = (
            pd.to_numeric(observed["readCount"], errors="coerce")
            if "readCount" in observed.columns
            else pd.Series(0, index=observed.index, dtype=float)
        )
        if "uniqueMoleculeCount" in observed.columns:
            observed["_count"] = pd.to_numeric(
                observed["uniqueMoleculeCount"], errors="coerce"
            ).fillna(read_counts)
        else:
            observed["_count"] = read_counts
        observed["_count"] = observed["_count"].fillna(0)

        properties = self.trees_properties(all_consensuses=False).loc[
            :, ["treeId", "v", "j", "consensus_CDR3"]
        ].rename(columns={"consensus_CDR3": "consensus_cdr3"})
        found_samples = observed["sample_id"].dropna().astype(str).unique().tolist()
        if self.metadata is not None and "sample_id" in self.metadata.columns:
            metadata_order = self.metadata["sample_id"].dropna().astype(str).tolist()
            sample_ids = [sample_id for sample_id in metadata_order if sample_id in found_samples]
            sample_ids.extend(sorted(set(found_samples) - set(sample_ids)))
        else:
            sample_ids = sorted(found_samples)

        abundance = observed.assign(sample_id=observed["sample_id"].astype(str)).groupby(
            ["treeId", "sample_id"], sort=False, dropna=False
        )["_count"].sum()
        count_table = properties.copy()
        for sample_id in sample_ids:
            sample_counts = abundance.xs(sample_id, level="sample_id", drop_level=True)
            values = pd.to_numeric(
                count_table["treeId"].map(sample_counts),
                errors="coerce",
            ).fillna(0)
            numeric_values = values.astype(float).to_numpy()
            if np.all(np.isclose(numeric_values, np.round(numeric_values))):
                values = pd.Series(
                    np.round(numeric_values),
                    index=values.index,
                    dtype="Int64",
                )
            count_table[sample_id] = values
        return count_table


__all__ = ["TreeAnalyzer"]
