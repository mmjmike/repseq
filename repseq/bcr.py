"""B-cell receptor sequence and somatic-hypermutation tree analysis."""

from __future__ import annotations

import os
import re
from collections import Counter
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

_MUTATION_REGIONS = [
    ("CDR1", 5, 6),
    ("FR2", 6, 7),
    ("CDR2", 7, 8),
    ("FR3", 8, 9),
    ("CDR3", 9, 18),
    ("FR4", 18, 19),
]


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


def _format_axis_value(value):
    if pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)) and np.isfinite(value) and value.is_integer():
        return str(int(value))
    return str(value)


def get_mutation_positions(mutations):
    """Parse MiXCR substitution, deletion, and insertion positions."""
    if not isinstance(mutations, str) or mutations == "":
        return [], [[], [], []]
    mutation_tokens = mutations.replace("D", "S").replace("I", "S").split("S")[1:]
    substitution_motif = re.compile(r"^[ATGC](\d+)[ATGC]$")
    deletion_motif = re.compile(r"^[ATGC](\d+)$")
    insertion_motif = re.compile(r"^(\d+)[ATGC]$")

    substitutions = []
    deletions = []
    insertions = []
    for mutation in mutation_tokens:
        substitution = substitution_motif.fullmatch(mutation)
        deletion = deletion_motif.fullmatch(mutation)
        insertion = insertion_motif.fullmatch(mutation)
        if substitution is not None:
            substitutions.append(int(substitution.group(1)))
        elif deletion is not None:
            deletions.append(int(deletion.group(1)))
        elif insertion is not None:
            insertions.append(int(insertion.group(1)))
    positions = substitutions + deletions + insertions
    return positions, [substitutions, deletions, insertions]


def _absolute_mutation_positions(row, segment):
    alignment = row.get(f"all{segment.upper()}Alignments")
    if not isinstance(alignment, str) or alignment == "":
        return []
    alignment_fields = alignment.split(";", 1)[0].split("|")
    if len(alignment_fields) < 6:
        return []
    try:
        target_from = int(alignment_fields[0])
        query_from = int(alignment_fields[3])
    except (TypeError, ValueError):
        return []
    mutation_positions = get_mutation_positions(alignment_fields[5])[0]
    segment_shift = query_from - target_from
    return [position + segment_shift for position in mutation_positions]


def _parse_refpoint_regions(ref_points):
    if not isinstance(ref_points, str) or ref_points == "":
        return None
    required_indices = sorted(
        {index for _, start, end in _MUTATION_REGIONS for index in (start, end)}
    )
    candidates = []
    for target_ref_points in ref_points.split(","):
        values = target_ref_points.split(":")
        parsed = []
        for value in values:
            try:
                parsed.append(int(value) if value != "" else None)
            except ValueError:
                parsed.append(None)
        coverage = sum(
            index < len(parsed) and parsed[index] is not None
            for index in required_indices
        )
        candidates.append((coverage, parsed))
    if not candidates:
        return None
    _, points = max(candidates, key=lambda candidate: candidate[0])
    if any(index >= len(points) or points[index] is None for index in required_indices):
        return None

    cdr1_begin = points[5]
    regions = []
    for region, start_index, end_index in _MUTATION_REGIONS:
        start = points[start_index] - cdr1_begin
        end = points[end_index] - cdr1_begin
        if start < 0 or end <= start:
            return None
        regions.append((region, start, end))
    return cdr1_begin, regions


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
        for column in ["uniqueMoleculeCount", "uniqueMoleculeFraction"]:
            if column not in data.columns:
                data[column] = np.nan
            else:
                data[column] = pd.to_numeric(data[column], errors="coerce")
        required = {"fileName", "cloneId"}
        missing = _observed_mask(data) & (
            data["uniqueMoleculeCount"].isna()
            | data["uniqueMoleculeFraction"].isna()
        )
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
        required_umi_columns = {"sample_id", "cloneId", "uniqueMoleculeCount"}
        if not required_umi_columns.issubset(pooled.columns):
            return data
        value_columns = [
            column
            for column in ["uniqueMoleculeCount", "uniqueMoleculeFraction"]
            if column in pooled.columns
        ]
        umi_df = pooled.loc[:, ["sample_id", "cloneId", *value_columns]].copy()
        umi_df["sample_id"] = umi_df["sample_id"].astype(str)
        umi_df["_clone_id"] = umi_df["cloneId"].map(_id_key)
        for column in value_columns:
            umi_df[column] = pd.to_numeric(umi_df[column], errors="coerce")
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
        for column in value_columns:
            pooled_column = f"{column}_pooled"
            data[column] = data[column].fillna(data.pop(pooled_column))
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

    def get_tree_clonotypes(self, treeId):
        """Return full clonoset rows for all observed clonotypes in one tree."""
        if self.trees_df is None:
            raise ValueError("Read a trees table before retrieving tree clonotypes")
        required_columns = {"treeId", "isObserved", "fileName", "cloneId", "sample_id"}
        missing_columns = sorted(required_columns - set(self.trees_df.columns))
        if missing_columns:
            raise ValueError(
                f"trees table is missing required columns: {', '.join(missing_columns)}"
            )

        tree_id_key = _id_key(treeId)
        tree_rows = self.trees_df.loc[
            (self.trees_df["treeId"].map(_id_key) == tree_id_key)
            & _observed_mask(self.trees_df)
        ].copy()
        if tree_rows.empty:
            raise ValueError(f"treeId {treeId!r} does not contain observed clonotypes")

        folders, resolved_filenames = self._source_folders_and_filenames(tree_rows)
        if not folders:
            raise ValueError(f"Unable to locate source clonoset folders for treeId {treeId!r}")
        clonosets_df = clonosets.find_all_mixcr_clonosets(folders)
        if clonosets_df.empty:
            raise ValueError(f"No exported clonoset tables found for treeId {treeId!r}")

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
        if filename_to_sample:
            discovered_ids = tree_rows["fileName"].astype(str).map(filename_to_sample)
            tree_rows.loc[discovered_ids.notna(), "sample_id"] = discovered_ids.dropna()

        sample_ids = tree_rows["sample_id"].dropna().astype(str).unique()
        clonosets_df = clonosets_df.loc[
            clonosets_df["sample_id"].astype(str).isin(sample_ids)
        ].reset_index(drop=True)
        if clonosets_df.empty:
            raise ValueError(f"No source clonosets matched treeId {treeId!r} samples")

        pooled = clonosets.pool_clonotypes_from_clonosets_df(
            clonosets_df,
            cl_filter=Filter(convert=False),
        )
        if not {"sample_id", "cloneId"}.issubset(pooled.columns):
            raise ValueError("Source clonosets must contain sample_id and cloneId columns")

        mutation_column = _first_existing(
            tree_rows.columns,
            ["nMutationsRate", "nMutationRate", "mutationRate", "mutation_rate"],
        )
        tree_rows["nMutationsRate"] = (
            pd.to_numeric(tree_rows[mutation_column], errors="coerce")
            if mutation_column is not None else np.nan
        )
        if "isotype" not in tree_rows.columns:
            tree_rows["isotype"] = pd.NA

        metadata_columns = []
        if self.metadata is not None:
            metadata_columns = [
                column
                for column in self.metadata.columns
                if column != "sample_id" and column in tree_rows.columns
            ]
        tree_columns = [
            "treeId",
            "sample_id",
            "cloneId",
            "nMutationsRate",
            *metadata_columns,
            "isotype",
        ]
        tree_info = tree_rows.loc[:, list(dict.fromkeys(tree_columns))].copy()
        tree_info["sample_id"] = tree_info["sample_id"].astype(str)
        tree_info["_clone_id"] = tree_info["cloneId"].map(_id_key)
        pooled = pooled.copy()
        pooled["sample_id"] = pooled["sample_id"].astype(str)
        pooled["_clone_id"] = pooled["cloneId"].map(_id_key)

        merged = tree_info.merge(
            pooled,
            on=["sample_id", "_clone_id"],
            how="left",
            suffixes=("_tree", ""),
            sort=False,
        )
        if "cloneId" not in merged.columns or merged["cloneId"].isna().any():
            missing = merged.loc[
                merged.get("cloneId", pd.Series(index=merged.index, dtype="object")).isna(),
                ["sample_id", "cloneId_tree"],
            ]
            raise ValueError(
                "Unable to find source clonotypes for: "
                + ", ".join(
                    f"{row.sample_id}/cloneId={row.cloneId_tree}"
                    for row in missing.itertuples(index=False)
                )
            )

        output = pd.DataFrame(index=merged.index)
        output["treeId"] = merged[
            "treeId_tree" if "treeId_tree" in merged.columns else "treeId"
        ]
        output["sample_id"] = merged["sample_id"]
        output["cloneId"] = merged["cloneId"]
        requested_tree_columns = ["nMutationsRate", *metadata_columns, "isotype"]
        for column in requested_tree_columns:
            tree_column = f"{column}_tree" if f"{column}_tree" in merged.columns else column
            output[column] = merged[tree_column]

        reserved_columns = set(output.columns) | {"_clone_id"}
        for column in pooled.columns:
            if column not in reserved_columns:
                output[column] = merged[column]
        return output.reset_index(drop=True)

    def _get_mutation_rate_df(self, treeId):
        if self.trees_df is None:
            raise ValueError("Read a trees table before plotting mutation rates")
        required_columns = {
            "treeId",
            "isObserved",
            "allVAlignments",
            "allDAlignments",
            "allJAlignments",
            "refPoints",
        }
        missing_columns = sorted(required_columns - set(self.trees_df.columns))
        if missing_columns:
            raise ValueError(
                "trees table is missing required mutation columns: "
                + ", ".join(missing_columns)
            )

        tree_id_key = _id_key(treeId)
        tree_rows = self.trees_df.loc[
            (self.trees_df["treeId"].map(_id_key) == tree_id_key)
            & _observed_mask(self.trees_df)
        ]
        if tree_rows.empty:
            raise ValueError(f"treeId {treeId!r} does not contain observed nodes")

        parsed_rows = []
        for _, row in tree_rows.iterrows():
            parsed_ref_points = _parse_refpoint_regions(row["refPoints"])
            if parsed_ref_points is None:
                continue
            cdr1_begin, regions = parsed_ref_points
            region_lengths = tuple(end - start for _, start, end in regions)
            mutation_positions = set()
            for segment in ["v", "d", "j"]:
                mutation_positions.update(_absolute_mutation_positions(row, segment))
            parsed_rows.append((cdr1_begin, region_lengths, mutation_positions))

        if not parsed_rows:
            raise ValueError(
                f"treeId {treeId!r} does not contain complete CDR1-to-FR4 refPoints"
            )

        canonical_lengths = Counter(
            region_lengths for _, region_lengths, _ in parsed_rows
        ).most_common(1)[0][0]
        regions = []
        region_start = 0
        for (region, _, _), length in zip(_MUTATION_REGIONS, canonical_lengths):
            regions.append((region, region_start, region_start + length))
            region_start += length
        total_length = region_start

        mutation_counts = np.zeros(total_length, dtype=float)
        for cdr1_begin, _, mutation_positions in parsed_rows:
            relative_positions = {
                position - cdr1_begin
                for position in mutation_positions
                if 0 <= position - cdr1_begin < total_length
            }
            for position in relative_positions:
                mutation_counts[position] += 1

        rates = mutation_counts / len(parsed_rows)
        mutation_rate_df = pd.DataFrame(
            {
                "position": np.arange(total_length, dtype=int),
                "rate": rates,
            }
        )
        mutation_rate_df["region"] = ""
        for region, start, end in regions:
            mutation_rate_df.loc[
                mutation_rate_df["position"].between(start, end - 1),
                "region",
            ] = region
        return mutation_rate_df

    def plot_mutations_rate(self, treeId, ax=None):
        """Plot observed-node mutation frequencies across CDR1 through FR4."""
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.patches import Patch

        mutation_rate_df = self._get_mutation_rate_df(treeId)
        region_order = [region for region, _, _ in _MUTATION_REGIONS]
        palette = dict(
            zip(region_order, sns.color_palette("Set2", n_colors=len(region_order)))
        )
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 4))

        colors = mutation_rate_df["region"].map(palette)
        ax.bar(
            mutation_rate_df["position"],
            mutation_rate_df["rate"],
            width=0.7,
            color=colors,
            edgecolor="none",
        )
        region_ends = (
            mutation_rate_df.groupby("region", sort=False)["position"].max().tolist()
        )
        for border in region_ends[:-1]:
            ax.axvline(border + 0.5, color="0.45", linestyle="--", linewidth=1)

        handles = [
            Patch(facecolor=palette[region], edgecolor="none", label=region)
            for region in region_order
            if region in set(mutation_rate_df["region"])
        ]
        ax.legend(
            handles=handles,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.18),
            ncol=len(handles),
        )
        ax.set_xlabel("Position")
        ax.set_ylabel("Mutation frequency")
        ax.set_title(f"Tree {treeId} mutation frequencies")
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.7, mutation_rate_df["position"].max() + 0.7)
        return ax

    def _get_timepoint_trajectory_df(
        self,
        treeId,
        timepoint_feature=None,
        by_freq=True,
    ):
        feature = "timepoint" if timepoint_feature is None else timepoint_feature
        if self.trees_df is None or self.metadata is None:
            print(
                "Cannot plot timepoint trajectory: first read the trees table and "
                "run ta.read_metadata(metadata)"
            )
            return None
        if feature not in self.trees_df.columns:
            print(
                f"Cannot plot timepoint trajectory: trees_df does not contain "
                f"'{feature}'. Run ta.read_metadata(metadata) with this column, "
                "or specify timepoint_feature='column_name'"
            )
            return None

        value_column = "uniqueMoleculeFraction" if by_freq else "uniqueMoleculeCount"
        data = self._ensure_umi_values()
        if value_column not in data.columns:
            print(
                f"Cannot plot timepoint trajectory: trees_df does not contain "
                f"'{value_column}'"
            )
            return None
        tree_id_key = _id_key(treeId)
        tree_rows = data.loc[
            (data["treeId"].map(_id_key) == tree_id_key)
            & _observed_mask(data)
        ].copy()
        if tree_rows.empty:
            raise ValueError(f"treeId {treeId!r} does not contain observed nodes")

        tree_rows[value_column] = pd.to_numeric(
            tree_rows[value_column], errors="coerce"
        )
        sample_values = (
            tree_rows.dropna(subset=["sample_id", feature])
            .groupby(["sample_id", feature], sort=False, dropna=False)[value_column]
            .sum(min_count=1)
            .dropna()
            .reset_index(name="value")
        )
        all_timepoints = data.loc[
            _observed_mask(data), feature
        ].dropna()
        if all_timepoints.empty:
            print(
                f"Cannot plot timepoint trajectory: trees_df has no '{feature}' values"
            )
            return None

        timepoint_values = all_timepoints.drop_duplicates().tolist()
        if isinstance(all_timepoints.dtype, pd.CategoricalDtype):
            present = set(timepoint_values)
            timepoint_order = [
                value for value in all_timepoints.cat.categories
                if value in present
            ]
        else:
            try:
                timepoint_order = sorted(timepoint_values)
            except TypeError:
                timepoint_order = sorted(timepoint_values, key=lambda value: str(value))

        if sample_values.empty:
            trajectory = pd.DataFrame(
                columns=[feature, "mean", "minimum", "maximum", "dispersion", "samples"]
            )
        else:
            trajectory = (
                sample_values.groupby(feature, sort=False, dropna=False)["value"]
                .agg(
                    mean="mean",
                    minimum="min",
                    maximum="max",
                    dispersion="std",
                    samples="size",
                )
                .reset_index()
            )
        trajectory = pd.DataFrame({feature: timepoint_order}).merge(
            trajectory,
            on=feature,
            how="left",
            sort=False,
        )
        value_columns = ["mean", "minimum", "maximum", "dispersion"]
        trajectory[value_columns] = trajectory[value_columns].fillna(0)
        trajectory["samples"] = trajectory["samples"].fillna(0).astype(int)
        return trajectory

    def timepoint_trajectory(
        self,
        treeId,
        timepoint_feature=None,
        by_freq=True,
        ax=None,
    ):
        """Plot mean lineage abundance and sample dispersion over time."""
        import matplotlib.pyplot as plt

        feature = "timepoint" if timepoint_feature is None else timepoint_feature
        trajectory = self._get_timepoint_trajectory_df(
            treeId,
            timepoint_feature=timepoint_feature,
            by_freq=by_freq,
        )
        if trajectory is None:
            return None
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))

        x_positions = np.arange(len(trajectory))
        means = trajectory["mean"].to_numpy(dtype=float)
        lower_errors = means - trajectory["minimum"].to_numpy(dtype=float)
        upper_errors = trajectory["maximum"].to_numpy(dtype=float) - means
        ax.plot(x_positions, means, color="#F05670", linewidth=1.5, zorder=2)
        ax.scatter(x_positions, means, color="0.15", s=32, zorder=3)
        ax.errorbar(
            x_positions,
            means,
            yerr=np.vstack([lower_errors, upper_errors]),
            fmt="none",
            ecolor="0.25",
            elinewidth=1,
            capsize=3,
            zorder=1,
        )
        ax.set_xticks(x_positions)
        ax.set_xticklabels(trajectory[feature].map(_format_axis_value))
        ax.set_xlabel(feature)
        ax.set_ylabel(
            "Lineage fraction in repertoire by UMI count"
            if by_freq else "Lineage UMI count"
        )
        ax.set_title(f"Tree {treeId} timepoint trajectory")
        ax.grid(axis="y", color="0.9", linewidth=0.8)
        return ax

    def _get_timepoint_isotypes_df(
        self,
        treeId,
        timepoint_feature=None,
        by_freq=True,
    ):
        feature = "timepoint" if timepoint_feature is None else timepoint_feature
        if self.trees_df is None or self.metadata is None:
            print(
                "Cannot plot timepoint isotypes: first read the trees table and "
                "run ta.read_metadata(metadata)"
            )
            return None
        if feature not in self.trees_df.columns:
            print(
                f"Cannot plot timepoint isotypes: trees_df does not contain "
                f"'{feature}'. Run ta.read_metadata(metadata) with this column, "
                "or specify timepoint_feature='column_name'"
            )
            return None

        value_column = "uniqueMoleculeFraction" if by_freq else "uniqueMoleculeCount"
        data = self._ensure_umi_values()
        if value_column not in data.columns:
            print(
                f"Cannot plot timepoint isotypes: trees_df does not contain "
                f"'{value_column}'"
            )
            return None
        if "isotype" not in data.columns:
            print("Cannot plot timepoint isotypes: trees_df does not contain 'isotype'")
            return None

        observed = _observed_mask(data)
        all_timepoints = data.loc[observed, feature].dropna()
        if all_timepoints.empty:
            print(f"Cannot plot timepoint isotypes: trees_df has no '{feature}' values")
            return None
        timepoint_values = all_timepoints.drop_duplicates().tolist()
        if isinstance(all_timepoints.dtype, pd.CategoricalDtype):
            present = set(timepoint_values)
            timepoint_order = [
                value for value in all_timepoints.cat.categories if value in present
            ]
        else:
            try:
                timepoint_order = sorted(timepoint_values)
            except TypeError:
                timepoint_order = sorted(timepoint_values, key=lambda value: str(value))

        tree_id_key = _id_key(treeId)
        tree_rows = data.loc[
            (data["treeId"].map(_id_key) == tree_id_key) & observed
        ].copy()
        if tree_rows.empty:
            raise ValueError(f"treeId {treeId!r} does not contain observed nodes")
        tree_rows = tree_rows.dropna(subset=["sample_id", feature])
        if tree_rows.empty:
            raise ValueError(
                f"treeId {treeId!r} has no observed nodes with sample and timepoint values"
            )

        tree_rows["_isotype"] = tree_rows["isotype"].map(rsplot._recode_isotype)
        tree_rows[value_column] = pd.to_numeric(
            tree_rows[value_column], errors="coerce"
        )
        isotype_order = rsplot._isotype_order(tree_rows["_isotype"])
        sample_timepoints = tree_rows.loc[:, ["sample_id", feature]].drop_duplicates()
        conflicting_samples = sample_timepoints["sample_id"].duplicated(keep=False)
        if conflicting_samples.any():
            samples = sample_timepoints.loc[
                conflicting_samples, "sample_id"
            ].drop_duplicates().tolist()
            raise ValueError(f"Samples have multiple timepoint values: {samples}")

        sample_values = (
            tree_rows.groupby(
                ["sample_id", "_isotype"], sort=False, dropna=False
            )[value_column]
            .sum(min_count=1)
            .reset_index(name="value")
        )
        sample_grid = sample_timepoints.assign(_key=1).merge(
            pd.DataFrame({"_isotype": isotype_order, "_key": 1}),
            on="_key",
            how="inner",
        ).drop(columns="_key")
        sample_grid = sample_grid.merge(
            sample_values,
            on=["sample_id", "_isotype"],
            how="left",
        )
        sample_grid["value"] = sample_grid["value"].fillna(0)
        means = (
            sample_grid.groupby([feature, "_isotype"], sort=False, dropna=False)[
                "value"
            ]
            .mean()
            .reset_index(name="mean")
        )

        trajectory = pd.DataFrame({feature: timepoint_order}).assign(_key=1).merge(
            pd.DataFrame({"_isotype": isotype_order, "_key": 1}),
            on="_key",
            how="inner",
        ).drop(columns="_key")
        trajectory = trajectory.merge(
            means,
            on=[feature, "_isotype"],
            how="left",
            sort=False,
        )
        trajectory["mean"] = trajectory["mean"].fillna(0)
        trajectory["_isotype"] = pd.Categorical(
            trajectory["_isotype"],
            categories=isotype_order,
            ordered=True,
        )
        return trajectory.rename(columns={"_isotype": "isotype"})

    def timepoint_isotypes(
        self,
        treeId,
        timepoint_feature=None,
        by_freq=True,
        ax=None,
    ):
        """Plot mean sample-level isotype abundance over time without error bars."""
        import matplotlib.pyplot as plt

        feature = "timepoint" if timepoint_feature is None else timepoint_feature
        trajectory = self._get_timepoint_isotypes_df(
            treeId,
            timepoint_feature=timepoint_feature,
            by_freq=by_freq,
        )
        if trajectory is None:
            return None
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))

        isotype_order = list(trajectory["isotype"].cat.categories)
        palette = rsplot._isotype_colors(isotype_order)
        timepoint_order = trajectory[feature].drop_duplicates().tolist()
        x_positions = np.arange(len(timepoint_order))
        for isotype in isotype_order:
            values = trajectory.loc[
                trajectory["isotype"] == isotype, "mean"
            ].to_numpy(dtype=float)
            ax.plot(
                x_positions,
                values,
                color=palette[isotype],
                marker="o",
                linewidth=1.8,
                markersize=5,
                label=isotype,
            )

        ax.set_xticks(x_positions)
        ax.set_xticklabels([_format_axis_value(value) for value in timepoint_order])
        ax.set_xlabel(feature)
        ax.set_ylabel(
            "Mean isotype fraction by UMI count"
            if by_freq else "Mean isotype UMI count"
        )
        ax.set_title(f"Tree {treeId} isotypes by timepoint")
        ax.grid(axis="y", color="0.9", linewidth=0.8)
        ax.legend(
            title="Isotype",
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.18),
            ncol=min(5, len(isotype_order)),
        )
        return ax

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
