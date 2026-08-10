"""Stateful differential-enrichment analysis workflow."""

from copy import deepcopy

import pandas as pd

from .clone_filter import Filter


CHAIN_SPECIFIC_PARAMETERS = {
    "cl_filter",
    "overlap_type",
    "mismatches",
    "clustering",
    "count_by_freq",
    "min_samples",
    "min_count",
    "min_total_count",
    "method",
    "p_adjust_method",
    "presence_threshold",
    "sample_totals",
    "hurdle_combine_method",
    "cpm_scale",
    "pseudocount",
    "n_permutations",
    "negative_binomial_alpha",
    "max_p_adj",
    "max_p_val",
    "min_group_mean",
    "min_logfc",
    "sort",
}
KNOWN_CHAINS = {"TRA", "TRB", "TRG", "TRD", "TRAD", "IGH", "IGK", "IGL", "IGKL", "XCR"}
COUNT_PARAMETERS = {"cl_filter", "overlap_type", "mismatches", "clustering", "count_by_freq"}
PREFILTER_PARAMETERS = {"min_samples", "min_count", "min_total_count"}
STATISTICS_PARAMETERS = {
    "method",
    "simplify",
    "p_adjust_method",
    "log2fc_zero_value",
    "presence_threshold",
    "sample_totals",
    "hurdle_combine_method",
    "cpm_scale",
    "pseudocount",
    "n_permutations",
    "negative_binomial_alpha",
    "sort",
}
POSTFILTER_PARAMETERS = {
    "max_p_adj",
    "max_p_val",
    "min_group_mean",
    "min_logfc",
    "groups",
    "groups_exclude",
    "sort",
}
RUNTIME_PARAMETERS = {"cpu", "verbose"}
STAGES = ("count_table", "prefiltered", "statistics_df", "postfiltered")
PARAMETER_STAGE = {
    **{name: 0 for name in COUNT_PARAMETERS},
    **{name: 1 for name in PREFILTER_PARAMETERS},
    **{name: 2 for name in STATISTICS_PARAMETERS},
    **{name: 3 for name in POSTFILTER_PARAMETERS},
}


def _default_parameters(default_sort_columns):
    return {
        "samples_df": None,
        "cl_filter": Filter(functionality="f", by_umi=True),
        "overlap_type": "aaV",
        "mismatches": 1,
        "clustering": False,
        "count_by_freq": False,
        "min_samples": 3,
        "min_count": 2,
        "min_total_count": 10,
        "method": "mann_whitney",
        "simplify": True,
        "p_adjust_method": "fdr_bh",
        "log2fc_zero_value": 100,
        "presence_threshold": 2,
        "sample_totals": None,
        "hurdle_combine_method": "fisher",
        "cpm_scale": 1_000_000,
        "pseudocount": 0.5,
        "n_permutations": 10_000,
        "negative_binomial_alpha": 1.0,
        "max_p_adj": None,
        "max_p_val": None,
        "min_group_mean": 2,
        "min_logfc": 1,
        "groups": None,
        "groups_exclude": None,
        "sort": list(default_sort_columns),
        "cpu": None,
        "verbose": True,
        "pair_chains": True,
        "pairing_method": "jsd",
    }


def _values_equal(left, right):
    if isinstance(left, pd.DataFrame) and isinstance(right, pd.DataFrame):
        return left.equals(right)
    if isinstance(left, pd.Series) and isinstance(right, pd.Series):
        return left.equals(right)
    if type(left) is type(right) and hasattr(left, "__dict__"):
        return _values_equal(vars(left), vars(right))
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            _values_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return len(left) == len(right) and all(
            _values_equal(a, b) for a, b in zip(left, right)
        )
    try:
        result = left == right
        return bool(result) if not hasattr(result, "all") else bool(result.all())
    except (TypeError, ValueError):
        return False


class _AnalyzerTable(pd.DataFrame):
    _metadata = ["_analyzer", "_result_name"]

    @property
    def _constructor(self):
        return pd.DataFrame

    def __call__(self, chain=None):
        return self._analyzer._get_result(self._result_name, chain=chain)


class Analyzer:
    """Run and cache a complete differential-enrichment analysis workflow."""

    def __init__(self, parameters=None, **kwargs):
        from .diff_enrichment import DEFAULT_SORT_COLUMNS

        if parameters is not None and not isinstance(parameters, dict):
            raise TypeError("parameters must be a dictionary")
        self._parameters = _default_parameters(DEFAULT_SORT_COLUMNS)
        self._samples_df = None
        self._source_samples_df = None
        self._chains = []
        self._active_chain = None
        self._results = {}
        self._signatures = {}
        self._pairing_matrix = None
        updates = dict(parameters or {})
        updates.update(kwargs)
        samples_df = updates.pop("samples_df", None)
        if updates:
            self.update_parameters(updates)
        if samples_df is not None:
            self.read_samples(samples_df)

    def _split_parameter_name(self, name):
        if name in self._parameters:
            return name, None
        if "_" in name:
            base, chain = name.rsplit("_", 1)
            if base in CHAIN_SPECIFIC_PARAMETERS and chain in KNOWN_CHAINS:
                return base, chain
        raise ValueError(f"Unknown Analyzer parameter: {name!r}")

    def get_parameters(self):
        parameters = deepcopy(self._parameters)
        parameters["samples_df"] = None if self._samples_df is None else self._samples_df.copy()
        return parameters

    def update_parameters(self, parameters=None, **kwargs):
        if parameters is not None and not isinstance(parameters, dict):
            raise TypeError("parameters must be a dictionary")
        updates = dict(parameters or {})
        updates.update(kwargs)
        if "samples_df" in updates:
            self.read_samples(updates.pop("samples_df"))
        for name, value in updates.items():
            base, chain = self._split_parameter_name(name)
            old_value = self._parameters.get(name)
            if name in self._parameters and _values_equal(old_value, value):
                continue
            if chain is not None and chain not in self._chains and self._chains:
                raise ValueError(f"Chain {chain!r} is not present in samples_df")
            if base == "pair_chains" and value and self._samples_df is not None:
                if len(self._chains) == 2 and "sample" not in self._samples_df.columns:
                    raise ValueError("samples_df must contain a 'sample' column when pairing two chains")
            self._parameters[name] = deepcopy(value)
            if base in RUNTIME_PARAMETERS or base == "sort":
                continue
            if base in PARAMETER_STAGE:
                self._invalidate(PARAMETER_STAGE[base], chains=[chain] if chain else None)
            elif base in {"pair_chains", "pairing_method"}:
                self._pairing_matrix = None
        return self

    @property
    def samples_df(self):
        return self._samples_df

    @property
    def chains(self):
        if not self._chains:
            print("No samples_df has been read.")
            return []
        print(f"Active chain: {self._active_chain}")
        return list(self._chains)

    def read_samples(self, samples_df):
        if samples_df is None:
            changed = self._source_samples_df is not None
            self._samples_df = self._source_samples_df = None
            self._chains = []
            self._active_chain = None
            if changed:
                self._reset_results()
            return self
        if not isinstance(samples_df, pd.DataFrame):
            raise TypeError("samples_df must be a pandas DataFrame")
        if self._source_samples_df is not None and self._source_samples_df.equals(samples_df):
            return self
        effective = self._validate_and_normalize_samples(samples_df)
        self._source_samples_df = samples_df.copy()
        self._samples_df = effective
        self._chains = effective["chain"].drop_duplicates().tolist()
        self._active_chain = self._chains[0]
        self._parameters["samples_df"] = None
        self._reset_results()
        return self

    def _validate_and_normalize_samples(self, samples_df):
        required = {"sample_id", "filename", "group"}
        missing = sorted(required.difference(samples_df.columns))
        if missing:
            raise ValueError(f"samples_df is missing required columns: {missing}")
        effective = samples_df.copy()
        if "chain" not in effective.columns:
            if effective["sample_id"].duplicated().any():
                raise ValueError("sample_id values must be unique")
            effective["chain"] = "XCR"
            return effective
        if effective["chain"].isna().any():
            raise ValueError("chain values must not be missing")
        effective["chain"] = effective["chain"].astype(str).str.upper()
        if effective.duplicated(["chain", "sample_id"]).any():
            raise ValueError("sample_id values must be unique within each chain")
        chains = set(effective["chain"])
        if "IGH" in chains and {"IGK", "IGL"}.issubset(chains):
            raise ValueError("IGH with both IGK and IGL is not supported; merge IGK and IGL manually into IGKL")
        if len(chains) > 2:
            raise ValueError("samples_df may contain at most two chains")
        rename = {}
        if "TRAD" in chains:
            if "TRB" in chains:
                rename["TRAD"] = "TRA"
            elif "TRG" in chains:
                rename["TRAD"] = "TRD"
        if "IGH" in chains and "IGK" in chains:
            rename["IGK"] = "IGKL"
        if "IGH" in chains and "IGL" in chains:
            rename["IGL"] = "IGKL"
        if rename:
            for old, new in rename.items():
                if self._parameters["verbose"]:
                    print(f"Renaming chain {old} to {new} for differential-enrichment analysis.")
            effective["chain"] = effective["chain"].replace(rename)
        chains = effective["chain"].drop_duplicates().tolist()
        if len(chains) == 2:
            pair = frozenset(chains)
            supported = {frozenset(("TRA", "TRB")), frozenset(("TRG", "TRD")), frozenset(("IGH", "IGKL"))}
            if pair not in supported:
                raise ValueError(
                    f"Unsupported chain combination {chains}. Supported pairs are TRA-TRB, TRG-TRD, and IGH-IGKL."
                )
            if self._parameters["pair_chains"] and "sample" not in effective.columns:
                raise ValueError("samples_df must contain a 'sample' column when pairing two chains")
        return effective

    def select_chain(self, chain):
        if not self._chains:
            raise RuntimeError("Read samples_df before selecting a chain")
        if len(self._chains) == 1:
            print(f"Only one chain ({self._chains[0]}) is available; active chain was not changed.")
            return self
        if chain not in self._chains:
            raise ValueError(f"Unknown chain {chain!r}; available chains: {self._chains}")
        self._active_chain = chain
        if self._parameters["verbose"]:
            print(f"Active chain selected: {chain}")
        return self

    def _reset_results(self):
        self._results = {chain: {stage: None for stage in STAGES} for chain in self._chains}
        self._signatures = {chain: {} for chain in self._chains}
        self._pairing_matrix = None

    def _invalidate(self, stage_index, chains=None):
        targets = self._chains if chains is None else [chain for chain in chains if chain in self._chains]
        for chain in targets:
            for stage in STAGES[stage_index:]:
                self._results[chain][stage] = None
                self._signatures[chain].pop(stage, None)
        if targets:
            self._pairing_matrix = None

    def _effective_parameter(self, name, chain=None):
        chain = chain or self._active_chain
        specific = f"{name}_{chain}"
        return self._parameters.get(specific, self._parameters[name])

    def _branch_samples(self, chain=None):
        chain = chain or self._active_chain
        return self._samples_df.loc[self._samples_df["chain"] == chain].copy()

    def _require_samples(self):
        if self._samples_df is None:
            raise RuntimeError("Read samples_df before running the analysis")

    def _update_stage_parameters(self, allowed, kwargs):
        unknown = set(kwargs).difference(allowed | RUNTIME_PARAMETERS)
        if unknown:
            raise TypeError(f"Unsupported parameters: {sorted(unknown)}")
        updates = {}
        for name, value in kwargs.items():
            if len(self._chains) == 2 and name in CHAIN_SPECIFIC_PARAMETERS:
                updates[f"{name}_{self._active_chain}"] = value
            else:
                updates[name] = value
        if updates:
            self.update_parameters(updates)

    def _signature(self, names, chain=None):
        return {name: deepcopy(self._effective_parameter(name, chain)) for name in names}

    def _already_done(self, stage, signature, chain):
        done = self._results[chain][stage] is not None and _values_equal(
            self._signatures[chain].get(stage), signature
        )
        if done:
            print(f"{stage.replace('_', ' ').title()} for chain {chain} is already calculated with these parameters; nothing to do.")
        return done

    def _store(self, stage, table, signature, chain):
        wrapped = _AnalyzerTable(table)
        wrapped.attrs = table.attrs.copy()
        wrapped._analyzer = self
        wrapped._result_name = stage
        self._results[chain][stage] = wrapped
        self._signatures[chain][stage] = deepcopy(signature)
        return wrapped

    def run_count_table(self, **kwargs):
        from . import intersections

        self._require_samples()
        self._update_stage_parameters(COUNT_PARAMETERS, kwargs)
        chain = self._active_chain
        signature = self._signature(COUNT_PARAMETERS, chain)
        if self._already_done("count_table", signature, chain):
            return self._results[chain]["count_table"]
        params = {name: self._effective_parameter(name, chain) for name in COUNT_PARAMETERS}
        samples = self._branch_samples(chain)
        cpu, verbose = self._parameters["cpu"], self._parameters["verbose"]
        if params["clustering"]:
            from . import clustering

            clusters = clustering.Clusters()
            clusters.read_from_clonosets_df(samples, cl_filter=params["cl_filter"], verbose=verbose)
            clusters.create_clusters(
                overlap_type=params["overlap_type"],
                mismatches=params["mismatches"],
                cpu=cpu,
                verbose=verbose,
            )
            table = clusters.to_count_table(by_freq=params["count_by_freq"])
        else:
            table = intersections.count_table(
                samples,
                cl_filter=params["cl_filter"],
                overlap_type=params["overlap_type"],
                mismatches=0,
                by_freq=params["count_by_freq"],
                cpu=cpu,
                verbose=verbose,
            )
        self._invalidate(1, chains=[chain])
        return self._store("count_table", table, signature, chain)

    def run_prefilter(self, **kwargs):
        from .diff_enrichment import prefilter

        self._require_samples()
        self._update_stage_parameters(PREFILTER_PARAMETERS, kwargs)
        chain = self._active_chain
        if self._results[chain]["count_table"] is None:
            raise RuntimeError("Run count-table calculation before prefiltering")
        signature = self._signature(PREFILTER_PARAMETERS, chain)
        if self._already_done("prefiltered", signature, chain):
            return self._results[chain]["prefiltered"]
        table = prefilter(
            self._results[chain]["count_table"],
            **{name: self._effective_parameter(name, chain) for name in PREFILTER_PARAMETERS},
            verbose=self._parameters["verbose"],
        )
        self._invalidate(2, chains=[chain])
        return self._store("prefiltered", table, signature, chain)

    def _statistics_signature_names(self, chain):
        method = self._effective_parameter("method", chain)
        names = {"method", "simplify", "p_adjust_method", "log2fc_zero_value", "presence_threshold", "sort"}
        if method in {"fisher_count", "hurdle", "quasi_binomial", "negative_binomial", "permutation"}:
            names.add("sample_totals")
        if method == "hurdle":
            names.update({"hurdle_combine_method", "cpm_scale", "pseudocount"})
        elif method in {"quasi_binomial", "permutation"}:
            names.update({"cpm_scale", "pseudocount"})
        if method == "permutation":
            names.add("n_permutations")
        if method == "negative_binomial":
            names.add("negative_binomial_alpha")
        return names

    def run_statistics(self, **kwargs):
        from .diff_enrichment import calc_statistics

        self._require_samples()
        self._update_stage_parameters(STATISTICS_PARAMETERS, kwargs)
        chain = self._active_chain
        if self._results[chain]["prefiltered"] is None:
            raise RuntimeError("Run prefiltering before calculating statistics")
        signature_names = self._statistics_signature_names(chain)
        signature = self._signature(signature_names, chain)
        if self._already_done("statistics_df", signature, chain):
            return self._results[chain]["statistics_df"]
        arguments = {name: self._effective_parameter(name, chain) for name in STATISTICS_PARAMETERS}
        table = calc_statistics(
            self._results[chain]["prefiltered"],
            self._branch_samples(chain),
            cpu=self._parameters["cpu"],
            verbose=self._parameters["verbose"],
            **arguments,
        )
        self._invalidate(3, chains=[chain])
        return self._store("statistics_df", table, signature, chain)

    def run_postfilter(self, **kwargs):
        from .diff_enrichment import postfilter

        self._require_samples()
        self._update_stage_parameters(POSTFILTER_PARAMETERS, kwargs)
        chain = self._active_chain
        if self._results[chain]["statistics_df"] is None:
            raise RuntimeError("Run statistics before postfiltering")
        signature = self._signature(POSTFILTER_PARAMETERS, chain)
        if self._already_done("postfiltered", signature, chain):
            return self._results[chain]["postfiltered"]
        table = postfilter(
            self._results[chain]["statistics_df"],
            verbose=self._parameters["verbose"],
            **{name: self._effective_parameter(name, chain) for name in POSTFILTER_PARAMETERS},
        )
        self._pairing_matrix = None
        return self._store("postfiltered", table, signature, chain)

    def pair_chains(self, pairing_method=None):
        from .diff_enrichment import POSTFILTER_COLUMN, pair_chains

        self._require_samples()
        if pairing_method is not None:
            self.update_parameters(pairing_method=pairing_method)
        if len(self._chains) == 1:
            print(f"Only one chain ({self._chains[0]}) is present; chain pairing was skipped.")
            return None
        missing = [chain for chain in self._chains if self._results[chain]["postfiltered"] is None]
        if missing:
            raise RuntimeError(f"Run postfiltering for both chains before pairing; missing: {missing}")
        if self._pairing_matrix is not None:
            print("Chain pairing is already calculated; nothing to do.")
            return self._pairing_matrix
        first, second = self._chains
        table1, table2 = self._results[first]["postfiltered"], self._results[second]["postfiltered"]
        ids1 = table1.loc[table1[POSTFILTER_COLUMN], table1.columns[0]].tolist()
        ids2 = table2.loc[table2[POSTFILTER_COLUMN], table2.columns[0]].tolist()
        pairing_metadata = self._samples_df
        if pairing_metadata["sample_id"].duplicated().any():
            pairing_metadata = pairing_metadata.copy()
            pairing_tables = []
            for chain, table in ((first, table1), (second, table2)):
                chain_rows = pairing_metadata["chain"].eq(chain)
                sample_ids = pairing_metadata.loc[chain_rows, "sample_id"].tolist()
                renamed_ids = {sample_id: f"{sample_id}__{chain}" for sample_id in sample_ids}
                pairing_metadata.loc[chain_rows, "sample_id"] = [
                    renamed_ids[sample_id] for sample_id in sample_ids
                ]
                pairing_tables.append(table.rename(columns=renamed_ids))
            table1, table2 = pairing_tables
        self._pairing_matrix = pair_chains(
            table1,
            table2,
            pairing_metadata,
            method=self._parameters["pairing_method"],
            filter_ids1=ids1,
            filter_ids2=ids2,
        )
        return self._pairing_matrix

    def run(self):
        self._require_samples()
        original_chain = self._active_chain
        for chain in self._chains:
            self._active_chain = chain
            for label, runner in (
                ("count table", self.run_count_table),
                ("prefilter", self.run_prefilter),
                ("statistics", self.run_statistics),
                ("postfilter", self.run_postfilter),
            ):
                if self._parameters["verbose"]:
                    print(f"Running {label} for chain {chain}.")
                runner()
        self._active_chain = original_chain
        if len(self._chains) == 2 and self._parameters["pair_chains"]:
            if self._parameters["verbose"]:
                print("Pairing chains.")
            self.pair_chains()
        return self

    def _get_result(self, name, chain=None):
        self._require_samples()
        selected = chain or self._active_chain
        if selected not in self._chains:
            raise ValueError(f"Unknown chain {selected!r}; available chains: {self._chains}")
        table = self._results[selected][name]
        if table is None:
            print(f"{name.replace('_', ' ').title()} for chain {selected} has not been calculated yet.")
        return table

    @property
    def count_table(self):
        return self._get_result("count_table")

    @property
    def prefiltered(self):
        return self._get_result("prefiltered")

    @property
    def statistics_df(self):
        return self._get_result("statistics_df")

    @property
    def postfiltered(self):
        return self._get_result("postfiltered")

    @property
    def pairing_matrix(self):
        if self._pairing_matrix is None:
            print("Pairing matrix has not been calculated yet.")
        return self._pairing_matrix

    def plot_volcano(self, chain=None, postfiltered=True, **kwargs):
        from .plot import de_volcano

        self._require_samples()
        selected = chain or self._active_chain
        if selected not in self._chains:
            raise ValueError(f"Unknown chain {selected!r}; available chains: {self._chains}")
        table = None
        if postfiltered:
            table = self._results[selected]["postfiltered"]
        if table is None:
            table = self._get_result("statistics_df", chain=selected)
        if table is None:
            return None
        return de_volcano(table, **kwargs)

    def plot_heatmap(self, chain=None, **kwargs):
        from .plot import de_heatmap

        selected = chain or self._active_chain
        table = self._get_result("postfiltered", chain=selected)
        if table is None:
            return None
        return de_heatmap(table, self._branch_samples(selected), **kwargs)

    def __repr__(self):
        if self._samples_df is None:
            return "Analyzer(samples_df not read; no analysis results)"
        sample_column = "sample" if "sample" in self._samples_df.columns else "sample_id"
        lines = [
            "Differential enrichment Analyzer",
            f"Samples: {self._samples_df[sample_column].nunique()}",
            f"Groups: {', '.join(map(str, self._samples_df['group'].drop_duplicates()))}",
            f"Chains: {', '.join(self._chains)} (active: {self._active_chain})",
        ]
        for chain in self._chains:
            state = []
            for stage in STAGES:
                table = self._results[chain][stage]
                if table is not None:
                    detail = f" ({len(table)} rows)"
                    if stage == "prefiltered" and "prefilter_pass" in table:
                        detail = f" ({int(table['prefilter_pass'].sum())} of {len(table)} passed)"
                    if stage == "postfiltered" and "postfilter_pass" in table:
                        detail = f" ({int(table['postfilter_pass'].sum())} of {len(table)} passed)"
                    state.append(stage.replace("_df", "") + detail)
            lines.append(f"{chain}: " + (", ".join(state) if state else "no steps calculated"))
        if self._pairing_matrix is not None:
            lines.append(f"Pairing matrix: {self._pairing_matrix.shape[0]} x {self._pairing_matrix.shape[1]}")
        return "\n".join(lines)

    __str__ = __repr__
