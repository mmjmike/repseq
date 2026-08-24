import pandas as pd
from .common_functions import run_parallel_calculation, print_progress_bar
from .logo import create_motif_dict, sum_motif_dicts, get_consensus_from_motif_dict, get_logo_for_list_of_clonotypes
from .clone_filter import Filter
from .io import read_clonoset
from .plot import _isotype_colors, _isotype_order, _recode_isotype
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from .common_functions import overlap_type_to_flags, overlap_type_uses_sequence
import numpy as np
import networkx as nx
from networkx.algorithms import community
from scipy.sparse import csr_matrix
import os
import json
import copy
import functools
import math
import operator
import re
from collections.abc import Iterable
from numbers import Real

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D


_NODE_PROPERTY_ALIASES = {
    "cdr3aa": "seq_aa",
    "cdr3nt": "seq_nt",
}


def _recode_cluster_isotype(constant_call):
    """Recode IGH constants and return None for non-IGH chains."""
    if pd.isna(constant_call):
        return None
    normalized_call = str(constant_call).strip().upper()
    normalized_call = re.split(r"[,;|]", normalized_call, maxsplit=1)[0]
    normalized_call = normalized_call.split("(", 1)[0].split("*", 1)[0]
    if not normalized_call.startswith("IGH"):
        return None
    try:
        return _recode_isotype(constant_call)
    except ValueError:
        return None


def _node_property_value(node, property_name):
    if not isinstance(property_name, str) or not property_name:
        raise TypeError("Node property names must be non-empty strings.")
    attribute_name = _NODE_PROPERTY_ALIASES.get(property_name, property_name)
    if hasattr(node, attribute_name):
        return getattr(node, attribute_name)
    if property_name in node.additional_properties:
        return node.additional_properties[property_name]
    if attribute_name in node.additional_properties:
        return node.additional_properties[attribute_name]
    raise ValueError(
        f"Node property '{property_name}' was not found in node attributes "
        "or additional_properties."
    )


def _normalize_match_values(values):
    if isinstance(values, (str, bytes)) or not isinstance(values, Iterable):
        return (values,)
    if isinstance(values, (set, frozenset)):
        return tuple(sorted(values, key=str))
    return tuple(values)


def _node_matches(node, property_name, values):
    return _node_property_value(node, property_name) in values


def _node_weight(node, weight):
    if weight == "nodes":
        return 1
    value = _node_property_value(node, weight)
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(
            f"Node property '{weight}' must contain numeric values to be used "
            "as a weight."
        )
    if not np.isfinite(value) or value < 0:
        raise ValueError(
            f"Node property '{weight}' weights must be finite and non-negative."
        )
    return value


def _expression_name_part(value):
    text = re.sub(r"[^0-9A-Za-z]+", "_", str(value)).strip("_")
    return text or "value"


def _selection_expression_name(prefix, property_name, values, weight=None):
    value_name = "_or_".join(_expression_name_part(value) for value in values)
    name = f"{prefix}_{_expression_name_part(property_name)}_{value_name}"
    if weight is not None:
        name += f"_by_{_expression_name_part(weight)}"
    return name


class ClusterExpression:
    """Lazy cluster-level value used for filtering and property calculation."""

    is_boolean = False

    def __init__(self, name):
        self.name = name

    def evaluate(self, cluster):
        raise NotImplementedError

    def alias(self, name):
        if not isinstance(name, str) or not name:
            raise TypeError("Expression aliases must be non-empty strings.")
        return _AliasedClusterExpression(self, name)

    def _compare(self, other, comparison, symbol):
        return _ComparisonExpression(self, other, comparison, symbol)

    def __lt__(self, other):
        return self._compare(other, operator.lt, "<")

    def __le__(self, other):
        return self._compare(other, operator.le, "<=")

    def __eq__(self, other):
        return self._compare(other, operator.eq, "==")

    def __ne__(self, other):
        return self._compare(other, operator.ne, "!=")

    def __ge__(self, other):
        return self._compare(other, operator.ge, ">=")

    def __gt__(self, other):
        return self._compare(other, operator.gt, ">")

    def __and__(self, other):
        return _LogicalExpression(self, other, operator.and_, "&")

    def __or__(self, other):
        return _LogicalExpression(self, other, operator.or_, "|")

    def __invert__(self):
        return _NotExpression(self)

    def __bool__(self):
        raise TypeError(
            "Cluster expressions cannot be converted to bool directly. Use "
            "comparison operators and combine predicates with &, |, and ~."
        )

    def __repr__(self):
        return self.name


class _AliasedClusterExpression(ClusterExpression):
    def __init__(self, expression, name):
        super().__init__(name)
        self.expression = expression
        self.is_boolean = expression.is_boolean

    def evaluate(self, cluster):
        return self.expression.evaluate(cluster)


class _FunctionClusterExpression(ClusterExpression):
    def __init__(self, name, evaluator, is_boolean=False):
        super().__init__(name)
        self.evaluator = evaluator
        self.is_boolean = is_boolean

    def evaluate(self, cluster):
        return self.evaluator(cluster)


class _ComparisonExpression(ClusterExpression):
    is_boolean = True

    def __init__(self, left, right, comparison, symbol):
        if not isinstance(left, ClusterExpression):
            raise TypeError("The left comparison value must be a cluster expression.")
        right_name = right.name if isinstance(right, ClusterExpression) else repr(right)
        super().__init__(f"({left.name} {symbol} {right_name})")
        self.left = left
        self.right = right
        self.comparison = comparison

    def evaluate(self, cluster):
        left_value = self.left.evaluate(cluster)
        right_value = (
            self.right.evaluate(cluster)
            if isinstance(self.right, ClusterExpression)
            else self.right
        )
        return bool(self.comparison(left_value, right_value))


class _LogicalExpression(ClusterExpression):
    is_boolean = True

    def __init__(self, left, right, logical_operator, symbol):
        if not isinstance(left, ClusterExpression) or not left.is_boolean:
            raise TypeError(f"The left operand of {symbol} must be a predicate.")
        if not isinstance(right, ClusterExpression) or not right.is_boolean:
            raise TypeError(f"The right operand of {symbol} must be a predicate.")
        super().__init__(f"({left.name} {symbol} {right.name})")
        self.left = left
        self.right = right
        self.logical_operator = logical_operator

    def evaluate(self, cluster):
        left_value = bool(self.left.evaluate(cluster))
        if self.logical_operator is operator.and_:
            return left_value and bool(self.right.evaluate(cluster))
        return left_value or bool(self.right.evaluate(cluster))


class _NotExpression(ClusterExpression):
    is_boolean = True

    def __init__(self, expression):
        if not isinstance(expression, ClusterExpression) or not expression.is_boolean:
            raise TypeError("The operand of ~ must be a predicate.")
        super().__init__(f"~({expression.name})")
        self.expression = expression

    def evaluate(self, cluster):
        return not bool(self.expression.evaluate(cluster))


class _TotalCountExpression(ClusterExpression):
    def __init__(self, property_name=None, values=None, weight="count"):
        if property_name is None:
            name = "total_count" if weight == "count" else f"total_{weight}"
            normalized_values = None
        else:
            normalized_values = _normalize_match_values(values)
            name = _selection_expression_name(
                "total_count", property_name, normalized_values, weight
            )
        super().__init__(name)
        self.property_name = property_name
        self.values = normalized_values
        self.weight = weight

    def __call__(self, property_name, values, weight="count"):
        return _TotalCountExpression(property_name, values, weight=weight)

    def evaluate(self, cluster):
        return sum(
            _node_weight(node, self.weight)
            for node in cluster
            if self.property_name is None
            or _node_matches(node, self.property_name, self.values)
        )


cluster_size = _FunctionClusterExpression("cluster_size", len)
total_count = _TotalCountExpression()


def proportion(property_name, values, weight="nodes"):
    normalized_values = _normalize_match_values(values)
    name = _selection_expression_name(
        "proportion", property_name, normalized_values, weight
    )

    def evaluate(cluster):
        denominator = sum(_node_weight(node, weight) for node in cluster)
        if denominator == 0:
            return 0.0
        numerator = sum(
            _node_weight(node, weight)
            for node in cluster
            if _node_matches(node, property_name, normalized_values)
        )
        return numerator / denominator

    return _FunctionClusterExpression(name, evaluate)


def all_nodes(property_name, values):
    normalized_values = _normalize_match_values(values)
    name = _selection_expression_name("all_nodes", property_name, normalized_values)

    def evaluate(cluster):
        nodes = list(cluster)
        return bool(nodes) and all(
            _node_matches(node, property_name, normalized_values) for node in nodes
        )

    return _FunctionClusterExpression(name, evaluate, is_boolean=True)


def any_nodes(property_name, values):
    normalized_values = _normalize_match_values(values)
    name = _selection_expression_name("any_nodes", property_name, normalized_values)

    def evaluate(cluster):
        return any(
            _node_matches(node, property_name, normalized_values) for node in cluster
        )

    return _FunctionClusterExpression(name, evaluate, is_boolean=True)


_CLUSTER_PROPERTY_COLUMNS = [
    "cluster_no",
    "cluster_id",
    "nodes",
    "edges",
    "total_count",
    "diameter",
    "density",
    "eccentricity",
    "concensus_cdr3aa",
    "concensus_cdr3nt",
    "concensus_v",
    "concensus_j",
]


def _cluster_properties_worker(args):
    cluster_no, cluster, use_first_v, use_first_j = args
    nodes = list(cluster)
    if not nodes:
        raise ValueError(f"Cluster {cluster_no} does not contain any nodes.")

    first_node = nodes[0]
    node_count = len(nodes)
    edge_count = cluster.number_of_edges()
    total_node_count = sum(node.count for node in nodes)
    if node_count == 1:
        diameter = 0
        density = 0
        average_eccentricity = 0
        aa_consensus = first_node.seq_aa
        nt_consensus = first_node.seq_nt
        v_consensus = first_node.v
        j_consensus = first_node.j
    else:
        diameter = nx.diameter(cluster)
        density = nx.density(cluster)
        average_eccentricity = np.mean(
            list(nx.eccentricity(cluster).values())
        )
        aa_consensus = cluster.calc_cluster_consensus(
            seq_type="prot", weigh_by=None
        )
        nt_consensus = cluster.calc_cluster_consensus(
            seq_type="dna", weigh_by=None
        )
        v_consensus = (
            first_node.v
            if use_first_v
            else cluster.calc_cluster_consensus_segment(
                segment_type="v", weigh_by=None
            )
        )
        j_consensus = (
            first_node.j
            if use_first_j
            else cluster.calc_cluster_consensus_segment(
                segment_type="j", weigh_by=None
            )
        )

    return (
        cluster_no,
        f"cluster_{cluster_no}",
        node_count,
        edge_count,
        total_node_count,
        diameter,
        density,
        average_eccentricity,
        aa_consensus,
        nt_consensus,
        v_consensus,
        j_consensus,
    )


def _intersect_clusters_with_clonoset_worker(args):
    (
        sample_id,
        filename,
        cl_filter,
        cluster_dicts,
        overlap_type,
        mismatches,
        by_freq,
    ) = args
    from .intersections import (
        _similarity_pair_matches,
        _similarity_total_count,
        prepare_clonoset_for_intersection,
    )

    clonoset = read_clonoset(filename)
    clonoset = cl_filter.apply(clonoset)
    uses_sequence = overlap_type_uses_sequence(overlap_type)
    target_dict = prepare_clonoset_for_intersection(
        clonoset,
        overlap_type=overlap_type,
        by_freq=False,
        len_vj_format=uses_sequence,
    )
    target_total = _similarity_total_count(target_dict, uses_sequence)

    values = {}
    for cluster_no, comparison_dict in cluster_dicts.items():
        matched_target_counts = {}
        matches = _similarity_pair_matches(
            target_dict,
            comparison_dict,
            overlap_type,
            mismatches,
        )
        for target_clone, _, target_count, _, _ in matches:
            matched_target_counts.setdefault(target_clone, target_count)
        value = float(sum(matched_target_counts.values()))
        if by_freq:
            value = value / target_total if target_total else 0.0
        values[cluster_no] = value
    return sample_id, values


class ClusterCollectionSelector:
    """Collection-level selector accepted by :meth:`Clusters.filter`."""

    def select(self, clusters):
        raise NotImplementedError


class _TopClustersSelector(ClusterCollectionSelector):
    def __init__(self, number_of_clusters):
        if (
            not isinstance(number_of_clusters, (int, np.integer))
            or isinstance(number_of_clusters, bool)
            or number_of_clusters < 0
        ):
            raise ValueError("number_of_clusters must be a non-negative integer.")
        self.number_of_clusters = int(number_of_clusters)

    def select(self, clusters):
        ranked_clusters = sorted(
            clusters.clusters,
            key=lambda cluster: (
                -len(cluster),
                -sum(node.count for node in cluster),
                str(
                    cluster.calc_cluster_consensus(
                        seq_type="prot", weigh_by=None
                    )
                ),
            ),
        )
        return ranked_clusters[: self.number_of_clusters]

    def __repr__(self):
        return f"top_clusters({self.number_of_clusters})"


def top_clusters(number_of_clusters):
    """Select the largest clusters using size, count, and AA consensus order."""
    return _TopClustersSelector(number_of_clusters)


# ? add freq, count
class Node:
    def __init__(self, node_id, seq_nt, seq_aa, v, j, sample_id, freq, count):
        self.id = node_id
        self.v = v
        self.j = j
        self.seq_aa = seq_aa
        self.seq_nt = seq_nt
        self.sample_id = sample_id
        self.freq = freq
        self.count = count
        self.additional_properties = {}
        
    def is_neighbour_of(self, other, mismatches=1, aa=True, check_v=False, check_j=False):
        """function compares two strings and returns
        True if their are equal
            or if they have one mismatch and equal length
        False in all other conditions
        """
        
        if aa:
            string1 = self.seq_aa
            string2 = other.seq_aa
        else:
            string1 = self.seq_nt
            string2 = other.seq_nt
        if len(string1) != len(string2):
             return False
        if check_v and self.v != other.v:
            return False
        if check_j and self.j != other.j:
            return False
        hamm_dist = sum([a != b for a,b in zip(string1,string2)]) 
        if hamm_dist > mismatches:
            return False     
        return True
    

# !!! id, cluster_no, sample_id, v, j, seq_aa, seq_nt, size: count and freq
    def __str__(self):
        return (
            f'id={self.id} | cluster_no={self.additional_properties["cluster_no"]} | V={self.v} | J={self.j} | '
            f'aa={self.seq_aa} | nt={self.seq_nt} | '
            f'count={self.count} | freq={self.freq:.3e}')
    


    def add_properties(self, metadata):
        if self.sample_id not in metadata:
            for property in list(metadata[list(metadata)[0]].keys()):
                self.additional_properties[property] = None
        else:
            self.additional_properties.update(metadata[self.sample_id])


class Cluster(nx.Graph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.id = None
    

    def get_node_by_index(self, index):
        return list(self.nodes)[index]


    def plot_cluster_logo(self, seq_type="prot", weigh_by=None):
        """
        Plots cluster sequence logo.

        Args:
            seq_type (str): "prot" or "dna".
            weigh_by (str | None): "freq", "count", or None.

        Returns:
            None
        """
        list_of_clonotypes = []
        seq_types = ["prot", "dna"]
        if seq_type not in seq_types:
            raise ValueError(f"Wrong 'seq_type'! Possible values: {', '.join(seq_types)}")
        for node in self.nodes:
            if seq_type == "dna":
                seq = node.seq_nt
                if seq == "-":
                    raise ValueError(f"Node '{node.id}' does not have specified 'cdr3nt' value. Unable to create Logo for 'dna' seq_type")
            else:
                seq = node.seq_aa
            if weigh_by:
                if weigh_by == 'freq':
                    weight = node.freq
                    clone = (seq, weight)
                elif weigh_by == 'count':
                    weight = node.count
                    clone = (seq, weight)
            else:
                clone = (seq,)
            list_of_clonotypes.append(clone)
        get_logo_for_list_of_clonotypes(list_of_clonotypes, seq_type, plot=True)   
    

    def calc_cluster_consensus(self, seq_type="dna", weigh_by=None):
        """
        Calculates cluster consensus.

        Args:
            seq_type (str): "prot" or "dna".
            weigh_by (str | None): "freq", "count", or None.

        Returns:
            consensus_seq: Consensus sequence
        """
        motif_dicts = []
        for node in self.nodes:
            weight = 1
            if weigh_by == 'freq':
                weight = node.freq
            elif weigh_by == 'count':
                weight = node.count
            if seq_type == "dna":
                seq = node.seq_nt
                if seq == "-":
                    return "-"
            else:
                seq = node.seq_aa
            node_motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
            motif_dicts.append(node_motif_dict)
        motif_dict = sum_motif_dicts(motif_dicts)
        consensus_seq = get_consensus_from_motif_dict(motif_dict)
        return consensus_seq


    def calc_cluster_consensus_segment(self, segment_type="v", weigh_by=None):
        segments = {}
        for node in self.nodes:
            if segment_type == "v":
                segment = node.v
            else:
                segment = node.j
            weight = 1
            if weigh_by == 'freq':
                weight = node.freq
            elif weigh_by == 'count':
                weight = node.count
            if segment not in segments:
                segments[segment] = weight
            else:
                segments[segment] += weight
        best_segment = ""
        best_score = 0
        for segment in segments:
            score = segments[segment]
            if segments[segment] > best_score:
                best_score = score
                best_segment = segment
        return best_segment


class Clusters(list):
    def __init__(self):
        super().__init__() 
        self.clonosets = None
        self.clonotypes = None
        self.cl_filter = None
        self.mismatches = None
        self.overlap_type = None
        self.is_pooled = None
        self.clusters = []
        self.cluster_communities_louvain = None
        self.alice_results = None
        self.tcrdist_radius = None
        self.state = {
            "empty": True,
            "clonotypes_read": False,
            "clusters_created": False,
            "metadata_added": False,
            "node_pgen_calculated": False,
            "alice_calculated": False,
        }
        self.state_parameters = {
            "clonotypes_read": None,
            "clusters_created": None,
            "metadata_added": {"columns": []},
            "node_pgen_calculated": None,
            "alice_calculated": None,
        }
        self._metadata_frames = []
        self._clonotypes_revision = 0
        self._cluster_cache_signature = None
        self._properties_cache = None
        self.TCRDIST_BLOSUM = {('A', 'A'): 0,  ('A', 'C'): 4,  ('A', 'D'): 4,  ('A', 'E'): 4,  ('A', 'F'): 4,  ('A', 'G'): 4,  ('A', 'H'): 4,  ('A', 'I'): 4,  ('A', 'K'): 4,  ('A', 'L'): 4,  ('A', 'M'): 4,  ('A', 'N'): 4,  ('A', 'P'): 4,  ('A', 'Q'): 4,  ('A', 'R'): 4,  ('A', 'S'): 3,  ('A', 'T'): 4,  ('A', 'V'): 4,  ('A', 'W'): 4,  ('A', 'Y'): 4,  ('C', 'A'): 4,  ('C', 'C'): 0,  ('C', 'D'): 4,  ('C', 'E'): 4,  ('C', 'F'): 4,  ('C', 'G'): 4,  ('C', 'H'): 4,  ('C', 'I'): 4,  ('C', 'K'): 4,  ('C', 'L'): 4,  ('C', 'M'): 4,  ('C', 'N'): 4,  ('C', 'P'): 4,  ('C', 'Q'): 4,  ('C', 'R'): 4,  ('C', 'S'): 4,  ('C', 'T'): 4,  ('C', 'V'): 4,  ('C', 'W'): 4,  ('C', 'Y'): 4,  ('D', 'A'): 4,  ('D', 'C'): 4,  ('D', 'D'): 0,  ('D', 'E'): 2,  ('D', 'F'): 4,  ('D', 'G'): 4,  ('D', 'H'): 4,  ('D', 'I'): 4,  ('D', 'K'): 4,  ('D', 'L'): 4,  ('D', 'M'): 4,  ('D', 'N'): 3,  ('D', 'P'): 4,  ('D', 'Q'): 4,  ('D', 'R'): 4,  ('D', 'S'): 4,  ('D', 'T'): 4,  ('D', 'V'): 4,  ('D', 'W'): 4,  ('D', 'Y'): 4,  ('E', 'A'): 4,  ('E', 'C'): 4,  ('E', 'D'): 2,  ('E', 'E'): 0,  ('E', 'F'): 4,  ('E', 'G'): 4,  ('E', 'H'): 4,  ('E', 'I'): 4,  ('E', 'K'): 3,  ('E', 'L'): 4,  ('E', 'M'): 4,  ('E', 'N'): 4,  ('E', 'P'): 4,  ('E', 'Q'): 2,  ('E', 'R'): 4,  ('E', 'S'): 4,  ('E', 'T'): 4,  ('E', 'V'): 4,  ('E', 'W'): 4,  ('E', 'Y'): 4,  ('F', 'A'): 4,  ('F', 'C'): 4,  ('F', 'D'): 4,  ('F', 'E'): 4,  ('F', 'F'): 0,  ('F', 'G'): 4,  ('F', 'H'): 4,  ('F', 'I'): 4,  ('F', 'K'): 4,  ('F', 'L'): 4,  ('F', 'M'): 4,  ('F', 'N'): 4,  ('F', 'P'): 4,  ('F', 'Q'): 4,  ('F', 'R'): 4,  ('F', 'S'): 4,  ('F', 'T'): 4,  ('F', 'V'): 4,  ('F', 'W'): 3,  ('F', 'Y'): 1,  ('G', 'A'): 4,  ('G', 'C'): 4,  ('G', 'D'): 4,  ('G', 'E'): 4,  ('G', 'F'): 4,  ('G', 'G'): 0,  ('G', 'H'): 4,  ('G', 'I'): 4,  ('G', 'K'): 4,  ('G', 'L'): 4,  ('G', 'M'): 4,  ('G', 'N'): 4,  ('G', 'P'): 4,  ('G', 'Q'): 4,  ('G', 'R'): 4,  ('G', 'S'): 4,  ('G', 'T'): 4,  ('G', 'V'): 4,  ('G', 'W'): 4,  ('G', 'Y'): 4,  ('H', 'A'): 4,  ('H', 'C'): 4,  ('H', 'D'): 4,  ('H', 'E'): 4,  ('H', 'F'): 4,  ('H', 'G'): 4,  ('H', 'H'): 0,  ('H', 'I'): 4,  ('H', 'K'): 4,  ('H', 'L'): 4,  ('H', 'M'): 4,  ('H', 'N'): 3,  ('H', 'P'): 4,  ('H', 'Q'): 4,  ('H', 'R'): 4,  ('H', 'S'): 4,  ('H', 'T'): 4,  ('H', 'V'): 4,  ('H', 'W'): 4,  ('H', 'Y'): 2,  ('I', 'A'): 4,  ('I', 'C'): 4,  ('I', 'D'): 4,  ('I', 'E'): 4,  ('I', 'F'): 4,  ('I', 'G'): 4,  ('I', 'H'): 4,  ('I', 'I'): 0,  ('I', 'K'): 4,  ('I', 'L'): 2,  ('I', 'M'): 3,  ('I', 'N'): 4,  ('I', 'P'): 4,  ('I', 'Q'): 4,  ('I', 'R'): 4,  ('I', 'S'): 4,  ('I', 'T'): 4,  ('I', 'V'): 1,  ('I', 'W'): 4,  ('I', 'Y'): 4,  ('K', 'A'): 4,  ('K', 'C'): 4,  ('K', 'D'): 4,  ('K', 'E'): 3,  ('K', 'F'): 4,  ('K', 'G'): 4,  ('K', 'H'): 4,  ('K', 'I'): 4,  ('K', 'K'): 0,  ('K', 'L'): 4,  ('K', 'M'): 4,  ('K', 'N'): 4,  ('K', 'P'): 4,  ('K', 'Q'): 3,  ('K', 'R'): 2,  ('K', 'S'): 4,  ('K', 'T'): 4,  ('K', 'V'): 4,  ('K', 'W'): 4,  ('K', 'Y'): 4,  ('L', 'A'): 4,  ('L', 'C'): 4,  ('L', 'D'): 4,  ('L', 'E'): 4,  ('L', 'F'): 4,  ('L', 'G'): 4,  ('L', 'H'): 4,  ('L', 'I'): 2,  ('L', 'K'): 4,  ('L', 'L'): 0,  ('L', 'M'): 2,  ('L', 'N'): 4,  ('L', 'P'): 4,  ('L', 'Q'): 4,  ('L', 'R'): 4,  ('L', 'S'): 4,  ('L', 'T'): 4,  ('L', 'V'): 3,  ('L', 'W'): 4,  ('L', 'Y'): 4,  ('M', 'A'): 4,  ('M', 'C'): 4,  ('M', 'D'): 4,  ('M', 'E'): 4,  ('M', 'F'): 4,  ('M', 'G'): 4,  ('M', 'H'): 4,  ('M', 'I'): 3,  ('M', 'K'): 4,  ('M', 'L'): 2,  ('M', 'M'): 0,  ('M', 'N'): 4,  ('M', 'P'): 4,  ('M', 'Q'): 4,  ('M', 'R'): 4,  ('M', 'S'): 4,  ('M', 'T'): 4,  ('M', 'V'): 3,  ('M', 'W'): 4,  ('M', 'Y'): 4,  ('N', 'A'): 4,  ('N', 'C'): 4,  ('N', 'D'): 3,  ('N', 'E'): 4,  ('N', 'F'): 4,  ('N', 'G'): 4,  ('N', 'H'): 3,  ('N', 'I'): 4,  ('N', 'K'): 4,  ('N', 'L'): 4,  ('N', 'M'): 4,  ('N', 'N'): 0,  ('N', 'P'): 4,  ('N', 'Q'): 4,  ('N', 'R'): 4,  ('N', 'S'): 3,  ('N', 'T'): 4,  ('N', 'V'): 4,  ('N', 'W'): 4,  ('N', 'Y'): 4,  ('P', 'A'): 4,  ('P', 'C'): 4,  ('P', 'D'): 4,  ('P', 'E'): 4,  ('P', 'F'): 4,  ('P', 'G'): 4,  ('P', 'H'): 4,  ('P', 'I'): 4,  ('P', 'K'): 4,  ('P', 'L'): 4,  ('P', 'M'): 4,  ('P', 'N'): 4,  ('P', 'P'): 0,  ('P', 'Q'): 4,  ('P', 'R'): 4,  ('P', 'S'): 4,  ('P', 'T'): 4,  ('P', 'V'): 4,  ('P', 'W'): 4,  ('P', 'Y'): 4,  ('Q', 'A'): 4,  ('Q', 'C'): 4,  ('Q', 'D'): 4,  ('Q', 'E'): 2,  ('Q', 'F'): 4,  ('Q', 'G'): 4,  ('Q', 'H'): 4,  ('Q', 'I'): 4,  ('Q', 'K'): 3,  ('Q', 'L'): 4,  ('Q', 'M'): 4,  ('Q', 'N'): 4,  ('Q', 'P'): 4,  ('Q', 'Q'): 0,  ('Q', 'R'): 3,  ('Q', 'S'): 4,  ('Q', 'T'): 4,  ('Q', 'V'): 4,  ('Q', 'W'): 4,  ('Q', 'Y'): 4,  ('R', 'A'): 4,  ('R', 'C'): 4,  ('R', 'D'): 4,  ('R', 'E'): 4,  ('R', 'F'): 4,  ('R', 'G'): 4,  ('R', 'H'): 4,  ('R', 'I'): 4,  ('R', 'K'): 2,  ('R', 'L'): 4,  ('R', 'M'): 4,  ('R', 'N'): 4,  ('R', 'P'): 4,  ('R', 'Q'): 3,  ('R', 'R'): 0,  ('R', 'S'): 4,  ('R', 'T'): 4,  ('R', 'V'): 4,  ('R', 'W'): 4,  ('R', 'Y'): 4,  ('S', 'A'): 3,  ('S', 'C'): 4,  ('S', 'D'): 4,  ('S', 'E'): 4,  ('S', 'F'): 4,  ('S', 'G'): 4,  ('S', 'H'): 4,  ('S', 'I'): 4,  ('S', 'K'): 4,  ('S', 'L'): 4,  ('S', 'M'): 4,  ('S', 'N'): 3,  ('S', 'P'): 4,  ('S', 'Q'): 4,  ('S', 'R'): 4,  ('S', 'S'): 0,  ('S', 'T'): 3,  ('S', 'V'): 4,  ('S', 'W'): 4,  ('S', 'Y'): 4,  ('T', 'A'): 4,  ('T', 'C'): 4,  ('T', 'D'): 4,  ('T', 'E'): 4,  ('T', 'F'): 4,  ('T', 'G'): 4,  ('T', 'H'): 4,  ('T', 'I'): 4,  ('T', 'K'): 4,  ('T', 'L'): 4,  ('T', 'M'): 4,  ('T', 'N'): 4,  ('T', 'P'): 4,  ('T', 'Q'): 4,  ('T', 'R'): 4,  ('T', 'S'): 3,  ('T', 'T'): 0,  ('T', 'V'): 4,  ('T', 'W'): 4,  ('T', 'Y'): 4,  ('V', 'A'): 4,  ('V', 'C'): 4,  ('V', 'D'): 4,  ('V', 'E'): 4,  ('V', 'F'): 4,  ('V', 'G'): 4,  ('V', 'H'): 4,  ('V', 'I'): 1,  ('V', 'K'): 4,  ('V', 'L'): 3,  ('V', 'M'): 3,  ('V', 'N'): 4,  ('V', 'P'): 4,  ('V', 'Q'): 4,  ('V', 'R'): 4,  ('V', 'S'): 4,  ('V', 'T'): 4,  ('V', 'V'): 0,  ('V', 'W'): 4,  ('V', 'Y'): 4,  ('W', 'A'): 4,  ('W', 'C'): 4,  ('W', 'D'): 4,  ('W', 'E'): 4,  ('W', 'F'): 3,  ('W', 'G'): 4,  ('W', 'H'): 4,  ('W', 'I'): 4,  ('W', 'K'): 4,  ('W', 'L'): 4,  ('W', 'M'): 4,  ('W', 'N'): 4,  ('W', 'P'): 4,  ('W', 'Q'): 4,  ('W', 'R'): 4,  ('W', 'S'): 4,  ('W', 'T'): 4,  ('W', 'V'): 4,  ('W', 'W'): 0,  ('W', 'Y'): 2,  ('Y', 'A'): 4,  ('Y', 'C'): 4,  ('Y', 'D'): 4,  ('Y', 'E'): 4,  ('Y', 'F'): 1,  ('Y', 'G'): 4,  ('Y', 'H'): 2,  ('Y', 'I'): 4,  ('Y', 'K'): 4,  ('Y', 'L'): 4,  ('Y', 'M'): 4,  ('Y', 'N'): 4,  ('Y', 'P'): 4,  ('Y', 'Q'): 4,  ('Y', 'R'): 4,  ('Y', 'S'): 4,  ('Y', 'T'): 4,  ('Y', 'V'): 4,  ('Y', 'W'): 2,  ('Y', 'Y'): 0}
        self.TCRDIST_CDR3_N_CUT = 3
        self.TCRDIST_CDR3_C_CUT = 2
        self.TCRDIST_CDR3_SCORE_MULTIPLIER = 3
        
    # to enable list-like behaviour 
    def __getitem__(self, index):
        return self.clusters[index]

    def __len__(self):
        return len(self.clusters)

    def __iter__(self):
       return iter(self.clusters)

    def _sync_inferred_state(self):
        if self.clonotypes is not None:
            self.state["empty"] = False
            self.state["clonotypes_read"] = True
        if self.clusters:
            self.state["empty"] = False
            self.state["clusters_created"] = True
        if self.alice_results is not None:
            self.state["node_pgen_calculated"] = True
            self.state["alice_calculated"] = True


    def _has_clonotypes(self):
        self._sync_inferred_state()
        return self.state["clonotypes_read"]


    def _has_clusters(self):
        self._sync_inferred_state()
        return self.state["clusters_created"]


    def _require_clonotypes(self, action):
        if not self._has_clonotypes():
            raise RuntimeError(
                f"Cannot {action}: clonotypes have not been read. Run "
                "read_from_pooled_clonoset(...) or "
                "read_from_clonosets_df(...) first."
            )


    def _require_clusters(self, action):
        if not self._has_clusters():
            raise RuntimeError(
                f"Cannot {action}: clusters have not been created. Read "
                "clonotypes and run create_clusters(...) first."
            )


    def _require_alice(self, action):
        if not self.state["alice_calculated"] and self.alice_results is None:
            raise RuntimeError(
                f"Cannot {action}: ALICE has not been calculated. Run "
                "alice(...) after creating clusters first."
            )


    def _invalidate_properties_cache(self):
        self._properties_cache = None


    def _invalidate_cluster_dependent_analysis(self):
        self.cluster_communities_louvain = None
        self.alice_results = None
        self.state["node_pgen_calculated"] = False
        self.state["alice_calculated"] = False
        self.state_parameters["node_pgen_calculated"] = None
        self.state_parameters["alice_calculated"] = None


    @staticmethod
    def _filter_parameters(cl_filter):
        return {
            name: value
            for name, value in vars(cl_filter).items()
            if not name.startswith("_")
        }


    def _reset_after_clonotype_read(self):
        self.clusters = []
        self.cluster_communities_louvain = None
        self.alice_results = None
        self._metadata_frames = []
        self._cluster_cache_signature = None
        self._invalidate_properties_cache()
        self.state.update(
            {
                "empty": False,
                "clonotypes_read": True,
                "clusters_created": False,
                "metadata_added": False,
                "node_pgen_calculated": False,
                "alice_calculated": False,
            }
        )
        self.state_parameters.update(
            {
                "clusters_created": None,
                "metadata_added": {"columns": []},
                "node_pgen_calculated": None,
                "alice_calculated": None,
            }
        )
        self._clonotypes_revision += 1


    def _mark_clonotypes_read(self, source, extra_parameters=None):
        self._reset_after_clonotype_read()
        sample_count = (
            self.clonotypes["sample_id"].nunique()
            if "sample_id" in self.clonotypes.columns
            else None
        )
        parameters = {
            "source": source,
            "clonotypes": len(self.clonotypes),
            "samples": sample_count,
        }
        if extra_parameters:
            parameters.update(extra_parameters)
        self.state_parameters["clonotypes_read"] = parameters


    @staticmethod
    def _metadata_columns(metadata):
        return [column for column in metadata.columns if column != "sample_id"]


    @staticmethod
    def _apply_metadata_to_clusters(clusters, metadata):
        metadata_dict = metadata.set_index("sample_id").to_dict("index")
        for cluster in clusters:
            for node in cluster:
                node.add_properties(metadata_dict)


    def _apply_stored_metadata(self):
        for metadata in self._metadata_frames:
            self._apply_metadata_to_clusters(self.clusters, metadata)


    def _cluster_summary(self):
        total_clusters = len(self.clusters)
        multi_node_clusters = sum(len(cluster) > 1 for cluster in self.clusters)
        single_node_clusters = total_clusters - multi_node_clusters
        total_nodes = sum(len(cluster) for cluster in self.clusters)
        total_edges = sum(cluster.number_of_edges() for cluster in self.clusters)
        return {
            "total_clusters": total_clusters,
            "multi_node_clusters": multi_node_clusters,
            "single_node_clusters": single_node_clusters,
            "total_nodes": total_nodes,
            "total_edges": total_edges,
        }


    @staticmethod
    def _format_state_parameters(parameters, indent="    "):
        if not parameters:
            return []
        lines = []
        for name, value in parameters.items():
            if isinstance(value, dict):
                lines.append(f"{indent}{name}:")
                lines.extend(
                    Clusters._format_state_parameters(value, indent=indent + "  ")
                )
            else:
                lines.append(f"{indent}{name}: {value}")
        return lines


    def __str__(self):
        self._sync_inferred_state()
        lines = ["Clusters state:"]
        if self.state["empty"]:
            lines.append("  State: empty")
        else:
            completed = [
                label
                for flag, label in (
                    ("clonotypes_read", "clonotypes read"),
                    ("clusters_created", "clusters created"),
                    ("metadata_added", "metadata added"),
                    ("node_pgen_calculated", "node Pgen calculated"),
                    ("alice_calculated", "ALICE calculated"),
                )
                if self.state[flag]
            ]
            lines.append(f"  Completed: {', '.join(completed)}")

        read_parameters = self.state_parameters.get("clonotypes_read")
        if self.state["clonotypes_read"]:
            lines.append("Clonotypes:")
            lines.extend(self._format_state_parameters(read_parameters))
        else:
            lines.append("Clonotypes: not read")

        if self.state["clusters_created"]:
            summary = self._cluster_summary()
            lines.append(
                f"Graph: {summary['total_nodes']} nodes and "
                f"{summary['total_edges']} edges"
            )
            cluster_word = (
                "cluster" if summary["multi_node_clusters"] == 1 else "clusters"
            )
            single_node_word = (
                "single node"
                if summary["single_node_clusters"] == 1
                else "single nodes"
            )
            lines.append(
                f"Clusters: {summary['multi_node_clusters']} {cluster_word} "
                f"(2 or more nodes) and {summary['single_node_clusters']} "
                f"{single_node_word}. Total: {summary['total_clusters']}"
            )
            lines.append("Clustering parameters:")
            lines.extend(
                self._format_state_parameters(
                    self.state_parameters.get("clusters_created")
                )
            )
        else:
            lines.append("Clusters: not created")

        metadata_columns = self.state_parameters["metadata_added"]["columns"]
        lines.append(
            "Metadata columns: " + ", ".join(metadata_columns)
            if metadata_columns
            else "Metadata: not added"
        )
        lines.append(
            "Node Pgen: calculated"
            if self.state["node_pgen_calculated"]
            else "Node Pgen: not calculated"
        )
        lines.append(
            "ALICE: calculated"
            if self.state["alice_calculated"]
            else "ALICE: not calculated"
        )
        if self.state["alice_calculated"]:
            lines.append("ALICE parameters:")
            lines.extend(
                self._format_state_parameters(
                    self.state_parameters.get("alice_calculated")
                )
            )
        return "\n".join(lines)


    def __repr__(self):
        return self.__str__()


    def write_cluster_no_to_nodes(self):
        cluster_no = 0
        for cluster in self.clusters:
            for node in cluster.nodes():
                cluster.nodes[node]['cluster_no'] = cluster_no
            cluster_no += 1


    def add_metadata(self, metadata):
        self._require_clonotypes("add metadata")
        if not isinstance(metadata, pd.DataFrame):
            raise TypeError("metadata must be a pandas DataFrame.")
        if "sample_id" not in metadata.columns:
            raise ValueError("metadata must contain a 'sample_id' column.")
        metadata = metadata.copy()
        self._metadata_frames.append(metadata)
        if self._has_clusters():
            self._apply_metadata_to_clusters(self.clusters, metadata)
        columns = self.state_parameters["metadata_added"]["columns"]
        for column in self._metadata_columns(metadata):
            if column not in columns:
                columns.append(column)
        self.state["metadata_added"] = True
        self.state["empty"] = False


    @staticmethod
    def _cluster_number(cluster, fallback):
        if cluster.id is not None:
            return cluster.id
        for node in cluster:
            if "cluster_no" in node.additional_properties:
                return node.additional_properties["cluster_no"]
            if "cluster_no" in cluster.nodes[node]:
                return cluster.nodes[node]["cluster_no"]
            break
        return fallback


    def _copy_with_clusters(self, selected_clusters):
        result = Clusters()
        for attribute_name, value in self.__dict__.items():
            if attribute_name not in {
                "clusters",
                "state",
                "state_parameters",
                "_properties_cache",
                "_cluster_cache_signature",
            }:
                setattr(result, attribute_name, value)
        result.clusters = list(selected_clusters)
        result.state = copy.deepcopy(self.state)
        result.state_parameters = copy.deepcopy(self.state_parameters)
        result.state["clusters_created"] = True
        result.state["empty"] = False
        result._properties_cache = None
        result._cluster_cache_signature = None
        result._invalidate_cluster_dependent_analysis()
        return result


    def _resolve_cluster_identifiers(self, cluster_identifiers):
        if isinstance(cluster_identifiers, (str, int, np.integer)) and not isinstance(
            cluster_identifiers, bool
        ):
            identifiers = [cluster_identifiers]
        else:
            try:
                identifiers = list(cluster_identifiers)
            except TypeError as error:
                raise TypeError(
                    "cluster identifiers must be an identifier or iterable of "
                    "identifiers."
                ) from error

        clusters_by_number = {}
        for fallback_number, cluster in enumerate(self.clusters):
            cluster_no = self._cluster_number(cluster, fallback_number)
            if cluster_no in clusters_by_number:
                raise RuntimeError(f"Duplicate cluster_no found: {cluster_no}")
            clusters_by_number[cluster_no] = cluster

        requested_numbers = []
        invalid_identifiers = []
        for identifier in identifiers:
            if isinstance(identifier, (int, np.integer)) and not isinstance(
                identifier, bool
            ):
                cluster_no = int(identifier)
            elif isinstance(identifier, str):
                match = re.fullmatch(r"cluster_(\d+)", identifier)
                if match is None:
                    invalid_identifiers.append(identifier)
                    continue
                cluster_no = int(match.group(1))
            else:
                invalid_identifiers.append(identifier)
                continue
            if cluster_no not in requested_numbers:
                requested_numbers.append(cluster_no)

        missing_numbers = [
            cluster_no
            for cluster_no in requested_numbers
            if cluster_no not in clusters_by_number
        ]
        if invalid_identifiers or missing_numbers:
            missing = [repr(value) for value in invalid_identifiers]
            missing.extend(f"cluster_{number}" for number in missing_numbers)
            raise KeyError("Unknown cluster identifiers: " + ", ".join(missing))
        return [
            (cluster_no, f"cluster_{cluster_no}", clusters_by_number[cluster_no])
            for cluster_no in requested_numbers
        ]


    def select(self, cluster_identifiers):
        """Select clusters by integer cluster numbers or ``cluster_N`` IDs.

        Input order is preserved and duplicate identifiers are returned once.
        Scalars, ranges, pandas Series, and other iterables are accepted.
        """
        self._require_clusters("select clusters")
        resolved_clusters = self._resolve_cluster_identifiers(cluster_identifiers)
        return self._copy_with_clusters(
            [cluster for _, _, cluster in resolved_clusters]
        )


    def top_clusters(self, number_of_clusters):
        """Return the top-ranked clusters as a new collection."""
        return self.filter(top_clusters(number_of_clusters))


    def filter(self, condition, inplace=False):
        """Select clusters satisfying a composable cluster predicate.

        Args:
            condition (ClusterExpression | ClusterCollectionSelector):
                Boolean expression created with comparisons, ``all_nodes``,
                ``any_nodes``, and ``&``, ``|``, or ``~`` operators, or a
                collection selector such as ``top_clusters``.
            inplace (bool): Replace this collection when ``True``.

        Returns:
            Clusters: Filtered collection. Original cluster identifiers are
            preserved.
        """
        self._require_clusters("filter clusters")
        if not isinstance(inplace, (bool, np.bool_)):
            raise TypeError("inplace must be a boolean.")
        if isinstance(condition, ClusterCollectionSelector):
            selected_clusters = condition.select(self)
        else:
            if (
                not isinstance(condition, ClusterExpression)
                or not condition.is_boolean
            ):
                raise TypeError(
                    "condition must be a boolean cluster expression or "
                    "collection selector."
                )
            selected_clusters = [
                cluster for cluster in self.clusters if condition.evaluate(cluster)
            ]
        if inplace:
            self.clusters = selected_clusters
            self._cluster_cache_signature = None
            self._invalidate_properties_cache()
            self._invalidate_cluster_dependent_analysis()
            return None
        return self._copy_with_clusters(selected_clusters)


    def custom_properties(self, expressions):
        """Calculate custom cluster-level expressions as a dataframe.

        ``cluster_no`` and ``cluster_id`` are always included before the
        requested expression columns. Use ``expression.alias(name)`` to set a
        custom output column name.
        """
        self._require_clusters("calculate custom cluster properties")
        if isinstance(expressions, ClusterExpression):
            expressions = [expressions]
        else:
            try:
                expressions = list(expressions)
            except TypeError as error:
                raise TypeError(
                    "expressions must be a cluster expression or iterable of "
                    "cluster expressions."
                ) from error
        if not expressions:
            raise ValueError("At least one cluster expression must be provided.")
        if any(
            not isinstance(expression, ClusterExpression)
            for expression in expressions
        ):
            raise TypeError("Every custom property must be a cluster expression.")

        expression_names = [expression.name for expression in expressions]
        reserved_names = {"cluster_no", "cluster_id"}
        duplicate_names = {
            name for name in expression_names if expression_names.count(name) > 1
        }
        invalid_names = reserved_names.intersection(expression_names)
        if duplicate_names or invalid_names:
            names = sorted(duplicate_names.union(invalid_names))
            raise ValueError(
                "Custom property names must be unique and cannot use reserved "
                f"columns: {', '.join(names)}"
            )

        rows = []
        for fallback_number, cluster in enumerate(self.clusters):
            cluster_no = self._cluster_number(cluster, fallback_number)
            row = {
                "cluster_no": cluster_no,
                "cluster_id": f"cluster_{cluster_no}",
            }
            for expression in expressions:
                row[expression.name] = expression.evaluate(cluster)
            rows.append(row)
        return pd.DataFrame(
            rows,
            columns=["cluster_no", "cluster_id", *expression_names],
        )


    @staticmethod
    def _plot_node_property(node, property_name):
        value = _node_property_value(node, property_name)

        if value is None:
            return "NA"
        missing = pd.isna(value)
        if isinstance(missing, (bool, np.bool_)) and missing:
            return "NA"
        try:
            hash(value)
        except TypeError:
            return str(value)
        return value


    @staticmethod
    def _plot_color_map(levels, palette):
        if isinstance(palette, dict):
            missing_levels = [level for level in levels if level not in palette]
            if missing_levels:
                missing_text = ", ".join(str(level) for level in missing_levels)
                raise ValueError(
                    f"Palette does not define colors for: {missing_text}"
                )
            color_map = {level: palette[level] for level in levels}
        else:
            if palette is None:
                palette = "tab10" if len(levels) <= 10 else "husl"
            if isinstance(palette, str):
                colors = sns.color_palette(palette, n_colors=len(levels))
            else:
                colors = list(palette)
                if len(colors) < len(levels):
                    raise ValueError(
                        "Palette must contain at least as many colors as color levels."
                    )
            color_map = dict(zip(levels, colors))

        for color_value in color_map.values():
            to_rgba(color_value)
        return color_map


    @staticmethod
    def _plot_cluster_layout(cluster, layout, seed):
        layout_name = layout.lower().replace("-", "_")
        if layout_name in {"spring", "fruchterman_reingold"}:
            return nx.spring_layout(cluster, seed=seed)
        if layout_name in {"kamada_kawai", "kk"}:
            return nx.kamada_kawai_layout(cluster)
        if layout_name == "circular":
            return nx.circular_layout(cluster)
        if layout_name == "shell":
            return nx.shell_layout(cluster)
        if layout_name == "spectral":
            return nx.spectral_layout(cluster)
        raise ValueError(
            "Unknown layout. Possible values: spring, kamada_kawai, circular, "
            "shell, spectral."
        )


    def plot_logo(
        self,
        cluster_no,
        seq_type="prot",
        weight="count",
        plot=True,
    ):
        """Create weighted sequence logos selected by cluster number or ID.

        ``cluster_no`` accepts an integer cluster number, a ``cluster_N`` ID,
        or an iterable mixing both forms. For multiple clusters, ``plot=False``
        returns a dictionary keyed by cluster ID; ``plot=True`` draws each logo
        and returns ``None``.
        """
        self._require_clusters("plot a cluster logo")
        resolved_clusters = self._resolve_cluster_identifiers(cluster_no)
        if not resolved_clusters:
            raise ValueError("At least one cluster identifier must be provided.")
        if seq_type not in {"prot", "dna"}:
            raise ValueError("seq_type must be either 'prot' or 'dna'.")
        if not isinstance(weight, str) or not weight:
            raise TypeError("weight must be a non-empty node property name.")
        if not isinstance(plot, (bool, np.bool_)):
            raise TypeError("plot must be a boolean.")

        results = {}
        for _, cluster_id, cluster in resolved_clusters:
            clonotypes = []
            for node in cluster:
                sequence = node.seq_aa if seq_type == "prot" else node.seq_nt
                if seq_type == "dna" and sequence == "-":
                    raise ValueError(
                        f"Node '{node.id}' does not have a cdr3nt sequence."
                    )
                clonotypes.append((sequence, _node_weight(node, weight)))
            results[cluster_id] = get_logo_for_list_of_clonotypes(
                clonotypes, seq_type, plot=bool(plot)
            )

        if len(results) == 1:
            return next(iter(results.values()))
        if plot:
            return None
        return results


    def plot_cluster(
        self,
        cluster_no,
        layout="spring",
        color=None,
        palette=None,
        label=None,
        shape=None,
        ncols=None,
        figsize=None,
        seed=1,
        size="count",
        min_size=50,
        log_scaled=False,
        linear_scale=1,
        max_clusters=50,
        height=4,
        aspect=1,
    ):
        """Plot one, selected, or all clusters as network facets.

        Node area is uniform when ``size=None``. Otherwise, ``size`` names a
        numeric :class:`Node` attribute or ``node.additional_properties`` value.
        Log scaling uses ``min_size + linear_scale * log2(value + 1)``; linear
        scaling uses ``min_size + linear_scale * value``.

        Args:
            cluster_no (int | str | iterable | None): Cluster numbers,
                ``cluster_N`` IDs, mixed iterables, or ``None`` for all clusters.
            layout (str): ``spring`` (default), ``kamada_kawai``, ``circular``,
                ``shell``, or ``spectral``.
            color (str | None): Node property used for color grouping.
            palette (dict | list | str | None): Custom level-to-color mapping,
                color sequence, or seaborn palette name.
            label (str | None): Node property displayed in the node center.
            shape (str | None): Node property used for shape grouping. At most
                five levels are supported.
            ncols (int | None): Number of facet columns.
            figsize (tuple | None): Matplotlib figure size.
            seed (int): Seed used by the spring layout.
            size (str | None): Numeric node property used for node area. Use
                ``None`` for uniform node sizes.
            min_size (float): Minimum matplotlib node area.
            log_scaled (bool): Apply logarithmic scaling when ``True``.
            linear_scale (float): Multiplier for linear or log2-transformed
                values.
            max_clusters (int): Maximum number of clusters automatically plotted
                when ``cluster_no=None``. Ignored for explicit identifiers.
            height (float): Height of each facet in inches.
            aspect (float): Facet width divided by facet height.

        Returns:
            matplotlib.figure.Figure: The generated figure.
        """
        self._require_clusters("plot clusters")
        plot_all_clusters = cluster_no is None
        if plot_all_clusters:
            if (
                not isinstance(max_clusters, (int, np.integer))
                or isinstance(max_clusters, bool)
                or max_clusters < 1
            ):
                raise ValueError(
                    "max_clusters must be a positive integer when "
                    "cluster_no=None."
                )
            cluster_count = len(self.clusters)
            if cluster_count > max_clusters:
                raise ValueError(
                    f"plot_cluster(cluster_no=None) would plot all "
                    f"{cluster_count} clusters, which exceeds "
                    f"max_clusters={max_clusters}. No plot was created. "
                    "Increase max_clusters if plotting every cluster is "
                    "intentional, or pass cluster numbers/IDs for a smaller "
                    "selection."
                )
            identifiers = [
                self._cluster_number(cluster, fallback_number)
                for fallback_number, cluster in enumerate(self.clusters)
            ]
        else:
            if isinstance(cluster_no, (str, int, np.integer)) and not isinstance(
                cluster_no, bool
            ):
                identifiers = cluster_no
            else:
                try:
                    identifiers = list(cluster_no)
                except TypeError as error:
                    raise TypeError(
                        "cluster_no must be None, a cluster number/ID, or an "
                        "iterable of cluster numbers/IDs."
                    ) from error
                if len(identifiers) > 50:
                    raise ValueError("At most 50 clusters can be plotted at once.")

        resolved_clusters = self._resolve_cluster_identifiers(identifiers)
        if not resolved_clusters:
            raise ValueError("At least one cluster identifier must be provided.")
        cluster_numbers = [number for number, _, _ in resolved_clusters]
        cluster_ids = [cluster_id for _, cluster_id, _ in resolved_clusters]
        selected_clusters = [cluster for _, _, cluster in resolved_clusters]
        if palette is not None and color is None:
            raise ValueError("palette requires a color property.")
        if ncols is not None and (not isinstance(ncols, int) or ncols < 1):
            raise ValueError("ncols must be a positive integer.")
        if min_size < 0 or linear_scale < 0:
            raise ValueError("Node size parameters must be non-negative.")
        if (
            not isinstance(height, Real)
            or isinstance(height, (bool, np.bool_))
            or height <= 0
        ):
            raise ValueError("height must be a positive number.")
        if (
            not isinstance(aspect, Real)
            or isinstance(aspect, (bool, np.bool_))
            or aspect <= 0
        ):
            raise ValueError("aspect must be a positive number.")
        if not isinstance(log_scaled, (bool, np.bool_)):
            raise TypeError("log_scaled must be a boolean.")

        selected_nodes = [
            node for cluster in selected_clusters for node in cluster.nodes
        ]

        color_values = {}
        color_levels = []
        if color is not None:
            for node in selected_nodes:
                value = self._plot_node_property(node, color)
                color_values[node] = value
                if value not in color_levels:
                    color_levels.append(value)
            if color == "isotype":
                color_levels = _isotype_order(color_levels)
                if palette is None:
                    palette = _isotype_colors(color_levels)
            color_map = self._plot_color_map(color_levels, palette)
        else:
            color_map = {None: "#4C78A8"}

        shape_markers = ["o", "^", "D", "h", "s"]
        shape_values = {}
        shape_levels = []
        if shape is not None:
            for node in selected_nodes:
                value = self._plot_node_property(node, shape)
                shape_values[node] = value
                if value not in shape_levels:
                    shape_levels.append(value)
            if len(shape_levels) > len(shape_markers):
                raise ValueError(
                    "shape supports at most five levels: circle, triangle, "
                    "rhombus, hexagon, and square."
                )
            shape_map = dict(zip(shape_levels, shape_markers))
        else:
            shape_map = {None: shape_markers[0]}

        node_sizes = {}
        for node in selected_nodes:
            if size is None:
                node_sizes[node] = min_size
                continue
            size_value = self._plot_node_property(node, size)
            try:
                size_value = float(size_value)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Node '{node.id}' property '{size}' must be numeric."
                ) from error
            if not np.isfinite(size_value) or size_value < 0:
                raise ValueError(
                    f"Node '{node.id}' property '{size}' must be finite and "
                    "non-negative."
                )
            if log_scaled:
                node_sizes[node] = (
                    min_size + linear_scale * np.log2(size_value + 1)
                )
            else:
                node_sizes[node] = min_size + linear_scale * size_value

        facet_count = len(selected_clusters)
        if ncols is None:
            ncols = min(5, math.ceil(math.sqrt(facet_count)))
        ncols = min(ncols, facet_count)
        nrows = math.ceil(facet_count / ncols)
        if figsize is None:
            figsize = (float(height) * float(aspect) * ncols, float(height) * nrows)
        figure, axes = plt.subplots(
            nrows, ncols, figsize=figsize, squeeze=False
        )
        axes_flat = axes.ravel()

        for axis, cluster_id, cluster in zip(
            axes_flat, cluster_ids, selected_clusters
        ):
            positions = self._plot_cluster_layout(cluster, layout, seed)
            nx.draw_networkx_edges(
                cluster,
                positions,
                ax=axis,
                edge_color="#A9A9A9",
                width=0.6,
                alpha=0.75,
            )
            levels_to_draw = shape_levels if shape is not None else [None]
            for shape_level in levels_to_draw:
                nodes_for_shape = [
                    node
                    for node in cluster.nodes
                    if shape is None or shape_values[node] == shape_level
                ]
                node_colors = [
                    color_map[color_values[node]] if color is not None
                    else color_map[None]
                    for node in nodes_for_shape
                ]
                nx.draw_networkx_nodes(
                    cluster,
                    positions,
                    nodelist=nodes_for_shape,
                    node_size=[node_sizes[node] for node in nodes_for_shape],
                    node_color=node_colors,
                    node_shape=shape_map[shape_level],
                    edgecolors="white",
                    linewidths=0.8,
                    ax=axis,
                )
            if label is not None:
                for node, position in positions.items():
                    label_value = self._plot_node_property(node, label)
                    axis.text(
                        position[0],
                        position[1],
                        str(label_value),
                        ha="center",
                        va="center",
                        fontsize=8,
                        zorder=4,
                    )
            axis.set_title(cluster_id)
            axis.set_axis_off()

        for axis in axes_flat[facet_count:]:
            axis.set_visible(False)

        legend_rows = 0
        if color is not None:
            color_handles = [
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="",
                    markerfacecolor=color_map[level],
                    markeredgecolor="white",
                    label=str(level),
                    markersize=8,
                )
                for level in color_levels
            ]
            figure.legend(
                handles=color_handles,
                title=color,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.01),
                ncol=min(8, len(color_handles)),
                frameon=False,
            )
            legend_rows += 1
        if shape is not None:
            shape_handles = [
                Line2D(
                    [],
                    [],
                    marker=shape_map[level],
                    linestyle="",
                    markerfacecolor="#808080",
                    markeredgecolor="white",
                    label=str(level),
                    markersize=8,
                )
                for level in shape_levels
            ]
            shape_y = 0.08 if color is not None else 0.01
            figure.legend(
                handles=shape_handles,
                title=shape,
                loc="lower center",
                bbox_to_anchor=(0.5, shape_y),
                ncol=min(5, len(shape_handles)),
                frameon=False,
            )
            legend_rows += 1

        bottom_margin = 0.04 + 0.08 * legend_rows
        figure.tight_layout(rect=(0, bottom_margin, 1, 1))
        plt.close(figure)
        return figure


    def properties(self, cpu=None):
        """Return cached basic cluster properties calculated in parallel.

        Args:
            cpu (int | None): Number of worker processes. Use ``1`` for a
                sequential calculation; ``None`` uses the executor default.
        """
        self._require_clusters("calculate cluster properties")
        if self._properties_cache is None:
            overlap_type = self.overlap_type
            if overlap_type is None:
                cluster_parameters = self.state_parameters.get(
                    "clusters_created"
                )
                if cluster_parameters:
                    overlap_type = cluster_parameters.get("overlap_type")
            if overlap_type is None:
                check_v = False
                check_j = False
            else:
                _, check_v, check_j = overlap_type_to_flags(overlap_type)

            tasks = [
                (
                    self._cluster_number(cluster, fallback_number),
                    cluster,
                    check_v,
                    check_j,
                )
                for fallback_number, cluster in enumerate(self.clusters)
            ]
            results = run_parallel_calculation(
                _cluster_properties_worker,
                tasks,
                "Calculating cluster properties",
                object_name="clusters",
                cpu=cpu,
            )
            self._properties_cache = pd.DataFrame(
                results, columns=_CLUSTER_PROPERTY_COLUMNS
            )
        return self._properties_cache.copy()


    def to_count_table(self, by_freq=False):
        """Create a wide table of cluster abundance by sample.

        Args:
            by_freq (bool): If ``True``, sum node frequencies. Otherwise, sum
                node counts.

        Returns:
            pd.DataFrame: Cluster identifiers and consensus properties followed
                by one abundance column per sample.
        """
        self._require_clusters("create a cluster count table")
        property_columns = [
            "cluster_id",
            "concensus_cdr3aa",
            "concensus_v",
            "concensus_j",
        ]
        count_table = self.properties()[property_columns].copy()
        count_table.insert(
            1,
            "consensus",
            count_table[["concensus_cdr3aa", "concensus_v", "concensus_j"]]
            .astype(str)
            .agg("|".join, axis=1),
        )

        sample_ids = []
        clonosets_df = getattr(self, "clonosets_df", None)
        if isinstance(clonosets_df, pd.DataFrame) and "sample_id" in clonosets_df.columns:
            sample_ids.extend(clonosets_df["sample_id"].drop_duplicates().tolist())
        elif isinstance(self.clonotypes, pd.DataFrame) and "sample_id" in self.clonotypes.columns:
            sample_ids.extend(self.clonotypes["sample_id"].drop_duplicates().tolist())

        for cluster in self.clusters:
            for node in cluster:
                if node.sample_id not in sample_ids:
                    sample_ids.append(node.sample_id)

        values_by_sample = {sample_id: [] for sample_id in sample_ids}
        for cluster in self.clusters:
            cluster_values = {sample_id: 0 for sample_id in sample_ids}
            for node in cluster:
                value = node.freq if by_freq else node.count
                cluster_values[node.sample_id] += value
            for sample_id in sample_ids:
                values_by_sample[sample_id].append(cluster_values[sample_id])

        for sample_id, values in values_by_sample.items():
            count_table[sample_id] = values

        return count_table


    def intersect_with_clonosets(
        self,
        clonosets_df,
        cl_filter=None,
        overlap_type="aaVJ",
        by_freq=False,
        mismatches=1,
        cpu=None,
    ):
        """Measure target-sample clonotypes similar to each cluster.

        Each target clonotype contributes at most once to a cluster, even when
        it is similar to several nodes in that cluster. Frequency values use
        the directional similarity definition: matched target counts divided
        by the filtered target clonoset's total count.
        """
        from .intersections import prepare_clonoset_for_intersection

        self._require_clusters("intersect clusters with clonosets")
        if not isinstance(clonosets_df, pd.DataFrame):
            raise TypeError("clonosets_df must be a pandas DataFrame.")
        required_columns = {"sample_id", "filename"}
        missing_columns = required_columns.difference(clonosets_df.columns)
        if missing_columns:
            raise ValueError(
                "clonosets_df is missing columns: "
                + ", ".join(sorted(missing_columns))
            )
        if clonosets_df["sample_id"].duplicated().any():
            raise ValueError("clonosets_df sample_id values must be unique.")
        if not isinstance(by_freq, (bool, np.bool_)):
            raise TypeError("by_freq must be a boolean.")
        overlap_type_to_flags(overlap_type)
        if (
            not isinstance(mismatches, (int, np.integer))
            or isinstance(mismatches, bool)
        ):
            raise TypeError("mismatches must be a non-negative integer.")
        if mismatches < 0:
            raise ValueError("mismatches must be a non-negative integer.")
        uses_sequence = overlap_type_uses_sequence(overlap_type)
        effective_mismatches = int(mismatches) if uses_sequence else 0

        cluster_dicts = {}
        for fallback_number, cluster in enumerate(self.clusters):
            cluster_no = self._cluster_number(cluster, fallback_number)
            cluster_clonotypes = pd.DataFrame(
                [
                    {
                        "cdr3aa": node.seq_aa,
                        "cdr3nt": node.seq_nt,
                        "v": node.v,
                        "j": node.j,
                        "count": 1,
                        "freq": 1,
                    }
                    for node in cluster
                ]
            )
            cluster_dicts[cluster_no] = prepare_clonoset_for_intersection(
                cluster_clonotypes,
                overlap_type=overlap_type,
                by_freq=False,
                len_vj_format=uses_sequence,
            )

        if cl_filter is None:
            cl_filter = Filter()
        tasks = [
            (
                row["sample_id"],
                row["filename"],
                cl_filter.spawn(),
                cluster_dicts,
                overlap_type,
                effective_mismatches,
                bool(by_freq),
            )
            for _, row in clonosets_df.iterrows()
        ]
        sample_results = run_parallel_calculation(
            _intersect_clusters_with_clonoset_worker,
            tasks,
            "Intersecting clusters with clonosets",
            object_name="samples",
            verbose=False,
            cpu=cpu,
        )
        values_by_sample = dict(sample_results)

        cluster_rows = []
        cluster_numbers = []
        for fallback_number, cluster in enumerate(self.clusters):
            cluster_no = self._cluster_number(cluster, fallback_number)
            cluster_numbers.append(cluster_no)
            aa_consensus = cluster.calc_cluster_consensus(
                seq_type="prot", weigh_by=None
            )
            v_consensus = cluster.calc_cluster_consensus_segment(
                segment_type="v", weigh_by=None
            )
            j_consensus = cluster.calc_cluster_consensus_segment(
                segment_type="j", weigh_by=None
            )
            cluster_rows.append(
                {
                    "cluster_id": f"cluster_{cluster_no}",
                    "consensus": (
                        f"{aa_consensus}|{v_consensus}|{j_consensus}"
                    ),
                    "concensus_cdr3aa": aa_consensus,
                    "concensus_v": v_consensus,
                    "concensus_j": j_consensus,
                }
            )
        count_table = pd.DataFrame(
            cluster_rows,
            columns=[
                "cluster_id",
                "consensus",
                "concensus_cdr3aa",
                "concensus_v",
                "concensus_j",
            ],
        )
        for sample_id in clonosets_df["sample_id"]:
            sample_values = values_by_sample.get(sample_id, {})
            count_table[sample_id] = [
                sample_values.get(cluster_no, 0.0)
                for cluster_no in cluster_numbers
            ]
        return count_table


    def split(self, method="leiden", resolution=0.5, threshold=1e-07, seed=1):
        """
        Performs a community detection on pre-calculated clusters. Available methods are `louvain` and `leiden`. 
        The latter uses `leidenalg` implementation. The result is saved as either `cluster_communities_leiden` or `cluster_communities_louvain`.

        Args:
            method (str): method for detecting communities. Possible options are `leiden` and `louvain`.
            resolution (float): Default is 0.5
            threshold (float): Default is 1e-07.
            seed (int): default is 1.
        
        Returns:
            None
        """
        self._require_clusters("split clusters into communities")
        if method == 'louvain':
            self.find_cluster_communities_louvain(resolution=resolution, 
                                                  threshold=threshold, 
                                                  seed=seed)
        elif method == 'leiden':
            self.find_cluster_communities_leiden(resolution=resolution, 
                                                  seed=seed)
        else:
            raise ValueError(f'Unknown method: {method}')


    def check_compulsory_columns(self, clonoset, compulsory_columns):
        return all(c in clonoset.columns for c in compulsory_columns)


    def set_pooled(self, pooled: bool):
        self.is_pooled = pooled


    @staticmethod
    def _resolve_verbosity(verbose=True, verbosity=None):
        if verbosity is not None:
            if not isinstance(verbosity, (bool, np.bool_)):
                raise TypeError("verbosity must be a boolean.")
            return bool(verbosity)
        if not isinstance(verbose, (bool, np.bool_)):
            raise TypeError("verbose must be a boolean.")
        return bool(verbose)


    def read_from_clonosets_df(
        self,
        clonosets_df: "pd.DataFrame",
        cl_filter=Filter(),
        verbose=True,
        verbosity=None,
    ):
        """Read and pool clonotypes described by a sample/file dataframe."""
        verbose = self._resolve_verbosity(verbose, verbosity)
        if not isinstance(clonosets_df, pd.DataFrame):
            raise TypeError("clonosets_df must be a pandas DataFrame.")
        required_columns = {"sample_id", "filename"}
        missing_columns = required_columns.difference(clonosets_df.columns)
        if missing_columns:
            raise ValueError(
                "clonosets_df is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        self.cl_filter = cl_filter
        self.clonosets_df = clonosets_df.copy()
        clonotypes_dfs = []
        for _, row in clonosets_df.iterrows():
            sample_id = row["sample_id"]
            clonoset = read_clonoset(row["filename"])
            clonoset = cl_filter.apply(clonoset)
            clonoset["sample_id"] = sample_id
            clonotypes_dfs.append(clonoset)

        self.clonotypes = pd.concat(clonotypes_dfs).reset_index(drop=True)
        self.set_pooled(False)
        self._mark_clonotypes_read(
            "clonosets_df",
            {
                "input_samples": len(clonosets_df),
                "filter": self._filter_parameters(cl_filter),
            },
        )
        if verbose:
            print(
                f"Pooled {len(self.clonotypes)} clonotypes from "
                f"{len(clonosets_df)} samples"
            )


    def read_from_pooled_clonoset(
        self, pooled_clonoset: "pd.DataFrame", verbose=True, verbosity=None
    ):
        """Read clonotypes from an already pooled dataframe."""
        verbose = self._resolve_verbosity(verbose, verbosity)
        if not isinstance(pooled_clonoset, pd.DataFrame):
            raise TypeError("pooled_clonoset must be a pandas DataFrame.")
        compulsory_columns = [
            "freq",
            "count",
            "v",
            "j",
            "cdr3aa",
            "cdr3nt",
            "sample_id",
        ]
        self.set_pooled(True)
        self.clonosets_df = None
        self.cl_filter = None
        converted = False
        if self.check_compulsory_columns(pooled_clonoset, compulsory_columns):
            self.clonotypes = pooled_clonoset.copy()
        else:
            if verbose:
                print("Trying to convert pooled clonoset...")
            pooled_df = Filter(
                by_umi=True, convert=False, recount_fractions=False
            ).apply(pooled_clonoset)
            if not self.check_compulsory_columns(pooled_df, compulsory_columns):
                error_message = (
                    "Couldn't find at least one of compulsory columns in "
                    "pooled_df: "
                    + ", ".join(compulsory_columns)
                )
                if "sample_id" not in pooled_clonoset.columns:
                    error_message += (
                        "\nNo `sample_id` column. Ensure this column is present "
                        "or add it manually before proceeding."
                    )
                raise ValueError(error_message)
            self.clonotypes = pooled_df
            converted = True
        self._mark_clonotypes_read(
            "pooled_clonoset", {"converted": converted}
        )


    @staticmethod
    def _canonicalize_edge_nodes(edges, nodes):
        """Replace deserialized edge nodes with the original node instances."""
        nodes_by_id = {}
        for node in nodes:
            if node.id in nodes_by_id:
                raise RuntimeError(f"Duplicate node_id generated: {node.id}")
            nodes_by_id[node.id] = node

        canonical_edges = []
        for node_1, node_2, edge_value in edges:
            try:
                canonical_node_1 = nodes_by_id[node_1.id]
                canonical_node_2 = nodes_by_id[node_2.id]
            except KeyError as error:
                raise RuntimeError(
                    f"Edge references unknown node_id: {error.args[0]}"
                ) from error
            canonical_edges.append(
                (canonical_node_1, canonical_node_2, edge_value)
            )
        return canonical_edges


    def find_nodes_and_edges(self, mismatches, overlap_type, cpu=None, verbose=True):
        """
        Builds nodes from clonotypes and finds edges between similar sequences. 
        Parallel calculations are implemented.

        Args:
            self
            mismatches: Number of allowed mismatches for clustering.
            overlap_type: Rules for sequence comparison (e.g., "aaV", "ntVJ").

        Returns:
            tuple: (single nodes, edges)
        """
        clonoset = self.clonotypes
        aa, check_v, check_j = overlap_type_to_flags(overlap_type)
        clonoset = Filter(by_umi=True).apply(clonoset)
        
        igh = "c" in clonoset.columns
        is_freq = "freq" in clonoset.columns
        is_count = "count" in clonoset.columns
        if igh:
            clonoset["isotype"] = clonoset["c"].apply(_recode_cluster_isotype)
        # if is_count:
        #     clonoset["count"] = 1

        nodes_by_comparison_group = {}
        list_of_all_nodes = []

        
        for node_id, (_, row) in enumerate(clonoset.iterrows()):
            v = row["v"]
            j = row["j"]
            cdr3aa = row["cdr3aa"]
            if "cdr3nt" in clonoset.columns:
                cdr3nt = row["cdr3nt"]
            else:
                cdr3nt = "-"
            if is_freq:
                freq = row["freq"]
            else:
                freq = 1 / len(clonoset)
            count = row["count"]
            sample_id = row["sample_id"]
            len_cdr3aa = len(cdr3aa)
            comparison_group = [len_cdr3aa]
            if check_v:
                comparison_group.append(v)
            if check_j:
                comparison_group.append(j)
            comparison_group = tuple(comparison_group)
            if comparison_group not in nodes_by_comparison_group:
                nodes_by_comparison_group[comparison_group] = []
            node = Node(node_id=node_id,
                        seq_nt=cdr3nt, 
                        seq_aa=cdr3aa, 
                        v=v, 
                        j=j, 
                        sample_id=sample_id, 
                        freq=freq,
                        count=count)
            if igh:
                isotype = row["isotype"]
                node.additional_properties["isotype"] = (
                    None if pd.isna(isotype) else isotype
                )
            nodes_by_comparison_group[comparison_group].append(node)
            list_of_all_nodes.append(node)
        if verbose:
            print("Nodes list created: {} nodes".format(len(list_of_all_nodes)))
        
        tasks = []
        for nodes_list in nodes_by_comparison_group.values():
            task = (nodes_list, mismatches, aa, check_v, check_j)
            tasks.append(task)
        
        program_name = "Find neighbour clonotypes"
        result_list = run_parallel_calculation(
            self.find_edges_in_nodes_set_mp,
            tasks,
            program_name,
            verbose=verbose,
            cpu=cpu,
        )
        edges = [edge for result in result_list for edge in result]
        edges = self._canonicalize_edge_nodes(edges, list_of_all_nodes)
        if verbose:
            print("found {} edges".format(len(edges)))
        nodes_set = set(list_of_all_nodes)
        connected_nodes_set = set()
        for edge in edges:
            connected_nodes_set.add(edge[0])
            connected_nodes_set.add(edge[1])
        single_nodes_set = nodes_set.difference(connected_nodes_set)
        
        return list(single_nodes_set), edges


    def find_edges_in_nodes_set_mp(self, args):
        (nodes_list, mismatches, aa, check_v, check_j) = args
        edges = []
        
        for i in range(len(nodes_list)-1):
            for j in range(i+1,len(nodes_list)):
                node_1=nodes_list[i]
                node_2=nodes_list[j]
                if node_1.is_neighbour_of(node_2, mismatches, aa, check_v, check_j):
                    edges.append((node_1, node_2, 1)) #save unique codes
        return edges


    def find_nodes_and_edges_tcrdist_no_gaps(self, radius=16, cpu=None, verbose=True):
        
        with open(os.path.join(os.path.dirname(__file__), 'tcrdist_ab_v_segments.json')) as f:
            TCRDIST_V_DIST = json.load(f)

        rename_some_v_segments_dict = {"TRAV14/DV4": "TRAV14DV4",
                                    "TRAV23/DV6": "TRAV23DV6",
                                    "TRAV29/DV5": "TRAV29DV5",
                                    "TRAV36/DV7": "TRAV36DV7",
                                    "TRAV38-2/DV8": "TRAV38-2DV8"}
        
        # if isinstance(self.clonotypes, str):
        #     clonoset = pd.read_csv(self.clonotypes,sep="\t")
        # else:
        #     clonoset = self.clonotypes
        clonoset = self.clonotypes
        clonoset = clonoset.rename(columns={"bestVGene": "v",
                                            "bestJGene": "j",
                                            "CDR3.amino.acid.sequence": "cdr3aa",
                                            "CDR3.nucleotide.sequence": "cdr3nt",
                                            "allVHitsWithScore": "v",
                                            "allJHitsWithScore": "j",
                                            "allCHitsWithScore": "c",
                                            "aaSeqCDR3": "cdr3aa",
                                            "nSeqCDR3": "cdr3nt",
                                            "Sample":"sample_id",
                                            "cloneFraction":"freq",
                                            "Read.count": "count",
                                            "cloneCount": "count",
                                            "uniqueUMICount":"uniqueMoleculeCount",
                                            "uniqueUMIFraction":"uniqueMoleculeFraction"})
        
        count_column = "count"
        fraction_column = "freq"
        if "uniqueMoleculeCount" in clonoset.columns:
            count_column = "uniqueMoleculeCount"
            fraction_column = "uniqueMoleculeFraction"
        elif "readCount" in clonoset.columns:
            count_column = "readCount"
            fraction_column = "readFraction"
        clonoset = clonoset.rename(columns={count_column: "count", fraction_column: "freq"})
        
        clonoset["v"] = clonoset["v"].apply(lambda x: x.split("*")[0])
        clonoset["v"] = clonoset["v"].apply(lambda x: rename_some_v_segments_dict[x] if x in rename_some_v_segments_dict else x)
        clonoset["j"] = clonoset["j"].apply(lambda x: x.split("*")[0])

        igh = "c" in clonoset.columns 

        if igh:
            clonoset["isotype"] = clonoset["c"].apply(_recode_cluster_isotype)


        nodes_by_len = {}
        list_of_all_nodes = []

        is_freq = "freq" in clonoset.columns
        is_count = "count" in clonoset.columns
        if not is_count:
            clonoset["count"] = 1

        for node_id, (_, row) in enumerate(clonoset.iterrows()):
            v = row["v"]
            if v not in TCRDIST_V_DIST:
                if verbose:
                    print(
                        f"Warning! {v} gene was not recognized in reference "
                        "db; no CDR sequence could be inferred. The clone was "
                        "skipped"
                    )
                continue
            j = row["j"]
            cdr3aa = row["cdr3aa"]
            if "cdr3nt" in clonoset.columns:
                cdr3nt = row["cdr3nt"]
            else:
                cdr3nt = "-"
            if is_freq:
                freq = row["freq"]
            else:
                freq = 1 / len(clonoset)
            count = row["count"]
            sample_id = row["sample_id"]
            len_cdr3aa = len(cdr3aa)
            if len_cdr3aa not in nodes_by_len:
                nodes_by_len[len_cdr3aa] = []
            node = Node(node_id, cdr3nt, cdr3aa, v, j, sample_id, freq=freq, count=count)
            if igh:
                isotype = row["isotype"]
                node.additional_properties["isotype"] = (
                    None if pd.isna(isotype) else isotype
                )
            nodes_by_len[len_cdr3aa].append(node)
            list_of_all_nodes.append(node)
        if verbose:
            print("Nodes list created: {} nodes".format(len(list_of_all_nodes)))

        tasks = []
        for aa_len in nodes_by_len:
            nodes_list=nodes_by_len[aa_len]
            task = (nodes_list, radius)
            tasks.append(task)
        
        program_name = "Find neighbour clonotypes (TCRdist with no CDR3 gaps)"
        result_list = run_parallel_calculation(
            self.find_nodes_and_edges_tcrdist_no_gaps_mp,
            tasks,
            program_name,
            verbose=verbose,
            cpu=cpu,
        )
        edges = [edge for result in result_list for edge in result]
        edges = self._canonicalize_edge_nodes(edges, list_of_all_nodes)
        if verbose:
            print("found {} edges".format(len(edges)))
        nodes_set = set(list_of_all_nodes)
        connected_nodes_set = set()
        for edge in edges:
            connected_nodes_set.add(edge[0])
            connected_nodes_set.add(edge[1])
        single_nodes_set = nodes_set.difference(connected_nodes_set)
        
        return list(single_nodes_set), edges


    def find_nodes_and_edges_tcrdist_no_gaps_mp(self, args):
        (nodes_list, radius) = args
        edges = []
        
        for i in range(len(nodes_list)-1):
            for j in range(i+1,len(nodes_list)):
                node_1=nodes_list[i]
                node_2=nodes_list[j]
                dist = self.tcr_dist(node_1, node_2, radius)
                if dist >= 0:
                    edges.append((node_1, node_2, dist)) #save unique codes
        return edges


    def tcr_dist(self, node_1, node_2, radius):

        dist=self.TCRDIST_V_DIST[node_1.v][node_2.v]
        if dist > radius:
            return -1
        dist += self.TCRDIST_CDR3_SCORE_MULTIPLIER*sum([self.TCRDIST_BLOSUM[(a,b)] for a,b in zip(node_1.seq_aa[self.TCRDIST_CDR3_N_CUT:-self.TCRDIST_CDR3_C_CUT],
                                                        node_2.seq_aa[self.TCRDIST_CDR3_N_CUT:-self.TCRDIST_CDR3_C_CUT])])
        if dist <= radius:
            return dist
        return -1
    

    def filter_one_node_clusters(self, inplace=False):
        if inplace:
            self.clusters = [c for c in self.clusters if len(c) > 1]
            self._cluster_cache_signature = None
            self._invalidate_properties_cache()
            self._invalidate_cluster_dependent_analysis()
        else:
            return [c for c in self.clusters if len(c) > 1]


# !!! add igh check inside the function - only if there is a c column 
    def create_clusters(self,
                        overlap_type='aaVJ',
                        mismatches=1,
                        tcrdist_radius=None,
                        cpu=None,
                        verbose=True,
                        verbosity=None):
        """
        Creates clusters of clonotypes using either mismatch-based or distance-based methods.

        If `tcrdist_radius` is provided, clustering uses TCRdist algorithm 
        without CDR3 gaps. Otherwise, clustering is based on sequence similarity using 
        the specified `overlap_type` and allowed number of `mismatches`.

        Args:
            self
            overlap_type (str): Sequence comparison rule (possible values are: "aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ").
            mismatches (int): Maximum allowed mismatches for clustering.
            tcrdist_radius (int, optional): TCRdist radius.

        Returns:
            list: List of clusters (nx.Graph objects).
        """
        verbose = self._resolve_verbosity(verbose, verbosity)
        self._require_clonotypes("create clusters")
        possible_overlap_types = ["aa", "aaV", "aaVJ", "nt", "ntV", "ntVJ", "VJ", "VJlen"]
        compulsory_columns = ["freq", "count", "v", "j", "cdr3aa", "cdr3nt", "sample_id"]
        tcr_dist = isinstance(tcrdist_radius, int)

        if tcrdist_radius is None: 

            if overlap_type not in possible_overlap_types:
                error_message = f'Incorrect overlap type  {overlap_type}. Possible values: {", ".join(possible_overlap_types)}'
                raise ValueError(error_message)

            if not isinstance(mismatches, int):
                error_message = f"Incorrect value for mismatches: {mismatches}. Expected a non-negative integer (e.g., 0, 1, 2)."
                raise ValueError(error_message)

            self.overlap_type = overlap_type
            self.mismatches = mismatches

        if tcr_dist:
            cluster_parameters = {
                "method": "tcrdist",
                "tcrdist_radius": tcrdist_radius,
            }
            cache_parameters = ("tcrdist", tcrdist_radius)
        else:
            cluster_parameters = {
                "method": "mismatches",
                "overlap_type": overlap_type,
                "mismatches": mismatches,
            }
            cache_parameters = ("mismatches", overlap_type, mismatches)
        cache_signature = (
            self._clonotypes_revision,
            id(self.clonotypes),
            len(self.clonotypes),
            cache_parameters,
        )
        if (
            self._cluster_cache_signature == cache_signature
            and self.state["clusters_created"]
        ):
            if verbose:
                summary = self._cluster_summary()
                print(
                    "Clusters were already created from the same clonotypes with "
                    "the same parameters. "
                    f"Graph: {summary['total_nodes']} nodes and "
                    f"{summary['total_edges']} edges. "
                    f"{summary['multi_node_clusters']} clusters (2 or more nodes) "
                    f"and {summary['single_node_clusters']} single nodes. Total: "
                    f"{summary['total_clusters']}"
                )
            return None

        if tcr_dist:
            self.tcrdist_radius = tcrdist_radius
            nodes, edges = self.find_nodes_and_edges_tcrdist_no_gaps(
                radius=tcrdist_radius, cpu=cpu, verbose=verbose
            )
        else:
            nodes, edges = self.find_nodes_and_edges(
                mismatches, overlap_type, cpu=cpu, verbose=verbose
            )
        
        graph_nodes_by_id = {}
        for node in [*nodes, *(node for edge in edges for node in edge[:2])]:
            existing_node = graph_nodes_by_id.get(node.id)
            if existing_node is not None and existing_node is not node:
                raise RuntimeError(
                    f"Multiple Node objects were produced for node_id {node.id}."
                )
            graph_nodes_by_id[node.id] = node

        main_graph = Cluster()
        main_graph.add_nodes_from(graph_nodes_by_id.values())
        if verbose:
            print("-----------------------------\nNexworkX graph created")

        program_name = "Adding edges..."
        edges_done = 0
        edges_total = len(edges)
        node_id_dict = {node.id: node for node in main_graph}
        for edge in edges:
            node1 = node_id_dict[edge[0].id]
            node2 = node_id_dict[edge[1].id]
            length = edge[2]
            main_graph.add_edge(node1, node2, length=length)
            edges_done += 1

        self.clusters = [main_graph.subgraph(c).copy() for c in nx.connected_components(main_graph)]
        total_clusters = len(self.clusters)
        cluster_num = len(self.filter_one_node_clusters(inplace=False))
        singletons = total_clusters - cluster_num
        
        if verbose:
            print(f"Found {cluster_num} clusters (2 or more nodes) and {singletons} single nodes. Total: {total_clusters}")
        
        self.clusters.sort(key=lambda x: (-len(x), x.calc_cluster_consensus(seq_type="prot", weigh_by=None)))
        self.write_cluster_no_to_nodes()

        for i, cluster in enumerate(self.clusters):
            cluster.id = i
            for j, node in enumerate(cluster):
                node.additional_properties["cluster_no"] = i
                node.additional_properties['n_neighbours'] = cluster.degree(node)


        self._apply_stored_metadata()
        self._cluster_cache_signature = cache_signature
        self._invalidate_properties_cache()
        self.alice_results = None
        self.state.update(
            {
                "empty": False,
                "clusters_created": True,
                "node_pgen_calculated": False,
                "alice_calculated": False,
            }
        )
        self.state_parameters["clusters_created"] = cluster_parameters
        self.state_parameters["node_pgen_calculated"] = None
        self.state_parameters["alice_calculated"] = None


# !!! add check_progress
# !!! add cluster_no before communities and visa versa
# !!! add wrapper for louvain and leiden
# !!! remove networkx representation from the user
    def find_cluster_communities_louvain(self, resolution=1, threshold=1e-07, seed=1):
        """
        Apply Louvain community detection to each cluster.
        
        Args:
            resolution (float): Resolution parameter for Louvain algorithm.
            threshold (float): Convergence threshold.
            seed (int): Random seed for reproducibility.
        
        Returns:
            List of NetworkX Graphs corresponding to detected communities.
        """

        self._require_clusters("calculate Louvain communities")
        total_communities = 0
        self.cluster_communities_louvain = ClusterCommunities()
        self.cluster_communities_louvain.resolution = resolution
        self.cluster_communities_louvain.threshold = threshold
        self.cluster_communities_louvain.seed = seed 

        for cluster in self.clusters:
            if len(cluster) < 2:
                for node in cluster:
                    node.additional_properties["community"] = total_communities
                total_communities += 1 
                self.cluster_communities_louvain.communities.append(cluster)
            else:
                cluster_communities = community.louvain_communities(cluster, resolution=resolution, threshold=threshold, seed=seed)
                for com in cluster_communities:
                    com_nodes = []
                    for node in com:
                        node.additional_properties["community"] = total_communities
                        com_nodes.append(node)
                    total_communities += 1  
                    self.cluster_communities_louvain.communities.append(cluster.subgraph(com_nodes))
        self.cluster_communities_louvain.communities.sort(key=lambda x: len(x), reverse=True)


    def find_cluster_communities_leiden(self, resolution=1, seed=1):
        """
        Apply Leiden community detection to each cluster.
        
        Args:
            resolution (float): Resolution parameter for Leiden algorithm.
            seed (int): Random seed for reproducibility.
        
        Returns:
            List of NetworkX Graphs corresponding to detected communities.
        """
        self._require_clusters("calculate Leiden communities")
        try:
            import igraph as ig
            import leidenalg
        except ImportError as exc:
            raise ImportError(
                "Leiden community detection requires optional clustering "
                "dependencies. Install them with `pip install repseq[clustering]`."
            ) from exc

        total_communities = 0
        self.cluster_communities_leiden = ClusterCommunities()
        self.cluster_communities_leiden.resolution = resolution
        self.cluster_communities_leiden.seed = seed 

        for cluster in self.clusters:
            if len(cluster) < 2:
                for node in cluster:
                    node.additional_properties["leiden_community"] = total_communities
                total_communities += 1 
                self.cluster_communities_leiden.communities.append(cluster)
            else:
                # leidenalg requires an igraph Graph, hence the conversion
                cluster_igraph = ig.Graph.from_networkx(cluster)
                nodes = cluster_igraph.vs["_nx_name"]

                partition = leidenalg.find_partition(
                    cluster_igraph,
                    leidenalg.RBConfigurationVertexPartition,
                    resolution_parameter=resolution,
                    seed=seed)

                for com in partition:
                    com_nodes = []
                    for i in com:
                        node = nodes[i]
                        node.additional_properties["leiden_community"] = total_communities
                        com_nodes.append(node)
                    total_communities += 1  
                    self.cluster_communities_leiden.communities.append(cluster.subgraph(com_nodes))

        self.cluster_communities_leiden.communities.sort(key=lambda x: len(x), reverse=True)


    @staticmethod
    def weight_function(length):
        return 1 /(1 + length)


    def save_to_cytoscape(self, output_prefix, sample_metadata=None):
        self._require_clusters("save clusters to Cytoscape")
        sif_filename = output_prefix + ".sif"
        properties_metadata_filename = output_prefix + ".prop.metadata.tsv"
        edges = []
        nodes = []

        additional_properties=[]
        for node in self.clusters[0]:
            additional_properties = list(node.additional_properties.keys())
            break
        for cluster in self.clusters:
            attributes = nx.get_edge_attributes(cluster,'length')
            for u,v in cluster.edges():
                node1_id = str(u.id)
                node2_id = str(v.id)
                try:
                    length = attributes[(u,v)]
                    edges.append(f"{node1_id}\t{length}\t{node2_id}")
                except KeyError:
                    try:
                        length = attributes[(v,u)]
                        edges.append(f"{node1_id}\t{length}\t{node2_id}")
                    except KeyError:
                        edges.append(f"{node1_id}\ttneighbour\t{node2_id}")
            if len(cluster) == 1:
                list(cluster.nodes())[0].id
                edges.append(str(list(cluster.nodes())[0]))
            for node in cluster:
                add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
                nodes.append((str(node.id), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.freq, node.count, *add_properties_values))

        with open(sif_filename, "w") as f:
            f.write("\n".join(edges))
        print("Saved edges to: {}".format(sif_filename))
        
        properties_names = ["code", "cdr3aa", "v", "j", "cdr3nt", "sample_id", "freq", "count"] + additional_properties
        properties_df = pd.DataFrame(nodes, columns=properties_names)
        if sample_metadata is not None:
            properties_df = properties_df.merge(sample_metadata)
        properties_df.to_csv(properties_metadata_filename, index=False, sep="\t")
        print("Saved node properties and metadata to: {}".format(properties_metadata_filename))


    def as_dataframe(self, filter_one_node_clusters=False):
        self._require_clusters("convert clusters to a dataframe")
        additional_properties=[]
        for node in self.clusters[0]:
            additional_properties = list(node.additional_properties.keys())
            break
        nodes = []
        if not filter_one_node_clusters:
            for cluster in self.clusters:
                for node in cluster:
                    add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
                    nodes.append((node.id, node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.freq, node.count, *add_properties_values))
        else:
            for cluster in self.filter_one_node_clusters(inplace=False):
                for node in cluster:
                    add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
                    nodes.append((node.id, node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.freq, node.count, *add_properties_values))
        properties_names = ["node_id", "cdr3aa", "v", "j", "cdr3nt", "sample_id", "freq", "count"] + additional_properties
        df = pd.DataFrame(nodes, columns=properties_names)
        first_columns = ["cluster_no", "node_id"]
        df = df[first_columns + [c for c in df.columns if c not in first_columns]]
        return df


    # def pool_clonotypes_to_df(self, folders, samples_list=None, top=0, functional=True, exclude_singletons=False, cdr3aa_len_range=[], metadata_filename="vdjtools_metadata.txt"):
    #     all_metadata = combine_metadata_from_folders(folders, metadata_filename=metadata_filename)
    #     #pool_metadata(folders, metadata_filename,samples_list)     
    #     clonotypes_dfs = []
    #     for index, row in all_metadata.iterrows():
    #         sample_id = row["sample.id"]
    #         if samples_list is not None:
    #             if sample_id not in samples_list:
    #                 continue
    #         clonoset_data=pd.read_csv(row["#file.name"],sep="\t")
    #         clonoset_data = clonoset_data.rename(columns={"bestVGene": "v",
    #                                         "bestJGene": "j",
    #                                         "CDR3.amino.acid.sequence": "cdr3aa",
    #                                         "CDR3.nucleotide.sequence": "cdr3nt",
    #                                         "allVHitsWithScore": "v",
    #                                         "allJHitsWithScore": "j",
    #                                         "aaSeqCDR3": "cdr3aa",
    #                                         "nSeqCDR3": "cdr3nt",
    #                                         "Sample":"sample_id",
    #                                         "cloneFraction":"freq",
    #                                         "Read.count": "count",
    #                                         "cloneCount": "count"})
        
    #         clonoset_data["v"] = clonoset_data["v"].apply(lambda x: x.split("*")[0])
    #         clonoset_data["j"] = clonoset_data["j"].apply(lambda x: x.split("*")[0])

    #         if exclude_singletons:
    #             clonoset_data=clonoset_data.loc[clonoset_data["count"]>1]
    #         if functional:
    #             clonoset_data=clonoset_data.loc[~clonoset_data["cdr3aa"].str.contains("\*|_")]
    #             clonoset_data=clonoset_data.sample(frac=1, random_state=1) #shuffle
    #             clonoset_data=clonoset_data.sort_values(by="count", ascending=False) #sort by counts "back" 
    #         if top > 0:
    #             clonoset_data=clonoset_data.iloc[:top]
    #         if cdr3aa_len_range:
    #             clonoset_data=clonoset_data.loc[(clonoset_data["cdr3aa"].str.len() <= cdr3aa_len_range[-1]) 
    #                                             & (clonoset_data["cdr3aa"].str.len() >= cdr3aa_len_range[0])]
    #         clonoset_data["freq"]=clonoset_data["count"]/clonoset_data["count"].sum()
    #         sample_id = row["sample.id"]
    #         # clonotypes_num = clonoset_data.shape[0]
    #         clonoset_data["sample_id"] = sample_id
    #         clonotypes_dfs.append(clonoset_data)
    # #         print("Added {} clonotypes from {}".format(clonotypes_num, sample_id))
    #     result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)
    #     self.pooled_clonosets = result_df
    #     clonotypes_number = len(result_df)
    #     samples_number = len(result_df["sample_id"].unique())
    #     print("Pooled {} clonotypes from {} samples".format(clonotypes_number, samples_number))
    #     return result_df


    def alice(self, 
        cl_filter=None, 
        overlap_type=None, 
        mismatches=None,
        generation_model='human_T_beta', 
        Q=9.41, 
        alpha=0.05, 
        olga_warnings=False,
        skip_single_nodes=False,
        method='bonferroni'):

        self._require_clusters("run ALICE")
        if overlap_type is None:
            overlap_type = self.overlap_type
        if mismatches is None:
            mismatches = self.mismatches
        if method not in ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by','fdr_tsbh', 'fdr_tsbky']:
            raise ValueError("P-value adjustment method is not one on the list. Possible values are: ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by','fdr_tsbh', 'fdr_tsbky']")
        if mismatches > 1:
            print(f'Using {mismatches} may increase runtime.')

        try:
            from .pgen_calculation import calculate_clonotypes_pgen
        except ImportError as exc:
            raise ImportError(
                "ALICE pgen calculation requires optional pgen dependencies. "
                "Install them with `pip install repseq[pgen]`."
            ) from exc
        
        if skip_single_nodes:
            clusters_all = self.as_dataframe()
        clusters = self.as_dataframe(filter_one_node_clusters=skip_single_nodes)
        clusters_all_pgen = calculate_clonotypes_pgen(clonosets_df=clusters, 
                                                                    cl_filter=cl_filter, 
                                                                    overlap_type=overlap_type, 
                                                                    mismatches=mismatches, 
                                                                    generation_model=generation_model, 
                                                                    olga_warnings=olga_warnings)

        aa, check_v, check_j = overlap_type_to_flags(overlap_type)
        if not check_v and not check_j:
            clusters_all_pgen['n'] = len(clusters_all_pgen)
        else:
            columns_to_check = ['v']
            if check_j:
                columns_to_check = ['v', 'j']
            group_counts = clusters_all_pgen[columns_to_check + ['cdr3aa']].groupby(by=columns_to_check).count().rename(columns={'cdr3aa': 'n'}).reset_index()
            clusters_all_pgen = clusters_all_pgen.merge(group_counts)
            
        alice_lambda = clusters_all_pgen['n'] * Q * clusters_all_pgen['pgen'].to_numpy()
        d = clusters_all_pgen['n_neighbours'].to_numpy()
        p = poisson.pmf(d, mu=alice_lambda)
        reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(p, alpha=alpha, method=method)
        clusters_all_pgen['p_value'] = p
        clusters_all_pgen['p_value_adj'] = pvals_corrected
        clusters_all_pgen['is_alice_hit'] = pvals_corrected < alpha
        if skip_single_nodes:
            self.alice_results = clusters_all.merge(clusters_all_pgen, how='left')
        else:
            self.alice_results = clusters_all_pgen.sort_values(by='p_value_adj').reset_index(drop=True)

        alice_parameters = {
            "overlap_type": overlap_type,
            "mismatches": mismatches,
            "generation_model": generation_model,
            "Q": Q,
            "alpha": alpha,
            "olga_warnings": olga_warnings,
            "skip_single_nodes": skip_single_nodes,
            "method": method,
        }
        if cl_filter is not None:
            alice_parameters["filter"] = self._filter_parameters(cl_filter)
        self.state["node_pgen_calculated"] = True
        self.state["alice_calculated"] = True
        self.state_parameters["node_pgen_calculated"] = {
            "overlap_type": overlap_type,
            "mismatches": mismatches,
            "generation_model": generation_model,
        }
        self.state_parameters["alice_calculated"] = alice_parameters
        return self.alice_results


    def add_alice_hits_to_clusters(self):
        self._require_alice("add ALICE hits to clusters")
        node_lookup = {}
        for cluster in self.clusters:
            for node in cluster:
                node_lookup[node.id] = node

        for _, row in self.alice_results.iterrows():
            node_id = row['node_id']
            node = node_lookup[node_id]
            node.additional_properties['pgen'] = row['pgen']
            node.additional_properties['alice_neighbour_count'] = row['alice_neighbour_count']
            node.additional_properties['p_value'] = row['p_value']
            node.additional_properties['p_value_adj'] = row['p_value_adj']
            node.additional_properties['is_alice_hit'] = row['is_alice_hit']


    # def export_clusters_to_gae(self):
        
    #     # additional_properties=[]
    #     # for node in clusters[0]:
    #     #     additional_properties = list(node.additional_properties.keys())
    #     #     break
        
    #     node_id = 0
    #     clone_id_dict = {}
    #     weights = []
    #     rows = []
    #     columns = []
    #     clone_list = []
    #     for cluster in self.clusters:
    #         for node in cluster:
    #             # add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
    #             # clone_list.append((str(node), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.size, *add_properties_values))
    #             clone_list.append((node_id, node.v, node.j, node.seq_aa, node.seq_nt, node.sample_id, node.size))
    #             clone_id_dict[node.id] = node_id
    #             weights.append(1)
    #             rows.append(node_id)
    #             columns.append(node_id)
    #             node_id += 1

    #         attributes = nx.get_edge_attributes(cluster,'length')
    #         for u,v in cluster.edges():
    #             length = attributes[(u,v)]
    #             weight = self.weight_function(length)
    #             u_id = clone_id_dict[u.id]
    #             v_id = clone_id_dict[v.id]
    #             weights.append(weight)
    #             weights.append(weight)
    #             rows.append(u_id)
    #             rows.append(v_id)
    #             columns.append(v_id)
    #             columns.append(u_id)

    #     col_names = ["node_id", "v", "j", "cdr3aa", "cdr3nt", "sample_id", "size"]
    #     # col_names += additional_properties
    #     clonoset = pd.DataFrame(clone_list, columns=col_names)
    #     clone_count = len(clonoset)

    #     # create sparse adjacency matrix with floats
    #     adjacency_matrix = csr_matrix((np.array(weights), (np.array(rows), np.array(columns))),
    #                         shape = (clone_count, clone_count), 
    #                         dtype = float).toarray()
    #     return (adjacency_matrix, clonoset, weights, rows, columns)

    #     def add_alice_hits_to_clusters(self, alice_hits_df, check_samples=True):
    #     nt_seq_colname = "CDR3.nucleotide.sequence"
    #     if nt_seq_colname not in alice_hits_df.columns:
    #         nt_seq_colname = "cdr3nt"
    #     if check_samples:
    #         samples = list(alice_hits_df["sample_id"].unique())
    #         alice_hits_dict = {}
    #         for sample in samples:
    #             hits = set()
    #             for index, row in alice_hits_df.loc[alice_hits_df["sample_id"] == sample].iterrows():
    #                 v = row["bestVGene"]
    #                 j = row["bestJGene"]
    #                 cdr3nt = row[nt_seq_colname]
    #                 hit = (cdr3nt, v, j)
    #                 hits.add(hit)
    #             alice_hits_dict[sample] = hits
    #         for cluster in self.clusters:
    #             for node in cluster:
    #                 if node.sample_id in samples and (node.seq_nt, node.v, node.j) in alice_hits_dict[node.sample_id]:
    #                     node.additional_properties["alice_hit"] = True
    #                 else:
    #                     node.additional_properties["alice_hit"] = False
    #     else:
    #         hits = set()
    #         for index, row in alice_hits_df.iterrows():
    #             v = row["bestVGene"]
    #             j = row["bestJGene"]
    #             cdr3nt = row[nt_seq_colname]
    #             hit = (cdr3nt, v, j)
    #             hits.add(hit)
    #             for cluster in self.clusters:
    #                 for node in cluster:
    #                     if (node.seq_nt, node.v, node.j) in hits:
    #                         node.additional_properties["alice_hit"] = True
    #                     else:
    #                         node.additional_properties["alice_hit"] = False


    # def filter_clusters_with_alice_hits(self):
    #     filtered_clusters = []
    #     for cluster in self.clusters:
    #         good_node = False
    #         for node in cluster:
    #             try:
    #                 a = node.additional_properties["alice_hit"]
    #             except KeyError:
    #                 print ("ALICE hits are not specified for these clusters. Specify this by ... [still not written]")
    #                 return None
    #             if node.additional_properties["alice_hit"]:
    #                 good_node = True
    #                 break
    #         if good_node:
    #             filtered_clusters.append(cluster)
    #             continue
    #     return filtered_clusters
    

    #    def pool_alice_hits_to_df(self, folders, samples_list=None, metadata_filename="vdjtools_metadata.txt"):
    #     all_metadata = combine_metadata_from_folders(folders, metadata_filename=metadata_filename)
    #     #pool_metadata(folders, metadata_filename,samples_list)
    #     clonotypes_dfs = []
    #     for index, row in all_metadata.iterrows():
    #         sample_id = row["sample.id"]
    #         if samples_list is not None:
    #             if sample_id not in samples_list:
    #                 continue
    #         clonoset_data=pd.read_csv(row["#file.name"],sep="\t")
    #         clonoset_data["freq"]=clonoset_data["Read.count"]/clonoset_data["Read.count"].sum()
    #         sample_id = row["sample.id"]
    #         # clonotypes_num = clonoset_data.shape[0]
    #         clonoset_data["sample_id"] = sample_id
    #         clonotypes_dfs.append(clonoset_data)
    # #         print("Added {} clonotypes from {}".format(clonotypes_num, sample_id))
    #     result_df = pd.concat(clonotypes_dfs).reset_index(drop=True)
    #     clonotypes_number = len(result_df)
    #     samples_number = len(result_df["sample_id"].unique())
    #     print("Pooled {} clonotypes from {} samples".format(clonotypes_number, samples_number))
    #     return result_df


class ClusterCommunities(list):

    def __init__(self):
        super().__init__() 
        self.communities = []
        self.resolution = None
        self.threshold = None
        self.seed = None
        

    # to enable list-like behaviour 
    def __getitem__(self, index):
        return self.communities[index]


    def __len__(self):
        return len(self.communities)


    def __iter__(self):
       return iter(self.communities)


    def __str__(self):
        total_communities = len(self.communities)
        possible_params = ['resolution', 'threshold', 'seed']
        params_used = {param: getattr(self, param, None)
                        for param in possible_params
                        if getattr(self, param, None) is not None}
        params_used = ''.join(f'{k}: {v}\n' for k, v in params_used.items())
        single_nodes = len(self.filter_one_node_communities(inplace=False))
        return f'Communities with {total_communities} nodes, of which {single_nodes} are single nodes.\nParameters:\n{params_used}'

    def __repr__(self):
        return self.__str__()
    

    @property
    def properties(self, weigh_by=None):
        properties_list = ["cluster_no", "cluster_id", "nodes", "edges", "diameter", "density", "eccentricity",
                       "concensus_cdr3aa", "concensus_cdr3nt", "concensus_v", "concensus_j"]
        results = []
        for cluster in self.communities:
            for node in cluster:
                break
            cluster_no = node.additional_properties["cluster_no"]
            cluster_id = f"cluster_{cluster_no}"
            average_eccentricity = np.mean(list(nx.eccentricity(cluster).values()))
            aa_consensus = cluster.calc_cluster_consensus(seq_type="prot", weigh_by=weigh_by)
            nt_consensus = cluster.calc_cluster_consensus(seq_type="dna", weigh_by=weigh_by)
            v_consensus = cluster.calc_cluster_consensus_segment(segment_type="v", weigh_by=weigh_by)
            j_consensus = cluster.calc_cluster_consensus_segment(segment_type="j", weigh_by=weigh_by)
            result = (cluster_no,
                    cluster_id,
                    len(cluster), 
                    nx.number_of_edges(cluster), 
                    nx.diameter(cluster),
                    nx.density(cluster), 
                    average_eccentricity,
                    aa_consensus,
                    nt_consensus,
                    v_consensus,
                    j_consensus)
            results.append(result)
        return pd.DataFrame(results, columns=properties_list)


    def filter_one_node_communities(self, inplace=False):
        if inplace:
            self.communities = [c for c in self.communities if len(c) > 1]
        else:
            return [c for c in self.communities if len(c) > 1]


def save_to_cytoscape(self, output_prefix, sample_metadata=None):
        sif_filename = output_prefix + ".sif"
        properties_metadata_filename = output_prefix + ".prop.metadata.tsv"
        edges = []
        nodes = []

        additional_properties=[]
        for node in self.communities[0]:
            additional_properties = list(node.additional_properties.keys())
            break
        for cluster in self.communities:
            attributes = nx.get_edge_attributes(cluster,'length')
            for u,v in cluster.edges():
                node1_id = str(u.id)
                node2_id = str(v.id)
                try:
                    length = attributes[(u,v)]
                    edges.append(f"{node1_id}\t{length}\t{node2_id}")
                except KeyError:
                    try:
                        length = attributes[(v,u)]
                        edges.append(f"{node1_id}\t{length}\t{node2_id}")
                    except KeyError:
                        edges.append(f"{node1_id}\ttneighbour\t{node2_id}")
            if len(cluster) == 1:
                list(cluster.nodes())[0].id
                edges.append(str(list(cluster.nodes())[0]))
            for node in cluster:
                add_properties_values = [node.additional_properties[add_property] for add_property in additional_properties]
                nodes.append((str(node.id), node.seq_aa, node.v, node.j, node.seq_nt, node.sample_id, node.freq, node.count, *add_properties_values))

        with open(sif_filename, "w") as f:
            f.write("\n".join(edges))
        print("Saved edges to: {}".format(sif_filename))
        
        properties_names = ["code", "cdr3aa", "v", "j", "cdr3nt", "sample_id", "freq", "count"] + additional_properties
        properties_df = pd.DataFrame(nodes, columns=properties_names)
        if sample_metadata is not None:
            properties_df = properties_df.merge(sample_metadata)
        properties_df.to_csv(properties_metadata_filename, index=False, sep="\t")
        print("Saved node properties and metadata to: {}".format(properties_metadata_filename))
        
    
