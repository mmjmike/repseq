"""Sparse feature-by-sample count tables."""

import pandas as pd
from scipy.sparse import csr_matrix


class CountTable:
    """A CSR count matrix with aligned per-feature metadata and sample IDs.

    Rows are features, columns are samples. ``to_pandas()`` materializes the
    full wide table only when explicitly requested.
    """

    def __init__(self, matrix, features, sample_ids):
        self.matrix = csr_matrix(matrix)
        self.features = features.reset_index(drop=True).copy()
        self.index = features.index.copy()
        self.sample_ids = list(sample_ids)
        if self.matrix.shape != (len(self.features), len(self.sample_ids)):
            raise ValueError("matrix shape must match features and sample_ids")
        if self.features.columns.intersection(self.sample_ids).size:
            raise ValueError("feature columns and sample IDs must not overlap")
        self.attrs = {}

    def __len__(self):
        return self.matrix.shape[0]

    def __contains__(self, column):
        return column in self.columns

    @property
    def shape(self):
        return self.matrix.shape

    @property
    def columns(self):
        return self.features.columns.append(pd.Index(self.sample_ids))

    def __getitem__(self, column):
        if column in self.features:
            return self.features[column].set_axis(self.index)
        if column in self.sample_ids:
            return pd.Series(self.matrix[:, self.sample_ids.index(column)].toarray().ravel(), index=self.index)
        raise KeyError(column)

    def copy(self):
        result = CountTable(self.matrix.copy(), self.features, self.sample_ids)
        result.index = self.index.copy()
        result.attrs = self.attrs.copy()
        return result

    def to_pandas(self):
        counts = pd.DataFrame(self.matrix.toarray(), columns=self.sample_ids)
        result = pd.concat([self.features.reset_index(drop=True), counts], axis=1)
        result.index = self.index.copy()
        result.attrs = self.attrs.copy()
        return result

    def __call__(self, chain=None):
        if not hasattr(self, "_analyzer"):
            raise TypeError("Only Analyzer results support chain selection")
        return self._analyzer._get_result(self._result_name, chain=chain)
