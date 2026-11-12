"""RM19: per-subclass sums from a CSC 10x matrix, and gene counts over distinct peaks."""

import numpy as np
import polars as pl
from scipy import sparse

from entropy_age.experiments.atac_hippocampus import barcode_of, gene_counts, subclass_peak_sums


def test_subclass_sums_match_dense_column_sums():
    rng = np.random.default_rng(0)
    dense = rng.poisson(0.4, (8, 20)).astype(float)  # features x barcodes
    m = sparse.csc_matrix(dense)
    groups = {"a": np.array([0, 3, 7]), "b": np.array([5]), "empty": np.array([], dtype=int)}
    out = subclass_peak_sums(m.indptr, m.indices, m.data, groups, 8)
    assert np.allclose(out["a"], dense[:, [0, 3, 7]].sum(axis=1))
    assert np.allclose(out["b"], dense[:, 5])
    assert np.allclose(out["empty"], 0)


def test_gene_counts_add_each_peak_once_per_gene():
    peaks = np.array([10.0, 20.0, 30.0])
    pairs = pl.DataFrame({"peak": [0, 0, 1, 2], "gene": ["A", "A", "A", "B"]})  # peak 0 holds two CpGs of A
    assert gene_counts(peaks, pairs, ["A", "B", "C"]).tolist() == [30.0, 30.0, 0.0]


def test_barcode_of_drops_donor_and_library_tags():
    assert barcode_of("hc77_deep_AAACAGCC-1") == "AAACAGCC-1"
    assert barcode_of("hc1134_AAACAGCC-1") == "AAACAGCC-1"
