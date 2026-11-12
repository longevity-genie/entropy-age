"""RM19: per-cell-type pseudobulk from a CSR count matrix, and the donor ids of the demographics table."""

import numpy as np
from scipy import sparse

from entropy_age.experiments.atac_cortex import donor_id, group_sums, parse_peaks


def test_parse_peaks():
    p = parse_peaks(["chr1:100-600", "chrX:5-10"])
    assert (
        p["chr"].to_list() == ["chr1", "chrX"] and p["start"].to_list() == [100, 5] and p["end"].to_list() == [600, 10]
    )


def test_hbcc_donors_are_bare_numbers_in_the_table():
    assert donor_id(887.0, "HBCC") == "HBCC-887"
    assert donor_id("KEN-845", "UKY") == "KEN-845"


def test_group_sums_match_dense_sums_over_kept_columns_and_totals_over_all():
    rng = np.random.default_rng(0)
    dense = rng.poisson(0.5, (30, 12)).astype(float)
    m = sparse.csr_matrix(dense)
    groups = rng.integers(0, 4, 30)
    keep = np.array([1, 4, 7, 11])
    keep_col = np.full(12, -1)
    keep_col[keep] = np.arange(keep.size)
    sums, totals = group_sums(m.indptr, m.indices, m.data, groups, keep_col, 4, keep.size)
    for g in range(4):
        assert np.allclose(sums[g], dense[groups == g][:, keep].sum(axis=0))
        assert np.isclose(totals[g], dense[groups == g].sum())
