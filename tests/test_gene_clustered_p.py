import numpy as np



def test_gene_clustered_p_counts_genes_not_cpgs():
    """Twenty genes whose five CpGs share one shift: the CpG permutation calls it, the gene-level p is calibrated."""
    from entropy_age.stochasticity import gene_clustered_p

    rng = np.random.default_rng(3)
    false_calls = 0
    for _ in range(200):
        genes = np.repeat(np.arange(20), 5).astype(str)
        in_set = np.repeat(rng.normal(0, 1, 20), 5) + rng.normal(0, 0.2, 100)
        values = np.r_[in_set, rng.normal(0, 1, 2000)]
        labels = np.r_[np.ones(100, bool), np.zeros(2000, bool)]
        strata = rng.integers(0, 5, values.size)
        p, n = gene_clustered_p(values, labels, strata, genes)
        assert n == 20
        false_calls += p < 0.05
    assert false_calls / 200 < 0.09
