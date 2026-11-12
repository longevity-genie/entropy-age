"""RM10b: shared-lineage covariance and the phylogenetic GLS."""

import numpy as np
import pytest

from entropy_age.data.taxonomy import binomial, lineage_covariance, shared_depth
from entropy_age.experiments.species_phylo import fit_lambda, gls, jackknife


def test_binomial_drops_subspecies():
    assert binomial("Ceratotherium simum simum") == "Ceratotherium simum"
    assert binomial("Mus musculus") == "Mus musculus"


def test_shared_depth_counts_leading_ancestors_only():
    assert shared_depth(["r", "a", "b", "x"], ["r", "a", "c", "x"]) == 2
    assert shared_depth(["r", "a"], ["q", "a"]) == 0


def test_covariance_orders_relatedness_and_has_unit_diagonal():
    mouse = ["r", "mam", "glires", "rodent", "mus"]
    rat = ["r", "mam", "glires", "rodent", "rattus"]
    whale = ["r", "mam", "cetacea", "balaena"]
    c = lineage_covariance([mouse, rat, whale])
    assert np.allclose(np.diag(c), 1.0)
    assert c[0, 1] > c[0, 2] > 0
    assert c[0, 2] == c[1, 2]
    assert np.allclose(c, c.T)


def _clades(n_clades: int, per: int):
    lineages = [["root", f"c{k}", f"s{k}_{i}"] for k in range(n_clades) for i in range(per)]
    return lineages, lineage_covariance(lineages)


def test_independent_noise_gives_lambda_near_zero_and_ols_coefficients():
    rng = np.random.default_rng(0)
    lineages, c = _clades(10, 6)
    n = len(lineages)
    x = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = 0.5 * x[:, 1] + rng.normal(0, 1, n)
    assert fit_lambda(y, x, c) < 0.4
    coef, *_ = gls(y, x, c, 0.0)
    assert coef == pytest.approx(np.linalg.lstsq(x, y, rcond=None)[0])


def test_clade_shared_noise_is_found_and_widens_the_lifespan_se():
    # a trait shared within clades and a predictor that also varies by clade: independence understates the SE
    rng = np.random.default_rng(1)
    lineages, c = _clades(8, 8)
    clade = np.repeat(np.arange(8), 8)
    n = clade.size
    pred = rng.normal(size=8)[clade]  # a clade-level trait, as lifespan nearly is
    x = np.column_stack([np.ones(n), pred])
    y = 2.0 * rng.normal(size=8)[clade] + 0.3 * rng.normal(size=n)  # no true effect of pred
    lam = fit_lambda(y, x, c)
    assert lam > 0.6
    _, cov0, _, _ = gls(y, x, c, 0.0)
    _, cov1, _, _ = gls(y, x, c, lam)
    # unit branch lengths let siblings share at most half their path below the deepest shared node, so the
    # taxonomy covariance captures only part of a strong clade effect: the SE widens, but not to the clade count
    assert np.sqrt(cov1[1, 1]) > 1.2 * np.sqrt(cov0[1, 1])


def test_a_species_in_several_tissues_does_not_make_the_fit_singular():
    rng = np.random.default_rng(2)
    lineages, _ = _clades(6, 4)
    rows = lineages + lineages  # every species twice, as in two tissues
    c = lineage_covariance(rows)
    n = len(rows)
    x = np.column_stack([np.ones(n), rng.normal(size=n), np.repeat([0.0, 1.0], n // 2)])
    y = rng.normal(size=n)
    lam = fit_lambda(y, x, c)
    assert 0.0 <= lam < 1.0
    assert np.isfinite(gls(y, x, c, lam)[0]).all()


def test_jackknife_names_the_species_that_carries_a_slope():
    rng = np.random.default_rng(4)
    species = [f"s{i}" for i in range(20)]
    x = np.column_stack([np.ones(20), np.linspace(0, 1, 20)])
    y = rng.normal(0, 0.01, 20)
    y[-1] += 1.0  # one species alone makes the slope
    j = jackknife(y, x, species)
    assert j["influential_species"] == "s19"
    assert abs(j["coef_without"]) < 0.05
    assert j["coef_min"] <= j["coef_without"] <= j["coef_max"]
