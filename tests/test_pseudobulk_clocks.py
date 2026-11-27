import numpy as np
import polars as pl
import pytest

from entropy_age.experiments import pseudobulk_clocks as pc


def test_canonical_key_moves_bottom_strand_to_top_c():
    assert pc.canonical_key("chr1:100", "+") == "1:100"
    assert pc.canonical_key("chrX:101", "-") == "X:100"
    with pytest.raises(ValueError):
        pc.canonical_key("chr1:100", "")


def test_to_cpgs_pools_both_strands_onto_the_top_c():
    keys = pl.DataFrame({"key": [pc.key_of("1:100")]})
    cov = pl.DataFrame(
        {"key": [pc.key_of("1:100"), pc.key_of("1:101"), pc.key_of("1:102")], "meth": [1, 2, 5], "total": [1, 3, 5]}
    )
    out = pc.to_cpgs(cov, keys)
    assert out.height == 1
    assert out["meth"][0] == 3 and out["total"][0] == 4


def test_binarise_drops_ties_and_empty():
    out = pc.binarise(np.array([2.0, 0.0, 1.0, 0.0]), np.array([2.0, 2.0, 2.0, 0.0]))
    assert out[0] == 1 and out[1] == 0
    assert np.isnan(out[2]) and np.isnan(out[3])


def test_scage_recovers_a_planted_age():
    rng = np.random.default_rng(0)
    n = 3000
    slope = rng.uniform(-0.02, 0.02, n)
    intercept = rng.uniform(0.3, 0.7, n)
    r = np.sign(slope) * np.abs(slope) / 0.02
    truth = 18.0
    p = np.clip(intercept + slope * truth, 0.001, 0.999)
    state = (rng.random(n) < p).astype(float)
    age, used = pc.scage_cell(state, intercept, slope, r, top=2000)
    assert used == 2000
    assert abs(age - truth) < 4


def test_scage_reference_fits_a_line():
    age = np.array([2.0, 2.0, 10.0, 10.0, 20.0, 20.0])
    beta = np.column_stack([0.1 + 0.01 * age, 0.9 - 0.02 * age])
    beta[0, 1] = np.nan
    a, b, r = pc.scage_reference(beta, age)
    assert np.allclose(a, [0.1, 0.9]) and np.allclose(b, [0.01, -0.02])
    assert np.allclose(r, [1.0, -1.0])


def test_within_group_link_uses_every_relabelling():
    x = np.array([1.0, 2.0, 3.0, 11.0, 12.0, 13.0])
    y = np.array([0.1, 0.2, 0.3, 5.1, 5.2, 5.3])
    group = np.array([0, 0, 0, 1, 1, 1])
    out = pc.within_group_link(x, y, group)
    assert out["permutations"] == 36
    assert out["r"] == pytest.approx(1.0)
    # only the identity and nothing else reach |r| = 1 within groups: two of 36 (r = +1 and r = -1 relabellings)
    assert out["perm_p"] == pytest.approx(2 / 36)


def test_fill_missing_takes_the_first_finite_source():
    beta = np.array([[np.nan, 0.2, np.nan]])
    out = pc.fill_missing(beta, [np.array([0.7, np.nan, np.nan]), np.array([np.nan, np.nan, 0.9])])
    assert np.allclose(out, [[0.7, 0.2, 0.9]])
    assert np.allclose(pc.fill_missing(np.array([[np.nan]]), []), [[0.5]])


def test_leave_one_out_mean():
    beta = np.array([[1.0, np.nan], [3.0, 2.0], [5.0, np.nan]])
    out = pc.leave_one_out_mean(beta)
    assert np.allclose(out[:, 0], [4.0, 3.0, 2.0])
    assert out[0, 1] == pytest.approx(2.0) and np.isnan(out[1, 1])


def test_weight_share_counts_absolute_weight():
    clock = pc.Clock("x", np.array(["a", "b", "c"]), np.array([1.0, -3.0, 0.0]), 0.0, np.full(3, np.nan))
    share = pc.weight_share(clock, np.array([[True, False, True], [False, True, False]]))
    assert np.allclose(share, [0.25, 0.75])


def test_quantile_to_target_averages_ties():
    target = np.array([0.0, 1.0, 2.0, 3.0])
    out = pc.quantile_to_target(np.array([[0.5, 0.1, 0.5, 0.9]]), target)
    assert np.allclose(out, [[1.5, 0.0, 1.5, 3.0]])


def test_tage_prepare_centres_the_reference():
    rng = np.random.default_rng(3)
    symbols = np.array([f"G{i}" for i in range(200)])
    lookup = {f"G{i}": str(i) for i in range(150)}
    counts = rng.poisson(rng.lognormal(4, 1, 200), size=(8, 200)).astype(float)
    reference = np.array([True] * 4 + [False] * 4)
    x, ids = pc.tage_prepare(counts, symbols, lookup, reference)
    assert ids.size <= 150
    assert np.allclose(np.median(x[reference], axis=0), 0.0)


def test_residual_z_removes_group_means():
    age = np.array([3, 3, 3, 24, 24, 24])
    sex = np.array(["f", "m", "f", "m", "f", "m"])
    y = np.array([1.0, 2.0, 3.0, 11.0, 12.0, 13.0])
    z = pc.residual_z(y, age, sex)
    assert np.isfinite(z).all()
    assert abs(z.sum()) < 1e-9
    assert np.isnan(pc.residual_z(np.array([1.0, 2.0]), np.array([3, 24]), np.array(["f", "f"]))).all()


def test_scage_refuses_a_cell_with_too_few_cpgs():
    state = np.array([1.0, 0.0, np.nan])
    age, used = pc.scage_cell(state, np.full(3, 0.5), np.full(3, 0.01), np.full(3, 0.5), top=3)
    assert np.isnan(age) and used == 2


def test_matrix_fills_both_inputs_that_name_one_cpg():
    long = pl.DataFrame({"sample": ["a", "b"], "key": [7, 7], "beta": [0.2, 0.4]})
    out = pc._matrix(long, ["a", "b"], "sample", np.array([7, 9, 7]), "beta")
    assert np.allclose(out[:, 0], [0.2, 0.4]) and np.allclose(out[:, 2], [0.2, 0.4])
    assert np.isnan(out[:, 1]).all()


def test_tage_prepare_refuses_a_flat_sample():
    counts = np.full((3, 5), 50.0)
    counts[1] = [10, 20, 30, 40, 50]
    with pytest.raises(ValueError):
        pc.tage_prepare(counts, np.array(list("abcde")), {c.upper(): c for c in "abcde"}, np.array([True, True, False]))


def test_residual_z_of_a_constant_is_nan():
    age = np.array([3, 3, 3, 24, 24, 24])
    sex = np.array(["f", "m", "f", "m", "f", "m"])
    assert np.isnan(pc.residual_z(np.full(6, 2.0), age, sex)).all()
