import numpy as np
import pytest

from entropy_age.clock import ClockConfig, evaluate, train_clock
from entropy_age.simulate.bulk import Bound, age_grid, random_ground_state, simulate_bulk
from entropy_age.simulate.methylation import MaintenanceRates, binomial_variance, simulate_methylation


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(7)


def test_age_grid_and_output_order(rng):
    ages = np.array([5, 0, 3, 5, 1])
    x = simulate_bulk(np.full(4, 0.5), ages, sigma=0.0, ground_sigma=0.0, rng=rng)
    assert x.shape == (5, 4)
    assert np.allclose(x, 0.5)
    assert age_grid(3, 2).tolist() == [0, 0, 1, 1, 2, 2]


def test_rows_follow_input_order(rng):
    # with a per-step drift only through clipping, older samples move further from the ground
    ground = np.full(200, 0.98)
    ages = np.array([40, 0, 20])
    x = simulate_bulk(ground, ages, sigma=0.05, ground_sigma=0.0, bound=Bound.CLIP, rng=rng)
    drift = np.abs(x - ground).mean(axis=1)
    assert drift[1] < drift[2] < drift[0]


def test_clip_keeps_unit_interval(rng):
    x = simulate_bulk(random_ground_state(300, rng), age_grid(30, 2), sigma=0.3, bound=Bound.CLIP, rng=rng)
    assert x.min() >= 0 and x.max() <= 1


def test_logit_every_step_is_the_same_model_as_logit():
    ground = np.random.default_rng(0).uniform(0.05, 0.95, 100)
    ages = age_grid(20, 2)
    a = simulate_bulk(ground, ages, sigma=0.2, bound=Bound.LOGIT, rng=np.random.default_rng(1))
    b = simulate_bulk(ground, ages, sigma=0.2, bound=Bound.LOGIT_EVERY_STEP, rng=np.random.default_rng(1))
    assert np.allclose(a, b)


def test_per_feature_sigma(rng):
    sigma = np.r_[np.zeros(50), np.full(50, 0.1)]
    x = simulate_bulk(np.full(100, 0.5), np.full(20, 30), sigma=sigma, ground_sigma=0.0, rng=rng)
    assert np.allclose(x[:, :50], 0.5)
    assert x[:, 50:].std() > 0.1


def test_bound_decides_whether_a_clock_exists(rng):
    ground = random_ground_state(400, rng)
    ages = age_grid(60, 2)
    results = {}
    for bound in (Bound.NONE, Bound.CLIP):
        clock = train_clock(
            simulate_bulk(ground, ages, sigma=0.05, bound=bound, rng=rng), ages, config=ClockConfig.fast()
        )
        results[bound] = evaluate(ages, clock.predict(simulate_bulk(ground, ages, sigma=0.05, bound=bound, rng=rng)))
    assert results[Bound.CLIP].pearson_r > 0.9
    assert not (results[Bound.NONE].pearson_r > 0.5)  # nan (constant fit) or weak


def test_maintenance_rate_algebra(rng):
    rates = MaintenanceRates.random(500, rng, em_limit=0.95, ed_limit=0.2)
    z = rates.equilibrium
    assert np.all((z >= 0) & (z <= 1))
    back = MaintenanceRates.from_rate_and_equilibrium(rates.relaxation_rate, z)
    assert np.allclose(back.em, rates.em) and np.allclose(back.ed, rates.ed)
    assert np.allclose(rates.expected(np.full(500, 0.3), 10_000), z, atol=1e-6)


def test_from_equilibrium_respects_limits_and_equilibrium(rng):
    eq = np.r_[0.0, 0.01, 0.2, 0.5, 0.9, 0.999, 1.0, rng.uniform(0, 1, 200)]
    rates = MaintenanceRates.from_equilibrium(eq, rng, em_limit=0.97, ed_limit=0.05)
    assert np.all(rates.em >= 0.97 - 1e-12) and np.all(rates.ed <= 0.05 + 1e-12)
    assert np.allclose(rates.equilibrium, np.clip(eq, 1e-4, 1 - 1e-4), atol=1e-9)


def test_binomial_matches_explicit_cells():
    rng = np.random.default_rng(3)
    ground = rng.uniform(0, 1, 30)
    rates = MaintenanceRates.random(30, rng, em_limit=0.9, ed_limit=0.1)
    steps, n_cells = 12, 100
    fast = simulate_methylation(ground, np.full(600, steps), rates, n_cells=n_cells, jitter=0, rng=rng)
    slow = simulate_methylation(ground, np.full(300, steps), rates, n_cells=n_cells, jitter=0, rng=rng, method="cells")
    start = np.floor(n_cells * ground) / n_cells
    expected_mean = rates.expected(start, steps)
    expected_var = binomial_variance(ground, steps, rates, n_cells)
    for sim, n in ((fast, 600), (slow, 300)):
        se = np.sqrt(expected_var / n)
        assert np.all(np.abs(sim.mean(axis=0) - expected_mean) < 5 * se + 1e-9)
        assert 0.8 < np.median(sim.var(axis=0) / expected_var) < 1.2


def test_universal_symmetric_rates_relax_to_half(rng):
    x = simulate_methylation(
        np.r_[np.zeros(20), np.ones(20)], np.full(5, 2000), MaintenanceRates.universal(40, 0.99), n_cells=500, rng=rng
    )
    assert abs(x.mean() - 0.5) < 0.02
