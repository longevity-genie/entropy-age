import numpy as np
import polars as pl
import pytest

from entropy_age.clock import HorvathAge, LinearRescale, LogLogRelativeAge, evaluate
from entropy_age.simulate.methylation import MaintenanceRates
from entropy_age.stochasticity import (
    binary_entropy,
    compare_groups,
    drift_metrics,
    fit_relaxation,
    paired_comparison,
    regress_out,
)
from entropy_age.transcriptome import binarize_by_median, correct_biological_age


@pytest.mark.parametrize(
    "transform", [HorvathAge(), LogLogRelativeAge(67), LinearRescale(1, 73)], ids=["horvath", "loglog", "rescale"]
)
def test_transforms_round_trip(transform):
    ages = np.array([1.0, 5.0, 19.0, 21.0, 40.0, 66.0])
    back = transform.inverse(transform.forward(ages))
    expected = ages / 67 if isinstance(transform, LogLogRelativeAge) else ages
    assert np.allclose(back, expected)


def test_evaluate_handles_constant_predictions():
    m = evaluate([1, 2, 3, 4], [2, 2, 2, 2])
    assert m.n == 4 and np.isnan(m.pearson_r)
    assert evaluate([1, 2, 3, 4], [1, 2, 3, 5]).pearson_r > 0.95


def test_binary_entropy_bounds():
    assert binary_entropy(np.array([0.0, 1.0])).max() < 1e-6
    assert binary_entropy(np.array([0.5]))[0] == pytest.approx(1.0)


def _cohort(rng, rates, n_people=300, n_cells=2000):
    ages = rng.uniform(0, 80, n_people)
    ground = np.full(len(rates), 0.05)
    mean = rates.expected(ground, ages)
    x = rng.binomial(n_cells, mean) / n_cells
    return x, ages


def test_fit_relaxation_recovers_rates():
    rng = np.random.default_rng(5)
    true_rate = np.r_[np.full(50, 0.01), np.full(50, 0.05)]
    rates = MaintenanceRates.from_rate_and_equilibrium(true_rate, np.full(100, 0.6))
    x, ages = _cohort(rng, rates)
    fit = fit_relaxation(x, ages, t0=0.0)
    slow, fast = fit["rate"].to_numpy()[:50], fit["rate"].to_numpy()[50:]
    assert np.median(slow) == pytest.approx(0.01, rel=0.3)
    assert np.median(fast) == pytest.approx(0.05, rel=0.3)
    assert np.median(fit["equilibrium"].to_numpy()[50:]) == pytest.approx(0.6, abs=0.05)


def test_drift_metrics_and_group_contrast():
    rng = np.random.default_rng(6)
    true_rate = np.r_[np.full(60, 0.005), np.full(240, 0.03)]
    rates = MaintenanceRates.from_rate_and_equilibrium(true_rate, np.full(300, 0.5))
    x, ages = _cohort(rng, rates)
    names = [f"cg{i}" for i in range(300)]
    metrics = drift_metrics(x, ages, names)
    assert (metrics["entropy_slope"] > 0).mean() > 0.95
    membership = pl.DataFrame({"feature": names[:60], "group": "protected"})
    # every site starts at 0.05, so an independent reference baseline puts all of them in one stratum
    metrics = metrics.with_columns(reference=pl.lit(0.05) + pl.Series(rng.normal(0, 1e-3, 300)))
    result = compare_groups(metrics, membership, metric="entropy_slope", stratify_by="reference", n_permutations=300)
    row = result.row(0, named=True)
    assert row["matched_diff"] < 0 and row["p_permutation"] < 0.01 and row["n_matched"] == 60

    null = pl.DataFrame({"feature": list(rng.choice(names, 60, replace=False)), "group": "random"})
    assert (
        compare_groups(metrics, null, metric="entropy_slope", stratify_by="reference", n_permutations=300)[
            "p_permutation"
        ][0]
        > 0.01
    )


def test_in_cohort_baseline_separates_fast_and_slow_sites():
    """Documents the pitfall: early drift moves fast sites into other baseline bins."""
    rng = np.random.default_rng(6)
    true_rate = np.r_[np.full(60, 0.005), np.full(240, 0.03)]
    x, ages = _cohort(rng, MaintenanceRates.from_rate_and_equilibrium(true_rate, np.full(300, 0.5)))
    names = [f"cg{i}" for i in range(300)]
    metrics = drift_metrics(x, ages, names)
    result = compare_groups(
        metrics,
        pl.DataFrame({"feature": names[:60], "group": "slow"}),
        metric="entropy_slope",
        stratify_by="baseline",
        n_permutations=100,
    )
    assert result["n_matched"][0] < 60


def test_paired_comparison():
    a = pl.DataFrame({"feature": [f"f{i}" for i in range(50)], "m": np.linspace(1, 2, 50)})
    b = a.with_columns(pl.col("m") - 0.5)
    result = paired_comparison(a, b, metric="m", names=("neuron", "blood"))
    assert result["fraction_neuron_greater"] == 1.0 and result["p_wilcoxon"] < 1e-6


def test_regress_out_removes_covariate_signal():
    rng = np.random.default_rng(8)
    cov = rng.uniform(0, 1, (200, 2))
    x = 0.3 + cov @ np.array([[0.2, -0.1], [0.05, 0.3]]) + 0.01 * rng.standard_normal((200, 2))
    adjusted = regress_out(x, cov)
    assert np.allclose(adjusted.mean(axis=0), x.mean(axis=0))
    assert np.abs(np.corrcoef(np.c_[adjusted, cov].T)[:2, 2:]).max() < 0.05


def test_transcriptome_helpers():
    x = np.array([[0.0, 0.2], [0.5, 0.4], [0.9, 0.6], [0.7, 0.0]])
    assert binarize_by_median(x).tolist() == [[0, 0], [0, 0], [1, 1], [0, 0]]
    assert np.isfinite(correct_biological_age(np.array([100.0, 300.0]))).all()


def test_batch_adjusted_slopes_remove_a_batch_age_confound():
    """Old donors run on another batch: unadjusted slopes see the batch shift, partial slopes do not."""
    rng = np.random.default_rng(9)
    n = 400
    batch = np.r_[np.zeros(200), np.ones(200)]
    ages = np.where(batch == 0, rng.uniform(16, 55, n), rng.uniform(45, 88, n))
    true_slope = np.r_[np.zeros(50), np.full(50, 0.002)]
    shift = 0.05  # batch offset on every feature
    x = 0.3 + np.outer(ages, true_slope) + shift * batch[:, None] + rng.normal(0, 0.01, (n, 100))
    entropy_shift = binary_entropy(np.full(1, 0.3 + shift))[0] - binary_entropy(np.full(1, 0.3))[0]
    naive = drift_metrics(x, ages)
    adjusted = drift_metrics(x, ages, covariates=batch)
    assert naive["slope"].to_numpy()[:50].mean() > 5e-4  # spurious drift from the batch
    assert abs(adjusted["slope"].to_numpy()[:50].mean()) < 1e-4
    assert np.allclose(adjusted["slope"].to_numpy()[50:].mean(), 0.002, rtol=0.05)
    assert entropy_shift > 0 and abs(adjusted["entropy_slope"].to_numpy()[:50].mean()) < 1e-4
    assert abs(adjusted["spearman"].to_numpy()[:50].mean()) < 0.05 < naive["spearman"].to_numpy()[:50].mean()

    # partial rank correlation for a binary batch = correlation of within-batch-centred ranks
    from scipy.stats import rankdata

    def centred(v):
        return v - np.where(batch == 1, v[batch == 1].mean(), v[batch == 0].mean())

    rt = centred(rankdata(ages))
    expected = [np.corrcoef(rt, centred(rankdata(x[:, j])))[0, 1] for j in (50, 75, 99)]
    assert np.allclose(adjusted["spearman"].to_numpy()[[50, 75, 99]], expected)


def test_composition_components_recover_a_shared_cell_mix_factor():
    from entropy_age.stochasticity import composition_components

    rng = np.random.default_rng(3)
    n, p = 240, 3000
    batch = rng.integers(0, 2, n).astype(float)
    ages = rng.uniform(20, 90, n)
    mix = rng.normal(size=n) + 0.02 * ages  # a cell-mix factor that partly tracks age
    loaded = rng.random(p) < 0.3  # cell-type-specific CpGs follow the mix
    x = 0.5 + 0.03 * rng.normal(size=(n, p)) + 0.05 * batch[:, None]
    x[:, loaded] += 0.04 * mix[:, None] * rng.uniform(0.5, 1.5, loaded.sum())

    scores, explained = composition_components(x, 3, covariates=batch[:, None])
    assert scores.shape == (n, 3)
    assert np.all(np.diff(explained) <= 0) and explained.sum() < 1
    design = np.column_stack([np.ones(n), batch])
    assert np.abs(design.T @ scores).max() < 1e-6 * np.abs(scores).sum()  # orthogonal to batch and intercept
    mix_free = mix - design @ np.linalg.lstsq(design, mix, rcond=None)[0]
    assert abs(np.corrcoef(scores[:, 0], mix_free)[0, 1]) > 0.95  # the first component is the cell mix

    # with age among the covariates, the components carry no age trend
    age_free, _ = composition_components(x, 3, covariates=np.column_stack([batch, ages]))
    assert np.abs(np.corrcoef(np.column_stack([age_free, ages]).T)[-1, :-1]).max() < 1e-6
    assert abs(np.corrcoef(scores[:, 0], ages)[0, 1]) > 0.2  # but by default they keep it
