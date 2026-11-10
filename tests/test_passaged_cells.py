import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.passaged_cells import (
    across_lines,
    arm_slopes,
    arrest_contrast,
    deviations,
    matched_scores,
    score_samples,
)
from entropy_age.stochasticity import compare_groups


def test_matched_scores_is_compare_groups_over_one_sample():
    rng = np.random.default_rng(0)
    d = rng.gamma(2.0, 1.0, size=(3, 400))
    in_group = rng.random(400) < 0.2
    strata = rng.integers(0, 6, 400)
    strata[in_group & (strata == 5)] = 4  # stratum 5 holds no set CpG and must not count
    p, b = matched_scores(d, in_group, strata)
    features = [f"c{j}" for j in range(400)]
    membership = pl.DataFrame({"feature": [f for f, g in zip(features, in_group, strict=True) if g], "group": "s"})
    for i in range(3):
        metrics = pl.DataFrame({"feature": features, "d": d[i], "stratum": strata})
        expected = compare_groups(metrics, membership, metric="d", stratify_by="stratum", n_permutations=0)
        assert p[i] == pytest.approx(expected["matched_diff"][0])
    usable = np.isin(strata, np.unique(strata[in_group]))
    w = np.bincount(strata[in_group], minlength=6)[:6] / in_group.sum()
    out_means = np.array([d[:, ~in_group & (strata == k)].mean(axis=1) if w[k] else np.zeros(3) for k in range(6)])
    assert b == pytest.approx((w[:, None] * out_means).sum(axis=0))
    assert usable.sum() < 400
    with pytest.raises(ValueError, match="no stratum"):
        matched_scores(d, np.zeros(400, dtype=bool), strata)


def _culture(selection: bool, rng: np.random.Generator):
    """Three lines passaged 0..60 PD, and an MMC experiment; set CpGs are features 0..299."""
    n_cpg, n_set = 3000, 300
    start = rng.uniform(0.05, 0.95, n_cpg)
    samples, columns = [], []
    for line in ("L1", "L2", "L3"):
        direction = rng.choice([-1.0, 1.0], n_cpg)
        for step, pd in enumerate(np.linspace(0, 60, 9)):
            drift = 0.004 * pd * np.ones(n_cpg)
            # chromatin: set CpGs always drift half as fast; selection: set drift is purged as PD accumulate
            drift[:n_set] = drift[:n_set] * (np.exp(-pd / 15) if selection else 0.5)
            x = np.clip(start + direction * drift * start * (1 - start) * 4 + rng.normal(0, 0.004, n_cpg), 0, 1)
            columns.append(x)
            samples.append(
                {
                    "sample": f"{line}_{step}",
                    "subexperiment": "Baseline profiling",
                    "donor": line,
                    "doublings": pd,
                    "days": 5.0 * step,
                    "oxygen": "Ambient",
                    "agent": "NA",
                    "serum": "Full",
                    "vector": "NA",
                }
            )
    for arm, pd_step in (("NA", 0.0), ("MMC", 0.0), ("control", 2.5)):
        for rep, day in enumerate((8, 18, 25, 25, 25) if arm != "NA" else (8,)):
            if arm != "NA" and day == 8:
                continue
            pd = pd_step * day / 5
            drift = 0.004 * pd * np.ones(n_cpg) + 0.002
            drift[:n_set] *= np.exp(-pd / 15) if selection else 0.5
            columns.append(np.clip(start + drift * start * (1 - start) * 4 + rng.normal(0, 0.004, n_cpg), 0, 1))
            samples.append(
                {
                    "sample": f"M_{arm}_{rep}",
                    "subexperiment": "Growth arrest: MMC",
                    "donor": "L1",
                    "doublings": 5 + pd,
                    "days": float(day),
                    "oxygen": "Ambient",
                    "agent": arm,
                    "serum": "Full",
                    "vector": "NA",
                }
            )
    features = [f"cg{j:05d}" for j in range(n_cpg)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": rng.random(n_cpg) < 0.5})
    membership = pl.DataFrame({"feature": features[:n_set], "group": "housekeeping", "gene": "G"})
    return pl.DataFrame(samples), np.vstack(columns), features, probes, membership


@pytest.mark.parametrize("selection", [True, False])
def test_relative_gap_grows_with_divisions_only_under_selection(selection):
    samples, x, features, probes, membership = _culture(selection, np.random.default_rng(3))
    scores = score_samples(samples, x, features, probes, membership)
    assert scores.height == 3 * 8 + 2 * 4  # every sample but each experiment's earliest
    assert (scores.filter(pl.col("subexperiment") == "Baseline profiling")["doublings_gained"] > 0).all()
    lines = across_lines(arm_slopes(scores))
    f = lines.filter((pl.col("measure") == "F") & (pl.col("group") == "housekeeping")).row(0, named=True)
    b = lines.filter((pl.col("measure") == "B") & (pl.col("group") == "housekeeping")).row(0, named=True)
    assert f["lines"] == 3 and b["ci_low"] > 0  # deviation accumulates with divisions either way
    if selection:
        assert f["ci_high"] < 0 and f["negative"] == 3
    else:
        assert abs(f["mean_slope_per_10_pd"]) < 0.02
    arrest = arrest_contrast(scores).filter(pl.col("group") == "housekeeping").row(0, named=True)
    assert arrest["pd_gained_mmc"] == pytest.approx(0.0) and arrest["pd_gained_control"] > 5
    assert arrest["B_control"] > arrest["B_mmc"]
    if selection:
        assert arrest["F_control"] < arrest["F_mmc"]


def test_deviations_scale_by_the_binomial_spread():
    d = deviations(np.array([[0.6, 0.0]]), np.array([0.5, 0.0]))
    assert d[0, 0] == pytest.approx(0.2) and d[0, 1] == 0.0
