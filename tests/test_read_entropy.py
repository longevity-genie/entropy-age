"""Read-level epiallele measures: entropy, epipolymorphism, MHL, relative entropy, fixed-effect slopes, matched gap."""

import numpy as np
import pytest

from entropy_age.experiments.read_entropy import (
    bh_q,
    epiallele_metrics,
    epiallele_relative_entropy,
    fe_age_slopes,
    matched_gap,
    methylated_runs,
    young_reference,
)


def _counts(*pairs):
    """One window per argument: a dict {epiallele bit string: reads}."""
    out = np.zeros((len(pairs), 16), dtype=int)
    for i, pair in enumerate(pairs):
        for pattern, n in pair.items():
            out[i, int(pattern, 2)] = n
    return out


def test_entropy_and_epipolymorphism_on_toy_reads():
    m = epiallele_metrics(_counts({"0000": 40}, {"0000": 20, "1111": 20}, {f"{c:04b}": 5 for c in range(16)}))
    assert m["entropy"][0] == 0 and m["epipolymorphism"][0] == 0
    assert m["entropy"][1] == pytest.approx(1 / 4)  # one bit over four CpGs
    assert m["epipolymorphism"][1] == pytest.approx(0.5)
    assert m["entropy"][2] == pytest.approx(1.0) and m["epipolymorphism"][2] == pytest.approx(15 / 16)
    assert m["mean"][2] == pytest.approx(0.5)


def test_mhl_on_toy_reads():
    # all methylated: every stretch is fully methylated, MHL 1; all unmethylated: 0
    m = epiallele_metrics(_counts({"1111": 12}, {"0000": 12}, {"1010": 12}, {"1100": 12}))
    assert m["mhl"][0] == pytest.approx(1.0) and m["mhl"][1] == 0
    # 1010: length 1 -> 2/4 methylated, longer stretches never fully methylated: (1 x 0.5) / 10
    assert m["mhl"][2] == pytest.approx(0.05)
    # 1100: length 1 -> 2/4, length 2 -> 1/3: (0.5 + 2/3) / 10; same mean as 1010, more linkage
    assert m["mhl"][3] == pytest.approx((0.5 + 2 / 3) / 10)
    assert m["mean"][2] == m["mean"][3]


def test_methylated_runs_counts_stretches():
    runs = methylated_runs(4)
    assert runs[0b1111].tolist() == [4, 3, 2, 1]
    assert runs[0b0110].tolist() == [2, 1, 0, 0]


def test_low_coverage_is_missing():
    m = epiallele_metrics(_counts({"0000": 3}))
    assert np.isnan(m["entropy"][0]) and np.isnan(m["mhl"][0])


def test_relative_entropy_zero_at_reference_and_grows_with_departure():
    ref = _counts({"0000": 30})
    same = epiallele_relative_entropy(_counts({"0000": 20}), ref)
    moved = epiallele_relative_entropy(_counts({"0000": 15, "1000": 5}), ref)
    far = epiallele_relative_entropy(_counts({"0000": 10, "1000": 10}), ref)
    assert same[0] == pytest.approx(0, abs=1e-12)
    assert 0 < moved[0] < far[0]


def test_fixed_effect_slope_ignores_between_group_age_differences():
    rng = np.random.default_rng(0)
    group = np.repeat([0, 1, 2], 10)
    age = np.concatenate([rng.uniform(20, 30, 10), rng.uniform(40, 50, 10), rng.uniform(60, 70, 10)])
    female = rng.integers(0, 2, 30)
    level = np.array([0.0, 5.0, -3.0])[group]  # group offsets that a pooled fit would read as age
    y = (level + 0.1 * age + 0.4 * female)[:, None] + rng.normal(0, 0.01, (30, 2))
    y[3, 1] = np.nan  # a missing sample in one window
    slope, se, df = fe_age_slopes(y, age, female, group, min_samples=5)
    assert slope == pytest.approx([0.1, 0.1], abs=0.005)
    assert np.all(se < 0.01) and df.tolist() == [30 - 3 - 2, 29 - 3 - 2]


def test_matched_gap_compares_within_strata():
    values = np.array([1.0, 1.0, 0.0, 0.0, 10.0, 10.0, 9.0, 9.0, 50.0])
    member = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=bool)
    strata = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2])  # stratum 2 holds members only: dropped
    # stratum 0: 1 - 1/3; stratum 1: 10 - 28/3; weighted by one member each
    assert matched_gap(values, member, strata) == pytest.approx(((1 - 1 / 3) + (10 - 28 / 3)) / 2)


def test_bh_q():
    q = bh_q(np.array([0.01, 0.04, np.nan, 0.03]))
    assert q[0] == pytest.approx(0.03) and q[3] == pytest.approx(0.04) and q[1] == pytest.approx(0.04)
    assert np.isnan(q[2])


def test_young_reference_leaves_the_sample_out():
    counts = np.zeros((3, 1, 16), dtype=int)
    counts[0, 0, 0], counts[1, 0, 15], counts[2, 0, 3] = 10, 20, 30
    age, group = np.array([20.0, 21.0, 70.0]), np.zeros(3, dtype=int)
    ref = young_reference(counts, age, group, q=0.5)  # young: the two samples at or below age 21
    assert ref[0, 0, 15] == 20 and ref[0, 0, 0] == 0  # sample 0 against sample 1 only
    assert ref[1, 0, 0] == 10 and ref[1, 0, 15] == 0
    assert ref[2, 0, 0] == 10 and ref[2, 0, 15] == 20  # an old sample against both
    boot = young_reference(counts[[0, 0, 2]], age[[0, 0, 2]], group, ids=np.array([0, 0, 2]), q=0.5)
    assert boot[0].sum() == 0  # both copies of sample 0 are left out of its own reference
