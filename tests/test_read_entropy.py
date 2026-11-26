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


# ----------------------------------------------------------------------------- review 7
def test_matched_sides_report_both_sides_of_the_gap():
    from entropy_age.experiments.read_entropy import matched_sides

    values = np.array([-1.0, -3.0, 2.0, 2.0, 5.0, 5.0, 6.0, 8.0])
    member = np.array([1, 1, 0, 0, 1, 0, 0, 0], dtype=bool)
    strata = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    inside, matched = matched_sides(values, member, strata)
    # the matched side weights each stratum by its set windows: two in stratum 0, one in stratum 1
    assert inside == pytest.approx((-1 - 3 + 5) / 3)
    assert matched == pytest.approx((2 * 2 + 19 / 3) / 3)
    assert matched_gap(values, member, strata) == pytest.approx(inside - matched)


def test_sign_tail_drops_nan_draws_and_is_never_zero():
    from entropy_age.experiments.read_entropy import sign_tail

    p, n = sign_tail(np.array([0.5] * 99 + [np.nan] * 50))
    assert n == 99 and p == pytest.approx(2 / 100)
    # NaN draws are not counted on the non-positive side
    p_nan, _ = sign_tail(np.array([1.0] * 10 + [np.nan] * 10))
    assert p_nan == pytest.approx(2 / 11)
    assert sign_tail(np.array([np.nan]))[1] == 0
    assert sign_tail(np.array([1.0, -1.0]))[0] == 1.0


def test_cov_age_slopes_matches_least_squares_with_missing_values():
    from entropy_age.experiments.read_entropy import cov_age_slopes

    rng = np.random.default_rng(4)
    n = 40
    age = rng.uniform(20, 80, n)
    cov = np.column_stack([rng.integers(0, 2, n), rng.normal(0, 1, n) + 0.02 * age])
    y = (0.05 * age + 0.3 * cov[:, 1])[:, None] + rng.normal(0, 0.2, (n, 3))
    y[[1, 5], 0] = np.nan
    slope, se, df = cov_age_slopes(y, age, cov, min_samples=10)
    for j in range(3):
        ok = np.isfinite(y[:, j])
        x = np.column_stack([np.ones(ok.sum()), age[ok], cov[ok]])
        coef, rss, *_ = np.linalg.lstsq(x, y[ok, j], rcond=None)
        assert slope[j] == pytest.approx(coef[1], rel=1e-8)
        var = rss[0] / (ok.sum() - 4) * np.linalg.inv(x.T @ x)[1, 1]
        assert se[j] == pytest.approx(np.sqrt(var), rel=1e-6) and df[j] == ok.sum() - 4
    # a covariate equal to age leaves age unidentifiable
    assert np.isnan(cov_age_slopes(y, age, age[:, None], min_samples=10)[0]).all()


def test_rarefy_draws_exact_depth_without_replacement():
    from entropy_age.experiments.read_entropy import rarefy

    counts = np.zeros((2, 3, 16), dtype=int)
    counts[0, 0, [0, 15]] = [30, 10]
    counts[0, 1, 3] = 5  # below depth: dropped
    counts[1, 2, [1, 2, 4]] = [4, 4, 4]
    out = rarefy(counts, 10, np.random.default_rng(0))
    assert out.sum(axis=-1).tolist() == [[10, 0, 0], [0, 0, 10]]
    assert (out <= counts).all()
    many = np.stack([rarefy(counts[:1, :1], 10, np.random.default_rng(s))[0, 0] for s in range(400)])
    assert many[:, 0].mean() == pytest.approx(7.5, abs=0.15)  # 3/4 of ten fragments


def test_balanced_windows_compare_younger_and_older_halves_within_groups():
    from entropy_age.experiments.read_entropy import balanced_windows

    age = np.array([20.0, 30.0, 60.0, 70.0, 25.0, 65.0])
    group = np.array([0, 0, 0, 0, 1, 1])
    reads = np.array([[30, 30], [30, 30], [30, 5], [30, 5], [30, 30], [30, 30]])
    # window 1: every young sample reaches 20 fragments, one old sample in three
    assert balanced_windows(reads, age, group, 20).tolist() == [True, False]


def test_fractions_recover_a_known_mixture():
    from entropy_age.experiments.read_entropy import estimate_fractions, marker_windows

    rng = np.random.default_rng(5)
    types = rng.uniform(0.05, 0.95, (4, 600))
    truth = np.array([[0.6, 0.2, 0.1, 0.1], [0.1, 0.1, 0.3, 0.5]])
    samples = truth @ types + rng.normal(0, 0.005, (2, 600))
    markers = marker_windows(types, per_side=20)
    assert 40 <= markers.size <= 160
    assert np.allclose(estimate_fractions(samples, types, markers), truth, atol=0.03)


def test_mixture_counts_follow_the_mixed_shares_and_depths():
    from entropy_age.experiments.read_entropy import mixture_counts

    shares = np.zeros((2, 2, 16))
    shares[0, :, 0] = 1.0  # type 0: all unmethylated
    shares[1, :, 15] = 1.0  # type 1: all methylated
    shares[:, 1] = np.nan  # window 1 missing in the references
    reads = np.array([[1000, 50], [1000, 50]])
    out = mixture_counts(shares, np.array([[0.9, 0.1], [0.5, 0.5]]), reads, np.random.default_rng(0))
    assert out[:, 0].sum(axis=-1).tolist() == [1000, 1000] and out[:, 1].sum() == 0
    assert out[0, 0, 15] / 1000 == pytest.approx(0.1, abs=0.03)
    assert out[1, 0, 15] / 1000 == pytest.approx(0.5, abs=0.05)


def test_putative_donors_join_samples_of_one_body():
    import polars as pl

    from entropy_age.experiments.read_entropy import putative_donors

    s = pl.DataFrame(
        {
            "age": [19.0, 19.0, 19.0, 25.0],
            "female": [False, False, False, False],
            "lab": ["UA", "UA", "KI", "KI"],
            "tissue": ["Pancreas", "Pancreas", "Brain", "Brain"],
        }
    )
    d = putative_donors(s)
    assert d[0] == d[1] and len(set(d)) == 3


def test_cache_key_depends_on_window_order_and_file(tmp_path):
    import os

    from entropy_age.experiments.read_entropy import cache_key

    pat = tmp_path / "GSM1_x.hg38.pat.gz"
    pat.write_bytes(b"abc")
    starts = np.array([5, 9, 12])
    key = cache_key(starts, pat)
    assert key == cache_key(starts.copy(), pat)
    assert key != cache_key(starts[::-1], pat)  # same size, other order: another panel
    assert key != cache_key(np.array([5, 9, 13]), pat)
    pat.write_bytes(b"abcd")
    os.utime(pat, (1, 1))
    assert key != cache_key(starts, pat)


def test_each_measure_keeps_its_own_windows():
    from entropy_age.experiments.read_entropy import MEASURES, _keeps

    fits = {m: (np.array([0.1, 0.2, 0.3]), None, None) for m in MEASURES}
    fits["relative_entropy"] = (np.array([np.nan, 0.2, np.nan]), None, None)
    keep = _keeps(fits, np.array([0.5, 0.5, np.nan]))
    assert keep["entropy"].tolist() == [True, True, False]  # KL missing does not drop a Shannon window
    assert keep["relative_entropy"].tolist() == [False, True, False]
    assert keep["log_depth"].tolist() == keep["entropy"].tolist()
