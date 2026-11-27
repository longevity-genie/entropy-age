import numpy as np
import polars as pl
import pytest
from scipy import sparse

from entropy_age.experiments.tms_cell_noise import (
    between_spread,
    distinct_mice_ok,
    eligible,
    gene_boot_diffs,
    matched_count,
    matched_diff,
    mouse_draw,
    mouse_weights,
    noise_residual,
    one_to_one_orthologs,
    plant_noise,
    pool_gaps,
    strata,
    thin_rows,
    trend_residual,
    within_moments,
)


def test_thinning_reaches_the_target_depth_and_refuses_shallow_rows():
    rng = np.random.default_rng(0)
    x = sparse.csr_matrix(rng.poisson(50, size=(200, 400)).astype(np.float32))
    out = thin_rows(x, 5000, rng)
    totals = np.asarray(out.sum(axis=1)).ravel()
    assert abs(totals.mean() - 5000) < 30
    assert np.all(out.toarray() <= x.toarray())
    with pytest.raises(ValueError):
        thin_rows(x, 10**7, rng)


def test_within_mouse_variance_ignores_differences_between_mouse_means():
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 1, size=(400, 3))
    mouse = np.repeat(["a", "b", "c", "d"], 100)
    shift = np.repeat([0.0, 10.0, -10.0, 20.0], 100)[:, None]
    mean, var, dropout = within_moments(noise + shift + 50, mouse)
    assert np.allclose(var, 1.0, atol=0.2)
    assert np.allclose(mean, 55, atol=0.3)
    assert np.all(dropout == 0)


def test_trend_residual_removes_a_mean_trend_and_keeps_a_planted_excess():
    rng = np.random.default_rng(2)
    log_mean = rng.uniform(-2, 6, 3000)
    log_cv2 = -log_mean + rng.normal(0, 0.1, 3000)
    log_cv2[:100] += 2.0
    res = trend_residual(log_mean, log_cv2)
    assert abs(np.median(res[100:])) < 0.05
    assert np.median(res[:100]) > 1.8
    assert np.isnan(trend_residual(np.r_[np.nan, log_mean], np.r_[0.0, log_cv2])[0])


def test_noise_residual_ranks_overdispersed_genes_above_poisson_ones():
    rng = np.random.default_rng(3)
    means = np.exp(rng.uniform(0, 4, 600))
    x = rng.poisson(means, size=(500, 600)).astype(float)
    x[:, :60] = rng.negative_binomial(1, 1 / (1 + means[:60]), size=(500, 60))
    res, _, _ = noise_residual(x, np.repeat(["m1", "m2"], 250))
    assert np.nanmean(res[:60]) > np.nanmean(res[60:]) + 1


def test_strata_cross_mean_deciles_with_dropout_tertiles():
    rng = np.random.default_rng(4)
    mean = rng.uniform(0, 1, 3000)
    dropout = rng.uniform(0, 1, 3000)
    s = strata(mean, dropout)
    assert np.unique(s).size == 30
    top = s // 3 == 9
    assert mean[top].min() >= np.quantile(mean, 0.89)
    within = top & (s % 3 == 0)
    assert dropout[within].max() <= dropout[top & (s % 3 == 2)].min()


def test_matched_diff_removes_a_difference_carried_by_the_matching_variable():
    # members sit in high strata where everyone has higher values; matching must cancel that
    stratum = np.repeat([0, 1], 100)
    values = np.where(stratum == 1, 5.0, 0.0)
    member = np.zeros(200, bool)
    member[150:170] = True
    member[10:12] = True
    assert matched_diff(values, member, ~member, stratum) == pytest.approx(0.0)
    values[member] -= 1.0
    assert matched_diff(values, member, ~member, stratum) == pytest.approx(-1.0)
    values[0] = np.nan
    assert np.isfinite(matched_diff(values, member, ~member, stratum))


def test_gene_bootstrap_brackets_a_planted_gap():
    rng = np.random.default_rng(5)
    stratum = rng.integers(0, 5, 2000)
    member = rng.random(2000) < 0.1
    values = rng.normal(0, 1, 2000) + stratum - 0.5 * member
    boot = gene_boot_diffs(values, member, ~member, stratum, 200, 6)
    lo, hi = np.quantile(boot, [0.025, 0.975])
    assert lo < -0.5 < hi and hi < 0


def test_between_spread_is_unbiased_for_the_mean_absolute_deviation_scale():
    rng = np.random.default_rng(7)
    x = rng.normal(0, 1, size=(4, 20000))
    # E|X - Xbar| * sqrt(n/(n-1)) = sqrt(2/pi) for unit normals
    assert between_spread(x).mean() == pytest.approx(np.sqrt(2 / np.pi), rel=0.02)
    with pytest.raises(ValueError):
        between_spread(x[:1])


def test_mouse_draw_keeps_group_sizes_and_draws_only_its_own_mice():
    draw = mouse_draw({"young": ["a", "b", "c"], "old": ["d", "e"]}, np.random.default_rng(8))
    assert len(draw["young"]) == 3 and set(draw["young"]) <= {"a", "b", "c"}
    assert len(draw["old"]) == 2 and set(draw["old"]) <= {"d", "e"}


def test_eligible_needs_cells_and_mice_in_both_age_groups():
    def cells(ct, age, mice, per):
        return [
            {"tissue": "T", "cell_type": ct, "age_group": age, "mouse": f"{age}{m}"}
            for m in range(mice)
            for _ in range(per)
        ]

    rows = (
        cells("ok", "young", 3, 20)
        + cells("ok", "old", 4, 15)
        + cells("two_mice", "young", 2, 40)
        + cells("two_mice", "old", 4, 40)
        + cells("few_cells", "young", 3, 12)
        + cells("few_cells", "old", 3, 30)
        + cells("one_age", "young", 5, 30)
    )
    out = eligible(pl.DataFrame(rows))
    assert out["cell_type"].to_list() == ["ok"]


def test_thinning_at_the_rate_of_all_reads_keeps_the_share_of_the_given_columns():
    rng = np.random.default_rng(10)
    x = sparse.csr_matrix(rng.poisson(50, size=(300, 100)).astype(np.float32))
    all_reads = np.asarray(x.sum(axis=1)).ravel() * 2  # the columns hold half of each cell's reads
    out = thin_rows(x, int(all_reads.min() // 2), rng, totals=all_reads)
    kept = np.asarray(out.sum(axis=1)).ravel() / np.asarray(x.sum(axis=1)).ravel()
    assert kept.mean() == pytest.approx((all_reads.min() // 2) / all_reads.mean(), rel=0.05)


def test_mouse_weights_reproduce_and_drop_mice():
    rng = np.random.default_rng(11)
    x = rng.poisson(5, size=(300, 50)).astype(float)
    mouse = np.repeat(["a", "b", "c"], 100)
    plain = within_moments(x, mouse)
    ones = within_moments(x, mouse, {"a": 1.0, "b": 1.0, "c": 1.0})
    for p, o in zip(plain, ones, strict=True):
        assert np.allclose(p, o)
    no_c = within_moments(x, mouse, {"a": 1.0, "b": 1.0, "c": 0.0})
    two = within_moments(x[:200], mouse[:200])
    for p, o in zip(two, no_c, strict=True):
        assert np.allclose(p, o)


def test_mouse_weights_are_positive_with_mean_one_per_age_group():
    w = mouse_weights({"young": ["a", "b", "c"], "old": ["d", "e"]}, np.random.default_rng(12))
    assert set(w) == {"a", "b", "c", "d", "e"} and min(w.values()) > 0
    assert sum(w[m] for m in "abc") == pytest.approx(3) and sum(w[m] for m in "de") == pytest.approx(2)


def test_a_mouse_drawn_twice_counts_once():
    assert not distinct_mice_ok(["a", "a", "a"], {"a", "b"})
    assert distinct_mice_ok(["a", "a", "b"], {"a", "b"})
    assert not distinct_mice_ok(["a", "c", "c"], {"a", "b"})  # c gave no cells here


def test_strata_with_ties_keep_equal_dropout_in_one_bin():
    mean = np.arange(300, dtype=float)
    dropout = np.zeros(300)
    dropout[::2] = 0.5
    s = strata(mean, dropout, n_mean=1, ties=True)
    assert np.unique(s[dropout == 0]).size == 1 and np.unique(s[dropout == 0.5]).size == 1
    assert np.unique(strata(mean, dropout, n_mean=1)[dropout == 0]).size > 1  # input order splits them


def test_planted_noise_raises_log_cv2_by_the_shift_and_keeps_means_and_zeros():
    rng = np.random.default_rng(13)
    means = np.exp(rng.uniform(1, 4, 400))
    x = rng.negative_binomial(5, 5 / (5 + means), size=(2000, 400)).astype(float)
    mouse = np.repeat(["a", "b", "c", "d"], 500)
    rows = np.zeros(2000, bool)
    rows[:1000] = True
    cols = np.zeros(400, bool)
    cols[:200] = True
    out = plant_noise(x, rows, cols, 0.3, mouse, np.random.default_rng(14))
    assert np.array_equal(out == 0, x == 0)
    assert np.array_equal(out[~rows], x[~rows]) and np.array_equal(out[:, ~cols], x[:, ~cols])

    def log_cv2(a, m):
        mean, var, _ = within_moments(a, m)
        return np.log(var / mean**2)

    rise = log_cv2(out[rows], mouse[rows]) - log_cv2(x[rows], mouse[rows])
    assert np.mean(rise[:200]) == pytest.approx(0.3, abs=0.03)
    assert np.abs(rise[200:]).max() == 0
    assert np.allclose(out[rows][:, :200].mean(axis=0), x[rows][:, :200].mean(axis=0), rtol=0.05)


def test_one_to_one_orthologs_drop_classes_with_paralogs():
    rows = [
        ("1", "mouse", "Aaa", 1),
        ("1", "human", "AAA", 11),
        ("2", "mouse", "Bbb1", 2),
        ("2", "mouse", "Bbb2", 3),
        ("2", "human", "BBB", 12),
        ("3", "mouse", "Ccc", 4),
    ]
    frame = pl.DataFrame(rows, schema=["class_key", "organism", "symbol", "entrez"], orient="row")
    out = one_to_one_orthologs(frame)
    assert out["mouse_symbol"].to_list() == ["Aaa"] and out["human_entrez"].to_list() == [11]


def test_matched_count_skips_members_without_background_in_their_stratum():
    stratum = np.array([0, 0, 1, 1, 2])
    member = np.array([True, False, True, False, True])
    background = ~member
    values = np.array([1.0, 2.0, 3.0, np.nan, 5.0])
    assert matched_count(values, member, background, stratum) == 1


def test_pooling_excludes_groups_without_a_finite_gap_and_needs_them_all_in_each_resample():
    # regression (Codex review 8): a group with no retained non-essential genes made the pooled row NaN
    gaps = pl.DataFrame(
        {
            "set": ["s"] * 3,
            "tissue": ["T"] * 3,
            "cell_type": ["a", "b", "c"],
            "sd": [1.0, 1.0, 1.0],
            "gap_sd": [0.1, 0.3, float("nan")],
        }
    )
    boot = pl.DataFrame(
        {
            "set": ["s"] * 6,
            "tissue": ["T"] * 6,
            "cell_type": ["a", "b", "c", "a", "b", "a"],
            "boot": [0, 0, 0, 1, 1, 2],
            "matched_diff": [0.0, 0.2, float("nan"), 0.4, float("nan"), 0.1],
        }
    )
    out = pool_gaps(gaps, boot, ["set"]).row(0, named=True)
    assert out["cell_types"] == 2 and out["cell_types_excluded"] == 1
    assert out["mean_gap_sd"] == pytest.approx(0.2)
    assert out["complete_resamples"] == 1 and out["ci_low"] == pytest.approx(0.1)
