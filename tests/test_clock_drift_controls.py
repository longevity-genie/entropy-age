import numpy as np
import polars as pl

from entropy_age.experiments.clock_drift_controls import (
    corr,
    flag_primary,
    matched_arrays,
    matched_placebo,
    matched_subset,
    paired_weights,
    panel_indices,
    partial_r,
    partial_r_checked,
    placebo_pair,
    placebo_pool,
    primary_groups,
    seeded_panel,
    summarise,
)


def test_panel_indices_level_compression_and_spread():
    x = np.array([[0.5, 0.5, 0.5, 0.5], [0.0, 1.0, 0.0, 1.0], [0.1, 0.1, 0.9, 0.9]])
    ind = panel_indices(x)
    assert np.allclose(ind["level"], [0.5, 0.5, 0.5])
    assert np.allclose(ind["compression"], [0.0, 0.5, 0.4])  # mean |beta - 0.5|: lower is more compressed
    assert np.allclose(ind["sd"], [0.0, 0.5, 0.4])


def test_compression_falls_when_betas_shrink_toward_half():
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1, (1, 5000))
    shrunk = 0.5 + 0.7 * (x - 0.5)
    assert panel_indices(shrunk)["compression"][0] < panel_indices(x)["compression"][0]
    assert np.isclose(panel_indices(shrunk)["level"][0], 0.5 + 0.7 * (panel_indices(x)["level"][0] - 0.5))


def test_partial_r_removes_a_shared_driver():
    rng = np.random.default_rng(2)
    z = rng.normal(size=2000)
    x = z + rng.normal(0, 0.3, z.size)
    y = z + rng.normal(0, 0.3, z.size)
    assert corr(x, y) > 0.8
    assert abs(partial_r(x, y, z[:, None])) < 0.1


def test_partial_r_keeps_a_direct_link():
    rng = np.random.default_rng(3)
    z = rng.normal(size=2000)
    x = rng.normal(size=z.size)
    y = 0.8 * x + z + rng.normal(0, 0.3, z.size)
    assert partial_r(x, y, z[:, None]) > 0.9


def test_matched_placebo_reproduces_stratum_counts():
    rng = np.random.default_rng(4)
    strata = rng.choice(np.array(["a", "b", "c"]), 600)
    real = np.zeros(600, dtype=bool)
    real[:60] = True
    pool = [np.arange(j, min(j + 4, 600)) for j in range(60, 600, 4)]
    mask, genes, fill = matched_placebo(strata, real, pool, np.random.default_rng(0))
    assert fill == 1.0
    assert not (mask & real).any()
    for s in "abc":
        assert (mask & (strata == s)).sum() == (real & (strata == s)).sum()
    assert sum(g.size for g in genes) == mask.sum()


def test_placebo_pool_skips_set_genes_and_members():
    kept = ["c1", "c2", "c3", "c4"]
    links = pl.DataFrame(
        {
            "cpg": ["c1", "c2", "c3", "c4", "c4"],
            "symbol": ["HK", "G1", "G1", "G2", "HK"],
            "promoter": [True, True, True, True, True],
        }
    )
    member = np.array([True, False, True, False])
    pool = placebo_pool(kept, member, links, {"HK"})
    assert [g.tolist() for g in pool] == [[1], [3]]


def test_loo_drift_matches_loo_scores():
    from entropy_age.experiments.clock_drift_controls import loo_drift
    from entropy_age.experiments.health_states import loo_scores

    rng = np.random.default_rng(5)
    x = np.clip(rng.beta(2, 5, (12, 400)) + rng.normal(0, 0.02, (12, 400)), 0, 1)
    island = rng.random(400) < 0.5
    a = np.zeros(400, dtype=bool)
    a[:40] = True
    b = np.zeros(400, dtype=bool)
    b[40:90] = True
    sets = {"hk": (a, [np.arange(40)]), "other": (b, [np.arange(40, 90)])}
    ours, theirs = loo_drift(x, island, sets), loo_scores(x, island, sets)
    for key in ("P_hk", "B_hk", "F_hk", "F_other"):
        assert np.allclose(ours[key], theirs[key])


def test_loo_drift_drop_leaves_cpgs_out_of_the_contrast():
    from entropy_age.experiments.clock_drift_controls import loo_drift

    rng = np.random.default_rng(6)
    x = rng.beta(2, 5, (10, 300))
    island = np.zeros(300, dtype=bool)
    a = np.zeros(300, dtype=bool)
    a[:30] = True
    hk = np.zeros(300, dtype=bool)
    hk[30:80] = True
    x[:, 30:80] += rng.normal(0, 0.2, (10, 50))  # noisy CpGs: in the look-alikes they raise B
    sets = {"pl": (a, [np.arange(30)]), "hk": (hk, [np.arange(30, 80)])}
    kept = loo_drift(x, island, sets)
    dropped = loo_drift(x, island, sets, {"pl": hk})
    assert not np.allclose(kept["B_pl"], dropped["B_pl"])
    assert np.allclose(kept["F_hk"], dropped["F_hk"])
    assert np.allclose(loo_drift(x, island, sets, {"pl": np.zeros(300, dtype=bool)})["F_pl"], kept["F_pl"])


# ----------------------------------------------------------------------------- review 4, finding 1
def test_matched_subset_takes_the_placebos_attained_quotas():
    rng = np.random.default_rng(7)
    strata = rng.choice(np.array(["a", "b", "c", "d"]), 800)
    real = np.zeros(800, dtype=bool)
    real[:200] = True
    pool = [np.arange(j, min(j + 3, 800)) for j in range(200, 260, 3)]  # a small pool: some strata stay short
    placebo, _, fill = matched_placebo(strata, real, pool, np.random.default_rng(1))
    assert fill < 1.0
    sub = matched_subset(strata, real, placebo, np.random.default_rng(2))
    assert not (sub & ~real).any()
    for s in "abcd":
        assert (sub & (strata == s)).sum() == (placebo & (strata == s)).sum()
    assert sub.sum() == placebo.sum() < real.sum()


def test_paired_weights_share_one_background_and_one_set_of_stratum_weights():
    rng = np.random.default_rng(8)
    n = 600
    strata = rng.integers(0, 5, n)
    hk = np.zeros(n, dtype=bool)
    hk[:150] = True
    a = hk & (rng.random(n) < 0.4)  # the matched housekeeping subset
    b = np.zeros(n, dtype=bool)
    b[150 : 150 + int(a.sum())] = True  # the placebo
    strata[150 : 150 + int(a.sum())] = strata[np.flatnonzero(a)]  # the same count per stratum
    background = ~(hk | b)
    pa, pb, bw, bal = paired_weights(a, b, background, strata)
    assert np.all(bw[hk | b] == 0)  # the rest of the housekeeping set is in neither the sets nor the background
    assert np.all(pa[b] == 0) and np.all(pb[a] == 0)
    assert np.allclose(pa[background], pb[background])  # one background with the same weights for both
    for k in range(5):
        assert np.isclose(pa[a & (strata == k)].sum(), pb[b & (strata == k)].sum())  # identical stratum weights
    assert bal["tvd"] == 0 and bal["unusable"] == 0
    d = rng.random(n)
    d[a] = d[b] = 0.7  # the two sets deviate alike: their F must agree, whatever the background
    assert np.isclose(d @ pa / (d @ bw), d @ pb / (d @ bw))


def test_placebo_pair_background_leaves_out_every_housekeeping_and_placebo_cpg():
    rng = np.random.default_rng(13)
    strata = rng.choice(np.array(["a", "b"]), 400)
    hk = np.zeros(400, dtype=bool)
    hk[:80] = True
    pool = [np.arange(j, j + 2) for j in range(100, 160, 2)]
    sub, placebo, background, fill = placebo_pair(strata, hk, pool, np.random.default_rng(3))
    assert 0 < fill < 1 and sub.sum() == placebo.sum() and not (sub & ~hk).any()
    assert not (background & (hk | placebo)).any()
    assert (background | hk | placebo).all()


def test_paired_weights_pool_unequal_stratum_counts():
    strata = np.array([0] * 10 + [1] * 5 + [0] * 5 + [1] * 10 + [0] * 20 + [1] * 20)
    a = np.zeros(70, dtype=bool)
    a[:15] = True
    b = np.zeros(70, dtype=bool)
    b[15:30] = True
    pa, pb, bw, bal = paired_weights(a, b, ~(a | b), strata)
    for k in (0, 1):  # each set's weight in a stratum is the pair's pooled share there, 15 of 30 in both strata
        assert np.isclose(pa[a & (strata == k)].sum(), 0.5) and np.isclose(pb[b & (strata == k)].sum(), 0.5)
        assert np.isclose(bw[strata == k].sum(), 0.5)
    assert np.isclose(bal["tvd"], 1 / 3)


def test_loo_drift_pairs_share_b():
    from entropy_age.experiments.clock_drift_controls import loo_drift

    rng = np.random.default_rng(9)
    x = rng.beta(2, 5, (10, 300))
    island = rng.random(300) < 0.5
    hk = np.zeros(300, dtype=bool)
    hk[:60] = True
    sub = np.zeros(300, dtype=bool)
    sub[:30] = True
    pl_ = np.zeros(300, dtype=bool)
    pl_[100:130] = True
    out = loo_drift(x, island, {"hk": (hk, [np.arange(60)])}, pairs={"p": (sub, pl_, ~(hk | pl_))})
    assert np.allclose(out["F_p_a"] * out["B_p"], out["P_p_a"])
    assert np.allclose(out["F_p_b"] * out["B_p"], out["P_p_b"])
    assert not np.allclose(out["F_p_a"], out["F_p_b"])


# ----------------------------------------------------------------------------- finding 3
def test_seeded_panel_keeps_drawing_until_full():
    pool = np.array([f"cg{i:05d}" for i in range(10_000)])
    complete = lambda batch: {c for c in batch if int(c[2:]) % 4 == 0}  # noqa: E731  only a quarter complete
    picked = seeded_panel(pool, np.random.default_rng(0), complete, size=2000)
    assert len(picked) == 2000 and all(int(c[2:]) % 4 == 0 for c in picked)
    order = pool[np.random.default_rng(0).permutation(pool.size)]
    assert picked == [c for c in order if int(c[2:]) % 4 == 0][:2000]  # the seeded order, first come first kept
    assert len(seeded_panel(pool, np.random.default_rng(0), complete, size=5000)) == 2500  # pool exhausted


# ----------------------------------------------------------------------------- finding 5
def test_primary_group_per_accession_by_eligible_donors():
    t = pl.DataFrame(
        {
            "accession": ["G1", "G1", "G1", "G2", "G2", "G3", "G3"],
            "series": ["s1", "s1", "s1b", "s2", "s2_no_field_d10", "s3", "s3"],
            "group": ["all", "controls", "fibrosis", "x", "x", "every donor", "cases"],
            "donors": [100, 40, 120, 50, 60, 70, 70],
            "n": [100, 40, 30, 50, 60, 70, 70],  # clock-scored n must not decide
            "sensitivity": [False, False, False, False, True, False, False],
        }
    )
    assert sorted(primary_groups(t).rows()) == [("s1b", "fibrosis"), ("s2", "x"), ("s3", "every donor")]


def test_summary_over_primary_groups_and_accession_medians():
    t = pl.DataFrame(
        {
            "accession": ["G1", "G1", "G2"],
            "series": ["s1", "s1b", "s2"],
            "group": ["a", "b", "c"],
            "donors": [100, 50, 80],
            "sensitivity": [False, False, False],
            "brain": [False, False, False],
            "clock": ["epiTOC2"] * 3,
            "mitotic": [True] * 3,
            "r_F_raw": [-0.6, -0.2, -0.4],
            "r_F_adj": [-0.5, -0.1, float("nan")],
            "hk_matched_adj": [-0.4, -0.1, -0.3],
            "placebo_abs_p95_adj": [0.2, 0.2, 0.2],
            "exceeds_p95_raw": [True, False, True],
            "exceeds_p95_adj": [True, False, True],
        }
    )
    flagged = flag_primary(t)
    assert flagged.height == 3 and flagged["primary"].to_list() == [True, False, True]
    s = summarise(flagged)
    prim = s.filter(pl.col("level") == "primary groups").row(0, named=True)
    assert prim["groups"] == 2 and prim["exceeding_adj"] == 2 and prim["undefined_adj"] == 1
    assert np.isclose(prim["r_F_raw"], -0.5)
    acc = s.filter(pl.col("level") == "accession medians").row(0, named=True)
    assert acc["groups"] == 2 and np.isclose(acc["r_F_raw"], np.median([-0.4, -0.4]))


# ----------------------------------------------------------------------------- finding 6
def test_partial_r_undefined_when_the_design_explains_a_vector():
    rng = np.random.default_rng(10)
    z = rng.normal(size=(60, 3)) * np.array([50.0, 1.0, 1e-3])  # columns on very different scales
    x = z @ np.array([0.02, 1.0, 300.0])
    y = z @ np.array([-0.01, 2.0, 100.0])
    r, why = partial_r_checked(x, y, z, ("acceleration", "F"))
    assert np.isnan(r) and "explained by the design" in why
    assert np.isnan(partial_r(x, y, z))


def test_partial_r_undefined_with_too_few_residual_df():
    rng = np.random.default_rng(11)
    z = rng.normal(size=(14, 5))
    x, y = rng.normal(size=14), rng.normal(size=14)
    r, why = partial_r_checked(x, y, z)
    assert np.isnan(r) and why.startswith("residual df 8 < 10")


def test_partial_r_counts_the_rank_of_a_degenerate_design():
    rng = np.random.default_rng(12)
    chip = np.repeat(np.eye(3), 10, axis=0)  # dummies of three chips on 30 rows, with the intercept: one too many
    z = np.column_stack([chip, chip[:, :1]])  # and a duplicated column
    x = rng.normal(size=30)
    y = x + rng.normal(0, 0.5, 30)
    r, why = partial_r_checked(x, y, z)
    assert why == "" and np.isfinite(r)  # rank 3, df 27: defined despite the duplicated columns


# ----------------------------------------------------------------------------- finding 7
def test_matched_arrays_drops_samples_read_from_another_array():
    series_qc = pl.DataFrame(
        {
            "prefix": ["GSM1_100_R01C01", "GSM2_100_R02C01", "u3_noid", "GSM4_100_R04C01"],
            "sample": ["GSM1", "GSM2", "u3_noid", "GSM4"],
            "chip": ["100", "100", "NA", "100"],
            "position": ["R01C01", "R02C01", "NA", "R04C01"],
        }
    )
    clock_qc = pl.DataFrame(
        {
            "prefix": ["GSM1_100_R01C01", "GSM2_200_R05C01", "u3_noid"],
            "sample": ["GSM1", "GSM2", "u3_noid"],
            "chip": ["100", "200", "NA"],
            "position": ["R01C01", "R05C01", "NA"],
        }
    )
    same, reasons = matched_arrays(series_qc, clock_qc)
    assert same == {"GSM1", "u3_noid"}  # GSM4 has no clock-track QC row, so it is unavailable too
    assert len(reasons) == 1 and reasons[0].startswith("GSM2: drift array 100_R02C01, clock array 200_R05C01")
