import numpy as np

from entropy_age.experiments.clock_drift_classes import (
    class_draw,
    class_scores,
    clusters,
    contrasts,
    loo_codes,
    partial_rs,
    random_effects,
    shared_weights,
    thin,
    units,
)
from entropy_age.experiments.clock_drift_controls import paired_weights, partial_r_checked, per_gene_counts
from entropy_age.protection import fit_reference


def test_clusters_join_neighbours_under_the_gap():
    chrom = np.array(["1", "1", "1", "2", "2", "1"])
    pos = np.array([100, 900, 5000, 100, 200, 5500])
    lab = clusters(chrom, pos, gap=1000)
    assert lab[0] == lab[1] != lab[2] == lab[5]
    assert lab[3] == lab[4] and lab[3] not in (lab[0], lab[2])


def test_units_group_members_by_label():
    out = units(np.array([7, 7, 3, 3, 9]), np.array([0, 1, 3, 4]))
    assert [u.tolist() for u in out] == [[3], [0, 1], [4]]


def test_thin_reaches_the_wanted_counts():
    strata = np.array([0, 0, 0, 1, 1, 2])
    mask = np.array([True, True, True, True, True, False])
    out = thin(strata, mask, np.array([2, 1, 0]), np.random.default_rng(0))
    assert out.sum() == 3 and not (out & ~mask).any() and (out & (strata == 0)).sum() == 2


def test_class_draw_gives_three_disjoint_classes_with_identical_strata_one_per_unit():
    rng = np.random.default_rng(1)
    n = 3000
    strata = rng.choice(np.array(["a", "b", "c", "d"]), n)
    hk_genes = [np.arange(j, j + 5) for j in range(0, 300, 5)]  # 60 genes
    op_genes = [np.arange(j, j + 4) for j in range(300, 1500, 4)]
    np_units = [np.arange(j, j + 3) for j in range(1500, 3000, 3)]
    m = class_draw(strata, hk_genes, op_genes, np_units, np.random.default_rng(2))
    assert not (m["hk"] & m["op"]).any() and not (m["op"] & m["np"]).any() and not (m["hk"] & m["np"]).any()
    for s in "abcd":
        counts = {c: int((m[c] & (strata == s)).sum()) for c in m}
        assert len(set(counts.values())) == 1
    assert per_gene_counts(m["hk"], hk_genes).max() == 1
    assert per_gene_counts(m["op"], op_genes).max() == 1 and per_gene_counts(m["np"], np_units).max() == 1
    assert 0 < m["hk"].sum() <= 60 and not m["op"][:300].any() and m["np"][1500:].sum() == m["np"].sum()


def test_shared_weights_match_paired_weights_for_two_sets():
    rng = np.random.default_rng(3)
    codes = rng.integers(0, 4, 500)
    a = np.zeros(500, dtype=bool)
    a[:40] = True
    b = np.zeros(500, dtype=bool)
    b[40:90] = True
    bg = ~(a | b)
    (wa, wb), bw, bal = shared_weights([a, b], bg, codes)
    pa, pb, bw2, bal2 = paired_weights(a, b, bg, codes)
    assert np.allclose(bw, bw2) and np.allclose(wa, pa + bw2) and np.allclose(wb, pb + bw2)
    assert np.isclose(bal["tvd"], bal2["tvd"]) and np.isclose(wa.sum(), 1) and np.isclose(bw.sum(), 1)


def test_shared_weights_equal_scores_for_classes_that_deviate_alike():
    rng = np.random.default_rng(4)
    codes = np.repeat(np.arange(3), 100)
    sets = [np.zeros(300, dtype=bool) for _ in range(3)]
    for i, s in enumerate(sets):
        s[i * 10 : i * 10 + 10] = s[100 + i * 10 : 100 + i * 10 + 10] = s[200 + i * 10 : 200 + i * 10 + 10] = True
    bg = ~np.any(sets, axis=0)
    ws, bw, bal = shared_weights(sets, bg, codes)
    d = rng.random(300)
    d[np.any(sets, axis=0)] = np.tile([0.2, 0.2, 0.5], 1)[codes[np.any(sets, axis=0)]]  # deviation set by stratum
    scores = [d @ w for w in ws]
    assert np.allclose(scores, scores[0]) and bal["tvd"] == 0 and bal["unusable"] == 0


def test_loo_codes_give_fit_references_partition():
    rng = np.random.default_rng(5)
    x = rng.beta(2, 5, (8, 400))
    island = rng.random(400) < 0.4
    ref = fit_reference(x[1:], island, np.ones(400, dtype=bool))
    ours = loo_codes(x[1:].mean(axis=0), island)
    pairs = set(zip(ref.strata.tolist(), ours.tolist(), strict=True))
    assert len(pairs) == len(set(ref.strata.tolist())) == len(set(ours.tolist()))


def test_class_scores_use_the_leave_one_out_mean():
    rng = np.random.default_rng(6)
    x = rng.beta(2, 5, (10, 300))
    island = rng.random(300) < 0.5
    s1 = np.zeros(300, dtype=bool)
    s1[:30] = True
    s2 = np.zeros(300, dtype=bool)
    s2[30:60] = True
    bg = ~(s1 | s2)
    a, b, bal = class_scores(x, island, [([s1, s2], bg)])
    from entropy_age.protection import deviations

    level = x[1:].mean(axis=0)
    (w1, _), bw, _ = shared_weights([s1, s2], bg, loo_codes(level, island))
    d = deviations(x[0], level)
    assert np.isclose(a[0, 0, 0], d @ w1) and np.isclose(b[0, 0], d @ bw)


def test_partial_rs_agree_with_partial_r_checked():
    rng = np.random.default_rng(7)
    g = rng.normal(size=(80, 4))
    v = g @ rng.normal(size=4) + rng.normal(size=80)
    t = np.column_stack([v + rng.normal(size=80), rng.normal(size=80), g[:, 0] * 2.0])
    r, note = partial_rs(v, t, g)
    assert note == ""
    for j in range(2):
        assert np.isclose(r[j], partial_r_checked(v, t[:, j], g)[0])
    assert np.isnan(r[2])  # a target the design explains


def test_partial_rs_undefined_when_the_clock_is_explained():
    rng = np.random.default_rng(8)
    age = rng.uniform(20, 80, 50)
    r, note = partial_rs(2 * age + 1, rng.normal(size=(50, 3)), age[:, None])
    assert np.isnan(r).all() and "explained" in note


def test_contrasts_order():
    assert np.allclose(contrasts(np.array([-0.6, -0.5, -0.2])), [-0.3, -0.1])


def test_random_effects_pools_and_widens_with_heterogeneity():
    mu, se, tau2 = random_effects(np.array([0.1, 0.1, 0.1]), np.array([0.05, 0.05, 0.05]))
    assert np.isclose(mu, 0.1) and tau2 == 0 and np.isclose(se, 0.05 / np.sqrt(3))
    mu, se2, tau2 = random_effects(np.array([-0.3, 0.0, 0.3]), np.array([0.05, 0.05, 0.05]))
    assert np.isclose(mu, 0.0) and tau2 > 0 and se2 > se
