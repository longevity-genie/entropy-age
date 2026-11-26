import numpy as np
import polars as pl

from entropy_age.experiments.clock_drift_controls import (
    clock_stats,
    corr,
    export_check,
    flag_primary,
    gallery_cpgs,
    gallery_inputs,
    gene_balanced_pair,
    matched_arrays,
    matched_placebo,
    matched_subset,
    missing_estimates,
    paired_weights,
    panel_indices,
    panel_overlap,
    partial_r,
    partial_r_checked,
    per_gene_counts,
    placebo_band,
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


# ============================================================================= review 5
# ----------------------------------------------------------------------------- finding 3: both targets by gene
def _genes(start: int, sizes: list[int]) -> list[np.ndarray]:
    out, j = [], start
    for k in sizes:
        out.append(np.arange(j, j + k))
        j += k
    return out


def test_gene_balanced_pair_takes_one_cpg_per_gene_on_both_sides():
    rng = np.random.default_rng(20)
    hk_genes = _genes(0, [6] * 40)  # 240 housekeeping CpGs in 40 clustered genes
    pool = _genes(240, [5] * 200)  # 1000 placebo candidates in 200 genes
    strata = rng.choice(np.array(["a", "b", "c"]), 1240)
    sub, placebo, background, fill = gene_balanced_pair(strata, hk_genes, pool, np.random.default_rng(1))
    assert per_gene_counts(sub, hk_genes).max() == 1 and per_gene_counts(placebo, pool).max() == 1
    assert per_gene_counts(sub, hk_genes).size == sub.sum() == placebo.sum() == per_gene_counts(placebo, pool).size
    for s_ in "abc":  # the common support: equal counts in every stratum
        assert (sub & (strata == s_)).sum() == (placebo & (strata == s_)).sum()
    assert sub.sum() <= 40 and 0 < fill <= 1
    hk_all = np.zeros(1240, dtype=bool)
    hk_all[:240] = True
    assert not (background & (hk_all | placebo)).any() and (background | hk_all | placebo).all()


def test_placebo_band_reports_the_margin_on_common_support_per_scheme():
    p = pl.DataFrame(
        {
            "series": ["s"] * 4,
            "group": ["g"] * 4,
            "clock": ["epiTOC2"] * 4,
            "hk_r_raw": [-0.5, -0.6, -0.4, float("nan")],
            "hk_r_adj": [-0.5, -0.6, -0.4, -0.5],
            "r_raw": [-0.4, -0.5, -0.5, -0.3],
            "r_adj": [-0.4, -0.5, -0.5, float("nan")],
            "balance_tvd": [0.0] * 4,
            "balance_unusable": [0.0] * 4,
            "adj_note": ["", "", "", "x"],
        }
    )
    band = placebo_band(p, suffix="_gb").row(0, named=True)
    assert np.isclose(band["margin_adj_gb"], 0.1)  # median of 0.1, 0.1, -0.1: housekeeping ahead by 0.1
    assert "hk_matched_adj_gb" in band and "margin_adj" not in band


# ----------------------------------------------------------------------------- finding 4: preprocessing inputs
def _fake_biolearn(tmp_path):
    pkg = tmp_path / "biolearn"
    (pkg / "data").mkdir(parents=True)
    (pkg / "model.py").write_text(
        "from biolearn.norm import normalise\n"
        "model_definitions = {\n"
        '    "Plain": {\n        "model": {"type": "LinearMethylationModel", "file": "Plain.csv"},\n    },\n'
        '    "Paced": {\n        "model": {\n            "type": "LinearMethylationModel",\n'
        '            "file": "Paced.csv",\n            "preprocess": normalise,\n        },\n    },\n'
        "}\n"
    )
    (pkg / "norm.py").write_text(
        "from biolearn.util import get_data_file\n"
        "def normalise(df):\n    return _gold(df)\n"
        "def _gold(df):\n    return get_data_file('Gold.csv')\n"
    )
    (pkg / "data" / "Plain.csv").write_text("CpGmarker,CoefficientTraining\ncg001,1\nintercept,2\n")
    (pkg / "data" / "Paced.csv").write_text("CpGmarker,CoefficientTraining\ncg002,1\n")
    (pkg / "data" / "Gold.csv").write_text('CpGmarker,mean\n"cg100",0.5\n"cg101",0.4\n')
    return pkg


def test_gallery_exclusion_holds_preprocessing_probe_lists(tmp_path, monkeypatch):
    pkg = _fake_biolearn(tmp_path)
    monkeypatch.setenv("ENTROPY_AGE_BIOLEARN_SRC", str(pkg))
    assert gallery_inputs(pkg)["Paced"] == ["Gold.csv", "Paced.csv"]
    assert gallery_cpgs(tmp_path, {"Plain", "Paced"}) == {"cg001", "cg002", "cg100", "cg101"}
    panel = pl.DataFrame({"series": ["s"] * 3, "group": ["g"] * 3, "rank": [0, 1, 2], "cpg": ["cg100", "cg5", "cg6"]})
    assert panel_overlap(panel, gallery_cpgs(tmp_path, {"Paced"})) == 1


# ----------------------------------------------------------------------------- finding 5: registry before predictions
def test_primary_group_fixed_by_the_registry_and_missing_clocks_reported():
    registry = pl.DataFrame(
        {
            "accession": ["G1", "G1", "G2"],
            "series": ["s1", "s1", "s2"],
            "group": ["all", "sub", "x"],
            "donors": [100, 50, 20],
            "eligible": [True, True, False],
            "sensitivity": [False, False, False],
            "brain": [False, False, False],
        }
    )
    table = pl.DataFrame(  # the largest group emitted no clock row
        {
            "accession": ["G1"],
            "series": ["s1"],
            "group": ["sub"],
            "donors": [50],
            "sensitivity": [False],
            "clock": ["epiTOC2"],
        }
    )
    assert flag_primary(table, registry)["primary"].to_list() == [False]
    missing = missing_estimates(registry, table, ["epiTOC2", "HypoClock"], {("s1", "all", "*"): "no clock scores"})
    assert missing.select("series", "group", "clock", "reason").rows() == [
        ("s1", "all", "epiTOC2", "no clock scores"),
        ("s1", "all", "HypoClock", "no clock scores"),
    ]  # G2 is not eligible, so it has no primary group


# ----------------------------------------------------------------------------- finding 6: one scaled design
def test_partial_r_fits_with_the_scaled_design_it_ranks():
    rng = np.random.default_rng(21)
    t = rng.normal(size=80)
    given = np.column_stack([rng.normal(size=80), 1e-16 * t])  # a tiny column the scaled rank counts
    r, why = partial_r_checked(t, -3 * t + 1e-3 * rng.normal(size=80) * 1e-14, given)
    assert np.isnan(r) and "explained by the design" in why


def test_clock_linear_in_age_gives_no_correlation():
    rng = np.random.default_rng(22)
    age = rng.uniform(20, 80, 60)
    v = 0.7 * age + 3.0  # a clock that is age, exactly
    f = rng.normal(size=60)
    design = np.column_stack([age, rng.normal(size=(60, 2))])
    ind = {k: rng.normal(size=60) for k in ("level", "compression", "sd")}
    st = clock_stats(v, age, design, f, rng.normal(size=60), rng.normal(size=60), ind)
    assert st["age_note"] == "clock explained by age" and np.isnan(st["r_F_raw"])
    assert np.isnan(st["r_F_adj"]) and "clock explained" in st["r_F_adj_note"]
    noisy = clock_stats(v + rng.normal(0, 1, 60), age, design, f, f, rng.normal(size=60), ind)
    assert noisy["age_note"] == "" and np.isfinite(noisy["r_F_adj"]) and np.isfinite(noisy["r_F_raw"])


# ----------------------------------------------------------------------------- finding 7: stale exports
def test_export_check_flags_an_export_older_than_the_rebuild():
    ids = {"a", "b", "c"}
    kw = dict(export_samples=ids, clock_samples=ids, passing=ids)
    assert export_check(qc_mtime=10, export_mtime=20, clocks_mtime=30, **kw) == (
        set(),
        "no manifest; times and ids agree",
        False,
    )
    drop, note, stale = export_check(qc_mtime=25, export_mtime=20, clocks_mtime=30, **kw)
    assert stale and "older than the clock-track QC" in note and not drop
    _, note, stale = export_check(
        qc_mtime=10, export_mtime=20, clocks_mtime=30, export_samples=ids, clock_samples=ids, passing=ids | {"d"}
    )
    assert stale and "1 passing arrays absent" in note


def test_export_check_drops_samples_whose_exported_array_is_not_the_current_one():
    current = pl.DataFrame(
        {"sample": ["a", "b"], "chip": ["100", "100"], "position": ["R01C01", "R02C01"], "prefix": ["pa", "pb"]}
    )
    manifest = pl.DataFrame({"sample": ["a", "b"], "array": ["100_R01C01", "200_R07C01"]})  # b: an earlier repeat
    drop, note, stale = export_check(
        qc_mtime=99,
        export_mtime=1,
        clocks_mtime=2,
        export_samples={"a", "b"},
        clock_samples={"a", "b"},
        passing={"a", "b"},
        manifest=manifest,
        current=current,
    )
    assert drop == {"b"} and not stale and note.startswith("manifest: 1 arrays differ")
