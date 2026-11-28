import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.solo_wcgw_culture import (
    island_class,
    joint_rates,
    matched_weights,
    reference_rows,
    score_unit,
    t_interval,
    within_slope,
)


def test_island_class_drops_the_side():
    assert [island_class(r) for r in ("N_Shore", "S_Shelf", "OpenSea", "Island")] == [
        "Shore",
        "Shelf",
        "OpenSea",
        "Island",
    ]


def test_matched_weights_reweight_the_pool_to_the_set_strata():
    strata = np.array(["a", "a", "a", "b", "b", "b"])
    count = np.array([1.0, 0, 0, 1, 1, 0])
    pool = np.array([True, True, True, True, True, True])  # a set CpG is never its own match
    s, m = matched_weights(count, pool, strata)
    assert s.sum() == pytest.approx(1) and m.sum() == pytest.approx(1)
    assert m[0] == 0 and m[3] == 0 and m[4] == 0
    # stratum a carries a third of the set, b two thirds: a's two pool CpGs share 1/3, b's one gets 2/3
    assert m[1] == pytest.approx(1 / 6) and m[2] == pytest.approx(1 / 6) and m[5] == pytest.approx(2 / 3)
    x = np.array([[5.0, 1, 1, 9, 9, 3]])
    assert (x @ s)[0] == pytest.approx((5 + 9 + 9) / 3)
    assert (x @ m)[0] == pytest.approx(1 / 3 + 2 / 3 * 3)


def test_matched_weights_drop_strata_without_pool_and_count_resampled_genes():
    strata = np.array(["a", "a", "b", "b"])
    s, m = matched_weights(np.array([2.0, 0, 1, 0]), np.array([False, True, False, False]), strata)
    assert s.tolist() == [1.0, 0, 0, 0] and m.tolist() == [0, 1.0, 0, 0]  # b has no pool: left out
    with pytest.raises(ValueError):
        matched_weights(np.array([1.0, 0]), np.array([False, False]), np.array(["a", "a"]))


def test_within_and_joint_slopes_recover_planted_rates():
    rng = np.random.default_rng(1)
    units = np.repeat(np.array(["u1", "u2", "u3"]), 20)
    doublings = rng.uniform(0, 5, 60)
    days = rng.uniform(0, 8, 60)
    offset = {"u1": 0.0, "u2": 3.0, "u3": -2.0}
    y = np.array([offset[u] for u in units]) + 0.7 * doublings + 0.0 * days + rng.normal(0, 0.01, 60)
    assert within_slope(y, doublings, units) == pytest.approx(0.7, abs=0.01)
    chips = np.tile(np.array(["c1", "c2"]), 30)
    assert within_slope(y, doublings, units, chips) == pytest.approx(0.7, abs=0.01)
    pd_rate, day_rate = joint_rates(y, doublings, days, units)
    assert pd_rate == pytest.approx(0.7, abs=0.01) and day_rate == pytest.approx(0.0, abs=0.01)


def test_t_interval_needs_three():
    assert np.isnan(t_interval(np.array([1.0, 2.0]))["low"])
    r = t_interval(np.array([1.0, 2.0, 3.0]))
    assert r["mean"] == 2.0 and r["low"] < 2.0 < r["high"] and r["n"] == 3


def test_reference_rows_use_days_then_doublings():
    with_days = pl.DataFrame({"days": [3.0, 3.0, 10.0], "doublings": [5.0, 4.0, 9.0]})
    assert reference_rows(with_days).tolist() == [0, 1]
    without = pl.DataFrame(
        {"days": [None, None, None], "doublings": [5.0, 4.0, 9.0]}, schema_overrides={"days": pl.Float64}
    )
    assert reference_rows(without).tolist() == [1]


def test_score_unit_finds_a_planted_shield():
    rng = np.random.default_rng(7)
    n = 3000
    hk = np.zeros(n, dtype=bool)
    hk[:300] = True
    promoter = np.zeros(n, dtype=bool)
    promoter[:1200] = True
    table = pl.DataFrame(
        {
            "cpg": [f"cg{i}" for i in range(n)],
            "island": rng.choice(["OpenSea", "Shore"], n),
            "hk": hk,
            "hk_gene": [f"g{i // 3}" if i < 300 else None for i in range(n)],
            "promoter": promoter,
            "non_essential": np.zeros(n, dtype=bool),
            "pmd": np.zeros(n, dtype=bool),
        }
    )
    reference = rng.uniform(0.3, 0.9, n)
    steps = np.arange(1, 5.0)
    per_step = np.where(hk, 0.5, 1.0)  # housekeeping CpGs lose half as much per step
    loss = np.outer(steps, per_step) + rng.normal(0, 0.05, (4, n))
    s = score_unit(loss, reference, table, 0.0)
    assert np.allclose(s["loss_hk_other"], 0.5 * steps, atol=0.05)
    assert np.allclose(s["matched_other"], steps, atol=0.05)
    assert np.allclose(s["gap_promoter"], -0.5 * steps, atol=0.05)
    assert np.isnan(s["loss_pmd"]).all() and np.isnan(s["gap_non_essential"]).all()
    high = score_unit(loss, reference, table, 0.6)  # the band keeps only CpGs starting at 0.6 or above
    assert np.allclose(high["gap_other"], -0.5 * steps, atol=0.05)
    counts = np.where(hk, 2.0, 0.0)  # every gene drawn twice: the same means
    assert np.allclose(score_unit(loss, reference, table, 0.0, counts)["gap_other"], s["gap_other"])


def test_confounded_exposure_gives_no_slope():
    units = np.array(["u"] * 6)
    chips = np.array(["a", "a", "b", "b", "c", "c"])
    t = np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0])  # doublings fixed by chip
    y = t + np.array([0.0, 0.1, 0.0, 0.1, 0.0, 0.1])
    assert np.isfinite(within_slope(y, t, units))
    assert np.isnan(within_slope(y, t, units, chips))
    pd_rate, day_rate = joint_rates(y, t, 2 * t, units)  # days a multiple of doublings: neither identified
    assert np.isnan(pd_rate) and np.isnan(day_rate)


def test_an_undrawn_housekeeping_gene_stays_out_of_both_pools():
    strata = np.array(["a"] * 6)
    hk = np.array([True, True, False, False, False, False])
    counts = np.array([1.0, 0, 0, 0, 0, 0])  # the second gene was not drawn
    s, m = matched_weights(counts, ~hk, strata)
    assert s[1] == 0 and m[1] == 0 and m[2:].sum() == pytest.approx(1)


def test_arrest_pairs_on_oxygen_and_skips_regrowth():
    from entropy_age.experiments.solo_wcgw_culture import arrest_contrast

    def row(arm, arrested, slope):
        return {
            "series": "s",
            "unit": f"c | {arm}",
            "culture": "c",
            "arm": arm,
            "arrested": arrested,
            "doublings_in_arm": 1.0,
            "axis": "days",
            "metric": "loss_all",
            "arrays": 5,
            "span": 100.0,
            "slope_per_10": slope,
            "se": 0.1,
        }

    slopes = pl.DataFrame(
        [
            row("oxygen=3%, agent=Contact_Inhibition", True, 0.0),
            row("agent=Contact_Inhibition", True, 0.1),
            row("oxygen=3%, agent=control", False, 1.0),
            row("agent=control", False, 2.0),
            row("agent=Contact_Inhibition_Regrowth", False, 9.0),
        ]
    )
    out = arrest_contrast(slopes)
    assert out["dividing_per_10_days"].to_list() == [1.0, 2.0]
    assert out["dividing_arm"].to_list() == ["oxygen=3%, agent=control", "agent=control"]


def test_presplit_reference_only_when_the_culture_starts_in_the_standard_arm():
    from entropy_age.experiments.solo_wcgw_culture import with_presplit

    culture = pl.DataFrame(
        {
            "arm": ["standard", "agent=MMC", "agent=MMC", "agent=control", "agent=control"],
            "days": [8.0, 18.0, 25.0, 18.0, 25.0],
            "doublings": [7.0, 8.0, 8.2, 9.0, 11.0],
        }
    )
    mmc = culture.filter(pl.col("arm") == "agent=MMC")
    assert with_presplit(mmc, culture)["days"].to_list() == [8.0, 18.0, 25.0]
    split = culture.filter(pl.col("arm") != "standard")  # no common start: each arm keeps its own
    assert with_presplit(split.filter(pl.col("arm") == "agent=MMC"), split).height == 2


def test_flank_context_sorts_the_three_contexts():
    from entropy_age.experiments.solo_wcgw_culture import flank_context

    pad = "A" * 40
    assert flank_context(pad + "[CG]" + pad) == "solo_wcgw"
    assert flank_context(pad[:-10] + "CG" + "A" * 8 + "[CG]" + pad) == "wcgw"  # a CpG neighbour within 35 bp
    assert flank_context("A" * 39 + "C[CG]G" + "A" * 39) == "scgs"
    assert flank_context("A" * 39 + "C[CG]A" + "A" * 39) == "other"
    assert flank_context("AA[CG]AA") is None


def test_timing_classes_and_fixed_strata():
    from entropy_age.experiments.solo_wcgw_culture import fixed_strata, timing_classes

    t = np.array([10.0, 50.0, 90.0, np.nan])
    assert timing_classes(t, np.array([30.0, 70.0])).tolist() == ["late", "mid", "early", "unknown"]
    assert fixed_strata(t, np.array([30.0, 70.0])).tolist() == [0, 1, 2, 3]


def _toy_unit(rng, n=600):
    hk = np.zeros(n, dtype=bool)
    hk[:60] = True
    table = pl.DataFrame(
        {
            "cpg": [f"cg{i}" for i in range(n)],
            "context": ["solo_wcgw"] * n,
            "island": rng.choice(["OpenSea", "Shore"], n),
            "hk": hk,
            "hk_gene": [f"g{i // 3}" if i < 60 else None for i in range(n)],
            "promoter": rng.random(n) < 0.5,
            "timing": rng.uniform(0, 100, n),
            "timing_class": rng.choice(["early", "late"], n),
            "timing_stratum": rng.integers(0, 5, n),
            "domain": rng.choice(["commonPMD", "commonHMD", "neither"], n),
        }
    )
    return {"loss": rng.normal(0, 1, (4, n)), "reference": rng.uniform(0.2, 0.9, n), "table": table}


def test_design_scores_equal_matched_weights_and_drop_undrawn_genes():
    from entropy_age.experiments.solo_wcgw_culture import design_parts, restriction, score_parts

    rng = np.random.default_rng(3)
    u = _toy_unit(rng)
    genes = {f"g{i}": i for i in range(20)}
    parts = design_parts(u, ("solo_wcgw", "early", True), genes)
    counts = np.ones(20)
    counts[0] = 0  # gene g0 not drawn
    counts[1] = 2  # g1 drawn twice
    hk_loss, matched, _ = score_parts(parts, "other", counts)
    t = u["table"]
    keep = restriction(t, "early")
    from entropy_age.stochasticity import baseline_strata

    level = np.full(t.height, -1)
    level[keep] = baseline_strata(u["reference"][keep], 10)
    strata = np.char.add(
        np.char.add(np.char.add(np.char.add(t["island"].to_numpy().astype(str), "|"), level.astype(str)), "|"),
        t["timing_stratum"].to_numpy().astype(str),
    )
    gene_index = np.array([genes[g] if g is not None else 0 for g in t["hk_gene"].to_list()])
    count_in = np.where(t["hk"].to_numpy() & keep, counts[gene_index], 0.0)
    s, m = matched_weights(count_in, keep & ~t["hk"].to_numpy(), strata)
    assert np.allclose(hk_loss, u["loss"] @ s) and np.allclose(matched, u["loss"] @ m)
    assert s[:3].sum() == 0  # g0's CpGs carry no weight


def test_slope_frame_matches_within_slope_and_leaves_lines_out():
    from entropy_age.experiments.solo_wcgw_culture import SlopeFrame

    rng = np.random.default_rng(5)
    rows = []
    for line in ("a", "b", "c", "d"):
        for unit in (1, 2):
            for d in range(5):
                rows.append({"donor": line, "unit": f"{line}{unit}", "doublings_gained": float(d * 3 + unit)})
    meta = pl.DataFrame(rows)
    y = rng.normal(0, 1, meta.height) + 0.4 * meta["doublings_gained"].to_numpy() / 10
    frame = SlopeFrame(meta)
    t, units = meta["doublings_gained"].to_numpy() / 10, meta["unit"].to_numpy()
    assert frame.slope(y) == pytest.approx(within_slope(y, t, units))
    keep = meta["donor"].to_numpy() != "b"
    assert frame.leave_one_out(y)[1] == pytest.approx(within_slope(y[keep], t[keep], units[keep]))
    assert frame.line_boot(y, np.array([[0, 1, 2, 3]]))[0] == pytest.approx(frame.slope(y))


def test_slope_frame_applies_the_identifiability_rule():
    from entropy_age.experiments.solo_wcgw_culture import SlopeFrame

    # nearly all doubling variance lies between units: within_slope gives no slope, nor may the fast path
    meta = pl.DataFrame(
        {
            "donor": ["a"] * 4 + ["b"] * 4,
            "unit": ["a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2"],
            "doublings_gained": [0.0, 0.01, 100.0, 100.01, 200.0, 200.01, 300.0, 300.01],
        }
    )
    y = np.arange(8.0)
    t, units = meta["doublings_gained"].to_numpy() / 10, meta["unit"].to_numpy()
    assert np.isnan(within_slope(y, t, units))
    frame = SlopeFrame(meta)
    assert np.isnan(frame.slope(y)) and np.isnan(frame.leave_one_out(y)).all()


def test_outside_hmd_leaves_unknown_domains_out():
    from entropy_age.experiments.solo_wcgw_culture import restriction

    t = pl.DataFrame({"domain": ["commonPMD", "commonHMD", "neither", "unknown"]})
    assert restriction(t, "notHMD").tolist() == [True, False, True, False]


def test_non_essential_design_repeats_score_unit_and_keeps_both_sets_out_of_the_pool():
    from entropy_age.experiments.solo_wcgw_culture import design_parts, score_parts

    rng = np.random.default_rng(11)
    u = _toy_unit(rng)
    n = u["table"].height
    ne_gene = [f"n{i // 4}" if 60 <= i < 140 else None for i in range(n)]
    ne_gene[0] = "n_both"  # a CpG linked to a housekeeping and a non-essential gene stays out of the set
    u["table"] = u["table"].with_columns(
        non_essential_gene=pl.Series(ne_gene, dtype=pl.String),
        non_essential=pl.Series([g is not None for g in ne_gene]),
        pmd=pl.lit(False),
    )
    genes = {g: i for i, g in enumerate(sorted({g for g in ne_gene if g is not None}))}
    parts = design_parts(u, ("solo_wcgw", None, False), genes, "non_essential")
    assert parts["hk_loss"].shape[1] == 80  # CpGs 60-139, not CpG 0
    set_loss, matched, _ = score_parts(parts, "other", np.ones(len(genes)))
    expected = score_unit(u["loss"], u["reference"], u["table"], 0.0)["gap_non_essential"]
    assert np.allclose(set_loss - matched, expected)
    pool = parts["pools"]["other"]["n"].sum()
    assert pool == n - 60 - 80  # neither housekeeping nor non-essential CpGs are controls


def test_non_essential_design_needs_twenty_cpgs_as_score_unit_does():
    from entropy_age.experiments.solo_wcgw_culture import design_parts

    rng = np.random.default_rng(12)
    u = _toy_unit(rng)
    n = u["table"].height
    ne_gene = [f"n{i}" if 60 <= i < 70 else None for i in range(n)]
    u["table"] = u["table"].with_columns(
        non_essential_gene=pl.Series(ne_gene, dtype=pl.String),
        non_essential=pl.Series([g is not None for g in ne_gene]),
        pmd=pl.lit(False),
    )
    genes = {g: i for i, g in enumerate(g for g in ne_gene if g is not None)}
    assert design_parts(u, ("solo_wcgw", None, False), genes, "non_essential") is None
    assert np.isnan(score_unit(u["loss"], u["reference"], u["table"], 0.0)["gap_non_essential"]).all()


def test_sturm_series_end_to_end_from_metadata_qc_and_betas(tmp_path):
    """Metadata -> QC -> rebuilt betas -> scores -> fits on a tiny Sturm fixture (Codex review 11): shuffled sample
    and probe order, a drug arm and a SURF1 line left out, a QC failure, an arrested arm on its own reference, a CpG
    missing in one array, and a chip offset that only the chip-held slope removes."""
    from entropy_age.data.cohorts import layout
    from entropy_age.experiments import solo_wcgw_culture as swc
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    rng = np.random.default_rng(11)
    n = 400
    cpgs = [f"cg{i:04d}" for i in range(n)]
    hk = np.arange(n) % 10 == 0  # 40 housekeeping CpGs, in every island class and starting decile
    table = pl.DataFrame(
        {
            "cpg": cpgs,
            "island": np.array(["OpenSea", "Shore", "Island", "Shelf"])[np.arange(n) % 4],
            "hk": hk,
            "hk_gene": [f"g{i // 20}" if h else None for i, h in enumerate(hk)],
            "promoter": hk | (np.arange(n) % 3 == 0),
            "non_essential": np.zeros(n, dtype=bool),
            "pmd": np.zeros(n, dtype=bool),
        }
    )
    base = rng.uniform(0.3, 0.9, n)
    rate = np.where(hk, 0.25, 0.5)  # points lost per doubling: housekeeping CpGs half as much

    rows, betas, qc = [], {}, []

    def array(sample, donor, sub, agent, doublings, days, chip="A", offset=0.0, passed=True):
        rows.append(
            {
                "sample": sample,
                "donor": donor,
                "subexperiment": sub,
                "oxygen": "Ambient",
                "agent": agent,
                "serum": "Full",
                "vector": "NA",
                "doublings": doublings,
                "days": days,
            }
        )
        betas[sample] = base - (rate * doublings + offset) / 100.0
        qc.append({"sample": sample, "chip": chip, "frac_dt": 0.99 if passed else 0.5, "platform": "EPIC"})

    for donor, chip_b in (("D1", 2.0), ("D2", 6.0)):
        for k, d in enumerate((0.0, 2.0, 4.0, 6.0, 8.0)):
            on_b = d == chip_b
            array(
                f"GSM{donor}{k}", donor, "Baseline profiling", "control", d, 10 * d, "B" if on_b else "A", float(on_b)
            )
    for k, day in enumerate((50.0, 60.0, 70.0, 80.0)):  # arrested: days pass, no doublings, no further loss
        array(f"GSMCI{k}", "D1", "Study part 3", "Contact_Inhibition", 8.0, day)
    array("GSMCI9", "D1", "Study part 3", "Contact_Inhibition", 8.0, 90.0, passed=False)  # fails QC
    array("GSMAZA", "D1", "Baseline profiling", "5-azacytidine", 9.0, 100.0)  # a drug arm: left out
    array("GSMSURF", "M1", "SURF1 mutation", "control", 4.0, 40.0)  # a genetic disease line: left out
    betas["GSMD13"][7] = np.nan  # one CpG missing in one array of D1: out of D1's complete matrix
    samples = pl.DataFrame(rows).sample(fraction=1.0, shuffle=True, seed=3)
    out = layout("sturm", paths).samples
    out.parent.mkdir(parents=True)
    samples.write_parquet(out)
    folder = tmp_path / "data" / "interim" / "sesame" / swc.STURM
    folder.mkdir(parents=True)
    pl.DataFrame(qc).sample(fraction=1.0, shuffle=True, seed=4).write_csv(folder / f"{swc.STURM}.qc.csv")
    order = rng.permutation(n)
    names = list(betas)[::-1]
    chunk = pl.DataFrame({"probe": [cpgs[i] for i in order]} | {s: betas[s][order] for s in names})
    chunk.write_parquet(folder / f"{swc.STURM}.betas.chunk000.parquet")

    s = swc.cultures(paths)
    assert set(s["sample"]) == {f"GSMD{d}{k}" for d in (1, 2) for k in range(5)} | {f"GSMCI{k}" for k in range(4)}
    assert s.filter(pl.col("arrested"))["sample"].sort().to_list() == [f"GSMCI{k}" for k in range(4)]
    chip = dict(s.select("sample", "chip").iter_rows())
    assert chip["GSMD11"] == "B" and chip["GSMD23"] == "B" and chip["GSMD10"] == "A"

    units = swc.prepare(paths, s, table)
    meta = {u["meta"]["unit"][0]: u for u in units}
    d1 = next(u for name, u in meta.items() if "Baseline profiling | D1" in name)
    ci = next(u for name, u in meta.items() if "Contact_Inhibition" in name)
    assert d1["meta"]["cpgs"][0] == n - 1 and "cg0007" not in d1["table"]["cpg"].to_list()
    assert d1["meta"]["reference_samples"][0] == "GSMD10" and ci["meta"]["reference_samples"][0] == "GSMCI0"
    # every loss is bound to its own CpG and array: reference minus beta, in points
    cols = d1["table"]["cpg"].to_list()
    j = cols.index("cg0010")
    later = d1["meta"]["sample"].to_list()
    assert later == ["GSMD11", "GSMD12", "GSMD13", "GSMD14"]
    assert d1["loss"][later.index("GSMD14"), j] == pytest.approx(0.25 * 8.0)
    assert d1["loss"][later.index("GSMD11"), j] == pytest.approx(0.25 * 2.0 + 1.0)  # the chip-B offset

    sc = swc.scores(units, 0.0)
    slopes = swc.unit_slopes(sc)
    ci_rows = slopes.filter(pl.col("arrested"))
    assert set(ci_rows["axis"]) == {"days"}  # no doublings in the arrested arm: no slope on them
    assert ci_rows.filter(pl.col("metric") == "loss_all")["slope_per_10"][0] == pytest.approx(0.0, abs=1e-9)

    fits = swc.series_fits(sc, slopes, 0, np.random.default_rng(0))
    sturm = {r["metric"]: r for r in fits.filter(pl.col("selection") == "Sturm, 5 lines").iter_rows(named=True)}
    mean_rate = float(np.mean(rate[np.arange(n) != 7]))  # roughly, D1 lacks one CpG
    assert sturm["loss_all"]["slope_chip_held"] == pytest.approx(10 * mean_rate, abs=0.01)
    assert abs(sturm["loss_all"]["slope_per_10_pd"] - sturm["loss_all"]["slope_chip_held"]) > 0.05  # chip bias
    for pool in ("other", "promoter"):  # the uniform chip offset cancels in the gap
        assert sturm[f"gap_{pool}"]["slope_per_10_pd"] == pytest.approx(10 * (0.25 - 0.5), abs=1e-6)
        assert sturm[f"loss_hk_{pool}"]["slope_chip_held"] == pytest.approx(2.5, abs=1e-6)
    assert sturm["loss_all"]["arrays"] == 8 and sturm["loss_all"]["units"] == 2
