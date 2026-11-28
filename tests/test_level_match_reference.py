"""P6: matching variables from another cohort's young donors (``experiments.level_match_reference``)."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.level_match_reference import SAME, impute, reference_rows, restrict, summarise, tables
from entropy_age.experiments.tissue_contrast import Prepared


def _prepared(n: int, p: int, seed: int, label: str, features: list[str]) -> Prepared:
    rng = np.random.default_rng(seed)
    ages = np.sort(rng.uniform(20, 80, n))
    x = np.clip(rng.uniform(0.05, 0.95, p) + rng.normal(0, 0.02, (n, p)), 0, 1).astype(np.float32)
    return Prepared(label, x, ages, None, features, np.arange(n).astype(str), None)


def _probes(features):
    return pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(len(features))], "probe_type": "I"})


def test_reference_rows_take_the_youngest_fifth_or_all():
    p = _prepared(50, 5, 0, "a", [f"c{i}" for i in range(5)])
    rows = np.arange(50)
    young = reference_rows(p, rows, True)
    assert young.size == 10 and p.ages[young].max() <= np.quantile(p.ages, 0.2)
    assert reference_rows(p, rows, False).tolist() == rows.tolist()


def test_tables_take_level_and_sd_from_the_reference_and_drift_from_the_drift_rows():
    from entropy_age.experiments.tissue_contrast import _lifespan_metrics

    features = [f"c{i}" for i in range(40)]
    drift = _prepared(200, 40, 1, "drift", features)
    rng = np.random.default_rng(9)
    drift = Prepared(
        "drift",
        drift.x + (drift.ages[:, None] * rng.uniform(0, 2e-3, 40)).astype(np.float32),
        drift.ages,
        None,
        features,
        drift.donors,
        None,
    )
    ref = _prepared(30, 40, 2, "ref", features)
    rows = np.arange(5, 200, 2)  # a non-trivial subset of the drift cohort
    ref_rows = np.arange(10, 30)
    out = tables(
        drift,
        rows,
        {SAME: (drift, drift.x, reference_rows(drift, rows, True)), "other": (ref, ref.x, ref_rows)},
        _probes(features),
    )
    other, same = out["other"], out[SAME]
    assert np.allclose(other["baseline"].to_numpy(), ref.x[ref_rows].astype(float).mean(axis=0))
    assert np.allclose(other["young_sd"].to_numpy(), ref.x[ref_rows].astype(float).std(axis=0, ddof=1))
    young = reference_rows(drift, rows, True)
    assert np.allclose(same["young_sd"].to_numpy(), drift.x[young].astype(float).std(axis=0, ddof=1))
    expected = _lifespan_metrics(drift, rows)
    for m in ("entropy_slope", "dispersion_slope", "abs_slope", "kl_slope"):
        assert other["feature"].to_list() == expected["feature"].to_list()
        assert np.allclose(other[m].to_numpy(), expected[m].to_numpy())
        assert np.abs(expected[m].to_numpy()).max() > 0
    whole = _lifespan_metrics(drift)  # the rows matter: the subset's slopes are not the whole cohort's
    assert not np.allclose(whole["abs_slope"].to_numpy(), expected["abs_slope"].to_numpy())
    with pytest.raises(ValueError, match="CpG order"):
        tables(drift, rows, {"x": (restrict(ref, features[::-1]), ref.x, ref_rows)}, _probes(features))


def test_reference_stats_skip_gaps_and_need_enough_measured_rows():
    from entropy_age.experiments.level_match_reference import MIN_OBSERVED, reference_stats

    raw = np.full((30, 3), 0.5)
    raw[:, 0] = np.linspace(0.1, 0.4, 30)
    raw[::2, 1] = np.nan
    raw[1::2, 1] = 0.9  # measured values only: the gaps are not filled from other donors
    raw[MIN_OBSERVED - 1 :, 2] = np.nan
    mean, sd = reference_stats(raw, np.arange(30))
    assert mean[0] == pytest.approx(raw[:, 0].mean()) and sd[0] == pytest.approx(raw[:, 0].std(ddof=1))
    assert mean[1] == pytest.approx(0.9) and sd[1] == pytest.approx(0.0)
    assert np.isnan(mean[2]) and np.isnan(sd[2])


def test_impute_drops_gappy_cpgs_and_fills_the_rest_with_medians():
    x = np.array([[0.1, np.nan, 0.5], [0.3, np.nan, np.nan], [0.2, np.nan, 0.6]] + [[0.2, 0.4, 0.6]] * 37)
    out, raw, kept = impute(x, ["a", "b", "c"])
    assert kept == ["a", "c"] and out[1, 1] == pytest.approx(0.6) and np.isnan(raw[1, 1])


def test_summarise_pairs_each_reference_with_the_same_donor_gap():
    base = {"tissue": "t", "pool": "rest", "scope": "all CpGs", "matching": "m", "group": "g", "metric": "kl_slope"}
    point = pl.DataFrame(
        [base | {"reference": SAME, "matched_diff": -2.0}, base | {"reference": "r", "matched_diff": -1.0}]
    )
    boot = pl.DataFrame(
        [
            base | {"reference": ref, "matched_diff": v + b, "boot": b}
            for b in range(10)
            for ref, v in ((SAME, -2.0), ("r", -1.0))
        ]
    )
    spread = pl.DataFrame({"tissue": ["t"], "metric": ["kl_slope"], "sd": [2.0]})
    s = summarise(point, boot, spread).sort("reference")
    r = s.filter(pl.col("reference") == "r").row(0, named=True)
    assert r["gap_sd"] == pytest.approx(-0.5) and r["delta_sd"] == pytest.approx(0.5)
    assert r["delta_lo"] == pytest.approx(0.5) and r["delta_hi"] == pytest.approx(0.5)  # paired in each resample


def test_deltas_only_on_the_whole_scope_and_gated_on_their_own_resamples():
    base = {"tissue": "t", "pool": "rest", "matching": "m", "group": "g", "metric": "kl_slope"}
    point, boot = [], []
    for scope in ("all CpGs", "lowest decile"):
        point += [
            base | {"scope": scope, "reference": SAME, "matched_diff": -2.0},
            base | {"scope": scope, "reference": "r", "matched_diff": -1.0},
        ]
        for b in range(10):
            same = -2.0 if b == 0 else float("nan")  # the same-donor gap is finite in one resample only
            boot += [
                base | {"scope": scope, "reference": SAME, "matched_diff": same, "boot": b},
                base | {"scope": scope, "reference": "r", "matched_diff": -1.0, "boot": b},
            ]
    spread = pl.DataFrame({"tissue": ["t"], "metric": ["kl_slope"], "sd": [1.0]})
    s = summarise(pl.DataFrame(point), pl.DataFrame(boot), spread)
    r = s.filter((pl.col("reference") == "r") & (pl.col("scope") == "all CpGs")).row(0, named=True)
    assert r["delta_sd"] == pytest.approx(1.0) and r["delta_lo"] is None and r["gap_lo"] == pytest.approx(-1.0)
    low = s.filter((pl.col("reference") == "r") & (pl.col("scope") == "lowest decile")).row(0, named=True)
    assert low["delta_sd"] is None


def test_one_processing_per_run_and_rebuilt_runs_need_every_rebuild(monkeypatch, tmp_path):
    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    monkeypatch.setenv("ENTROPY_AGE_BETAS", "methylprep")
    with pytest.raises(RuntimeError, match="deposited or sesame"):
        lmr.processing()
    monkeypatch.setenv("ENTROPY_AGE_BETAS", "sesame")
    monkeypatch.setattr(lmr.reprocessed, "available", lambda paths, gse: gse != lmr.ERISK)
    with pytest.raises(RuntimeError, match="GSE105018"):
        lmr.processing(ProjectPaths(tmp_path))
    monkeypatch.setattr(lmr.reprocessed, "available", lambda paths, gse: True)
    monkeypatch.setattr(lmr, "recipe", lambda paths, gse: "QCDPB")
    assert lmr.processing(ProjectPaths(tmp_path)) == "sesame"
    monkeypatch.delenv("ENTROPY_AGE_BETAS")
    assert lmr.processing(ProjectPaths(tmp_path)) == "deposited"


def test_rebuilt_runs_take_one_recipe_on_hm450(monkeypatch, tmp_path):
    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    folder = tmp_path / "data" / "interim" / "sesame" / "GSE1"
    folder.mkdir(parents=True)
    qc = {"sample": ["GSM1", "GSM2"], "chip": ["1", "1"], "frac_dt": [0.99, 0.99], "platform": ["HM450", "HM450"]}
    pl.DataFrame(qc).write_csv(folder / "GSE1.qc.csv")
    assert lmr.recipe(paths, "GSE1").startswith("QCDPB")  # no prep column: the first script's fixed recipe
    pl.DataFrame(qc | {"prep": ["QCDPB", "CDPB"]}).write_csv(folder / "GSE1.qc.csv")
    with pytest.raises(RuntimeError, match="CDPB"):
        lmr.recipe(paths, "GSE1")
    pl.DataFrame(qc | {"prep": ["QCDPB", "QCDPB"]}).write_csv(folder / "GSE1.qc.csv")
    assert lmr.recipe(paths, "GSE1") == "QCDPB"
    pl.DataFrame(qc | {"platform": ["HM450", "EPIC"]}).write_csv(folder / "GSE1.qc.csv")
    with pytest.raises(RuntimeError, match="EPIC"):
        lmr.recipe(paths, "GSE1")
    monkeypatch.setenv("ENTROPY_AGE_SESAME_DIR", "sesame_clock")
    with pytest.raises(RuntimeError, match="default rebuild"):
        lmr.recipe(paths, "GSE1")


def test_paired_mask_keeps_cells_both_versions_measure():
    from entropy_age.experiments.level_match_reference import paired_mask

    a = np.array([[0.1, np.nan], [0.2, 0.3]], dtype=np.float32)
    b = np.array([[0.1, 0.4], [np.inf, 0.3]], dtype=np.float32)
    out = paired_mask({("c", "deposited"): a, ("c", "sesame"): b, ("d", "sesame"): a})
    assert (
        np.isnan(out[("c", "deposited")][[0, 1], [1, 0]]).all() and np.isnan(out[("c", "sesame")][[0, 1], [1, 0]]).all()
    )
    assert out[("c", "sesame")][1, 1] == pytest.approx(0.3) and np.isnan(out[("d", "sesame")][0, 1])


def test_rebuilt_plan_has_no_hannum():
    from entropy_age.experiments.level_match_reference import PLAN, SESAME_PLAN, active_plan

    assert active_plan("deposited") is PLAN and active_plan("sesame") is SESAME_PLAN
    assert list(SESAME_PLAN) == ["epic_italy"]
    refs = {ref for ref, _ in SESAME_PLAN["epic_italy"].values()}
    assert refs == {"epic_italy", "charge", "young_adults", "erisk"}


def test_one_per_family_keeps_the_drawn_twin_when_usable_else_the_co_twin(tmp_path):
    from entropy_age.experiments.level_match_reference import erisk_people, one_per_family

    matrix = tmp_path / "GSE105018_series_matrix.txt.gz"
    _erisk_matrix(matrix, ["11", "12", "21", "22", "31"], ["1", "1", "2", "2", "3"], ["1", "2", "1", "2", "1"])
    people = erisk_people(matrix)
    for seed in range(20):
        whole = one_per_family(people, np.random.default_rng(seed))
        drawn = dict(whole.select("family", "sample").iter_rows())
        dropped = drawn["1"]  # family 1's drawn twin fails; family 3's only array fails
        usable = [s for s in people["sample"].to_list() if s not in (dropped, drawn["3"])]
        cut = one_per_family(people, np.random.default_rng(seed), usable)
        picked = dict(cut.select("family", "sample").iter_rows())
        assert set(picked) == {"1", "2"} and picked["2"] == drawn["2"]
        assert picked["1"] != dropped and picked["1"] in usable


def test_young_sd_contrast_pairs_versions_and_separates_cohorts():
    from entropy_age.experiments.level_match_reference import processing_table

    rng = np.random.default_rng(3)
    p = 400
    features = [f"c{i}" for i in range(p)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(p)]})
    house = features[::4]
    membership = pl.DataFrame({"feature": house, "group": "housekeeping"})
    membership = pl.concat([membership, pl.DataFrame({"feature": features[1::4], "group": "core_essential"})])
    membership = pl.concat([membership, pl.DataFrame({"feature": features[2::8], "group": "non_essential"})])
    level = rng.uniform(0.05, 0.9, p)
    spread = np.where(np.isin(features, house), 0.01, 0.03)  # housekeeping CpGs vary less between donors

    def cohort(n, scale):
        return (level + rng.normal(0, 1, (n, p)) * spread * scale).astype(np.float32)

    a, b = cohort(60, 1.0), cohort(80, 0.5)
    b_dep = b.copy()
    b_dep[:, 5] = np.nan  # a CpG the deposited E-Risk lacks drops out of every version
    versions = {
        ("young_adults", "deposited"): a,
        ("young_adults", "sesame"): a,
        ("erisk", "deposited"): b_dep,
        ("erisk", "sesame"): b,
    }
    t = processing_table(versions, features, probes, membership, 30, np.random.default_rng(0))
    assert set(t["cpgs"].to_list()) == {p - 1}
    assert set(t["strata"].to_list()) == {"own", "deposited"}
    hk = t.filter((pl.col("group") == "housekeeping") & (pl.col("strata") == "own"))
    row = {(r["cohort"], r["processing"]): r for r in hk.iter_rows(named=True)}
    ya, er = row[("GSE73103 young adults", "sesame")], row[("E-Risk one twin per family", "sesame")]
    assert ya["sd_diff_pp"] < er["sd_diff_pp"] < 0  # half the spread, half the difference
    assert ya["donors"] == 60 and er["donors"] == 80
    same = row[("GSE73103 young adults", "rebuilt minus deposited")]
    assert same["sd_diff_pp"] == pytest.approx(0.0, abs=1e-9) and same["lo_pp"] == pytest.approx(0.0, abs=1e-9)
    gap = row[("GSE73103 minus E-Risk", "sesame")]
    assert gap["sd_diff_pp"] == pytest.approx(ya["sd_diff_pp"] - er["sd_diff_pp"])
    assert gap["lo_pp"] <= gap["sd_diff_pp"] <= gap["hi_pp"] < 0
    fixed = t.filter((pl.col("group") == "housekeeping") & (pl.col("strata") == "deposited"))
    row = {(r["cohort"], r["processing"]): r for r in fixed.iter_rows(named=True)}
    assert row[("GSE73103 young adults", "deposited")]["sd_diff_pp"] == pytest.approx(ya["sd_diff_pp"])


def test_every_reference_is_read_on_the_same_cpgs():
    features = [f"c{i}" for i in range(30)]
    drift = _prepared(100, 30, 3, "drift", features)
    ref = _prepared(40, 30, 4, "ref", features)
    raw = ref.x.copy()
    raw[:35, 7] = np.nan  # c7 measured in 5 reference rows only: no matching variables there
    rows = np.arange(100)
    out = tables(
        drift,
        rows,
        {SAME: (drift, drift.x, reference_rows(drift, rows, True)), "other": (ref, raw, np.arange(40))},
        _probes(features),
    )
    assert "c7" not in out[SAME]["feature"].to_list()
    assert out[SAME]["feature"].to_list() == out["other"]["feature"].to_list()


def test_young_adult_reference_keeps_ages_18_to_30(tmp_path):
    import gzip

    from entropy_age.experiments.level_match_reference import young_adult_people

    rows = [
        "!Sample_geo_accession\t" + "\t".join(f'"GSM{i}"' for i in range(5)),
        "!Sample_characteristics_ch1\t"
        + "\t".join(f'"Sex: {s}"' for s in ("Male", "Female", "Male", "Female", "Male")),
        "!Sample_characteristics_ch1\t" + "\t".join('"tissue: blood"' for _ in range(5)),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"age: {a}"' for a in (15, 18, 25, 30, 31)),
        "!Sample_characteristics_ch1\t" + "\t".join('"weight category: Normal-weight"' for _ in range(5)),
    ]
    matrix = tmp_path / "GSE73103_series_matrix.txt.gz"
    with gzip.open(matrix, "wt") as handle:
        handle.write("\n".join(rows) + "\n")
    people = young_adult_people(matrix)
    assert people["sample"].to_list() == ["GSM1", "GSM2", "GSM3"]
    assert people["age"].to_list() == [18.0, 25.0, 30.0] and people["sex"].to_list() == ["Female", "Male", "Female"]


def test_young_adult_betas_bind_by_probe_and_sample_name(tmp_path, monkeypatch):
    import gzip

    from entropy_age.experiments import level_match_reference as lmr

    ages = [15, 30, 18, 25]
    header = [
        "!Sample_geo_accession\t" + "\t".join(f'"GSM{i}"' for i in range(4)),
        "!Sample_characteristics_ch1\t" + "\t".join('"Sex: Male"' for _ in range(4)),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"age: {a}"' for a in ages),
        "!Sample_characteristics_ch1\t" + "\t".join('"weight category: Normal-weight"' for _ in range(4)),
        "!series_matrix_table_begin",
        '"ID_REF"\t"GSM3"\t"GSM1"\t"GSM0"\t"GSM2"',  # samples out of accession order
        '"cg2"\t0.23\t0.21\t0.20\t0.22',  # probes out of the requested order; cg9 absent
        '"cg1"\t0.13\t0.11\t0.10\t0.12',
        "!series_matrix_table_end",
    ]
    geo = tmp_path / "input" / "geo" / "GSE73103"
    geo.mkdir(parents=True)
    with gzip.open(geo / "GSE73103_series_matrix.txt.gz", "wt") as handle:
        handle.write("\n".join(header) + "\n")
    monkeypatch.setattr(lmr, "matrix_path", lambda paths, gse: geo / f"{gse}_series_matrix.txt.gz")
    monkeypatch.setattr(lmr, "usable", lambda matrix, samples: samples)
    p, raw = lmr._young_adults(None, ["cg1", "cg9", "cg2"])
    assert p.donors.tolist() == ["GSM2", "GSM3", "GSM1"] and p.ages.tolist() == [18.0, 25.0, 30.0]
    assert p.features == ["cg1", "cg2"]  # the absent CpG is dropped
    assert np.allclose(p.x, [[0.12, 0.22], [0.13, 0.23], [0.11, 0.21]]) and np.allclose(raw, p.x)


def test_cohort_summary_reports_the_rows_each_reference_uses():
    from entropy_age.experiments.level_match_reference import PLAN, cohort_summary

    features = ["c0", "c1"]
    cohorts = {k: (p := _prepared(50, 2, i, k, features), p.x) for i, k in enumerate(("hannum", "epic_italy"))}
    cohorts |= {k: (p := _prepared(30, 2, 9, k, features), p.x) for k in ("charge", "young_adults", "erisk")}
    own = {k: {"c0", "c1"} | ({"c9"} if k != "erisk" else set()) for k in cohorts}
    rows = {r["key"]: r for r in cohort_summary(cohorts, 2, own)}
    assert rows["erisk"]["own_cpgs"] == 2 and rows["erisk"]["common_cpgs_without"] == 3
    assert rows["hannum"]["common_cpgs_without"] == 2 and "dasen" in rows["erisk"]["processing"]
    assert rows["erisk"]["reference_rule"] == "all donors" and rows["erisk"]["reference_arrays"] == 30
    assert set(rows) == {ref for refs in PLAN.values() for ref, _ in refs.values()}
    assert rows["hannum"]["reference_rule"] == "youngest fifth" and rows["hannum"]["reference_arrays"] == 10
    assert rows["young_adults"]["reference_rule"] == "all donors" and rows["young_adults"]["reference_arrays"] == 30
    assert "BMIQ" in rows["young_adults"]["processing"]


def _erisk_matrix(path, titles, families, twins):
    import gzip

    n = len(titles)
    rows = [
        "!Sample_title\t" + "\t".join(f'"whole blood taken from sample {t}"' for t in titles),
        "!Sample_geo_accession\t" + "\t".join(f'"GSM{i}"' for i in range(n)),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"familyid: {f}"' for f in families),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"twinid: {t}"' for t in twins),
        "!Sample_characteristics_ch1\t" + "\t".join('"gender: F"' for _ in range(n)),
        "!Sample_characteristics_ch1\t" + "\t".join('"zygosity: MZ"' for _ in range(n)),
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(rows) + "\n")


def test_erisk_keeps_one_twin_per_family_drawn_by_the_seed(tmp_path):
    from entropy_age.experiments.level_match_reference import erisk_people, one_per_family

    matrix = tmp_path / "GSE105018_series_matrix.txt.gz"
    _erisk_matrix(matrix, ["11", "12", "21", "22", "31"], ["1", "1", "2", "2", "3"], ["1", "2", "1", "2", "1"])
    people = erisk_people(matrix)
    assert people["column"].to_list() == ["11", "12", "21", "22", "31"]
    chosen = [one_per_family(people, np.random.default_rng(s)) for s in range(40)]
    assert all(c["family"].sort().to_list() == ["1", "2", "3"] for c in chosen)
    picks = {tuple(c["sample"].to_list()) for c in chosen}
    assert len(picks) > 1  # random among twins
    assert one_per_family(people, np.random.default_rng(5)).equals(
        one_per_family(people.reverse(), np.random.default_rng(5))
    )


def test_erisk_betas_stream_into_parts_and_bind_by_cpg_and_sample_number(tmp_path, monkeypatch):
    import gzip

    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    geo = paths.input / "geo" / "GSE105018"
    geo.mkdir(parents=True)
    _erisk_matrix(geo / "GSE105018_series_matrix.txt.gz", ["11", "12", "21"], ["1", "1", "2"], ["1", "2", "1"])
    lines = ['"","12","21","11"'] + [f'"cg{i}",{i / 100 + 0.002},{i / 100 + 0.003},{i / 100 + 0.001}' for i in range(7)]
    lines[3] = '"cg2",NA,0.023,0.021'
    lines[4] = '"cg3",0.032,0.033,NA'  # a gap in a selected column stays NaN
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    monkeypatch.setattr(lmr, "ERISK_ROWS_PER_PART", 3)  # several parts
    parts = lmr.erisk_parts(paths)
    assert len(list(parts.glob("*.parquet"))) == 3
    assert not list(parts.parent.glob("erisk_betas.partial*"))  # the temporary folder became the cache
    assert lmr.erisk_parts(paths) == parts  # a finished cache is reused
    x = lmr.read_erisk_betas(parts, ["cg5", "cg99", "cg2", "cg3"], ["11", "21"])
    assert np.allclose(x[:, 0], [0.051, 0.053]) and np.isnan(x[:, 1]).all()
    assert np.allclose(x[:, 2], [0.021, 0.023]) and np.isnan(x[0, 3]) and x[1, 3] == pytest.approx(0.033)
    monkeypatch.setattr(lmr, "matrix_path", lambda paths, gse: geo / f"{gse}_series_matrix.txt.gz")
    rng = np.random.default_rng(0)
    p, raw = lmr._erisk(paths, ["cg1", "cg2", "cg3"], rng)
    assert p.x.shape[0] == 2 and set(p.ages) == {18.0}
    first = "11" if "GSM0" in p.donors else "12"
    assert raw[list(p.donors).index("GSM2"), 0] == pytest.approx(0.013)
    assert raw[list(p.donors).index("GSM0" if first == "11" else "GSM1"), 0] == pytest.approx(
        0.011 if first == "11" else 0.012
    )


def test_erisk_cache_refuses_a_header_only_csv(tmp_path):
    import gzip

    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    geo = paths.input / "geo" / "GSE105018"
    geo.mkdir(parents=True)
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write('"","11"\n')
    with pytest.raises(ValueError, match="no CpG rows"):
        lmr.erisk_parts(paths)
    assert not list((paths.interim / "cohorts").iterdir())


def test_equal_size_references_draw_from_the_reference_rows(monkeypatch):
    from entropy_age.experiments import level_match_reference as lmr

    features = [f"c{i}" for i in range(60)]
    cohorts = {k: _prepared(n, 60, s, k, features) for k, n, s in (("d", 150, 1), ("r", 40, 2))}
    cohorts = {k: (p, p.x) for k, p in cohorts.items()}
    seen = {}

    def fake_tables(drift, rows, references, probes, drift_raw=None):
        seen.update({name: rows for name, (_, _, rows) in references.items()})
        return {}

    monkeypatch.setattr(lmr, "PLAN", {"d": {SAME: ("d", True), "other": ("r", False)}})
    monkeypatch.setattr(lmr, "tables", fake_tables)
    monkeypatch.setattr(lmr.pl, "concat", lambda frames: pl.DataFrame())
    full = {k: np.arange(p.x.shape[0]) for k, (p, _) in cohorts.items()}
    lmr.sample_rows("d", cohorts, None, None, full, 60, np.random.default_rng(1))
    assert seen["other"].size == 60 and np.unique(seen["other"]).size < 40  # with replacement, beyond the 40 rows
    young = reference_rows(cohorts["d"][0], full["d"], True)
    assert seen[SAME].size == 60 and set(seen[SAME]) <= set(young)


def test_subsample_summary_reports_spread_over_draws():
    from entropy_age.experiments.level_match_reference import summarise_subsamples

    base = {"tissue": "t", "reference": "r", "pool": "rest", "scope": "all CpGs", "matching": "m", "group": "g"}
    draws = pl.DataFrame(
        [base | {"metric": "kl_slope", "matched_diff": float(v), "sd_diff": v / 1000, "draw": v} for v in range(11)]
    )
    spread = pl.DataFrame({"tissue": ["t"], "metric": ["kl_slope"], "sd": [2.0]})
    r = summarise_subsamples(draws, spread).row(0, named=True)
    assert r["n_draws"] == 11 and r["gap_mean"] == pytest.approx(2.5)
    assert r["sd_diff_mean"] == pytest.approx(0.005) and r["sd_diff_sd_over_draws"] > 0
    assert r["gap_lo"] == pytest.approx(0.0) and r["gap_hi"] == pytest.approx(5.0)  # the shared nearest-rank quantile
    # undefined draws: statistics over the finite ones, counted, and null when fewer than half are finite
    nan = draws.with_columns(
        matched_diff=pl.when(pl.col("draw") < 3).then(float("nan")).otherwise(pl.col("matched_diff"))
    )
    r = summarise_subsamples(nan, spread).row(0, named=True)
    assert r["gap_valid"] == 8 and r["gap_mean"] == pytest.approx(3.25) and r["sd_diff_valid"] == 11
    mostly = draws.with_columns(
        matched_diff=pl.when(pl.col("draw") < 7).then(float("nan")).otherwise(pl.col("matched_diff"))
    )
    r = summarise_subsamples(mostly, spread).row(0, named=True)
    assert r["gap_mean"] is None and r["gap_lo"] is None and r["sd_diff_mean"] is not None


def test_erisk_reader_end_to_end_binds_sample_numbers_to_accessions(tmp_path, monkeypatch):
    """The deposited E-Risk path from a tiny fixture: series-matrix titles map the CSV's sample-number columns (in
    another order than the accessions) to GSMs, CpG rows come shuffled and are read back in the requested order, a
    failed array is replaced by its co-twin, and the rebuilt-run's deposited twin reads the same cells."""
    import gzip

    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    monkeypatch.delenv("ENTROPY_AGE_BETAS", raising=False)
    paths = ProjectPaths(tmp_path)
    geo = paths.input / "geo" / "GSE105018"
    geo.mkdir(parents=True)
    # GSM0..GSM5 carry sample numbers 305, 101, 102, 201, 306, 202 (families 3, 1, 1, 2, 3, 2)
    numbers = ["305", "101", "102", "201", "306", "202"]
    _erisk_matrix(geo / "GSE105018_series_matrix.txt.gz", numbers, ["3", "1", "1", "2", "3", "2"], list("112212"))
    value = {num: {f"cg{c}": int(num) / 1000 + c / 100_000 for c in range(8)} for num in numbers}
    value["201"]["cg4"] = None  # a gap in a read column stays a gap
    columns = ["202", "101", "306", "201", "305", "102"]  # the CSV's own column order
    cpg_rows = [5, 2, 7, 0, 3, 6, 1, 4]  # and its own CpG order
    lines = ['"",' + ",".join(f'"{c}"' for c in columns)]
    for c in cpg_rows:
        cells = ["NA" if value[n][f"cg{c}"] is None else repr(value[n][f"cg{c}"]) for n in columns]
        lines.append(f'"cg{c}",' + ",".join(cells))
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    monkeypatch.setattr(lmr, "ERISK_ROWS_PER_PART", 3)
    monkeypatch.setattr(lmr, "matrix_path", lambda paths, gse: geo / f"{gse}_series_matrix.txt.gz")
    failed = ("GSM1", "GSM5")  # 101 and 202 failed: families 1 and 2 keep 102 (GSM2) and 201 (GSM3) whatever the draw
    monkeypatch.setattr(lmr, "usable", lambda matrix, samples: [s for s in samples if s not in failed])
    features = ["cg6", "cg4", "cg0", "cg9", "cg3"]  # not sorted; cg9 absent from the CSV
    p, raw = lmr._erisk(paths, features, np.random.default_rng(1))
    gsm_to_number = dict(zip([f"GSM{i}" for i in range(6)], numbers, strict=True))
    assert len(p.donors) == 3 and {"GSM2", "GSM3"} <= set(p.donors) and not set(failed) & set(p.donors)
    assert sorted(gsm_to_number[d][0] for d in p.donors) == ["1", "2", "3"]  # one twin per family
    # the absent CpG and cg4 (missing in a third of the arrays, above MAX_MISSING) dropped, the requested order kept
    assert p.features == ["cg6", "cg0", "cg3"]
    for i, donor in enumerate(p.donors):
        num = gsm_to_number[donor]
        for j, f in enumerate(p.features):
            expected = value[num][f]
            if expected is None:
                assert np.isnan(raw[i, j])
            else:
                assert raw[i, j] == pytest.approx(expected, rel=1e-6)
    twin = lmr.deposited_twin(paths, "erisk", dataclasses_replace(p, features=["cg3", "cg9", "cg4", "cg6"]))
    for i, donor in enumerate(p.donors):
        num = gsm_to_number[donor]
        assert twin[i, 0] == pytest.approx(value[num]["cg3"], rel=1e-6) and np.isnan(twin[i, 1])
        assert twin[i, 3] == pytest.approx(value[num]["cg6"], rel=1e-6)
        gap = value[num]["cg4"] is None  # the gap stays a gap, bound to its own array
        assert np.isnan(twin[i, 2]) if gap else twin[i, 2] == pytest.approx(value[num]["cg4"], rel=1e-6)
    assert np.isnan(twin[:, 2]).sum() == 1


def dataclasses_replace(p, **kwargs):
    import dataclasses

    return dataclasses.replace(p, **kwargs)


@pytest.mark.parametrize("content", ["", '"","11"\n\n', '"","11"\n,0.5\n'], ids=["empty", "blank", "no-id"])
def test_erisk_cache_refuses_a_csv_without_cpgs_and_leaves_no_partial_folder(tmp_path, content):
    import gzip

    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    geo = paths.input / "geo" / "GSE105018"
    geo.mkdir(parents=True)
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write(content)  # no header at all, a header and a blank line, or a row with no CpG id
    with pytest.raises(ValueError, match="no CpG rows"):
        lmr.erisk_parts(paths)
    assert not list((paths.interim / "cohorts").iterdir())


def test_erisk_cache_cleans_up_a_failed_build_and_rebuilds_an_empty_cache(tmp_path):
    import gzip

    from entropy_age.experiments import level_match_reference as lmr
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    geo = paths.input / "geo" / "GSE105018"
    geo.mkdir(parents=True)
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write('"","11"\n"cg1",0.1\n"cg2",not-a-number\n')
    with pytest.raises(pl.exceptions.ComputeError):  # a value that is not a number
        lmr.erisk_parts(paths)
    assert not list((paths.interim / "cohorts").iterdir())  # the writer's folder removed
    (paths.interim / "cohorts" / "erisk_betas").mkdir()  # an empty cache left by a killed writer is not reused
    with gzip.open(geo / lmr.ERISK_CSV, "wt") as handle:
        handle.write('"","11"\n"cg1",0.1\n"cg2",0.2\n')
    parts = lmr.erisk_parts(paths)
    assert lmr.read_erisk_betas(parts, ["cg2", "cg1"], ["11"]).ravel().tolist() == pytest.approx([0.2, 0.1])


def test_erisk_reader_refuses_a_cpg_in_two_parts(tmp_path):
    from entropy_age.experiments.level_match_reference import read_erisk_betas

    for k, v in enumerate((0.1, 0.2)):
        pl.DataFrame({"cpg": ["cg1"], "11": [v]}).write_parquet(tmp_path / f"part-{k:04d}.parquet")
    with pytest.raises(ValueError, match="twice"):
        read_erisk_betas(tmp_path, ["cg1"], ["11"])


def _resampling_state():
    features = [f"c{i}" for i in range(240)]
    rng = np.random.default_rng(21)
    drift = _prepared(90, 240, 5, "drift", features)
    drift = Prepared(
        "drift",
        drift.x + (drift.ages[:, None] * rng.uniform(0, 2e-3, 240)).astype(np.float32),
        drift.ages,
        None,
        features,
        drift.donors,
        None,
    )
    ref = _prepared(40, 240, 6, "ref", features)
    raw = ref.x.copy()
    raw[::7, 3] = np.nan  # gaps in the reference's measured matrix
    membership = pl.concat(
        [
            pl.DataFrame({"feature": features[0::6], "group": "housekeeping"}),
            pl.DataFrame({"feature": features[1::6], "group": "core_essential"}),
            pl.DataFrame({"feature": features[2::6], "group": "non_essential"}),
        ]
    )
    plan = {"d": {SAME: ("d", True), "other": ("r", False)}}
    return {
        "cohorts": {"d": (drift, drift.x), "r": (ref, raw)},
        "probes": _probes(features),
        "membership": membership,
        "plan": plan,
    }


@pytest.mark.parametrize("task", ["_boot", "_subsample"])
def test_serial_and_spawned_resampling_agree(task):
    from polars.testing import assert_frame_equal

    from entropy_age.experiments import level_match_reference as lmr

    state = _resampling_state()
    jobs = [("d", s) for s in (3, 17, 2**30 + 5)]
    fn = getattr(lmr, task)
    serial = pl.DataFrame(lmr._parallel(jobs, 1, state, fn))
    spawned = pl.DataFrame(lmr._parallel(jobs, 2, state, fn))
    assert serial.height > 0
    key = "boot" if task == "_boot" else "draw"
    assert sorted(set(serial[key])) == [3, 17, 2**30 + 5]
    assert_frame_equal(serial, spawned)
    again = pl.DataFrame(lmr._parallel(jobs[::-1], 1, state, fn))  # a draw depends on its seed, not its turn
    cols = ["reference", "scope", "matching", "group", "metric", key]
    assert_frame_equal(serial.sort(cols), again.sort(cols))
