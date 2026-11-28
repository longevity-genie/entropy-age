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


def test_sesame_mode_is_refused(monkeypatch):
    from entropy_age.experiments.level_match_reference import deposited_only

    monkeypatch.setenv("ENTROPY_AGE_BETAS", "sesame")
    with pytest.raises(RuntimeError, match="deposited"):
        deposited_only()


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
    cohorts |= {k: (p := _prepared(30, 2, 9, k, features), p.x) for k in ("charge", "young_adults")}
    rows = {r["key"]: r for r in cohort_summary(cohorts, 2)}
    assert set(rows) == {ref for refs in PLAN.values() for ref, _ in refs.values()}
    assert rows["hannum"]["reference_rule"] == "youngest fifth" and rows["hannum"]["reference_arrays"] == 10
    assert rows["young_adults"]["reference_rule"] == "all donors" and rows["young_adults"]["reference_arrays"] == 30
    assert "BMIQ" in rows["young_adults"]["processing"]
