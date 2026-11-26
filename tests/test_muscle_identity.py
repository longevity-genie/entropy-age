import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.muscle_identity import (
    axis_projection,
    band_gaps,
    differential,
    identity_gene_cpgs,
    matched_partners,
    null_projection,
    slope_interval,
)


def test_axis_projection_places_the_ends_at_zero_and_one():
    origin = np.array([0.1, 0.9, 0.5, 0.2])
    target = np.array([0.8, 0.2, 0.5, 0.6])
    x = np.vstack([origin, target, (origin + target) / 2, origin + 0.25 * (target - origin)])
    assert np.allclose(axis_projection(x, origin, target), [0.0, 1.0, 0.5, 0.25])


def test_axis_projection_ignores_missing_values():
    origin, target = np.array([0.1, 0.9, 0.3]), np.array([0.8, 0.2, np.nan])
    x = np.array([[0.8, 0.2, 0.0], [np.nan, 0.2, 0.3]])
    assert np.allclose(axis_projection(x, origin, target), [1.0, 1.0])


def test_differential_needs_threshold_and_confirmation_of_sign():
    a = np.array([0.9, 0.9, 0.9, 0.5, np.nan])
    b = np.array([0.1, 0.1, 0.1, 0.4, 0.1])
    confirm = np.array([0.5, -0.5, 0.05, 0.5, 0.5])
    assert differential(a, b, 0.3).tolist() == [True, True, True, False, False]
    assert differential(a, b, 0.3, confirm, 0.1).tolist() == [True, False, False, False, False]


def test_identity_gene_cpgs_keeps_promoters_and_enhancers_of_listed_genes():
    links = pl.DataFrame(
        {
            "cpg": ["c1", "c2", "c3", "c4", "c1"],
            "symbol": ["DES", "DES", "DES", "GAPDH", "DES"],
            "region": ["TSS200", "Body", "Body", "TSS200", "Body"],
        }
    )
    out = identity_gene_cpgs(links, {"c2", "c4"})
    assert out["cpg"].to_list() == ["c1", "c2"]
    assert dict(out.select("cpg", "kind").iter_rows()) == {"c1": "promoter", "c2": "enhancer"}


def test_slope_interval_recovers_a_known_slope():
    rng = np.random.default_rng(1)
    age = rng.uniform(20, 80, 200)
    sex = rng.integers(0, 2, 200).astype(float)
    y = 0.003 * age + 0.1 * sex + rng.normal(0, 0.02, 200)
    r = slope_interval(y, age, sex[:, None], 200, np.random.default_rng(2))
    assert abs(r["slope"] - 0.003) < 0.0005
    assert r["ci_low"] < 0.003 < r["ci_high"] and r["p"] < 1e-10


def test_null_projection_follows_the_matched_partners():
    origin = np.array([0.1, 0.9, 0.1, 0.9])
    target = np.array([0.6, 0.4, 0.1, 0.9])
    x = np.array([[0.2, 0.9, 0.6, 0.9]])  # the first background CpG moved fully toward the target direction
    strata = np.array(["a", "b", "a", "b"])
    in_axis = np.array([True, True, False, False])
    partners = matched_partners(strata, in_axis, ~in_axis, 5, np.random.default_rng(0))
    assert (partners == [[2, 3]] * 5).all()
    assert np.allclose(null_projection(x, origin, target, partners, np.flatnonzero(in_axis)), [0.5])


def test_band_gaps_sees_a_set_that_drifts_faster_than_matched_cpgs():
    rng = np.random.default_rng(3)
    n = 4000
    baseline = rng.uniform(0.05, 0.95, n)
    in_set = np.zeros(n, dtype=bool)
    in_set[:400] = True
    metrics = pl.DataFrame(
        {m: rng.normal(0, 1, n) + 2.0 * in_set for m in ("kl_slope", "dispersion_slope", "entropy_slope")}
    )
    metrics = metrics.with_columns(baseline=pl.Series(baseline))
    rows = band_gaps(metrics, rng.integers(0, 2, n).astype(bool), {"set": in_set})
    gap = next(r for r in rows if r["band"] == "all" and r["metric"] == "kl_slope")["gap_sd"]
    assert gap > 1.0


# ---------------------------------------------------------------- review 7
def test_relaxed_qc_applies_the_intensity_rule_after_the_lower_cut():
    from entropy_age.experiments.muscle_identity import relaxed_qc

    qc = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d", "d"],
            "chip": [1, 1, 2, 2, 2],
            "frac_dt": [0.95, 0.85, 0.82, 0.70, 0.99],
            "mean_intensity": [5000.0, 5100.0, 900.0, 5000.0, 5000.0],
        }
    )
    out = relaxed_qc(qc, 0.80)
    assert dict(out.select("sample", "passed").iter_rows()) == {"a": True, "b": True, "c": False, "d": False}
    assert out.filter(pl.col("sample") == "c")["low_intensity"].item()  # dim array dropped, as in the 90 % set


def test_chip_identifiability_flags_age_aliased_with_chip():
    from entropy_age.experiments.muscle_identity import chip_identifiability

    ages = np.array([20.0, 21, 22, 60, 61, 62] * 4)
    held = np.zeros((ages.size, 0))
    aliased = np.array(["c1", "c1", "c1", "c2", "c2", "c2"] * 4)
    r = chip_identifiability(ages, np.char.add(aliased, np.repeat(["a", "b", "c", "d"], 6)), held)
    assert r["age_identifiable_with_chip"] == "yes" and r["chip_age_r2"] > 0.9
    per_age = np.array([f"chip{a:.0f}" for a in ages])  # one chip per age: age is the chip
    assert chip_identifiability(ages, per_age, held)["age_identifiable_with_chip"] == "no"
    assert chip_identifiability(ages, np.array(["NA"] * ages.size), held)["age_identifiable_with_chip"] == "unknown"


def test_shared_panel_needs_every_platform_and_one_sign():
    from entropy_age.experiments.muscle_identity import shared_panel

    origin = np.array([0.1, 0.1, 0.1, 0.9, 0.5])
    t450 = np.array([0.8, 0.8, 0.2, 0.1, 0.9])
    tepic = np.array([0.7, 0.3, 0.8, 0.2, 0.1])
    confirm = np.array([-0.5, -0.5, -0.5, 0.5, -0.5])
    assert shared_panel(origin, [t450, tepic], 0.3, confirm, 0.1).tolist() == [True, False, False, True, False]


def test_holdout_rows_score_held_out_bodies_and_calibrate_admixture():
    from entropy_age.experiments.muscle_identity import holdout_rows

    rng = np.random.default_rng(0)
    n_cpg = 200
    origin = rng.uniform(0.05, 0.3, n_cpg)
    adipo = origin + 0.6
    muscle = origin + rng.normal(0, 0.01, (6, n_cpg))
    fat = adipo - 0.1 + rng.normal(0, 0.01, (6, n_cpg))
    rows = holdout_rows(muscle, fat, origin, {"450k": adipo, "EPIC": adipo}, "450k", 2, np.random.default_rng(1))
    assert len(rows) == 2 and all(r["panel_cpgs"] == n_cpg for r in rows)
    r = rows[0]
    assert abs(r["muscle_position"]) < 0.05 and r["fat_position"] > 0.7 and r["fat_above_muscle"] == 1.0
    step = r["position_10pct_fat"] - r["position_0pct_fat"]
    assert step == pytest.approx(0.1 * (r["fat_position"] - r["muscle_position"]), rel=1e-6)


def test_crossfit_kl_keeps_a_young_sample_out_of_its_reference():
    from entropy_age.experiments.muscle_identity import crossfit_kl
    from entropy_age.stochasticity import relative_entropy

    x = np.array([[0.2], [0.4], [0.9], [0.9], [0.9]])
    ages = np.array([20.0, 21.0, 60.0, 61.0, 62.0])
    kl = crossfit_kl(x, ages, q=0.4)  # young: the two youngest
    assert kl[0, 0] == pytest.approx(relative_entropy(np.array([0.2]), np.array([0.4]))[0])
    assert kl[2, 0] == pytest.approx(relative_entropy(np.array([0.9]), np.array([0.3]))[0])
    assert kl[0, 0] > 0  # an in-sample reference would have pulled it toward its own value
    boot = crossfit_kl(x[[0, 0, 1, 2, 3]], ages[[0, 0, 1, 2, 3]], ids=np.array([0, 0, 1, 2, 3]), q=0.6)
    assert boot[0, 0] == pytest.approx(kl[0, 0])  # both copies of sample 0 left out


def test_weighted_matched_gap_weights_genes_and_matches_within_class():
    from entropy_age.experiments.muscle_identity import weighted_matched_gap

    values = np.array([3.0, 3.0, 3.0, 1.0, 0.0, 0.0, 2.0, 2.0])
    member = np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)
    strata = np.array(["p", "p", "p", "e", "p", "p", "e", "e"])
    # three promoter CpGs of one gene (+3 each) and one enhancer CpG of another gene (-1)
    per_cpg = np.where(member, 1.0, 0.0)
    per_gene = np.where(member, [1 / 3, 1 / 3, 1 / 3, 1.0, 0, 0, 0, 0], 0.0)
    assert weighted_matched_gap(values, member, per_cpg, strata) == pytest.approx((3 * 3 - 1) / 4)
    assert weighted_matched_gap(values, member, per_gene, strata) == pytest.approx((3 - 1) / 2)


def test_regulatory_class_of_cpgs():
    from entropy_age.experiments.muscle_identity import regulatory_class

    links = pl.DataFrame({"cpg": ["a", "b", "c"], "symbol": ["X", "Y", "Z"], "region": ["TSS200", "Body", "Body"]})
    assert regulatory_class(["a", "b", "c", "d"], links, {"a", "b"}).tolist() == [
        "promoter",
        "enhancer",
        "other",
        "other",
    ]
