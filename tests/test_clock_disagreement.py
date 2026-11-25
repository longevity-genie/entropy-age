import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.clock_disagreement import (
    PANELS,
    age_clocks,
    aligned,
    breast_rows,
    disagreement,
    fit_clogit,
    fit_logit,
    partial_r,
    quality_terms,
    read_clock_scores,
    summarise,
)


def _clocks(rng, n=200):
    age = rng.uniform(20, 80, n)
    shared = rng.normal(0, 3, n)
    noise = rng.normal(0, 2, (n, 4))
    preds = age[:, None] + shared[:, None] + noise
    return age, preds


def test_disagreement_is_sd_of_z_scored_age_residuals():
    rng = np.random.default_rng(0)
    age, preds = _clocks(rng)
    dis, mean, use = disagreement(preds, age)
    assert use.all()
    d = np.column_stack([np.ones(age.size), age])
    acc = preds - d @ np.linalg.lstsq(d, preds, rcond=None)[0]
    z = (acc - acc.mean(axis=0)) / acc.std(axis=0)
    np.testing.assert_allclose(dis, z.std(axis=1, ddof=1))
    np.testing.assert_allclose(mean, z.mean(axis=1))


def test_disagreement_ignores_each_clocks_offset_scale_and_age_slope():
    rng = np.random.default_rng(1)
    age, preds = _clocks(rng)
    dis, mean, _ = disagreement(preds, age)
    moved = preds * np.array([1.0, 0.1, 5.0, 2.0]) + np.array([0.0, 40.0, -3.0, 0.0]) + np.outer(age, [0, 0, 0, 0.5])
    dis2, mean2, _ = disagreement(moved, age)
    np.testing.assert_allclose(dis, dis2, atol=1e-10)
    np.testing.assert_allclose(mean, mean2, atol=1e-10)


def test_agreeing_clocks_give_zero_and_one_outlier_clock_gives_spread():
    rng = np.random.default_rng(2)
    age = rng.uniform(20, 80, 100)
    shared = rng.normal(0, 3, 100)
    preds = np.column_stack([age + shared, 2 * (age + shared) + 1, age + shared - 5])
    dis, _, _ = disagreement(preds, age)
    np.testing.assert_allclose(dis, 0, atol=1e-10)
    preds[0, 0] += 30  # one clock reads person 0 much older
    dis, mean, _ = disagreement(preds, age)
    assert dis[0] > 1 and np.argmax(dis) == 0 and mean[0] > 0


def test_clock_against_age_is_left_out():
    rng = np.random.default_rng(3)
    age, preds = _clocks(rng)
    preds[:, 3] = -age + rng.normal(0, 2, age.size)
    dis, _, use = disagreement(preds, age)
    assert list(use) == [True, True, True, False]
    full, _, use_all = disagreement(preds, age, min_age_r=None)
    assert use_all.all() and not np.allclose(dis, full)


def test_too_few_clocks_gives_nan():
    rng = np.random.default_rng(4)
    age, preds = _clocks(rng)
    dis, mean, _ = disagreement(preds[:, :1], age)
    assert np.isnan(dis).all() and np.isnan(mean).all()


OUTPUTS = {
    "Horvathv1": "Age (Years)",
    "Hannum": "Age (Years)",
    "PhenoAge": "Age (Years)",
    "GrimAgeV2": "Mortality Adjusted Age (Years)",
    "Zhang_10": "Mortality Risk",
    "DunedinPACE": "Aging Rate (Years/Year)",
    "EpiTOC1": "Stem Cell Division Rate",
    "StocH": "Age (Years)",
    "BMI_McCartney": "BMI",
    "DNAmClockCortical": "Human Cortex Age (Years)",
    "Horvathv2": "Age (Years)",
    "Bohlin": "Age (days)",
}


def test_panel_keeps_chronological_year_clocks_only():
    # risk, pace, mitotic, stochastic, trait, phenotypic-target and day-unit outputs are out of the primary panel
    assert age_clocks(OUTPUTS, "blood") == ["Horvathv1", "Hannum", "Horvathv2"]
    assert age_clocks(OUTPUTS, "prefrontal cortex") == ["Horvathv1", "DNAmClockCortical", "Horvathv2"]
    assert age_clocks(OUTPUTS, "kidney", PANELS[1]) == ["Horvathv1", "Hannum", "Horvathv2"]
    assert age_clocks(OUTPUTS, "kidney", PANELS[2]) == ["Horvathv1", "Horvathv2"]
    assert age_clocks(OUTPUTS, "blood", PANELS[3]) == ["Horvathv1", "Hannum", "PhenoAge", "GrimAgeV2", "Horvathv2"]
    assert age_clocks({"Horvathv1": "Age (Years)"}, "blood", PANELS[2]) == []  # the fixed panel or nothing


def test_near_constant_residual_is_left_out():
    rng = np.random.default_rng(5)
    age, preds = _clocks(rng)
    preds[:, 2] = 3 * age + 1 + 1e-12 * rng.normal(size=age.size)  # a clock that is exactly linear in age
    _, _, use = disagreement(preds, age)
    assert list(use) == [True, True, False, True]


def test_quality_covariates_change_the_accelerations():
    rng = np.random.default_rng(6)
    age, preds = _clocks(rng)
    batch = rng.normal(0, 1, age.size)
    preds[:, 0] += 5 * batch  # one clock moved by array quality
    plain, _, _ = disagreement(preds, age)
    adjusted, _, _ = disagreement(preds, age, covariates=batch[:, None])
    assert np.corrcoef(plain, np.abs(batch))[0, 1] > np.corrcoef(adjusted, np.abs(batch))[0, 1] + 0.2


def test_quality_terms_report_unknown_intensity_and_keep_order():
    samples = [f"s{i}" for i in range(40)]
    qc = pl.DataFrame(
        {
            "sample": samples[::-1],
            "chip": [str(i % 4) for i in range(40)],
            "frac_dt": np.linspace(0.95, 0.99, 40),
            "mean_intensity": [None] + [5000.0] * 39,
        }
    )
    x, terms, intensity = quality_terms(samples, qc)
    assert intensity == "unknown" and "intensity unknown" in terms and "chip (4 chips)" in terms
    np.testing.assert_allclose(x[:, 0], np.linspace(0.95, 0.99, 40)[::-1])  # joined in the order of samples
    known = qc.with_columns(mean_intensity=pl.Series(np.linspace(4000, 6000, 40)))
    x, terms, intensity = quality_terms(samples, known)
    assert intensity == "known" and terms.startswith("chip") and x.shape[1] == 2 + 3
    x, terms, intensity = quality_terms(samples, None)
    assert x.shape[1] == 0 and "no QC table" in terms


def test_chip_dropped_when_it_leaves_too_few_degrees_of_freedom():
    samples = [f"s{i}" for i in range(20)]
    qc = pl.DataFrame({"sample": samples, "chip": [str(i // 2) for i in range(20)], "frac_dt": np.linspace(0.9, 1, 20)})
    x, terms, _ = quality_terms(samples, qc)
    assert "chip dropped" in terms and x.shape[1] == 1


def test_aligned_keeps_sample_order_and_refuses_repeats(tmp_path):
    long = pl.DataFrame(
        {
            "sample": ["b", "a", "b", "a"],
            "clock": ["H1", "H1", "H2", "H2"],
            "output": ["Age (Years)"] * 4,
            "prediction": [2.0, 1.0, 20.0, 10.0],
        }
    )
    np.testing.assert_allclose(aligned(["a", "b"], long, ["H2", "H1"]), [[10.0, 1.0], [20.0, 2.0]])
    path = tmp_path / "x_clocks.csv"
    pl.concat([long, long.head(1)]).with_columns(pl.col("sample").replace({"a": "007"})).write_csv(path)
    with pytest.raises(ValueError, match="repeated"):
        read_clock_scores(path)
    long.with_columns(sample=pl.Series(["007", "1", "007", "1"])).write_csv(path)
    read, _ = read_clock_scores(path)
    assert read["sample"].dtype == pl.String and "007" in read["sample"].to_list()  # ids stay strings


def test_partial_r_matches_residual_correlation():
    rng = np.random.default_rng(7)
    z = rng.normal(size=300)
    x, y = z + rng.normal(size=300), z + rng.normal(size=300)
    assert np.corrcoef(x, y)[0, 1] > 0.3 and abs(partial_r(x, y, z)) < 0.15


def test_logit_diagnostics():
    rng = np.random.default_rng(8)
    x = rng.normal(size=400)
    y = (rng.uniform(size=400) < 1 / (1 + np.exp(-(0.5 * x - 1)))).astype(float)
    design = np.column_stack([np.ones(400), x])
    beta, se, status = fit_logit(y, design)
    assert status == "ok" and abs(beta[1] - 0.5) < 0.3
    separated = (x > 0).astype(float)
    beta, se, status = fit_logit(separated, design)
    assert status != "ok" and np.isnan(beta).all() and np.isnan(se).all()
    beta, _, status = fit_logit(y, np.column_stack([design, 2 * x]))
    assert status == "design not of full rank" and np.isnan(beta).all()
    assert fit_logit(np.zeros(400), design)[2] == "one outcome only"


def test_clogit_diagnostics():
    rng = np.random.default_rng(9)
    sets = [[3 * i, 3 * i + 1, 3 * i + 2] for i in range(60)]
    x = rng.normal(size=(180, 1))
    beta, _, status = fit_clogit(x, sets)
    assert status == "ok" and np.isfinite(beta).all()
    for s in sets:  # the case always has the highest value: separation
        x[s[0], 0] = x[s, 0].max() + 1
    beta, _, status = fit_clogit(x, sets)
    assert status != "ok" and np.isnan(beta).all()


def _women(rng, n=600):
    chips = np.array([f"c{i % 30}" for i in range(n)])
    icd = ["C50.9" if u < 0.15 else None for u in rng.uniform(size=n)]
    frame = pl.DataFrame(
        {
            "sample": [f"GSM{i}" for i in range(n)],
            "age": rng.uniform(40, 70, n),
            "chip": chips,
            "icd": icd,
            "years_to_diagnosis": rng.uniform(0, 10, n),
            "mean_beta": rng.normal(0.5, 0.01, n),
            "mean_intensity": rng.uniform(3000, 6000, n),
            "frac_dt": rng.uniform(0.95, 0.99, n),
            "disagreement": rng.gamma(2, 0.4, n),
            "mean_acceleration": rng.normal(size=n),
            **{f"fraction_{i}": rng.uniform(size=n) for i in range(1, 3)},
        }
    )
    return frame


def test_breast_rows_report_every_model_with_retention():
    rng = np.random.default_rng(10)
    women = _women(rng)
    rows = pl.DataFrame(breast_rows(women, ["fraction_1", "fraction_2"], ["Horvathv1", "Hannum"]))
    models = set(rows["model"])
    assert {"covariates, chip effects", "matched sets, same chip only, covariates", "covariates + quality"} <= models
    exact = rows.filter(pl.col("model") == "matched sets, same chip only, covariates")
    assert (exact["controls_on_case_chip"] == 1.0).all()
    assert (exact["cases"] <= exact["eligible_cases"]).all()
    chip = rows.filter(pl.col("model") == "covariates, chip effects")
    assert (chip["cases"] <= chip["eligible_cases"]).all() and (chip["chips"] <= 30).all()
    assert (rows["fit"] == "ok").all() and rows["odds_ratio"].is_finite().all()
    assert set(rows["estimand"]) == {rows["estimand"][0]} and "incident breast cancer" in rows["estimand"][0]


def test_summary_leaves_out_sensitivity_and_subset_rows():
    drift = pl.DataFrame(
        {
            "panel": ["valid"] * 4,
            "independent": [True, True, False, False],
            "r_dis_B": [0.1, 0.3, 0.9, 0.9],
            "r_dis_F": [0.0] * 4,
            "partial_r_dis_B": [0.1, 0.3, 0.9, 0.9],
            "partial_r_dis_B_age": [0.1, 0.3, 0.9, 0.9],
            "partial_r_dis_B_quality": [0.05, float("nan"), 0.9, 0.9],
            "partial_r_dis_F_quality": [0.0] * 4,
        }
    )
    row = summarise(drift).row(0, named=True)
    assert row["groups"] == 2 and np.isclose(row["median_r_dis_B"], 0.2)
    assert row["groups_with_quality"] == 1 and np.isclose(row["median_partial_r_dis_B_quality"], 0.05)
    assert np.isclose(row["median_partial_r_dis_B_age_same_groups"], 0.1)
