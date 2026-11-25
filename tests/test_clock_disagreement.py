import numpy as np

from entropy_age.experiments.clock_disagreement import age_clocks, disagreement


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


def test_panel_keeps_valid_age_clocks_only():
    outputs = {
        "Horvathv1": "Age (Years)",
        "Hannum": "Age (Years)",
        "GrimAgeV2": "Age (Years)",
        "Zhang_10": "Mortality Risk",
        "DunedinPACE": "Aging Rate (Years/Year)",
        "EpiTOC1": "Stem Cell Division Rate",
        "StocH": "Age (Years)",
        "BMI_McCartney": "BMI",
        "DNAmClockCortical": "Age (Years)",
    }
    assert age_clocks(outputs, "blood") == ["Horvathv1", "Hannum", "GrimAgeV2", "Zhang_10"]
    assert age_clocks(outputs, "prefrontal cortex") == ["Horvathv1", "DNAmClockCortical"]
    assert age_clocks(outputs, "kidney", outside_tissue=True) == ["Horvathv1", "Hannum", "GrimAgeV2", "Zhang_10"]
