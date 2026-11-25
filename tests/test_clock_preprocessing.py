"""Idea 5: clock predictions on the same arrays under two preprocessing tracks."""

import numpy as np
import polars as pl

from entropy_age.experiments.clock_preprocessing import compare, read_scores, read_status, summarise


def _long(samples, clock, output, values):
    return pl.DataFrame({"sample": samples, "clock": clock, "output": output, "prediction": values})


def _tracks(seed=0, n=40):
    rng = np.random.default_rng(seed)
    samples = [f"s{i}" for i in range(n)]
    age = rng.uniform(20, 80, n)
    horvath_m = age + rng.normal(0, 4, n)
    horvath_u = horvath_m + 1.5 + rng.normal(0, 1, n)  # shifted by 1.5 years, noisy
    pace_m = rng.normal(1, 0.1, n)
    pace_u = pace_m.copy()  # identical on both tracks
    masked = pl.concat([_long(samples, "Horvathv1", "Age (Years)", horvath_m), _long(samples, "Pace", "Rate", pace_m)])
    unmasked = pl.concat(
        [
            # one array only on the unmasked track, one array with a non-finite prediction
            _long(samples + ["extra"], "Horvathv1", "Age (Years)", np.append(horvath_u, 50.0)),
            _long(samples, "Pace", "Rate", np.where(np.arange(n) == 0, np.nan, pace_u)),
            _long(samples, "GrimAgeV1", "Mortality Adjusted Age (Years)", age + 3),
        ]
    )
    ages = pl.DataFrame({"sample": samples[:-5], "age": age[:-5]})  # five arrays without an age
    return masked, unmasked, ages, (age, horvath_m, horvath_u, pace_m)


def test_compare_matches_direct_computation():
    masked, unmasked, ages, (age, hm, hu, pm) = _tracks()
    rows = {r["clock"]: r for r in compare(masked, unmasked, ages, {"GrimAgeV1": "coverage"}, {})}
    h = rows["Horvathv1"]
    d = hu - hm
    sd = (np.std(hm, ddof=1) + np.std(hu, ddof=1)) / 2
    assert h["arrays"] == len(hm)  # the array scored on one track only is left out
    assert np.isclose(h["r"], np.corrcoef(hm, hu)[0, 1])
    assert np.isclose(h["mean_diff"], d.mean()) and h["mean_diff"] > 1  # unmasked minus masked
    assert np.isclose(h["sd_diff"], np.std(d, ddof=1))
    assert np.isclose(h["mean_diff_per_sd"], d.mean() / sd)
    assert np.isclose(h["sd_diff_per_sd"], np.std(d, ddof=1) / sd)
    assert h["arrays_age"] == len(hm) - 5
    assert np.isclose(h["mae_masked"], np.mean(np.abs(hm[:-5] - age[:-5])))
    assert np.isclose(h["mae_unmasked"], np.mean(np.abs(hu[:-5] - age[:-5])))
    assert not h["one_track_only"]

    p = rows["Pace"]
    assert p["arrays"] == len(pm) - 1  # the NaN array drops out
    assert np.isclose(p["r"], 1.0) and np.isclose(p["mean_diff"], 0.0) and np.isclose(p["sd_diff_per_sd"], 0.0)
    assert p["mae_masked"] is None  # not an age clock

    g = rows["GrimAgeV1"]
    assert g["arrays"] == 0 and g["r"] is None
    assert (g["status_masked"], g["status_unmasked"]) == ("coverage", "kept") and g["one_track_only"]


def test_not_scored_is_not_a_coverage_difference():
    masked, unmasked, ages, _ = _tracks()
    rows = {r["clock"]: r for r in compare(masked, unmasked, ages, {"GrimAgeV1": "GrimageModel needs META"}, {})}
    assert not rows["GrimAgeV1"]["one_track_only"]


def test_summarise_per_clock():
    rows = []
    for seed in range(3):
        masked, unmasked, ages, _ = _tracks(seed)
        rows += [{"series": f"S{seed}"} | r for r in compare(masked, unmasked, ages, {"GrimAgeV1": "coverage"}, {})]
    pairs = pl.DataFrame(rows, infer_schema_length=None)
    clocks = {r["clock"]: r for r in summarise(pairs).iter_rows(named=True)}
    h = pairs.filter(pl.col("clock") == "Horvathv1")
    assert clocks["Horvathv1"]["series_compared"] == 3
    assert np.isclose(clocks["Horvathv1"]["median_r"], float(np.median(h["r"].to_numpy())))
    assert np.isclose(clocks["Horvathv1"]["median_abs_mean_diff"], float(np.median(np.abs(h["mean_diff"].to_numpy()))))
    assert clocks["GrimAgeV1"]["unmasked_only"] == 3 and clocks["GrimAgeV1"]["series_compared"] == 0


def test_read_scores_empty(tmp_path):
    empty = tmp_path / "X_clocks.csv"
    empty.write_text("")
    assert read_scores(empty) is None
    assert read_scores(tmp_path / "missing_clocks.csv") is None


def test_coverage_on_common_arrays_replaces_each_tracks_own_panel():
    masked, unmasked, ages, _ = _tracks()
    own = {"GrimAgeV1": "coverage", "Pace": "kept"}
    rows = {
        r["clock"]: r
        for r in compare(
            masked,
            unmasked,
            ages,
            own,
            {"Horvathv1": "kept"},
            arrays={f"s{i}" for i in range(20)},
            common_masked={"GrimAgeV1": "kept", "Horvathv1": "kept", "Pace": "kept"},
            common_unmasked={"GrimAgeV1": "kept", "Horvathv1": "kept", "Pace": "coverage"},
        )
    }
    g = rows["GrimAgeV1"]  # one track only on the panels' own arrays, kept on both on the common arrays
    assert g["one_track_only_own_arrays"] and not g["one_track_only"]
    assert (g["status_masked"], g["status_masked_own_arrays"]) == ("kept", "coverage")
    assert rows["Pace"]["one_track_only"] and not rows["Pace"]["one_track_only_own_arrays"]
    assert rows["Horvathv1"]["arrays"] == 20 and rows["Horvathv1"]["status_basis"] == "common arrays"


def test_not_scored_reason_survives_common_coverage():
    masked, unmasked, ages, _ = _tracks()
    rows = {
        r["clock"]: r
        for r in compare(
            masked,
            unmasked,
            ages,
            {"GrimAgeV1": "GrimageModel needs META"},
            {},
            common_masked={"GrimAgeV1": "coverage"},
            common_unmasked={"GrimAgeV1": "kept"},
        )
    }
    assert rows["GrimAgeV1"]["status_masked"] == "GrimageModel needs META"
    assert not rows["GrimAgeV1"]["one_track_only"]


def test_read_status_finds_the_panel_beside_any_scores_name(tmp_path):
    for name in ("GSE1_clocks.csv", "GSE1_clocks_unmasked.csv", "GSE1_qc712_clocks.csv"):
        panel = tmp_path / name.replace(".csv", "_panel.csv")
        pl.DataFrame({"clock": ["Horvathv1"], "status": ["kept"]}).write_csv(panel)
        assert read_status(tmp_path / name) == {"Horvathv1": "kept"}
