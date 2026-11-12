"""RM19: accessibility spread with age, per peak, and peaks tied to 450k CpGs."""

import numpy as np
import polars as pl

from entropy_age.experiments.atac_variance import first_library, log_cpm, peak_cpgs, spread_slopes


def test_first_library_per_donor():
    cols = ["CR_027_2", "CR_027_1", "CR_028_2", "CR_100_1"]
    assert first_library(cols) == {"CR_027": "CR_027_1", "CR_028": "CR_028_2", "CR_100": "CR_100_1"}


def test_log_cpm_normalises_library_size():
    counts = np.array([[10.0, 90.0], [20.0, 180.0]])  # the second library is twice as deep
    y = log_cpm(counts)
    assert np.allclose(y[0], y[1])


def test_spread_slope_is_positive_where_people_spread_apart_with_age_and_not_otherwise():
    rng = np.random.default_rng(0)
    age = np.linspace(20, 80, 4000)
    flat = rng.normal(0, 1, age.size)
    growing = rng.normal(0, 1, age.size) * (age / 20)
    shifted = 0.05 * age + rng.normal(0, 1, age.size)  # mean trend only: removed with age
    slopes = spread_slopes(np.column_stack([flat, growing, shifted]), age, np.empty((age.size, 0)))
    assert slopes[1] > 0.005
    assert abs(slopes[0]) < 0.002 and abs(slopes[2]) < 0.002


def test_peak_cpgs_uses_half_open_intervals_per_chromosome():
    peaks = pl.DataFrame({"chr": ["chr1", "chr1", "chr2"], "start": [100, 500, 100], "end": [200, 600, 200]})
    manifest = pl.DataFrame(
        {"chr": ["chr1", "chr1", "chr1", "chr2", "chr3"], "pos": [100, 199, 200, 150, 150], "cpg": list("abcde")}
    )
    out = peak_cpgs(peaks, manifest).sort("cpg")
    assert dict(zip(out["cpg"], out["peak"].to_list(), strict=True)) == {"a": 0, "b": 0, "d": 2}
