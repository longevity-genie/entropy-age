import numpy as np
import polars as pl

from entropy_age.experiments.syndrome_concordance import noise_concordance, site_classes


def test_classes_split_noisy_from_directed():
    n = 1000
    rng = np.random.default_rng(0)
    m = pl.DataFrame(
        {
            "feature": [str(i) for i in range(n)],
            "baseline": rng.uniform(size=n),
            "slope": rng.normal(size=n),
            "dispersion_slope": rng.normal(size=n),
        }
    )
    m = m.with_columns(abs_slope=pl.col("slope").abs())
    c = site_classes(m, top=0.05)
    counts = dict(c.group_by("cls").len().iter_rows())
    assert counts["noisy"] == 50 and 45 <= counts["directed"] <= 55
    d = c.join(m.select("feature", "dispersion_slope"), on="feature")
    assert d.filter(pl.col("cls") == "directed")["dispersion_slope"].max() <= np.median(m["dispersion_slope"])


def test_noise_concordance_is_matched_within_strata():
    log_vr = np.array([1.0, 0.0, 3.0, 2.0])
    noisy = np.array([True, False, True, False])
    strata = np.array(["a", "a", "b", "b"])
    assert noise_concordance(log_vr, noisy, strata) == 1.0
