"""RM12: divergence after a transplant bottleneck."""

import numpy as np
import polars as pl

from entropy_age.experiments.hsct_bottleneck import KINDS, comparisons, divergence_scores, kind_effect, parse_ids


def test_parse_and_compare_within_batch():
    ids = ["1_D_1_1", "1_R_1_1", "1_R_1_2", "1_R_2_1", "2_D_1_1", "2_D_1_2"]
    s = parse_ids(ids).with_columns(years=pl.Series([0.0, 1.0, 3.0, 5.0, 0.0, 0.0]))
    c = comparisons(s)
    kinds = c.group_by("kind").len().sort("kind")
    got = dict(kinds.iter_rows())
    assert got[KINDS[0]] == 2 and got[KINDS[1]] == 1 and got[KINDS[2]] == 1
    rr = c.filter(pl.col("kind") == KINDS[1]).row(0, named=True)
    assert (rr["first"], rr["second"], rr["years"]) == (1, 2, 2.0)


def test_divergence_share_is_negative_when_set_diverges_less():
    rng = np.random.default_rng(1)
    x = rng.uniform(0.1, 0.9, size=(1, 200)).repeat(2, axis=0)
    in_set = np.zeros(200, bool)
    in_set[:50] = True
    step = rng.normal(scale=0.05, size=200)
    step[in_set] *= 0.3
    x[1] += step
    out = divergence_scores(x, np.array([0]), np.array([1]), {"housekeeping": (in_set, [])}, np.zeros(200).astype(str))
    assert out["F_housekeeping"][0] < -0.5


def test_kind_effect_reads_the_bottleneck_term():
    n = 40
    kind = np.array([KINDS[0], KINDS[1]] * (n // 2))
    years = np.tile(np.arange(n // 2, dtype=float), 2)
    f = -0.2 * (kind == KINDS[0]) + 0.01 * np.log1p(years)
    t = pl.DataFrame({"kind": kind, "years": years, "batch": ["1"] * n, "F": f})
    assert abs(kind_effect(t, "F")[0] + 0.2) < 1e-9
