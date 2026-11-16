"""RM35 item 1: loci held at one half by allele."""

import numpy as np
import polars as pl

from entropy_age.experiments.allelic_controls import half_methylated, strata


def test_half_methylated_keeps_candidates_near_one_half():
    m = pl.DataFrame({"feature": ["a", "b", "c", "d"], "baseline": [0.5, 0.9, 0.3, 0.45]})
    assert half_methylated(m, {"a", "b", "c"}) == {"a", "c"}


def test_strata_pair_island_relation_with_level():
    rng = np.random.default_rng(0)
    n = 400
    m = pl.DataFrame({"feature": [f"f{i}" for i in range(n)], "baseline": rng.random(n)})
    probes = pl.DataFrame({"cpg": [f"f{i}" for i in range(n)], "island_relation": ["Island", None] * (n // 2)})
    s = strata(m, probes)
    assert s["stratum"].n_unique() == 40
    assert s.filter(pl.col("island_relation").is_null())["stratum"].str.starts_with("OpenSea").all()


def test_sex_codes_read_every_spelling():
    from entropy_age.experiments.sex_drift import sex_codes

    assert sex_codes(["F", "FEMALE", "male", "Male", None, "unknown"]).tolist() == ["F", "F", "M", "M", "", ""]
