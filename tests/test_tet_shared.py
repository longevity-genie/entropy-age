import numpy as np
import polars as pl

from entropy_age.experiments.tet_shared import base_ids, complete, effects


def test_base_ids_strips_suffix_and_keeps_first():
    wide = pl.DataFrame({"cpg": ["cg01_TC11", "cg01_BC21", "cg02", "cg03_TC12"], "a": [0.1, 0.9, 0.2, 0.3]})
    got = base_ids(wide)
    assert got["cpg"].to_list() == ["cg01", "cg02", "cg03"]
    assert got["a"].to_list() == [0.1, 0.2, 0.3]


def test_complete_drops_cpgs_with_missing_values():
    wide = pl.DataFrame({"cpg": ["a", "b", "c"], "s1": [0.1, None, 0.3], "s2": [0.2, 0.5, 0.4]})
    x, feats = complete(wide, ["s1", "s2"])
    assert feats == ["a", "c"]
    assert x.shape == (2, 2)


def test_effects_scores_only_the_kept_cpgs():
    rng = np.random.default_rng(0)
    features = [f"cg{i:03d}" for i in range(400)]
    x = np.clip(rng.normal(0.2, 0.05, (16, 400)), 0.01, 0.99)
    treated = np.arange(16) >= 8
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(400)]})
    members = pl.DataFrame(
        {
            "feature": features[:80],
            "group": ["housekeeping"] * 80,
            "gene": [f"g{i // 2}" for i in range(80)],
        }
    )
    design = np.column_stack([np.ones(16), treated.astype(float)])
    keep = features[:300]
    rows = effects(x, features, keep, treated, design, probes, members)
    assert {r["measure"] for r in rows} == {"B", "F", "S"}
    assert all(r["cpgs"] == 300 and r["set_cpgs"] == 80 for r in rows)
    assert all(np.isfinite(r["effect"]) for r in rows)
