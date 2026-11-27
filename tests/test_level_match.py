"""P6: finer level matching in the lowest decile (``experiments.level_match``)."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.level_match import (
    FLOORS,
    cpg_table,
    gaps,
    held_out_split,
    lowest_decile,
    nn_gap,
    strata,
    stratum_gap,
    summarise,
    young_rows,
    young_sd,
)
from entropy_age.experiments.tissue_contrast import LIFESPAN_METRICS, Prepared
from entropy_age.stochasticity import compare_groups


def _frame(n: int = 4000, seed: int = 0) -> tuple[pl.DataFrame, pl.DataFrame]:
    """CpGs whose drift is a function of the young level alone; set CpGs sit lower inside the lowest decile."""
    rng = np.random.default_rng(seed)
    level = rng.uniform(0.0, 1.0, n)
    low = level < np.quantile(level, 0.1)
    in_set = np.zeros(n, dtype=bool)
    lows = np.flatnonzero(low)
    in_set[lows[np.argsort(level[lows])[: lows.size // 2]]] = True  # the lower half of the lowest decile
    in_set[rng.choice(np.flatnonzero(~low), size=n // 10, replace=False)] = True
    drift = np.sqrt(level) + rng.normal(0, 1e-3, n)  # pure floor: less drift the lower the start
    frame = pl.DataFrame(
        {
            "feature": [f"c{i}" for i in range(n)],
            "baseline": level,
            "young_sd": 0.01 + 0.1 * level,
            "cpg_island": rng.random(n) < 0.5,
            "probe_type": np.where(rng.random(n) < 0.5, "I", "II"),
            **{m: drift for m in LIFESPAN_METRICS},
        }
    )
    membership = pl.DataFrame({"feature": frame["feature"].filter(pl.Series(in_set)), "group": "housekeeping"})
    return frame, membership


def test_fine_bins_split_only_the_lowest_decile_and_refine_it():
    frame, _ = _frame()
    coarse = strata(frame, by_type=True)
    fine = strata(frame, by_type=True, fine=20)
    low = lowest_decile(frame["baseline"].to_numpy())
    # every fine stratum lies inside one coarse stratum
    pairs = pl.DataFrame({"c": coarse, "f": fine}).unique()
    assert pairs["f"].n_unique() == pairs.height
    # outside the lowest decile nothing changes; inside it each island x probe cell gets up to 20 bins
    assert len(np.unique(fine[~low])) == len(np.unique(coarse[~low]))
    assert len(np.unique(fine[low])) == 20 * len(np.unique(coarse[low]))


def test_finer_bins_remove_the_level_residual_and_the_floor_gap():
    frame, membership = _frame()
    in_group = frame["feature"].is_in(membership["feature"].implode()).to_numpy()
    cols = ("kl_slope", "baseline")
    coarse = stratum_gap(frame, in_group, strata(frame, by_type=True), cols)
    fine = stratum_gap(frame, in_group, strata(frame, by_type=True, fine=50), cols)
    assert coarse["baseline"] < 0 and coarse["kl_slope"] < 0  # matched on decile, the set still starts lower
    assert abs(fine["baseline"]) < abs(coarse["baseline"]) / 3
    assert abs(fine["kl_slope"]) < abs(coarse["kl_slope"]) / 3
    assert fine["n_matched"] <= coarse["n_matched"] == in_group.sum()


def test_one_sided_strata_are_dropped():
    frame = pl.DataFrame({"x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], "baseline": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]})
    in_group = np.array([True, True, False, False, True, True])
    codes = np.array([0, 0, 0, 0, 1, 1])  # stratum 1 has no controls
    out = stratum_gap(frame, in_group, codes, ("x",))
    assert out["n_matched"] == 2
    assert out["x"] == pytest.approx(1.5 - 3.5)
    assert out["abs_level_diff"] == pytest.approx(np.mean([abs(0.1 - 0.35), abs(0.2 - 0.35)]))


def test_stratum_gap_equals_compare_groups_on_unequal_strata():
    rng = np.random.default_rng(5)
    n = 600
    codes = rng.choice(7, size=n, p=[0.4, 0.2, 0.15, 0.1, 0.08, 0.05, 0.02])
    in_group = rng.random(n) < (0.1 + 0.1 * codes)  # set share differs by stratum
    frame = pl.DataFrame(
        {"feature": [f"c{i}" for i in range(n)], "v": rng.normal(codes * 0.3, 1), "baseline": rng.random(n)}
    )
    mine = stratum_gap(frame, in_group, codes, ("v",))
    membership = pl.DataFrame({"feature": frame["feature"].filter(pl.Series(in_group)), "group": "s"})
    ref = compare_groups(
        frame.with_columns(stratum=pl.Series(codes.astype(str))),
        membership,
        metric="v",
        stratify_by="stratum",
        n_permutations=0,
    )
    assert mine["v"] == pytest.approx(ref["matched_diff"][0])
    assert mine["n_matched"] == ref["n_matched"][0]


def test_nearest_neighbours_match_a_brute_force_search_on_shuffled_rows():
    rng = np.random.default_rng(3)
    n = 400
    frame = pl.DataFrame(
        {
            "baseline": rng.uniform(0, 0.1, n),
            "young_sd": rng.uniform(0.0, 0.05, n),
            "cpg_island": rng.random(n) < 0.5,
            "probe_type": np.where(rng.random(n) < 0.5, "I", "II"),
            "v": rng.normal(0, 1, n),
        }
    )
    in_group = rng.random(n) < 0.3
    scale = {"baseline": 0.02, "young_sd": 0.01}  # the SD weighs more than the level: the two distances disagree
    for on in (("baseline",), ("baseline", "young_sd")):
        out = nn_gap(frame, in_group, ("v",), k=3, scale=scale, on=on)
        z = np.column_stack([frame[c].to_numpy() / scale[c] for c in on])
        cell = (frame["cpg_island"].cast(pl.String) + frame["probe_type"]).to_numpy()
        v, level, expected, absdiff = frame["v"].to_numpy(), frame["baseline"].to_numpy(), [], []
        for i in np.flatnonzero(in_group):
            ctrl = np.flatnonzero((cell == cell[i]) & ~in_group)
            nearest = ctrl[np.argsort(((z[ctrl] - z[i]) ** 2).sum(axis=1), kind="stable")[:3]]
            expected.append(v[i] - v[nearest].mean())
            absdiff.append(abs(level[i] - level[nearest].mean()))
        assert out["n_matched"] == in_group.sum()
        assert out["v"] == pytest.approx(np.mean(expected))
        assert out["abs_level_diff"] == pytest.approx(np.mean(absdiff))
    level_only = nn_gap(frame, in_group, ("young_sd",), k=3, scale=scale, on=("baseline",))
    both = nn_gap(frame, in_group, ("young_sd",), k=3, scale=scale, on=("baseline", "young_sd"))
    assert abs(both["young_sd"]) < abs(level_only["young_sd"]) or abs(level_only["young_sd"]) < 1e-3


def test_nearest_neighbours_skip_cells_with_too_few_controls():
    frame = pl.DataFrame(
        {
            "baseline": [0.1, 0.2, 0.3, 0.4],
            "young_sd": [0.01, 0.02, 0.03, 0.04],
            "cpg_island": [True, True, False, False],
            "probe_type": ["I", "I", "I", "I"],
            "kl_slope": [1.0, 2.0, 3.0, 4.0],
        }
    )
    out = nn_gap(frame, np.array([True, False, True, False]), ("kl_slope",), k=2)
    assert out["n_matched"] == 0 and np.isnan(out["kl_slope"])


def test_gaps_cover_every_scope_and_floors_keep_only_cpgs_above_them():
    frame, membership = _frame()
    rows = pl.DataFrame(gaps(frame, membership)).filter(pl.col("group") == "housekeeping")
    scopes = set(rows["scope"])
    assert scopes == {"all CpGs", "lowest decile", *(f"young level > {f:.2f}" for f in FLOORS)}
    is_set = frame["feature"].is_in(membership["feature"].implode())
    for f in FLOORS:
        n = rows.filter(pl.col("scope") == f"young level > {f:.2f}")["n"].unique().to_list()
        assert n == [int((is_set & (frame["baseline"] > f)).sum())]
    low = rows.filter(pl.col("scope") == "lowest decile")
    assert low["n"].unique().to_list() == [int((is_set & pl.Series(lowest_decile(frame["baseline"].to_numpy()))).sum())]


def test_young_sd_uses_the_youngest_fifth_and_cpg_table_needs_probe_design():
    rng = np.random.default_rng(1)
    ages = np.arange(20.0, 70.0)
    x = rng.uniform(0, 1, (ages.size, 6))
    young = ages <= np.quantile(ages, 0.2)
    assert np.allclose(young_sd(x, ages), x[young].std(axis=0, ddof=1))
    p = Prepared("t", x, ages, None, [f"c{i}" for i in range(6)], np.arange(ages.size), None)
    probes = pl.DataFrame({"cpg": [f"c{i}" for i in range(6)], "cpg_island": [True] * 6, "probe_type": ["I", "II"] * 3})
    table = cpg_table(p, probes)
    assert np.allclose(table["baseline"].to_numpy(), x[young].mean(axis=0))
    with pytest.raises(ValueError):
        cpg_table(p, probes.with_columns(probe_type=pl.Series(["I", None, "I", "II", "I", "II"])))


def test_held_out_split_keeps_donors_apart_and_draws_reference_from_the_young():
    rng = np.random.default_rng(2)
    ages = np.repeat(np.arange(20.0, 80.0), 2)  # each donor twice, as in a resample
    donors = np.repeat(np.arange(60), 2)
    ref, drift = held_out_split(ages, donors, rng)
    assert set(donors[ref]).isdisjoint(donors[drift])
    assert set(ref) | set(drift) == set(range(ages.size))
    young = set(donors[young_rows(ages)])
    assert set(donors[ref]) <= young and len(set(donors[ref])) == len(young) // 2


def test_non_set_pool_drops_other_sets_from_the_controls():
    frame, membership = _frame()
    other = frame["feature"].filter(~frame["feature"].is_in(membership["feature"].implode())).head(300)
    membership = pl.concat([membership, pl.DataFrame({"feature": other, "group": "core_essential"})])
    rows = pl.DataFrame(gaps(frame, membership)).filter(
        (pl.col("group") == "housekeeping")
        & (pl.col("scope") == "all CpGs")
        & (pl.col("matching") == "island x decile")
    )
    rest, non_set = (rows.filter(pl.col("pool") == p) for p in ("rest", "non-set"))
    assert rest["n"].item() == non_set["n"].item()
    assert rest["kl_slope"].item() != non_set["kl_slope"].item()  # 300 controls fewer changes the contrast


def test_summarise_ignores_undefined_draws_and_withholds_thin_intervals():
    keys = {"tissue": "t", "reference": "same donors", "pool": "rest", "scope": "all CpGs", "group": "g"}
    point = pl.DataFrame(
        [
            keys | {"matching": m, "metric": "kl_slope", "matched_diff": d, "level_diff": 0.0}
            for m, d in (("island x decile", -1.0), ("lowest decile in 50 bins", -0.5))
        ]
    )
    draws = [-1.0, -2.0, -3.0, float("nan"), float("inf")]
    boot = pl.DataFrame(
        [
            keys | {"matching": m, "metric": "kl_slope", "matched_diff": d * f, "level_diff": 0.0, "boot": b}
            for b, d in enumerate(draws)
            for m, f in (("island x decile", 1.0), ("lowest decile in 50 bins", 0.5))
        ]
    )
    spread = pl.DataFrame({"tissue": ["t"], "reference": ["same donors"], "metric": ["kl_slope"], "sd": [2.0]})
    out = summarise(point, boot, spread)
    fine = out.filter(pl.col("matching") == "lowest decile in 50 bins").row(0, named=True)
    assert fine["gap_sd"] == pytest.approx(-0.25) and fine["delta_sd"] == pytest.approx(0.25)
    assert fine["n_valid"] == 3 and fine["n_boot"] == 5
    assert fine["gap_lo"] == pytest.approx(np.quantile([-0.25, -0.5, -0.75], 0.025, method="nearest"))
    thin = summarise(point, boot.filter(pl.col("boot") >= 2), spread).row(0, named=True)
    assert thin["n_valid"] == 1 and thin["gap_lo"] is None
