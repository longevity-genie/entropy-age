"""RM28 SATSA: consecutive pairs, divergence separated from the measurement floor, within-person slopes."""

import numpy as np
import polars as pl

from entropy_age.experiments.within_person import consecutive_pairs, pair_metrics, person_slopes, within_slope


def test_pairs_follow_each_person_in_age_order():
    person = np.array(["a", "b", "a", "a", "b", "c"])
    age = np.array([60.0, 50.0, 50.0, 70.0, 55.0, 40.0])
    first, second = consecutive_pairs(person, age)
    pairs = sorted(zip(first.tolist(), second.tolist(), strict=True))
    assert pairs == [(1, 4), (2, 0), (0, 3)] or set(pairs) == {(2, 0), (0, 3), (1, 4)}


def test_divergence_is_found_beside_a_measurement_floor():
    rng = np.random.default_rng(3)
    n_people, draws = 400, 4
    person = np.repeat(np.arange(n_people), draws)
    age = np.concatenate([np.sort(rng.uniform(50, 80, draws)) for _ in range(n_people)])
    q = np.array([0.0, 0.002])  # two CpGs: no walk, a random walk
    walk = np.zeros((person.size, 2))
    for p in range(n_people):
        rows = np.flatnonzero(person == p)
        steps = rng.normal(size=(draws - 1, 2)) * np.sqrt(np.diff(age[rows])[:, None] * q)
        walk[rows] = np.vstack([np.zeros(2), np.cumsum(steps, axis=0)])
    x = 0.5 + 0.003 * (age[:, None] - 65) + walk + rng.normal(scale=0.02, size=walk.shape)
    first, second = consecutive_pairs(person.astype(str), age)
    m = pair_metrics(x, age, np.zeros(age.size, bool), first, second)
    assert abs(m["directed"][0] - 0.003) < 5e-4
    assert abs(m["divergence"][0]) < 3e-4 and abs(m["divergence"][1] - 0.002) < 6e-4
    assert abs(m["floor"][0] - 2 * 0.02**2) < 3e-4


def test_within_slope_ignores_differences_between_people():
    rng = np.random.default_rng(5)
    person = np.repeat(np.arange(200), 3)
    age = rng.uniform(50, 90, person.size)
    level = rng.normal(scale=5, size=200)[person]  # people differ a lot, unrelated to age
    v = level + 0.1 * age + rng.normal(scale=0.1, size=person.size)
    assert abs(within_slope(v, age, person, np.zeros(person.size, bool)) - 0.1) < 0.01


def test_person_slopes_report_within_and_between():
    rng = np.random.default_rng(8)
    person = np.repeat(np.arange(100).astype(str), 2)
    age = np.concatenate([[a, a + 8] for a in rng.uniform(50, 80, 100)])
    people = pl.DataFrame({"person": person, "age": age, "epic": np.tile([False, True], 100)})
    rows = person_slopes(people, {"m": 0.05 * age + rng.normal(scale=0.2, size=age.size)}, rng)
    assert rows[0]["within_per_decade"] > 0 and rows[0]["low"] > 0 and rows[0]["people"] == 100


def test_read_sdrf_keeps_green_channel(tmp_path):
    from entropy_age.experiments.within_person import read_sdrf

    head = "Source Name\tCharacteristics[individual]\tCharacteristics[age]\tCharacteristics[sex]\tLabel\tAssay Name\n"
    rows = [
        "S1\tIND_1\t60.5\tfemale\tCy3\t97_R01C01_Grn\n",
        "S1\tIND_1\t60.5\tfemale\tCy5\t97_R01C01_Red\n",
        "S2\tIND_1\t70.1\tfemale\tCy3\t98_R02C01_Grn\n",
    ]
    path = tmp_path / "s.sdrf.txt"
    path.write_text(head + "".join(rows))
    s = read_sdrf(path)
    assert s["sample"].to_list() == ["97_R01C01", "98_R02C01"] and s["person"].to_list() == ["IND_1", "IND_1"]
    assert not s["epic"].any() and s["female"].all()
