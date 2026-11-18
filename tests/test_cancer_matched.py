import numpy as np
import polars as pl

from entropy_age.experiments.cancer_matched import clogit, match


def test_match_prefers_same_chip_and_uses_controls_once():
    cases = pl.DataFrame({"age": [50.0, 51.0], "chip": ["a", "b"]})
    controls = pl.DataFrame({"age": [50.5, 50.0, 51.0, 70.0], "chip": ["b", "a", "b", "a"]})
    sets = match(cases, controls, k=1, caliper=3.0)
    assert len(sets) == 2
    flat = [j for s in sets for j in s[1:]]
    assert len(flat) == len(set(flat)) and 2 + 3 not in flat  # the 70-year-old is out of the caliper
    by_case = {s[0]: s[1] for s in sets}
    assert by_case[1] == 2 + 2  # case on chip b, age 51: the same-chip 51-year-old
    assert by_case[0] == 2 + 1  # case on chip a: the same-chip control


def test_clogit_recovers_a_known_effect():
    rng = np.random.default_rng(1)
    rows, sets = [], []
    for _ in range(1500):
        x = rng.normal(size=4)
        w = np.exp(0.7 * x)
        case = rng.choice(4, p=w / w.sum())
        order = [case, *[j for j in range(4) if j != case]]
        sets.append(list(range(len(rows), len(rows) + 4)))
        rows += [x[j] for j in order]
    beta, se = clogit(np.array(rows)[:, None], sets)
    assert abs(beta[0] - 0.7) < 3 * se[0]
