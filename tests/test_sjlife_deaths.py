"""RM28 SJLIFE: first draws per person, and death models with the treatment covariates."""

import numpy as np
import polars as pl

from entropy_age.experiments.sjlife_deaths import EXPOSURES, death_rows, first_draws


def test_first_draw_is_the_youngest_sample_of_each_person():
    p = pl.DataFrame({"title": ["A1", "A1_1", "A1_2", "B2"], "age": [30.0, 25.0, 35.0, 40.0]})
    f = first_draws(p).sort("person")
    assert f["title"].to_list() == ["A1_1", "B2"]


def test_death_model_finds_a_marker_beyond_age_and_treatment():
    rng = np.random.default_rng(21)
    n = 3000
    age = rng.uniform(20, 60, n)
    chest = (rng.uniform(size=n) < 0.3).astype(float)
    marker = rng.normal(size=n)
    risk = -4 + 0.04 * (age - 40) + 1.0 * chest + 0.6 * marker
    died = rng.uniform(size=n) < 1 / (1 + np.exp(-risk))
    people = pl.DataFrame(
        {
            "age": age,
            "female": rng.uniform(size=n) < 0.5,
            "survivor": rng.uniform(size=n) < 0.9,
            "diagnosis": rng.choice(["ALL", "HL", "control"], n),
            **{e: (chest if e == EXPOSURES[0] else rng.integers(0, 2, n).astype(float)) for e in EXPOSURES},
            "age_at_death": np.where(died, age + rng.uniform(0.5, 10, n), np.nan),
        }
    ).with_columns(pl.col("age_at_death").fill_nan(None))
    rows = {
        (r["subset"], r["lag_over_2y"], r["measure"]): r
        for r in death_rows(people, {"m": marker, "null": rng.normal(size=n)})
    }
    assert rows[("everyone", False, "m")]["odds_ratio"] > 1.5 and rows[("everyone", False, "m")]["p"] < 1e-6
    assert rows[("everyone", False, "null")]["p"] > 0.001
    assert rows[("everyone", True, "m")]["deaths"] < rows[("everyone", False, "m")]["deaths"]


def test_streaming_reader_picks_rows_and_columns(tmp_path):
    import gzip

    from entropy_age.experiments.sjlife_deaths import stream_header, stream_rows

    path = tmp_path / "t.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write("CpG_ID,b1,b2,b3\ncg1,0.1,0.2,0.3\ncg2,0.4,NA,0.6\ncg3,0.7,0.8,0.9\n")
    assert stream_header(path) == ["CpG_ID", "b1", "b2", "b3"]
    names, x = stream_rows(path, {"cg3", "cg2"}, ["b3", "b2"])
    assert names == ["cg2", "cg3"] and x.shape == (2, 2)
    assert x[0].tolist() == [0.6, 0.9] and np.isnan(x[1, 0]) and x[1, 1] == 0.8
