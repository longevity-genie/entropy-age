"""series_table reads a samples CSV for arrays outside GEO (TCGA normal-adjacent)."""

import polars as pl

from entropy_age.experiments.series import series_table


def test_samples_csv_gives_the_sample_table(tmp_path):
    csv = tmp_path / "samples.csv"
    pl.DataFrame(
        {
            "sample": ["TCGA-A-11A", "TCGA-B-11A", "TCGA-C-11A"],
            "basename": ["u1_noid", "u2_noid", "u3_noid"],
            "case_id": ["c1", "c2", "c3"],
            "age": [61.2, None, 70.0],
            "sex": ["male", "female", "male"],
        }
    ).write_csv(csv)
    spec = {"samples_csv": str(csv), "age": "age", "sex": "sex", "donor": "case_id", "groups": {"kidney": {}}}
    table = series_table(None, spec)
    assert table["sample"].to_list() == ["u1_noid", "u3_noid"]  # array names of the rebuild; no age, dropped
    assert table["donor"].to_list() == ["c1", "c3"]
    assert series_table(None, {**spec, "samples_csv": str(tmp_path / "missing.csv")}) is None
