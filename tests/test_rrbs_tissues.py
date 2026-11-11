"""RM2: GSE120132's CGfinal tables and sample sheet, and the six-tissue gaps."""

import gzip

import numpy as np
import polars as pl
import pytest

from entropy_age.data.rrbs import read_cgfinal, read_thompson_samples
from entropy_age.experiments.common import scaled_intervals
from entropy_age.experiments.mouse_rrbs_tissues import tissue_gaps


def test_cgfinal_keeps_chromosomes_and_covered_cpgs(tmp_path):
    path = tmp_path / "x.CGfinal.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write("chr1 cG 3014928 3 4 0.750 0.7\n")
        h.write("chrUn_JH584304 CG 100 1 2 0.5 0.7\n")
        h.write("chrX Cg 500 0 0 0 0\n")
        h.write("chrM CG 10 5 10 0.5 0.7\n")
    frame = read_cgfinal(path)
    assert frame.rows() == [("1:3014928", 0.75, 4)]


def test_ages_must_be_in_months(tmp_path):
    def matrix(age: str):
        path = tmp_path / f"m{len(age)}.txt.gz"
        lines = [
            '!Sample_title\t"Mouse_Blood_1"',
            '!Sample_geo_accession\t"GSM1"',
            '!Sample_characteristics_ch1\t"strain: C57BL/6J"',
            '!Sample_characteristics_ch1\t"tissue: Blood"',
            f'!Sample_characteristics_ch1\t"age: {age}"',
            '!Sample_characteristics_ch1\t"Sex: Male"',
            '!Sample_supplementary_file_1\t"ftp://ftp.ncbi.nlm.nih.gov/x.gz"',
        ]
        with gzip.open(path, "wt") as h:
            h.write("\n".join(lines) + "\n")
        return path

    assert read_thompson_samples(matrix("11.1111mo"))["age_months"].to_list() == [pytest.approx(11.1111)]
    with pytest.raises(ValueError, match="not in months"):
        read_thompson_samples(matrix("40 weeks"))


def test_gaps_find_a_set_that_drifts_less():
    rng = np.random.default_rng(0)
    n, p = 60, 800
    ages = rng.uniform(2, 20, n)
    in_set = np.arange(p) < 200
    rate = np.where(in_set, 0.3, 1.0) * 0.01
    x = np.clip(
        0.1 + np.abs(rng.normal(0, 1, (n, p))) * rate * ages[:, None], 0, 1
    )  # drift away from 0.1 raises entropy
    features = [f"c{j}" for j in range(p)]
    members = pl.DataFrame({"feature": features[:200], "group": "housekeeping"})
    cov = np.column_stack([rng.integers(0, 2, n), rng.integers(0, 2, n)]).astype(float)
    strata = np.array(["s"] * p)
    rows = {r["metric"]: r["matched_diff"] / r["sd"] for r in tissue_gaps(x, ages, cov, features, strata, members)}
    assert rows["dispersion_slope"] < -0.5 and rows["entropy_slope"] < 0


def test_intervals_are_scaled_by_the_full_samples_spread():
    point = pl.DataFrame({"metric": ["m"], "group": ["g"], "matched_diff": [-0.2], "sd": [1.0]})
    # resamples: the unscaled gap scatters around -0.2, while each resample's own spread is inflated to 2
    boot = pl.DataFrame({"metric": "m", "group": "g", "matched_diff": np.linspace(-0.3, -0.1, 41), "sd": 2.0})
    row = scaled_intervals(point, boot, ["metric", "group"]).row(0, named=True)
    assert row["gap_sd"] == pytest.approx(-0.2)
    assert row["ci_low"] < row["gap_sd"] < row["ci_high"]
    assert row["ci_low"] == pytest.approx(-0.295, abs=0.01)
