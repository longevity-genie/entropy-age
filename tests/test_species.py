"""RM10: the consortium's sample sheet and the per-group gap."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.species_shield import group_gap, read_samples


def test_the_sample_sheet_reads_species_lifespan_and_barcode_whatever_row_a_key_sits_in(tmp_path):
    path = tmp_path / "GSE223748_series_matrix.txt.gz"
    lines = [
        '!Sample_title\t"genomic DNA from Mouse Blood 202897220093_R01C01"'
        '\t"genomic DNA from Naked mole rat Skin 203000000001_R02C01"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"',
        '!Sample_source_name_ch1\t"Blood"\t"Skin"',
        '!Sample_characteristics_ch1\t"age: 0.5"\t"age: 12"',
        '!Sample_characteristics_ch1\t"speciescommonname: Mouse"\t"tissue: Skin"',  # GEO shifts a key to another row
        '!Sample_characteristics_ch1\t"maximumlifespanyears: 4"\t"speciescommonname: Naked mole rat"',
        '!Sample_characteristics_ch1\t"average_adultweight: 25"\t"maximumlifespanyears: 31"',
        '!Sample_characteristics_ch1\t""\t"average_adultweight: 35"',
    ]
    with gzip.open(path, "wt") as h:
        h.write("\n".join(lines) + "\n")
    s = read_samples(path).sort("sample")
    assert s["species"].to_list() == ["Mouse", "Naked mole rat"]
    assert s["max_lifespan"].to_list() == [4.0, 31.0]
    assert s["barcode"].to_list() == ["202897220093_R01C01", "203000000001_R02C01"]
    assert s["tissue"].to_list() == ["Blood", "Skin"]
    assert "experiment" in s.columns


def test_a_set_that_drifts_less_gives_a_negative_gap_with_an_interval():
    rng = np.random.default_rng(0)
    n, p = 40, 600
    ages = rng.uniform(0, 10, n)
    in_set = np.arange(p) < 150
    rate = np.where(in_set, 0.3, 1.0) * 0.01
    x = np.clip(0.1 + np.abs(rng.normal(0, 1, (n, p))) * rate * ages[:, None], 0, 1)
    genes = [np.arange(j, j + 3) for j in range(0, 150, 3)]
    rows = group_gap(
        x,
        ages,
        [f"c{j}" for j in range(p)],
        np.zeros(p, dtype=bool),
        {"housekeeping": (in_set, genes)},
        50,
        np.random.default_rng(1),
    )
    d = {r["metric"]: r for r in rows}
    assert d["dispersion_slope"]["gap_sd"] < -0.5
    assert d["dispersion_slope"]["ci_low"] < d["dispersion_slope"]["gap_sd"] < d["dispersion_slope"]["ci_high"]


def test_background_rate_is_the_lookalikes_drift_and_gap_is_unscaled():
    # set CpGs drift at 0.3 of the background; the background's dispersion slope sets background_rate
    rng = np.random.default_rng(2)
    n, p = 60, 800
    ages = rng.uniform(0, 10, n)
    in_set = np.arange(p) < 200
    genes = [np.arange(j, j + 4) for j in range(0, 200, 4)]
    out = {}
    for scale in (1.0, 2.0):
        rate = np.where(in_set, 0.3, 1.0) * 0.01 * scale
        x = np.clip(0.1 + np.abs(rng.normal(0, 1, (n, p))) * rate * ages[:, None], 0, 1)
        rows = group_gap(
            x,
            ages,
            [f"c{j}" for j in range(p)],
            np.zeros(p, dtype=bool),
            {"housekeeping": (in_set, genes)},
            10,
            np.random.default_rng(3),
        )
        out[scale] = {r["metric"]: r for r in rows}["dispersion_slope"]
    assert out[1.0]["background_rate"] > 0
    assert out[2.0]["background_rate"] / out[1.0]["background_rate"] == pytest.approx(2.0, rel=0.25)
    assert out[1.0]["gap"] < 0 and out[1.0]["gap"] / out[1.0]["background_rate"] == pytest.approx(-0.7, abs=0.2)


def test_outcomes_give_log_rate_and_relative_gap_and_drop_nonpositive_rates():
    import polars as pl

    from entropy_age.experiments.species_shield import outcomes

    sel = pl.DataFrame({"gap_sd": [-0.3, -0.2], "gap": [-0.02, -0.01], "background_rate": [0.1, -0.05]})
    out = outcomes(sel)
    assert out["log background rate"][0] == pytest.approx(np.log(0.1)) and np.isnan(out["log background rate"][1])
    assert out["relative gap"][0] == pytest.approx(-0.2) and np.isnan(out["relative gap"][1])
