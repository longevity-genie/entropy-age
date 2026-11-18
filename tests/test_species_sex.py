import polars as pl

from entropy_age.experiments.species_sex import sex_of, summarise


def test_sex_codes():
    assert [sex_of(v) for v in ("1", "0", "Female", "male", None, "NA")] == ["F", "M", "F", "M", None, None]


def test_summary_pairs_sexes_within_groups():
    rows = []
    for i in range(8):
        for sex, rate, gap in (("F", 2.0, -0.3), ("M", 1.0, -0.1)):
            rows.append(
                {
                    "species": f"s{i}",
                    "tissue": "Blood",
                    "metric": "kl_slope",
                    "group": "housekeeping",
                    "sex": sex,
                    "background_rate": rate,
                    "gap_pooled_sd": gap,
                }
            )
    s = summarise(pl.DataFrame(rows))
    r = s.row(0, named=True)
    assert r["groups"] == 8 and r["share_rate_higher_in_f"] == 1.0 and r["share_gap_more_negative_in_f"] == 1.0
    assert abs(r["median_log_rate_ratio_f_m"] - 0.6931) < 1e-3
