import polars as pl

from entropy_age.experiments.jena_rrbs import present_enough


def test_present_enough_keeps_cpgs_measured_in_nine_of_ten_libraries():
    cols = {f"s{i}": [0.1, 0.2, None if i < 2 else 0.3] for i in range(10)}
    betas = pl.DataFrame({"cpg": ["1:10", "1:20", "1:30"], **cols})
    assert present_enough(betas)["cpg"].to_list() == ["1:10", "1:20"]  # the third is missing in 2 of 10
