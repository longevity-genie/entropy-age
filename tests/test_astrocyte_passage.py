import polars as pl

from entropy_age.experiments.astrocyte_passage import bulk_cultures


def test_bulk_cultures_keeps_unsorted_cultures_with_doublings():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d"],
            "FACS.code": ["Bulk", "Sen", "Bulk", "Bulk"],
            "donor": ["1", "1", "2", "1"],
            "cpd": ["3.5", "9", "NA", "0"],
        }
    )
    kept = bulk_cultures(chars)
    assert kept["sample"].to_list() == ["d", "a"]
    assert kept["doublings"].to_list() == [0.0, 3.5]
