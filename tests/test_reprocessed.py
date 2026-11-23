def test_low_intensity_arrays_fail_qc():
    """An array passing detection at a tenth of the median signal is dropped and flagged (protocol change 4)."""
    import polars as pl

    from entropy_age.data.reprocessed import _intensity

    qc = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d"],
            "passed": [True, True, True, False],
            "mean_intensity": [5000.0, 4800.0, 450.0, 100.0],
        }
    )
    out = _intensity(qc)
    assert out["passed"].to_list() == [True, True, False, False]
    assert out["low_intensity"].to_list() == [False, False, True, False]
