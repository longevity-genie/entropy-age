import numpy as np

from entropy_age.experiments.mouse_tissues import MouseTissue, chip_overlap, shared_chips


def test_shared_chips_and_overlap_table():
    groups = np.array(["Aged_ctrl", "Aged_ctrl", "Young_ctrl", "Young_ctrl", "Aged_ctrl", "Young_ctrl"])
    chips = np.array(["c1", "c1", "c2", "c2", "c3", "c3"])
    aged = groups == "Aged_ctrl"
    assert shared_chips(chips, aged).tolist() == ["c3"]
    t = MouseTissue(
        "spleen",
        np.zeros((6, 2)),
        groups,
        np.array(list("abcdef")),
        ["x", "y"],
        chips,
        np.array(["R01C01", "R02C01", "R01C01", "R02C01", "R03C01", "R04C01"]),
    )
    row = chip_overlap(t)
    assert row["shared_chips"] == 1 and row["aged_on_shared_chips"] == 1 and row["young_on_shared_chips"] == 1
    assert 0 < row["r2_age_on_chips"] < 1
    full = MouseTissue("spleen", np.zeros((4, 2)), groups[:4], np.array(list("abcd")), ["x", "y"], chips[:4])
    assert chip_overlap(full)["r2_age_on_chips"] == 1.0  # age fully confounded with chip
