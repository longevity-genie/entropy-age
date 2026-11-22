from entropy_age.experiments.mouse_atlas_rrbs import parse_names


def test_parse_names_reads_lane_mouse_age_sex_and_tissue():
    names = [
        "GSM7040480_lane1_ATCACGTT_M04-20525_3_M_Liver_iPCRtag1_merged_L001_R1.cov.gz",
        "GSM7040481_lane1_GAGCCAAT_M00017913_15_M_Heart_iTAG_82_merged_L001_R1.cov.gz",
        "not_a_library.txt",
    ]
    frame = parse_names(names)
    assert frame["sample"].to_list() == ["GSM7040480", "GSM7040481"]
    assert frame["mouse"].to_list() == ["M04-20525", "M00017913"]
    assert frame["age_months"].to_list() == [3.0, 15.0]
    assert frame["tissue"].to_list() == ["Liver", "Heart"] and frame["lane"].to_list() == ["lane1", "lane1"]
