"""RM17: mice grouped into treated, old and young arms, with their beta columns."""

import gzip

from entropy_age.experiments.reprogramming import arm_of, read_mice

MATRIX = """!Sample_title\t"Sample1.Mouse.204529320080_R01C01"\t"Sample2.Mouse.204529320080_R02C01"
!Sample_geo_accession\t"GSM1"\t"GSM2"
!Sample_source_name_ch1\t"Liver"\t"Skin"
!Sample_characteristics_ch1\t"age: 1.75"\t"age: 0.25"
!Sample_characteristics_ch1\t"Sex: Male"\t"Sex: Female"
!Sample_characteristics_ch1\t"group: 3.4F+7mDox"\t"group: 7.4F-young"
"""


def test_mice_carry_barcode_column_tissue_and_group(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write(MATRIX)
    m = read_mice(path).sort("sample")
    assert m["barcode"].to_list() == ["X204529320080_R01C01", "X204529320080_R02C01"]
    assert m["tissue"].to_list() == ["Liver", "Skin"]
    assert m["group"].to_list() == ["3.4F+7mDox", "7.4F-young"]


def test_arms():
    assert arm_of("2.4F+1mDox") == "treated" and arm_of("4.4F+10mDox") == "treated"
    assert arm_of("5.4F-old") == "old" and arm_of("6.B6-old") == "old"
    assert arm_of("7.4F-young") == "young" and arm_of("1.B6+Dox") is None
