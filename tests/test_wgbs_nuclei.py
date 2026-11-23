from entropy_age.experiments.wgbs_nuclei import read_samples


def test_read_samples_parses_geo_text_records(tmp_path):
    text = """^SAMPLE = GSM1
!Sample_title = 1535_WGBS_NeuN
!Sample_characteristics_ch1 = cell type: NeuN
!Sample_characteristics_ch1 = subject id: 1535
!Sample_characteristics_ch1 = age: 52
!Sample_characteristics_ch1 = disease state: Control
!Sample_characteristics_ch1 = Sex: M
^SAMPLE = GSM2
!Sample_characteristics_ch1 = cell type: Olig2
!Sample_characteristics_ch1 = subject id: 1535
!Sample_characteristics_ch1 = age: 52
!Sample_characteristics_ch1 = disease state: Control
!Sample_characteristics_ch1 = Sex: M
"""
    path = tmp_path / "s.txt"
    path.write_text(text)
    frame = read_samples(path)
    assert frame["sample"].to_list() == ["GSM1", "GSM2"]
    assert frame["cell_type"].to_list() == ["NeuN", "Olig2"] and frame["donor"].to_list() == ["1535", "1535"]
    assert frame["age"].to_list() == [52.0, 52.0]
