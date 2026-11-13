"""RM2: Meer 2018 sample sheet."""

import gzip

from entropy_age.experiments.meer_tissues import read_mice


def test_read_mice_parses_months(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\n')
        h.write('!Sample_source_name_ch1\t"Heart"\t"Cortex"\n')
        h.write('!Sample_characteristics_ch1\t"strain: C57BL/6"\t"strain: C57BL/6"\n')
        h.write('!Sample_characteristics_ch1\t"age,months: 6"\t"age,months: 30"\n')
    m = read_mice(path)
    assert m["age_months"].to_list() == [6.0, 30.0] and m["tissue"].to_list() == ["Heart", "Cortex"]
