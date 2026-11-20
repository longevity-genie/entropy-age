import gzip

from entropy_age.experiments.orbitofrontal import groups, read_people


def test_read_people_drops_duplicates_and_groups_by_status(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\n')
        handle.write('!Sample_characteristics_ch1\t"age: 30"\t"age: 60"\t"age: 61"\t"age: 70"\n')
        handle.write('!Sample_characteristics_ch1\t"Sex: Male"\t"Sex: Female"\t"Sex: Male"\t"Sex: Male"\n')
        handle.write(
            '!Sample_characteristics_ch1\t"technical_duplicate: no"\t"technical_duplicate: no"\t'
            '"technical_duplicate: yes"\t"technical_duplicate: no"\n'
        )
        handle.write(
            '!Sample_characteristics_ch1\t"disease_status: 0"\t"disease_status: 1"\t"disease_status: 1"\t'
            '"disease_status: 1"\n'
        )
    people = read_people(path)
    assert people["sample"].to_list() == ["GSM1", "GSM2", "GSM4"]
    g = groups(people)
    assert g["status 0"]["sample"].to_list() == ["GSM1"]
    assert g["status 1"]["sample"].to_list() == ["GSM2", "GSM4"]
    assert g["everyone, status held"].height == 3
