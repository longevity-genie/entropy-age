"""RM17: one OSKM cycle, RRBS."""

import gzip

from entropy_age.experiments.reprogramming_rrbs import read_arms, read_cov


def test_read_arms_from_cell_type(tmp_path):
    path = tmp_path / "m.txt.gz"
    cells = ['"cell type: Old OSKM liver"', '"cell type: Young control liver"', '"cell type: Old control spleen"']
    with gzip.open(path, "wt") as h:
        h.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n')
        h.write('!Sample_source_name_ch1\t"Liver"\t"Liver"\t"Spleen"\n')
        h.write('!Sample_characteristics_ch1\t"tissue: Liver"\t"tissue: Liver"\t"tissue: Spleen"\n')
        h.write("!Sample_characteristics_ch1\t" + "\t".join(cells) + "\n")
    a = read_arms(path)
    assert a["arm"].to_list() == ["oskm", "young", "old"] and a["tissue"].to_list() == ["Liver", "Liver", "Spleen"]


def test_read_cov_sums_counts(tmp_path):
    path = tmp_path / "s.cov.gz"
    with gzip.open(path, "wt") as h:
        h.write("5\t100\t100\t75\t3\t1\nchrX\t7\t7\t0\t0\t12\n")
    c = read_cov(path)
    assert c["cpg"].to_list() == ["5:100", "X:7"] and c["coverage"].to_list() == [4, 12]
    assert abs(c["beta"][0] - 0.75) < 1e-6
