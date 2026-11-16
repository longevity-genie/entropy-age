"""RM4: Tabula Muris Senis bulk reader."""

import gzip

from entropy_age.experiments.rna_tms import read_tms


def test_read_tms_parses_tissue_age_and_counts(tmp_path):
    meta = tmp_path / "meta.csv.gz"
    with gzip.open(meta, "wt") as f:
        f.write("Sample name,title,source name,organism,characteristics: age,stage,characteristics: sex\n")
        f.write("A1_P1_S1,t,Liver_24,m,12,months,m\n")
        f.write("A2_P1_S2,t,Brain_3,m,1,months,f\n")  # too young: left out
        f.write("A3_P1_S3,t,Small_Intestine_7,m,24,months,f\n")
    counts = tmp_path / "counts.csv.gz"
    with gzip.open(counts, "wt") as f:
        f.write("gene,A1_P1_S1.gencode.vM19,A2_P1_S2.gencode.vM19,A3_P1_S3.gencode.vM19\n")
        f.write("Actb,10,20,30\nGapdh,1,2,3\n")
    s, genes, x = read_tms(meta, counts)
    assert s["tissue"].to_list() == ["Liver", "Small_Intestine"]
    assert s["age"].to_list() == [12.0, 24.0] and genes == ["Actb", "Gapdh"]
    assert x.shape == (2, 2) and x[1, 0] == 30
