import gzip

import polars as pl

from entropy_age.experiments.sturm_expression import control_cultures, control_rna, gene_means


def test_control_cultures_keeps_untreated_healthy_ambient_lines_with_enough_samples():
    rows = [("HC1", "control", "Ambient", float(i)) for i in range(6)]
    rows += [("HC2", "control", "Ambient", 1.0)]  # too few
    rows += [("HC1", "DEX", "Ambient", 3.0), ("HC1", "control", "3%", 3.0), ("SURF1_1", "control", "Ambient", 3.0)]
    rows += [("HC1", "control", "Ambient", None)]
    s = pl.DataFrame(
        {
            "sample": [f"s{i}" for i in range(len(rows))],
            "donor": [r[0] for r in rows],
            "agent": [r[1] for r in rows],
            "oxygen": [r[2] for r in rows],
            "doublings": [r[3] for r in rows],
        }
    )
    kept = control_cultures(s)
    assert kept["sample"].to_list() == [f"s{i}" for i in range(6)]


def test_control_rna_and_gene_means(tmp_path):
    chars = pl.DataFrame(
        {
            "rnaseq_sampleid": ["1", "2", "3"],
            "treatments": ["Control", "DEX", "Control"],
            "clinical_condition": ["Normal", "Normal", "SURF1_Mutation"],
            "percent_oxygen": ["21", "21", "21"],
        }
    )
    assert control_rna(chars) == ["1"]
    path = tmp_path / "tpm.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write('"","1_lenscaleTPM","2_lenscaleTPM","10_lenscaleTPM"\n"A",3,100,7\n"B",0,5,1\n')
    means = gene_means(path, ["1", "10"])
    assert means["A"] == (2.0 + 3.0) / 2  # log2(4), log2(8)
    assert means["B"] == 0.5  # log2(1), log2(2)
