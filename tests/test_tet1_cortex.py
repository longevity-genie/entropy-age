"""RM12: the long Tet1 table as bisulfite samples by probe, and genotypes from animal names."""

import polars as pl

from entropy_age.experiments.tet1_cortex import genotype, wide_bisulfite


def test_wide_keeps_bisulfite_only_and_names_animals():
    long = pl.DataFrame(
        {
            "ProbeID": ["cg1", "cg1", "cg2", "cg2", "cg1"],
            "Betas": [0.1, 0.9, 0.2, 0.8, 0.3],
            "Sample": ["WT1_BS", "WT1_bACE", "WT1_BS", "WT1_bACE", "KO1_BS"],
        }
    )
    w = wide_bisulfite(long).sort("cpg")
    assert set(w.columns) == {"cpg", "WT1", "KO1"}
    assert w["WT1"].to_list() == [0.1, 0.2]
    assert w["KO1"].to_list()[0] == 0.3


def test_genotype():
    assert genotype("WT3") == "WT" and genotype("Hxd7") == "Hxd" and genotype("KO12") == "KO"
