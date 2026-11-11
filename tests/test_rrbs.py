"""RM3: the GSE80672 RRBS tables, mm10 promoters and islands, and the intervention effect."""

import gzip

import numpy as np
import polars as pl
import pytest

from entropy_age.data.rrbs import cpg_genes, in_intervals, promoters, read_overlap, read_samples
from entropy_age.experiments.interventions import effect


def test_a_sample_table_maps_refseq_chromosomes_and_drops_other_contigs(tmp_path):
    path = tmp_path / "GSM1_x.overlap.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write("index\tX_A.Percentage\tX_A.Coverage\n")
        h.write("gi|149247747|ref|NT_166280.1|:23455\t81.8\t11\n")  # an unplaced contig
        h.write("gi|372099109|ref|NC_000067.6|:3000100\t25\t12\n")  # chromosome 1
        h.write("gi|372099090|ref|NC_000086.7|:500\t100\t30\n")  # X
    frame = read_overlap(path)
    assert frame["cpg"].to_list() == ["1:3000100", "X:500"]
    assert frame["beta"].to_list() == pytest.approx([0.25, 1.0])


def test_promoters_follow_the_strand(tmp_path):
    refseq = tmp_path / "ncbiRefSeq.txt.gz"
    with gzip.open(refseq, "wt") as h:
        h.write("0\tNM_1\tchr1\t+\t10000\t20000\n")
        h.write("0\tNM_2\tchr2\t-\t10000\t20000\n")
        h.write("0\tNM_3\tchrUn_x\t+\t10000\t20000\n")
    link = tmp_path / "ncbiRefSeqLink.txt.gz"
    with gzip.open(link, "wt") as h:
        h.write("NM_1\t\t\t\t\t\t111\n")
        h.write("NM_2\t\t\t\t\t\t222\n")
        h.write("NM_3\t\t\t\t\t\t333\n")
    prom = promoters(refseq, link).sort("entrez")
    assert prom.rows() == [("1", 8500, 10501, 111), ("2", 19499, 21500, 222)]


def test_intervals_may_overlap():
    # the second interval lies inside the first, as nested promoters of one gene's transcripts do
    intervals = pl.DataFrame({"chr": ["1", "1", "2"], "start": [10, 15, 0], "end": [40, 20, 5]})
    chroms = np.array(["1", "1", "1", "1", "2", "3"])
    pos = np.array([9, 10, 30, 40, 4, 4])
    assert in_intervals(chroms, pos, intervals).tolist() == [False, True, True, False, True, False]


def test_cpgs_take_every_promoter_that_holds_them():
    cpgs = pl.DataFrame({"cpg": ["1:101", "1:150", "1:400"], "chr": ["1", "1", "1"], "pos": [101, 150, 400]})
    prom = pl.DataFrame({"chr": ["1", "1"], "start": [100, 140], "end": [160, 200], "entrez": [7, 8]})
    links = cpg_genes(cpgs, prom).sort("feature", "entrez")
    assert links.rows() == [("1:101", 7), ("1:150", 7), ("1:150", 8)]


def test_ages_are_read_as_months_with_their_leading_number(tmp_path):
    path = tmp_path / "GSE80672_series_matrix.txt.gz"
    lines = [
        '!Sample_title\t"20D02"\t"F1"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"',
        '!Sample_characteristics_ch1\t"tissue: whole blood"\t"tissue: Lung fibroblasts"',
        '!Sample_characteristics_ch1\t"strain: C57BL/6"\t"strain: C57BL/6"',
        '!Sample_characteristics_ch1\t"age (years): 27"\t"age (years): 2.5 (at isolation)"',
        '!Sample_characteristics_ch1\t"gender: Male"\t"gender: Male"',
        '!Sample_characteristics_ch1\t"genetic condition: Wild Type"\t"genetic condition: Wild Type"',
        '!Sample_characteristics_ch1\t"diet: standard"\t"diet: standard"',
        '!Sample_supplementary_file_1\t"ftp://ftp.ncbi.nlm.nih.gov/a.gz"\t"ftp://ftp.ncbi.nlm.nih.gov/b.gz"',
    ]
    with gzip.open(path, "wt") as h:
        h.write("\n".join(lines) + "\n")
    samples = read_samples(path)
    assert samples["age_months"].to_list() == [27.0, 2.5]
    assert samples["url"][0] == "https://ftp.ncbi.nlm.nih.gov/a.gz"


def test_effect_recovers_a_treatment_shift_with_age_held_fixed():
    rng = np.random.default_rng(0)
    age = rng.uniform(10, 28, 60)
    treated = np.arange(60) < 30
    age[treated] += 3  # treated mice are older on average: age must be held fixed
    y = 0.02 * age - 0.3 * treated + rng.normal(0, 0.05, 60)
    row = effect(pl.DataFrame({"age_months": age, "y": y}), "y", treated, True)
    assert row["effect"] == pytest.approx(-0.3, abs=0.06) and row["p"] < 1e-6
