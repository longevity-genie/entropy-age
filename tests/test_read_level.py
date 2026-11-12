"""RM22: epialleles from PAT reads, and per-window entropy and discordance."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.read_level import epiallele_counts, window_metrics


def test_reads_covering_a_window_whole_are_counted_by_pattern(tmp_path):
    pat = tmp_path / "s.pat.gz"
    with gzip.open(pat, "wt") as h:
        h.write("chr1\t10\tCCCCCC\t3\n")  # covers windows 10..13, 11..14, 12..15
        h.write("chr1\t11\tTTTT\t2\n")  # covers 11..14 only
        h.write("chr1\t11\tCT.T\t5\n")  # unknown call: not counted
        h.write("chr1\t12\tCTCT\t1\n")  # covers 12..15
    c = epiallele_counts(pat, np.array([10, 11, 12, 50]))
    assert c[0, 15] == 3 and c[0].sum() == 3
    assert c[1, 15] == 3 and c[1, 0] == 2 and c[1].sum() == 5
    assert c[2, 0b1010] == 1 and c[2, 15] == 3
    assert c[3].sum() == 0


def test_uniform_windows_have_no_entropy_and_mixed_ones_do():
    counts = np.zeros((3, 16), dtype=int)
    counts[0, 0] = 40  # all unmethylated: uniform
    counts[1, 0], counts[1, 15] = 20, 20  # two clean epialleles: entropy, no discordance
    counts[2, :] = 5  # every pattern: maximal entropy, mostly discordant
    m = window_metrics(counts)
    assert m["entropy"][0] == 0 and m["pdr"][0] == 0 and m["mean"][0] == 0
    assert m["entropy"][1] == pytest.approx(0.25) and m["pdr"][1] == 0
    assert m["entropy"][2] == pytest.approx(1.0) and m["pdr"][2] == pytest.approx(14 / 16)


def test_low_coverage_windows_are_missing():
    counts = np.zeros((1, 16), dtype=int)
    counts[0, 0] = 3
    assert np.isnan(window_metrics(counts)["entropy"][0])
