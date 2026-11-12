"""wgbstools-style CpG index from a FASTA, and reading a .beta file."""

import numpy as np

from entropy_age.data.cpg_index import build_index, cpg_positions, read_beta


def test_cpg_positions_are_one_based_and_case_blind():
    assert cpg_positions("ACGTcgA").tolist() == [2, 5]


def test_index_runs_across_chromosomes_in_natural_order(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr2 x\nCGCG\n>chr1\nAC\nGT\n>chrUn_x\nCGCG\n>chr10\nCG\n")
    idx = build_index(fa, chromosomes=("chr1", "chr2", "chr10"))
    assert idx["chr"].to_list() == ["chr1", "chr2", "chr2", "chr10"]
    assert idx["pos"].to_list() == [2, 1, 3, 1]  # chr1's CG spans a line break
    assert idx["index"].to_list() == [1, 2, 3, 4]


def test_beta_holds_two_counts_per_cpg(tmp_path):
    path = tmp_path / "s.beta"
    np.array([[3, 4], [0, 0], [10, 10]], dtype=np.uint8).tofile(path)
    meth, cov = read_beta(path)
    assert meth.tolist() == [3, 0, 10] and cov.tolist() == [4, 0, 10]
