import gzip
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import probe_membership, read_gene_sets, read_manifest
from entropy_age.stochasticity import _matched_difference, _permute_within_strata, compare_groups

MANIFEST = """Illumina, Inc. 2007
[Heading],,
Loci Count,4,,
[Assay],,
IlmnID,Name,Chr,Gene_ID,Symbol,Accession,Distance_to_TSS,CPG_ISLAND,Product,,
cg01,cg01,1,GeneID:16,AARS,NM_001605.2,120,TRUE,"alanyl-tRNA synthetase, cytoplasmic",,
cg02,cg02,2,GeneID:999999,OLDNAME,NM_000001.1,800,FALSE,something,,
cg03,cg03,3,GeneID:7,ZZZ,NM_123456.3,50,TRUE,other,,
cg04,cg04,X,,,,,FALSE,,,
cg05,cg05,5,GeneID:55,OLDHK,NM_555555.1,30,TRUE,renamed housekeeping gene,,
[Controls],,
ctl,1,
"""


def test_manifest_and_membership(tmp_path: Path):
    path = tmp_path / "manifest.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(MANIFEST)
    probes = read_manifest(path)
    assert probes["cpg"].to_list() == ["cg01", "cg02", "cg03", "cg04", "cg05"]
    assert probes["cpg_island"].to_list() == [True, False, True, False, True]
    assert probes["refseq"].to_list()[:3] == ["NM_001605", "NM_000001", "NM_123456"]
    assert probes["entrez"].to_list() == [16, 999999, 7, None, 55]

    (tmp_path / "CEGv2.txt").write_text("GENE\tHGNC_ID\tENTREZ_ID\nAARS\tHGNC:20\t16\nRENAMED\tHGNC:1\t999999\n")
    (tmp_path / "NEGv1.txt").write_text("GENE\tHGNC_ID\tENTREZ_ID\nNOTONARRAY\tHGNC:2\t5\n")
    (tmp_path / "HK_genes.txt").write_text("ZZZ \tNM_999999\nNEWHK \tNM_555555\n")
    sets = read_gene_sets(tmp_path)
    membership = probe_membership(probes, sets)
    pairs = set(zip(membership["feature"].to_list(), membership["group"].to_list(), strict=True))
    # each key is the only route to one probe: Entrez -> cg02 (renamed), RefSeq -> cg05 (renamed), symbol -> cg03
    assert pairs == {
        ("cg01", "core_essential"),
        ("cg02", "core_essential"),
        ("cg03", "housekeeping"),
        ("cg05", "housekeeping"),
    }


def test_matched_difference_equals_per_stratum_loop():
    rng = np.random.default_rng(1)
    values = rng.normal(size=500)
    codes = rng.integers(0, 6, 500)
    labels = rng.random(500) < 0.3
    expected_num, expected_den = 0.0, 0
    for s in range(6):
        inside, outside = values[(codes == s) & labels], values[(codes == s) & ~labels]
        expected_num += inside.size * (inside.mean() - outside.mean())
        expected_den += inside.size
    assert np.isclose(_matched_difference(values, labels, codes, 6), expected_num / expected_den)


def test_categorical_strata_remove_a_confound():
    """CpG-island sites drift less; a gene set rich in islands looks protected unless matched on island status."""
    rng = np.random.default_rng(2)
    n = 2000
    island = rng.random(n) < 0.5
    drift = np.where(island, 0.1, 0.3) + rng.normal(0, 0.02, n)
    in_set = np.where(island, rng.random(n) < 0.3, rng.random(n) < 0.05)
    metrics = pl.DataFrame({"feature": [f"cg{i}" for i in range(n)], "drift": drift, "island": island})
    members = pl.DataFrame({"feature": metrics["feature"].filter(pl.Series(in_set)), "group": "set"})
    naive = compare_groups(metrics, members, metric="drift", stratify_by=None, n_permutations=200).row(0, named=True)
    matched = compare_groups(metrics, members, metric="drift", stratify_by="island", n_permutations=200).row(
        0, named=True
    )
    assert naive["matched_diff"] < -0.05 and naive["p_permutation"] < 0.01
    assert abs(matched["matched_diff"]) < 0.01 and matched["p_permutation"] > 0.05
    assert matched["n_matched"] == int(in_set.sum())


def test_permutation_keeps_every_stratum_composition():
    rng = np.random.default_rng(3)
    codes = rng.integers(0, 5, 1000)
    labels = rng.random(1000) < np.array([0.05, 0.2, 0.5, 0.8, 0.95])[codes]
    permuted = _permute_within_strata(labels, codes, rng)
    assert not np.array_equal(permuted, labels)
    for s in range(5):
        assert permuted[codes == s].sum() == labels[codes == s].sum()
