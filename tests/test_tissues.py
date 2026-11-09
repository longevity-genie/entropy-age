import gzip
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import link_membership, read_manifest_450k
from entropy_age.data.geo import read_supplementary_betas
from entropy_age.experiments.common import dummies
from entropy_age.stochasticity import dispersion_levels

MANIFEST_450K = """Illumina, Inc.
[Heading],,
[Assay],,
IlmnID,Name,CHR,UCSC_RefGene_Name,UCSC_RefGene_Accession,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island
cg01,cg01,1,AARS;AARS,NM_001605;NM_001605,TSS200;Body,Island
cg02,cg02,2,ZZZ;NEWHK,NR_1;NM_555555,Body;TSS1500,N_Shore
cg03,cg03,3,,,,
cg04,cg04,X,AARS,NM_001605,3'UTR,Island
[Controls],,
ctl,1
"""


def test_450k_manifest_links_and_promoters(tmp_path: Path):
    path = tmp_path / "m450.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(MANIFEST_450K)
    probes, links = read_manifest_450k(path)
    assert probes["cpg"].to_list() == ["cg01", "cg02", "cg03", "cg04"]
    assert probes["cpg_island"].to_list() == [True, False, False, True]
    assert probes["island_relation"].to_list()[2] == "OpenSea"
    # promoter if ANY of the CpG's genes has it in TSS1500/TSS200/5'UTR/1stExon
    assert probes["promoter"].to_list() == [True, True, False, False]
    pairs = set(zip(links["cpg"].to_list(), links["symbol"].to_list(), links["region"].to_list(), strict=True))
    assert pairs == {
        ("cg01", "AARS", "TSS200"),
        ("cg01", "AARS", "Body"),
        ("cg02", "ZZZ", "Body"),
        ("cg02", "NEWHK", "TSS1500"),
        ("cg04", "AARS", "3'UTR"),
    }

    sets = pl.DataFrame(
        {
            "gene_set": ["core_essential", "housekeeping", "housekeeping"],
            "symbol": ["AARS", "RENAMED", "ZZZ"],
            "entrez": [16, None, None],
            "refseq": [None, "NM_555555", None],
        }
    )
    member = link_membership(links, sets)
    got = set(zip(member["feature"].to_list(), member["group"].to_list(), strict=True))
    # cg01 by symbol in a promoter region; cg02 by RefSeq (its promoter gene); ZZZ is gene-body only;
    # cg04 is AARS but 3'UTR, so not a promoter link
    assert got == {("cg01", "core_essential"), ("cg02", "housekeeping")}
    everything = link_membership(links, sets, promoter_only=False)
    assert ("cg04", "core_essential") in set(
        zip(everything["feature"].to_list(), everything["group"].to_list(), strict=True)
    )


def test_supplementary_betas_suffix(tmp_path: Path):
    path = tmp_path / "proc.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write('"","S1_Beta","S1_detectP","S2_Beta","S2_detectP"\n"cg1",0.1,0.01,0.2,0.02\n"cg2",0.9,0,NA,0\n')
    frame = read_supplementary_betas(path, suffix="_Beta")
    assert frame.columns == ["cpg", "S1", "S2"]
    assert frame["cpg"].to_list() == ["cg1", "cg2"]
    assert frame["S2"][1] is None and abs(frame["S1"][1] - 0.9) < 1e-6


def test_dispersion_levels_adjust_and_scale():
    rng = np.random.default_rng(4)
    n = 400
    covariate = rng.normal(size=n)
    noise_sd = np.array([0.01, 0.05])
    x = np.column_stack([0.5 + 0.1 * covariate, 0.1 + 0.0 * covariate]) + rng.normal(size=(n, 2)) * noise_sd
    raw = dispersion_levels(x)
    adjusted = dispersion_levels(x, covariates=covariate[:, None])
    # adjusting for the covariate removes its spread from the first feature only
    assert raw["residual_sd"][0] > 5 * adjusted["residual_sd"][0]
    assert np.allclose(adjusted["residual_sd"].to_numpy(), noise_sd, rtol=0.1)
    m = x.mean(axis=0)
    assert np.allclose(adjusted["dispersion"].to_numpy(), adjusted["residual_sd"].to_numpy() / np.sqrt(m * (1 - m)))


def test_dummies_single_level_is_empty():
    assert dummies(np.array(["a", "a"])).shape == (2, 0)
