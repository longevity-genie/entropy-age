import gzip
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import link_membership, read_manifest_450k
from entropy_age.data.geo import read_supplementary_betas
from entropy_age.experiments.common import Context, Scale, dummies
from entropy_age.experiments.tissue_contrast import (
    ADJUSTMENTS,
    Prepared,
    _adjusted,
    _analyse,
    _boot_lunnon,
    _figure_cells,
    _run_parallel,
)
from entropy_age.paths import ProjectPaths
from entropy_age.stochasticity import dispersion_levels

MANIFEST_450K = """Illumina, Inc.
[Heading],,
[Assay],,
IlmnID,Name,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Accession,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island
cg01,cg01,1,15865,AARS;AARS,NM_001605;NM_001605,TSS200;Body,Island
cg02,cg02,2,534242,ZZZ;NEWHK,NR_1;NM_555555,Body;TSS1500,N_Shore
cg03,cg03,3,,,,,
cg04,cg04,X,710097,AARS,NM_001605,3'UTR,Island
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
    assert probes["pos"].to_list() == [15865, 534242, None, 710097]
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


def test_supplementary_betas_preamble_and_id_rows(tmp_path: Path):
    path = tmp_path / "lunnon.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "# row1: Sample barcodes\n# row2: New GSMs\n# row3: Old GSMs\n\n"
            ",bc1,bc2\n,GSM1,GSM2\n,GSMold1,GSMold2\ncg1,0.1,0.2\ncg2,0.3,0.4\n"
        )
    frame = read_supplementary_betas(path, preamble=4, after_header=2)
    assert frame.columns == ["cpg", "bc1", "bc2"]
    assert frame["cpg"].to_list() == ["cg1", "cg2"]
    assert np.allclose(frame["bc2"].to_numpy(), [0.2, 0.4])


def test_parallel_bootstrap_equals_serial():
    """Spawned workers get the state through the initializer and reproduce the serial results."""
    rng = np.random.default_rng(5)
    features = [f"cg{i}" for i in range(80)]
    donors = np.array([f"d{i}" for i in range(12)])
    tissues = {}
    for tissue in ("whole blood", "frontal cortex", "superior temporal gyrus", "entorhinal cortex", "cerebellum"):
        x = np.clip(rng.uniform(0.05, 0.95, 80) + rng.normal(0, 0.05, (12, 80)), 0, 1)
        cov = np.column_stack([rng.uniform(70, 90, 12), (rng.random(12) < 0.5).astype(float)])
        tissues[tissue] = Prepared(tissue, x, cov[:, 0], cov, features, donors, cov[:, 1:])
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(80)]})
    membership = pl.DataFrame({"feature": features[:20], "group": "housekeeping", "gene": "g"})
    state = {"lunnon": tissues, "probes": probes, "membership": membership}
    serial = pl.DataFrame(_run_parallel(_boot_lunnon, [1, 2, 3], 1, state)).sort("boot", "tissue")
    parallel = pl.DataFrame(_run_parallel(_boot_lunnon, [1, 2, 3], 2, state)).sort("boot", "tissue")
    assert serial.height == 15
    assert serial.equals(parallel)


def _synthetic_tissues(rng: np.random.Generator, n: int, features: list[str]) -> dict[str, Prepared]:
    donors = np.array([f"d{i}" for i in range(n)])
    out = {}
    for tissue in ("whole blood", "frontal cortex", "superior temporal gyrus", "entorhinal cortex", "cerebellum"):
        mix = rng.normal(size=n)
        x = rng.uniform(0.05, 0.95, len(features)) + rng.normal(0, 0.03, (n, len(features)))
        x += 0.03 * mix[:, None] * (rng.random(len(features)) < 0.3)
        sex = (rng.random(n) < 0.5).astype(float)[:, None]
        ages = rng.uniform(70, 90, n)
        out[tissue] = Prepared(tissue, np.clip(x, 0, 1), ages, np.column_stack([ages, sex]), features, donors, sex)
    return out


def test_composition_adjustment_uses_background_only_and_appends_components():
    rng = np.random.default_rng(8)
    features = [f"cg{i}" for i in range(300)]
    p = _synthetic_tissues(rng, 40, features)["whole blood"]
    members = set(features[:60])
    adjusted, info = _adjusted(p, members, 4, False)
    assert adjusted.covariates.shape == (40, p.covariates.shape[1] + 4)
    assert np.array_equal(adjusted.covariates[:, :2], p.covariates)
    assert info["background_cpgs"] == 240
    # the components ignore the gene-set CpGs: scrambling them changes nothing
    scrambled = Prepared(p.label, p.x.copy(), p.ages, p.covariates, p.features, p.donors, p.technical)
    scrambled.x[:, :60] = rng.permutation(scrambled.x[:, :60])
    again, _ = _adjusted(scrambled, members, 4, False)
    assert np.allclose(np.abs(again.covariates), np.abs(adjusted.covariates))
    age_free, info_free = _adjusted(p, members, 4, True)
    assert info_free["max_abs_r_age"] < 1e-8 < info["max_abs_r_age"]


def test_analyse_runs_with_every_adjustment(tmp_path: Path):
    rng = np.random.default_rng(9)
    features = [f"cg{i}" for i in range(200)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(200)]})
    membership = pl.DataFrame(
        {"feature": features[:60], "group": ["core_essential", "housekeeping", "non_essential"] * 20, "gene": "g"}
    )
    same = _synthetic_tissues(rng, 30, features)
    lifespan = {"blood": same["whole blood"], "cortex": same["frontal cortex"]}
    gaps = []
    for name, spec in ADJUSTMENTS.items():
        life, lun = lifespan, same
        if spec is not None:
            members = set(membership["feature"].to_list())
            life = {k: _adjusted(p, members, *spec)[0] for k, p in lifespan.items()}
            lun = {k: _adjusted(p, members, *spec)[0] for k, p in same.items()}
        point, boot_life, boot_same, contrasts = _analyse(life, lun, probes, membership, 3, 1)
        assert contrasts.height == 3 * (3 + 4)  # 3 sets x (3 lifespan metrics + 4 brain regions)
        assert set(boot_same["tissue"].unique()) == set(same)
        gaps.append(contrasts.with_columns(adjustment=pl.lit(name)))
    ctx = Context(ProjectPaths(tmp_path), Scale.of("quick"))
    assert _figure_cells(ctx, pl.concat(gaps)).exists()
