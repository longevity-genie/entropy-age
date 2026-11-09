import gzip
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.cohorts import _huang
from entropy_age.data.mouse import mouse_gene_sets, mouse_membership, read_orthologs, read_platform_table
from entropy_age.experiments.mouse_tissues import TURNOVER, MouseTissue, _boot, summarise, tissue_measures
from entropy_age.experiments.tissue_contrast import _matched_sd, _run_parallel

PLATFORM = (
    "^PLATFORM = GPL32543\n!Platform_title = HorvathMammalMethylChip320\n#ID = \n#SYMBOL = \n"
    "!platform_table_begin\n"
    "ID\tSPOT_ID\tSYMBOL\tmm10.Chr\tannotation\tPromoter\tCpG.island\tENTREZID\n"
    "cg01_BC21\tcg01\tGapdh\t6\tPromoter (<=1kb)\tTRUE\tyes\t14433\n"
    'cg02_BC21\tcg02\tZfr\t15\t"Exon (Zfr-201/Zfr\t exon 16 of 20)"\tFALSE\tno\t22763\n'
    "cg03_TC21\tcg03\tNA\tNA\tNA\tNA\tNA\tNA\n"
    "cg04_BC11\tcg04\tRbmx\tX\t3' UTR\tTRUE\tNA\t19655\n"
    "!platform_table_end\n"
)

ORTHOLOGS = (
    "DB Class Key\tCommon Organism Name\tNCBI Taxon ID\tSymbol\tEntrezGene ID\tMouse MGI ID\n"
    "1\tmouse, laboratory\t10090\tGapdh\t14433\tMGI:1\n"
    "1\thuman\t9606\tGAPDH\t2597\t\n"
    "2\tmouse, laboratory\t10090\tZfr\t22763\tMGI:2\n"
    "2\tmouse, laboratory\t10090\tZfr2\t99999\tMGI:3\n"
    "2\thuman\t9606\tZFR\t51663\t\n"
    "3\thuman\t9606\tHUMANONLY\t777\t\n"
    "4\trat\t10116\tGapdh\t24383\t\n"
)


def test_platform_table_keeps_mapped_probes_and_flags(tmp_path: Path):
    path = tmp_path / "GPL32543_table.txt"
    path.write_text(PLATFORM)
    probes = read_platform_table(path)
    assert probes["cpg"].to_list() == ["cg01_BC21", "cg02_BC21", "cg04_BC11"]  # cg03 has no mm10 position
    assert probes["promoter"].to_list() == [True, False, True]
    assert probes["cpg_island"].to_list() == [True, False, False]  # NA is not an island
    assert probes["entrez"].to_list() == [14433, 22763, 19655]
    assert probes["chr"].to_list() == ["6", "15", "X"]


def test_gene_sets_carry_to_every_mouse_ortholog(tmp_path: Path):
    path = tmp_path / "HOM.rpt"
    path.write_text(ORTHOLOGS)
    orthologs = read_orthologs(path)
    assert set(orthologs["organism"].to_list()) == {"human", "mouse"}  # rat rows dropped
    human = pl.DataFrame(
        {
            "gene_set": ["housekeeping", "core_essential", "core_essential", "non_essential"],
            "symbol": ["GAPDH", "ZFR", "WRONGSYM", "HUMANONLY"],
            "entrez": [None, None, 51663, 777],
        },
        schema={"gene_set": pl.String, "symbol": pl.String, "entrez": pl.Int64},
    )
    sets = mouse_gene_sets(human, orthologs)
    pairs = set(zip(sets["gene_set"].to_list(), sets["symbol"].to_list(), strict=True))
    # GAPDH by symbol; ZFR by symbol and (as WRONGSYM) by Entrez -> both mouse paralogs; HUMANONLY has none
    assert pairs == {("housekeeping", "Gapdh"), ("core_essential", "Zfr"), ("core_essential", "Zfr2")}

    platform = tmp_path / "GPL32543_table.txt"
    platform.write_text(PLATFORM)
    membership = mouse_membership(read_platform_table(platform), sets)
    # only promoter probes count: cg02 (Zfr) is an exon probe
    assert set(zip(membership["feature"].to_list(), membership["group"].to_list(), strict=True)) == {
        ("cg01_BC21", "housekeeping")
    }


def test_huang_sheet_maps_groups_and_columns(tmp_path: Path):
    meta = tmp_path / "meta.csv.gz"
    with gzip.open(meta, "wt") as handle:
        handle.write(
            "sname,tissue,treatment,animal\n"
            "208_R01C01,WBC,Young_ctrl,A1\n208_R02C01,FCTX,Aged_rapamycin,A2\n209_R01C01,Colon,Aged_ctrl_KD,A3\n"
        )
    betas = tmp_path / "betas.csv.gz"
    with gzip.open(betas, "wt") as handle:
        handle.write(
            "Probe,209_R01C01,208_R01C01,208_R02C01,extra\ncg1_BC21,0.1,0.2,0.3,0.9\ncg2_BC21,0.4,0.5,0.6,0.9\n"
        )
    samples, frame = _huang(meta, betas)
    assert samples["tissue"].to_list() == ["white blood cells", "frontal cortex", "colon"]
    assert samples["age_group"].to_list() == ["young", "aged", "aged"]
    assert samples["treatment"].to_list() == ["ctrl", "rapamycin", "ctrl"]
    assert samples["knockdown"].to_list() == [False, False, True]
    assert frame.columns == ["cpg", "208_R01C01", "208_R02C01", "209_R01C01"]  # sheet order, extra dropped
    assert np.allclose(frame["208_R02C01"].to_numpy(), [0.3, 0.6])


def _synthetic(rng: np.random.Generator, features: list[str], members: int) -> dict[str, MouseTissue]:
    """Aged spread grows everywhere; set members are protected in proportion to tissue turnover."""
    protect = {"renewing": 0.3, "slow": 0.6, "post-mitotic": 1.0}
    groups = np.array(["Young_ctrl"] * 8 + ["Aged_ctrl"] * 8 + ["Aged_exenatide"] * 6 + ["Aged_rapamycin"] * 5)
    donors = np.array([f"m{i}" for i in range(groups.size)])
    out = {}
    for tissue, turnover in TURNOVER.items():
        level = rng.uniform(0.1, 0.9, len(features))
        sd = np.where(groups[:, None] == "Young_ctrl", 0.02, 0.06) * np.ones((1, len(features)))
        sd[:, :members] = np.where(groups[:, None] == "Young_ctrl", 0.02, 0.02 + 0.04 * protect[turnover])
        x = np.clip(level + sd * rng.normal(size=sd.shape), 0, 1)
        out[tissue] = MouseTissue(tissue, x, groups, donors, features)
    return out


def test_turnover_gradient_is_recovered_and_bootstrap_is_reproducible():
    rng = np.random.default_rng(11)
    features = [f"cg{i}" for i in range(400)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(400)]})
    membership = pl.DataFrame({"feature": features[:80], "group": "housekeeping", "gene": "g"})
    tissues = _synthetic(rng, features, 80)

    measures = tissue_measures(tissues["colon"], np.arange(27), probes)
    assert set(measures) == {
        "aging",
        "aged dispersion",
        "dispersion Aged_ctrl",
        "dispersion Aged_exenatide",
        "dispersion Aged_rapamycin",
    }

    def aged_spread(tissue: str) -> float:
        metrics, names = tissue_measures(tissues[tissue], np.arange(27), probes)["aged dispersion"]
        return _matched_sd(metrics, membership, names, 0, 0)["matched_diff_sd"][0]

    assert aged_spread("colon") < aged_spread("liver") < aged_spread("heart") < 0.15

    animals = {
        g: np.unique(tissues["colon"].donors[tissues["colon"].groups == g]) for g in np.unique(tissues["colon"].groups)
    }
    state = {"tissues": tissues, "probes": probes, "membership": membership, "animals": animals}
    serial = pl.DataFrame(_run_parallel(_boot, list(range(12)), 1, state))
    parallel = pl.DataFrame(_run_parallel(_boot, list(range(12)), 2, state))
    keys = ["boot", "tissue", "measure", "metric", "group"]
    assert serial.sort(keys).equals(parallel.sort(keys))
    spread = serial.group_by(keys[1:]).agg(pl.col("matched_diff_sd").std())["matched_diff_sd"]
    assert spread.min() > 0  # every resample draws other animals

    summary = summarise(serial.head(0), serial)
    gradient = summary.filter(
        (pl.col("contrast") == "renewing - post-mitotic") & (pl.col("measure") == "aged dispersion")
    )
    assert gradient.height == 1 and gradient["ci_high"][0] < 0  # more protection where cells renew
    assert summary.filter(pl.col("contrast").str.contains("Aged_rapamycin - Aged_ctrl")).height == 3
