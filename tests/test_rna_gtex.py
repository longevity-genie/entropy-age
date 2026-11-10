import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from entropy_age.data.annotation import AnnotationLayout, gene_set_ensembl
from entropy_age.data.gtex_rna import ATTRIBUTES_ENV, COUNTS_ENV, SAMPLE_ATTRIBUTES, SUBJECT_PHENOTYPES
from entropy_age.experiments.common import Context, Scale
from entropy_age.experiments.rna_gtex import run
from entropy_age.paths import ProjectPaths
from entropy_age.stochasticity import expression_drift


def _hgnc(rows):
    return pl.DataFrame(
        rows,
        schema={
            "symbol": pl.String,
            "previous": pl.List(pl.String),
            "aliases": pl.List(pl.String),
            "entrez": pl.Int64,
            "ensembl": pl.String,
            "protein_coding": pl.Boolean,
        },
        orient="row",
    )


def test_gene_sets_reach_ensembl_by_entrez_then_symbol_then_unique_old_names():
    hgnc = _hgnc(
        [
            ("AARS1", ["AARS"], [], 16, "ENSG1", True),
            ("ACTB", [], ["PS1TP5BP1"], 60, "ENSG2", True),
            ("ONEX", [], ["SHARED"], 70, "ENSG3", True),
            ("TWOX", [], ["SHARED"], 80, "ENSG4", True),
            ("NOENS", [], [], 90, None, True),
        ]
    )
    sets = pl.DataFrame(
        {
            "gene_set": ["core_essential", "housekeeping", "housekeeping", "housekeeping", "housekeeping"],
            "symbol": ["WRONGNAME", "AARS", "PS1TP5BP1", "SHARED", "NOENS"],
            "entrez": [60, None, None, None, None],
            "refseq": [None] * 5,
        },
        schema_overrides={"entrez": pl.Int64, "refseq": pl.String},
    )
    out = gene_set_ensembl(sets, hgnc)
    got = {(r["group"], r["feature"], r["gene"], r["matched_by"]) for r in out.iter_rows(named=True)}
    assert got == {
        ("core_essential", "ENSG2", "ACTB", "entrez"),  # the Entrez ID wins over a wrong symbol
        ("housekeeping", "ENSG1", "AARS1", "previous"),  # renamed in 2019
        ("housekeeping", "ENSG2", "ACTB", "alias"),
    }  # SHARED names two genes and NOENS has no Ensembl ID: both left out


def test_expression_drift_finds_the_genes_whose_spread_grows():
    rng = np.random.default_rng(0)
    ages = np.repeat(np.arange(25.0, 80.0, 10.0), 60)
    sd = np.where(np.arange(40) < 20, 0.2 + 0.6 * (ages[:, None] - 25) / 50, 0.5)
    x = 5 + rng.normal(size=(ages.size, 40)) * sd + 0.01 * ages[:, None]
    out = expression_drift(x, ages, covariates=rng.normal(size=(ages.size, 2)))
    spreading, steady = out["dispersion_slope"].to_numpy()[:20], out["dispersion_slope"].to_numpy()[20:]
    assert spreading.min() > steady.max()
    assert np.allclose(out["slope"].to_numpy(), 0.01, atol=0.01)
    assert np.allclose(out["baseline"].to_numpy(), x.mean(axis=0))
    with pytest.raises(ValueError, match="one age per sample"):
        expression_drift(x, ages[:-1])


TISSUES = {"Whole Blood": True, "Liver": True, "Muscle - Skeletal": True, "Brain - Cortex": False}
GENES = 360
CORE, HOUSEKEEPING, NON_ESSENTIAL = range(0, 50), range(50, 110), range(110, 140)


def _toy_gtex(root: Path) -> tuple[Path, Path, ProjectPaths]:
    """Four 100-donor tissues and a small one; set genes resist growing spread except in brain."""
    rng = np.random.default_rng(1)
    genes = [f"ENSG{j:05d}" for j in range(GENES)]
    hgnc = _hgnc(
        [(f"G{j}", [], [], 1000 + j, genes[j], True) for j in range(GENES - 20)]
        + [(f"NC{j}", [], [], 5000 + j, genes[j], False) for j in range(GENES - 20, GENES)]
    )
    paths = ProjectPaths(root / "project").ensure()
    annotation = AnnotationLayout(paths)
    annotation.interim.mkdir(parents=True, exist_ok=True)
    hgnc.write_parquet(annotation.hgnc)
    sets = [("core_essential", j) for j in CORE] + [("housekeeping", j) for j in HOUSEKEEPING]
    sets += [("non_essential", j) for j in NON_ESSENTIAL]
    pl.DataFrame(
        {
            "gene_set": [s for s, _ in sets],
            "symbol": [f"G{j}" for _, j in sets],
            "entrez": [1000 + j for _, j in sets],
            "refseq": [None] * len(sets),
        },
        schema_overrides={"refseq": pl.String},
    ).write_parquet(annotation.gene_sets)

    counts_dir = root / "gtex" / "counts"
    counts_dir.mkdir(parents=True)
    level = rng.uniform(np.log(20), np.log(2000), GENES)
    members = np.zeros(GENES, dtype=bool)
    members[: NON_ESSENTIAL.stop] = True
    samples, attributes, columns = [], [], {}
    donors = [f"GTEX-{d:04d}" for d in range(110)]
    bins = rng.integers(2, 8, len(donors)) * 10
    for tissue, protected in {**TISSUES, "Bladder": True}.items():
        n = 100 if tissue != "Bladder" else 12
        for d in range(n):
            sample = f"{donors[d]}-{tissue[:3]}"
            age = bins[d] + 4.5
            gain = np.where(members & protected, 0.0, 0.5) * (age - 25) / 50
            mu = np.exp(level + rng.normal(size=GENES) * (0.15 + gain)) * 2.0
            columns[sample] = rng.poisson(mu).astype(np.uint32)
            samples.append((sample, donors[d], tissue, ["male", "female"][d % 2], float(bins[d]), bins[d] + 9.0))
            attributes.append(f"{sample}\t{rng.uniform(6, 9):.1f}\t{rng.uniform(100, 900):.0f}")
    attributes.append(f"{samples[0][0]}-extra\tNA\t")  # a sample the counts do not have
    pl.DataFrame(
        samples, schema=["sample_id", "donor_id", "tissue", "sex", "age_min", "age_max"], orient="row"
    ).with_columns(exclude_reason=pl.lit(None, dtype=pl.String)).write_parquet(root / "gtex" / "samples.parquet")
    pl.DataFrame({"gene_id": genes}).write_parquet(root / "gtex" / "genes.parquet")
    names = list(columns)
    index = {}
    for part, chunk in enumerate((names[:250], names[250:])):
        pl.DataFrame({"gene_id": genes, **{s: columns[s] for s in chunk}}).write_parquet(
            counts_dir / f"part-{part:05d}.parquet"
        )
        index.update(dict.fromkeys(chunk, f"part-{part:05d}.parquet"))
    (counts_dir / "index.json").write_text(json.dumps(index))
    attr = root / "attributes"
    attr.mkdir()
    (attr / SAMPLE_ATTRIBUTES).write_text("SAMPID\tSMRIN\tSMTSISCH\n" + "\n".join(attributes) + "\n")
    hardy = "\n".join(f"{d}\t{i % 5 if i % 7 else ''}" for i, d in enumerate(donors))
    (attr / SUBJECT_PHENOTYPES).write_text("SUBJID\tDTHHRDY\n" + hardy + "\n")
    return root / "gtex", attr, paths


def test_rna_gtex_refuses_without_its_settings(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv(COUNTS_ENV, raising=False)
    monkeypatch.setenv(ATTRIBUTES_ENV, str(tmp_path))
    ctx = Context(ProjectPaths(tmp_path).ensure(), Scale.of("quick"))
    assert run(ctx) == []
    assert f"{COUNTS_ENV} not set" in capsys.readouterr().err
    monkeypatch.setenv(COUNTS_ENV, str(tmp_path / "nowhere"))
    assert run(ctx) == []
    assert "GTEx RNA inputs missing" in capsys.readouterr().err


def test_rna_gtex_finds_the_missing_protection_in_brain(tmp_path, monkeypatch):
    counts, attributes, paths = _toy_gtex(tmp_path)
    monkeypatch.setenv(COUNTS_ENV, str(counts))
    monkeypatch.setenv(ATTRIBUTES_ENV, str(attributes))
    monkeypatch.setenv("ENTROPY_AGE_BOOT", "12")
    monkeypatch.setenv("ENTROPY_AGE_WORKERS", "1")
    figures = run(Context(paths, Scale.of("quick")))
    assert figures == [paths.figures / "rna_gtex.png"] and figures[0].exists()

    out = paths.output / "rna_gtex"
    tissues = pl.read_csv(out / "tissues.csv")
    used = dict(tissues.select("tissue", "used").iter_rows())
    assert used == {**dict.fromkeys(TISSUES, True), "Bladder": False}
    # non-coding genes never enter, so at most the protein-coding ones are kept
    assert tissues.filter(pl.col("used"))["genes"].max() <= GENES - 20
    assert tissues.filter(pl.col("used"))["core_essential"].to_list() == [len(CORE)] * len(TISSUES)

    contrasts = pl.read_csv(out / "contrasts.csv")
    pick = contrasts.filter((pl.col("metric") == "dispersion_slope") & (pl.col("group") == "core_essential")).sort(
        "adjustment", "contrast"
    )
    by = {(r["adjustment"], r["contrast"]): r for r in pick.iter_rows(named=True)}
    for adjustment in ("covariates", "pc5_age_free"):
        for faster in ("renewing", "slow", "heart & muscle"):
            row = by[(adjustment, f"{faster} - brain")]
            assert row["ci_high"] < 0, row  # protected tissues minus the unprotected brain
        assert abs(by[(adjustment, "renewing - slow")]["difference"]) < abs(
            by[(adjustment, "renewing - brain")]["difference"]
        )

    point = pl.read_csv(out / "point_estimates.csv")
    core = point.filter(
        (pl.col("metric") == "dispersion_slope")
        & (pl.col("group") == "core_essential")
        & (pl.col("adjustment") == "covariates")
    )
    values = dict(core.select("tissue", "matched_diff_sd").iter_rows())
    assert values["Brain - Cortex"] > -0.3 and max(v for t, v in values.items() if TISSUES[t]) < -0.5
    intervals = pl.read_csv(out / "tissue_intervals.csv")
    assert set(intervals["resamples"].to_list()) == {12}
