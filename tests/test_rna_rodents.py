import gzip
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import AnnotationLayout, symbol_ensembl
from entropy_age.experiments.common import Context, Scale
from entropy_age.experiments.rna_rodents import SPECIES, read_atlas, run
from entropy_age.paths import ProjectPaths


def _hgnc(rows):
    schema = {
        "symbol": pl.String,
        "previous": pl.List(pl.String),
        "aliases": pl.List(pl.String),
        "entrez": pl.Int64,
        "ensembl": pl.String,
        "protein_coding": pl.Boolean,
    }
    return pl.DataFrame(rows, schema=schema, orient="row")


def test_rodent_symbols_reach_human_genes_by_upper_case_and_unique_old_names():
    hgnc = _hgnc(
        [
            ("ACTB", [], [], 60, "E1", True),
            ("AARS1", ["AARS"], [], 16, "E2", True),
            ("ONEX", [], ["SHARED"], 70, "E3", True),
            ("TWOX", [], ["SHARED"], 80, "E4", True),
            ("GAPDH", [], ["G3PD"], 2597, "E5", True),
        ]
    )
    out = symbol_ensembl(["Actb", "Aars", "Shared", "Gapdh", "G3pd", "Nothing"], hgnc)
    got = {(r["symbol"], r["ensembl"], r["matched_by"]) for r in out.iter_rows(named=True)}
    # Gapdh and G3pd both reach E5, so neither stands for it; Shared names two genes; Nothing is unknown
    assert got == {("Actb", "E1", "symbol"), ("Aars", "E2", "previous")}


TISSUES = {"spleen": True, "liver": True, "left ventricle": True, "cortex": False}
GENES = 320
SET_GENES = range(0, 120)


def _gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _toy_atlas(root: Path, rng: np.random.Generator) -> ProjectPaths:
    paths = ProjectPaths(root).ensure()
    annotation = AnnotationLayout(paths)
    annotation.interim.mkdir(parents=True, exist_ok=True)
    _hgnc([(f"GENE{j}", [], [], 1000 + j, f"ENSG{j:05d}", j < GENES - 20) for j in range(GENES)]).write_parquet(
        annotation.hgnc
    )
    groups = ["core_essential"] * 40 + ["housekeeping"] * 50 + ["non_essential"] * 30
    pl.DataFrame(
        {"gene_set": groups, "symbol": [f"GENE{j}" for j in SET_GENES], "entrez": [1000 + j for j in SET_GENES]},
    ).with_columns(refseq=pl.lit(None, dtype=pl.String)).write_parquet(annotation.gene_sets)
    level = rng.uniform(np.log(20), np.log(2000), GENES)
    for species, (gse, counts_name) in SPECIES.items():
        folder = paths.input / "geo" / gse
        folder.mkdir(parents=True)
        titles, ids, tissues, sexes, ages, columns = [], [], [], [], [], {}
        n_tissues = {**TISSUES, "lens": True} if species == "rat" else TISSUES
        for tissue, protected in n_tissues.items():
            for a in range(40):
                sid = f"{len(ids) + 100000}"
                age = float([6, 12, 18, 21, 24, 27][a % 6])
                sex = ["male", "female"][a % 2]
                gain = np.where((np.arange(GENES) < SET_GENES.stop) & protected, 0.0, 0.6) * (age - 6) / 21
                columns[sid] = rng.poisson(np.exp(level + rng.normal(size=GENES) * (0.15 + gain)) * 2.0)
                ids.append(sid)
                titles.append(f"{sex}_{tissue}_{age:g} months_{sid}")
                tissues.append(tissue)
                sexes.append(sex)
                ages.append(f"{age:g} months")
        q = lambda v: "\t".join(f'"{x}"' for x in v)  # noqa: E731
        _gz(
            folder / f"{gse}_series_matrix.txt.gz",
            f"!Sample_title\t{q(titles)}\n!Sample_geo_accession\t{q(f'GSM{i}' for i in ids)}\n"
            f"!Sample_characteristics_ch1\t{q('tissue: ' + t for t in tissues)}\n"
            f"!Sample_characteristics_ch1\t{q('Sex: ' + s for s in sexes)}\n"
            f"!Sample_characteristics_ch1\t{q('age: ' + a for a in ages)}\n"
            '!series_matrix_table_begin\n"ID_REF"\n!series_matrix_table_end\n',
        )
        header = "GeneName\t" + "\t".join(ids[1:])  # the first GEO sample has no count column
        body = "\n".join(f"Gene{j}\t" + "\t".join(str(columns[s][j]) for s in ids[1:]) for j in range(GENES))
        _gz(folder / counts_name, header + "\n" + body + "\n")
    return paths


def test_read_atlas_maps_titles_to_count_columns(tmp_path, capsys):
    paths = _toy_atlas(tmp_path, np.random.default_rng(0))
    gse, counts = SPECIES["mouse"]
    samples, genes, x = read_atlas(
        paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz", paths.input / "geo" / gse / counts
    )
    assert "1 GEO samples have no count column" in capsys.readouterr().err
    assert samples.height == x.shape[0] == 4 * 40 - 1 and x.shape[1] == GENES == len(genes)
    assert samples["age"].min() == 6.0 and samples["age"].max() == 27.0
    assert samples["sample"][0] == "100001"


def test_rna_rodents_finds_the_missing_protection_in_brain(tmp_path, monkeypatch):
    paths = _toy_atlas(tmp_path, np.random.default_rng(1))
    monkeypatch.setenv("ENTROPY_AGE_BOOT", "12")
    monkeypatch.setenv("ENTROPY_AGE_WORKERS", "1")
    assert run(Context(paths, Scale.of("quick"))) == [paths.figures / "rna_rodents.png"]
    out = paths.output / "rna_rodents"
    tissues = pl.read_csv(out / "tissues.csv")
    lens = tissues.filter(pl.col("tissue") == "lens").row(0, named=True)
    assert lens["used"] is False and lens["species"] == "rat"
    assert tissues.filter(pl.col("used"))["genes"].max() <= GENES - 20  # non-coding genes never enter
    genes = pl.read_csv(out / "genes.csv")
    assert genes["protein_coding"].to_list() == [GENES - 20] * 2 and genes["subsample"].to_list() == [27, 27]
    contrasts = pl.read_csv(out / "contrasts.csv").filter(
        (pl.col("metric") == "dispersion_slope") & (pl.col("group") == "core_essential")
    )
    for species in SPECIES:
        for adjustment in ("covariates", "pc5_age_free"):
            for faster in ("renewing", "slow", "heart & muscle"):
                row = contrasts.filter(
                    (pl.col("species") == species)
                    & (pl.col("adjustment") == adjustment)
                    & (pl.col("contrast") == f"{faster} - brain")
                ).row(0, named=True)
                assert row["ci_high"] < 0, row
    assert set(pl.read_csv(out / "tissue_intervals.csv")["resamples"].to_list()) == {12}
