"""Mouse blood RRBS of Petkovich et al. 2017 (GSE80672), on promoter CpGs of GRCm38 (roadmap RM3).

GEO ships one table per sample, ``<GSM>_<title>.overlap.txt.gz``: a CpG per row as
``gi|...|ref|NC_0000xx.v|:pos``, its percentage methylated and its read coverage. The contigs are GRCm38
(mm10) RefSeq chromosomes. Only CpGs in promoters (1,500 bp upstream to 500 bp downstream of a RefSeq
transcription start, the 450k array's TSS1500 and TSS200) are kept, at a coverage of at least ``MIN_COVERAGE``,
and only those covered in at least ``MIN_PRESENT`` of the blood samples enter the matrix.

Promoters, Entrez IDs and CpG islands come from UCSC's mm10 tables (NCBI RefSeq genes and their links, and
the CpG-island track), so gene sets are matched by Entrez ID, as the mouse array's are.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.download import download
from entropy_age.data.geo import _split
from entropy_age.paths import ProjectPaths

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/{name}"
TRACKS = ("ncbiRefSeq.txt.gz", "ncbiRefSeqLink.txt.gz", "cpgIslandExt.txt.gz")
UPSTREAM, DOWNSTREAM = 1500, 500
MIN_COVERAGE = 10
MIN_PRESENT = 0.9
#: GRCm38 RefSeq chromosome accessions: NC_000067 is chromosome 1 ... NC_000085 is 19, then X and Y
_CHROMOSOMES = {f"NC_{67 + k:06d}": str(k + 1) for k in range(19)} | {"NC_000086": "X", "NC_000087": "Y"}
_CONTIG = re.compile(r"ref\|(NC_\d{6})\.\d+\|:(\d+)$")


@dataclass(frozen=True, slots=True)
class RrbsLayout:
    paths: ProjectPaths

    @property
    def input(self) -> Path:
        return self.paths.input / "geo" / "GSE80672"

    @property
    def matrix(self) -> Path:
        return self.input / "GSE80672_series_matrix.txt.gz"

    def track(self, name: str) -> Path:
        return self.paths.input / "annotation" / "mm10" / name

    @property
    def samples(self) -> Path:
        return self.paths.interim / "cohorts" / "petkovich_samples.parquet"

    @property
    def betas(self) -> Path:
        return self.paths.interim / "cohorts" / "petkovich_betas.parquet"

    @property
    def cpgs(self) -> Path:
        return self.paths.interim / "cohorts" / "petkovich_cpgs.parquet"


def read_samples(path: Path) -> pl.DataFrame:
    """One row per sample: accession, title, tissue, strain, age in months, sex, condition, diet, file URL.

    GEO's field says ``age (years)``, but the values are months (Petkovich et al. 2017 profile mice of
    0.67-35 months); values such as ``2.5 (at isolation)`` keep their leading number.
    """
    meta: dict[str, list[str]] = {}
    keyed: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1"):
                cells = _split(line)[1:]
                key = cells[0].split(":", 1)[0].strip()
                keyed[key] = [c.split(":", 1)[1].strip() if ":" in c else "" for c in cells]
            elif line.startswith(("!Sample_geo_accession", "!Sample_title", "!Sample_supplementary_file_1")):
                cells = _split(line)
                meta[cells[0].removeprefix("!Sample_")] = cells[1:]
    needed = ("tissue", "strain", "age (years)", "gender", "genetic condition", "diet")
    missing = [k for k in needed if k not in keyed]
    if missing:
        raise ValueError(f"{path.name}: no characteristics {missing}")
    age = [re.match(r"[\d.]+", a) for a in keyed["age (years)"]]
    return pl.DataFrame(
        {
            "sample": meta["geo_accession"],
            "title": meta["title"],
            "tissue": keyed["tissue"],
            "strain": keyed["strain"],
            "age_months": [float(a.group()) if a else None for a in age],
            "sex": keyed["gender"],
            "condition": keyed["genetic condition"],
            "diet": keyed["diet"],
            "url": [u.replace("ftp://", "https://") for u in meta["supplementary_file_1"]],
        }
    )


def fetch_rrbs(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    """The series matrix, the mm10 tracks and every whole-blood sample's table (about 3 GB)."""
    lay = RrbsLayout(paths)
    gse = "GSE80672"
    files = [
        download(
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/{gse}/matrix/{gse}_series_matrix.txt.gz", lay.matrix
        )
    ]
    files += [download(UCSC.format(name=n), lay.track(n), overwrite=overwrite, timeout=600) for n in TRACKS]
    blood = read_samples(lay.matrix).filter(pl.col("tissue") == "whole blood")
    for url in blood["url"]:
        files.append(download(url, lay.input / url.rsplit("/", 1)[1], overwrite=overwrite, timeout=1800))
    return files


def _columns(path: Path, index: list[int], names: list[str], **kwargs) -> pl.DataFrame:
    """Columns ``index`` of a headerless UCSC table, named ``names``, every column read as text first."""
    frame = pl.read_csv(path, separator="\t", has_header=False, infer_schema=False, quote_char=None, **kwargs)
    return frame.select([frame.columns[i] for i in index]).rename(
        dict(zip([frame.columns[i] for i in index], names, strict=True))
    )


def promoters(refseq: Path, link: Path) -> pl.DataFrame:
    """``chr``, ``start``, ``end`` (0-based, half-open) and ``entrez`` of every RefSeq transcript's promoter."""
    genes = _columns(refseq, [1, 2, 3, 4, 5], ["name", "chrom", "strand", "tx_start", "tx_end"]).with_columns(
        pl.col("tx_start").cast(pl.Int64), pl.col("tx_end").cast(pl.Int64)
    )
    links = _columns(link, [0, 6], ["name", "entrez"], truncate_ragged_lines=True).with_columns(
        pl.col("entrez").cast(pl.Int64, strict=False)
    )
    tss = pl.when(pl.col("strand") == "+").then(pl.col("tx_start")).otherwise(pl.col("tx_end") - 1)
    up = pl.when(pl.col("strand") == "+").then(UPSTREAM).otherwise(DOWNSTREAM)
    down = pl.when(pl.col("strand") == "+").then(DOWNSTREAM).otherwise(UPSTREAM)
    return (
        genes.join(links.drop_nulls("entrez"), on="name", how="inner")
        .with_columns(chr=pl.col("chrom").str.strip_prefix("chr"), tss=tss)
        .filter(pl.col("chr").is_in(list(_CHROMOSOMES.values())))
        .select("chr", (pl.col("tss") - up).clip(0).alias("start"), (pl.col("tss") + down + 1).alias("end"), "entrez")
        .unique()
    )


def islands(path: Path) -> pl.DataFrame:
    frame = _columns(path, [1, 2, 3], ["chrom", "start", "end"])
    return frame.select(
        pl.col("chrom").str.strip_prefix("chr").alias("chr"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )


def in_intervals(chroms: np.ndarray, positions: np.ndarray, intervals: pl.DataFrame) -> np.ndarray:
    """Whether each 0-based position falls in any interval of its chromosome (intervals may overlap)."""
    out = np.zeros(positions.size, dtype=bool)
    for (chrom,), sel in intervals.group_by("chr"):
        rows = np.flatnonzero(chroms == chrom)
        if rows.size == 0:
            continue
        merged = _merge(sel["start"].to_numpy(), sel["end"].to_numpy())
        i = np.searchsorted(merged[:, 0], positions[rows], side="right") - 1
        ok = i >= 0
        out[rows[ok]] = positions[rows[ok]] < merged[i[ok], 1]
    return out


def _merge(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    order = np.argsort(starts)
    s, e = starts[order], ends[order]
    out = []
    for a, b in zip(s, e, strict=True):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def read_overlap(path: Path, keep: pl.DataFrame | None = None) -> pl.DataFrame:
    """``cpg`` (``chr:pos``, 1-based), ``beta`` and ``coverage`` of one sample's table, chromosomes only.

    With ``keep`` (a frame of ``cpg``), only those CpGs are returned.
    """
    frame = pl.read_csv(path, separator="\t", has_header=True, new_columns=["index", "percent", "coverage"])
    parts = frame["index"].str.extract_groups(_CONTIG.pattern)
    frame = frame.with_columns(acc=parts.struct.field("1"), pos=parts.struct.field("2").cast(pl.Int64))
    frame = frame.filter(pl.col("acc").is_in(list(_CHROMOSOMES)))
    out = frame.select(
        pl.format("{}:{}", pl.col("acc").replace_strict(_CHROMOSOMES), pl.col("pos")).alias("cpg"),
        (pl.col("percent") / 100).cast(pl.Float32).alias("beta"),
        pl.col("coverage").cast(pl.Int32),
    )
    return out if keep is None else out.join(keep, on="cpg", how="semi")


def prepare_rrbs(paths: ProjectPaths) -> RrbsLayout:
    """Samples, the promoter CpG table (with island flag and Entrez IDs) and the wide beta matrix of blood."""
    lay = RrbsLayout(paths)
    samples = read_samples(lay.matrix).filter(pl.col("tissue") == "whole blood")
    prom = promoters(lay.track(TRACKS[0]), lay.track(TRACKS[1]))
    isl = islands(lay.track(TRACKS[2]))
    long = []
    for row in samples.iter_rows(named=True):
        table = read_overlap(lay.input / row["url"].rsplit("/", 1)[1]).filter(pl.col("coverage") >= MIN_COVERAGE)
        chrom = table["cpg"].str.split(":").list.get(0).to_numpy()
        pos = table["cpg"].str.split(":").list.get(1).cast(pl.Int64).to_numpy() - 1
        table = table.filter(pl.Series(in_intervals(chrom, pos, prom)))
        long.append(table.select("cpg", "beta").with_columns(sample=pl.lit(row["sample"])))
    long = pl.concat(long)
    counts = long.group_by("cpg").len()
    kept = counts.filter(pl.col("len") >= MIN_PRESENT * samples.height).select("cpg")
    wide = long.join(kept, on="cpg", how="semi").pivot(on="sample", index="cpg", values="beta").sort("cpg")
    chrom = wide["cpg"].str.split(":").list.get(0)
    pos = wide["cpg"].str.split(":").list.get(1).cast(pl.Int64)
    cpgs = pl.DataFrame({"cpg": wide["cpg"], "chr": chrom, "pos": pos}).with_columns(
        cpg_island=pl.Series(in_intervals(chrom.to_numpy(), pos.to_numpy() - 1, isl))
    )
    lay.samples.parent.mkdir(parents=True, exist_ok=True)
    samples.drop("url").write_parquet(lay.samples)
    wide.write_parquet(lay.betas)
    cpgs.write_parquet(lay.cpgs)
    return lay


def cpg_genes(cpgs: pl.DataFrame, prom: pl.DataFrame) -> pl.DataFrame:
    """``feature`` (CpG) and ``entrez`` for every promoter that holds it."""
    rows = []
    for (chrom,), sel in prom.group_by("chr"):
        c = cpgs.filter(pl.col("chr") == chrom)
        if c.height == 0:
            continue
        pos = c["pos"].to_numpy() - 1
        starts, ends, entrez = sel["start"].to_numpy(), sel["end"].to_numpy(), sel["entrez"].to_numpy()
        order = np.argsort(pos)
        sorted_pos = pos[order]
        names = c["cpg"].to_numpy()[order]
        lo = np.searchsorted(sorted_pos, starts, side="left")
        hi = np.searchsorted(sorted_pos, ends, side="left")
        for a, b, e in zip(lo, hi, entrez, strict=True):
            rows.extend((n, int(e)) for n in names[a:b])
    return pl.DataFrame(rows, schema={"feature": pl.String, "entrez": pl.Int64}, orient="row").unique()


# ----------------------------------------------------------------------------- Thompson et al. 2018 (GSE120132)
THOMPSON_MATRIX = "GSE120132-GPL21103_series_matrix.txt.gz"


@dataclass(frozen=True, slots=True)
class ThompsonLayout:
    """Mouse RRBS of six tissues (blood, liver, muscle, lung, adipose, kidney), two strains, 2-21 months."""

    paths: ProjectPaths

    @property
    def input(self) -> Path:
        return self.paths.input / "geo" / "GSE120132"

    @property
    def matrix(self) -> Path:
        return self.input / THOMPSON_MATRIX

    @property
    def samples(self) -> Path:
        return self.paths.interim / "cohorts" / "thompson_samples.parquet"

    def betas(self, tissue: str) -> Path:
        return self.paths.interim / "cohorts" / f"thompson_{tissue.lower()}_betas.parquet"

    def cpgs(self, tissue: str) -> Path:
        return self.paths.interim / "cohorts" / f"thompson_{tissue.lower()}_cpgs.parquet"


def read_thompson_samples(path: Path) -> pl.DataFrame:
    """Accession, title, tissue, strain, age in months (``age: 10mo``), sex and the table's URL."""
    meta: dict[str, list[str]] = {}
    keyed: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1"):
                cells = _split(line)[1:]
                keyed[cells[0].split(":", 1)[0].strip()] = [c.split(":", 1)[1].strip() for c in cells]
            elif line.startswith(("!Sample_geo_accession", "!Sample_title", "!Sample_supplementary_file_1")):
                cells = _split(line)
                meta[cells[0].removeprefix("!Sample_")] = cells[1:]
    missing = [k for k in ("tissue", "strain", "age", "Sex") if k not in keyed]
    if missing:
        raise ValueError(f"{path.name}: no characteristics {missing}")
    ages = [re.fullmatch(r"([\d.]+)\s*mo", a) for a in keyed["age"]]
    if not all(ages):
        raise ValueError(
            f"{path.name}: ages not in months, e.g. {[a for a, m in zip(keyed['age'], ages, strict=True) if not m][:3]}"
        )
    return pl.DataFrame(
        {
            "sample": meta["geo_accession"],
            "title": meta["title"],
            "tissue": keyed["tissue"],
            "strain": keyed["strain"],
            "age_months": [float(m.group(1)) for m in ages],
            "sex": keyed["Sex"],
            "url": [u.replace("ftp://", "https://") for u in meta["supplementary_file_1"]],
        }
    )


def read_cgfinal(path: Path) -> pl.DataFrame:
    """``cpg`` (``chr:pos``), ``beta`` and ``coverage`` of a BS_Seeker2 CGfinal table, chromosomes only.

    The table pools both strands of a CpG.
    """
    frame = pl.read_csv(
        path,
        separator=" ",
        has_header=False,
        new_columns=["chrom", "context", "pos", "methylated", "total", "frequency", "ci"],
        schema_overrides={"chrom": pl.String, "pos": pl.Int64, "methylated": pl.Int64, "total": pl.Int64},
    )
    chroms = [f"chr{c}" for c in _CHROMOSOMES.values()]
    return frame.filter(pl.col("chrom").is_in(chroms) & (pl.col("total") > 0)).select(
        pl.format("{}:{}", pl.col("chrom").str.strip_prefix("chr"), pl.col("pos")).alias("cpg"),
        (pl.col("methylated") / pl.col("total")).cast(pl.Float32).alias("beta"),
        pl.col("total").cast(pl.Int32).alias("coverage"),
    )


def fetch_thompson(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    """The series matrix and every sample's CGfinal table (about 6.4 GB); the mm10 tracks as for GSE80672."""
    lay = ThompsonLayout(paths)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120132/matrix/{THOMPSON_MATRIX}"
    files = [download(url, lay.matrix)]
    files += [
        download(UCSC.format(name=n), RrbsLayout(paths).track(n), overwrite=overwrite, timeout=600) for n in TRACKS
    ]
    for u in read_thompson_samples(lay.matrix)["url"]:
        files.append(download(u, lay.input / u.rsplit("/", 1)[1], overwrite=overwrite, timeout=1800))
    return files


def prepare_thompson(paths: ProjectPaths) -> ThompsonLayout:
    """Per tissue, the wide beta matrix of promoter CpGs covered in ``MIN_PRESENT`` of its samples."""
    lay = ThompsonLayout(paths)
    rr = RrbsLayout(paths)
    samples = read_thompson_samples(lay.matrix)
    prom = promoters(rr.track(TRACKS[0]), rr.track(TRACKS[1]))
    isl = islands(rr.track(TRACKS[2]))
    lay.samples.parent.mkdir(parents=True, exist_ok=True)
    samples.drop("url").write_parquet(lay.samples)
    for (tissue,), sel in samples.group_by("tissue", maintain_order=True):
        long = []
        for row in sel.iter_rows(named=True):
            table = read_cgfinal(lay.input / row["url"].rsplit("/", 1)[1]).filter(pl.col("coverage") >= MIN_COVERAGE)
            parts = table["cpg"].str.split(":")
            keep = in_intervals(parts.list.get(0).to_numpy(), parts.list.get(1).cast(pl.Int64).to_numpy() - 1, prom)
            long.append(table.filter(pl.Series(keep)).select("cpg", "beta").with_columns(sample=pl.lit(row["sample"])))
        long = pl.concat(long)
        kept = long.group_by("cpg").len().filter(pl.col("len") >= MIN_PRESENT * sel.height).select("cpg")
        wide = long.join(kept, on="cpg", how="semi").pivot(on="sample", index="cpg", values="beta").sort("cpg")
        chrom = wide["cpg"].str.split(":").list.get(0)
        pos = wide["cpg"].str.split(":").list.get(1).cast(pl.Int64)
        cpgs = pl.DataFrame({"cpg": wide["cpg"], "chr": chrom, "pos": pos}).with_columns(
            cpg_island=pl.Series(in_intervals(chrom.to_numpy(), pos.to_numpy() - 1, isl))
        )
        wide.write_parquet(lay.betas(tissue))
        cpgs.write_parquet(lay.cpgs(tissue))
    return lay
