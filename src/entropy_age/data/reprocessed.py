"""Betas rebuilt from raw IDATs with one SeSAMe recipe (roadmap RM41; ``scripts/sesame_preprocess.R``).

Each series is written to ``data/interim/sesame/<GSE>/<GSE>.betas.chunkNNN.parquet`` (``probe`` plus one column per
sample, named by GSM where the IDAT names carry one), ``<GSE>.qc.csv`` (per-sample SeSAMe QC with ``chip`` and
``position``) and ``<GSE>.session.txt``. This module reads them back.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.paths import ProjectPaths

MIN_DETECTED = 0.90  # samples with fewer probes passing pOOBAH are dropped (SeSAMe's usual cut)
#: an array whose mean intensity is below this share of the median of the arrays passing detection is dropped too
#: (protocol change 4, 2026-11-24): in EPIC-Italy 127 arrays at about a tenth of the median signal passed detection,
#: carried breast cases twice as often as controls, and drove a spurious drift association (exact-chip review)
MIN_INTENSITY_SHARE = 0.25


def series_dir(paths: ProjectPaths, gse: str, sesame_dir: str | None = None) -> Path:
    """``data/interim/sesame/<GSE>``; ``ENTROPY_AGE_SESAME_DIR`` (or ``sesame_dir``) names another rebuild (the clock
    track, ``sesame_clock``: SeSAMe without the design-flag mask, so clocks keep their CpGs)."""
    return paths.interim / (sesame_dir or os.environ.get("ENTROPY_AGE_SESAME_DIR", "sesame")) / gse


def available(paths: ProjectPaths, gse: str) -> bool:
    d = series_dir(paths, gse)
    return (d / f"{gse}.qc.csv").exists() and any(d.glob(f"{gse}.betas.chunk*.parquet"))


COLLAPSE_PLATFORMS = {"EPICv2"}  # replicate probes of one CpG, averaged to the cg id
MULTI_SPECIES = ("Mammal40",)  # probes for many species: in one species a quarter never hybridise


def species_detection(paths: ProjectPaths, gse: str) -> pl.DataFrame:
    """``sample`` and ``frac_dt_species``: each sample's detection over the probes this series measures (unmasked in
    at least ``MIN_DETECTED`` of its samples), for multi-species arrays where the whole-array rate reads the
    probes made for other species as failures (all of a mouse series fell below 90 %, 2026-11-18)."""

    frames = [pl.read_parquet(c) for c in sorted(series_dir(paths, gse).glob(f"{gse}.betas.chunk*.parquet"))]
    rows = []
    x = np.hstack([f.drop("probe").to_numpy() for f in frames])
    names = [c for f in frames for c in f.columns[1:]]
    measured = (~np.isnan(x)).mean(axis=1) >= MIN_DETECTED
    rate = (~np.isnan(x[measured])).mean(axis=0)
    rows = {"sample": names, "frac_dt_species": rate.tolist()}
    return pl.DataFrame(rows)


def read_qc(paths: ProjectPaths, gse: str, sesame_dir: str | None = None) -> pl.DataFrame:
    """Per-sample QC with ``passed``: detection at or above ``MIN_DETECTED`` over the whole array, or, on a
    multi-species array, over the probes SeSAMe's species mask keeps (``frac_dt_species``, written by the R script
    since RM43); rebuilds made before that fall back to the probes the series measures (circular, kept for them)."""
    qc = pl.read_csv(series_dir(paths, gse, sesame_dir) / f"{gse}.qc.csv", infer_schema_length=None)
    # one physical array per sample, the one the rebuild keeps in the betas (the first in file order): betas, QC,
    # detection and chip then describe the same array; a sample whose kept array fails detection is left out
    qc = qc.unique(subset="sample", keep="first", maintain_order=True)
    if "frac_dt_species" in qc.columns:  # SeSAMe's species mask decided the expected probes (RM43): not circular
        return _intensity(
            qc.with_columns(
                passed=pl.col("frac_dt_species").cast(pl.Float64) >= MIN_DETECTED, chip=pl.col("chip").cast(pl.String)
            )
        )
    if "platform" in qc.columns and set(qc["platform"].unique().to_list()) <= set(MULTI_SPECIES):
        qc = qc.join(species_detection(paths, gse), on="sample", how="left")
        passed = pl.col("frac_dt_species") >= MIN_DETECTED
    else:
        if "frac_dt" not in qc.columns:  # never let every array pass for want of a detection column
            raise ValueError(f"{gse}: the QC file has no detection column (frac_dt); the 90 % cut cannot be applied")
        passed = pl.col("frac_dt") >= MIN_DETECTED
    return _intensity(qc.with_columns(passed=passed, chip=pl.col("chip").cast(pl.String)))


def _intensity(qc: pl.DataFrame) -> pl.DataFrame:
    """``passed`` also requires a mean intensity of at least ``MIN_INTENSITY_SHARE`` of the passing arrays' median;
    ``low_intensity`` says which arrays this rule removed."""
    if "mean_intensity" not in qc.columns:
        return qc.with_columns(low_intensity=pl.lit(False))
    level = pl.col("mean_intensity").cast(pl.Float64, strict=False)
    median = qc.filter(pl.col("passed")).select(level.median()).item()
    low = level < MIN_INTENSITY_SHARE * median if median is not None else pl.lit(False)
    return qc.with_columns(low_intensity=pl.col("passed") & low.fill_null(False)).with_columns(
        passed=pl.col("passed") & ~pl.col("low_intensity")
    )


def collapse_probes(frame: pl.DataFrame) -> pl.DataFrame:
    """Probe ids to their CpG prefix (``cg00000029_TC21`` -> ``cg00000029``), replicates averaged (EPICv2, MM285)."""
    keyed = frame.with_columns(pl.col("probe").str.replace(r"_.*$", "").alias("probe"))
    if keyed["probe"].n_unique() == keyed.height:
        return keyed
    return keyed.group_by("probe", maintain_order=True).agg(pl.exclude("probe").mean())


def read_betas(paths: ProjectPaths, gse: str, *, passed_only: bool = True, collapse: bool = True) -> pl.DataFrame:
    """``probe`` plus one column per sample, chunks joined; failing samples left out unless ``passed_only`` is False."""
    chunks = sorted(series_dir(paths, gse).glob(f"{gse}.betas.chunk*.parquet"))
    frame = pl.read_parquet(chunks[0])
    for c in chunks[1:]:
        more = pl.read_parquet(c)
        more = more.select("probe", *(n for n in more.columns[1:] if n not in frame.columns))  # one array per sample
        frame = frame.join(more, on="probe", how="full", coalesce=True)
    if passed_only:
        keep = set(read_qc(paths, gse).filter(pl.col("passed"))["sample"].to_list())
        frame = frame.select("probe", *(c for c in frame.columns[1:] if c in keep))
    if not collapse:
        return frame
    qc = read_qc(paths, gse)
    platforms = set(qc["platform"].unique().to_list()) if "platform" in qc.columns else set()
    # only EPIC v2 repeats a CpG under suffixed ids; MM285's suffixes are its probe names, as the annotation has them
    return collapse_probes(frame) if platforms & COLLAPSE_PLATFORMS else frame


def in_use() -> bool:
    """True when the run reads rebuilt betas (``ENTROPY_AGE_BETAS=sesame``)."""
    return os.environ.get("ENTROPY_AGE_BETAS", "deposited") == "sesame"


def by_accession(paths: ProjectPaths, gse: str, names: dict[str, str]) -> pl.DataFrame:
    """``cpg`` plus one column per passing array, renamed from its GSM by ``names`` (GSM -> the experiment's sample
    name); arrays not in ``names`` are left out. Raises when the series was not rebuilt."""
    if not available(paths, gse):
        where = series_dir(paths, gse)
        raise FileNotFoundError(f"{gse} has no R SeSAMe rebuild in {where}; run scripts/sesame_preprocess.R")
    frame = read_betas(paths, gse)
    keep = [c for c in frame.columns[1:] if c in names]
    return frame.select(pl.col("probe").alias("cpg"), *(pl.col(c).alias(names[c]) for c in keep))


def matrix(frame: pl.DataFrame, features: list[str], samples: list[str]) -> np.ndarray:
    """``frame`` (``cpg`` plus sample columns) as samples by ``features``, NaN where a CpG is absent."""
    import numpy as np

    aligned = pl.DataFrame({"cpg": features}).join(frame.select("cpg", *samples), on="cpg", how="left")
    return aligned.select(samples).to_numpy().T.astype(float) if samples else np.empty((0, len(features)))


def cohort_betas(paths: ProjectPaths, gse: str, name: str, samples_path: Path) -> Path:
    """``data/interim/cohorts/<name>_betas_rsesame.parquet``: the rebuilt betas as the cohort table expects them
    (``cpg`` plus one column per sample id of ``<name>_samples.parquet``, GSM or chip_position), passing samples
    only; written on first use. Raises when the series was not rebuilt."""
    out = paths.interim / "cohorts" / f"{name}_betas_rsesame.parquet"  # not *_sesame: that name is R61's methylprep
    if out.exists():
        return out
    if not available(paths, gse):
        raise FileNotFoundError(
            f"{gse} has no R SeSAMe rebuild in {series_dir(paths, gse)}; run scripts/sesame_preprocess.R"
        )
    qc = read_qc(paths, gse).filter(pl.col("passed"))
    frame = read_betas(paths, gse)
    ids = set(pl.read_parquet(samples_path)["sample"].to_list())
    titles = {k: _fibre_id(v, ids) for k, v in _titles(paths, gse).items()} if name == "fibres" else _titles(paths, gse)
    rename = {}
    for row in qc.iter_rows(named=True):
        position = f"{row['chip']}_{row['position']}"
        if row["sample"] in ids:
            rename[row["sample"]] = row["sample"]
        elif position in ids:
            rename[row["sample"]] = position
        elif titles.get(row["sample"]) in ids:
            rename[row["sample"]] = titles[row["sample"]]
    keep = [c for c in frame.columns[1:] if c in rename]
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.select(pl.col("probe").alias("cpg"), *(pl.col(c).alias(rename[c]) for c in keep)).write_parquet(out)
    return out


def _titles(paths: ProjectPaths, gse: str) -> dict[str, str]:
    """GSM accession -> GEO sample title, from the series matrix (empty when there is none)."""
    import gzip

    matrix = paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
    if not matrix.exists():
        return {}
    acc = titles = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_geo_accession"):
                acc = cells
            if acc and titles:
                break
    return dict(zip(acc or [], titles or [], strict=False))


FIBRE = {"WM": "WM", "TII": "T2", "TI": "T1"}


def _fibre_id(title: str, ids: set[str] | frozenset[str] = frozenset()) -> str:
    """GEO's fibre title (``s203_pre_WM_rep1``) as the authors' matrix names the sample: ``WM.s203-preWM`` for a single
    array, ``WM.s201-preWM-2`` for technical replicate 2 where the table keeps replicates apart (``ids`` says which).
    A title that matches no table name is returned unchanged, so it matches nothing."""
    import re

    m = re.fullmatch(r"([a-z]\d+)_(pre|4w|12w)_(WM|TII|TI)_rep(\d+)", title)
    if not m:
        return title
    fibre = FIBRE[m.group(3)]
    base = f"{fibre}.{m.group(1)}-{m.group(2)}{fibre}"
    for candidate in (f"{base}-{m.group(4)}", base if m.group(4) == "1" else None):
        if candidate is not None and (candidate in ids or (not ids and candidate == base)):
            return candidate
    return title


def conversion_index(paths: ProjectPaths, gse: str, samples: list[str]) -> np.ndarray:
    """SeSAMe's bisulfite-conversion index (GCT) of ``samples`` (the ids of the QC table), in that order."""
    qc = read_qc(paths, gse)
    gct = dict(qc.select("sample", pl.col("GCT").cast(pl.Float64, strict=False)).iter_rows())
    return np.array([gct[s] for s in samples], dtype=float)


def conversion_by_group(gct: np.ndarray, carrier: np.ndarray) -> dict:
    """Mean conversion index of carriers and controls, their difference in pooled SDs and Welch's p."""
    from scipy import stats

    a, b = gct[carrier], gct[~carrier]
    pooled = float(np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2))
    return {
        "gct_carriers": float(a.mean()),
        "gct_controls": float(b.mean()),
        "difference_sd": float((a.mean() - b.mean()) / pooled) if pooled > 0 else float("nan"),
        "p_welch": float(stats.ttest_ind(a, b, equal_var=False).pvalue),
    }
