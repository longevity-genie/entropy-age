"""Betas rebuilt from raw IDATs with one SeSAMe recipe (roadmap RM41; ``scripts/sesame_preprocess.R``).

Each series is written to ``data/interim/sesame/<GSE>/<GSE>.betas.chunkNNN.parquet`` (``probe`` plus one column per
sample, named by GSM where the IDAT names carry one), ``<GSE>.qc.csv`` (per-sample SeSAMe QC with ``chip`` and
``position``) and ``<GSE>.session.txt``. This module reads them back.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from entropy_age.paths import ProjectPaths

MIN_DETECTED = 0.90  # samples with fewer probes passing pOOBAH are dropped (SeSAMe's usual cut)


def series_dir(paths: ProjectPaths, gse: str) -> Path:
    return paths.interim / "sesame" / gse


def available(paths: ProjectPaths, gse: str) -> bool:
    d = series_dir(paths, gse)
    return (d / f"{gse}.qc.csv").exists() and any(d.glob(f"{gse}.betas.chunk*.parquet"))


def read_qc(paths: ProjectPaths, gse: str) -> pl.DataFrame:
    """Per-sample QC with ``passed`` (detection at or above ``MIN_DETECTED``)."""
    qc = pl.read_csv(series_dir(paths, gse) / f"{gse}.qc.csv", infer_schema_length=None)
    detected = "frac_dt" if "frac_dt" in qc.columns else None
    passed = pl.col(detected) >= MIN_DETECTED if detected else pl.lit(True)
    return qc.with_columns(passed=passed, chip=pl.col("chip").cast(pl.String))


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
        frame = frame.join(pl.read_parquet(c), on="probe", how="full", coalesce=True)
    if passed_only:
        keep = set(read_qc(paths, gse).filter(pl.col("passed"))["sample"].to_list())
        frame = frame.select("probe", *(c for c in frame.columns[1:] if c in keep))
    return collapse_probes(frame) if collapse else frame
