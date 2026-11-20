"""Lunnon's five tissues normalised together from the raw signal intensities (roadmap RM43 item 7).

The deposited betas of GSE59685 were normalised with wateRmelon's dasen "for each of the tissues separately", so a
difference between blood and a brain region of the same donors carries whatever the per-tissue normalisation did.
GEO also holds the raw methylated and unmethylated intensities of every array
(``GSE59685_signal_intensities.csv.gz``, from GenomeStudio, no background correction). Here all arrays are
normalised in one pass: within each probe design (Infinium I, II), the methylated and the unmethylated intensities
are each quantile-normalised across all arrays, and beta = M / (M + U + 100). This is dasen's core step applied
jointly; it is not a rebuild from IDATs (none exist), and slide still goes with tissue, so a same-donor contrast
between tissues stays batch-confounded whatever it shows.

``build`` writes ``data/interim/lunnon_joint/betas.parquet`` (``cpg`` plus one column per GEO sample). With
``ENTROPY_AGE_LUNNON=joint`` the cohort table of ``lunnon`` reads it and outputs carry ``_lunnonjoint``.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.paths import ProjectPaths

INTENSITIES = "GSE59685_signal_intensities.csv.gz"
OFFSET = 100.0
MAX_DETECTION_P = 0.01


def in_use() -> bool:
    return os.environ.get("ENTROPY_AGE_LUNNON", "deposited") == "joint"


def betas_path(paths: ProjectPaths) -> Path:
    return paths.interim / "lunnon_joint" / "betas.parquet"


def quantile_normalise(x: np.ndarray) -> np.ndarray:
    """Columns of ``x`` (probes x arrays) forced to one distribution: the mean of the sorted columns. NaN stays NaN
    and takes no part (a column with missing values is mapped by its ranks among its own finite values)."""
    out = np.full(x.shape, np.nan)
    finite = np.isfinite(x)
    full = finite.all(axis=0)
    reference = (
        np.sort(x[:, full], axis=0).mean(axis=1) if full.any() else np.sort(np.nan_to_num(x), axis=0).mean(axis=1)
    )
    grid = (np.arange(reference.size) + 0.5) / reference.size
    for j in range(x.shape[1]):
        rows = np.flatnonzero(finite[:, j])
        ranks = np.argsort(np.argsort(x[rows, j], kind="stable"), kind="stable")
        out[rows, j] = np.interp((ranks + 0.5) / rows.size, grid, reference)
    return out


def joint_betas(m: np.ndarray, u: np.ndarray, design: np.ndarray, detection: np.ndarray | None = None) -> np.ndarray:
    """Betas from methylated ``m`` and unmethylated ``u`` intensities (probes x arrays), each quantile-normalised
    across arrays within each probe ``design``; probes failing detection in an array are NaN there."""
    mn, un = np.full(m.shape, np.nan), np.full(u.shape, np.nan)
    for kind in np.unique(design):
        rows = design == kind
        mn[rows], un[rows] = quantile_normalise(m[rows]), quantile_normalise(u[rows])
    beta = mn / (mn + un + OFFSET)
    if detection is not None:
        beta[detection > MAX_DETECTION_P] = np.nan
    return beta


def build(paths: ProjectPaths, *, overwrite: bool = False) -> Path:
    """Read the intensities, normalise all arrays jointly and write the beta table named by GEO sample."""
    from entropy_age.data.annotation import AnnotationLayout
    from entropy_age.data.geo import read_series_matrix

    out = betas_path(paths)
    if out.exists() and not overwrite:
        return out
    folder = paths.input / "geo" / "GSE59685"
    meta = read_series_matrix(folder / "GSE59685_series_matrix.txt.gz").samples
    sample_of = dict(zip(meta["barcode"].to_list(), meta["sample"].to_list(), strict=True))
    path = folder / INTENSITIES
    with gzip.open(path, "rt") as handle:
        skip = 0
        for line in handle:
            if line.startswith(","):
                break
            skip += 1
    frame = pl.read_csv(path, skip_rows=skip, null_values=["NA", ""], infer_schema_length=1000)
    frame = frame.rename({frame.columns[0]: "cpg"})
    barcodes = [b for b in sample_of if f"{b} Methylated Signal" in frame.columns]
    probes = pl.read_parquet(AnnotationLayout(paths).probes_450k).select("cpg", "probe_type")
    frame = frame.join(probes, on="cpg", how="inner")

    def block(suffix: str) -> np.ndarray:
        return frame.select([f"{b} {suffix}" for b in barcodes]).to_numpy().astype(float)

    beta = joint_betas(
        block("Methylated Signal"),
        block("Unmethylated Signal"),
        frame["probe_type"].to_numpy(),
        block("Detection Pval"),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    table = pl.DataFrame({"cpg": frame["cpg"], **{sample_of[b]: beta[:, j] for j, b in enumerate(barcodes)}})
    table.write_parquet(out)
    return out
