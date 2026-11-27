"""Published mouse clocks on single-cell data pooled per mouse: do the mice whose cells drift more look older?

Anton asked whether the clocks were checked at pseudobulk. Two single-cell sets carry the per-cell drift readings:

**Methylation: muscle stem cells (GSE121437, ``musc_single_cell``).** scBS-seq of 275 cells from three 2-month and
three 24-month mice. Each cell's Bismark coverage lists each strand's C apart; the reads at the two Cs of a CpG are
pooled (the CpG is keyed by its top-strand C, 1-based, mm10). Per mouse, every cell's reads are summed into a
pseudobulk beta per CpG. Scored with the mouse RRBS clocks pyaging ships (``scripts/export_mouse_clocks.py``):
Thompson 2018 (multi-tissue, months), Meer 2018 (multi-tissue, days), Petkovich 2017 (blood) and Stubbs 2017
(multi-tissue, quantile-normalised to its gold standard), each with its published conversion to months. A CpG counts
as measured in a mouse with at least ``MIN_PB_READS`` pooled reads; a clock's prediction is valid only where the
measured CpGs carry at least ``MIN_WEIGHT_SHARE`` of its absolute weight (the project's coverage rule). Missing CpGs
are filled with the mean of the other mice, else the Thompson bulk muscle mean, else the clock's own reference (or
0.5); a valid prediction still leans on fills for up to a tenth of the weight (``filled_weight_share``), and for
Stubbs every fill also moves the quantile normalisation. Two clock inputs (one per strand) can name one CpG; both get
its reads. The rule is repeated with any pooled read counting as measured (``<clock>_reads1``). As a positive
control every clock also scores the bulk muscle RRBS of Thompson et al. (GSE120132, 60 mice of two strains, 2-20
months; Thompson's clock was trained on these data, so for it the control is circular).

Single cells are scored with scAge (Trapp et al. 2021): per CpG, a linear model of bulk methylation on age (here
fitted on the 30 C57BL/6J muscle samples of GSE120132, CpGs with ``SCAGE_MIN_COV`` reads in ``SCAGE_MIN_SHARE`` of
them); each cell's CpGs are binarised (methylated share above 0.5 -> 1, below -> 0, exactly 0.5 dropped); the
``SCAGE_TOP`` CpGs of that cell with the largest absolute age correlation are kept, and the cell's age is the
maximum of the Bernoulli likelihood over ``SCAGE_GRID`` (-20 to 60 months by 0.1, as published), with the predicted
level clipped to [0.001, 0.999]; a cell with fewer than ``SCAGE_TOP`` usable CpGs is not scored. Differences from
Trapp's run on these cells: fixed CpG counts (250, 500, 1000) instead of the top 1 %, binarisation after pooling the
two strands, no 500,000-CpG cell threshold (these cells have far fewer) and no reference outlier screen. Bulk muscle is
mostly myonuclei, a different cell type from the stem cells, so scAge here cannot give calibrated stem-cell ages.

The per-mouse shield gap: per mouse, the between-cell variance of each promoter's methylated share over that mouse's
cells (``musc_single_cell.between_cell_variance``, promoters with ``MIN_READS`` reads in at least
``MIN_CELLS_MOUSE`` of its cells), housekeeping and core-essential promoters against promoters of the same island
status, young-level decile and coverage fifth (the pooled analysis's strata), in SD of that mouse's variances.
Negative = set promoters differ less between the mouse's cells than matched promoters. The link: each mouse's
acceleration (prediction minus its age-group mean) against its gap minus its age-group mean, Pearson over the six
mice, with the exact permutation p over the 36 relabellings within age groups.

**Expression: Tabula Muris Senis FACS (``tms_cell_noise``).** Per tissue x cell type eligible there (3 months against
18-24 months), every mouse with ``MIN_CELLS_PER_MOUSE`` cells gets a pseudobulk: raw Smart-seq2 counts summed over its
cells, all genes but ERCC spike-ins. Scored with tAge (Tyshkovskiy et al. 2026, scaled-difference multi-tissue
model, through pyaging's coefficients): genes with 10 counts in 20 % of the cell type's mice, mapped to mouse Entrez,
RLE-normalised, log10(x + 1), z-scored per mouse, centred on the median of the 3-month mice; the prediction is an
age difference in months against the 3-month mice. Measured = a clock gene that passes the filter; the coverage rule
as above. Per cell type: Spearman of prediction on age, and old (18-24 m) minus young mean prediction; pooled over
valid cell types, with an interval from resampling mice within age group (one draw shared by all cell types, the
whole tAge preparation redone per draw; for valid units a resample counts only when every unit is present and still
passes the coverage rule).

Within one cell type many of tAge's weighted genes are not expressed (the smoke run kept 34-65 % of its weight), so
the same scoring is also run on a whole-tissue pseudobulk per mouse (every FACS cell of the tissue, tissues with
``MIN_MICE`` young and old mice of ``MIN_CELLS_PER_MOUSE`` cells), closer to the bulk tissues tAge was trained on but
carrying each mouse's sorted cell mix.

The per-mouse noise gap: per cell type, cells thinned to one depth as in ``tms_cell_noise`` (the
``DEPTH_QUANTILE`` quantile of their protein-coding read totals; shallower cells dropped), genes detected in
``MIN_DETECTED`` of cells kept; per mouse, residual log CV2 on its mean trend over that mouse's cells, housekeeping
(and core-essential) genes against background genes of the same mean decile and dropout tertile
(``tms_cell_noise.matched_diff``), in SD of the mouse's residuals. Acceleration and gap are each residualised on age
group and sex within the cell type (or tissue) and z-scored, then averaged per mouse; the link is Spearman over mice
of mean acceleration (from cell types or from tissues) against mean noise gap (from cell types), with a mouse
bootstrap interval (fitted scores held fixed) and a permutation p within exact ages (3, 18, 21, 24 months; the 21-month
mice are all female and the 24-month mice all male, so age and sex are partly confounded).

With six and nineteen mice the links can show a direction and little more. Output ``pseudobulk_clocks``; run as a job
with ``uv run python -c`` and a ``Context`` (``ENTROPY_AGE_WORKERS`` processes, ``ENTROPY_AGE_BOOT`` resamples).
"""

from __future__ import annotations

import itertools
import json
import multiprocessing
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, rankdata, spearmanr

from entropy_age.data.rrbs import ThompsonLayout, read_thompson_samples
from entropy_age.experiments import musc_single_cell as musc
from entropy_age.experiments import tms_cell_noise as tms
from entropy_age.experiments.common import Context, write_table
from entropy_age.protection import matched_scores

METH_CLOCKS = ("thompson", "meer", "petkovich", "stubbs")
SHIELD_SETS = ("housekeeping", "core_essential")
MIN_WEIGHT_SHARE = 0.9
MIN_PB_READS = 5  # pooled reads for a CpG to count as measured in a mouse's pseudobulk
MIN_BULK_READS = 5  # the same for a bulk RRBS sample
MIN_CELLS_MOUSE = 10  # cells of one mouse a promoter needs for its per-mouse variance
MIN_PROMOTERS_MOUSE = 1000  # promoters a mouse needs for its gap (one mouse's shallow cells pass only a few dozen)
PB_READS_SENSITIVITY = 1  # the coverage rule again with any pooled read counting as measured
MUSC_AGE = {"Y": 2.0, "O": 24.0}  # months, from the series matrix
SCAGE_MIN_COV = 5
SCAGE_MIN_SHARE = 0.9
SCAGE_TOP = 500
SCAGE_TOPS = (250, 500, 1000)
SCAGE_GRID = np.round(np.arange(-20.0, 60.0 + 1e-9, 0.1), 1)  # months, as Trapp et al. 2021
SCAGE_CLIP = 0.001
TAGE_COUNTS = 10
TAGE_PERCENT = 20
TMS_AGE = {"3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0}
CHROM_CODE = {str(i): i for i in range(1, 20)} | {"X": 20, "Y": 21, "M": 22, "MT": 22}
KEY_BASE = 1_000_000_000  # a CpG's integer key: chromosome code * KEY_BASE + top-strand C position


def key_expr(chrom: str = "chr", pos: str = "pos") -> pl.Expr:
    """Integer CpG key from a chromosome column (no ``chr`` prefix) and a position; null off the main chromosomes."""
    code = pl.col(chrom).replace_strict(CHROM_CODE, default=None, return_dtype=pl.Int64)
    return (code * KEY_BASE + pl.col(pos)).alias("key")


def key_of(text: str) -> int:
    """``N:pos`` -> the integer key."""
    chrom, pos = text.split(":")
    return CHROM_CODE[chrom] * KEY_BASE + int(pos)


# ----------------------------------------------------------------------------- clocks
@dataclass(frozen=True)
class Clock:
    """A linear clock from ``scripts/export_mouse_clocks.py``: input keys, weights and the conversion."""

    name: str
    keys: np.ndarray  # methylation: "chr:pos" of the CpG's top-strand C (no "chr"); tAge: mouse Entrez as str
    weight: np.ndarray
    intercept: float
    reference: np.ndarray  # pyaging's fill value per input, NaN where none
    sd: np.ndarray | None = None  # Stubbs: gold-standard SDs
    used: np.ndarray | None = None  # Stubbs: inputs that carry a weight

    @property
    def weighted(self) -> np.ndarray:
        return self.weight != 0


def canonical_key(feature: str, strand: str) -> str:
    """``chrN:pos`` (1-based, either strand's C) -> ``N:pos`` of the CpG's top-strand C."""
    chrom, pos = feature.split(":")
    p = int(pos) - (1 if strand == "-" else 0)
    if strand not in {"+", "-"}:
        raise ValueError(f"{feature}: not a CpG on mm10 (strand {strand!r})")
    return f"{chrom.removeprefix('chr')}:{p}"


def load_clock(folder: Path, name: str) -> Clock:
    frame = pl.read_csv(folder / f"{name}.csv", schema_overrides={"feature": pl.String})
    meta = json.loads((folder / f"{name}.json").read_text())
    if name == "tage":
        keys = frame["feature"].to_numpy().astype(str)
    else:
        keys = np.array([canonical_key(f, s or "") for f, s in frame.select("feature", "strand").iter_rows()])
    reference = frame["reference"].cast(pl.Float64).fill_null(np.nan).to_numpy()
    sd = frame["sd"].to_numpy().astype(float) if "sd" in frame.columns else None
    used = np.asarray(meta["used_index"], dtype=np.int64) if "used_index" in meta else None
    return Clock(name, keys, frame["weight"].to_numpy().astype(float), float(meta["intercept"]), reference, sd, used)


def quantile_to_target(x: np.ndarray, target_sorted: np.ndarray) -> np.ndarray:
    """Each row given the sorted ``target`` by rank, ties as preprocessCore (mean of the targets either side of
    the tie group's average rank), as pyaging's Stubbs preprocessing."""
    out = np.empty_like(x, dtype=float)
    for i, row in enumerate(x):
        a = rankdata(row, method="average") - 1.0
        lo, hi = np.floor(a).astype(np.int64), np.ceil(a).astype(np.int64)
        out[i] = 0.5 * (target_sorted[lo] + target_sorted[hi])
    return out


def linear_predictor(clock: Clock, x: np.ndarray) -> np.ndarray:
    """The clock's linear predictor for samples x inputs ``x`` (no missing values left)."""
    if np.isnan(x).any():
        raise ValueError(f"{clock.name}: missing inputs; fill them first")
    if clock.name == "stubbs":
        gold = clock.reference
        target = np.sort(gold)[np.linspace(0, gold.size - 1, x.shape[1]).astype(np.int64)]
        sd = np.where(np.abs(clock.sd) < 1e-9, 1.0, clock.sd)
        z = (quantile_to_target(x, target) - gold) / sd
        return z[:, clock.used] @ clock.weight[clock.used] + clock.intercept
    return x @ clock.weight + clock.intercept


def to_months(name: str, x: np.ndarray) -> np.ndarray:
    """Each clock's published conversion of its linear predictor (pyaging's post-processing) to months."""
    x = np.asarray(x, dtype=float)
    if name == "thompson":
        return x
    if name == "meer":
        return x / 30.5
    if name == "petkovich":
        a, b, c = 0.1666, 0.4185, -1.712
        with np.errstate(invalid="ignore"):
            return ((x - c) / a) ** (1 / b) / 30.5
    if name == "stubbs":
        return (np.exp(0.1207 * x**2 + 1.2424 * x + 2.5440) - 3) * 7 / 30.5
    raise ValueError(f"no conversion for {name}")


def weight_share(clock: Clock, measured: np.ndarray) -> np.ndarray:
    """Per sample (rows of ``measured``, samples x inputs): share of the clock's absolute weight on measured
    inputs."""
    w = np.abs(clock.weight)
    return measured.astype(float) @ w / w.sum()


def fill_missing(beta: np.ndarray, fills: list[np.ndarray]) -> np.ndarray:
    """NaN entries of ``beta`` (samples x inputs) taken from the first finite value in ``fills`` (each per input),
    then 0.5."""
    out = beta.copy()
    for f in [*fills, np.full(beta.shape[1], 0.5)]:
        hole = np.isnan(out)
        if not hole.any():
            break
        out = np.where(hole, np.broadcast_to(f, out.shape), out)
    return out


def leave_one_out_mean(beta: np.ndarray) -> np.ndarray:
    """Per sample and input: the mean of the other samples' finite values (NaN where none)."""
    ok = np.isfinite(beta)
    total = np.where(ok, beta, 0.0).sum(axis=0)
    count = ok.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (total - np.where(ok, beta, 0.0)) / (count - ok.astype(int))


# ----------------------------------------------------------------------------- scAge
def scage_reference(beta: np.ndarray, age: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per CpG (columns of samples x CpGs ``beta``, NaN = not measured): intercept, slope and Pearson r of the
    methylated share on age over the samples that measured it."""
    ok = np.isfinite(beta)
    a = np.where(ok, age[:, None], 0.0)
    b = np.where(ok, beta, 0.0)
    n = ok.sum(axis=0).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma, mb = a.sum(axis=0) / n, b.sum(axis=0) / n
        saa = (a**2).sum(axis=0) - n * ma**2
        sbb = (b**2).sum(axis=0) - n * mb**2
        sab = (a * b).sum(axis=0) - n * ma * mb
        slope = sab / saa
        r = sab / np.sqrt(saa * sbb)
    return mb - slope * ma, slope, r


def scage_cell(
    state: np.ndarray, intercept: np.ndarray, slope: np.ndarray, r: np.ndarray, top: int = SCAGE_TOP
) -> tuple[float, int]:
    """Maximum-likelihood age (months, over ``SCAGE_GRID``) of one cell from its binary CpG states and the
    reference models of the same CpGs; the ``top`` CpGs by absolute r are used. Returns the age and CpGs used; the
    age is NaN for a cell with fewer than ``top`` usable CpGs."""
    ok = np.isfinite(r) & np.isfinite(state)
    if ok.sum() < top:
        return float("nan"), int(ok.sum())
    idx = np.flatnonzero(ok)
    idx = idx[np.argsort(-np.abs(r[idx]), kind="stable")[:top]]
    p = np.clip(intercept[idx, None] + slope[idx, None] * SCAGE_GRID[None, :], SCAGE_CLIP, 1 - SCAGE_CLIP)
    s = state[idx, None]
    ll = (s * np.log(p) + (1 - s) * np.log(1 - p)).sum(axis=0)
    return float(SCAGE_GRID[int(np.argmax(ll))]), int(idx.size)


def binarise(meth: np.ndarray, total: np.ndarray) -> np.ndarray:
    """1 where the methylated share is above 0.5, 0 below, NaN at exactly 0.5 or without reads."""
    with np.errstate(invalid="ignore", divide="ignore"):
        share = meth / total
    out = np.where(share > 0.5, 1.0, np.where(share < 0.5, 0.0, np.nan))
    return np.where(total > 0, out, np.nan)


# ----------------------------------------------------------------------------- links
def within_group_link(x: np.ndarray, y: np.ndarray, group: np.ndarray) -> dict:
    """Pearson of ``x`` and ``y`` after each is centred on its group mean, with the exact permutation p (two-sided)
    over every relabelling of ``y`` within groups (all of them when there are at most 50,000)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, group = x[ok], y[ok], np.asarray(group)[ok]
    xc, yc = x.copy(), y.copy()
    for g in np.unique(group):
        sel = group == g
        xc[sel] -= x[sel].mean()
        yc[sel] -= y[sel].mean()
    out = {"n": int(ok.sum()), "r": float("nan"), "perm_p": float("nan"), "permutations": 0}
    if ok.sum() < 3 or np.std(xc) == 0 or np.std(yc) == 0:
        return out
    r = float(np.corrcoef(xc, yc)[0, 1])
    blocks = [np.flatnonzero(group == g) for g in np.unique(group)]
    sizes = [int(np.prod(range(1, b.size + 1))) for b in blocks]
    if int(np.prod(sizes)) > 50_000:
        raise ValueError("within_group_link: too many permutations for the exact test")
    null = []
    for perm in itertools.product(*[itertools.permutations(b) for b in blocks]):
        yp = yc.copy()
        for b, p in zip(blocks, perm, strict=True):
            yp[b] = yc[list(p)]
        null.append(np.corrcoef(xc, yp)[0, 1])
    null = np.array(null)
    out |= {"r": r, "perm_p": float(np.mean(np.abs(null) >= abs(r) - 1e-12)), "permutations": int(null.size)}
    return out


def direction(r: float) -> str:
    if not np.isfinite(r):
        return "not estimable"
    return "older-looking mice have the larger gap" if r > 0 else "older-looking mice have the smaller gap"


# ----------------------------------------------------------------------------- methylation inputs
def _cov(path: Path) -> pl.DataFrame:
    return (
        pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=["chr", "pos", "end", "pct", "meth", "unmeth"],
            schema_overrides={"chr": pl.String, "pct": pl.Float64},
        )
        .with_columns(pl.col("chr").str.strip_prefix("chr"))
        .select(
            key_expr(),
            pl.col("meth").cast(pl.Int32),
            (pl.col("meth") + pl.col("unmeth")).cast(pl.Int32).alias("total"),
        )
        .drop_nulls("key")
    )


def to_cpgs(cov: pl.DataFrame, keys: pl.DataFrame) -> pl.DataFrame:
    """Strand rows of a coverage table (``key``, ``meth``, ``total``; the key built from each row's own position)
    pooled onto the CpGs in ``keys`` (``key`` of the top-strand C): a row at a CpG's top C carries its key, a row
    at its bottom C (one base on) the key plus one."""
    top = cov.join(keys, on="key", how="inner")
    bottom = cov.with_columns(pl.col("key") - 1).join(keys, on="key", how="inner")
    return pl.concat([top, bottom]).group_by("key").agg(pl.col("meth").sum(), pl.col("total").sum())


def read_musc_cells(folder: Path, keys: pl.DataFrame) -> pl.DataFrame:
    """``cell``, ``mouse``, ``key``, ``meth``, ``total`` over the CpGs in ``keys`` for every cell."""
    frames = []
    for path in sorted(folder.glob("*.cov.gz")):
        found = musc.CELL.search(path.name)
        if found is None:
            continue
        cell = f"{found.group(1)}_TA_Hi_{found.group(2)}"
        frames.append(
            to_cpgs(_cov(path), keys).with_columns(pl.lit(cell).alias("cell"), pl.lit(found.group(1)).alias("mouse"))
        )
    return pl.concat(frames)


def _bulk_samples(ctx: Context) -> pl.DataFrame:
    lay = ThompsonLayout(ctx.paths)
    samples = read_thompson_samples(lay.matrix).filter(pl.col("tissue") == "Muscle")
    files = {p.name.split("_")[0]: str(p) for p in lay.input.glob("GSM*_Mouse_Muscle_*.CGfinal.txt.gz")}
    return samples.with_columns(path=pl.col("sample").replace_strict(files, default=None)).drop_nulls("path")


def read_bulk(path: Path, min_reads: int) -> pl.DataFrame:
    """``key`` (top-strand C, as the CGfinal table keys a pooled CpG) and ``beta`` of CpGs with ``min_reads``
    reads."""
    frame = pl.read_csv(
        path,
        separator=" ",
        has_header=False,
        new_columns=["chr", "context", "pos", "meth", "total", "frequency", "ci"],
        schema_overrides={"chr": pl.String, "pos": pl.Int64, "meth": pl.Int64, "total": pl.Int64},
    )
    return (
        frame.filter(pl.col("total") >= min_reads)
        .with_columns(pl.col("chr").str.strip_prefix("chr"))
        .select(key_expr(), (pl.col("meth") / pl.col("total")).cast(pl.Float32).alias("beta"))
        .drop_nulls("key")
    )


def _matrix(long: pl.DataFrame, rows: list[str], row_col: str, keys: np.ndarray, value: str) -> np.ndarray:
    """Rows x ``keys`` (integer CpG keys) of ``value`` from a long table with ``key``, NaN where absent."""
    keys = np.asarray(keys, dtype=np.int64)
    # two clock inputs (one per strand) can name the same CpG: fill the unique keys, then expand (Codex review 13)
    unique, inverse = np.unique(keys, return_inverse=True)
    frame = long.filter(
        pl.col("key").is_in(pl.Series(unique).implode()) & pl.col(row_col).is_in(pl.Series(rows).implode())
    )
    row = {r: i for i, r in enumerate(rows)}
    out = np.full((len(rows), unique.size), np.nan, dtype=np.float32 if value == "beta" else np.float64)
    i = frame[row_col].replace_strict(row, return_dtype=pl.Int64).to_numpy()
    j = np.searchsorted(unique, frame["key"].to_numpy())
    out[i, j] = frame[value].to_numpy()
    if unique.size != keys.size:
        out = out[:, inverse]
    return out.astype(np.float64) if out.size < 50_000_000 else out


def score_methylation(
    clocks: dict[str, Clock],
    beta: dict[str, np.ndarray],
    measured: dict[str, np.ndarray],
    fills: dict[str, list[np.ndarray]],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per clock: (prediction in months, weight share measured) per sample."""
    out = {}
    for name, clock in clocks.items():
        x = fill_missing(beta[name], fills[name])
        out[name] = (to_months(name, linear_predictor(clock, x)), weight_share(clock, measured[name]))
    return out


def per_mouse_gaps(data: musc.Promoters) -> pl.DataFrame:
    """Per mouse and set: the between-cell variance gap (set minus matched promoters, SD units) over that mouse's
    cells."""
    point = musc.readings(data.k, data.n, data.old)
    keep, variants = musc.strata_of(point, data)
    strata = variants["island x level"]
    rows = []
    for m in np.unique(data.mouse):
        sel = data.mouse == m
        _, variance, cells, reads = musc.between_cell_variance(data.k[sel], data.n[sel])
        ok = keep & (cells >= MIN_CELLS_MOUSE) & np.isfinite(variance)
        sd = float(np.std(variance[ok])) if ok.any() else float("nan")
        for s in SHIELD_SETS:
            gap = float("nan")
            if ok.sum() >= MIN_PROMOTERS_MOUSE:
                gap = float(matched_scores(variance[ok][None, :], data.sets[s][ok], strata[ok])[0][0]) / sd
            rows.append(
                {
                    "mouse": str(m),
                    "age_months": MUSC_AGE[str(m)[0]],
                    "set": s,
                    "heterogeneity_gap_sd": gap,
                    "promoters": int(ok.sum()),
                    "promoters_in_set": int((ok & data.sets[s]).sum()),
                    "cells": int(sel.sum()),
                    "mean_reads": float(np.nanmean(reads[ok])),
                }
            )
    return pl.DataFrame(rows)


def run_methylation(ctx: Context) -> list[Path]:
    clock_dir = ctx.paths.input / "clocks" / "mouse"
    folder = musc.folder_of(ctx)
    lay = ThompsonLayout(ctx.paths)
    if not folder.exists() or not (clock_dir / "thompson.csv").exists() or not lay.matrix.exists():
        print("pseudobulk_clocks: GSE121437 cells, GSE120132 or the exported clocks are missing", file=sys.stderr)
        return []
    out = ctx.out_dir("pseudobulk_clocks")
    clocks = {c: load_clock(clock_dir, c) for c in METH_CLOCKS}

    # bulk muscle: the scAge reference (C57BL/6J) and the clocks' positive control (both strains)
    samples = _bulk_samples(ctx)
    clock_ints = {name: np.array([key_of(k) for k in c.keys], dtype=np.int64) for name, c in clocks.items()}
    clock_keys = pl.DataFrame({"key": np.unique(np.concatenate(list(clock_ints.values())))})
    bl6_ids = set(samples.filter(pl.col("strain") == "C57BL/6J")["sample"].to_list())
    bulk_parts, ref_parts = [], []
    for s, p in samples.select("sample", "path").iter_rows():
        frame = read_bulk(p, min(SCAGE_MIN_COV, MIN_BULK_READS))
        bulk_parts.append(frame.join(clock_keys, on="key").with_columns(pl.lit(s).alias("sample")))
        if s in bl6_ids:
            ref_parts.append(frame.filter(pl.col("beta").is_not_null()).with_columns(pl.lit(s).alias("sample")))
    bulk = pl.concat(bulk_parts)
    ref_long = pl.concat(ref_parts)
    del ref_parts
    bl6 = samples.filter(pl.col("strain") == "C57BL/6J")
    need = int(np.ceil(SCAGE_MIN_SHARE * bl6.height))
    ref_keys = ref_long.group_by("key").len().filter(pl.col("len") >= need).select("key").sort("key")
    ref_long = ref_long.join(ref_keys, on="key")
    ref_beta = _matrix(ref_long, bl6["sample"].to_list(), "sample", ref_keys["key"].to_numpy(), "beta")
    del ref_long
    intercept, slope, r = scage_reference(ref_beta, bl6["age_months"].to_numpy().astype(float))
    del ref_beta
    print(f"[pseudobulk_clocks] scAge reference: {ref_keys.height} CpGs over {bl6.height} mice", file=sys.stderr)

    keys = pl.concat([clock_keys, ref_keys]).unique()
    cells = read_musc_cells(folder, keys)

    # bulk control for every clock (bulk tables keep CpGs with MIN_BULK_READS reads)
    bulk_ids = samples["sample"].to_list()
    bulk_mean = {}
    control = []
    for name, clock in clocks.items():
        b = _matrix(bulk, bulk_ids, "sample", clock_ints[name], "beta")
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            bulk_mean[name] = np.nanmean(b, axis=0)
        fills = [leave_one_out_mean(b), clock.reference]
        pred = to_months(name, linear_predictor(clock, fill_missing(b, fills)))
        share = weight_share(clock, np.isfinite(b))
        age = samples["age_months"].to_numpy().astype(float)
        valid = share >= MIN_WEIGHT_SHARE
        for strain in ("all", "C57BL/6J", "BALB/cByJ"):
            sel = valid & ((samples["strain"] == strain).to_numpy() if strain != "all" else True)
            rr = pearsonr(pred[sel], age[sel])[0] if sel.sum() > 2 and np.std(pred[sel]) > 0 else np.nan
            control.append(
                {
                    "clock": name,
                    "strain": strain,
                    "samples": int(sel.sum()),
                    "median_weight_share": float(np.median(share)),
                    "pearson_r": float(rr),
                    "median_abs_error_months": float(np.median(np.abs(pred[sel] - age[sel]))) if sel.any() else np.nan,
                    **{
                        f"mean_pred_{int(a)}m": float(np.mean(pred[sel & (age == a)]))
                        if (sel & (age == a)).any()
                        else np.nan
                        for a in (2.0, 10.0, 20.0)
                    },
                }
            )
    write_table(pl.DataFrame(control), out / "meth_bulk_control")

    # pseudobulk per mouse
    mice = sorted(cells["mouse"].unique().to_list())
    pb = cells.group_by("mouse", "key").agg(pl.col("meth").sum(), pl.col("total").sum())
    pb = pb.with_columns(beta=pl.col("meth") / pl.col("total"))
    mouse_age = np.array([MUSC_AGE[m[0]] for m in mice])
    rows = []
    for (name, clock), min_reads in itertools.product(clocks.items(), (MIN_PB_READS, PB_READS_SENSITIVITY)):
        label = name if min_reads == MIN_PB_READS else f"{name}_reads{min_reads}"
        total = np.nan_to_num(_matrix(pb, mice, "mouse", clock_ints[name], "total"))
        beta = _matrix(pb, mice, "mouse", clock_ints[name], "beta").astype(np.float64)
        measured = total >= min_reads
        beta = np.where(measured, beta, np.nan)
        fills = [leave_one_out_mean(beta), bulk_mean[name], clock.reference]
        pred = to_months(name, linear_predictor(clock, fill_missing(beta, fills)))
        share = weight_share(clock, measured)
        share_any = weight_share(clock, total > 0)
        for i, m in enumerate(mice):
            rows.append(
                {
                    "mouse": m,
                    "age_months": mouse_age[i],
                    "clock": label,
                    "weight_share": float(share[i]),
                    "filled_weight_share": float(1 - share[i]),
                    "weight_share_any_read": float(share_any[i]),
                    "valid": bool(share[i] >= MIN_WEIGHT_SHARE),
                    "predicted_months": float(pred[i]),
                    "cpgs_measured": int(measured[i][clock.weighted].sum()),
                    "cpgs_weighted": int(clock.weighted.sum()),
                }
            )

    # scAge per cell
    ref_sorted = ref_keys["key"].to_numpy()
    cell_rows = []
    sc = cells.join(ref_keys, on="key")
    for (cell, mouse), part in sc.group_by("cell", "mouse"):
        state = np.full(ref_sorted.size, np.nan)
        j = np.searchsorted(ref_sorted, part["key"].to_numpy())
        state[j] = binarise(part["meth"].to_numpy().astype(float), part["total"].to_numpy().astype(float))
        row = {"cell": cell, "mouse": mouse, "age_months": MUSC_AGE[mouse[0]], "cpgs": int(np.isfinite(state).sum())}
        for top in SCAGE_TOPS:
            age_hat, used = scage_cell(state, intercept, slope, r, top)
            row[f"scage_top{top}"] = age_hat
            row[f"used_top{top}"] = used
        cell_rows.append(row)
    scage = pl.DataFrame(cell_rows).sort("mouse", "cell")
    write_table(scage, out / "scage_cells")
    for top in SCAGE_TOPS:
        per = (
            scage.group_by("mouse")
            .agg(
                pl.col(f"scage_top{top}").drop_nans().mean().alias("pred"), pl.col(f"scage_top{top}").is_not_nan().sum()
            )
            .sort("mouse")
        )
        for m, value, scored in per.iter_rows():
            rows.append(
                {
                    "mouse": m,
                    "age_months": MUSC_AGE[m[0]],
                    "clock": f"scage_top{top}",
                    "weight_share": float("nan"),
                    "weight_share_any_read": float("nan"),
                    "valid": bool(scored >= MIN_CELLS_MOUSE),
                    "predicted_months": float(value) if value is not None else float("nan"),
                    "cpgs_measured": int(scored),  # cells scored, for scAge
                    "cpgs_weighted": top,
                }
            )
    mice_table = pl.DataFrame(rows).with_columns(
        acceleration=pl.col("predicted_months")
        - pl.when(pl.col("valid")).then(pl.col("predicted_months")).mean().over("clock", "age_months")
    )
    write_table(mice_table, out / "meth_clock_mice")

    # scAge, old minus young: cells, and mouse means with a bootstrap over cells within mouse
    rng = ctx.rng(281_100)
    summary = []
    for top in SCAGE_TOPS:
        col = f"scage_top{top}"
        scored = scage.filter(pl.col(col).is_not_nan())
        v = scored[col].to_numpy()
        mouse_of = scored["mouse"].to_numpy()
        kept = [m for m in mice if (mouse_of == m).sum() >= MIN_CELLS_MOUSE]
        means = {m: float(v[mouse_of == m].mean()) for m in kept}

        def contrast(values: dict[str, float], drop: str | None = None) -> float:
            old = [x for m, x in values.items() if m.startswith("O") and m != drop]
            young = [x for m, x in values.items() if m.startswith("Y") and m != drop]
            return float(np.mean(old) - np.mean(young)) if old and young else float("nan")

        boot = []
        for _ in range(2000):
            boot.append(contrast({m: rng.choice(v[mouse_of == m], (mouse_of == m).sum()).mean() for m in kept}))
        lomo = [contrast(means, m) for m in kept]
        old_cells = np.isin(mouse_of, [m for m in kept if m.startswith("O")])
        young_cells = np.isin(mouse_of, [m for m in kept if m.startswith("Y")])
        summary.append(
            {
                "top_cpgs": top,
                "mice": len(kept),
                "cells_young": int(young_cells.sum()),
                "cells_old": int(old_cells.sum()),
                "cells_unscored": int(scage.height - scored.height),
                "mean_young_cells": float(v[young_cells].mean()),
                "mean_old_cells": float(v[old_cells].mean()),
                "old_minus_young_mouse_means": contrast(means),
                "cell_boot_low": float(np.quantile(boot, 0.025)),
                "cell_boot_high": float(np.quantile(boot, 0.975)),
                "leave_one_mouse_out_min": float(np.nanmin(lomo)),
                "leave_one_mouse_out_max": float(np.nanmax(lomo)),
                "mouse_means": json.dumps({m: round(means[m], 2) for m in kept}),
            }
        )
    write_table(pl.DataFrame(summary), out / "scage_summary")

    # per-mouse shield gap and the link
    gaps = per_mouse_gaps(musc.load_promoters(ctx))
    write_table(gaps, out / "musc_mouse_gaps")
    links = []
    for clock in mice_table["clock"].unique().sort().to_list():
        t = mice_table.filter(pl.col("clock") == clock).sort("mouse")
        for s in SHIELD_SETS:
            g = gaps.filter(pl.col("set") == s).sort("mouse")
            joined = t.join(g.select("mouse", "heterogeneity_gap_sd"), on="mouse")
            valid = joined.filter(pl.col("valid"))
            res = within_group_link(
                valid["acceleration"].to_numpy(),
                valid["heterogeneity_gap_sd"].to_numpy(),
                valid["age_months"].to_numpy(),
            )
            links.append({"clock": clock, "set": s, **res, "direction": direction(res["r"])})
    write_table(pl.DataFrame(links), out / "meth_links")
    print(f"[pseudobulk_clocks] methylation: {len(mice)} mice, {scage.height} cells", file=sys.stderr, flush=True)
    return []


# ----------------------------------------------------------------------------- tAge
def tage_prepare(
    counts: np.ndarray, symbols: np.ndarray, lookup: dict[str, str], reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """tAge's cohort preparation of raw counts (samples x genes named by ``symbols``) as pyaging's ``_prepare_tage``
    (Tyshkovskiy 2026, scaled difference): genes with ``TAGE_COUNTS`` counts in ``TAGE_PERCENT`` % of samples,
    mapped to mouse Entrez (several genes on one id summed), RLE-normalised, log10(x + 1), z-scored per sample and
    centred on the median of the ``reference`` rows. Returns the matrix and its Entrez ids."""
    counts = np.asarray(counts, dtype=float)
    passing = (counts >= TAGE_COUNTS).sum(axis=0) >= counts.shape[0] * (TAGE_PERCENT / 100)
    target = np.array([lookup.get(str(s).upper(), "") for s in symbols])
    keep = passing & (target != "")
    ids, inverse = np.unique(target[keep], return_inverse=True)
    mapped = np.zeros((counts.shape[0], ids.size))
    np.add.at(mapped.T, inverse, counts[:, keep].T)
    lib = mapped.sum(axis=1)
    with np.errstate(divide="ignore"):
        geo = np.exp(np.log(mapped).mean(axis=0))
    usable = geo > 0
    if not usable.any():
        raise ValueError("tage_prepare: no gene expressed in every sample")
    factors = np.median(mapped[:, usable] / geo[usable], axis=1) / lib
    factors /= np.exp(np.log(factors).mean())
    x = np.log10(mapped / (lib * factors)[:, None] * 1e7 + 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        x = (x - x.mean(axis=1, keepdims=True)) / x.std(axis=1, ddof=1, keepdims=True)
    if not np.isfinite(x).all():
        raise ValueError("tage_prepare: a sample without variation across genes")
    return x - np.median(x[reference], axis=0), ids


def tage_score(clock: Clock, x: np.ndarray, ids: np.ndarray) -> tuple[np.ndarray, float]:
    """Predictions (months against the reference) and the share of absolute weight on genes present."""
    col = {k: j for j, k in enumerate(ids)}
    j = np.array([col.get(k, -1) for k in clock.keys])
    present = j >= 0
    inputs = np.tile(clock.reference, (x.shape[0], 1))
    inputs[:, present] = x[:, j[present]]
    inputs = np.where(np.isnan(inputs), np.nan_to_num(clock.reference), inputs)
    w = np.abs(clock.weight)
    return inputs @ clock.weight + clock.intercept, float(w[present].sum() / w.sum())


def tage_lookup(path: Path) -> dict[str, str]:
    """Upper-case source id (symbol, Ensembl or Entrez) -> mouse Entrez, mouse rows in file order (later wins)."""
    frame = pl.read_csv(path, schema_overrides={"source_id": pl.String, "mouse_entrez": pl.String})
    frame = frame.filter(pl.col("species") == "mouse")
    return dict(zip(frame["source_id"].str.to_uppercase().to_list(), frame["mouse_entrez"].to_list(), strict=True))


def residual_z(y: np.ndarray, age: np.ndarray, sex: np.ndarray) -> np.ndarray:
    """``y`` residualised on age group and sex dummies (OLS), divided by the residual SD; NaN when no degrees of
    freedom are left."""
    ok = np.isfinite(y)
    design = [np.ones(y.size)]
    for values in (age, sex):
        for level in np.unique(values)[1:]:
            design.append((values == level).astype(float))
    x = np.column_stack(design)
    out = np.full(y.size, np.nan)
    if ok.sum() <= np.linalg.matrix_rank(x[ok]) + 1:
        return out
    beta, *_ = np.linalg.lstsq(x[ok], y[ok], rcond=None)
    res = y[ok] - x[ok] @ beta
    sd = res.std(ddof=np.linalg.matrix_rank(x[ok]))
    out[ok] = res / sd if sd > 1e-9 * max(1.0, float(np.abs(y[ok]).max())) else np.nan
    return out


@dataclass
class ExprTask:
    path: Path
    tissue: str
    cell_type: str
    cells: pl.DataFrame  # row, mouse, age, sex
    all_index: np.ndarray  # every non-ERCC gene
    all_symbols: np.ndarray
    gene_index: np.ndarray  # protein-coding, mapped (positions in X)
    sets: dict[str, np.ndarray]  # over gene_index
    clock_dir: Path
    draws: list[dict[str, list[str]]]
    seed: int
    noise: bool = True  # False for a whole-tissue pseudobulk: tAge only


def analyse_expression(task: ExprTask) -> tuple[list[dict], list[dict], list[dict]]:
    """Per mouse rows (tAge and noise gaps), mouse-resample rows, and a cell-type row."""
    print(f"[pseudobulk_clocks] {task.tissue} / {task.cell_type}", file=sys.stderr, flush=True)
    clock = load_clock(task.clock_dir, "tage")
    lookup = tage_lookup(task.clock_dir / "tage_gene_mapping.csv.gz")
    cells = task.cells.sort("row")
    raw_all = tms.read_rows(task.path, cells["row"].to_numpy(), tms._n_genes(task.path))
    counts = raw_all[:, task.all_index].tocsr()
    coding = raw_all[:, task.gene_index].tocsr() if task.noise else None
    del raw_all
    mouse = cells["mouse"].to_numpy()
    per_mouse = cells.group_by("mouse").agg(pl.len().alias("n"), pl.col("age").first(), pl.col("sex").first())
    per_mouse = per_mouse.filter(pl.col("n") >= tms.MIN_CELLS_PER_MOUSE).sort("mouse")
    mice = per_mouse["mouse"].to_numpy()
    age = per_mouse["age"].replace_strict(TMS_AGE, return_dtype=pl.Float64).to_numpy()
    sex = per_mouse["sex"].to_numpy()
    young = age == 3.0
    if young.sum() < 2 or (~young).sum() < 2:
        return [], [], []
    pseudo = np.vstack([np.asarray(counts[mouse == m].sum(axis=0)).ravel() for m in mice])
    x, ids = tage_prepare(pseudo, task.all_symbols, lookup, young)
    pred, share = tage_score(clock, x, ids)

    gap = {s: np.full(mice.size, np.nan) for s in SHIELD_SETS}
    depth, genes_noise = 0, 0
    if task.noise:
        gap, depth, genes_noise = _noise_gaps(coding, mouse, mice, task)
    accel = residual_z(pred, age, sex)
    gap_z = {s: residual_z(gap[s], age, sex) for s in SHIELD_SETS}
    rows = []
    for i, m in enumerate(mice):
        rows.append(
            {
                "tissue": task.tissue,
                "cell_type": task.cell_type,
                "mouse": str(m),
                "age_months": float(age[i]),
                "sex": str(sex[i]),
                "cells": int(per_mouse["n"][i]),
                "tage_months": float(pred[i]),
                "acceleration_z": float(accel[i]),
                **{f"noise_gap_{s}": float(gap[s][i]) for s in SHIELD_SETS},
                **{f"noise_gap_{s}_z": float(gap_z[s][i]) for s in SHIELD_SETS},
            }
        )
    rho = spearmanr(age, pred)[0] if np.std(pred) > 0 else np.nan
    summary = {
        "tissue": task.tissue,
        "cell_type": task.cell_type,
        "mice": int(mice.size),
        "mice_young": int(young.sum()),
        "mice_old": int((~young).sum()),
        "genes_after_filter": int(ids.size),
        "weight_share": share,
        "valid": bool(share >= MIN_WEIGHT_SHARE),
        "spearman_age": float(rho),
        "old_minus_young_months": float(pred[~young].mean() - pred[young].mean()),
        "slope_per_month": float(np.polyfit(age, pred, 1)[0]),
        "depth": depth,
        "genes_noise": genes_noise,
    }
    # mouse resamples: the cohort redrawn within age group, the tAge preparation redone
    boot = []
    index = {m: i for i, m in enumerate(mice)}
    for b, draw in enumerate(task.draws):
        pick = [index[m] for g in ("young", "old") for m in draw[g] if m in index]
        yb = young[pick]
        if yb.sum() < 2 or (~yb).sum() < 2:
            continue
        try:
            xb, idb = tage_prepare(pseudo[pick], task.all_symbols, lookup, yb)
        except ValueError:
            continue
        pb, share_b = tage_score(clock, xb, idb)
        boot.append(
            {
                "tissue": task.tissue,
                "cell_type": task.cell_type,
                "boot": b,
                "old_minus_young_months": float(pb[~yb].mean() - pb[yb].mean()),
                "weight_share": share_b,
            }
        )
    return rows, boot, [summary]


def _noise_gaps(coding, mouse: np.ndarray, mice: np.ndarray, task: ExprTask):
    """Per mouse and set: residual log CV2 gap over the mouse's depth-matched cells (SD of its residuals)."""
    rng = np.random.default_rng(task.seed)
    totals = np.asarray(coding.sum(axis=1)).ravel()
    depth = int(np.quantile(totals, tms.DEPTH_QUANTILE))
    deep = totals >= depth
    thinned = tms.thin_rows(coding[deep], depth, rng).toarray()
    mouse_d = mouse[deep]
    detected = (thinned > 0).mean(axis=0) >= tms.MIN_DETECTED
    xg = thinned[:, detected]
    sets = {s: m[detected] for s, m in task.sets.items()}
    background = ~np.any(np.column_stack(list(sets.values())), axis=1)
    stratum = tms.strata(xg.mean(axis=0), (xg == 0).mean(axis=0))
    gap = {s: np.full(mice.size, np.nan) for s in SHIELD_SETS}
    for i, m in enumerate(mice):
        sel = mouse_d == m
        if sel.sum() < tms.MIN_CELLS_PER_MOUSE:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            res, _, _ = tms.noise_residual(xg[sel], mouse_d[sel])
        sd = float(np.nanstd(res))
        for s in SHIELD_SETS:
            gap[s][i] = tms.matched_diff(res, sets[s], background, stratum) / sd
    return gap, depth, int(detected.sum())


def eligible_tissues(obs: pl.DataFrame) -> list[str]:
    """Tissues with ``MIN_MICE`` young and old mice of ``MIN_CELLS_PER_MOUSE`` cells each (all cell types)."""
    per = obs.group_by("tissue", "age_group", "mouse").len().filter(pl.col("len") >= tms.MIN_CELLS_PER_MOUSE)
    mice = per.group_by("tissue", "age_group").len().filter(pl.col("len") >= tms.MIN_MICE)
    return sorted(mice.group_by("tissue").len().filter(pl.col("len") == 2)["tissue"].to_list())


def _link(accel: pl.DataFrame, gaps: pl.DataFrame, s: str, rng: np.random.Generator) -> dict:
    """Spearman over mice of mean acceleration z (rows of ``accel``) and mean noise-gap z (rows of ``gaps``)."""
    a_m = accel.group_by("mouse").agg(
        pl.col("acceleration_z").drop_nans().mean(), pl.col("age_months").first(), pl.len().alias("units")
    )
    g_m = gaps.group_by("mouse").agg(pl.col(f"noise_gap_{s}_z").drop_nans().mean().alias("gap_z"))
    joined = a_m.join(g_m, on="mouse").drop_nulls().drop_nans().sort("mouse")
    if joined.height < 4:
        return {"mice": joined.height}
    a, g = joined["acceleration_z"].to_numpy(), joined["gap_z"].to_numpy()
    group = joined["age_months"].to_numpy()  # exact age (the 21 m mice are all female, the 24 m all male)
    rho, p = _perm_spearman(a, g, group, rng)
    low, high = _spearman_boot(a, g, rng)
    return {
        "mice": joined.height,
        "spearman": rho,
        "mouse_boot_low": low,
        "mouse_boot_high": high,
        "perm_p_within_age": p,
        "direction": direction(rho),
    }


def _spearman_boot(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, n: int = 2000) -> tuple[float, float]:
    vals = []
    for _ in range(n):
        i = rng.integers(0, x.size, x.size)
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            vals.append(spearmanr(x[i], y[i])[0])
    return float(np.nanquantile(vals, 0.025)), float(np.nanquantile(vals, 0.975))


def _perm_spearman(x: np.ndarray, y: np.ndarray, group: np.ndarray, rng: np.random.Generator, n: int = 10_000):
    r = spearmanr(x, y)[0]
    null = np.empty(n)
    for k in range(n):
        yp = y.copy()
        for g in np.unique(group):
            sel = np.flatnonzero(group == g)
            yp[sel] = y[rng.permutation(sel)]
        null[k] = spearmanr(x, yp)[0]
    return float(r), float(np.mean(np.abs(null) >= abs(r) - 1e-12))


def run_expression(ctx: Context) -> list[Path]:
    clock_dir = ctx.paths.input / "clocks" / "mouse"
    path, needed = tms._inputs(ctx)
    if any(not p.exists() for p in needed) or not (clock_dir / "tage.csv").exists():
        print("pseudobulk_clocks: TMS FACS object, annotation or tAge coefficients missing", file=sys.stderr)
        return []
    out = ctx.out_dir("pseudobulk_clocks")
    workers, n_boot = tms._settings()
    data = tms._load(ctx)
    obs = tms._age_groups(data.obs, tms.YOUNG, tms.OLD)
    groups = tms.eligible(obs)
    if limit := int(os.environ.get("ENTROPY_AGE_TMS_CELL_TYPES", "0")):
        groups = groups.head(limit)
    all_index = np.flatnonzero(~np.char.startswith(data.symbols.astype(str), "ERCC-"))
    mice = {g: sorted(obs.filter(pl.col("age_group") == g)["mouse"].unique().to_list()) for g in ("young", "old")}
    rng = ctx.rng(281_200)
    draws = [tms.mouse_draw(mice, rng) for _ in range(n_boot)]
    tasks = [
        ExprTask(
            path,
            tissue,
            cell_type,
            obs.filter((pl.col("tissue") == tissue) & (pl.col("cell_type") == cell_type)).select(
                "row", "mouse", "age", "sex"
            ),
            all_index,
            data.symbols[all_index],
            data.gene_index,
            data.sets,
            clock_dir,
            draws,
            2000 + i,
        )
        for i, (tissue, cell_type) in enumerate(groups.iter_rows())
    ]
    tissues = eligible_tissues(obs)
    if limit:
        tissues = tissues[:1]
    tasks += [
        ExprTask(
            path,
            tissue,
            "all cells",
            obs.filter(pl.col("tissue") == tissue).select("row", "mouse", "age", "sex"),
            all_index,
            data.symbols[all_index],
            data.gene_index,
            data.sets,
            clock_dir,
            draws,
            3000 + i,
            noise=False,
        )
        for i, tissue in enumerate(tissues)
    ]
    print(f"[pseudobulk_clocks] {len(tasks)} cell types, {n_boot} mouse resamples", file=sys.stderr, flush=True)
    if workers > 1:
        # spawn: a forked child inherits polars' locked thread pool (tms_cell_noise, 2026-11-26)
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            results = list(pool.map(analyse_expression, tasks))
    else:
        results = [analyse_expression(t) for t in tasks]
    per_mouse = pl.DataFrame([r for rows, _, _ in results for r in rows], infer_schema_length=None)
    boot = pl.DataFrame([r for _, b, _ in results for r in b], infer_schema_length=None)
    units = pl.DataFrame([r for _, _, s in results for r in s], infer_schema_length=None)
    units = units.with_columns(
        level=pl.when(pl.col("cell_type") == "all cells").then(pl.lit("tissue")).otherwise(pl.lit("cell type"))
    )
    cell_types = units
    write_table(per_mouse, out / "tage_mice")
    write_table(cell_types, out / "tage_cell_types")
    boot.write_parquet(out / "tage_resamples.parquet")

    rng = ctx.rng(281_300)
    pooled = []
    for (level, label), which in _subsets(units):
        keys = which.select("tissue", "cell_type")
        point = float(which["old_minus_young_months"].mean()) if which.height else float("nan")
        # a resample counts only with every unit present and, for valid units, every one still above the rule
        per_boot = (
            boot.join(keys, on=["tissue", "cell_type"])
            .group_by("boot")
            .agg(
                pl.col("old_minus_young_months").mean(),
                pl.len().alias("k"),
                (pl.col("weight_share") >= MIN_WEIGHT_SHARE).all().alias("all_pass"),
            )
            .filter((pl.col("k") == which.height) & (pl.col("all_pass") | pl.lit(label == "all")))
        )
        pooled.append(
            {
                "level": level,
                "units": label,
                "n": which.height,
                "old_minus_young_months": point,
                "mouse_low": float(per_boot["old_minus_young_months"].quantile(0.025)) if per_boot.height else np.nan,
                "mouse_high": float(per_boot["old_minus_young_months"].quantile(0.975)) if per_boot.height else np.nan,
                "resamples": per_boot.height,
                "median_spearman_age": float(which["spearman_age"].median()) if which.height else np.nan,
                "share_positive_spearman": float((which["spearman_age"] > 0).mean()) if which.height else np.nan,
                "median_weight_share": float(units.filter(pl.col("level") == level)["weight_share"].median()),
            }
        )
    write_table(pl.DataFrame(pooled), out / "tage_pooled")

    gaps = per_mouse.join(
        units.filter(pl.col("level") == "cell type").select("tissue", "cell_type"), on=["tissue", "cell_type"]
    )
    links = []
    for (level, label), which in _subsets(units):
        accel = per_mouse.join(which.select("tissue", "cell_type"), on=["tissue", "cell_type"])
        for s in SHIELD_SETS:
            links.append({"level": level, "units": label, "set": s, **_link(accel, gaps, s, rng)})
    write_table(pl.DataFrame(links, infer_schema_length=None), out / "tage_links")
    return []


def _subsets(units: pl.DataFrame):
    """((level, valid|all), rows) for cell types and whole tissues."""
    for level in ("cell type", "tissue"):
        rows = units.filter(pl.col("level") == level)
        yield (level, "valid"), rows.filter(pl.col("valid"))
        yield (level, "all"), rows


def run(ctx: Context) -> list[Path]:
    return run_methylation(ctx) + run_expression(ctx)
