"""The shield away from the 0 % floor (roadmap RM38 item 1; a simulated referee's concern on round v0.1).

Essential island promoters sit near 0 % methylation, where a CpG cannot drift down and its spread is bounded. The
matching strata (island x young-level decile) hold the level, but the lowest decile still spans 0-3 % in blood, so
part of the gap could be the floor. This reruns the shield's matched gap (``tissue_contrast``'s features,
covariates, metrics and ``_matched_sd``) on CpGs whose young level is inside a band away from both bounds:
10-90 % and, stricter, 20-80 %. R59 did the same for the species array. Hannum blood and Jaffe cortex; a donor
bootstrap gives the intervals.

``entropy-age figures floor_check`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.experiments.common import Context, rescale, write_table
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    Prepared,
    _features,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _matched_sd,
    _settings,
    _with_strata,
)

BANDS = {"all levels": (0.0, 1.0), "10-90 %": (0.10, 0.90), "20-80 %": (0.20, 0.80)}
N_BOOT_MAX = 100


def in_band(metrics: pl.DataFrame, band: tuple[float, float]) -> pl.DataFrame:
    lo, hi = band
    return metrics.filter(pl.col("baseline").is_between(lo, hi))


def gaps(p: Prepared, rows: np.ndarray, probes, membership, permutations: int, seed: int) -> pl.DataFrame:
    cov = None if p.covariates is None else p.covariates[rows]
    q = dataclasses.replace(p, x=p.x[rows], ages=p.ages[rows], covariates=cov, donors=p.donors[rows])
    m = _lifespan_metrics(q)
    out = []
    for name, band in BANDS.items():
        t = _matched_sd(_with_strata(in_band(m, band), probes), membership, LIFESPAN_METRICS, permutations, seed)
        out.append(t.with_columns(band=pl.lit(name)))
    return pl.concat(out, how="diagonal")


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("floor_check")
    _, n_boot = _settings(ctx)
    n_boot = min(n_boot, N_BOOT_MAX)
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    rng = ctx.rng(38_001)
    point, boot = [], []
    for loader in (_hannum, _jaffe):
        p = loader(ctx.paths, features)
        p = dataclasses.replace(p, x=p.x.astype(np.float32))
        idx = np.arange(p.x.shape[0])
        point.append(gaps(p, idx, probes, membership, 1000, 1).with_columns(tissue=pl.lit(p.label)))
        for b in range(n_boot):
            draw = np.sort(rng.choice(idx, size=idx.size, replace=True))
            t = gaps(p, draw, probes, membership, 0, b)
            boot.append(t.select("band", "metric", "group", "matched_diff").with_columns(tissue=pl.lit(p.label)))
        print(f"floor_check: {p.label}", file=sys.stderr, flush=True)
    point = pl.concat(point, how="diagonal")
    keys = ["tissue", "band", "metric", "group"]
    if boot:
        ci = (
            rescale(pl.concat(boot), point, keys)  # every resample on the full sample's sd (RM9; audit 2026-11-18)
            .group_by(keys)
            .agg(lo=pl.col("matched_diff_sd").quantile(0.025), hi=pl.col("matched_diff_sd").quantile(0.975))
        )
        point = point.join(ci, on=keys, how="left")
    write_table(point.sort(keys), out / "gaps_by_band")
    return []
