"""Faster loss beyond the pooled fit: per-CpG gaps with a donor interval (roadmap RM38 item 2).

R67 (``equilibrium_drift``) found that essential promoters lose methylation faster than look-alikes from pooled
rates per set, with donor intervals, and from per-CpG bounded relaxation fits, with a permutation p only. A
simulated referee on round v0.1 asked whether the faster loss rests on the pooled fit. This reruns the per-CpG
reading (``equilibrium_drift._per_cpg``: each CpG's relaxation fit, the set against look-alikes of the same island
status and young-level decile, and again within its own chromatin state) on donor resamples, each resample's gap
divided by the full sample's SD (``common.rescale``), so every per-CpG gap gets a 95 % interval.

``entropy-age figures per_cpg_interval`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.experiments.common import Context, rescale, write_table
from entropy_age.experiments.equilibrium_drift import OWN_STATE, PERMUTATIONS, _per_cpg, _states, adjust
from entropy_age.experiments.tissue_contrast import _features, _hannum, _jaffe, _settings

N_BOOT_MAX = 100
KEYS = ["tissue", "matching", "metric", "group"]


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("per_cpg_interval")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))  # equilibrium_drift's CpGs
    _, n_boot = _settings(ctx)
    n_boot = min(n_boot, N_BOOT_MAX)
    scale = str(ctx.scale.name)
    rng = ctx.rng(380_002)
    point, boot = [], []
    for p in (_hannum(ctx.paths, features), _jaffe(ctx.paths, features)):
        x = np.clip(adjust(p.x, p.ages.astype(float), p.covariates), 0, 1).astype(np.float32)
        state = _states(ctx, p.features, probes, OWN_STATE[p.label])
        _, g, _ = _per_cpg(p, x, probes, membership, PERMUTATIONS[scale], ctx.seed, state)
        point.append(g.with_columns(tissue=pl.lit(p.label), donors=pl.lit(x.shape[0])))
        idx = np.arange(x.shape[0])
        for b in range(n_boot):
            draw = np.sort(rng.choice(idx, size=idx.size, replace=True))
            q = dataclasses.replace(p, ages=p.ages[draw])
            _, gb, _ = _per_cpg(q, x[draw], probes, membership, 0, b, state)
            boot.append(gb.select("matching", "metric", "group", "matched_diff").with_columns(tissue=pl.lit(p.label)))
        print(f"per_cpg_interval: {p.label}: {idx.size} donors, {n_boot} resamples", file=sys.stderr, flush=True)
    point = pl.concat(point, how="diagonal")
    if boot:
        ci = (
            rescale(pl.concat(boot), point, KEYS)
            .group_by(KEYS)
            .agg(lo=pl.col("matched_diff_sd").quantile(0.025), hi=pl.col("matched_diff_sd").quantile(0.975))
        )
        point = point.join(ci, on=KEYS, how="left")
    write_table(point.sort(KEYS), out / "per_cpg_gaps")
    return []
