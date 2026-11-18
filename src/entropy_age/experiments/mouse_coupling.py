"""Do the organs of one mouse drift together, at one age? (roadmap RM39 item 1)

R41 (``organ_coupling``) correlated age acceleration across tissues of the same human donors (Lunnon, GTEx). The
literature check (docs/LITERATURE_COVERAGE.md, idea 1) left open clock-free drift coupling, and mice of one fixed
age, where no spread in age can couple two tissues. Huang et al. (GSE279437) profiled up to 11 tissues of the
same mice; its untreated controls are two fixed-age groups (young, aged). Per tissue and age group, every control
mouse is scored against the other controls of its group and tissue (``consortium_interventions.arm_scores``'s
leave-one-out: drift ``B``, shield share ``F``, signed shift ``S`` per gene set), standardised within the group.
For each pair of tissues and score: Pearson r over the mice that gave both, the groups pooled (each centred).
The null shuffles mice within age group in one tissue of every pair (``N_PERM`` draws) for the mean r over pairs.

``entropy-age figures mouse_coupling`` as a Slurm job.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.health_states import MEASURES, loo_scores
from entropy_age.experiments.mouse_tissues import _features, _tissues
from entropy_age.experiments.per_person import _islands, _sets

CONTROLS = ("Young_ctrl", "Aged_ctrl")
MIN_SHARED = 10
N_PERM = {"paper": 1000, "quick": 100}


def pair_r(a: pl.DataFrame, b: pl.DataFrame, key: str) -> tuple[float, int]:
    j = a.select("donor", pl.col(key).alias("a")).join(b.select("donor", pl.col(key).alias("b")), on="donor")
    if j.height < MIN_SHARED:
        return float("nan"), j.height
    return float(np.corrcoef(j["a"], j["b"])[0, 1]), j.height


def shuffled(t: pl.DataFrame, rng) -> pl.DataFrame:
    """The same scores with mouse labels permuted within age group."""
    parts = []
    for _, g in t.group_by("group", maintain_order=True):
        parts.append(g.with_columns(donor=pl.Series(rng.permutation(g["donor"].to_numpy()))))
    return pl.concat(parts)


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("mouse_coupling")
    features, probes, membership = _features(ctx.paths, ctx.rng(130_000))
    per_tissue: dict[str, pl.DataFrame] = {}
    for name, t in _tissues(ctx.paths, features).items():
        island = _islands(t.features, probes)
        sets = _sets(t.features, membership)
        parts = []
        for group in CONTROLS:
            rows = np.flatnonzero(t.groups == group)
            if rows.size < 5:
                continue
            sc = loo_scores(t.x[rows], island, sets)
            keys = [k for k in sc if k.split("_")[0] in MEASURES]
            frame = pl.DataFrame({"donor": t.donors[rows].astype(str), "group": group, **{k: sc[k] for k in keys}})
            parts.append(frame.with_columns([((pl.col(k) - pl.col(k).mean()) / pl.col(k).std()) for k in keys]))
        if parts:
            per_tissue[name] = pl.concat(parts)
            print(f"mouse_coupling: {name}: {per_tissue[name].height} control mice", file=sys.stderr, flush=True)
    keys = [c for c in next(iter(per_tissue.values())).columns if c.split("_")[0] in MEASURES]
    rows = []
    for a, b in combinations(sorted(per_tissue), 2):
        for k in keys:
            r, n = pair_r(per_tissue[a], per_tissue[b], k)
            if n < MIN_SHARED or not np.isfinite(r):
                continue
            z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
            measure, _, group = k.partition("_")
            rows.append(
                {
                    "tissue_a": a,
                    "tissue_b": b,
                    "measure": measure,
                    "group": group,
                    "mice": n,
                    "r": r,
                    "low": float(np.tanh(z - 1.96 * se)),
                    "high": float(np.tanh(z + 1.96 * se)),
                    "p": float(2 * stats.norm.sf(abs(z) / se)),
                }
            )
    pairs = pl.DataFrame(rows)
    write_table(pairs, out / "pairs")
    rng = ctx.rng(390_001)
    summary = []
    for (measure, group), g in pairs.group_by("measure", "group", maintain_order=True):
        key = f"{measure}_{group}"
        observed = float(g["r"].mean())
        null = []
        for _ in range(N_PERM[str(ctx.scale.name)]):
            rs = [
                pair_r(per_tissue[a], shuffled(per_tissue[b], rng), key)[0]
                for a, b in g.select("tissue_a", "tissue_b").iter_rows()
            ]
            null.append(np.nanmean(rs))
        null = np.array(null)
        summary.append(
            {
                "measure": measure,
                "group": group,
                "pairs": g.height,
                "mean_r": observed,
                "share_positive": float((g["r"] > 0).mean()),
                "null_mean": float(null.mean()),
                "null_q95": float(np.quantile(null, 0.95)),
                "p_perm": float((1 + (null >= observed).sum()) / (1 + null.size)),
            }
        )
    write_table(pl.DataFrame(summary), out / "summary")
    return []
