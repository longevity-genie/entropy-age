"""Clocks of our own: does a clock need unprotected CpGs? (roadmap RM11)

Meyer and Schumacher built a clock from noise alone. If the shield is real, CpGs at essential promoters carry
little of that noise, and a clock trained on them alone should read age worse than one trained on matched
look-alike CpGs. Per cohort (Hannum blood, Jaffe cortex), three pools of equal size are drawn (``pools``):

* ``shielded``: CpGs of the housekeeping and core-essential promoters;
* ``look-alike``: background CpGs of the same island x young-level strata, as many per stratum (where background
  CpGs are scarcer, both pools take that number);
* ``random``: as many background CpGs drawn at random.

Each pool trains an elastic-net clock (``clock.train_clock``, fast search) with 5-fold outer cross-validation;
predictions are out of fold. Accuracy (Pearson r, median absolute error) per pool and draw, and the correlation
between pools' acceleration (does a shielded-CpG clock run fast in the same people as a look-alike clock?).

``entropy-age figures own_clocks`` (after the cohorts of ``figures per_person``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.model_selection import KFold

from entropy_age import plots
from entropy_age.clock import ClockConfig, evaluate, train_clock
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import _features, _hannum, _inputs_missing, _jaffe
from entropy_age.stochasticity import baseline_strata

DRAWS = 3
FOLDS = 5
YOUNG = (20.0, 35.0)
SHIELDED = ("housekeeping", "core_essential")


def pools(features: list[str], members: set[str], strata: np.ndarray, rng) -> dict[str, np.ndarray]:
    """Column indices of a shielded pool, a stratum-matched look-alike pool and a random pool of equal size.

    Per stratum both the shielded and the look-alike pools take as many CpGs as the scarcer of the two has there,
    so the pools match stratum by stratum; the random pool has their total size.
    """
    is_set = np.array([f in members for f in features])
    background = np.flatnonzero(~is_set)
    set_cols = np.flatnonzero(is_set)
    shielded, look = [], []
    for s in np.unique(strata[set_cols]):
        mine = set_cols[strata[set_cols] == s]
        candidates = background[strata[background] == s]
        k = min(mine.size, candidates.size)
        if k == 0:
            continue
        shielded.append(rng.choice(mine, size=k, replace=False))
        look.append(rng.choice(candidates, size=k, replace=False))
    shielded, look = np.sort(np.concatenate(shielded)), np.sort(np.concatenate(look))
    random = np.sort(rng.choice(background, size=shielded.size, replace=False))
    return {"shielded": shielded, "look-alike": look, "random": random}


def out_of_fold(x: np.ndarray, ages: np.ndarray, seed: int) -> np.ndarray:
    pred = np.empty(ages.size)
    for train, test in KFold(FOLDS, shuffle=True, random_state=seed).split(x):
        pred[test] = train_clock(x[train], ages[train], config=ClockConfig.fast()).predict(x[test])
    return pred


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "own_clocks"):
        return []
    plots.apply_style()
    out = ctx.out_dir("own_clocks")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    members = set(membership.filter(pl.col("group").is_in(list(SHIELDED)))["feature"].to_list())
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    metrics, accel = [], []
    for cohort, loader in (("blood", _hannum), ("cortex", _jaffe)):
        p = loader(ctx.paths, features)
        x = np.where(np.isnan(p.x), np.nanmean(p.x, axis=0), p.x)
        young = (p.ages >= YOUNG[0]) & (p.ages <= YOUNG[1])
        level = x[young].mean(axis=0) if young.sum() >= 10 else x.mean(axis=0)
        strata = np.char.add(
            np.array([str(bool(island.get(f, False))) for f in p.features]), baseline_strata(level, 10).astype(str)
        )
        for draw in range(DRAWS):
            rng = ctx.rng(160_000 + draw)
            preds = {}
            for name, cols in pools(p.features, members, strata, rng).items():
                pred = out_of_fold(x[:, cols], p.ages, seed=draw)
                preds[name] = pred
                m = evaluate(p.ages, pred)
                metrics.append({"cohort": cohort, "draw": draw, "pool": name, "cpgs": int(cols.size), **m.as_dict()})
            design = np.column_stack([np.ones(p.ages.size), p.ages])
            res = {k: v - design @ np.linalg.lstsq(design, v, rcond=None)[0] for k, v in preds.items()}
            names = list(res)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    accel.append(
                        {
                            "cohort": cohort,
                            "draw": draw,
                            "a": names[i],
                            "b": names[j],
                            "r": float(np.corrcoef(res[names[i]], res[names[j]])[0, 1]),
                        }
                    )
            print(f"own_clocks: {cohort}, draw {draw} done", file=sys.stderr, flush=True)
    metrics, accel = pl.DataFrame(metrics), pl.DataFrame(accel)
    write_table(metrics, out / "accuracy")
    write_table(accel, out / "acceleration_overlap")
    return [_figure(ctx, metrics)]


def _figure(ctx: Context, metrics: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.0)
    order = ["shielded", "look-alike", "random"]
    for i, (col, label) in enumerate(
        (("pearson_r", "Pearson r, predicted vs age"), ("median_abs_error", "Median absolute error (years)"))
    ):
        ax = axes[i]
        for c, cohort in enumerate(("blood", "cortex")):
            sub = metrics.filter(pl.col("cohort") == cohort)
            for k, pool in enumerate(order):
                v = sub.filter(pl.col("pool") == pool)[col].to_numpy()
                ax.scatter(
                    np.full(v.size, k + 0.15 * (2 * c - 1)),
                    v,
                    color=plots.PALETTE[c],
                    s=20,
                    label=cohort if k == 0 else None,
                )
        ax.set_xticks(range(len(order)), order)
        ax.set(ylabel=label, title=label.split(",")[0])
        if i == 0:
            ax.legend(fontsize=7, frameon=False)
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Clocks from shielded, look-alike and random CpGs of equal number (out-of-fold, 3 draws)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("own_clocks"))
