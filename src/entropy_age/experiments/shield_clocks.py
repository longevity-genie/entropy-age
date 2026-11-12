"""Do people whose clocks run fast have a weaker shield? (roadmap RM11)

The per-person shield score ``P`` (``per_person``) is how far a person's essential-gene CpGs have drifted beyond
look-alike CpGs since youth: ``P < 0`` is a shield, and the less negative, the weaker. Epigenetic clocks read a
person's age from methylation; their age acceleration is the clock's reading beyond the person's age.

Per cohort (Hannum blood, Jaffe cortex) and person, this computes ``P`` as ``per_person`` does, and each published
clock's linear predictor (``data.clocks``; the intercept and any output transform are monotone and drop out of a
residual). Acceleration is the predictor's residual after age and the cohort's technical covariates; ``P`` is
residualised alike. Their partial correlation answers the question: a positive one means people with faster
clocks have a weaker shield. Clocks include Tong et al.'s stochastic clocks (StocH, StocP, StocZ), built from
age-associated noise.

``entropy-age figures shield_clocks`` (after ``fetch clocks`` and the cohorts of ``figures per_person``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.clocks import CLOCKS, ClockLayout, read_clock
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import COMPONENTS, lifespan_scores
from entropy_age.experiments.tissue_contrast import SETS, _features, _hannum, _inputs_missing, _jaffe, _matrix

SAMPLE_ORDER = {"blood": "hannum", "cortex": "jaffe"}


def linear_score(x: np.ndarray, features: list[str], weights: pl.DataFrame) -> tuple[np.ndarray, int]:
    """Samples' ``x @ w`` over the clock CpGs present (a missing value takes the CpG's mean); CpGs used."""
    index = {f: i for i, f in enumerate(features)}
    used = [(index[c], w) for c, w in weights.iter_rows() if c in index]
    if not used:
        return np.full(x.shape[0], np.nan), 0
    cols, w = np.array([u[0] for u in used]), np.array([u[1] for u in used])
    sub = x[:, cols]
    sub = np.where(np.isnan(sub), np.nanmean(sub, axis=0), sub)
    return sub @ w, len(used)


def residual(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def partial_r(a: np.ndarray, b: np.ndarray, design: np.ndarray) -> tuple[float, float, int]:
    """Correlation of ``a`` and ``b`` after both are residualised on ``design``; two-sided p; n."""
    ra, rb = residual(a, design), residual(b, design)
    r, _ = stats.pearsonr(ra, rb)
    n, k = a.size, np.linalg.matrix_rank(design)
    t = r * np.sqrt((n - k - 1) / max(1e-12, 1 - r**2))
    return float(r), float(2 * stats.t.sf(abs(t), n - k - 1)), n


def _cohort_samples(cohort: str, paths) -> list[str]:
    name = SAMPLE_ORDER[cohort]
    s = pl.read_parquet(cohort_layout(name, paths).samples)
    if name == "hannum":
        return s.drop_nulls("age").sort("age", "sample")["sample"].to_list()
    return (
        s.filter((pl.col("group") == "Control") & ~pl.col("drop_sample") & (pl.col("age") >= 20))
        .sort("best_qc", descending=True)
        .unique(subset="donor", keep="first", maintain_order=True)
        .drop_nulls(["age", "neuron_fraction"])
        .sort("age", "sample")["sample"]
        .to_list()
    )


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "shield_clocks"):
        return []
    clock_layout = ClockLayout(ctx.paths)
    clocks = {name: read_clock(clock_layout.file(name)) for name in CLOCKS if clock_layout.file(name).exists()}
    if not clocks:
        print("shield_clocks: run `entropy-age fetch clocks` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("shield_clocks")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    rows = []
    for cohort, loader in (("blood", _hannum), ("cortex", _jaffe)):
        p = loader(ctx.paths, features)
        samples = _cohort_samples(cohort, ctx.paths)
        if len(samples) != p.x.shape[0]:
            raise ValueError(f"shield_clocks: {cohort}: {len(samples)} samples against {p.x.shape[0]} rows")
        cm = _matrix(SAMPLE_ORDER[cohort], ctx.paths, samples, clock_cpgs)
        x_clock = cm.values
        for k in COMPONENTS:
            scores = lifespan_scores(p, probes, membership, k, 0, ctx.rng(130_000 + k))
            r_idx = scores["row"].to_numpy()
            age = p.ages[r_idx]
            design = np.column_stack([np.ones(r_idx.size), age, p.covariates[r_idx]])
            design = design[:, np.r_[0, 1, 2 + np.flatnonzero(design[:, 2:].std(axis=0) > 0)]]
            for clock, weights in clocks.items():
                score, used = linear_score(x_clock[r_idx], cm.features, weights)
                if used < 10:
                    continue
                for group in SETS:
                    for measure in ("P", "F"):
                        r, pval, n = partial_r(scores[f"{measure}_{group}"].to_numpy(), score, design)
                        rows.append(
                            {
                                "cohort": cohort,
                                "components": k,
                                "clock": clock,
                                "kind": CLOCKS[clock][2],
                                "cpgs_used": used,
                                "cpgs_total": weights.height,
                                "group": group,
                                "measure": measure,
                                "r": r,
                                "p": pval,
                                "n": n,
                            }
                        )
            print(f"shield_clocks: {cohort}, {k} components done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "correlations")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    sel = table.filter((pl.col("measure") == "P") & (pl.col("group") == "housekeeping"))
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    for i, cohort in enumerate(("blood", "cortex")):
        ax = axes[i]
        sub = sel.filter(pl.col("cohort") == cohort)
        clocks = sorted(sub["clock"].unique().to_list(), key=lambda c: (CLOCKS[c][2], c))
        for j, k in enumerate(COMPONENTS):
            s = sub.filter(pl.col("components") == k)
            ys = [clocks.index(c) + 0.25 * j for c in s["clock"].to_list()]
            ax.scatter(s["r"], ys, marker="os"[j], color=plots.PALETTE[j], s=22, label=f"{k} components")
            for xv, yv, pv in zip(s["r"], ys, s["p"], strict=True):
                if pv < 0.05:
                    ax.annotate("*", (xv, yv), xytext=(4, -3), textcoords="offset points", fontsize=9)
        ax.set_yticks([j + 0.12 for j in range(len(clocks))], [f"{c} ({CLOCKS[c][2]})" for c in clocks], fontsize=7)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.legend(fontsize=7, frameon=False)
        ax.set(
            xlabel="Partial r: shield score P with clock acceleration",
            title={"blood": "Hannum blood", "cortex": "Jaffe cortex"}[cohort],
        )
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Do people with faster clocks have a weaker shield? (> 0: yes; housekeeping P; * p < 0.05)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("shield_clocks"))
