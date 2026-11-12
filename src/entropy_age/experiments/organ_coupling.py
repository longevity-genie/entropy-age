"""Do all organs of one person age at the same rate? (roadmap RM26)

Idea 1 of docs/BEYOND_SHIELD.md, asked by antonkulaga (2026-11-12) for questions beyond the shield. Where one
donor gave several tissues, each tissue's age acceleration is measured per sample and the accelerations of two
tissues are correlated over their shared donors. A clear positive correlation means a person-wide pace (a blood
test could speak for the brain); none means each organ keeps its own.

- Lunnon (GSE59685): whole blood and four brain regions of the same donors. Blood was drawn before death, so
  each sample keeps its own age. Covariates within each tissue: age, sex, Alzheimer's status (it may age several
  tissues at once, a shared cause rather than a shared pace).
- GTEx: up to 8 tissues of the same donors (methylation of 9 tissues); covariates age, sex.

Measures per sample, each residualised on the covariates within its tissue:

- clock acceleration: the linear predictor of published clocks (``CLOCKS_USED``);
- clock-free drift: ``B``, the drift of look-alike promoter CpGs against the rest of the tissue's samples
  (leave-one-out, as in ``health_states``), and the shield's share ``F``.

Then for each pair of tissues and measure: Pearson r of the residuals over shared donors, with a Fisher interval.
A shared batch could couple two tissues processed together: plate and slide differ by tissue in both cohorts.

``entropy-age figures organ_coupling`` (as a Slurm job).
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.shield_clocks import linear_score, residual
from entropy_age.experiments.tissue_contrast import _features, _matrix

CLOCKS_USED = ("Horvath 2013", "Hannum", "PhenoAge", "Cortical", "DunedinPACE", "Zhang 2019")
GROUP = "housekeeping"
MIN_SHARED = 20


def coupling(a: pl.DataFrame, b: pl.DataFrame, measure: str) -> dict | None:
    """Pearson r of two tissues' residual ``measure`` over shared donors, Fisher 95 % interval, p, n."""
    j = a.select("donor", measure).join(b.select("donor", measure), on="donor", suffix="_b").drop_nulls()
    n = j.height
    if n < MIN_SHARED:
        return None
    r, p = stats.pearsonr(j[measure], j[f"{measure}_b"])
    z, half = np.arctanh(r), 1.96 / np.sqrt(n - 3)
    return {"r": float(r), "low": float(np.tanh(z - half)), "high": float(np.tanh(z + half)), "p": float(p), "n": n}


def tissue_measures(
    x_clock: np.ndarray,
    clock_features: list[str],
    clocks: dict[str, pl.DataFrame],
    x_set: np.ndarray,
    island: np.ndarray,
    sets: dict,
    design: np.ndarray,
) -> dict[str, np.ndarray]:
    """Per sample: each clock's and the drift scores' residuals on ``design``."""
    out = {}
    for name, weights in clocks.items():
        score, used = linear_score(x_clock, clock_features, weights)
        if used:
            out[name] = residual(score, design)
    scores = loo_scores(x_set, island, sets)
    for measure in ("B", "F"):
        key = f"{measure}_{GROUP}"
        if key in scores:
            out[f"drift {measure}"] = residual(scores[key], design)
    return out


def _design(t: pl.DataFrame, categorical: list[str]) -> np.ndarray:
    parts = [np.ones((t.height, 1)), t["age"].to_numpy().astype(float)[:, None]]
    parts += [dummies(t[c].to_numpy()) for c in categorical if t[c].n_unique() > 1]
    return np.column_stack(parts)


def _cohort(
    name: str, samples: pl.DataFrame, categorical: list[str], ctx, features, clock_cpgs, probes, membership, clocks
) -> pl.DataFrame:
    frames = []
    for (tissue,), t in samples.group_by("tissue", maintain_order=True):
        t = t.sort("sample")
        if t.height < MIN_SHARED:
            continue
        m_set = _matrix(name.lower(), ctx.paths, t["sample"].to_list(), features)
        m_clock = _matrix(name.lower(), ctx.paths, t["sample"].to_list(), clock_cpgs)
        island = _islands(m_set.features, probes)
        sets = _sets(m_set.features, membership)
        res = tissue_measures(
            m_clock.values, m_clock.features, clocks, m_set.values, island, sets, _design(t, categorical)
        )
        frames.append(pl.DataFrame({"cohort": name, "tissue": tissue, "donor": t["donor"], **res}))
        print(f"organ_coupling: {name} {tissue}: {t.height} samples", file=sys.stderr, flush=True)
    return pl.concat(frames, how="diagonal")


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("organ_coupling")
    lay = ClockLayout(ctx.paths)
    clocks = {c: read_clock(lay.file(c)) for c in CLOCKS_USED if lay.file(c).exists()}
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    lunnon = (
        pl.read_parquet(cohort_layout("lunnon", ctx.paths).samples)
        .drop_nulls("age")
        .filter(pl.col("ad_status").is_in(["AD", "C"]))
    )
    gtex = pl.read_parquet(cohort_layout("gtex", ctx.paths).samples)
    per_sample = pl.concat(
        [
            _cohort("Lunnon", lunnon, ["sex", "ad_status"], ctx, features, clock_cpgs, probes, membership, clocks),
            _cohort("GTEx", gtex, ["sex"], ctx, features, clock_cpgs, probes, membership, clocks),
        ],
        how="diagonal",
    )
    measures = [c for c in per_sample.columns if c not in ("cohort", "tissue", "donor")]
    rows = []
    for (cohort,), c in per_sample.group_by("cohort", maintain_order=True):
        tissues = c["tissue"].unique(maintain_order=True).to_list()
        for ta, tb in combinations(tissues, 2):
            a, b = c.filter(pl.col("tissue") == ta), c.filter(pl.col("tissue") == tb)
            for measure in measures:
                r = coupling(a, b, measure)
                if r:
                    rows.append({"cohort": cohort, "tissue_a": ta, "tissue_b": tb, "measure": measure, **r})
    table = pl.DataFrame(rows)
    write_table(table, out / "couplings")
    write_table(per_sample, out / "per_sample")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    """Per measure, the distribution of tissue-pair correlations, Lunnon and GTEx side by side."""
    measures = table["measure"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(2, ncols=2, size=4.8)
    for i, cohort in enumerate(("Lunnon", "GTEx")):
        ax = axes[i]
        sub = table.filter(pl.col("cohort") == cohort)
        for y, measure in enumerate(measures):
            m = sub.filter(pl.col("measure") == measure)
            if m.is_empty():
                continue
            jitter = np.linspace(-0.25, 0.25, m.height) if m.height > 1 else np.zeros(1)
            ax.scatter(m["r"], y + jitter, s=10, color=plots.PALETTE[i], alpha=0.7)
            ax.plot([m["r"].median()] * 2, [y - 0.35, y + 0.35], color=plots.MUTED, lw=1.5)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(range(len(measures)), measures if i == 0 else [""] * len(measures), fontsize=8)
        ax.invert_yaxis()
        ax.set(xlabel="r between two tissues of the same donors", title=f"{cohort}: one dot per tissue pair")
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Do organs of one person age together? Acceleration correlated across tissues (bar: median)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("organ_coupling"))
