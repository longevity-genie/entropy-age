"""Does the shield predict death? The Vitality 90+ nonagenarians (roadmap RM28, first cohort)

Asked by antonkulaga (2026-11-12): cohorts with age and later health. GSE68194 (Tampere Vitality 90+, 450k
array, PBMC) gives 115 people born around 1920, alive or dead about 2.5 years after the blood draw, with the follow-
up in days, sex, methylation batch, BMI, a cognitive score (MMSE) and four flow-cytometry cell proportions. All
are nonagenarians, so age is near constant and needs no covariate.

Per person: ``B``, ``F`` and ``S`` of the gene sets against everyone else (leave-one-out, ``health_states``), and
the linear predictor of published clocks (``CLOCKS_USED``), each standardised. Each enters a Cox proportional-
hazards model alone with sex and batch (``cox``, fitted here by Newton steps on the partial likelihood, Breslow
ties), then with the CD8+CD28- share added (senescent T cells, the strongest composition effect in old blood).
Hazard ratios are per SD. MMSE is tested the same way by linear regression. 39 deaths: power is for large effects
only.

``entropy-age figures survival`` (as a Slurm job).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.shield_clocks import linear_score
from entropy_age.experiments.syndromes import matrix_path, read_betas
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE68194"
CLOCKS_USED = ("Horvath 2013", "Hannum", "PhenoAge", "DunedinPACE", "Zhang 2019", "StocH", "StocP")
GROUPS = ("core_essential", "housekeeping")
FIELDS = {
    "sex (1=woman; 0=man)": "woman",
    "methylation batch": "batch",
    "status (0=alive; 1=died)": "died",
    "follow-up time in days": "days",
    "body mass index (bmi)": "bmi",
    "mini mental state examination (mmse)": "mmse",
    "cd8+cd28- cell proportion": "cd8_cd28neg",
}


def cox(time: np.ndarray, event: np.ndarray, x: np.ndarray, *, steps: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Cox proportional-hazards coefficients and their SEs (Breslow ties), by Newton steps from zero."""
    order = np.argsort(-time, kind="stable")
    t, d, z = time[order], event[order].astype(bool), x[order]
    beta = np.zeros(z.shape[1])
    info = np.eye(z.shape[1])
    for _ in range(steps):
        w = np.exp(z @ beta)
        s0 = np.cumsum(w)
        s1 = np.cumsum(w[:, None] * z, axis=0)
        s2 = np.cumsum(w[:, None, None] * z[:, :, None] * z[:, None, :], axis=0)
        # risk sets end at the last of tied times: use the cumulative sums at each tie group's last row
        last = np.r_[np.flatnonzero(t[1:] != t[:-1]), t.size - 1]
        idx = last[np.searchsorted(last, np.arange(t.size))]
        s0, s1, s2 = s0[idx], s1[idx], s2[idx]
        mean = s1 / s0[:, None]
        grad = (z - mean)[d].sum(axis=0)
        info = (s2 / s0[:, None, None] - mean[:, :, None] * mean[:, None, :])[d].sum(axis=0)
        step = np.linalg.solve(info, grad)
        beta = beta + step
        if np.abs(step).max() < 1e-9:
            break
    return beta, np.sqrt(np.diag(np.linalg.inv(info)))


def read_people(matrix: Path) -> pl.DataFrame:
    chars = _characteristics(matrix)
    return chars.select("sample", *(pl.col(k).cast(pl.Float64, strict=False).alias(v) for k, v in FIELDS.items()))


def _z(v: np.ndarray) -> np.ndarray:
    return (v - v.mean()) / v.std(ddof=1)


def run(ctx: Context) -> list[Path]:
    matrix = matrix_path(ctx.paths, GSE)
    if not matrix.exists():
        print(f"survival: {matrix} missing (fetch {GSE}). Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("survival_shield")
    people = read_people(matrix).drop_nulls(["died", "days", "woman", "batch"])
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    lay = ClockLayout(ctx.paths)
    clocks = {c: read_clock(lay.file(c)) for c in CLOCKS_USED if lay.file(c).exists()}
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    samples = people["sample"].to_list()
    x = read_betas(matrix, features, samples)
    ok = np.isnan(x).mean(axis=0) == 0
    kept = [f for f, k in zip(features, ok, strict=True) if k]
    scores = loo_scores(x[:, ok], _islands(kept, probes), _sets(kept, membership))
    measures = {f"{m} {g}": scores[f"{m}_{g}"] for g in GROUPS for m in ("F", "B", "S") if f"{m}_{g}" in scores}
    xc = read_betas(matrix, clock_cpgs, samples)
    for name, weights in clocks.items():
        score, used = linear_score(xc, clock_cpgs, weights)
        if used:
            measures[name] = score
    print(f"survival: {people.height} people, {int(people['died'].sum())} deaths", file=sys.stderr, flush=True)
    time, event = people["days"].to_numpy(), people["died"].to_numpy()
    base = np.column_stack([people["woman"].to_numpy(), people["batch"].to_numpy()])
    cd8 = people["cd8_cd28neg"].to_numpy()
    rows = []
    for name, v in measures.items():
        z = _z(v)
        for adjust in ("sex, batch", "+ CD8+CD28-"):
            cov = base if adjust == "sex, batch" else np.column_stack([base, np.nan_to_num(cd8, nan=np.nanmean(cd8))])
            beta, se = cox(time, event, np.column_stack([z, cov]))
            rows.append(
                {
                    "measure": name,
                    "outcome": "death",
                    "adjustment": adjust,
                    "estimate": float(np.exp(beta[0])),
                    "low": float(np.exp(beta[0] - 1.96 * se[0])),
                    "high": float(np.exp(beta[0] + 1.96 * se[0])),
                    "p": float(2 * stats.norm.sf(abs(beta[0] / se[0]))),
                    "n": int(z.size),
                    "events": int(event.sum()),
                }
            )
        mmse = people["mmse"].to_numpy()
        has = ~np.isnan(mmse)
        coef, se_m, p = ols(mmse[has], np.column_stack([np.ones(has.sum()), z[has], base[has]]), 1)
        rows.append(
            {
                "measure": name,
                "outcome": "MMSE",
                "adjustment": "sex, batch",
                "estimate": coef,
                "low": coef - 1.96 * se_m,
                "high": coef + 1.96 * se_m,
                "p": p,
                "n": int(has.sum()),
                "events": None,
            }
        )
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    write_table(people.with_columns(**{k: pl.Series(v) for k, v in measures.items()}), out / "per_person")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    measures = table["measure"].unique(maintain_order=True).to_list()
    for i, outcome in enumerate(("death", "MMSE")):
        ax = axes[i]
        for y, m in enumerate(measures):
            sub = table.filter((pl.col("measure") == m) & (pl.col("outcome") == outcome))
            for j, r in enumerate(sub.iter_rows(named=True)):
                ax.errorbar(
                    r["estimate"],
                    y + 0.2 * j,
                    xerr=[[r["estimate"] - r["low"]], [r["high"] - r["estimate"]]],
                    fmt="os"[j],
                    ms=4,
                    color=plots.PALETTE[j],
                )
        ax.axvline(1.0 if outcome == "death" else 0.0, color=plots.MUTED, lw=0.8)
        if outcome == "death":
            ax.set_xscale("log")
        ax.set_yticks(range(len(measures)), measures if i == 0 else [""] * len(measures), fontsize=7.5)
        ax.invert_yaxis()
        ax.set(
            xlabel="hazard ratio per SD (95 % CI)" if outcome == "death" else "MMSE points per SD (95 % CI)",
            title="Death within ~2.5 years" if outcome == "death" else "Cognition (MMSE)",
        )
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Vitality 90+ (GSE68194): shield, drift and clocks against death and cognition; squares: + CD8+CD28-",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("survival_shield"))
