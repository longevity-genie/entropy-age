"""Do published clocks weight protected CpGs up to make up for their slower change? (roadmap RM8; R64)

Asked by antonkulaga (2026-11-13): an age model fits protected and unprotected CpGs alike, so the shield may not
matter to it; do existing clocks' coefficients show that the fit simply scaled up the weights of protected CpGs,
which still change with age, only more slowly?

For each published clock (``entropy-age fetch clocks``) and two tissues (Hannum blood; Jaffe's adult cortex
controls, the samples of ``figures tissues``), each clock CpG is described by its age slope (beta per decade), its
spread slope (how fast people come to differ at it: the slope on age of the absolute residual from the age fit,
per decade) and its level in the youngest quarter of donors. Protected CpGs are those in promoters of
core-essential or housekeeping genes. Per clock: the share of protected CpGs, their median age and spread slopes
and weight against the clock's other CpGs, their share of the clock's age signal (sum of ``|weight × age slope|``),
and the weight they get beyond what their slope and starting level earn: ``log |weight| ~ log |age slope| +
young level + protected``, the protected term as a ratio, with a 95 % interval from resampling the clock's CpGs.

``entropy-age figures clock_weights`` as a Slurm job (after ``fetch clocks`` and the Hannum and Jaffe cohorts).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.clocks import CLOCKS, ClockLayout
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.clock_cpgs import clock_table
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.reprocessed_jaffe import adult_controls

PROTECTED = ("core_essential", "housekeeping")
MIN_PROTECTED = 5
N_BOOT = {"paper": 1000, "quick": 100}
TISSUES = {"blood (Hannum)": "hannum", "cortex (Jaffe)": "jaffe"}


def cpg_stats(x: np.ndarray, ages: np.ndarray) -> dict[str, np.ndarray]:
    """Per column of ``x`` (samples × CpGs, NaN allowed): age slope and spread slope per decade, young level."""
    out = {"slope": [], "spread_slope": [], "young": []}
    young = ages <= np.quantile(ages, 0.25)
    for j in range(x.shape[1]):
        v = x[:, j]
        ok = ~np.isnan(v)
        a, y = ages[ok], v[ok]
        if ok.sum() < 10 or np.ptp(a) == 0:
            for k in out:
                out[k].append(np.nan)
            continue
        b = np.polyfit(a, y, 1)
        resid = np.abs(y - np.polyval(b, a))
        out["slope"].append(b[0] * 10)
        out["spread_slope"].append(np.polyfit(a, resid, 1)[0] * 10)
        out["young"].append(float(np.mean(v[young & ok])) if (young & ok).any() else np.nan)
    return {k: np.asarray(v, dtype=float) for k, v in out.items()}


def _beyond_slope(absw: np.ndarray, abss: np.ndarray, young: np.ndarray, prot: np.ndarray) -> tuple[float, float]:
    """The protected coefficient of ``log|w| ~ log|slope| + young + protected`` and its OLS p."""
    X = np.column_stack([np.ones(absw.size), np.log(abss), young, prot.astype(float)])
    y = np.log(absw)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = y.size - X.shape[1]
    r = y - X @ b
    cov = (r @ r / dof) * np.linalg.pinv(X.T @ X)
    t = b[3] / np.sqrt(cov[3, 3])
    return float(b[3]), float(2 * stats.t.sf(abs(t), dof))


def clock_summary(g: pl.DataFrame, n_boot: int, rng: np.random.Generator) -> dict[str, float] | None:
    """One clock in one tissue: columns ``weight``, ``slope``, ``spread_slope``, ``young``, ``protected``."""
    g = g.drop_nulls(["slope", "young"]).filter(pl.col("weight") != 0)
    prot = g["protected"].to_numpy()
    if prot.sum() < MIN_PROTECTED or (~prot).sum() < MIN_PROTECTED:
        return None
    absw, abss = np.abs(g["weight"].to_numpy()), np.abs(g["slope"].to_numpy())
    young, spread = g["young"].to_numpy(), g["spread_slope"].to_numpy()
    keep = abss > 0
    coef, p = _beyond_slope(absw[keep], abss[keep], young[keep], prot[keep])
    draws = []
    idx = np.flatnonzero(keep)
    for _ in range(n_boot):
        s = rng.choice(idx, size=idx.size, replace=True)
        if prot[s].sum() < 2 or (~prot[s]).sum() < 2:
            continue
        draws.append(_beyond_slope(absw[s], abss[s], young[s], prot[s])[0])
    signal = absw * abss
    return {
        "cpgs": int(g.height),
        "protected_cpgs": int(prot.sum()),
        "share_cpgs": float(prot.mean()),
        "share_signal": float(signal[prot].sum() / signal.sum()),
        "slope_ratio": float(np.median(abss[prot]) / np.median(abss[~prot])),
        "spread_ratio": float(np.median(spread[prot]) / np.median(spread[~prot])),
        "weight_ratio": float(np.median(absw[prot]) / np.median(absw[~prot])),
        "young_protected": float(np.median(young[prot])),
        "young_other": float(np.median(young[~prot])),
        "weight_beyond_slope": float(np.exp(coef)),
        "ci_low": float(np.exp(np.quantile(draws, 0.025))),
        "ci_high": float(np.exp(np.quantile(draws, 0.975))),
        "p": p,
    }


def _samples(name: str, paths) -> pl.DataFrame:
    s = pl.read_parquet(cohort_layout(name, paths).samples)
    if name == "jaffe":
        return adult_controls(s)
    return s.drop_nulls("age").sort("age", "sample")


def run(ctx: Context) -> list[Path]:
    lay = ClockLayout(ctx.paths)
    if any(not lay.file(c).exists() for c in CLOCKS):
        print("clock_weights: run `entropy-age fetch clocks` first; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("clock_weights")
    ann = AnnotationLayout(ctx.paths)
    mem = link_membership(pl.read_parquet(ann.links_450k), pl.read_parquet(ann.gene_sets))
    protected = set(mem.filter(pl.col("group").is_in(PROTECTED))["feature"].to_list())
    clocks = clock_table(ctx.paths).with_columns(protected=pl.col("cpg").is_in(list(protected)))
    cpgs = sorted(set(clocks["cpg"].to_list()))
    rng = ctx.rng(640_000)
    rows, per_cpg = [], []
    for tissue, name in TISSUES.items():
        s = _samples(name, ctx.paths)
        frame = (
            pl.scan_parquet(cohort_layout(name, ctx.paths).betas)
            .select("cpg", *s["sample"].to_list())
            .filter(pl.col("cpg").is_in(cpgs))
            .collect()
        )
        st = cpg_stats(frame.drop("cpg").to_numpy().T.astype(float), s["age"].to_numpy().astype(float))
        stats_frame = pl.DataFrame({"cpg": frame["cpg"], **st})
        t = clocks.join(stats_frame, on="cpg")
        per_cpg.append(t.with_columns(tissue=pl.lit(tissue)))
        for (clock, kind), g in t.group_by("clock", "kind", maintain_order=True):
            r = clock_summary(g, N_BOOT[str(ctx.scale.name)], rng)
            if r is not None:
                rows.append({"tissue": tissue, "clock": clock, "kind": kind, "donors": s.height, **r})
        print(f"clock_weights: {tissue}: {s.height} donors, {frame.height} clock CpGs", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "summary")
    write_table(pl.concat(per_cpg), out / "per_cpg")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    ax = axes[0]
    clocks = table["clock"].unique(maintain_order=True).to_list()
    blood = pl.DataFrame({"clock": clocks}).join(
        table.filter(pl.col("tissue") == "blood (Hannum)"), on="clock", how="left"
    )
    y = np.arange(len(clocks))
    ax.barh(y - 0.2, blood["share_cpgs"], height=0.38, color=plots.MUTED, label="share of the clock's CpGs")
    ax.barh(y + 0.2, blood["share_signal"], height=0.38, color=plots.PALETTE[0], label="share of its age signal")
    ax.set_yticks(y, clocks, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="protected CpGs (share)", title="blood")
    ax.legend(fontsize=7, loc="lower right")
    plots.label_panel(ax, "a")
    ax = axes[1]
    for i, tissue in enumerate(TISSUES):
        t = pl.DataFrame({"clock": clocks}).join(table.filter(pl.col("tissue") == tissue), on="clock", how="left")
        y = np.arange(len(clocks)) + 0.2 * (i - 0.5)
        ax.errorbar(
            t["weight_beyond_slope"],
            y,
            xerr=[t["weight_beyond_slope"] - t["ci_low"], t["ci_high"] - t["weight_beyond_slope"]],
            fmt="o",
            ms=4,
            lw=1,
            color=plots.PALETTE[i],
            label=tissue,
        )
    ax.set_xscale("log")
    ticks = [0.25, 0.5, 1, 2, 4]
    ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.axvline(1, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(len(clocks)), [""] * len(clocks))
    ax.invert_yaxis()
    ax.set(xlabel="weight of protected CpGs, slope and level held (×)", title="extra weight")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "b")
    fig.suptitle("Published clocks' weights on protected CpGs", x=0.02, y=1.03, ha="left", fontsize=11)
    fig.subplots_adjust(wspace=0.08)
    return plots.save(fig, ctx.figure_path("clock_weights"))
