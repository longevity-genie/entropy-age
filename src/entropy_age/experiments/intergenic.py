"""Directed demethylation outside promoters: is aging methylation all stochastic? (roadmap RM15)

Asked by SvenDP: beyond variability growing at single loci, do specific intergenic regions lose methylation
together with age, each CpG too little for a per-CpG test, but visibly in a PCA? The known candidate is the loss at
solo-WCGW CpGs (A or T on both sides, no other CpG within 35 bp) in late-replicating, CpG-poor domains, which
tracks cell divisions (Zhou et al. 2018; Endicott et al. 2022).

On the 450k open-sea CpGs outside promoters (autosomes), in Hannum blood and Jaffe cortex:

* **directed loss** - each sample's mean methylation of solo-WCGW CpGs and of the other open-sea CpGs, fitted
  against age with the cohort's technical covariates; a stochastic process moves the mean little, a directed
  one moves it steadily;
* **coordination** - the leading principal components of the open-sea CpGs (technical covariates removed,
  age kept): their correlation with age, the share of variance each explains, and whether a component's
  loadings are larger on solo-WCGW CpGs than on the others;
* **blood against cortex** - neurons do not divide, so a division-linked loss should be weak in cortex.

The composition components of the other figures are leading PCs of background CpGs; if an age-tracking
component here is the solo-WCGW axis, those corrections have been removing a biological programme as cell mix.

Run as a job after ``entropy-age fetch annotation`` (it adds the solo-WCGW tags): ``entropy-age figures intergenic``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import _hannum, _inputs_missing, _jaffe
from entropy_age.stochasticity import regress_out

N_COMPONENTS = 10
MAX_PCA_CPGS = 40_000


def directed_loss(x: np.ndarray, ages: np.ndarray, covariates: np.ndarray | None, solo: np.ndarray) -> list[dict]:
    """Per-sample mean methylation of solo-WCGW and other CpGs, and their slope on age per decade."""
    rows = []
    cols = [np.ones(ages.size), ages / 10]
    if covariates is not None:
        cov = covariates[:, covariates.std(axis=0) > 0]
        cols.append(cov)
    design = np.column_stack(cols)
    for name, mask in (("solo-WCGW", solo), ("other open sea", ~solo)):
        y = x[:, mask].mean(axis=1)
        coef, se, p = ols(y, design, 1)
        rows.append(
            {"cpgs": name, "n_cpgs": int(mask.sum()), "mean": float(y.mean()), "per_decade": coef, "se": se, "p": p}
        )
    y = x[:, solo].mean(axis=1) - x[:, ~solo].mean(axis=1)
    coef, se, p = ols(y, design, 1)
    rows.append(
        {
            "cpgs": "solo minus other",
            "n_cpgs": int(solo.sum()),
            "mean": float(y.mean()),
            "per_decade": coef,
            "se": se,
            "p": p,
        }
    )
    return rows


def components(x: np.ndarray, ages: np.ndarray, solo: np.ndarray, k: int) -> tuple[list[dict], np.ndarray]:
    """Leading PCs of standardised CpGs: share of variance, r with age, and how much more the component weighs
    solo-WCGW CpGs (mean |loading|, tested) and in which direction they move with age (signed loading)."""
    z = x - x.mean(axis=0)
    sd = z.std(axis=0)
    keep = sd > 0
    z, s_ok = z[:, keep] / sd[keep], solo[keep]
    u, sv, vt = np.linalg.svd(z, full_matrices=False)
    share = sv**2 / (sv**2).sum()
    rows = []
    for j in range(k):
        score = u[:, j] * sv[j]
        r = stats.pearsonr(score, ages)
        load = vt[j]
        sign = np.sign(r.statistic) or 1.0
        signed = load * sign  # loadings oriented so that a positive value rises with age (a loss is negative)
        t = stats.ttest_ind(np.abs(load[s_ok]), np.abs(load[~s_ok]), equal_var=False)
        rows.append(
            {
                "component": j + 1,
                "variance_share": float(share[j]),
                "r_age": float(r.statistic),
                "p_age": float(r.pvalue),
                "solo_loading": float(signed[s_ok].mean()),
                "other_loading": float(signed[~s_ok].mean()),
                "solo_weight": float(np.abs(load[s_ok]).mean()),
                "other_weight": float(np.abs(load[~s_ok]).mean()),
                "p_loading": float(t.pvalue),
            }
        )
    return rows, u[:, :k] * sv[:k]


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "intergenic"):
        return []
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k)
    if "solo_wcgw" not in probes.columns:
        print(
            "intergenic: the probe table has no solo-WCGW tags; rerun `entropy-age fetch annotation`", file=sys.stderr
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("intergenic")
    sea = probes.filter(
        (pl.col("island_relation") == "OpenSea") & ~pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"])
    )
    features = sorted(sea["cpg"].to_list())
    solo_of = dict(sea.select("cpg", "solo_wcgw").iter_rows())
    rng = ctx.rng(200_000)
    loss_rows, pc_rows, scores = [], [], {}
    for label, loader in (("blood (Hannum)", _hannum), ("cortex (Jaffe)", _jaffe)):
        p = loader(ctx.paths, features)
        solo = np.array([bool(solo_of[f]) for f in p.features])
        loss_rows += [{"tissue": label, **r} for r in directed_loss(p.x, p.ages, p.technical, solo)]
        cov = p.technical[:, p.technical.std(axis=0) > 0] if p.technical is not None else None
        resid = regress_out(p.x, cov) if cov is not None and cov.shape[1] else p.x
        pick = np.sort(rng.choice(resid.shape[1], size=min(MAX_PCA_CPGS, resid.shape[1]), replace=False))
        rows, sc = components(resid[:, pick], p.ages, solo[pick], N_COMPONENTS)
        pc_rows += [{"tissue": label, **r} for r in rows]
        scores[label] = (p.ages, resid[:, solo].mean(axis=1), resid[:, ~solo].mean(axis=1))
        print(
            f"intergenic: {label}: {len(p.features)} open-sea CpGs, {int(solo.sum())} solo-WCGW",
            file=sys.stderr,
            flush=True,
        )
    loss, pcs = pl.DataFrame(loss_rows), pl.DataFrame(pc_rows)
    write_table(loss, out / "directed_loss")
    write_table(pcs, out / "components")
    return [_figure(ctx, loss, pcs, scores)]


def _figure(ctx: Context, loss: pl.DataFrame, pcs: pl.DataFrame, scores: dict) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    for i, (label, (ages, solo, other)) in enumerate(scores.items()):
        ax = axes[i]
        ax.scatter(ages, solo * 100, s=7, color=plots.PALETTE[1], label="solo-WCGW", alpha=0.7)
        ax.scatter(ages, other * 100, s=7, color=plots.PALETTE[0], label="other open sea", alpha=0.7)
        for y, c in ((solo, plots.PALETTE[1]), (other, plots.PALETTE[0])):
            f = stats.linregress(ages, y * 100)
            g = np.array([ages.min(), ages.max()])
            ax.plot(g, f.intercept + f.slope * g, color=c, lw=1.5)
        ax.set(
            xlabel="Age (years)", ylabel="Mean methylation, % (covariates removed)", title=f"{label}: directed loss?"
        )
        ax.legend(fontsize=7)
        plots.label_panel(ax, "ab"[i])
    ax = axes[2]
    for j, t in enumerate(scores):
        sub = pcs.filter(pl.col("tissue") == t)
        ax.bar(sub["component"] + (j - 0.5) * 0.4, sub["r_age"].abs(), width=0.4, color=plots.PALETTE[j], label=t)
    ax.set(xlabel="Principal component", ylabel="|r| with age", title="Which components track age")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "c")
    ax = axes[3]
    for j, t in enumerate(scores):
        sub = pcs.filter(pl.col("tissue") == t)
        ax.bar(
            sub["component"] + (j - 0.5) * 0.4,
            sub["solo_weight"] / sub["other_weight"],
            width=0.4,
            color=plots.PALETTE[j],
            label=t,
        )
    ax.axhline(1, color=plots.MUTED, lw=0.8)
    ax.set(
        xlabel="Principal component",
        ylabel="solo-WCGW / other, mean |loading|",
        title="Is a component the solo-WCGW axis?",
    )
    ax.legend(fontsize=7)
    plots.label_panel(ax, "d")
    fig.suptitle(
        "Outside promoters: coordinated, directed methylation loss with age?",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("intergenic"))
