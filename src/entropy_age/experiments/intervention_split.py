"""Do interventions reverse the directed part of aging or the noise? (roadmap RM27, RM12)

Idea 3 of docs/BEYOND_SHIELD.md. Perevoshchikova et al. 2024 (bioRxiv, doi 10.1101/2024.02.25.581928) report that in
mouse, caloric restriction and parabiosis move the "dynamic" part of methylation aging and not the "entropic" part.
Huang et al. (GSE279437) treated aged mice with exenatide (GLP-1 agonist, 13-14 per tissue) or rapamycin (5),
against young and aged controls, in the same 11 tissues: renewing (blood, colon, spleen) and not (heart, muscle,
brain). That also serves RM12: selection predicts an effect only where cells renew.

Per tissue and treatment, over promoter CpGs:

1. Aged controls are split at random into a selection half and a reference half (``SPLITS`` times, averaged), so
   that choosing CpGs on the aged controls does not make a treatment look rejuvenating by regression to the mean.
2. On the selection half against young controls: **directed** CpGs, the top ``TOP`` by |Welch t| of the mean
   change; **noisy** CpGs, the top ``TOP`` by the rise of log variance among those whose mean barely moves (|t| < 1).
3. **Directed reversal**: the slope, through the origin, of treated minus reference mean on young minus reference
   mean over directed CpGs (0: no reversal, 1: back to young). **Noise reversal**: the median over noisy CpGs of
   (log var reference - log var treated) / (log var reference - log var young).
4. A null from permuting treated and reference animals (``PERMUTATIONS``), with the same CpG sets; each statistic
   is reported with the null's mean subtracted and a one-sided p.

``entropy-age figures intervention_split`` (as a Slurm job).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.mouse_tissues import TURNOVER, _features, _tissues

TREATMENTS = ("Aged_exenatide", "Aged_rapamycin")
YOUNG, AGED = "Young_ctrl", "Aged_ctrl"
TOP = 0.02
SPLITS = 20
PERMUTATIONS = 200
MIN_TREATED = 5


def welch_t(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    return (a.mean(axis=0) - b.mean(axis=0)) / np.sqrt(va / a.shape[0] + vb / b.shape[0] + 1e-12)


def choose_sets(selection: np.ndarray, young: np.ndarray, top: float = TOP) -> tuple[np.ndarray, np.ndarray]:
    """Indices of directed CpGs (largest |t|) and noisy CpGs (largest log-variance rise with |t| < 1)."""
    t = welch_t(selection, young)
    k = max(10, int(top * t.size))
    directed = np.argsort(-np.abs(t))[:k]
    lv = np.log(selection.var(axis=0, ddof=1) + 1e-8) - np.log(young.var(axis=0, ddof=1) + 1e-8)
    calm = np.flatnonzero(np.abs(t) < 1)
    noisy = calm[np.argsort(-lv[calm])[:k]]
    return directed, noisy


def reversals(treated, reference, young, directed, noisy) -> tuple[float, float]:
    """Directed reversal (slope) and noise reversal (median fraction of the variance rise undone)."""
    dy = young[:, directed].mean(axis=0) - reference[:, directed].mean(axis=0)
    dt = treated[:, directed].mean(axis=0) - reference[:, directed].mean(axis=0)
    slope = float((dt * dy).sum() / (dy * dy).sum())

    def lv(x):
        return np.log(x[:, noisy].var(axis=0, ddof=1) + 1e-8)

    rise = lv(reference) - lv(young)
    ok = rise > 0
    frac = (lv(reference) - lv(treated))[ok] / rise[ok]
    return slope, float(np.median(frac)) if ok.any() else float("nan")


def split_test(x, groups, treatment, rng, splits=SPLITS, permutations=PERMUTATIONS) -> dict:
    """Both reversals for ``treatment`` against aged controls, averaged over random halvings, with a null."""
    young, treated = x[groups == YOUNG], x[groups == treatment]
    aged_idx = np.flatnonzero(groups == AGED)
    obs, null = [], []
    for _ in range(splits):
        perm = rng.permutation(aged_idx)
        half = perm.size // 2
        sel, ref = x[perm[:half]], x[perm[half:]]
        directed, noisy = choose_sets(sel, young)
        obs.append(reversals(treated, ref, young, directed, noisy))
        pool = np.vstack([treated, ref])
        for _ in range(max(1, permutations // splits)):
            p = rng.permutation(pool.shape[0])
            null.append(reversals(pool[p[: treated.shape[0]]], pool[p[treated.shape[0] :]], young, directed, noisy))
    obs, null = np.array(obs), np.array(null)
    o = np.nanmean(obs, axis=0)
    return {
        "directed_reversal": float(o[0]),
        "directed_null": float(np.nanmean(null[:, 0])),
        "directed_p": float((1 + (null[:, 0] >= o[0]).sum()) / (1 + null.shape[0])),
        "noise_reversal": float(o[1]),
        "noise_null": float(np.nanmean(null[:, 1])),
        "noise_p": float((1 + (null[:, 1] >= o[1]).sum()) / (1 + null.shape[0])),
        "treated": int(treated.shape[0]),
        "aged": int(aged_idx.size),
        "young": int(young.shape[0]),
    }


def run(ctx: Context) -> list[Path]:
    if not cohort_layout("huang", ctx.paths).betas.exists():
        print("intervention_split: Huang cohort missing (`fetch cohort huang`). Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("intervention_split")
    rng = ctx.rng(270_000)
    features, _, _ = _features(ctx.paths, ctx.rng(80_000))
    rows = []
    for name, t in _tissues(ctx.paths, features).items():
        for treatment in TREATMENTS:
            if (t.groups == treatment).sum() < MIN_TREATED or (t.groups == YOUNG).sum() < MIN_TREATED:
                continue
            r = split_test(t.x, t.groups, treatment, rng)
            rows.append({"tissue": name, "turnover": TURNOVER[name], "treatment": treatment, **r})
            print(f"intervention_split: {name} {treatment} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "reversals")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=1, size=4.4)
    for i, treatment in enumerate(TREATMENTS):
        ax = axes[i]
        t = table.filter(pl.col("treatment") == treatment).sort("turnover", "tissue")
        y = np.arange(t.height)
        ax.barh(
            y - 0.2,
            t["directed_reversal"] - t["directed_null"],
            height=0.38,
            color=plots.PALETTE[0],
            label="directed part",
        )
        ax.barh(y + 0.2, t["noise_reversal"] - t["noise_null"], height=0.38, color=plots.PALETTE[1], label="noise part")
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(y, [f"{a} ({b})" for a, b in t.select("tissue", "turnover").iter_rows()], fontsize=7.5)
        ax.invert_yaxis()
        ax.set(
            xlabel="share of aging reversed, minus permutation null",
            title=f"{treatment.replace('Aged_', '')} ({int(t['treated'].max() or 0)} mice per tissue)",
        )
        if i == 0:
            ax.legend(fontsize=7)
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "What interventions reverse: the directed or the noisy part of methylation aging (GSE279437)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("intervention_split"))
