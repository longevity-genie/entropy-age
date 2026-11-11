"""Six tissues of mice with numeric ages: the shield and its growth by turnover (roadmap RM2).

Thompson et al. 2018 (GSE120132, ``entropy-age fetch rrbs-tissues``) profiled blood, liver, muscle, lung,
adipose and kidney of C57BL/6J and BALB/cByJ mice of both sexes at about 2, 10 and 20 months by RRBS. The
promoter CpGs of each tissue (covered at least 10 times in 90 % of its samples) are tested as in
``figures interventions``:

* **aging gaps** - entropy and dispersion slopes with strain and sex as covariates, each gene set against
  promoter CpGs of the same CpG-island status and young (about 2 months) level decile, bootstrap over mice;
* **per-mouse scores** - each mouse against the young mice of its own strain and tissue, and the shield's
  share ``F`` and drift ``B`` fitted as ``~ age + sex + strain``, without and with five composition
  components of the young reference.

Selection predicts the shield and its growth with age to follow turnover: largest in blood, smallest in
muscle. Chromatin predicts no such order. Run as a job: ``entropy-age fetch rrbs-tissues`` (about 6.4 GB),
then ``entropy-age figures rrbs_tissues``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.rrbs import ThompsonLayout
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, scaled_intervals, write_table
from entropy_age.experiments.interventions import _sets, membership
from entropy_age.experiments.per_person import score_against
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.protection import fit_reference
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics

#: tissue -> turnover class, in plotting order
TURNOVER = {
    "Blood": "renewing",
    "Lung": "slow",
    "Kidney": "slow",
    "Liver": "slow",
    "Adipose": "slow",
    "Muscle": "post-mitotic",
}
YOUNG_MONTHS = 3.0
MIN_YOUNG = 5
N_BOOT = {"paper": 100, "quick": 20}


def covariates(frame: pl.DataFrame) -> np.ndarray:
    return np.column_stack([dummies(frame["strain"].to_numpy()), dummies(frame["sex"].to_numpy())])


def tissue_gaps(x, ages, cov, features, strata, members) -> list[dict]:
    """Set minus matched CpGs (SD of the metric) for the entropy and dispersion slopes."""
    keep = cov.std(axis=0) > 0
    metrics = drift_metrics(x, ages, features, covariates=cov[:, keep] if keep.any() else None, spearman=False)
    metrics = metrics.with_columns(stratum=pl.Series(strata))
    out = []
    for metric in ("entropy_slope", "dispersion_slope"):
        table = compare_groups(metrics, members, metric=metric, stratify_by="stratum", n_permutations=0)
        sd = float(metrics[metric].std())
        out += [
            {"metric": metric, "group": r["group"], "matched_diff": r["matched_diff"], "sd": sd}
            for r in table.to_dicts()
        ]
    return out


def run(ctx: Context) -> list[Path]:
    lay = ThompsonLayout(ctx.paths)
    if not lay.samples.exists():
        print("rrbs_tissues: run `entropy-age fetch rrbs-tissues` first; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("mouse_rrbs_tissues")
    samples_all = pl.read_parquet(lay.samples)
    gap_rows, effect_rows, sizes = [], [], []
    rng = ctx.rng(180_000)
    for tissue, turnover in TURNOVER.items():
        if not lay.betas(tissue).exists():
            raise ValueError(f"rrbs_tissues: {lay.betas(tissue)} missing; rerun `entropy-age fetch rrbs-tissues`")
        samples = samples_all.filter(pl.col("tissue") == tissue).sort("sample")
        m = FeatureMatrix.from_frame(
            pl.read_parquet(lay.betas(tissue)), feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.1
        )
        cpgs = pl.read_parquet(lay.cpgs(tissue))
        members = membership(ctx.paths, m.features, cpgs)
        island = dict(cpgs.select("cpg", "cpg_island").iter_rows())
        isl = np.array([bool(island[f]) for f in m.features])
        x, ages = m.values, samples["age_months"].to_numpy()
        young = ages <= YOUNG_MONTHS
        strata = np.char.add(isl.astype(str), baseline_strata(x[young].mean(axis=0), 10).astype(str))
        cov = covariates(samples)
        sizes.append({"tissue": tissue, "mice": x.shape[0], "young": int(young.sum()), "cpgs": len(m.features)})
        print(
            f"rrbs_tissues: {tissue}: {x.shape[0]} mice, {len(m.features)} promoter CpGs", file=sys.stderr, flush=True
        )

        point = tissue_gaps(x, ages, cov, m.features, strata, members)
        boot = []
        for b in range(N_BOOT[str(ctx.scale.name)]):
            rows = np.sort(rng.choice(x.shape[0], size=x.shape[0], replace=True))
            boot += [{"boot": b, **r} for r in tissue_gaps(x[rows], ages[rows], cov[rows], m.features, strata, members)]
        boot = pl.DataFrame(boot)
        gap_rows.append(
            scaled_intervals(pl.DataFrame(point), boot, ["metric", "group"]).with_columns(
                tissue=pl.lit(tissue), turnover=pl.lit(turnover)
            )
        )

        sets = _sets(m.features, members)
        member = np.any([s for s, _ in sets.values()], axis=0)
        for k in (0, 5):
            frames = []
            for (strain,), sel in samples.with_row_index("row").group_by("strain", maintain_order=True):
                rows = sel["row"].to_numpy()
                ref = rows[ages[rows] <= YOUNG_MONTHS]
                if ref.size < max(MIN_YOUNG, k + 2):
                    print(f"rrbs_tissues: {tissue}, {strain}: {ref.size} young mice; strain left out", file=sys.stderr)
                    continue
                reference = fit_reference(x[ref], isl, ~member, components=k)
                older = np.setdiff1d(rows, ref)
                frames.append(
                    pl.DataFrame(
                        {"sample": samples["sample"].to_numpy()[older], **score_against(x[older], reference, sets)}
                    )
                )
            scores = pl.concat(frames).join(samples, on="sample")
            design = np.column_stack([np.ones(scores.height), scores["age_months"].to_numpy() / 10, covariates(scores)])
            design = design[:, np.r_[True, True, design[:, 2:].std(axis=0) > 0]]
            for group in sets:
                for measure in ("B", "F"):
                    y = scores[f"{measure}_{group}"].to_numpy()
                    coef, se, p = ols(y, design, 1)
                    effect_rows.append(
                        {
                            "tissue": tissue,
                            "turnover": turnover,
                            "components": k,
                            "group": group,
                            "measure": measure,
                            "per_10_months": coef,
                            "se": se,
                            "p": p,
                            "mean": float(y.mean()),
                            "mice": int(y.size),
                        }
                    )
    gaps, effects = pl.concat(gap_rows), pl.DataFrame(effect_rows)
    write_table(pl.DataFrame(sizes), out / "tissues")
    write_table(gaps, out / "aging_gaps")
    write_table(effects, out / "age_effects")
    return [_figure(ctx, gaps, effects)]


def _figure(ctx: Context, gaps: pl.DataFrame, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=4.0)
    tissues = list(TURNOVER)
    for i, metric in enumerate(("entropy_slope", "dispersion_slope")):
        ax = axes[i]
        for j, group in enumerate(SETS):
            sub = {
                r["tissue"]: r
                for r in gaps.filter((pl.col("metric") == metric) & (pl.col("group") == group)).iter_rows(named=True)
            }
            xs = np.arange(len(tissues)) + (j - 1) * 0.25
            g = np.array([sub[t]["gap_sd"] if t in sub else np.nan for t in tissues])
            lo = np.array([sub[t]["ci_low"] if t in sub else np.nan for t in tissues])
            hi = np.array([sub[t]["ci_high"] if t in sub else np.nan for t in tissues])
            ax.errorbar(
                xs,
                g,
                yerr=[np.clip(g - lo, 0, None), np.clip(hi - g, 0, None)],
                fmt="o",
                ms=4,
                color=plots.PALETTE[j],
                label=SET_LABELS[group],
            )
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set_xticks(range(len(tissues)), [f"{t}\n({TURNOVER[t]})" for t in tissues], fontsize=7)
        ax.set(ylabel="Set minus matched CpGs (SD), 95 % CI", title=metric.replace("_", " ").capitalize())
        ax.legend(fontsize=6)
        plots.label_panel(ax, "ab"[i])
    ax = axes[2]
    for j, k in enumerate((0, 5)):
        sub = {
            r["tissue"]: r
            for r in effects.filter(
                (pl.col("group") == "housekeeping") & (pl.col("measure") == "F") & (pl.col("components") == k)
            ).iter_rows(named=True)
        }
        xs = np.arange(len(tissues)) + (j - 0.5) * 0.3
        v = np.array([sub[t]["per_10_months"] for t in tissues])
        e = np.array([1.96 * sub[t]["se"] for t in tissues])
        ax.errorbar(xs, v, yerr=e, fmt="os"[j], color=plots.PALETTE[3 + j], label=f"{k} components")
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set_xticks(range(len(tissues)), tissues, fontsize=7)
    ax.set(ylabel="Change in F per 10 months (95 % CI)", title="Housekeeping shield share with age")
    ax.legend(fontsize=6)
    plots.label_panel(ax, "c")
    fig.suptitle(
        "Six mouse tissues with numeric ages (GSE120132, RRBS)", x=0.02, ha="left", fontsize=11, fontweight="bold"
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("mouse_rrbs_tissues"))
