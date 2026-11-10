"""The per-person protection score, validated on cohorts already used here (roadmap RM6).

``entropy_age.protection`` gives each sample a score ``P`` per gene set: how far its set CpGs sit from
a reference state, minus how far CpGs of the same island status and reference level sit, in the same
sample (``docs/PER_PERSON_SCORE.md``). Three of the design's validation steps need no new data:

1. **It reproduces the cohort result.** Hannum blood and Jaffe cortex, the young (20-35 y) of each
   cohort as its reference; every older sample is scored, and ``P`` is regressed on age (per decade)
   with the cohort's technical covariates (plate, sex, ancestry, in blood the recruitment site, in
   cortex the neuron fraction). The cohort test (``figures tissues``) found the dispersion slope of
   set CpGs below that of matched CpGs in both tissues, more in blood; ``P`` should fall with age
   in both, more in blood.
2. **A sample's own error against the spread between people.** Each ``P`` carries an interval from
   200 resamples of the set's genes. The share of the between-person variance at a given age (the
   residual variance of the age fit) not explained by that sampling error is the score's
   reliability; technical replicates, the design's real test, are not in these cohorts.
3. **Turnover inside one person.** Lunnon's donors gave blood and four brain regions. Here the
   reference is the other donors' mean in the same tissue (leave one out): everyone is 70-100 y, so
   ``d`` measures how far a donor has moved from the old population's centre, the per-person form of
   the cohort's dispersion test. The relative gap ``F = P / B`` takes out that blood and brain differ in
   total dispersion; selection predicts blood's ``F`` below the same donor's brain ``F`` for most
   donors, not only on average.

Step 4 (division count in culture) is the culture test, ``figures passage`` and ``passage_sturm``,
which scores each culture sample the same way. Every step is run without and with five composition
components, fitted on the reference's background CpGs and projected into each sample (``Reference``).

Run as a job (about 10 GB): ``entropy-age figures per_person``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import (
    LUNNON_BRAIN,
    SET_LABELS,
    SETS,
    Prepared,
    _features,
    _hannum,
    _inputs_missing,
    _jaffe,
    _lunnon,
)
from entropy_age.protection import (
    deviations,
    entropy_gain,
    fit_reference,
    gene_counts,
    matched_scores,
    resampled_gaps,
)

REFERENCE_AGES = (20.0, 35.0)
MIN_REFERENCE = 20
COMPONENTS = (0, 5)
N_RESAMPLES = {"paper": 200, "quick": 40}


def _sets(features: list[str], membership: pl.DataFrame) -> dict[str, tuple[np.ndarray, list[np.ndarray]]]:
    out = {}
    for group in SETS:
        rows = membership.filter(pl.col("group") == group)
        members = {g: sel["feature"].to_list() for (g,), sel in rows.group_by("gene", maintain_order=True)}
        in_set, genes = gene_counts(features, members)
        if in_set.sum() >= 2:
            out[group] = (in_set, genes)
    return out


def _islands(features: list[str], probes: pl.DataFrame) -> np.ndarray:
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    return np.array([bool(island.get(f, False)) for f in features])


def score_against(
    x: np.ndarray, reference, sets: dict, *, n_resamples: int = 0, rng: np.random.Generator | None = None
) -> dict[str, np.ndarray]:
    """Every set's ``P``, ``B``, ``F`` and entropy gap ``E`` for the rows of ``x``, and ``P``'s resampling SE."""
    adjusted = reference.adjust(x)
    d = deviations(adjusted, reference.level)
    h = entropy_gain(np.clip(adjusted, 0.0, 1.0), reference.level)
    scores = {}
    for group, (in_set, genes) in sets.items():
        p, b = matched_scores(d, in_set, reference.strata)
        e, _ = matched_scores(h, in_set, reference.strata)
        scores |= {f"P_{group}": p, f"B_{group}": b, f"F_{group}": p / b, f"E_{group}": e}
        if n_resamples:
            draws = resampled_gaps(d, genes, reference.strata, n_resamples, rng)
            scores[f"se_{group}"] = draws.std(axis=0, ddof=1)
    return scores


# ----------------------------------------------------------------------------- step 1 and 2
def lifespan_scores(p: Prepared, probes, membership, components: int, n_resamples: int, rng) -> pl.DataFrame:
    lo, hi = REFERENCE_AGES
    young = (p.ages >= lo) & (p.ages <= hi)
    if young.sum() < MIN_REFERENCE:
        raise ValueError(f"{p.label}: {young.sum()} samples aged {lo:g}-{hi:g}, fewer than {MIN_REFERENCE}")
    sets = _sets(p.features, membership)
    member = np.any([s for s, _ in sets.values()], axis=0)
    reference = fit_reference(p.x[young], _islands(p.features, probes), ~member, components=components)
    scored = ~young
    scores = score_against(p.x[scored], reference, sets, n_resamples=n_resamples, rng=rng)
    return pl.DataFrame(
        {
            "tissue": p.label,
            "components": components,
            "row": np.flatnonzero(scored),
            "age": p.ages[scored],
            "reference_samples": int(young.sum()),
            **scores,
        }
    )


def age_fits(scores: pl.DataFrame, covariates: np.ndarray | None) -> list[dict]:
    """OLS of each score on age per decade with the sample's technical covariates; reliability of ``P``."""
    rows = scores["row"].to_numpy()
    age = scores["age"].to_numpy() / 10.0
    columns = [np.ones(age.size), age]
    if covariates is not None:
        cov = covariates[rows]
        columns.append(cov[:, cov.std(axis=0) > 0])
    design = np.column_stack(columns)
    rank = np.linalg.matrix_rank(design)
    out = []
    for group in SETS:
        if f"P_{group}" not in scores.columns:
            continue
        for measure in ("P", "F", "E", "B"):
            y = scores[f"{measure}_{group}"].to_numpy()
            coef, *_ = np.linalg.lstsq(design, y, rcond=None)
            resid = y - design @ coef
            dof = y.size - rank
            sigma2 = resid @ resid / dof
            se = float(np.sqrt(sigma2 * np.linalg.pinv(design.T @ design)[1, 1]))
            row = {
                "tissue": scores["tissue"][0],
                "components": int(scores["components"][0]),
                "group": group,
                "measure": measure,
                "samples": y.size,
                "slope_per_decade": float(coef[1]),
                "se": se,
                "p": float(2 * stats.t.sf(abs(coef[1] / se), dof)),
                "mean": float(y.mean()),
                "between_sd": float(np.sqrt(sigma2)),
            }
            if measure == "P" and f"se_{group}" in scores.columns:
                noise = float(np.mean(scores[f"se_{group}"].to_numpy() ** 2))
                row |= {"within_sd": float(np.sqrt(noise)), "reliability": float(1 - noise / sigma2)}
            out.append(row)
    return out


# ----------------------------------------------------------------------------- step 3
def leave_one_out_scores(p: Prepared, probes, membership, components: int) -> pl.DataFrame:
    """Each donor scored against the other donors of the same tissue."""
    sets = _sets(p.features, membership)
    member = np.any([s for s, _ in sets.values()], axis=0)
    islands = _islands(p.features, probes)
    rows = []
    for i in range(p.x.shape[0]):
        others = np.arange(p.x.shape[0]) != i
        reference = fit_reference(p.x[others], islands, ~member, components=components)
        scores = score_against(p.x[i : i + 1], reference, sets)
        rows.append({"donor": p.donors[i], "age": float(p.ages[i]), **{k: float(v[0]) for k, v in scores.items()}})
    return pl.DataFrame(rows).with_columns(tissue=pl.lit(p.label), components=pl.lit(components))


def paired_turnover(scores: pl.DataFrame) -> pl.DataFrame:
    """Per region and set: blood minus region in each donor; share of donors with blood below, sign test."""
    out = []
    for (components,), sel in scores.group_by("components", maintain_order=True):
        blood = sel.filter(pl.col("tissue") == "whole blood")
        for region in LUNNON_BRAIN:
            other = sel.filter(pl.col("tissue") == region)
            both = blood.join(other, on="donor", suffix="_r")
            if both.height == 0:
                continue
            for group in SETS:
                for measure in ("F", "P", "E"):
                    column = f"{measure}_{group}"
                    if column not in both.columns:
                        continue
                    diff = (both[column] - both[f"{column}_r"]).to_numpy()
                    below = int((diff < 0).sum())
                    out.append(
                        {
                            "components": components,
                            "region": region,
                            "group": group,
                            "measure": measure,
                            "donors": diff.size,
                            "blood_below": below,
                            "share_below": below / diff.size,
                            "median_diff": float(np.median(diff)),
                            "mean_diff": float(diff.mean()),
                            "p_sign": float(stats.binomtest(below, diff.size).pvalue),
                        }
                    )
    return pl.DataFrame(out)


# ----------------------------------------------------------------------------- run
def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "per_person"):
        return []
    plots.apply_style()
    out = ctx.out_dir("per_person")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    lifespan = {"blood": _hannum(ctx.paths, features), "cortex": _jaffe(ctx.paths, features)}
    life_scores, fits = [], []
    for key, p in lifespan.items():
        for k in COMPONENTS:
            s = lifespan_scores(p, probes, membership, k, n_resamples, ctx.rng(130_000 + 10 * k + len(fits)))
            life_scores.append(s.drop("row").with_columns(cohort=pl.lit(key)))
            fits += age_fits(s, p.technical)
            print(f"per_person: {p.label}, {k} components: {s.height} samples scored", file=sys.stderr)
    life = pl.concat(life_scores, how="diagonal")
    fits = pl.DataFrame(fits)
    write_table(life, out / "lifespan_scores")
    write_table(fits, out / "age_fits")

    lunnon = _lunnon(ctx.paths, features)
    same = pl.concat(
        [leave_one_out_scores(p, probes, membership, k) for p in lunnon.values() for k in COMPONENTS], how="diagonal"
    )
    paired = paired_turnover(same)
    write_table(same, out / "same_donor_scores")
    write_table(paired, out / "same_donor_turnover")
    return [_figure(ctx, life, fits, same, paired)]


def _figure(ctx: Context, life: pl.DataFrame, fits: pl.DataFrame, same: pl.DataFrame, paired: pl.DataFrame) -> Path:
    fig, axes = plots.panels(6, ncols=3, size=3.6)
    group = "housekeeping"
    for i, (cohort, label) in enumerate((("blood", "blood (Hannum)"), ("cortex", "cortex (Jaffe)"))):
        ax = axes[i]
        sel = life.filter((pl.col("cohort") == cohort) & (pl.col("components") == 0))
        age, y = sel["age"].to_numpy(), sel[f"P_{group}"].to_numpy()
        ax.errorbar(age, y, yerr=1.96 * sel[f"se_{group}"].to_numpy(), fmt="none", ecolor=plots.MUTED, alpha=0.25)
        ax.scatter(age, y, s=7, color=plots.PALETTE[i])
        fit = stats.linregress(age, y)
        grid = np.array([age.min(), age.max()])
        ax.plot(grid, fit.intercept + fit.slope * grid, color="black", lw=1.2)
        ax.axhline(0, color=plots.MUTED, lw=0.8, ls=":")
        ax.set(
            xlabel="Age (years)",
            ylabel="P, housekeeping (set − matched d)",
            title=f"{label}: one score per person",
        )
        plots.label_panel(ax, "ab"[i])

    ax = axes[2]
    y = 0.0
    ticks, labels = [], []
    for g in SETS:
        for tissue in ("blood (Hannum)", "cortex (Jaffe)"):
            for j, k in enumerate(COMPONENTS):
                r = fits.filter(
                    (pl.col("group") == g)
                    & (pl.col("tissue") == tissue)
                    & (pl.col("components") == k)
                    & (pl.col("measure") == "P")
                )
                if r.height == 0:
                    continue
                r = r.row(0, named=True)
                ax.errorbar(
                    r["slope_per_decade"],
                    y + 0.25 * j,
                    xerr=1.96 * r["se"],
                    fmt="os"[j],
                    color=plots.PALETTE[0 if tissue.startswith("blood") else 1],
                    ms=4,
                    label=f"{k} components" if (g, tissue) == (SETS[0], "blood (Hannum)") else None,
                )
            ticks.append(y + 0.12)
            labels.append(f"{SET_LABELS[g]} · {tissue.split(' ')[0]}")
            y += 1
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(ticks, labels, fontsize=7)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.set(xlabel="Change in P per decade (95 % CI)", title="Step 1: the score falls with age")
    plots.label_panel(ax, "c")

    ax = axes[3]
    rel = fits.filter((pl.col("measure") == "P") & (pl.col("components") == 0)).sort("tissue", "group")
    xs = np.arange(rel.height)
    ax.bar(xs - 0.2, rel["between_sd"], width=0.4, color=plots.PALETTE[2], label="between people (age-adjusted)")
    ax.bar(xs + 0.2, rel["within_sd"], width=0.4, color=plots.MUTED, label="one sample (gene resampling)")
    ax.set_xticks(
        xs,
        [f"{SET_LABELS[g]}\n{t.split(' ')[0]}" for g, t in zip(rel["group"], rel["tissue"], strict=True)],
        fontsize=7,
        rotation=30,
        ha="right",
    )
    ax.legend(fontsize=7)
    ax.set(ylabel="SD of P", title="Step 2: people differ more than the error")
    plots.label_panel(ax, "d")

    ax = axes[4]
    sel = same.filter(pl.col("components") == 0)
    blood = sel.filter(pl.col("tissue") == "whole blood").select("donor", pl.col(f"F_{group}").alias("blood"))
    for j, region in enumerate(LUNNON_BRAIN):
        other = sel.filter(pl.col("tissue") == region).select("donor", pl.col(f"F_{group}").alias("region"))
        both = blood.join(other, on="donor")
        diff = (both["blood"] - both["region"]).to_numpy()
        ax.scatter(j + np.random.default_rng(j).uniform(-0.15, 0.15, diff.size), diff, s=8, color=plots.PALETTE[3])
        ax.hlines(np.median(diff), j - 0.3, j + 0.3, color="black", lw=1.5)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls=":")
    ax.set_xticks(range(len(LUNNON_BRAIN)), [r.replace(" ", "\n", 1) for r in LUNNON_BRAIN], fontsize=7)
    ax.set(ylabel="F blood − F region, same donor", title="Step 3: housekeeping, each donor")
    plots.label_panel(ax, "e")

    ax = axes[5]
    share = paired.filter((pl.col("measure") == "F") & (pl.col("components") == 0))
    for j, g in enumerate(SETS):
        sub = share.filter(pl.col("group") == g)
        ax.bar(
            np.arange(sub.height) + (j - 1) * 0.27,
            sub["share_below"],
            width=0.27,
            color=plots.PALETTE[j],
            label=SET_LABELS[g],
        )
    ax.axhline(0.5, color=plots.MUTED, lw=0.8, ls=":")
    ax.set_xticks(range(len(LUNNON_BRAIN)), [r.replace(" ", "\n", 1) for r in LUNNON_BRAIN], fontsize=7)
    ax.set(ylim=(0, 1), ylabel="Share of donors with blood's F below", title="Step 3: donors, by set")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "f")

    fig.suptitle(
        "A per-person protection score (450k promoters; young reference in a and b, other donors in e and f)",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("per_person"))
