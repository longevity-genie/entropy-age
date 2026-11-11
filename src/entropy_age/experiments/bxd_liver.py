"""The shield's share with age and diet in a slowly renewing tissue: BXD mouse liver (roadmap RM3).

In living mouse blood the shield grows as a share of drift with age, and caloric restriction weakens it
(``figures interventions``). If that growth comes from selection among dividing stem cells, a tissue whose
cells rarely divide should show little of it. Roy et al. (GSE199979) profiled the livers of 339 BXD mice
(321 female, 45 strains, 6-33 months) on chow or a high-fat diet (HFD) on the mammalian methylation array;
the series matrix carries the beta values.

* **Aging gaps** - entropy and dispersion slopes on chow-fed females, with strain as covariates, each gene set
  against promoter probes of the same island status and young-reference level decile, with a bootstrap over
  mice.
* **Per-mouse scores** against young chow females (at most ``REFERENCE_DAYS`` old): drift ``B`` and the
  shield's share ``F``, fitted as ``~ age + HFD + strain`` (strains with at least two mice; the rest pooled).
  HFD drives hepatocytes into cycle, so selection predicts a stronger shield on HFD; chromatin, none.

Run as a job: ``entropy-age figures bxd``. Its only input is the series matrix
(``data/input/geo/GSE199979/GSE199979_series_matrix.txt.gz``) and the mouse annotation.
"""

from __future__ import annotations

import gzip
import io
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.common import Context, dummies, scaled_intervals, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.per_person import score_against
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.protection import fit_reference
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics

REFERENCE_DAYS = 270
MIN_REFERENCE = 20
N_BOOT = {"paper": 200, "quick": 30}


def read_bxd(matrix: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Samples (age in months, strain, sex, diet) and the ``cpg`` x sample beta table of the series matrix."""
    samples = _characteristics(matrix).select(
        "sample",
        (pl.col("age (days)").cast(pl.Float64) / 30.44).alias("age_months"),
        "strain",
        pl.col("Sex").alias("sex"),
        "diet",
    )
    with gzip.open(matrix, "rt") as handle:
        text = handle.read()
    start = text.index("!series_matrix_table_begin") + len("!series_matrix_table_begin\n")
    end = text.index("!series_matrix_table_end")
    betas = pl.read_csv(io.StringIO(text[start:end]), separator="\t", null_values=["", "NA", "null"])
    return samples, betas.rename({"ID_REF": "cpg"})


def probes_without_suffix(probes: pl.DataFrame) -> pl.DataFrame:
    """The annotation's ``cg..._BC21`` probe names as the 40k array writes them (``cg...``), first kept."""
    return probes.with_columns(pl.col("cpg").str.replace(r"_[A-Z]{2}\d+$", "")).unique(
        subset="cpg", keep="first", maintain_order=True
    )


def ols(y: np.ndarray, design: np.ndarray, column: int) -> tuple[float, float, float]:
    """Coefficient, SE and two-sided p of ``design[:, column]``."""
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    dof = y.size - np.linalg.matrix_rank(design)
    se = float(np.sqrt(resid @ resid / dof * np.linalg.pinv(design.T @ design)[column, column]))
    return float(coef[column]), se, float(2 * stats.t.sf(abs(coef[column] / se), dof))


def strain_dummies(strains: np.ndarray) -> np.ndarray:
    """One-hot strains with at least two mice; the others pooled as one level."""
    values, counts = np.unique(strains, return_counts=True)
    keep = set(values[counts >= 2])
    return dummies(np.array([s if s in keep else "other" for s in strains]))


def run(ctx: Context) -> list[Path]:
    matrix = ctx.paths.input / "geo" / "GSE199979" / "GSE199979_series_matrix.txt.gz"
    if not matrix.exists():
        print(f"bxd: {matrix} missing; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("bxd_liver")
    samples, betas = read_bxd(matrix)
    samples = samples.filter(pl.col("sex") == "F").sort("sample")
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    betas = betas.filter(pl.col("cpg").is_in(probes["cpg"].implode()))
    m = FeatureMatrix.from_frame(betas, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.05)
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(m.features)), pl.read_parquet(lay.gene_sets))
    members = members.filter(pl.col("group").is_in(SETS))
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in m.features])
    x, ages = m.values, samples["age_months"].to_numpy()
    chow = (samples["diet"] == "CD").to_numpy()
    young = chow & (ages <= REFERENCE_DAYS / 30.44)
    if young.sum() < MIN_REFERENCE:
        raise ValueError(
            f"bxd: {young.sum()} chow females at most {REFERENCE_DAYS} days old, fewer than {MIN_REFERENCE}"
        )
    strains = samples["strain"].to_numpy()
    print(
        f"bxd: {len(m.features)} promoter probes, {x.shape[0]} females, {int(young.sum())} in the reference",
        file=sys.stderr,
    )

    # aging gaps on chow females, strain as covariates
    level = x[young].mean(axis=0)
    strata = np.char.add(isl.astype(str), baseline_strata(level, 10).astype(str))

    def gaps(rows: np.ndarray) -> list[dict]:
        cov = strain_dummies(strains[rows])
        cov = cov[:, cov.std(axis=0) > 0] if cov.shape[1] else None
        metrics = drift_metrics(x[rows], ages[rows], m.features, covariates=cov, spearman=False)
        metrics = metrics.with_columns(stratum=pl.Series(strata))
        rows_out = []
        for metric in ("entropy_slope", "dispersion_slope"):
            table = compare_groups(metrics, members, metric=metric, stratify_by="stratum", n_permutations=0)
            sd = float(metrics[metric].std())
            rows_out += [
                {"metric": metric, "group": r["group"], "matched_diff": r["matched_diff"], "sd": sd}
                for r in table.to_dicts()
            ]
        return rows_out

    chow_rows = np.flatnonzero(chow)
    point = pl.DataFrame(gaps(chow_rows))
    rng = ctx.rng(170_000)
    boot = pl.DataFrame(
        [
            {"boot": b, **r}
            for b in range(N_BOOT[str(ctx.scale.name)])
            for r in gaps(np.sort(rng.choice(chow_rows, size=chow_rows.size, replace=True)))
        ]
    )
    aging = scaled_intervals(point, boot, ["metric", "group"]).with_columns(mice=pl.lit(int(chow.sum())))
    write_table(aging, out / "aging_gaps")

    # per-mouse scores
    sets = _sets(m.features, members)
    member = np.any([s for s, _ in sets.values()], axis=0)
    scored = np.flatnonzero(~young)
    rows = []
    all_scores = []
    for k in (0, 5):
        reference = fit_reference(x[young], isl, ~member, components=k)
        scores = pl.DataFrame(
            {"sample": samples["sample"].to_numpy()[scored], **score_against(x[scored], reference, sets)}
        )
        scores = scores.join(samples, on="sample").with_columns(components=pl.lit(k))
        all_scores.append(scores)
        hfd = (scores["diet"] == "HFD").to_numpy().astype(float)
        a = scores["age_months"].to_numpy()
        design = np.column_stack([np.ones(a.size), a / 10, hfd, strain_dummies(scores["strain"].to_numpy())])
        for group in sets:
            for measure in ("B", "F"):
                y = scores[f"{measure}_{group}"].to_numpy()
                for term, col in (("age, per 10 months", 1), ("HFD", 2)):
                    coef, se, p = ols(y, design, col)
                    rows.append(
                        {
                            "components": k,
                            "group": group,
                            "measure": measure,
                            "term": term,
                            "effect": coef,
                            "se": se,
                            "p": p,
                            "mean": float(y.mean()),
                            "mice": int(y.size),
                        }
                    )
    scores = pl.concat(all_scores)
    effects = pl.DataFrame(rows)
    write_table(scores, out / "mouse_scores")
    write_table(effects, out / "effects")
    return [_figure(ctx, aging, scores, effects)]


def _figure(ctx: Context, aging: pl.DataFrame, scores: pl.DataFrame, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    ax = axes[0]
    for j, metric in enumerate(("entropy_slope", "dispersion_slope")):
        sub = {r["group"]: r for r in aging.filter(pl.col("metric") == metric).iter_rows(named=True)}
        groups = [g for g in SETS if g in sub]
        xs = np.arange(len(groups)) + (j - 0.5) * 0.3
        g = np.array([sub[k]["gap_sd"] for k in groups])
        lo = np.array([sub[k]["ci_low"] for k in groups])
        hi = np.array([sub[k]["ci_high"] for k in groups])
        ax.errorbar(
            xs,
            g,
            yerr=[np.clip(g - lo, 0, None), np.clip(hi - g, 0, None)],
            fmt="os"[j],
            color=plots.PALETTE[j],
            label=metric.replace("_", " "),
        )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set_xticks(range(len(groups)), [SET_LABELS[k] for k in groups])
    ax.set(ylabel="Set minus matched probes (SD), 95 % CI", title="Aging, chow-fed BXD liver")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "a")
    for i, measure in enumerate(("B", "F")):
        ax = axes[1 + i]
        sel = scores.filter(pl.col("components") == 0)
        for j, diet in enumerate(("CD", "HFD")):
            sub = sel.filter(pl.col("diet") == diet)
            ax.scatter(
                sub["age_months"],
                sub[f"{measure}_housekeeping"],
                s=9,
                alpha=0.7,
                color=plots.PALETTE[j],
                label={"CD": "chow", "HFD": "high-fat"}[diet],
            )
        ax.set(
            xlabel="Age (months)",
            ylabel={"B": "B: matched drift", "F": "F: shield as a share of drift"}[measure],
            title=f"Housekeeping {measure}, liver",
        )
        ax.legend(fontsize=7)
        plots.label_panel(ax, "bc"[i])
    ax = axes[3]
    sel = effects.filter(pl.col("group") == "housekeeping")
    labels = []
    y = 0
    for (measure, term), sub in sel.group_by("measure", "term", maintain_order=True):
        for j, r in enumerate(sub.sort("components").iter_rows(named=True)):
            rel = r["effect"] / abs(r["mean"])
            ax.errorbar(
                rel,
                y + 0.2 * j,
                xerr=1.96 * r["se"] / abs(r["mean"]),
                fmt="os"[j],
                color=plots.PALETTE[2 + j],
                label=f"{r['components']} components" if y == 0 else None,
            )
        labels.append(f"{measure} · {term}")
        y += 1
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(np.arange(len(labels)) + 0.1, labels, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="Effect relative to |mean| (95 % CI)", title="Age and diet, strain held fixed")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "d")
    fig.suptitle(
        "BXD mouse liver (GSE199979): the shield's share with age and diet",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("bxd_liver"))
