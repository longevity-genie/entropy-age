"""Prediction 2 along a turnover gradient, and lifespan interventions, in the same mice (GSE279437).

Huang et al. profiled 11 tissues of the same male mice on the mammalian methylation array: young
and aged controls, and aged mice treated with exenatide (a GLP-1 receptor agonist) or rapamycin.
If selection among cells protects essential genes from entropy, the protection should grow with
cell turnover: largest in white blood cells, colon and spleen, smaller in liver, kidney and
adipose, smallest in muscle, heart and brain. Chromatin protection alone predicts no such order.

Measures, all on promoter CpGs of autosomes, each set against background CpGs of the same CpG-island
status and young-control level decile (``compare_groups``), in SD of the metric:

* **aging drift** - aged minus young controls (``drift_metrics`` with a 0/1 age): entropy gain,
  change in spread, and absolute change of the mean;
* **dispersion among aged mice** - how far aged animals differ from each other, with treatment
  groups as covariates (the Lunnon measure, here with an independent young baseline);
* **per-group dispersion** for the interventions. Selection predicts that slowing proliferation
  (rapamycin) weakens the protection in renewing tissue. With 5 rapamycin mice this is exploratory.

GEO gives young/aged labels but no numeric ages, so there are no age slopes. Animals are resampled
within their treatment group for the bootstrap, the same animals in every tissue.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.common import Context, dummies, rescale, write_table
from entropy_age.experiments.tissue_contrast import _STATE, SET_LABELS, SETS, _interval, _matched_sd, _run_parallel
from entropy_age.stochasticity import baseline_strata, dispersion_levels, drift_metrics

#: tissue -> turnover class, in plotting order
TURNOVER: dict[str, str] = {
    "white blood cells": "renewing",
    "colon": "renewing",
    "spleen": "renewing",
    "liver": "slow",
    "kidney": "slow",
    "adipose": "slow",
    "skeletal muscle": "post-mitotic",
    "heart": "post-mitotic",
    "frontal cortex": "post-mitotic",
    "hippocampus": "post-mitotic",
    "hypothalamus": "post-mitotic",
}
CLASSES = ("renewing", "slow", "post-mitotic")
AGING_METRICS = ("entropy_slope", "dispersion_slope", "abs_slope", "kl_slope")
INTERVENTION_GROUPS = ("Aged_ctrl", "Aged_exenatide", "Aged_rapamycin")
MIN_GROUP = 4  # fewer animals than this give no per-group dispersion
BACKGROUND = 60_000
MAX_MISSING = 0.05


@dataclass(frozen=True, slots=True)
class MouseTissue:
    name: str
    x: np.ndarray
    groups: np.ndarray  # treatment group per row, e.g. "Aged_rapamycin"
    donors: np.ndarray
    features: list[str]


# ----------------------------------------------------------------------------- inputs
def _features(paths, rng: np.random.Generator) -> tuple[list[str], pl.DataFrame, pl.DataFrame]:
    lay = MouseAnnotationLayout(paths)
    probes = pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    membership = mouse_membership(probes, pl.read_parquet(lay.gene_sets)).filter(pl.col("group").is_in(SETS))
    members = set(membership["feature"].to_list())
    rest = sorted(set(probes["cpg"].to_list()) - members)
    background = sorted(rng.choice(rest, size=min(BACKGROUND, len(rest)), replace=False).tolist())
    return sorted(members | set(background)), probes, membership


def _tissues(paths, features: list[str]) -> dict[str, MouseTissue]:
    lay = cohort_layout("huang", paths)
    samples = pl.read_parquet(lay.samples)
    frame = pl.scan_parquet(lay.betas).filter(pl.col("cpg").is_in(features)).collect()
    out = {}
    for tissue in TURNOVER:
        s = samples.filter(pl.col("tissue") == tissue).sort("group", "donor")
        if s.height == 0:
            continue
        m = FeatureMatrix.from_frame(
            frame.select("cpg", *s["sample"].to_list()),
            feature_col="cpg",
            samples=s["sample"].to_list(),
            max_missing=MAX_MISSING,
        )
        out[tissue] = MouseTissue(tissue, m.values, s["group"].to_numpy(), s["donor"].to_numpy(), m.features)
    return out


# ----------------------------------------------------------------------------- metrics
def _strata(metrics: pl.DataFrame, probes: pl.DataFrame, baseline: np.ndarray) -> pl.DataFrame:
    """Replace ``baseline`` by the young-control mean and add the island x decile stratum."""
    metrics = metrics.with_columns(baseline=pl.Series(baseline))
    metrics = metrics.join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")
    decile = baseline_strata(baseline, 10)
    return metrics.with_columns(stratum=pl.format("{}|{}", pl.col("cpg_island").fill_null(False), pl.Series(decile)))


def tissue_measures(t: MouseTissue, rows: np.ndarray, probes: pl.DataFrame) -> dict[str, tuple[pl.DataFrame, tuple]]:
    """Every measure of one tissue on the animals in ``rows``: name -> (metrics, metric columns)."""
    x, g = t.x[rows], t.groups[rows]
    young = g == "Young_ctrl"
    if not young.any():
        return {}
    ref = x[young].mean(axis=0)
    out: dict[str, tuple[pl.DataFrame, tuple]] = {}
    pair = young | (g == "Aged_ctrl")
    if (g == "Aged_ctrl").sum() >= 2 and young.sum() >= 2:
        aged = (g[pair] == "Aged_ctrl").astype(float)
        drift = drift_metrics(x[pair], aged, t.features, spearman=False)
        out["aging"] = (_strata(drift, probes, ref), AGING_METRICS)
    old = np.array([s.startswith("Aged") for s in g])
    if old.sum() > 4:
        cov = dummies(g[old])
        disp = dispersion_levels(x[old], t.features, covariates=cov if cov.shape[1] else None)
        out["aged dispersion"] = (_strata(disp, probes, ref), ("dispersion",))
    for group in INTERVENTION_GROUPS:
        sel = g == group
        if sel.sum() >= MIN_GROUP:
            out[f"dispersion {group}"] = (_strata(dispersion_levels(x[sel], t.features), probes, ref), ("dispersion",))
    return out


def _rows(t: MouseTissue, measure: str, metrics: pl.DataFrame, names, membership, permutations: int, seed: int):
    table = _matched_sd(metrics, membership, names, permutations, seed)
    return table.with_columns(tissue=pl.lit(t.name), turnover=pl.lit(TURNOVER[t.name]), measure=pl.lit(measure))


# ----------------------------------------------------------------------------- bootstrap
def _boot(seed: int) -> list[dict]:
    """One resample; ``_run_parallel`` installs the state into ``tissue_contrast._STATE``."""
    state = _STATE
    tissues, probes, membership, animals = state["tissues"], state["probes"], state["membership"], state["animals"]
    rng = np.random.default_rng(seed)
    picked = np.concatenate(
        [rng.choice(group, size=group.size, replace=True) for group in animals.values()]  # within treatment group
    )
    out = []
    for t in tissues.values():
        index: dict[str, int] = {d: i for i, d in enumerate(t.donors)}
        rows = np.array([index[d] for d in picked if d in index])
        for measure, (metrics, names) in tissue_measures(t, rows, probes).items():
            table = _rows(t, measure, metrics, names, membership, 0, seed)
            out += [
                {"boot": seed, **r}
                for r in table.select("tissue", "measure", "metric", "group", "matched_diff").to_dicts()
            ]
    return out


# ----------------------------------------------------------------------------- run
def run(ctx: Context) -> list[Path]:
    lay = cohort_layout("huang", ctx.paths)
    needed = [lay.betas, lay.samples, MouseAnnotationLayout(ctx.paths).probes]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        print(
            f"mouse: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation`, "
            "`entropy-age fetch mouse-annotation` and `entropy-age fetch cohort huang` first. Nothing was run.",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("mouse_tissues")
    workers = int(os.environ.get("ENTROPY_AGE_WORKERS", "1"))
    default = 200 if ctx.scale.name == "paper" else 100
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", default))
    if n_boot != default:
        print(f"mouse: ENTROPY_AGE_BOOT={n_boot} replaces the default {default} resamples", file=sys.stderr)

    features, probes, membership = _features(ctx.paths, ctx.rng(130_000))
    tissues = _tissues(ctx.paths, features)
    write_table(
        pl.DataFrame(
            [
                {"tissue": t.name, "turnover": TURNOVER[t.name], "cpgs": len(t.features)}
                | {g: int((t.groups == g).sum()) for g in sorted(set(t.groups.tolist()))}
                for t in tissues.values()
            ],
            infer_schema_length=None,
        ).fill_null(0),
        out / "tissues",
    )
    write_table(
        membership.group_by("group").agg(pl.len().alias("cpgs"), pl.col("gene").n_unique().alias("genes")),
        out / "set_sizes",
    )

    point = []
    for t in tissues.values():
        for measure, (metrics, names) in tissue_measures(t, np.arange(t.x.shape[0]), probes).items():
            point.append(_rows(t, measure, metrics, names, membership, 1000, 3))
    point = pl.concat(point, how="diagonal")
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))
    write_table(point, out / "point_estimates")

    groups: dict[str, set[str]] = {}
    for t in tissues.values():
        for donor, group in zip(t.donors.tolist(), t.groups.tolist(), strict=True):
            groups.setdefault(group, set()).add(donor)
    animals = {g: np.array(sorted(d)) for g, d in sorted(groups.items())}
    state = {"tissues": tissues, "probes": probes, "membership": membership, "animals": animals}
    boot = pl.DataFrame(_run_parallel(_boot, [5000 + b for b in range(n_boot)], workers, state))
    boot = rescale(boot, point, ["tissue", "measure", "metric"])
    boot.write_parquet(out / "bootstrap.parquet")

    summary = summarise(point, boot)
    write_table(summary, out / "contrasts")
    return [_figure(ctx, point, boot, summary)]


def _mean_by_class(frame: pl.DataFrame, value: str) -> pl.DataFrame:
    return frame.group_by("boot", "measure", "metric", "group", "turnover").agg(pl.col(value).mean())


def summarise(point: pl.DataFrame, boot: pl.DataFrame) -> pl.DataFrame:
    """Turnover-class contrasts and intervention effects, each with a bootstrap interval.

    * ``renewing - post-mitotic`` and ``slow - post-mitotic``: the class means of a set's matched
      difference (negative = more protection in the faster-renewing class);
    * ``<treatment> - Aged_ctrl``: the change of a set's matched dispersion difference under the
      treatment, per turnover class (positive for an essential set = weaker protection).
    """
    with_class = boot.with_columns(turnover=pl.col("tissue").replace_strict(TURNOVER))
    by_class = _mean_by_class(with_class, "matched_diff_sd")
    rows = []
    for (measure, metric, group), sel in by_class.group_by("measure", "metric", "group"):
        if measure.startswith("dispersion Aged_"):
            continue
        wide = sel.pivot(on="turnover", index="boot", values="matched_diff_sd")
        for faster in ("renewing", "slow"):
            if faster not in wide.columns or "post-mitotic" not in wide.columns:
                continue
            diff = (wide[faster] - wide["post-mitotic"]).drop_nulls().to_numpy()
            if diff.size < 10:
                continue
            lo, hi, p = _interval(diff)
            rows.append(
                {
                    "contrast": f"{faster} - post-mitotic",
                    "measure": measure,
                    "metric": metric,
                    "group": group,
                    "difference": float(np.median(diff)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                }
            )
    treated = with_class.filter(pl.col("measure").str.starts_with("dispersion Aged_"))
    treated = treated.with_columns(arm=pl.col("measure").str.strip_prefix("dispersion "))
    for (group, turnover), sel in treated.group_by("group", "turnover"):
        wide = (
            sel.group_by("boot", "arm")
            .agg(pl.col("matched_diff_sd").mean())
            .pivot(on="arm", index="boot", values="matched_diff_sd")
        )
        for arm in INTERVENTION_GROUPS[1:]:
            if arm not in wide.columns or "Aged_ctrl" not in wide.columns:
                continue
            diff = (wide[arm] - wide["Aged_ctrl"]).drop_nulls().to_numpy()
            if diff.size < 10:
                continue
            lo, hi, p = _interval(diff)
            rows.append(
                {
                    "contrast": f"{arm} - Aged_ctrl ({turnover})",
                    "measure": "dispersion",
                    "metric": "dispersion",
                    "group": group,
                    "difference": float(np.median(diff)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                }
            )
    if not rows:
        return pl.DataFrame(
            schema={
                "contrast": pl.String,
                "measure": pl.String,
                "metric": pl.String,
                "group": pl.String,
                "difference": pl.Float64,
                "ci_low": pl.Float64,
                "ci_high": pl.Float64,
                "p_bootstrap": pl.Float64,
            }
        )
    return pl.DataFrame(rows).sort("contrast", "measure", "metric", "group")


# ----------------------------------------------------------------------------- figure
def _figure(ctx: Context, point: pl.DataFrame, boot: pl.DataFrame, summary: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=5.6)
    colors = {g: plots.PALETTE[i] for i, g in enumerate(SETS)}
    tissues = [t for t in TURNOVER if t in set(point["tissue"].to_list())]

    def by_tissue(ax, measure: str, metric: str, title: str) -> None:
        for j, group in enumerate(SETS):
            for i, tissue in enumerate(tissues):
                est = point.filter(
                    (pl.col("tissue") == tissue)
                    & (pl.col("measure") == measure)
                    & (pl.col("metric") == metric)
                    & (pl.col("group") == group)
                )
                if est.height == 0:
                    continue
                v = boot.filter(
                    (pl.col("tissue") == tissue)
                    & (pl.col("measure") == measure)
                    & (pl.col("metric") == metric)
                    & (pl.col("group") == group)
                )["matched_diff_sd"].to_numpy()
                xx = i + (j - 1) * 0.22
                if v.size:
                    lo, hi = np.quantile(v, [0.025, 0.975])
                    ax.plot([xx, xx], [lo, hi], color=colors[group], lw=1.1)
                ax.scatter(
                    [xx], est["matched_diff_sd"], color=colors[group], s=14, label=SET_LABELS[group] if i == 0 else None
                )
        for boundary in (2.5, 5.5):
            ax.axvline(boundary, color=plots.MUTED, lw=0.6, ls=":")
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set_xticks(range(len(tissues)), tissues, rotation=45, ha="right", fontsize=7)
        ax.set(ylabel="Set minus matched CpGs (SD)", title=title)

    by_tissue(axes[0], "aging", "entropy_slope", "Entropy gain with age (aged minus young)")
    axes[0].legend(fontsize=7, loc="upper left")
    plots.label_panel(axes[0], "a")
    by_tissue(axes[1], "aged dispersion", "dispersion", "Spread among aged mice")
    plots.label_panel(axes[1], "b")

    def contrast_panel(ax, frame: pl.DataFrame, title: str, xlabel: str) -> None:
        y, ticks, labels = 0.0, [], []
        for row in frame.iter_rows(named=True):
            color = plots.ACCENT if row["ci_high"] < 0 or row["ci_low"] > 0 else plots.MUTED
            ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=color, lw=1.2)
            ax.scatter([row["difference"]], [y], color=color, s=16)
            ticks.append(y)
            labels.append(f"{SET_LABELS[row['group']]} · {row['label']}")
            y += 1
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(ticks, labels, fontsize=6.5)
        ax.invert_yaxis()
        ax.set(xlabel=xlabel, title=title)

    gradient = summary.filter(pl.col("contrast") == "renewing - post-mitotic").with_columns(
        label=pl.format("{} {}", pl.col("measure"), pl.col("metric").str.replace_all("_", " "))
    )
    contrast_panel(
        axes[2],
        gradient.sort("group", "measure", "metric"),
        "Renewing minus post-mitotic tissues",
        "Difference (SD); < 0 = more protection where cells renew",
    )
    plots.label_panel(axes[2], "c")
    arms = summary.filter(pl.col("contrast").str.contains("- Aged_ctrl")).with_columns(
        label=pl.col("contrast").str.replace("Aged_", "").str.replace(" - Aged_ctrl", " vs control")
    )
    contrast_panel(
        axes[3],
        arms.sort("group", "contrast"),
        "Treatment effect on the protection (aged mice)",
        "Change of the set's dispersion difference (SD)",
    )
    plots.label_panel(axes[3], "d")
    fig.suptitle(
        "Does the protection follow cell turnover? 11 tissues of the same mice",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("mouse_tissues"))
