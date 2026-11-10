"""Prediction 2 in human tissues: GTEx (GSE213478) along a turnover gradient, and muscle fibres (GSE304045).

**GTEx.** Nine tissues of 424 post-mortem donors on the EPIC array, ages in decade bins. The mouse
test (``mouse_tissues``) found the essential-gene protection in every tissue whose cells divide,
even slowly, and missing in cortex and hippocampus: a step rather than a slope. GTEx has one
post-mitotic tissue, skeletal muscle, so the human test is whether muscle is the tissue without the
protection. Classes: renewing (whole blood, colon, testis), slow (lung, kidney, breast, prostate,
ovary), post-mitotic (skeletal muscle). Drift metrics are partial slopes on the bin midpoint,
adjusted for sex, autolysis score and smoking; a second version adds five composition components
of background CpGs, made orthogonal to age (the ``tissues_cells`` check), which also absorb plate.

Muscle has 47 donors and no one under 30. Fewer donors make the per-CpG slopes noisier, which
inflates the metric's SD and shrinks every SD-scaled difference toward zero: on its own that would
look like "less protection in muscle". So the intervals come from subsamples of the same size in
every tissue (``SUBSAMPLE`` donors, drawn without replacement), and the class contrasts are taken
over those subsamples. The full-data estimate of each tissue is written too, but tissues of
different size are not compared on it.

**Muscle fibres.** Pooled type I fibres, pooled type II fibres and whole muscle of the same ~12
adults, before and after 12 weeks of interval training. Fibre nuclei never divide; whole muscle
also holds satellite cells, fibroblasts and endothelium, which do. Without ages, the measure is
between-person dispersion (``dispersion_levels``, adjusted for sex). All three preparations share
the same people, so sample size cannot bias the comparison. Selection predicts more protection in
whole muscle than in pure fibres: a negative ``WM - T1`` and ``WM - T2`` difference for the
essential sets. Between-person spread in a dozen adults of unknown age is a weak signal, so this is
exploratory.

Both use 450k promoter CpGs on autosomes (``tissue_contrast._features``), present on both EPIC
versions, with island x baseline-decile strata.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import (
    _STATE,
    LIFESPAN_METRICS,
    SET_LABELS,
    SETS,
    Prepared,
    _adjusted,
    _features,
    _interval,
    _lifespan_metrics,
    _matched_sd,
    _matrix,
    _run_parallel,
    _settings,
    _with_strata,
)
from entropy_age.stochasticity import baseline_strata, dispersion_levels

#: tissue -> turnover class, in plotting order
TURNOVER: dict[str, str] = {
    "whole blood": "renewing",
    "colon": "renewing",
    "testis": "renewing",
    "lung": "slow",
    "kidney": "slow",
    "breast": "slow",
    "prostate": "slow",
    "ovary": "slow",
    "skeletal muscle": "post-mitotic",
}
CLASSES = ("renewing", "slow", "post-mitotic")
#: name -> (components, orthogonal to age); None = technical covariates only
ADJUSTMENTS: dict[str, tuple[int, bool] | None] = {"covariates": None, "pc5_age_free": (5, True)}
SUBSAMPLE = 32  # donors per tissue and resample; the smallest tissue has 47
TECHNICAL = ("sex", "autolysis", "smoker")
FIBRES = ("T1", "T2", "WM")
FIBRE_LABELS = {"T1": "type I fibres", "T2": "type II fibres", "WM": "whole muscle"}
TIMEPOINTS = ("pre", "12 weeks")


# ============================================================================= GTEx
def _gtex_tissues(paths, features: list[str], members: set[str]) -> dict[tuple[str, str], Prepared]:
    """(tissue, adjustment) -> the tissue's donors, ready for ``_lifespan_metrics``."""
    samples = pl.read_parquet(cohort_layout("gtex", paths).samples)
    out = {}
    for tissue in TURNOVER:
        t = samples.filter(pl.col("tissue") == tissue).sort("donor")
        if t.height == 0:
            continue
        m = _matrix("gtex", paths, t["sample"].to_list(), features)
        technical = np.column_stack([dummies(t[c].fill_null("NA").to_numpy()) for c in TECHNICAL])
        base = Prepared(tissue, m.values, t["age"].to_numpy(), technical, m.features, t["donor"].to_numpy(), technical)
        for name, adjustment in ADJUSTMENTS.items():
            out[(tissue, name)] = base if adjustment is None else _adjusted(base, members, *adjustment)[0]
    return out


def _tissue_rows(tissue: str, adjustment: str, p: Prepared, rows, probes, membership, permutations, seed):
    metrics = _with_strata(_lifespan_metrics(p, rows), probes)
    table = _matched_sd(metrics, membership, LIFESPAN_METRICS, permutations, seed)
    return table.with_columns(tissue=pl.lit(tissue), turnover=pl.lit(TURNOVER[tissue]), adjustment=pl.lit(adjustment))


def _subsample(rng: np.random.Generator, n: int, size: int) -> np.ndarray:
    """``size`` distinct rows of ``n`` in ascending order (every row when ``n <= size``)."""
    if n <= size:
        return np.arange(n)
    return np.sort(rng.choice(n, size=size, replace=False))


def _boot_gtex(seed: int) -> list[dict]:
    """One resample: ``SUBSAMPLE`` donors of every tissue, without replacement, the same in both adjustments."""
    tissues, probes, membership, size = _STATE["tissues"], _STATE["probes"], _STATE["membership"], _STATE["size"]
    rng = np.random.default_rng(seed)
    out = []
    for tissue in TURNOVER:
        if (tissue, "covariates") not in tissues:
            continue
        n = tissues[(tissue, "covariates")].x.shape[0]
        rows = _subsample(rng, n, size)
        for adjustment in ADJUSTMENTS:
            table = _tissue_rows(tissue, adjustment, tissues[(tissue, adjustment)], rows, probes, membership, 0, seed)
            out += [
                {"boot": seed, **r}
                for r in table.select("tissue", "adjustment", "metric", "group", "matched_diff_sd").to_dicts()
            ]
    return out


def summarise_gtex(boot: pl.DataFrame, point: pl.DataFrame) -> pl.DataFrame:
    """Class contrasts and every tissue against skeletal muscle, over the size-matched resamples.

    ``difference`` is the median over resamples (negative = more protection in the first tissue
    or class); ``full_data`` is the same difference of the full-data estimates, for reference only.
    """
    boot = boot.with_columns(turnover=pl.col("tissue").replace_strict(TURNOVER))
    point = point.with_columns(turnover=pl.col("tissue").replace_strict(TURNOVER))
    keys = ("adjustment", "metric", "group")
    rows = []

    def add(contrast: str, diff: np.ndarray, full: float, key: tuple) -> None:
        if diff.size < 10:
            return
        lo, hi, p = _interval(diff)
        rows.append(
            dict(zip(keys, key, strict=True))
            | {
                "contrast": contrast,
                "difference": float(np.median(diff)),
                "ci_low": lo,
                "ci_high": hi,
                "p_bootstrap": p,
                "full_data": full,
            }
        )

    for key, sel in boot.group_by(*keys):
        full = point.filter(*(pl.col(k) == v for k, v in zip(keys, key, strict=True)))
        classes = sel.group_by("boot", "turnover").agg(pl.col("matched_diff_sd").mean())
        wide = classes.pivot(on="turnover", index="boot", values="matched_diff_sd").sort("boot")
        full_class = dict(full.group_by("turnover").agg(pl.col("matched_diff_sd").mean()).iter_rows())
        for faster in ("renewing", "slow"):
            if faster in wide.columns and "post-mitotic" in wide.columns:
                diff = (wide[faster] - wide["post-mitotic"]).drop_nulls().to_numpy()
                add(f"{faster} - post-mitotic", diff, full_class[faster] - full_class["post-mitotic"], key)
        per_tissue = sel.pivot(on="tissue", index="boot", values="matched_diff_sd").sort("boot")
        full_tissue = dict(full.select("tissue", "matched_diff_sd").iter_rows())
        if "skeletal muscle" not in per_tissue.columns:
            continue
        for tissue in TURNOVER:
            if tissue == "skeletal muscle" or tissue not in per_tissue.columns:
                continue
            diff = (per_tissue[tissue] - per_tissue["skeletal muscle"]).drop_nulls().to_numpy()
            add(f"{tissue} - skeletal muscle", diff, full_tissue[tissue] - full_tissue["skeletal muscle"], key)
    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows).sort("adjustment", "contrast", "metric", "group")


def _missing(ctx: Context, figure: str, cohort: str) -> bool:
    needed = [cohort_layout(cohort, ctx.paths).betas, AnnotationLayout(ctx.paths).links_450k]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        print(
            f"{figure}: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation` and "
            f"`entropy-age fetch cohort {cohort}` first. Nothing was run.",
            file=sys.stderr,
        )
    return bool(missing)


def run(ctx: Context) -> list[Path]:
    if _missing(ctx, "gtex", "gtex"):
        return []
    plots.apply_style()
    out = ctx.out_dir("gtex_tissues")
    workers, n_boot = _settings(ctx)
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    members = set(membership["feature"].to_list())
    tissues = _gtex_tissues(ctx.paths, features, members)
    write_table(
        pl.DataFrame(
            [
                {
                    "tissue": tissue,
                    "turnover": TURNOVER[tissue],
                    "donors": p.x.shape[0],
                    "cpgs": len(p.features),
                    "age_min_bin_mid": float(p.ages.min()),
                    "age_max_bin_mid": float(p.ages.max()),
                    "age_mean": float(p.ages.mean()),
                }
                for (tissue, adjustment), p in tissues.items()
                if adjustment == "covariates"
            ]
        ),
        out / "tissues",
    )
    point = pl.concat(
        [
            _tissue_rows(tissue, adjustment, p, np.arange(p.x.shape[0]), probes, membership, 1000, 3)
            for (tissue, adjustment), p in tissues.items()
        ],
        how="diagonal",
    )
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))
    write_table(point, out / "point_estimates")

    state = {"tissues": tissues, "probes": probes, "membership": membership, "size": SUBSAMPLE}
    boot = pl.DataFrame(_run_parallel(_boot_gtex, [7000 + b for b in range(n_boot)], workers, state))
    boot.write_parquet(out / "bootstrap.parquet")
    summary = summarise_gtex(boot, point)
    write_table(summary, out / "contrasts")
    return [_figure_gtex(ctx, boot, summary)]


# ============================================================================= muscle fibres
@dataclass(frozen=True, slots=True)
class FibreSet:
    """One timepoint: the same people in every preparation, rows in ``donors`` order."""

    timepoint: str
    x: dict[str, np.ndarray]  # fibre -> (people, features)
    sex: np.ndarray
    donors: np.ndarray
    features: list[str]
    baseline: np.ndarray  # whole-muscle mean, the common strata for every preparation


def _fibre_sets(paths, features: list[str]) -> dict[str, FibreSet]:
    samples = pl.read_parquet(cohort_layout("fibres", paths).samples)
    first = samples.sort("replicate").unique(subset=["donor", "fibre", "timepoint"], keep="first")
    out = {}
    for timepoint in TIMEPOINTS:
        s = first.filter(pl.col("timepoint") == timepoint)
        complete = s.group_by("donor").agg(pl.col("fibre").n_unique().alias("k")).filter(pl.col("k") == len(FIBRES))
        donors = sorted(complete["donor"].to_list())
        if len(donors) < 4:
            continue
        s = s.filter(pl.col("donor").is_in(donors))
        columns = {f: s.filter(pl.col("fibre") == f).sort("donor")["sample"].to_list() for f in FIBRES}
        m = _matrix("fibres", paths, [c for f in FIBRES for c in columns[f]], features)
        n = len(donors)
        x = {f: m.values[i * n : (i + 1) * n] for i, f in enumerate(FIBRES)}
        sex = s.filter(pl.col("fibre") == "WM").sort("donor")["sex"].fill_null("NA").to_numpy()
        out[timepoint] = FibreSet(timepoint, x, sex, np.array(donors), m.features, x["WM"].mean(axis=0))
    return out


def fibre_measures(fs: FibreSet, rows: np.ndarray, probes: pl.DataFrame, membership, permutations: int, seed: int):
    """Matched dispersion difference of every set in every preparation, on the people in ``rows``."""
    strata = baseline_strata(fs.baseline, 10)
    island = (
        pl.DataFrame({"feature": fs.features})
        .join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")["cpg_island"]
        .fill_null(False)
    )
    cov = dummies(fs.sex[rows])
    cov = cov[:, cov.std(axis=0) > 0] if cov.shape[1] else None
    tables = []
    for fibre in FIBRES:
        metrics = dispersion_levels(
            fs.x[fibre][rows], fs.features, covariates=cov if cov is not None and cov.size else None
        )
        metrics = metrics.with_columns(stratum=pl.format("{}|{}", island, pl.Series(strata)))
        table = _matched_sd(metrics, membership, ("dispersion",), permutations, seed)
        tables.append(table.with_columns(fibre=pl.lit(fibre), timepoint=pl.lit(fs.timepoint)))
    return pl.concat(tables, how="diagonal")


def _boot_fibres(seed: int) -> list[dict]:
    sets, probes, membership = _STATE["sets"], _STATE["probes"], _STATE["membership"]
    rng = np.random.default_rng(seed)
    out = []
    for fs in sets.values():
        rows = np.sort(rng.integers(0, fs.donors.size, fs.donors.size))  # the same people in every preparation
        table = fibre_measures(fs, rows, probes, membership, 0, seed)
        out += [{"boot": seed, **r} for r in table.select("timepoint", "fibre", "group", "matched_diff_sd").to_dicts()]
    return out


def summarise_fibres(boot: pl.DataFrame) -> pl.DataFrame:
    """Whole muscle minus each fibre type (negative = more protection in whole muscle)."""
    rows = []
    for (timepoint, group), sel in boot.group_by("timepoint", "group"):
        wide = sel.pivot(on="fibre", index="boot", values="matched_diff_sd").sort("boot")
        for fibre in ("T1", "T2"):
            diff = (wide["WM"] - wide[fibre]).drop_nulls().to_numpy()
            if diff.size < 10:
                continue
            lo, hi, p = _interval(diff)
            rows.append(
                {
                    "timepoint": timepoint,
                    "group": group,
                    "contrast": f"WM - {fibre}",
                    "difference": float(np.median(diff)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                }
            )
    return pl.DataFrame(rows).sort("timepoint", "group", "contrast") if rows else pl.DataFrame()


def run_fibres(ctx: Context) -> list[Path]:
    if _missing(ctx, "fibres", "fibres"):
        return []
    plots.apply_style()
    out = ctx.out_dir("muscle_fibres")
    workers, n_boot = _settings(ctx)
    features, probes, membership = _features(ctx.paths, ctx.rng(150_000))
    sets = _fibre_sets(ctx.paths, features)
    write_table(
        pl.DataFrame(
            [{"timepoint": fs.timepoint, "people": fs.donors.size, "cpgs": len(fs.features)} for fs in sets.values()]
        ),
        out / "people",
    )
    point = pl.concat(
        [fibre_measures(fs, np.arange(fs.donors.size), probes, membership, 1000, 3) for fs in sets.values()],
        how="diagonal",
    )
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))
    write_table(point, out / "point_estimates")
    state = {"sets": sets, "probes": probes, "membership": membership}
    boot = pl.DataFrame(_run_parallel(_boot_fibres, [8000 + b for b in range(n_boot)], workers, state))
    boot.write_parquet(out / "bootstrap.parquet")
    summary = summarise_fibres(boot)
    write_table(summary, out / "contrasts")
    return [_figure_fibres(ctx, point, boot, summary)]


# ============================================================================= figures
def _contrast_panel(ax, frame: pl.DataFrame, title: str, xlabel: str) -> None:
    y, ticks, labels = 0.0, [], []
    for row in frame.iter_rows(named=True):
        color = plots.ACCENT if row["ci_high"] < 0 or row["ci_low"] > 0 else plots.MUTED
        ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=color, lw=1.2)
        ax.scatter([row["difference"]], [y], color=color, s=16)
        ticks.append(y)
        labels.append(row["label"])
        y += 1
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(ticks, labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set(xlabel=xlabel, title=title)


def _figure_gtex(ctx: Context, boot: pl.DataFrame, summary: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=5.6)
    colors = {g: plots.PALETTE[i] for i, g in enumerate(SETS)}
    tissues = [t for t in TURNOVER if t in set(boot["tissue"].to_list())]

    def by_tissue(ax, adjustment: str, title: str) -> None:
        for j, group in enumerate(SETS):
            for i, tissue in enumerate(tissues):
                v = boot.filter(
                    (pl.col("tissue") == tissue)
                    & (pl.col("adjustment") == adjustment)
                    & (pl.col("metric") == "entropy_slope")
                    & (pl.col("group") == group)
                )["matched_diff_sd"].to_numpy()
                if v.size == 0:
                    continue
                xx = i + (j - 1) * 0.22
                lo, hi = np.quantile(v, [0.025, 0.975])
                ax.plot([xx, xx], [lo, hi], color=colors[group], lw=1.1)
                ax.scatter([xx], [np.median(v)], color=colors[group], s=14, label=SET_LABELS[group] if i == 0 else None)
        for boundary in (2.5, 7.5):
            ax.axvline(boundary, color=plots.MUTED, lw=0.6, ls=":")
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set_xticks(range(len(tissues)), tissues, rotation=45, ha="right", fontsize=7)
        ax.set(ylabel="Set minus matched CpGs (SD)", title=title)

    by_tissue(axes[0], "covariates", f"Entropy gain with age, {SUBSAMPLE} donors per tissue")
    plots.label_panel(axes[0], "a")
    by_tissue(axes[1], "pc5_age_free", "Same, plus 5 age-free composition components")
    plots.label_panel(axes[1], "b")

    entropy = summary.filter(pl.col("metric") == "entropy_slope")
    classes = entropy.filter(pl.col("contrast").str.ends_with("post-mitotic")).with_columns(
        label=pl.format(
            "{} · {} · {}",
            pl.col("group").replace_strict(SET_LABELS),
            pl.col("contrast").str.replace(" - post-mitotic", ""),
            pl.col("adjustment").replace({"covariates": "cov", "pc5_age_free": "+PCs"}),
        )
    )
    _contrast_panel(
        axes[2],
        classes.sort("group", "contrast", "adjustment"),
        "Turnover class minus skeletal muscle",
        "Difference (SD); < 0 = more protection than in muscle",
    )
    plots.label_panel(axes[2], "c")
    each = entropy.filter(
        pl.col("contrast").str.ends_with("skeletal muscle")
        & (pl.col("adjustment") == "covariates")
        & pl.col("group").is_in(["core_essential", "housekeeping"])
    ).with_columns(
        label=pl.format(
            "{} · {}",
            pl.col("group").replace_strict(SET_LABELS),
            pl.col("contrast").str.replace(" - skeletal muscle", ""),
        )
    )
    _contrast_panel(
        axes[3], each.sort("group", "contrast"), "Each tissue minus skeletal muscle", "Difference (SD), entropy gain"
    )
    plots.label_panel(axes[3], "d")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(
        "Is muscle the human tissue without the protection? GTEx, 9 tissues",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return plots.save(fig, ctx.figure_path("gtex_tissues"))


def _figure_fibres(ctx: Context, point: pl.DataFrame, boot: pl.DataFrame, summary: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=5.6)
    colors = {g: plots.PALETTE[i] for i, g in enumerate(SETS)}
    ax = axes[0]
    slots = [(t, f) for t in TIMEPOINTS for f in FIBRES]
    for j, group in enumerate(SETS):
        for i, (timepoint, fibre) in enumerate(slots):
            est = point.filter(
                (pl.col("timepoint") == timepoint) & (pl.col("fibre") == fibre) & (pl.col("group") == group)
            )
            if est.height == 0:
                continue
            v = boot.filter(
                (pl.col("timepoint") == timepoint) & (pl.col("fibre") == fibre) & (pl.col("group") == group)
            )["matched_diff_sd"].to_numpy()
            xx = i + (j - 1) * 0.22
            if v.size:
                lo, hi = np.quantile(v, [0.025, 0.975])
                ax.plot([xx, xx], [lo, hi], color=colors[group], lw=1.1)
            ax.scatter(
                [xx], est["matched_diff_sd"], color=colors[group], s=14, label=SET_LABELS[group] if i == 0 else None
            )
    ax.axvline(len(FIBRES) - 0.5, color=plots.MUTED, lw=0.6, ls=":")
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set_xticks(range(len(slots)), [f"{FIBRE_LABELS[f]}\n{t}" for t, f in slots], rotation=45, ha="right", fontsize=7)
    ax.set(ylabel="Set minus matched CpGs (SD)", title="Spread between people, by preparation")
    ax.legend(fontsize=7, loc="best")
    plots.label_panel(ax, "a")
    if summary.height:
        frame = summary.with_columns(
            label=pl.format(
                "{} · {} · {}", pl.col("group").replace_strict(SET_LABELS), pl.col("contrast"), pl.col("timepoint")
            )
        )
        _contrast_panel(
            axes[1],
            frame.sort("group", "contrast", "timepoint"),
            "Whole muscle minus fibres",
            "Difference (SD); < 0 = more protection in whole muscle",
        )
    plots.label_panel(axes[1], "b")
    fig.suptitle(
        "Muscle fibres against whole muscle, same people (GSE304045)", x=0.02, ha="left", fontsize=12, fontweight="bold"
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("muscle_fibres"))
