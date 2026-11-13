"""Extension, prediction 2: is the essential-gene protection larger in renewing than in post-mitotic tissue?

If selection among cells protects essential genes from entropy, the protection needs cell
turnover: it should be clear in blood and weak or absent in neurons. If active chromatin alone
protects them, it should be similar in both. Two tests on the Illumina 450k array, restricted to
promoter CpGs (TSS1500, TSS200, 5'UTR, first exon) on autosomes:

* **Lifespan** - Hannum whole blood (GSE40279, 19-101 y) against Jaffe prefrontal cortex
  (GSE74193, adult controls 20-97 y, one array per donor). Drift metrics are partial age slopes
  adjusted for processing plate, sex and ancestry, in blood also the recruitment site, and in
  cortex also the estimated neuron fraction. Different people, wide age ranges.
* **Same donors** - Lunnon (GSE59685): blood and four brain regions of the same donors, mostly
  70-100 y. With so narrow an age range, the measure is between-person dispersion at each CpG
  (adjusted for age, sex and disease status), which is the divergence accumulated by old age.

In every tissue a gene set is compared with background promoter CpGs of the same island status
and starting-level decile (``compare_groups``); the difference is expressed in standard
deviations of the metric so tissues are comparable. Blood minus brain gets a bootstrap interval:
independent resamples per cohort for the lifespan test, the same resampled donors in every tissue
for the paired test. Bootstrap workers are spawned processes (``_run_parallel``).
The prediction is a negative blood-minus-brain difference for core-essential
and housekeeping CpGs (more protection in blood).

Run on a machine with tens of GB of memory: ``entropy-age fetch cohort hannum jaffe lunnon`` then
``entropy-age figures tissues``. ``ENTROPY_AGE_WORKERS`` sets the bootstrap processes.
"""

from __future__ import annotations

import dataclasses
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.experiments.common import Context, dummies, rescale, write_table
from entropy_age.stochasticity import (
    baseline_strata,
    compare_groups,
    composition_components,
    dispersion_levels,
    drift_metrics,
)

SETS = ("core_essential", "housekeeping", "non_essential")
SET_LABELS = {"core_essential": "core essential", "housekeeping": "housekeeping", "non_essential": "non-essential"}
LIFESPAN_METRICS = ("entropy_slope", "dispersion_slope", "abs_slope", "kl_slope")
BACKGROUND = 60_000
MAX_MISSING = 0.05
LUNNON_BRAIN = ("frontal cortex", "superior temporal gyrus", "entorhinal cortex", "cerebellum")
# composition check (`figures tissues_cells`): name -> (components, orthogonal to age), None = as in `tissues`
ADJUSTMENTS: dict[str, tuple[int, bool] | None] = {
    "none": None,
    "pc5": (5, False),
    "pc10": (10, False),
    "pc5_age_free": (5, True),
}


@dataclass(frozen=True, slots=True)
class Prepared:
    """One tissue of one cohort, ready for the metric: values, ages, covariates and features."""

    label: str
    x: np.ndarray
    ages: np.ndarray
    covariates: np.ndarray | None
    features: list[str]
    donors: np.ndarray
    technical: np.ndarray | None  # the covariates other than age (plate, sex, ...), for the composition check


# ----------------------------------------------------------------------------- inputs
def _features(paths, rng: np.random.Generator) -> tuple[list[str], pl.DataFrame, pl.DataFrame]:
    annotation = AnnotationLayout(paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]))
    membership = link_membership(pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.gene_sets))
    membership = membership.filter(pl.col("group").is_in(SETS) & pl.col("feature").is_in(probes["cpg"].implode()))
    members = set(membership["feature"].to_list())
    rest = sorted(set(probes["cpg"].to_list()) - members)
    background = sorted(rng.choice(rest, size=min(BACKGROUND, len(rest)), replace=False).tolist())
    return sorted(members | set(background)), probes, membership


def _matrix(name: str, paths, samples: list[str], features: list[str]) -> FeatureMatrix:
    lay = cohort_layout(name, paths)
    frame = pl.scan_parquet(lay.betas).select("cpg", *samples).filter(pl.col("cpg").is_in(features)).collect()
    return FeatureMatrix.from_frame(frame, feature_col="cpg", samples=samples, max_missing=MAX_MISSING)


def _hannum(paths, features) -> Prepared:
    s = pl.read_parquet(cohort_layout("hannum", paths).samples).drop_nulls("age").sort("age", "sample")
    m = _matrix("hannum", paths, s["sample"].to_list(), features)
    cov = np.column_stack([dummies(s[c].to_numpy()) for c in ("plate", "source", "sex", "ethnicity")])
    return Prepared("blood (Hannum)", m.values, s["age"].to_numpy(), cov, m.features, s["donor"].to_numpy(), cov)


def _jaffe(paths, features) -> Prepared:
    s = (
        pl.read_parquet(cohort_layout("jaffe", paths).samples)
        .filter((pl.col("group") == "Control") & ~pl.col("drop_sample") & (pl.col("age") >= 20))
        .sort("best_qc", descending=True)
        .unique(subset="donor", keep="first", maintain_order=True)
        .drop_nulls(["age", "neuron_fraction"])
        .sort("age", "sample")
    )
    m = _matrix("jaffe", paths, s["sample"].to_list(), features)
    cov = np.column_stack(
        [*(dummies(s[c].to_numpy()) for c in ("plate", "sex", "race")), s["neuron_fraction"].to_numpy()]
    )
    return Prepared("cortex (Jaffe)", m.values, s["age"].to_numpy(), cov, m.features, s["donor"].to_numpy(), cov)


def _lunnon(paths, features) -> dict[str, Prepared]:
    s = pl.read_parquet(cohort_layout("lunnon", paths).samples).drop_nulls("age")
    blood_donors = s.filter(pl.col("tissue") == "whole blood")["donor"].unique()
    s = s.filter(pl.col("donor").is_in(blood_donors.implode()))
    out = {}
    for tissue in ("whole blood", *LUNNON_BRAIN):
        t = s.filter(pl.col("tissue") == tissue).sort("donor")
        m = _matrix("lunnon", paths, t["sample"].to_list(), features)
        technical = np.column_stack([dummies(t[c].to_numpy()) for c in ("sex", "ad_status")])
        cov = np.column_stack([t["age"].to_numpy(), technical])
        out[tissue] = Prepared(tissue, m.values, t["age"].to_numpy(), cov, m.features, t["donor"].to_numpy(), technical)
    return out


def _adjusted(p: Prepared, members: set[str], k: int, age_free: bool) -> tuple[Prepared, dict]:
    """``p`` with ``k`` composition components of its background CpGs added to the covariates.

    The components come from the CpGs outside every gene set, after the technical covariates (and,
    with ``age_free``, age) are regressed out, so they describe the shared movement of thousands of
    unrelated CpGs: in blood mostly the cell mix, in brain mostly the neuron-glia balance.
    """
    background = np.array([f not in members for f in p.features])
    base = p.technical
    if age_free:
        base = p.ages[:, None] if base is None else np.column_stack([p.ages, base])
    scores, explained = composition_components(p.x[:, background], k, covariates=base)
    covariates = scores if p.covariates is None else np.column_stack([p.covariates, scores])
    r_age = [abs(float(np.corrcoef(scores[:, j], p.ages)[0, 1])) for j in range(k)]
    info = {
        "tissue": p.label,
        "components": k,
        "age_free": age_free,
        "background_cpgs": int(background.sum()),
        "variance_explained": float(explained.sum()),
        "first_component": float(explained[0]),
        "max_abs_r_age": max(r_age),
    }
    return dataclasses.replace(p, covariates=covariates), info


# ----------------------------------------------------------------------------- metrics
def _with_strata(metrics: pl.DataFrame, probes: pl.DataFrame) -> pl.DataFrame:
    metrics = metrics.join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")
    decile = baseline_strata(metrics["baseline"].to_numpy(), 10)
    return metrics.with_columns(stratum=pl.format("{}|{}", pl.col("cpg_island").fill_null(False), pl.Series(decile)))


def _lifespan_metrics(p: Prepared, rows: np.ndarray | None = None) -> pl.DataFrame:
    rows = np.arange(p.x.shape[0]) if rows is None else rows
    cov = None if p.covariates is None else p.covariates[rows]
    if cov is not None:  # a resample can drop a level entirely; keep the design full rank
        cov = cov[:, cov.std(axis=0) > 0]
    return drift_metrics(p.x[rows], p.ages[rows], p.features, covariates=cov, spearman=False)


def _dispersion_metrics(p: Prepared, rows: np.ndarray | None = None) -> pl.DataFrame:
    rows = np.arange(p.x.shape[0]) if rows is None else rows
    cov = p.covariates[rows]
    cov = cov[:, cov.std(axis=0) > 0]
    return dispersion_levels(p.x[rows], p.features, covariates=cov)


def _matched_sd(metrics: pl.DataFrame, membership: pl.DataFrame, metric_names, permutations: int, seed: int):
    rows = []
    for metric in metric_names:
        table = compare_groups(
            metrics, membership, metric=metric, stratify_by="stratum", n_permutations=permutations, seed=seed
        )
        spread = float(metrics[metric].std())
        rows += [
            {"metric": metric, **r, "sd": spread, "matched_diff_sd": r["matched_diff"] / spread}
            for r in table.to_dicts()
        ]
    return pl.DataFrame(rows)


# ----------------------------------------------------------------------------- bootstrap workers
_STATE: dict = {}


def _boot_lifespan(job: tuple[str, int]) -> list[dict]:
    key, seed = job
    p, probes, membership = _STATE["lifespan"][key], _STATE["probes"], _STATE["membership"]
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.integers(0, p.x.shape[0], p.x.shape[0]))
    table = _matched_sd(_with_strata(_lifespan_metrics(p, rows), probes), membership, LIFESPAN_METRICS, 0, seed)
    return [{"cohort": key, "boot": seed, **r} for r in table.select("group", "metric", "matched_diff").to_dicts()]


def _boot_lunnon(seed: int) -> list[dict]:
    tissues, probes, membership = _STATE["lunnon"], _STATE["probes"], _STATE["membership"]
    rng = np.random.default_rng(seed)
    donors = np.unique(tissues["whole blood"].donors)
    picked = rng.choice(donors, size=donors.size, replace=True)
    out = []
    for tissue, p in tissues.items():
        index = {d: i for i, d in enumerate(p.donors)}
        rows = np.array([index[d] for d in picked if d in index])
        table = _matched_sd(_with_strata(_dispersion_metrics(p, rows), probes), membership, ("dispersion",), 0, seed)
        out += [
            {"tissue": tissue, "boot": seed, **r} for r in table.select("group", "metric", "matched_diff").to_dicts()
        ]
    return out


def _init_worker(state: dict) -> None:
    _STATE.clear()
    _STATE.update(state)


def _run_parallel(fn, jobs, workers: int, state: dict) -> list[dict]:
    """Run ``fn`` over ``jobs`` with ``state`` installed as the module's ``_STATE``.

    Workers are spawned, never forked: forking a process whose polars thread pool is running
    deadlocks the children on locks held by threads that do not exist in them. That happened on
    2026-11-09 (16 workers asleep after 19 CPU-seconds each). Each worker receives ``state``
    once, through the initializer.
    """
    if workers <= 1:
        _init_worker(state)
        return [row for job in jobs for row in fn(job)]
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=context, initializer=_init_worker, initargs=(state,)
    ) as pool:
        return [row for rows in pool.map(fn, jobs, chunksize=4) for row in rows]


def _interval(values: np.ndarray) -> tuple[float, float, float]:
    lo, hi = np.quantile(values, [0.025, 0.975])
    p = 2 * min(np.mean(values <= 0), np.mean(values >= 0))
    return float(lo), float(hi), float(min(1.0, max(p, 1 / values.size)))


# ----------------------------------------------------------------------------- run
def _inputs_missing(ctx: Context, figure: str) -> bool:
    needed = [cohort_layout(n, ctx.paths).betas for n in ("hannum", "jaffe", "lunnon")]
    needed.append(AnnotationLayout(ctx.paths).links_450k)
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        print(
            f"{figure}: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation` and "
            "`entropy-age fetch cohort hannum jaffe lunnon` first. Nothing was run.",
            file=sys.stderr,
        )
    return bool(missing)


def _settings(ctx: Context) -> tuple[int, int]:
    """Bootstrap workers (``ENTROPY_AGE_WORKERS``) and resamples (``ENTROPY_AGE_BOOT``, for smoke tests)."""
    workers = int(os.environ.get("ENTROPY_AGE_WORKERS", "1"))
    default = 200 if ctx.scale.name == "paper" else 100
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", default))
    if n_boot != default:
        print(f"tissues: ENTROPY_AGE_BOOT={n_boot} replaces the default {default} resamples", file=sys.stderr)
    return workers, n_boot


def _load(ctx: Context):
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    lifespan = {"blood": _hannum(ctx.paths, features), "cortex": _jaffe(ctx.paths, features)}
    return lifespan, _lunnon(ctx.paths, features), probes, membership


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "tissues"):
        return []
    plots.apply_style()
    out = ctx.out_dir("tissue_contrast")
    workers, n_boot = _settings(ctx)
    lifespan, lunnon, probes, membership = _load(ctx)
    summary = [
        {
            "part": "lifespan",
            "tissue": p.label,
            "samples": p.x.shape[0],
            "cpgs": len(p.features),
            "age_min": float(p.ages.min()),
            "age_max": float(p.ages.max()),
        }
        for p in lifespan.values()
    ] + [
        {
            "part": "same donors",
            "tissue": p.label,
            "samples": p.x.shape[0],
            "cpgs": len(p.features),
            "age_min": float(p.ages.min()),
            "age_max": float(p.ages.max()),
        }
        for p in lunnon.values()
    ]
    write_table(pl.DataFrame(summary), out / "cohorts")
    write_table(
        membership.group_by("group").agg(pl.len().alias("cpgs"), pl.col("gene").n_unique().alias("genes")),
        out / "set_sizes",
    )
    point, boot_life, boot_same, contrasts = _analyse(lifespan, lunnon, probes, membership, n_boot, workers)
    write_table(point, out / "point_estimates")
    boot_life.write_parquet(out / "bootstrap_lifespan.parquet")
    boot_same.write_parquet(out / "bootstrap_same_donors.parquet")
    write_table(contrasts, out / "tissue_differences")
    return [_figure(ctx, point, boot_life, boot_same, contrasts)]


def _analyse(lifespan: dict, lunnon: dict, probes, membership, n_boot: int, workers: int):
    """Point estimates, bootstraps and blood-minus-brain differences for one set of covariates."""
    point = []
    for key, p in lifespan.items():
        t = _matched_sd(_with_strata(_lifespan_metrics(p), probes), membership, LIFESPAN_METRICS, 1000, 1)
        point.append(t.with_columns(part=pl.lit("lifespan"), tissue=pl.lit(key)))
    for tissue, p in lunnon.items():
        t = _matched_sd(_with_strata(_dispersion_metrics(p), probes), membership, ("dispersion",), 1000, 2)
        point.append(t.with_columns(part=pl.lit("same donors"), tissue=pl.lit(tissue)))
    point = pl.concat(point, how="diagonal")
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))

    # bootstraps
    shared = {"probes": probes, "membership": membership}
    boot_life = pl.DataFrame(
        _run_parallel(
            _boot_lifespan,
            [(k, 1000 + b) for k in lifespan for b in range(n_boot)],
            workers,
            {**shared, "lifespan": lifespan},
        )
    )
    boot_same = pl.DataFrame(
        _run_parallel(_boot_lunnon, [2000 + b for b in range(n_boot)], workers, {**shared, "lunnon": lunnon})
    )
    sds = point.select("part", pl.col("tissue"), "metric", "sd")
    boot_life = rescale(
        boot_life, sds.filter(pl.col("part") == "lifespan").rename({"tissue": "cohort"}), ["cohort", "metric"]
    )
    boot_same = rescale(boot_same, sds.filter(pl.col("part") == "same donors"), ["tissue", "metric"])

    rows = []
    for group in SETS:
        for metric in LIFESPAN_METRICS:
            sel = boot_life.filter((pl.col("group") == group) & (pl.col("metric") == metric))
            blood = sel.filter(pl.col("cohort") == "blood").sort("boot")["matched_diff_sd"].to_numpy()
            cortex = sel.filter(pl.col("cohort") == "cortex").sort("boot")["matched_diff_sd"].to_numpy()
            k = min(blood.size, cortex.size)
            lo, hi, p = _interval(blood[:k] - cortex[:k])
            rows.append(
                {
                    "part": "lifespan",
                    "comparison": "blood - cortex",
                    "group": group,
                    "metric": metric,
                    "difference": float(np.median(blood[:k] - cortex[:k])),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                }
            )
        sel = boot_same.filter(pl.col("group") == group)
        blood = sel.filter(pl.col("tissue") == "whole blood").sort("boot")["matched_diff_sd"].to_numpy()
        for region in LUNNON_BRAIN:
            other = sel.filter(pl.col("tissue") == region).sort("boot")["matched_diff_sd"].to_numpy()
            k = min(blood.size, other.size)
            lo, hi, p = _interval(blood[:k] - other[:k])
            rows.append(
                {
                    "part": "same donors",
                    "comparison": f"blood - {region}",
                    "group": group,
                    "metric": "dispersion",
                    "difference": float(np.median(blood[:k] - other[:k])),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                }
            )
    return point, boot_life, boot_same, pl.DataFrame(rows)


def run_cells(ctx: Context) -> list[Path]:
    """Does the blood-minus-brain gap survive adjusting both tissues for composition?

    Neither blood cohort ships cell counts, and Lunnon's brain regions ship no neuron fraction, so
    every tissue gets the same reference-free correction: the leading principal components of its
    background CpGs as extra covariates (``ADJUSTMENTS``). Components that keep their age trend
    also absorb composition change with age (the conservative variant); age-free ones absorb only
    the person-to-person spread of the mix.
    """
    if _inputs_missing(ctx, "tissues_cells"):
        return []
    plots.apply_style()
    out = ctx.out_dir("tissue_contrast_cells")
    workers, n_boot = _settings(ctx)
    lifespan, lunnon, probes, membership = _load(ctx)
    members = set(membership["feature"].to_list())
    gaps, points, components = [], [], []
    for name, spec in ADJUSTMENTS.items():
        if spec is None:
            life, same = lifespan, lunnon
        else:
            k, age_free = spec
            life, same = {}, {}
            for target, source in ((life, lifespan), (same, lunnon)):
                for key, p in source.items():
                    target[key], info = _adjusted(p, members, k, age_free)
                    components.append({"adjustment": name, **info})
        print(f"tissues_cells: {name}", file=sys.stderr, flush=True)
        point, _, _, contrasts = _analyse(life, same, probes, membership, n_boot, workers)
        points.append(point.with_columns(adjustment=pl.lit(name)))
        gaps.append(contrasts.with_columns(adjustment=pl.lit(name)))
    gaps = pl.concat(gaps)
    write_table(pl.concat(points, how="diagonal"), out / "point_estimates")
    write_table(gaps, out / "tissue_differences")
    if components:
        write_table(pl.DataFrame(components), out / "components")
    return [_figure_cells(ctx, gaps)]


def _figure_cells(ctx: Context, gaps: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=6.4)
    names = list(ADJUSTMENTS)
    offsets = np.linspace(-0.3, 0.3, len(names))
    for ax, part, title in (
        (axes[0], "lifespan", "Lifespan cohorts: blood minus cortex"),
        (axes[1], "same donors", "Same donors: blood minus brain region"),
    ):
        frame = gaps.filter(pl.col("part") == part)
        keys = frame.select("group", "metric", "comparison").unique(maintain_order=True).to_dicts()
        for i, key in enumerate(keys):
            for j, name in enumerate(names):
                row = frame.filter(
                    (pl.col("group") == key["group"])
                    & (pl.col("metric") == key["metric"])
                    & (pl.col("comparison") == key["comparison"])
                    & (pl.col("adjustment") == name)
                ).row(0, named=True)
                yy = i + offsets[j]
                color = plots.PALETTE[j % len(plots.PALETTE)]
                ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=color, lw=1.1)
                ax.scatter([row["difference"]], [yy], color=color, s=12, label=name if i == 0 else None)
        labels = [
            f"{SET_LABELS[k['group']]} · "
            + (k["metric"].replace("_", " ") if part == "lifespan" else k["comparison"].split(" - ")[1])
            for k in keys
        ]
        ax.set_yticks(range(len(keys)), labels, fontsize=6.5)
        ax.invert_yaxis()
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set(xlabel="Blood minus brain (SD); < 0 = more protection in blood", title=title)
        ax.title.set_fontsize(10)
    handles, names_shown = axes[0].get_legend_handles_labels()
    fig.legend(handles, names_shown, fontsize=8, title="adjustment", loc="lower center", ncols=len(names))
    plots.label_panel(axes[0], "a")
    plots.label_panel(axes[1], "b")
    fig.suptitle(
        "Does the blood-brain gap survive a composition adjustment?", x=0.02, ha="left", fontsize=12, fontweight="bold"
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    return plots.save(fig, ctx.figure_path("tissue_contrast_cells"))


def _figure(ctx: Context, point, boot_life, boot_same, contrasts) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    colors = {"blood": plots.PALETTE[1], "cortex": plots.PALETTE[0]}

    def ci(frame: pl.DataFrame) -> tuple[float, float]:
        v = frame["matched_diff_sd"].to_numpy()
        return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))

    ax, y, ticks, labels = axes[0], 0.0, [], []
    for group in SETS:
        for metric in LIFESPAN_METRICS:
            for k, cohort in enumerate(("blood", "cortex")):
                est = point.filter(
                    (pl.col("part") == "lifespan")
                    & (pl.col("tissue") == cohort)
                    & (pl.col("group") == group)
                    & (pl.col("metric") == metric)
                )
                lo, hi = ci(
                    boot_life.filter(
                        (pl.col("cohort") == cohort) & (pl.col("group") == group) & (pl.col("metric") == metric)
                    )
                )
                yy = y + (k - 0.5) * 0.3
                ax.plot([lo, hi], [yy, yy], color=colors[cohort], lw=1.2)
                ax.scatter(
                    est["matched_diff_sd"],
                    [yy],
                    color=colors[cohort],
                    s=18,
                    label=cohort if (group, metric) == (SETS[0], LIFESPAN_METRICS[0]) else None,
                )
            ticks.append(y)
            labels.append(f"{SET_LABELS[group]} · {metric.replace('_', ' ')}")
            y += 1
        y += 0.5
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(ticks, labels, fontsize=7)
    ax.invert_yaxis()
    ax.legend(fontsize=7, loc="lower right")
    ax.set(xlabel="Set minus matched CpGs (SD), 95 % bootstrap", title="Lifespan cohorts: blood vs cortex")
    plots.label_panel(ax, "a")

    def diff_panel(ax, frame: pl.DataFrame, title: str) -> None:
        y, ticks, labels = 0.0, [], []
        for row in frame.iter_rows(named=True):
            color = plots.ACCENT if row["ci_high"] < 0 or row["ci_low"] > 0 else plots.MUTED
            ax.plot([row["ci_low"], row["ci_high"]], [y, y], color=color, lw=1.2)
            ax.scatter([row["difference"]], [y], color=color, s=18)
            ticks.append(y)
            labels.append(
                f"{SET_LABELS[row['group']]} · {row['metric'].replace('_', ' ')}"
                + ("" if row["part"] == "lifespan" else f" · {row['comparison'].split(' - ')[1]}")
            )
            y += 1
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(ticks, labels, fontsize=6.5)
        ax.invert_yaxis()
        ax.set(xlabel="Blood minus brain (SD); < 0 = more protection in blood", title=title)

    diff_panel(axes[1], contrasts.filter(pl.col("part") == "lifespan"), "Prediction 2, lifespan cohorts")
    plots.label_panel(axes[1], "b")

    ax, y, ticks, labels = axes[2], 0.0, [], []
    tissues = ("whole blood", *LUNNON_BRAIN)
    for group in SETS:
        for k, tissue in enumerate(tissues):
            est = point.filter(
                (pl.col("part") == "same donors") & (pl.col("tissue") == tissue) & (pl.col("group") == group)
            )
            lo, hi = ci(boot_same.filter((pl.col("tissue") == tissue) & (pl.col("group") == group)))
            yy = y + (k - 2) * 0.17
            color = plots.PALETTE[1] if tissue == "whole blood" else plots.PALETTE[(k + 1) % len(plots.PALETTE)]
            ax.plot([lo, hi], [yy, yy], color=color, lw=1.2)
            ax.scatter(est["matched_diff_sd"], [yy], color=color, s=14, label=tissue if group == SETS[0] else None)
        ticks.append(y)
        labels.append(SET_LABELS[group])
        y += 1.4
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(ticks, labels)
    ax.invert_yaxis()
    ax.legend(fontsize=6.5, loc="lower right")
    ax.set(xlabel="Set minus matched CpGs, between-person dispersion (SD)", title="Same donors (Lunnon)")
    plots.label_panel(ax, "c")

    diff_panel(axes[3], contrasts.filter(pl.col("part") == "same donors"), "Prediction 2, same donors")
    plots.label_panel(axes[3], "d")
    fig.suptitle(
        "Is essential-gene protection larger in blood than in brain?", x=0.02, ha="left", fontsize=12, fontweight="bold"
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("tissue_contrast"))
