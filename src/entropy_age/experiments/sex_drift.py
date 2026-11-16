"""Sex differences in drift and in the shield (roadmap RM35 item 2; BEYOND_SHIELD idea 10).

Women live longer than men, and several clocks run slower in women. Does the stochastic part differ too? Per
tissue, the donors are split by sex and each half gets what the pooled analysis gets: the per-CpG age slopes of
spread and of relative entropy from the young state (``drift_metrics``, the technical covariates other than sex),
the shield's matched gaps for the essential sets (``tissue_contrast._matched_sd``, island x young-level strata,
promoter CpGs), and the drift rate itself, the median slope of the 60,000 background promoter CpGs. Women minus
men, with a donor bootstrap within each sex.

Tissues: Hannum blood, Jaffe cortex, Lunnon's blood and four brain regions (same donors), and the GTEx tissues
with at least 20 donors of each sex (prostate, testis and ovary drop out by design).

``entropy-age figures sex_drift`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.human_tissues import _gtex_tissues
from entropy_age.experiments.tissue_contrast import (
    Prepared,
    _features,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _lunnon,
    _matched_sd,
    _settings,
    _with_strata,
)

METRICS = ("dispersion_slope", "kl_slope")
MIN_PER_SEX = 20
N_BOOT_MAX = 100


def sex_codes(values) -> np.ndarray:
    """``F``/``M`` from the cohorts' spellings (F, FEMALE, female, ...); anything else is empty."""
    out = [str(v).strip().upper()[:1] if v is not None else "" for v in values]
    return np.array([v if v in ("F", "M") else "" for v in out])


def measures(p: Prepared, rows: np.ndarray, probes, membership, members: set[str]) -> list[dict]:
    """Matched gaps of the sets and the background's median slopes, for the samples ``rows``."""
    cov = None if p.covariates is None else p.covariates[rows]
    q = dataclasses.replace(p, x=p.x[rows], ages=p.ages[rows], covariates=cov, donors=p.donors[rows])
    m = _with_strata(_lifespan_metrics(q), probes)
    gaps = _matched_sd(m, membership, METRICS, 0, 0)
    out = [{"metric": r["metric"], "group": r["group"], "value": r["matched_diff_sd"]} for r in gaps.to_dicts()]
    background = m.filter(~pl.col("feature").is_in(list(members)))
    for metric in (*METRICS, "abs_slope"):
        out.append({"metric": metric, "group": "background median", "value": float(background[metric].median())})
    return out


def _sexes(cohort: str, ctx: Context, p: Prepared) -> np.ndarray:
    s = pl.read_parquet(cohort_layout(cohort, ctx.paths).samples)
    first = s.group_by("donor").agg(pl.col("sex").first())
    lookup = dict(zip(first["donor"].to_list(), first["sex"].to_list(), strict=True))
    return sex_codes([lookup.get(d) for d in p.donors])


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("sex_drift")
    _, n_boot = _settings(ctx)
    n_boot = min(n_boot, N_BOOT_MAX)
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    members = set(membership["feature"].to_list())
    tissues: list[tuple[str, Prepared, np.ndarray]] = []
    h = _hannum(ctx.paths, features)
    tissues.append((h.label, h, _sexes("hannum", ctx, h)))
    j = _jaffe(ctx.paths, features)
    tissues.append((j.label, j, _sexes("jaffe", ctx, j)))
    for tissue, p in _lunnon(ctx.paths, features).items():
        p = dataclasses.replace(p, covariates=p.technical)  # age is the regressor, not a covariate
        tissues.append((f"{tissue} (Lunnon)", p, _sexes("lunnon", ctx, p)))
    for (tissue, adjustment), p in _gtex_tissues(ctx.paths, features, members).items():
        if adjustment == "covariates":
            tissues.append((f"{tissue} (GTEx)", p, _sexes("gtex", ctx, p)))

    rng = ctx.rng(35_002)
    point, boot, counts = [], [], []
    for label, p, sex in tissues:
        rows = {s: np.flatnonzero(sex == s) for s in ("F", "M")}
        counts.append({"tissue": label, "women": rows["F"].size, "men": rows["M"].size})
        if min(r.size for r in rows.values()) < MIN_PER_SEX:
            print(f"sex_drift: {label}: {rows['F'].size} women, {rows['M'].size} men; skipped", file=sys.stderr)
            continue
        p = dataclasses.replace(p, x=p.x.astype(np.float32))
        for s, r in rows.items():
            point += [{**d, "tissue": label, "sex": s} for d in measures(p, r, probes, membership, members)]
        for b in range(n_boot):
            for s, r in rows.items():
                draw = np.sort(rng.choice(r, size=r.size, replace=True))
                ms = measures(p, draw, probes, membership, members)
                boot += [{**d, "tissue": label, "sex": s, "b": b} for d in ms]
        print(f"sex_drift: {label}: {rows['F'].size} women, {rows['M'].size} men", file=sys.stderr, flush=True)
    point = pl.DataFrame(point)
    write_table(pl.DataFrame(counts), out / "donors")
    keys = ["tissue", "metric", "group"]
    diff = _women_minus_men(point, keys)
    if boot:
        boot = pl.DataFrame(boot)
        boot.write_parquet(out / "bootstrap.parquet")
        bd = _women_minus_men(boot, [*keys, "b"])
        ci = bd.group_by(keys).agg(
            lo=pl.col("difference").quantile(0.025),
            hi=pl.col("difference").quantile(0.975),
            p_boot=2 * pl.min_horizontal((pl.col("difference") <= 0).mean(), (pl.col("difference") >= 0).mean()),
        )
        diff = diff.join(ci, on=keys, how="left")
    write_table(point, out / "by_sex")
    write_table(diff.sort(keys), out / "women_minus_men")
    return [_figure(ctx, diff)]


def _women_minus_men(frame: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    w = frame.filter(pl.col("sex") == "F").select(*keys, women=pl.col("value"))
    m = frame.filter(pl.col("sex") == "M").select(*keys, men=pl.col("value"))
    return w.join(m, on=keys).with_columns(difference=pl.col("women") - pl.col("men"))


def _figure(ctx: Context, diff: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    panels = (
        ("background median", "drift rate: women - men (slope of spread)"),
        ("housekeeping", "shield gap: women - men (SD)"),
    )
    for ax, (group, title), letter in zip(axes, panels, "ab", strict=False):
        d = diff.filter((pl.col("metric") == "dispersion_slope") & (pl.col("group") == group))
        tissues = d["tissue"].to_list()
        y = np.arange(len(tissues))
        ax.plot(d["difference"], y, "o", color=plots.PALETTE[0])
        if "lo" in d.columns:
            for k, (lo, hi) in enumerate(zip(d["lo"], d["hi"], strict=True)):
                ax.plot([lo, hi], [k, k], color=plots.PALETTE[0], lw=1.2)
        ax.set_yticks(y, tissues if letter == "a" else [""] * len(tissues), fontsize=7)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.set(title=title)
        plots.label_panel(ax, letter)
    return plots.save(fig, ctx.figure_path("sex_drift"))
