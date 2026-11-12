"""One cell type at a time: the shield, drift and clocks in sorted blood cells (roadmap RM31)

Asked by antonkulaga (2026-11-12): clocks carry a cell-composition part that may hide what is stochastic; with one
cell type, a changing mix cannot make drift, the shield or a clock's age signal. MESA (Reynolds et al. 2014)
measured CD14+ monocytes of about 1,200 people aged 44-83 (GSE56046) and CD4+ T cells of about 200 (GSE56581) on
the 450k array, with age, a race-sex-site code and expression-based contamination scores (B, T, NK, neutrophil,
monocyte) per sample; the monocyte series adds carotid plaque and coronary artery calcium (CAC).

Per cell type:

1. **The shield** as everywhere in this project (``tissue_contrast``): per-CpG age slopes of entropy, dispersion
   and |mean change| with covariates (race-sex-site, contamination scores), set CpGs against look-alikes matched on
   island status and starting level, in SDs of the metric.
2. **Clocks**: each clock's correlation with age, and the share of its acceleration (residual on age) that the
   residual contamination scores still explain.
3. **Health** (monocytes): plaque count and log(1 + CAC) on each person's ``F``, ``B``, ``S`` (leave-one-out,
   ``health_states``) and clock accelerations, with age and race-sex-site; effects per SD.

``entropy-age figures sorted_cells`` (as a Slurm job; series matrices and ``*_methylome_normalized.txt.gz``
under ``data/input/geo``).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.shield_clocks import linear_score, residual
from entropy_age.experiments.syndromes import matrix_path
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    Prepared,
    _features,
    _lifespan_metrics,
    _matched_sd,
    _with_strata,
)

SERIES = {"CD14+ monocytes": "GSE56046", "CD4+ T cells": "GSE56581"}
CONTAMINATION = ("bcell", "tcell", "nkcell", "neutro", "mono")
CLOCKS_USED = ("Horvath 2013", "Hannum", "PhenoAge", "DunedinPACE", "Zhang 2019", "StocH", "StocP", "StocZ")
GROUPS = ("core_essential", "housekeeping")
PERMUTATIONS = 1000


def normalized_path(paths, gse: str) -> Path:
    return paths.input / "geo" / gse / f"{gse}_methylome_normalized.txt.gz"


def read_people(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM), ``id`` (the participant id the methylation table uses), age, codes and scores."""
    chars = _characteristics(matrix)
    titles = []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [c.strip().strip('"').split("_")[0] for c in line.rstrip("\n").split("\t")[1:]]
                break
    age = next(c for c in chars.columns if c.startswith("age"))
    present = [c for c in (*CONTAMINATION, "plaque", "cac") if c in chars.columns]
    return chars.with_columns(id=pl.Series(titles)).select(
        "sample",
        "id",
        pl.col(age).cast(pl.Float64, strict=False).alias("age"),
        pl.col("racegendersite").alias("group"),
        *(pl.col(c).cast(pl.Float64, strict=False) for c in present),
    )


def read_betas(path: Path, features: list[str], ids: list[str]) -> np.ndarray:
    """Betas (from M-values) of ``features`` for participants ``ids`` from the normalized table; NaN where absent."""
    index = {f: j for j, f in enumerate(features)}
    x = np.full((len(ids), len(features)), np.nan)
    columns = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            if columns is None:
                header = [c.strip('"') for c in cells]
                where = {h.removesuffix(".Mvalue"): i for i, h in enumerate(header) if h.endswith(".Mvalue")}
                columns = np.array([where[i] for i in ids])
                continue
            j = index.get(cells[0].strip('"'))
            if j is None:
                continue
            m = np.array([float(cells[c]) if cells[c] not in ("", "NA") else np.nan for c in columns])
            x[:, j] = 2.0**m / (1.0 + 2.0**m)
    return x


def covariates(people: pl.DataFrame) -> np.ndarray:
    scores = [people[c].fill_null(people[c].median()).to_numpy()[:, None] for c in CONTAMINATION if c in people.columns]
    return np.column_stack([dummies(people["group"].to_numpy()), *scores])


def clock_rows(cell: str, xc, clock_cpgs, clocks, people: pl.DataFrame) -> tuple[list[dict], dict[str, np.ndarray]]:
    age = people["age"].to_numpy()
    contamination = np.column_stack(
        [people[c].fill_null(people[c].median()) for c in CONTAMINATION if c in people.columns]
    )
    base = np.column_stack([np.ones(age.size), age, dummies(people["group"].to_numpy())])
    rows, accel = [], {}
    for name, weights in clocks.items():
        score, used = linear_score(xc, clock_cpgs, weights)
        if not used:
            continue
        acc = residual(score, base)
        explained = 1 - np.var(residual(acc, np.column_stack([np.ones(age.size), contamination]))) / np.var(acc)
        rows.append(
            {
                "cell": cell,
                "clock": name,
                "cpgs_used": used,
                "r_age": float(np.corrcoef(score, age)[0, 1]),
                "contamination_r2": float(explained),
            }
        )
        accel[name] = acc
    return rows, accel


def health_rows(people: pl.DataFrame, measures: dict[str, np.ndarray]) -> list[dict]:
    rows = []
    design0 = np.column_stack([people["age"].to_numpy(), dummies(people["group"].to_numpy())])
    outcomes = {"plaque": people["plaque"].to_numpy(), "log CAC": np.log1p(people["cac"].to_numpy())}
    for outcome, y in outcomes.items():
        has = ~np.isnan(y)
        for name, v in measures.items():
            z = (v - v.mean()) / v.std(ddof=1)
            design = np.column_stack([np.ones(has.sum()), z[has], design0[has]])
            coef, se, p = ols(y[has], design, 1)
            rows.append(
                {"outcome": outcome, "measure": name, "effect_per_sd": coef, "se": se, "p": p, "n": int(has.sum())}
            )
    return rows


def run(ctx: Context) -> list[Path]:
    missing = [g for g in SERIES.values() if not normalized_path(ctx.paths, g).exists()]
    if missing:
        print(f"sorted_cells: normalized tables missing for {missing}. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("sorted_cells")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    lay = ClockLayout(ctx.paths)
    clocks = {c: read_clock(lay.file(c)) for c in CLOCKS_USED if lay.file(c).exists()}
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    shield, clock_table, health = [], [], []
    for cell, gse in SERIES.items():
        people = read_people(matrix_path(ctx.paths, gse)).drop_nulls("age")
        ids = people["id"].to_list()
        wanted = sorted(set(features) | set(clock_cpgs))
        x_all = read_betas(normalized_path(ctx.paths, gse), wanted, ids)
        col = {f: j for j, f in enumerate(wanted)}
        x = x_all[:, [col[f] for f in features]]
        ok = np.isnan(x).mean(axis=0) == 0
        kept = [f for f, k in zip(features, ok, strict=True) if k]
        cov = covariates(people)
        prepared = Prepared(cell, x[:, ok], people["age"].to_numpy(), cov, kept, np.array(ids), cov)
        metrics = _with_strata(_lifespan_metrics(prepared), probes)
        table = _matched_sd(
            metrics, membership.filter(pl.col("feature").is_in(kept)), LIFESPAN_METRICS, PERMUTATIONS, ctx.seed
        )
        shield.append(table.with_columns(cell=pl.lit(cell), people=pl.lit(people.height)))
        print(f"sorted_cells: {cell}: {people.height} people, shield done", file=sys.stderr, flush=True)
        xc = x_all[:, [col[f] for f in clock_cpgs]]
        rows, accel = clock_rows(cell, xc, clock_cpgs, clocks, people)
        clock_table += rows
        if "plaque" in people.columns:
            scores = loo_scores(x[:, ok], _islands(kept, probes), _sets(kept, membership))
            measures = {f"{m} {g}": scores[f"{m}_{g}"] for g in GROUPS for m in ("F", "B", "S") if f"{m}_{g}" in scores}
            health += [{"cell": cell, **r} for r in health_rows(people, measures | accel)]
            print(f"sorted_cells: {cell}: health done", file=sys.stderr, flush=True)
    write_table(pl.concat(shield), out / "shield")
    write_table(pl.DataFrame(clock_table), out / "clocks")
    if health:
        write_table(pl.DataFrame(health), out / "health")
    return [_figure(ctx, pl.concat(shield), pl.DataFrame(clock_table))]


def _figure(ctx: Context, shield: pl.DataFrame, clocks: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    ax = axes[0]
    cells = shield["cell"].unique(maintain_order=True).to_list()
    labels = []
    for y, (metric, group) in enumerate((m, g) for m in LIFESPAN_METRICS for g in GROUPS):
        for i, cell in enumerate(cells):
            r = shield.filter((pl.col("cell") == cell) & (pl.col("metric") == metric) & (pl.col("group") == group))
            if r.height:
                ax.plot(
                    r["matched_diff_sd"][0], y + 0.2 * i, "o", color=plots.PALETTE[i], label=cell if y == 0 else None
                )
        labels.append((y + 0.1, f"{metric.replace('_', ' ')} · {group.replace('_', ' ')}"))
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks([p for p, _ in labels], [t for _, t in labels], fontsize=7.5)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.set(xlabel="set minus look-alikes (SD); < 0: shielded", title="The shield in one cell type")
    plots.label_panel(ax, "a")
    ax = axes[1]
    for i, cell in enumerate(cells):
        c = clocks.filter(pl.col("cell") == cell)
        ax.scatter(c["r_age"], c["contamination_r2"], color=plots.PALETTE[i], s=18, label=cell)
        for name, xv, yv in c.select("clock", "r_age", "contamination_r2").iter_rows():
            ax.annotate(name, (xv, yv), fontsize=6)
    ax.set(
        xlabel="clock r with age",
        ylabel="acceleration share explained by contamination",
        title="Clocks in sorted cells",
    )
    plots.label_panel(ax, "b")
    fig.suptitle(
        "Sorted blood cells (MESA): one cell type, no changing mix",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("sorted_cells"))
