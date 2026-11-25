"""Preprocessing moves clocks (RESEARCH_PLAN, ideas to pivot to, idea 5): every scored clock on the same arrays under
two SeSAMe tracks, series by series.

Tracks, both scored by scripts/clocks_biolearn.py at the pinned biolearn commit with the same coverage rule (90 % of a
clock's absolute weight measured in 90 % of arrays):

* masked: the standard rebuild, SeSAMe with the design/quality mask, ``data/interim/clocks/series/<GSE>_clocks.csv``;
  for EPIC-Italy (GSE51032), which has no file there, ``data/interim/clocks/GSE51032_clocks.csv`` (R160, the same
  masked rebuild scored before the quarter-intensity QC rule, so its arrays are a superset);
* unmasked: the clock track of protocol change 4, SeSAMe CDPB without the Q mask,
  ``data/interim/clocks/series_unmasked/<GSE>_clocks.csv``.

No clock file scored from the deposited GEO betas exists yet (both EPIC-Italy files under ``data/interim/clocks`` are
rebuilds), so the deposited track is left out until one is scored.

Per series and clock, on the arrays scored on both tracks: the number of arrays, Pearson r between the tracks, the
mean difference (unmasked minus masked) in the clock's own units, the SD of the difference, and both divided by the
clock's SD across those arrays (the mean of the two tracks' SDs), so that clocks in years, rates and risks share one
scale. Each clock's panel status on each track (kept, coverage, or why it was not scored) is reported, and a clock that
passes coverage on one track and fails it on the other is flagged. Clocks whose output is years get the mean absolute
error against chronological age on each track, on the common arrays with an age; ages come from the series' config
tables (``config/series``) and, for EPIC-Italy, from ``data/interim/clocks/GSE51032_meta.csv``. Per clock across
series: the median r and the median absolute mean difference, raw and per SD.

Writes ``clock_preprocessing/pairs`` and ``clock_preprocessing/clocks`` and a figure: per clock, r across series and
the absolute mean difference per SD, clocks sorted by median r.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.series import series_table, specs

MASKED = Path("data/interim/clocks/series")
UNMASKED = Path("data/interim/clocks/series_unmasked")
#: series whose masked scores live outside ``MASKED`` (scored before the per-series runs)
MASKED_ELSEWHERE = {"GSE51032": Path("data/interim/clocks/GSE51032_clocks.csv")}
#: series without a config table whose ages sit beside their clock scores
AGES_ELSEWHERE = {"GSE51032": Path("data/interim/clocks/GSE51032_meta.csv")}
YEARS = "Age (Years)"
MIN_ARRAYS = 3


def read_scores(path: Path) -> pl.DataFrame | None:
    """A long clock table (sample, clock, output, prediction); None when the file is missing or empty."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        frame = pl.read_csv(path, schema_overrides={"sample": pl.String})
    except pl.exceptions.NoDataError:
        return None
    return frame.select("sample", "clock", "output", pl.col("prediction").cast(pl.Float64)) if frame.height else None


def read_status(scores_path: Path) -> dict[str, str]:
    """Clock -> panel status (kept, coverage, or why not scored); the panel may still carry a ``.tmp`` suffix."""
    base = scores_path.with_name(scores_path.name.replace("_clocks.csv", "_clocks_panel.csv"))
    for path in (base, base.with_name(base.name + ".tmp")):
        if path.exists() and path.stat().st_size:
            return dict(pl.read_csv(path).select("clock", "status").iter_rows())
    return {}


def compare(
    masked: pl.DataFrame,
    unmasked: pl.DataFrame,
    ages: pl.DataFrame | None = None,
    status_masked: dict[str, str] | None = None,
    status_unmasked: dict[str, str] | None = None,
) -> list[dict]:
    """One row per clock seen on either track: agreement between the tracks on the arrays both scored."""
    status_masked, status_unmasked = status_masked or {}, status_unmasked or {}
    output = dict(pl.concat([masked, unmasked]).select("clock", "output").unique("clock", keep="first").iter_rows())
    clocks = sorted(
        set(masked["clock"])
        | set(unmasked["clock"])
        | {c for c, s in status_masked.items() if s in ("kept", "coverage")}
        | {c for c, s in status_unmasked.items() if s in ("kept", "coverage")}
    )
    pair = masked.select("sample", "clock", pl.col("prediction").alias("m")).join(
        unmasked.select("sample", "clock", pl.col("prediction").alias("u")), on=["sample", "clock"], how="inner"
    )
    if ages is not None:
        pair = pair.join(
            ages.select("sample", pl.col("age").cast(pl.Float64)).unique("sample"), on="sample", how="left"
        )
    rows = []
    for clock in clocks:
        on_m, on_u = clock in set(masked["clock"]), clock in set(unmasked["clock"])
        sm = status_masked.get(clock, "kept" if on_m else "absent")
        su = status_unmasked.get(clock, "kept" if on_u else "absent")
        part = pair.filter((pl.col("clock") == clock) & pl.col("m").is_finite() & pl.col("u").is_finite())
        m, u = part["m"].to_numpy(), part["u"].to_numpy()
        row = {
            "clock": clock,
            "output": output.get(clock),
            "status_masked": sm,
            "status_unmasked": su,
            "one_track_only": {sm, su} == {"kept", "coverage"},
            "arrays": len(m),
            "r": None,
            "mean_diff": None,
            "sd_diff": None,
            "sd_clock": None,
            "mean_diff_per_sd": None,
            "sd_diff_per_sd": None,
            "arrays_age": None,
            "mae_masked": None,
            "mae_unmasked": None,
        }
        if len(m) >= MIN_ARRAYS:
            d = u - m
            sd = float((np.std(m, ddof=1) + np.std(u, ddof=1)) / 2)
            row |= {"mean_diff": float(d.mean()), "sd_diff": float(np.std(d, ddof=1)), "sd_clock": sd}
            if np.std(m) > 0 and np.std(u) > 0:
                row["r"] = float(np.corrcoef(m, u)[0, 1])
            if sd > 0:
                row |= {"mean_diff_per_sd": row["mean_diff"] / sd, "sd_diff_per_sd": row["sd_diff"] / sd}
        if output.get(clock) == YEARS and "age" in part.columns:
            aged = part.drop_nulls("age")
            if aged.height >= MIN_ARRAYS:
                age = aged["age"].to_numpy()
                row |= {
                    "arrays_age": aged.height,
                    "mae_masked": float(np.mean(np.abs(aged["m"].to_numpy() - age))),
                    "mae_unmasked": float(np.mean(np.abs(aged["u"].to_numpy() - age))),
                }
        rows.append(row)
    return rows


def summarise(pairs: pl.DataFrame) -> pl.DataFrame:
    """Per clock across series: how many compared, median and minimum r, median absolute mean difference."""
    return (
        pairs.group_by("clock")
        .agg(
            pl.col("output").drop_nulls().first(),
            pl.col("r").is_not_null().sum().alias("series_compared"),
            pl.col("r").median().alias("median_r"),
            pl.col("r").min().alias("min_r"),
            pl.col("mean_diff").abs().median().alias("median_abs_mean_diff"),
            pl.col("mean_diff_per_sd").abs().median().alias("median_abs_mean_diff_per_sd"),
            pl.col("sd_diff_per_sd").median().alias("median_sd_diff_per_sd"),
            (pl.col("one_track_only") & (pl.col("status_unmasked") == "kept")).sum().alias("unmasked_only"),
            (pl.col("one_track_only") & (pl.col("status_masked") == "kept")).sum().alias("masked_only"),
            (pl.col("mae_unmasked") - pl.col("mae_masked")).median().alias("median_mae_change"),
        )
        .sort("median_r", "clock", nulls_last=True)
    )


def ages_of(ctx: Context, gse: str, by_accession: dict[str, list[dict]]) -> pl.DataFrame | None:
    frames = []
    for spec in by_accession.get(gse, []):
        table = series_table(ctx.paths, spec)
        if table is not None:
            frames.append(table.select(pl.col("sample").cast(pl.String), "age"))
    if gse in AGES_ELSEWHERE and (path := ctx.paths.root / AGES_ELSEWHERE[gse]).exists():
        frames.append(pl.read_csv(path).select(pl.col("sample").cast(pl.String), pl.col("age").cast(pl.Float64)))
    return pl.concat(frames).drop_nulls().unique("sample", keep="first") if frames else None


def run(ctx: Context) -> list[Path]:
    by_accession: dict[str, list[dict]] = {}
    for spec in specs().values():
        if spec.get("annotation", "human") == "human":
            by_accession.setdefault(spec["accession"], []).append(spec)
    rows, skipped = [], {}
    for upath in sorted((ctx.paths.root / UNMASKED).glob("*_clocks.csv")):
        gse = upath.name.removesuffix("_clocks.csv")
        mpath = ctx.paths.root / MASKED_ELSEWHERE.get(gse, MASKED / upath.name)
        masked, unmasked = read_scores(mpath), read_scores(upath)
        if masked is None or unmasked is None:
            skipped.setdefault("no scores on one track", []).append(gse)
            continue
        ages = ages_of(ctx, gse, by_accession)
        for row in compare(masked, unmasked, ages, read_status(mpath), read_status(upath)):
            rows.append({"series": gse, "masked_file": str(mpath.relative_to(ctx.paths.root))} | row)
    for reason, names in sorted(skipped.items()):
        print(f"clock_preprocessing: {len(names)} series skipped, {reason}: {', '.join(names)}", file=sys.stderr)
    if not rows:
        print("clock_preprocessing: no series scored on both tracks", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_preprocessing")
    pairs = pl.DataFrame(rows, infer_schema_length=None).sort("series", "clock")
    clocks = summarise(pairs)
    write_table(pairs, out / "pairs", csv_max_rows=10_000)
    write_table(clocks, out / "clocks")
    return [draw(pairs, clocks, ctx.figure_path("clock_preprocessing"))]


def draw(pairs: pl.DataFrame, clocks: pl.DataFrame, path: Path) -> Path:
    order = clocks.filter(pl.col("series_compared") > 0).sort("median_r", descending=True)["clock"].to_list()
    plots.apply_style()
    fig, axes = plots.plt.subplots(1, 2, figsize=(9, 0.2 * len(order) + 1.6), sharey=True)
    for ax, column, label in (
        (axes[0], "r", "r between tracks (per series)"),
        (axes[1], "mean_diff_per_sd", "|mean difference| / clock SD"),
    ):
        for i, clock in enumerate(order):
            values = pairs.filter((pl.col("clock") == clock) & pl.col(column).is_not_null())[column].to_numpy()
            values = np.abs(values) if column != "r" else values
            jitter = (np.arange(len(values)) % 5 - 2) * 0.06
            ax.scatter(values, i + jitter, s=7, color=plots.POINT, alpha=0.6, lw=0)
            if len(values):
                ax.plot([np.median(values)] * 2, [i - 0.35, i + 0.35], color=plots.ACCENT, lw=1.4)
        ax.set_xlabel(label)
    axes[0].set_xlim(max(-0.2, float(np.nanmin(pairs["r"].drop_nulls().to_numpy())) - 0.05), 1.01)
    axes[1].set_xscale("symlog", linthresh=0.05)
    axes[1].set_xlim(left=0)
    axes[0].set_yticks(range(len(order)), order, fontsize=6)
    axes[0].invert_yaxis()
    fig.suptitle(
        "Clock predictions on the same arrays: SeSAMe with the design mask against without (one point per series)",
        fontsize=9,
        x=0.02,
        ha="left",
    )
    fig.tight_layout()
    return plots.save(fig, path)
