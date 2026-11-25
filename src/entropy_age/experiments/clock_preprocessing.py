"""Preprocessing moves clocks (RESEARCH_PLAN, ideas to pivot to, idea 5): every scored clock on the same arrays under
two SeSAMe tracks, series by series.

Tracks, both scored by scripts/clocks_biolearn.py at the pinned biolearn commit with the same coverage rule (90 % of a
clock's absolute weight measured in 90 % of arrays):

* masked: the standard rebuild, SeSAMe with the design/quality mask, ``data/interim/clocks/series/<GSE>_clocks.csv``;
  for EPIC-Italy (GSE51032), ``data/interim/clocks/GSE51032_qc712_clocks.csv``: the masked rebuild rescored on the
  712 arrays that pass the current QC (detection and the quarter-intensity rule), with the same metadata as the
  unmasked track, so GrimAge is scored on both;
* unmasked: the clock track of protocol change 4, SeSAMe CDPB without the Q mask,
  ``data/interim/clocks/series_unmasked/<GSE>_clocks.csv``; for EPIC-Italy ``GSE51032_clocks_unmasked.csv``, the
  same arrays scored with the metadata file (GrimAge).

No clock file scored from the deposited GEO betas exists yet, so the deposited track is left out until one is scored.

The arrays compared are, per series, those both tracks exported that pass the current QC (``reprocessed.read_qc``;
every array both exported when the series has no rebuild QC). Coverage is judged on exactly those arrays on both
tracks: ``write_manifest`` lists them, scripts/clock_coverage.py applies clocks_biolearn.py's coverage rule to each
track's betas on them, and a clock's status on a track is that result (``kept`` or ``coverage``); the status each
track's own panel gave on its own arrays is kept beside it (``*_own_arrays``). Reasons a clock was not scored at all
(model type, GrimAge without metadata, a failed prediction) come from the panels. Without the common-array coverage
file, the panels' statuses are used and ``status_basis`` says so.

Per series and clock, on the common arrays: the number of arrays, Pearson r between the tracks, the mean difference
(unmasked minus masked) in the clock's own units, the SD of the difference, and both divided by the clock's SD across
those arrays (the mean of the two tracks' SDs), so that clocks in years, rates and risks share one scale. A clock that
passes coverage on one track and fails it on the other is flagged. Clocks whose output is years get the mean absolute
error against chronological age on each track, on the common arrays with an age; ages come from the series' config
tables (``config/series``) and, for EPIC-Italy, from ``data/interim/clocks/GSE51032_meta.csv``. These MAEs are on all
compared arrays (repeat arrays of a donor and every tissue included); they describe the shift between tracks and are
no measure of a clock's age accuracy, which needs donor-balanced, tissue-valid subsets. Per clock across series: the
median r and the median absolute mean difference, raw and per SD.

Writes ``clock_preprocessing/pairs`` and ``clock_preprocessing/clocks`` and a figure: per clock, r across series and
the absolute mean difference per SD, clocks sorted by median r.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.series import series_table, specs

MASKED = Path("data/interim/clocks/series")
UNMASKED = Path("data/interim/clocks/series_unmasked")
#: series whose scores or betas live outside ``MASKED`` / ``UNMASKED``: EPIC-Italy, scored before the per-series runs
#: and rescored on the current QC-passing arrays with its metadata (GrimAge)
MASKED_ELSEWHERE = {"GSE51032": Path("data/interim/clocks/GSE51032_qc712_clocks.csv")}
UNMASKED_ELSEWHERE = {"GSE51032": Path("data/interim/clocks/GSE51032_clocks_unmasked.csv")}
MASKED_BETAS_ELSEWHERE = {"GSE51032": Path("data/interim/clocks/GSE51032_cpgs.parquet")}
#: coverage of both tracks on the common QC-passing arrays (scripts/clock_coverage.py) and the manifest it reads
COMMON = Path("data/interim/clocks/common")
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
    """Clock -> panel status (kept, coverage, or why not scored); the panel (``<stem>_panel.csv``, as
    clocks_biolearn.py names it) may still carry a ``.tmp`` suffix."""
    base = scores_path.with_name(scores_path.name.removesuffix(".csv") + "_panel.csv")
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
    *,
    arrays: set[str] | None = None,
    common_masked: dict[str, str] | None = None,
    common_unmasked: dict[str, str] | None = None,
) -> list[dict]:
    """One row per clock seen on either track: agreement between the tracks on the arrays both scored (only those in
    ``arrays`` when given). ``status_*`` are the panels' statuses on each track's own arrays; ``common_*`` the
    coverage status (kept or coverage) on the common arrays, which replaces a panel's kept or coverage when given."""
    status_masked, status_unmasked = status_masked or {}, status_unmasked or {}
    common = common_masked is not None and common_unmasked is not None
    output = dict(pl.concat([masked, unmasked]).select("clock", "output").unique("clock", keep="first").iter_rows())
    clocks = sorted(
        set(masked["clock"])
        | set(unmasked["clock"])
        | {c for c, s in status_masked.items() if s in ("kept", "coverage")}
        | {c for c, s in status_unmasked.items() if s in ("kept", "coverage")}
    )
    pair = masked.select("sample", "clock", pl.col("prediction").alias("m")).join(
        unmasked.select("sample", "clock", pl.col("prediction").alias("u")),
        on=["sample", "clock"],
        how="inner",
        validate="1:1",
    )
    if arrays is not None:
        pair = pair.filter(pl.col("sample").is_in(sorted(arrays)))
    if ages is not None:
        pair = pair.join(
            ages.select("sample", pl.col("age").cast(pl.Float64)).unique("sample"), on="sample", how="left"
        )
    rows = []

    def status(own: str, on_common: dict[str, str] | None, clock: str) -> str:
        if on_common is None or own not in ("kept", "coverage"):
            return own  # not scored for another reason: the panel says why
        return on_common.get(clock, own)

    for clock in clocks:
        on_m, on_u = clock in set(masked["clock"]), clock in set(unmasked["clock"])
        own_m = status_masked.get(clock, "kept" if on_m else "absent")
        own_u = status_unmasked.get(clock, "kept" if on_u else "absent")
        sm, su = status(own_m, common_masked, clock), status(own_u, common_unmasked, clock)
        part = pair.filter((pl.col("clock") == clock) & pl.col("m").is_finite() & pl.col("u").is_finite())
        m, u = part["m"].to_numpy(), part["u"].to_numpy()
        row = {
            "clock": clock,
            "output": output.get(clock),
            "status_masked": sm,
            "status_unmasked": su,
            "status_masked_own_arrays": own_m,
            "status_unmasked_own_arrays": own_u,
            "status_basis": "common arrays" if common else "each track's own arrays (no common coverage file)",
            "one_track_only": {sm, su} == {"kept", "coverage"},
            "one_track_only_own_arrays": {own_m, own_u} == {"kept", "coverage"},
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
            (pl.col("one_track_only_own_arrays") & (pl.col("status_unmasked_own_arrays") == "kept"))
            .sum()
            .alias("unmasked_only_own_arrays"),
            (pl.col("one_track_only_own_arrays") & (pl.col("status_masked_own_arrays") == "kept"))
            .sum()
            .alias("masked_only_own_arrays"),
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


def tracks(ctx: Context) -> list[tuple[str, Path, Path, Path, Path]]:
    """(series, masked scores, unmasked scores, masked betas, unmasked betas) for every series of the clock track."""
    out = []
    for upath in sorted((ctx.paths.root / UNMASKED).glob("*_clocks.csv")):
        gse = upath.name.removesuffix("_clocks.csv")
        mpath = ctx.paths.root / MASKED_ELSEWHERE.get(gse, MASKED / upath.name)
        upath = ctx.paths.root / UNMASKED_ELSEWHERE.get(gse, UNMASKED / upath.name)
        mbetas = ctx.paths.root / MASKED_BETAS_ELSEWHERE.get(gse, MASKED / f"{gse}_cpgs.parquet")
        out.append((gse, mpath, upath, mbetas, ctx.paths.root / UNMASKED / f"{gse}_cpgs.parquet"))
    return out


def common_arrays(ctx: Context, gse: str, mbetas: Path, ubetas: Path) -> tuple[set[str] | None, str]:
    """The arrays both tracks exported that pass the current QC, and how they were chosen; None without both
    exports."""
    if not (mbetas.exists() and ubetas.exists()):
        return None, "no betas export on one track"
    both = (set(pl.read_parquet_schema(mbetas)) & set(pl.read_parquet_schema(ubetas))) - {"cpg"}
    if not reprocessed.available(ctx.paths, gse):
        return both, "exported on both tracks (no rebuild QC)"
    passed = set(reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))["sample"].cast(pl.String).to_list())
    return both & passed, "exported on both tracks, passing the current QC"


def write_manifest(ctx: Context) -> Path:
    """The manifest scripts/clock_coverage.py reads: per series and track, the betas and the common arrays."""
    rows = []
    for gse, _, _, mbetas, ubetas in tracks(ctx):
        arrays, _ = common_arrays(ctx, gse, mbetas, ubetas)
        if not arrays:
            continue
        joined = ";".join(sorted(arrays))
        rows += [
            {"series": gse, "track": "masked", "parquet": str(mbetas), "samples": joined},
            {"series": gse, "track": "unmasked", "parquet": str(ubetas), "samples": joined},
        ]
    path = ctx.paths.root / COMMON / "manifest.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(path)
    return path


def read_common(ctx: Context) -> dict[tuple[str, str], tuple[dict[str, str], int]]:
    """(series, track) -> (clock -> kept or coverage on the common arrays, number of arrays), when computed."""
    path = ctx.paths.root / COMMON / "coverage.csv"
    if not path.exists():
        return {}
    out: dict[tuple[str, str], tuple[dict[str, str], int]] = {}
    for (gse, track), part in pl.read_csv(path, schema_overrides={"series": pl.String}).group_by(
        "series", "track", maintain_order=True
    ):
        out[(gse, track)] = (dict(part.select("clock", "status").iter_rows()), int(part["arrays"].max()))
    return out


def run(ctx: Context) -> list[Path]:
    by_accession: dict[str, list[dict]] = {}
    for spec in specs().values():
        if spec.get("annotation", "human") == "human":
            by_accession.setdefault(spec["accession"], []).append(spec)
    rows, skipped = [], {}
    coverage = read_common(ctx)
    for gse, mpath, upath, mbetas, ubetas in tracks(ctx):
        masked, unmasked = read_scores(mpath), read_scores(upath)
        if masked is None or unmasked is None:
            skipped.setdefault("no scores on one track", []).append(gse)
            continue
        arrays, how = common_arrays(ctx, gse, mbetas, ubetas)
        cm, cu = coverage.get((gse, "masked")), coverage.get((gse, "unmasked"))
        if arrays is not None and cm is not None and cm[1] != len(arrays):
            print(f"clock_preprocessing: {gse}: common coverage on {cm[1]} arrays, now {len(arrays)}", file=sys.stderr)
            cm = cu = None  # stale: computed on other arrays
        ages = ages_of(ctx, gse, by_accession)
        found = compare(
            masked,
            unmasked,
            ages,
            read_status(mpath),
            read_status(upath),
            arrays=arrays,
            common_masked=cm[0] if cm and cu else None,
            common_unmasked=cu[0] if cm and cu else None,
        )
        for row in found:
            rows.append(
                {
                    "series": gse,
                    "masked_file": str(mpath.relative_to(ctx.paths.root)),
                    "unmasked_file": str(upath.relative_to(ctx.paths.root)),
                    "arrays_common": len(arrays) if arrays is not None else None,
                    "arrays_chosen": how,
                }
                | row
            )
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


if __name__ == "__main__":  # ``python -m entropy_age.experiments.clock_preprocessing manifest``: before clock_coverage
    from entropy_age.experiments.common import Scale
    from entropy_age.paths import ProjectPaths

    if sys.argv[1:] == ["manifest"]:
        print(write_manifest(Context(ProjectPaths.default(), Scale.of("paper"))))
