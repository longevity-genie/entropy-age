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
(model type, GrimAge without metadata, a failed prediction) come from the panels. A coverage result is reused only
when it was computed on exactly the current common arrays (a SHA-256 of their sorted ids) and on both beta exports
as they are now (path, size and modification time): a rebuild that keeps the array ids still changes the export, so
its old result is stale. Stale or missing results are recomputed (``refresh_common`` runs scripts/clock_coverage.py
on the affected series); when that fails, the common-array status is ``unknown``, ``one_track_only`` is empty and
``status_basis`` gives the reason, so such rows never enter the one-track-only counts.

Per series and clock, on the common arrays: the number of arrays, Pearson r between the tracks, the mean difference
(unmasked minus masked) in the clock's own units, the SD of the difference, and both divided by the clock's SD across
those arrays (the mean of the two tracks' SDs), so that clocks in years, rates and risks share one scale. A clock that
passes coverage on one track and fails it on the other is flagged. Clocks whose output is years get the mean absolute
error against chronological age on each track, on the common arrays with an age (an output whose unit, compared
without case, ends in ``Age (Years)`` and whose target is no mortality-adjusted age: Bocklandt's ``Age (years)`` and
the cortical clock's ``Human Cortex Age (Years)`` count, GrimAge does not); ages come from the series' config
tables (``config/series``) and, for EPIC-Italy, from ``data/interim/clocks/GSE51032_meta.csv``. These MAEs are on all
compared arrays (repeat arrays of a donor and every tissue included); they describe the shift between tracks and are
no measure of a clock's age accuracy, which needs donor-balanced, tissue-valid subsets. Per clock across series: the
median r and the median absolute mean difference, raw and per SD.

Diagnostic (``clock_preprocessing/forced``): a clock that fails coverage on the masked track can be scored with the
coverage rule switched off (scripts/clocks_forced.py, files named ``*_below_coverage_clocks.csv``). Each such file is
set beside the unmasked track on the common arrays: paired r, the mean shift (unmasked minus masked), r after age and
sex where known, and both tracks' coverage results. These rows are labelled ``diagnostic`` and never enter ``pairs``.

Writes ``clock_preprocessing/pairs``, ``clock_preprocessing/clocks`` and ``clock_preprocessing/forced`` and a figure:
per clock, r across series and the absolute mean difference per SD, clocks sorted by median r.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
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
#: scores of a masked clock below coverage, scored for comparison only (scripts/clocks_forced.py)
FORCED = {"GSE51032": Path("data/interim/clocks/GSE51032_qc712_grimage_below_coverage_clocks.csv")}
FORCED_LABEL = "diagnostic: masked track scored below coverage (scripts/clocks_forced.py), never a panel score"
COVERAGE_SCRIPT = Path("scripts/clock_coverage.py")
MIN_ARRAYS = 3


def is_age_years(output: str | None) -> bool:
    """A clock output of chronological-age years: the unit, compared without case, ends in ``age (years)`` and the
    target is no mortality-adjusted age (GrimAge's ``Mortality Adjusted Age (Years)``)."""
    unit = " ".join(str(output or "").lower().split())
    return unit.endswith("age (years)") and "mortality" not in unit


def samples_digest(samples) -> str:
    """SHA-256 of the sorted array ids joined by ``;`` (as scripts/clock_coverage.py computes it)."""
    return hashlib.sha256(";".join(sorted(str(s) for s in samples)).encode()).hexdigest()


def fingerprint(path: Path) -> str:
    """Size and modification time of a betas export (as scripts/clock_coverage.py records it)."""
    st = Path(path).stat()
    return f"{st.st_size}:{st.st_mtime_ns}"


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
    unknown: str | None = None,
) -> list[dict]:
    """One row per clock seen on either track: agreement between the tracks on the arrays both scored (only those in
    ``arrays`` when given). ``status_*`` are the panels' statuses on each track's own arrays; ``common_*`` the
    coverage status (kept or coverage) on the common arrays, which replaces a panel's kept or coverage when given.
    ``unknown`` (a reason) says the common-array coverage could not be validated: a scored clock's status is then
    ``unknown`` and ``one_track_only`` empty."""
    status_masked, status_unmasked = status_masked or {}, status_unmasked or {}
    common = common_masked is not None and common_unmasked is not None and unknown is None
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
        if own not in ("kept", "coverage"):
            return own  # not scored for another reason: the panel says why
        if unknown is not None:
            return "unknown"
        if on_common is None:
            return own
        return on_common.get(clock, own)

    if unknown is not None:
        basis = f"unknown: {unknown}"
    elif common:
        basis = "common arrays"
    else:
        basis = "each track's own arrays (no common coverage file)"

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
            "status_basis": basis,
            "one_track_only": None if "unknown" in (sm, su) else {sm, su} == {"kept", "coverage"},
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
        if is_age_years(output.get(clock)) and "age" in part.columns:
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
            pl.col("status_basis").str.starts_with("unknown").sum().alias("series_coverage_unknown"),
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


def write_manifest(ctx: Context, only: set[str] | None = None, path: Path | None = None) -> Path:
    """The manifest scripts/clock_coverage.py reads: per series and track, the betas, the common arrays, their
    digest and the export's fingerprint (``only``: these series)."""
    rows = []
    for gse, _, _, mbetas, ubetas in tracks(ctx):
        if only is not None and gse not in only:
            continue
        arrays, _ = common_arrays(ctx, gse, mbetas, ubetas)
        if not arrays:
            continue
        joined, digest = ";".join(sorted(arrays)), samples_digest(arrays)
        for track, betas in (("masked", mbetas), ("unmasked", ubetas)):
            rows.append(
                {
                    "series": gse,
                    "track": track,
                    "parquet": str(betas),
                    "samples": joined,
                    "samples_sha256": digest,
                    "parquet_fingerprint": fingerprint(betas),
                }
            )
    path = path or ctx.paths.root / COMMON / "manifest.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(path)
    return path


def read_common(ctx: Context, path: Path | None = None) -> dict[tuple[str, str], dict]:
    """(series, track) -> the coverage record on the common arrays: ``status`` (clock -> kept or coverage),
    ``arrays``, and, when scripts/clock_coverage.py recorded them, ``samples_sha256``, ``parquet`` and
    ``parquet_fingerprint``."""
    path = path or ctx.paths.root / COMMON / "coverage.csv"
    if not path.exists() or path.stat().st_size == 0:
        return {}
    table = pl.read_csv(path, schema_overrides={"series": pl.String}, infer_schema_length=None)
    out: dict[tuple[str, str], dict] = {}
    for (gse, track), part in table.group_by("series", "track", maintain_order=True):
        record = {"status": dict(part.select("clock", "status").iter_rows()), "arrays": int(part["arrays"].max())}
        for column in ("samples_sha256", "parquet", "parquet_fingerprint"):
            values = part[column].unique().to_list() if column in part.columns else []
            record[column] = values[0] if len(values) == 1 else None
        out[(gse, track)] = record
    return out


def validate_common(record: dict | None, arrays: set[str] | None, betas: Path) -> str | None:
    """Why a coverage record cannot stand for ``arrays`` on ``betas`` as they are now, or None when it can."""
    if record is None:
        return "no coverage on the common arrays"
    if arrays is None:
        return "no common arrays"
    if not record.get("samples_sha256") or not record.get("parquet_fingerprint"):
        return "coverage recorded without array ids or export fingerprint"
    if record["samples_sha256"] != samples_digest(arrays):
        return f"coverage computed on other arrays ({record['arrays']} then, {len(arrays)} now)"
    if record.get("parquet") and Path(record["parquet"]).resolve() != Path(betas).resolve():
        return "coverage computed on another export"
    if not Path(betas).exists():
        return "betas export missing"
    if record["parquet_fingerprint"] != fingerprint(betas):
        return "betas export changed since the coverage was computed"
    return None


def series_common(
    coverage: dict[tuple[str, str], dict], gse: str, arrays: set[str] | None, mbetas: Path, ubetas: Path
) -> tuple[dict[str, str] | None, dict[str, str] | None, str | None]:
    """(masked statuses, unmasked statuses, None) when both tracks' records hold for the current arrays and
    exports, else (None, None, the reason)."""
    reasons = []
    for track, betas in (("masked", mbetas), ("unmasked", ubetas)):
        why = validate_common(coverage.get((gse, track)), arrays, betas)
        if why is not None:
            reasons.append(f"{track}: {why}")
    if reasons:
        return None, None, "; ".join(reasons)
    return coverage[(gse, "masked")]["status"], coverage[(gse, "unmasked")]["status"], None


def refresh_common(ctx: Context, stale: set[str]) -> bool:
    """Recompute the common-array coverage of the ``stale`` series with scripts/clock_coverage.py and merge it into
    the coverage file; False when the script cannot run (``ENTROPY_AGE_COVERAGE_REFRESH=0`` skips it)."""
    if not stale or os.environ.get("ENTROPY_AGE_COVERAGE_REFRESH", "1") == "0":
        return False
    common = ctx.paths.root / COMMON
    manifest = write_manifest(ctx, only=stale, path=common / "manifest_refresh.csv")
    fresh = common / "coverage_refresh.csv"
    print(f"clock_preprocessing: recomputing common coverage for {', '.join(sorted(stale))}", file=sys.stderr)
    cmd = ["uv", "run", "--script", str(ctx.paths.root / COVERAGE_SCRIPT), str(manifest), str(fresh)]
    try:
        subprocess.run(cmd, check=True, cwd=ctx.paths.root)
    except (OSError, subprocess.CalledProcessError) as err:
        print(f"clock_preprocessing: coverage refresh failed: {err}", file=sys.stderr)
        return False
    new = pl.read_csv(fresh, schema_overrides={"series": pl.String}, infer_schema_length=None)
    path = common / "coverage.csv"
    if path.exists() and path.stat().st_size:
        old = pl.read_csv(path, schema_overrides={"series": pl.String}, infer_schema_length=None)
        new = pl.concat([old.filter(~pl.col("series").is_in(new["series"].unique().to_list())), new], how="diagonal")
    new.sort("series", "track", maintain_order=True).write_csv(path)
    fresh.unlink()
    return True


def forced_rows(ctx: Context, common: dict[tuple[str, str], dict]) -> list[dict]:
    """The labelled below-coverage diagnostic: per forced scores file and clock, the masked track scored with the
    coverage rule off against the unmasked track on the common arrays."""
    rows = []
    by_gse = {gse: (mpath, upath, mbetas, ubetas) for gse, mpath, upath, mbetas, ubetas in tracks(ctx)}
    for gse, rel in sorted(FORCED.items()):
        path = ctx.paths.root / rel
        forced = read_scores(path)
        if forced is None or gse not in by_gse:
            continue
        mpath, upath, mbetas, ubetas = by_gse[gse]
        unmasked = read_scores(upath)
        if unmasked is None:
            continue
        arrays, how = common_arrays(ctx, gse, mbetas, ubetas)
        cm, cu, unknown = series_common(common, gse, arrays, mbetas, ubetas)
        own_m, own_u = read_status(mpath), read_status(upath)
        meta = None
        if gse in AGES_ELSEWHERE and (mp := ctx.paths.root / AGES_ELSEWHERE[gse]).exists():
            meta = pl.read_csv(mp).select(
                pl.col("sample").cast(pl.String), pl.col("age").cast(pl.Float64), pl.col("female").cast(pl.Float64)
            )
        share = None
        raw = pl.read_csv(path, schema_overrides={"sample": pl.String})
        for clock in sorted(set(forced["clock"])):
            pair = (
                forced.filter(pl.col("clock") == clock)
                .select("sample", pl.col("prediction").alias("m"))
                .join(
                    unmasked.filter(pl.col("clock") == clock).select("sample", pl.col("prediction").alias("u")),
                    on="sample",
                    validate="1:1",
                )
                .filter(pl.col("m").is_finite() & pl.col("u").is_finite())
            )
            if arrays is not None:
                pair = pair.filter(pl.col("sample").is_in(sorted(arrays)))
            m, u = pair["m"].to_numpy(), pair["u"].to_numpy()
            if "measured_weight_share" in raw.columns:
                part = raw.filter((pl.col("clock") == clock) & pl.col("sample").is_in(pair["sample"].to_list()))
                share = part["measured_weight_share"].median()
            row = {
                "series": gse,
                "clock": clock,
                "label": FORCED_LABEL,
                "forced_file": str(rel),
                "unmasked_file": str(upath.relative_to(ctx.paths.root)),
                "arrays_chosen": how,
                "arrays": len(m),
                "masked_status_common": "unknown" if unknown else (cm or {}).get(clock),
                "unmasked_status_common": "unknown" if unknown else (cu or {}).get(clock),
                "status_basis": f"unknown: {unknown}" if unknown else "common arrays",
                "masked_status_own_arrays": own_m.get(clock),
                "unmasked_status_own_arrays": own_u.get(clock),
                "masked_median_weight_share": share,
                "r": None,
                "mean_diff": None,
                "sd_diff": None,
                "partial_r_age_sex": None,
            }
            if len(m) >= MIN_ARRAYS and np.std(m) > 0 and np.std(u) > 0:
                d = u - m
                row |= {
                    "r": float(np.corrcoef(m, u)[0, 1]),
                    "mean_diff": float(d.mean()),
                    "sd_diff": float(np.std(d, ddof=1)),
                }
                if meta is not None:
                    held = pair.join(meta, on="sample", how="left").drop_nulls(["age", "female"])
                    if held.height >= MIN_ARRAYS + 3:
                        x = np.column_stack([np.ones(held.height), held["age"].to_numpy(), held["female"].to_numpy()])
                        rm = held["m"].to_numpy() - x @ np.linalg.lstsq(x, held["m"].to_numpy(), rcond=None)[0]
                        ru = held["u"].to_numpy() - x @ np.linalg.lstsq(x, held["u"].to_numpy(), rcond=None)[0]
                        if np.std(rm) > 0 and np.std(ru) > 0:
                            row["partial_r_age_sex"] = float(np.corrcoef(rm, ru)[0, 1])
            rows.append(row)
    return rows


def run(ctx: Context) -> list[Path]:
    by_accession: dict[str, list[dict]] = {}
    for spec in specs().values():
        if spec.get("annotation", "human") == "human":
            by_accession.setdefault(spec["accession"], []).append(spec)
    rows, skipped = [], {}
    coverage = read_common(ctx)
    found_tracks = tracks(ctx)
    stale = set()
    for gse, _, _, mbetas, ubetas in found_tracks:
        arrays, _ = common_arrays(ctx, gse, mbetas, ubetas)
        if arrays and series_common(coverage, gse, arrays, mbetas, ubetas)[2] is not None:
            stale.add(gse)
    if stale and refresh_common(ctx, stale):
        coverage = read_common(ctx)
    for gse, mpath, upath, mbetas, ubetas in found_tracks:
        masked, unmasked = read_scores(mpath), read_scores(upath)
        if masked is None or unmasked is None:
            skipped.setdefault("no scores on one track", []).append(gse)
            continue
        arrays, how = common_arrays(ctx, gse, mbetas, ubetas)
        cm, cu, unknown = series_common(coverage, gse, arrays, mbetas, ubetas)
        if unknown is not None:
            print(f"clock_preprocessing: {gse}: common-array coverage unknown ({unknown})", file=sys.stderr)
        ages = ages_of(ctx, gse, by_accession)
        found = compare(
            masked,
            unmasked,
            ages,
            read_status(mpath),
            read_status(upath),
            arrays=arrays,
            common_masked=cm,
            common_unmasked=cu,
            unknown=unknown,
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
    if forced := forced_rows(ctx, coverage):
        write_table(pl.DataFrame(forced, infer_schema_length=None), out / "forced")
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
