"""Read an atlas export back and check it: ``uv run python scripts/atlas_validate.py OUT_DIR [--privacy]``.

Checks the manifest (sizes, sha256, row counts), that every file the card names exists, unique stable ids, that
every row of the other tables joins to a series row, ages in range, licences present, public sample ids only, and
counts against the QC tables (``reprocessed.read_qc``) and against ``clock_drift.parquet``. ``--privacy`` also prints
how many included people are unique on series, group, sex and age at each age precision. Exit status 1 on a failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))

import atlas_export as ax  # noqa: E402

TABLES = ["series", "gaps", "people", "clock_atlas", "clock_drift", "clock_drift_classes"]
failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}{': ' + detail if detail and not ok else ''}")
    if not ok:
        failures.append(name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--privacy", action="store_true")
    out = ap.parse_args().out_dir
    a_privacy = ap.parse_args().privacy

    # manifest and card
    m = json.loads((out / "manifest.json").read_text())
    for f in m["files"]:
        p = out / f["file"]
        good = p.exists() and p.stat().st_size == f["bytes"] and hashlib.sha256(p.read_bytes()).hexdigest() == f["sha256"]
        if good and "rows" in f:
            good = pl.scan_parquet(p).select(pl.len()).collect().item() == f["rows"]
        check(f"manifest matches {f['file']}", good)
    card = (out / "README.md").read_text()
    named = ax.CARD_FILE.findall(card)
    check("card names every table", set(f"{t}.parquet" for t in TABLES) <= set(named), str(named))
    check("every file the card names exists", all((out / f).exists() for f in named))
    check("inputs pinned", len(m.get("inputs", [])) > 0 and all("sha256" in i or "mtime" in i for i in m["inputs"]))

    t = {name: pl.read_parquet(out / f"{name}.parquet") for name in TABLES}
    series = t["series"]
    names = set(series["series"].to_list())

    # ids, licences, publication guard
    for name, frame in t.items():
        check(f"{name}: row_id unique", frame["row_id"].n_unique() == frame.height and frame["row_id"].null_count() == 0)
        check(f"{name}: licence on every row", frame["licence_access"].null_count() == 0)
        for col in ("producing_module", "producing_commit", "export_commit"):
            check(f"{name}: {col} on every row", col in frame.columns and frame[col].null_count() == 0)
        try:
            ax.check_public(frame, name)
            check(f"{name}: public strings and ids only", True)
        except ValueError as e:
            check(f"{name}: public strings and ids only", False, str(e))
    people = t["people"]
    check(
        "people: one row per series, group, array and clock",
        people.select("series", "group", "record_id", "clock").n_unique() == people.height,
    )
    per = people.select("accession", "record_id", "donor_key").unique()
    check("people: one donor_key per array", per.height == per.select("accession", "record_id").n_unique())
    check("people: status set", set(people["status"].unique().to_list()) <= {"included", "excluded"})
    check(
        "people: a reason for every excluded row and none for included",
        people.filter((pl.col("status") == "excluded") & pl.col("status_reason").is_null()).height == 0
        and people.filter((pl.col("status") == "included") & pl.col("status_reason").is_not_null()).height == 0,
    )

    # joins to series
    for name in ("gaps", "people", "clock_atlas", "clock_drift"):
        missing = set(t[name]["series"].unique().to_list()) - names
        check(f"{name}: every row joins to a series row", not missing, str(sorted(missing))[:200])

    # ages
    unit = dict(series.select("series", "age_unit").iter_rows())
    ages = people.with_columns(pl.col("series").replace_strict(unit, default=None).alias("unit"))
    human = ages.filter(pl.col("unit") == "years")
    check("people: human ages within 0-120 years", human.filter(~pl.col("age").is_between(0, 120)).height == 0)
    days = ages.filter(pl.col("unit") == "days")
    check("people: mouse ages within 0-1500 days", days.filter(~pl.col("age").is_between(0, 1500)).height == 0)
    if people["age_precision"].unique().to_list() == ["exact"]:
        rng = series.select("series", "age_min", "age_max")
        inc = people.filter(pl.col("status") == "included").join(rng, on="series")
        bad = inc.filter((pl.col("age") < pl.col("age_min")) | (pl.col("age") > pl.col("age_max")))
        check("people: included ages inside the series' age range", bad.height == 0, str(bad.head(3).rows()))

    # counts: people against QC tables and against the other tables
    inc = people.filter(pl.col("status") == "included")
    sizes = inc.group_by("series", "group").agg(pl.col("record_id").n_unique().alias("n"), pl.col("group_arrays").first())
    check("people: included arrays per group equal group_arrays", sizes.filter(pl.col("n") != pl.col("group_arrays")).height == 0)
    cd = t["clock_drift"].select("series", "group", pl.col("arrays").alias("cd_arrays")).unique()
    both = sizes.join(cd, on=["series", "group"], how="full", coalesce=True)
    diff = both.filter(pl.col("n").is_null() | pl.col("cd_arrays").is_null() | (pl.col("n") != pl.col("cd_arrays")))
    check(
        "people groups and arrays equal clock_drift.parquet's",
        diff.height == 0,
        "; ".join(f"{r['series']} / {r['group']}: people {r['n']}, clock_drift {r['cd_arrays']}" for r in diff.to_dicts()),
    )
    paths = ax.ProjectPaths.default()
    status = people.select("series", "record_id", "status", "status_reason").unique()
    for name, spec in ax.specs().items():
        row = series.filter(pl.col("series") == name)
        if row.height != 1 or row["arrays_passing_qc"][0] is None:
            continue
        table = ax.series_table(paths, spec)
        qc = ax.reprocessed.read_qc(paths, spec["accession"]).filter(pl.col("sample").is_in(table["sample"].implode()))
        ok = qc.height == row["arrays_rebuilt"][0] and int(qc["passed"].sum()) == row["arrays_passing_qc"][0]
        ids = ax.public_ids(spec, table)
        failed = {ids[s] for s in qc.filter(~pl.col("passed"))["sample"].to_list()}
        s = status.filter(pl.col("series") == name)
        flagged_ok = s.filter(pl.col("record_id").is_in(list(failed)) & (pl.col("status") == "included")).height == 0
        check(f"{name}: QC counts match the QC table; no failed array included", ok and flagged_ok)

    if a_privacy:
        inc1 = people.filter(pl.col("status") == "included").select("series", "group", "sex", "age", "record_id").unique()
        for label, expr in (("exact", pl.col("age")), ("year", pl.col("age").floor()), ("5y", pl.col("age") // 5 * 5)):
            k = inc1.with_columns(expr.alias("a")).group_by("series", "group", "sex", "a").agg(pl.len().alias("k"))
            uniq = int(k.filter(pl.col("k") == 1)["k"].sum())
            print(f"privacy  {label}: {uniq} of {inc1.height} array-group rows unique on series, group, sex and age")

    print(f"{len(failures)} failed" if failures else "all checks passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
