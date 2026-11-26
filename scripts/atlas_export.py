"""Data export for a public explorer of the entropy-age results (atlas v1).

``uv run python scripts/atlas_export.py OUT_DIR [--part light|people|manifest|all]`` writes, into OUT_DIR:

* ``series.parquet``: one row per description file in ``config/series`` and per forest cohort outside them;
* ``gaps.parquet``: every gap of the series runner (``data/output/series_*/gaps.csv``) and the forest's other rows;
* ``people.parquet``: per analysed array, drift B, shield share F and P (``health_states.loo_scores``, chosen as
  ``clock_drift`` chooses them) with every scored clock's prediction and age-residual acceleration (the heavy part:
  run it on Slurm);
* ``clock_atlas.parquet``, ``clock_drift.parquet``, ``clock_drift_classes.parquet``: the clock summaries;
* ``README.md`` (every column) and ``manifest.json`` (commit, date, row counts, sizes, sha256).

Reads only; writes nothing outside OUT_DIR. Every string column is checked for paths, host or account names before
a table is written, and sample ids must be GEO (GSM...) or TCGA barcodes.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "output"
sys.path.insert(0, str(ROOT / "scripts"))

import forest_series  # noqa: E402  (SERIES_INFO, SKIP_GROUPS: the forest's own choice of rows)

from entropy_age.data import reprocessed  # noqa: E402
from entropy_age.experiments.series import SENSITIVITY, selection, series_table, specs  # noqa: E402
from entropy_age.paths import ProjectPaths  # noqa: E402

PLATFORM = {"HM450": "450k", "EPIC": "EPIC", "EPICv2": "EPIC v2", "Mammal40": "mammalian 40k", "MM285": "MM285"}
PUBLIC_ID = re.compile(r"^(GSM\d+|TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z])$")
#: anything that would identify a machine, an account or a file
LEAK = re.compile(
    r"(/data/|/home/|uq2347|ie2441|haumea|gpu1[123]|\bnode\d|itmz|uni-rostock|\.parquet|\.csv|\.idat|_noid|@)", re.I
)


# ----------------------------------------------------------------------------------------------------- helpers
def sex_of(v) -> str:
    s = str(v or "").strip().lower()
    if s.startswith("f"):
        return "female"
    if s.startswith("m"):
        return "male"
    return "unknown"


def species_of(spec: dict) -> str:
    return "human" if spec.get("annotation", "human") == "human" else "mouse"


def age_unit(spec: dict) -> str:
    if "age_values" in spec:
        return "group code (" + ", ".join(f"{k} = {v}" for k, v in spec["age_values"].items()) + ")"
    return "years" if species_of(spec) == "human" else "days"


def tissue_of(group: str) -> str:
    return re.split(r"[,(]", group)[0].replace("_", " ").strip()


def public_ids(spec: dict, table: pl.DataFrame) -> dict[str, str]:
    """Array name in the rebuild -> its public id (GSM for GEO; the TCGA aliquot barcode, cut to the sample
    portion as GDC lists it, for TCGA)."""
    if "samples_csv" in spec:  # sample_table keeps no barcode column, so read the GDC table itself
        csv = pl.read_csv(ROOT / spec["samples_csv"], infer_schema_length=None)
        return dict(csv.select("basename", "sample").iter_rows())
    return {s: s for s in table["sample"].to_list()}


def check_public(frame: pl.DataFrame, name: str) -> None:
    for col, dtype in frame.schema.items():
        if dtype == pl.String:
            values = frame[col].drop_nulls().unique().to_list()
        elif dtype == pl.List(pl.String):
            values = frame[col].explode().drop_nulls().unique().to_list()
        else:
            continue
        bad = [v for v in values if LEAK.search(v)]
        if bad:
            raise ValueError(f"{name}.{col}: values that are not for publication, e.g. {bad[:3]}")
    if "sample_id" in frame.columns:
        bad = [v for v in frame["sample_id"].unique().to_list() if v is None or not PUBLIC_ID.match(v)]
        if bad:
            raise ValueError(f"{name}.sample_id: ids that are not GEO or TCGA ids, e.g. {bad[:3]}")


LICENCE = {
    "GEO": "GEO public: NCBI places no restrictions on use; submitters may hold rights; cite the series",
    "TCGA": "TCGA open-access tier (GDC): open data under the NIH GDC data policy; cite TCGA",
    "other": "public deposit named by the producing module; cite the source",
}
#: forest rows outside the runner: the module whose table they come from
FOREST_MODULE = {
    "Hannum": "entropy_age.experiments.tissue_contrast",
    "Jaffe": "entropy_age.experiments.tissue_contrast",
    "Jena RRBS": "entropy_age.experiments.jena_rrbs",
    "GSE225166": "entropy_age.experiments.mouse_atlas_rrbs",
    "Thompson RRBS": "entropy_age.experiments.mouse_rrbs_tissues",
    "Meer RRBS": "entropy_age.experiments.meer_tissues",
}
BUILD_DATE = dt.datetime.now(dt.UTC).isoformat(timespec="seconds")


def licence_expr(acc: pl.Expr) -> pl.Expr:
    return (
        pl.when(acc.str.starts_with("TCGA-"))
        .then(pl.lit(LICENCE["TCGA"]))
        .when(acc.str.starts_with("GSE"))
        .then(pl.lit(LICENCE["GEO"]))
        .otherwise(pl.lit(LICENCE["other"]))
    )


def git_last(path: Path) -> tuple[str, str]:
    """(commit, ISO date) of the last commit that touched ``path``; ('unknown', '') without git."""
    import subprocess

    try:
        r = subprocess.run(
            ["git", "-C", str(ROOT), "log", "-1", "--format=%H %cI", "--", str(path.relative_to(ROOT))],
            capture_output=True, text=True, check=True,
        )
        commit, date = r.stdout.split()
        return commit, date
    except (OSError, subprocess.CalledProcessError, ValueError):
        return "unknown", ""


def derived(frame: pl.DataFrame, module: str | pl.Expr, csv: Path | None) -> pl.DataFrame:
    """Mark a table derived: the module that computes it and the commit (and its date) that last wrote the input
    table; without an input table (computed here), this build's commit and date."""
    commit, date = git_last(csv) if csv is not None else (commit_head(), BUILD_DATE)
    return frame.with_columns(
        pl.lit("derived").alias("table_role"),
        (module if isinstance(module, pl.Expr) else pl.lit(module)).alias("producing_module"),
        pl.lit(commit).alias("producing_commit"),
        pl.lit(date).alias("produced_at"),
        pl.lit(commit_head()).alias("export_commit"),
        pl.lit(BUILD_DATE).alias("export_date"),
    )


def row_id(cols: list[str]) -> pl.Expr:
    """Stable row id: a hash of the row's natural key, unchanged across builds while the key is."""
    return pl.concat_str([pl.col(c).cast(pl.String).fill_null("") for c in cols], separator="|").map_elements(
        lambda s: hashlib.sha256(s.encode()).hexdigest()[:16], return_dtype=pl.String
    )


def donor_key(acc: str, donor: str) -> str:
    """Join key of one person within an accession: a hash of the deposited donor label, so arrays of one donor (other
    tissues, other platform files of the accession) share it without the label itself being exported."""
    return hashlib.sha256(f"{acc}|{donor}".encode()).hexdigest()[:16]


def write(frame: pl.DataFrame, out: Path, name: str) -> None:
    check_public(frame, name)
    frame.write_parquet(out / f"{name}.parquet", compression="zstd")
    print(f"{name}: {frame.height} rows", file=sys.stderr, flush=True)


# ------------------------------------------------------------------------------------------------- series
def series_rows(paths: ProjectPaths) -> pl.DataFrame:
    rows = []
    for name, spec in specs().items():
        acc = spec["accession"]
        table = series_table(paths, spec)
        rebuilt = reprocessed.available(paths, acc)
        row = {
            "series": name,
            "source": "config/series",
            "accession": acc,
            "repository": "GDC (TCGA)" if acc.startswith("TCGA-") else "GEO",
            "species": species_of(spec),
            "platform": None,
            "processing": "rebuilt" if rebuilt else "no rebuild",
            "sensitivity": bool(SENSITIVITY.search(name)),
            "age_unit": age_unit(spec),
            "arrays_with_age": None if table is None else table.height,
            "arrays_rebuilt": None,
            "arrays_passing_qc": None,
            "arrays_analysed": None,
            "groups": [],
            "group_tissues": [],
            "donors_per_group": [],
            "group_age_min": [],
            "group_age_max": [],
            "age_min": None,
            "age_max": None,
            "n_female": None,
            "n_male": None,
            "n_sex_unknown": None,
            "gaps_available": (OUT / f"series_{name}" / "gaps.csv").exists(),
            "preprocessing_track": None,
            "arrays_failed_detection": None,
            "arrays_low_intensity": None,
        }
        if table is not None and rebuilt:
            qc = reprocessed.read_qc(paths, acc).filter(pl.col("sample").is_in(table["sample"].implode()))
            row["arrays_rebuilt"] = qc.height
            row["arrays_passing_qc"] = int(qc["passed"].sum())
            row["arrays_failed_detection"] = int((~qc["passed"] & ~qc["low_intensity"]).sum())
            row["arrays_low_intensity"] = int(qc["low_intensity"].sum())
            if "prep" in qc.columns:
                row["preprocessing_track"] = "SeSAMe " + ", ".join(sorted(qc["prep"].cast(pl.String).unique().to_list()))
            if "platform" in qc.columns:
                row["platform"] = ", ".join(sorted({PLATFORM.get(p, p) for p in qc["platform"].unique().to_list()}))
            chosen = selection(paths, name, spec)
            if chosen is not None:
                groups, _ = chosen
                groups = {k: v.drop_nulls("age") for k, v in groups.items()}
                every = pl.concat(list(groups.values())).unique(subset="sample") if groups else None
                sizes_path = OUT / f"series_{name}" / "sizes.csv"
                sizes = pl.read_csv(sizes_path) if sizes_path.exists() else None
                for label, s in groups.items():
                    sz = sizes.filter(pl.col("group") == label) if sizes is not None else None
                    if sz is not None and sz.height:  # the gap runner's own count (it also drops incomplete rows)
                        r = sz.row(0, named=True)
                        donors, lo, hi = r["donors"], r["age_min"], r["age_max"]
                    else:
                        donors = s.height
                        lo = float(s["age"].min()) if s.height else None
                        hi = float(s["age"].max()) if s.height else None
                    row["groups"].append(label)
                    row["group_tissues"].append(tissue_of(label))
                    row["donors_per_group"].append(donors)
                    row["group_age_min"].append(lo)
                    row["group_age_max"].append(hi)
                if every is not None and every.height:
                    row["arrays_analysed"] = every.height
                    row["age_min"], row["age_max"] = float(every["age"].min()), float(every["age"].max())
                    sexes = [sex_of(v) for v in every["sex"].to_list()]
                    row["n_female"], row["n_male"] = sexes.count("female"), sexes.count("male")
                    row["n_sex_unknown"] = sexes.count("unknown")
        rows.append(row)
    schema = {
        "groups": pl.List(pl.String),
        "group_tissues": pl.List(pl.String),
        "donors_per_group": pl.List(pl.Int64),
        "group_age_min": pl.List(pl.Float64),
        "group_age_max": pl.List(pl.Float64),
        "age_min": pl.Float64,
        "age_max": pl.Float64,
        "platform": pl.String,
    }
    frame = pl.DataFrame(rows, schema_overrides=schema, infer_schema_length=None)
    forest = forest_only()
    if forest.height:
        extra = forest.select(
            pl.col("series"),
            pl.lit("forest plot cohort").alias("source"),
            pl.col("series").str.extract(r"(GSE\d+)").alias("accession"),
            pl.when(pl.col("series").str.contains("GSE")).then(pl.lit("GEO")).otherwise(None).alias("repository"),
            "species",
            "platform",
            "processing",
            pl.lit(False).alias("sensitivity"),
            pl.when(pl.col("species") == "human").then(pl.lit("years")).otherwise(None).alias("age_unit"),
            pl.concat_list(pl.col("series")).alias("groups"),
            pl.concat_list(pl.col("series").map_elements(tissue_of, return_dtype=pl.String)).alias("group_tissues"),
            pl.concat_list(pl.col("n").cast(pl.Int64)).alias("donors_per_group"),
            pl.lit(True).alias("gaps_available"),
        )
        frame = pl.concat([frame, extra], how="diagonal_relaxed")
    frame = frame.with_columns(
        pl.col("source").replace({"config/series": "series runner"}),
        pl.col("accession").alias("source_accession"),
        licence_expr(pl.col("accession").fill_null("")).alias("licence_access"),
        pl.when(pl.col("preprocessing_track").is_null())
        .then(pl.col("processing"))
        .otherwise(pl.col("preprocessing_track"))
        .alias("preprocessing_track"),
    )
    module = (
        pl.when(pl.col("source") == "series runner")
        .then(pl.lit("entropy_age.experiments.series (selection) and entropy_age.data.reprocessed (read_qc)"))
        .otherwise(pl.lit("scripts/forest_series.py"))
    )
    return derived(frame, module, None).with_columns(
        pl.lit("metadata").alias("table_role"), pl.col("series").alias("row_id")
    )


def forest_only() -> pl.DataFrame:
    """Forest rows that come from outside ``config/series``: the headline cohorts and the sequencing modules."""
    f = pl.read_csv(ROOT / "docs" / "forest_series.csv", infer_schema_length=None)
    return f.filter(pl.col("interval") != "donors or animals")


# --------------------------------------------------------------------------------------------------- gaps
def gap_rows() -> pl.DataFrame:
    parts = []
    for name, spec in specs().items():
        path = OUT / f"series_{name}" / "gaps.csv"
        if not path.exists():
            continue
        g = pl.read_csv(path, infer_schema_length=None).rename({"cohort": "group", "group": "gene_set"})
        g = g.with_columns(
            pl.lit(name).alias("series"),
            pl.lit(spec["accession"]).alias("accession"),
            pl.lit(bool(SENSITIVITY.search(name))).alias("sensitivity"),
        )
        # the forest's row: kl_slope, housekeeping, all CpGs, held + chip where its donor interval exists, else held
        g = g.with_columns(in_forest=pl.lit(False))
        if name in forest_series.SERIES_INFO:
            key = (pl.col("metric") == "kl_slope") & (pl.col("gene_set") == "housekeeping") & (pl.col("band") == "all")
            flags = []
            for label in g["group"].unique(maintain_order=True).to_list():
                if label in forest_series.SKIP_GROUPS:
                    continue
                sub = g.filter(key & (pl.col("group") == label))
                chip = sub.filter((pl.col("covariates") == "held + chip") & pl.col("donor_ci_low").is_not_nan())
                flags.append((label, "held + chip" if chip.height else "held"))
            for label, fit in flags:
                g = g.with_columns(
                    in_forest=pl.col("in_forest") | (key & (pl.col("group") == label) & (pl.col("covariates") == fit))
                )
        g = g.with_columns(
            pl.lit("series runner").alias("source"),
            pl.lit("donors").alias("interval_kind"),
            pl.lit("SeSAMe rebuild (drift track)").alias("preprocessing_track"),
        )
        parts.append(derived(g, "entropy_age.experiments.series", path))
    gaps = pl.concat(parts, how="diagonal_relaxed")
    forest = forest_only().select(
        "series",
        pl.col("series").str.extract(r"(GSE\d+)").alias("accession"),
        pl.col("series").alias("group"),
        pl.col("fit").alias("covariates"),
        pl.lit("all").alias("band"),
        pl.lit("kl_slope").alias("metric"),
        pl.lit("housekeeping").alias("gene_set"),
        pl.col("gap").alias("gap_sd"),
        pl.col("ci_low").alias("donor_ci_low"),
        pl.col("ci_high").alias("donor_ci_high"),
        pl.col("interval").alias("interval_kind"),
        pl.lit(False).alias("sensitivity"),
        pl.lit(True).alias("in_forest"),
        pl.lit("forest (module tables)").alias("source"),
        pl.col("processing").alias("preprocessing_track"),
    )
    module = pl.lit("unknown")
    for word, mod in FOREST_MODULE.items():
        module = pl.when(pl.col("series").str.contains(word, literal=True)).then(pl.lit(mod)).otherwise(module)
    forest = derived(forest, module, ROOT / "docs" / "forest_series.csv")
    gaps = pl.concat([gaps, forest], how="diagonal_relaxed").with_columns(
        pl.col("accession").alias("source_accession"),
        licence_expr(pl.col("accession").fill_null("")).alias("licence_access"),
        row_id(["series", "group", "gene_set", "metric", "band", "covariates"]).alias("row_id"),
    )
    first = ["series", "accession", "source", "group", "gene_set", "metric", "band", "covariates", "matching"]
    return gaps.select(*first, *[c for c in gaps.columns if c not in first]).with_columns(
        pl.col("interval_kind").replace({"donors or animals": "donors"})
    )


# ------------------------------------------------------------------------------------------------- clocks
def clock_tables(out: Path) -> None:
    acc_of = {name: spec["accession"] for name, spec in specs().items()}
    for name, csv, module, key in (
        ("clock_atlas", OUT / "clock_atlas" / "clocks.csv", "entropy_age.experiments.clock_atlas", ["series", "group", "clock"]),
        ("clock_drift", OUT / "clock_drift" / "correlations.csv", "entropy_age.experiments.clock_drift", ["series", "group", "clock"]),
        ("clock_drift_classes", OUT / "clock_drift_classes" / "summary.csv", "entropy_age.experiments.clock_drift_classes", ["clock"]),
    ):
        t = pl.read_csv(csv, infer_schema_length=None)
        if "series" in t.columns:
            acc = pl.col("series").replace_strict(acc_of, default=None)
            t = t.with_columns(acc.alias("source_accession"), licence_expr(acc.fill_null("")).alias("licence_access"))
        else:  # summaries over the primary groups of many accessions
            t = t.with_columns(
                pl.lit(None, pl.String).alias("source_accession"),
                pl.lit(LICENCE["GEO"] + " | " + LICENCE["TCGA"]).alias("licence_access"),
            )
        t = t.with_columns(pl.lit("SeSAMe rebuilds (drift and clock tracks)").alias("preprocessing_track"))
        write(derived(t, module, csv).with_columns(row_id(key).alias("row_id")), out, name)


# ------------------------------------------------------------------------------------------------- people
def people(out: Path) -> None:
    """clock_drift's per-group loop, kept per array: same features, arrays, exclusions and scores."""
    from entropy_age.data.annotation import AnnotationLayout
    from entropy_age.experiments import Context, Scale, mitotic_clocks
    from entropy_age.experiments.clock_atlas import valid
    from entropy_age.experiments.clock_drift import CLOCK_SCORES, CLOCK_TRACK, GROUP, MIN_ARRAYS, SERIES, residual
    from entropy_age.experiments.clock_increment import GENERATION
    from entropy_age.experiments.health_states import loo_scores
    from entropy_age.experiments.per_person import _islands, _sets
    from entropy_age.experiments.tissue_contrast import _features

    ctx = Context(ProjectPaths.default(), Scale.of("paper"), 42)
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # the stream clock_drift uses
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    mitotic = mitotic_clocks.load(ctx.paths.root)
    if mitotic is None:
        raise FileNotFoundError("atlas_export: no mitotic coefficients; run scripts/fetch_mitotic_clocks.py")
    clock_cpgs = mitotic_clocks.cpgs(mitotic)
    frames, skipped = [], []
    for name, spec in specs().items():
        if species_of(spec) != "human":
            skipped.append((name, "mouse: clock_drift scores human series only"))
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv"
        unmasked = ctx.paths.root / CLOCK_SCORES / f"{gse}_cpgs.parquet"
        if not (betas_path.exists() and clocks_path.exists()):
            skipped.append((name, "no clock-panel betas or clock scores"))
            continue
        long = pl.read_csv(clocks_path, schema_overrides={"sample": pl.String, "prediction": pl.Float64})
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(clock_cpgs))
        available = set(betas.columns) & set(long["sample"].to_list())
        mb = None
        if unmasked.exists():
            mb = pl.scan_parquet(unmasked).filter(pl.col("cpg").is_in(clock_cpgs)).collect()
            available &= set(mb.columns)
            available &= set(
                reprocessed.read_qc(ctx.paths, gse, sesame_dir=CLOCK_TRACK).filter(pl.col("passed"))["sample"]
            )
        chosen = selection(ctx.paths, name, spec, available)
        if chosen is None:
            skipped.append((name, "no rebuild or sample table"))
            continue
        groups, dropped = chosen
        if dropped:
            betas = betas.filter(~pl.col("cpg").is_in(dropped))
        ids = public_ids(spec, series_table(ctx.paths, spec))
        for label, s in groups.items():
            s = s.drop_nulls("age")
            if s.height < MIN_ARRAYS:
                skipped.append((f"{name} / {label}", f"{s.height} arrays, under {MIN_ARRAYS}"))
                continue
            samples = s["sample"].to_list()
            scored = long.filter(pl.col("sample").is_in(samples)).select(
                "sample", "clock", "output", "prediction", "measured_weight_share", pl.lit("biolearn panel").alias("clock_source")
            )
            if mb is not None:
                extra, _ = mitotic_clocks.score(mb, samples, mitotic)
                if extra.height:
                    scored = pl.concat(
                        [scored, extra.select("sample", "clock", pl.col("prediction").cast(pl.Float64)).with_columns(
                            pl.lit("mitotic_clocks").alias("clock_source"))],
                        how="diagonal_relaxed",
                    )
            sub = betas.select("cpg", *samples)
            x = sub.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                skipped.append((f"{name} / {label}", "no housekeeping CpGs measured"))
                continue
            print(f"atlas_export: {name} / {label}: {len(samples)} arrays, {len(kept)} CpGs", file=sys.stderr, flush=True)
            per = loo_scores(x[:, ok], _islands(kept, probes), sets)
            age = s["age"].to_numpy().astype(float)
            person = pl.DataFrame(
                {
                    "sample": samples,
                    "age": age,
                    "sex": [sex_of(v) for v in s["sex"].to_list()],
                    "B": per[f"B_{GROUP}"],
                    "F": per[f"F_{GROUP}"],
                    "P": per[f"P_{GROUP}"],
                }
            )
            # acceleration per clock: residual of the prediction on age within the group, on the arrays it scored,
            # where at least MIN_ARRAYS arrays carry a score and the scores vary (clock_drift's rule)
            accs = []
            age_of = dict(zip(samples, age, strict=True))
            for (clock,), c in scored.group_by("clock", maintain_order=True):
                c = c.filter(pl.col("prediction").is_finite())
                acc = np.full(c.height, np.nan)
                v = c["prediction"].to_numpy()
                if c.height >= MIN_ARRAYS and np.std(v) > 0:
                    acc = residual(v, np.array([age_of[t] for t in c["sample"].to_list()]))
                accs.append(c.with_columns(pl.Series("acceleration", acc)))
            scored = pl.concat(accs, how="diagonal_relaxed") if accs else scored
            frames.append(
                person.join(scored, on="sample", how="left")
                .with_columns(
                    pl.col("sample").replace_strict(ids, default=None).alias("sample_id"),
                    pl.lit(name).alias("series"),
                    pl.lit(gse).alias("accession"),
                    pl.lit(label).alias("group"),
                    pl.lit(bool(SENSITIVITY.search(name))).alias("sensitivity"),
                    pl.lit(len(samples)).alias("group_arrays"),
                    pl.lit(len(kept)).alias("drift_cpgs"),
                    pl.col("clock").replace_strict(GENERATION, default="other").alias("generation"),
                    pl.col("clock").map_elements(lambda k, g=label: valid(k, g), return_dtype=pl.Boolean).alias("valid"),
                )
            )
    table = pl.concat(frames, how="diagonal_relaxed").select(
        "series", "accession", "group", "sensitivity", "group_arrays", "drift_cpgs", "sample_id", "age", "sex",
        "B", "F", "P", "clock", "clock_source", "generation", "valid", "output", "measured_weight_share",
        "prediction", "acceleration",
    )
    check_public(table, "people")
    table.write_parquet(core_path(out))  # the heavy result, kept beside the export; enrich() builds people.parquet
    for what, why in skipped:
        print(f"atlas_export: people skipped {what}: {why}", file=sys.stderr)


# ----------------------------------------------------------------------------------------------- manifest
def commit_head() -> str:
    """HEAD's commit from .git without the git binary (compute nodes have none)."""
    git = ROOT / ".git"
    head = (git / "HEAD").read_text().strip()
    if not head.startswith("ref: "):
        return head
    ref = head.removeprefix("ref: ")
    if (git / ref).exists():
        return (git / ref).read_text().strip()
    for line in (git / "packed-refs").read_text().splitlines():
        if line.endswith(" " + ref):
            return line.split()[0]
    return "unknown"


def core_path(out: Path) -> Path:
    """The heavy part's per-array scores, kept outside the published folder."""
    return out.parent / f".{out.name}_people_core.parquet"


CORE = [
    "series", "accession", "group", "sensitivity", "group_arrays", "drift_cpgs", "sample_id", "age", "sex",
    "B", "F", "P", "clock", "clock_source", "generation", "valid", "output", "measured_weight_share",
    "prediction", "acceleration",
]


def enrich(out: Path, precision: str) -> None:
    """people.parquet from the heavy scores: keys, a status row for every array a group's filters select that the
    scores leave out (with the reason), per-row provenance and licence; light (sample tables and QC files only)."""
    from entropy_age.data.cohorts import layout as cohort_layout
    from entropy_age.experiments.clock_drift import CLOCK_SCORES, CLOCK_TRACK, MIN_ARRAYS, SERIES
    from entropy_age.experiments.series import group_filter

    src = core_path(out)
    if not src.exists() and (out / "people.parquet").exists():  # a first build wrote the scores to people.parquet
        pl.read_parquet(out / "people.parquet").select(CORE).write_parquet(src)
    core = pl.read_parquet(src).select(CORE)
    built = dt.datetime.fromtimestamp(src.stat().st_mtime, dt.UTC).isoformat(timespec="seconds")
    paths = ProjectPaths.default()
    parts, extra = [], []
    for name, spec in specs().items():
        acc = spec["accession"]
        table = series_table(paths, spec)
        if table is None:
            print(f"atlas_export: enrich: {name}: no sample table", file=sys.stderr)
            continue
        ids = public_ids(spec, table)
        donor_of = {r[0]: donor_key(acc, r[1]) for r in table.select("sample", "donor").iter_rows()}
        rebuilt = reprocessed.available(paths, acc)
        qc = reprocessed.read_qc(paths, acc) if rebuilt else None
        prep = dict(qc.select("sample", pl.col("prep").cast(pl.String)).iter_rows()) if qc is not None and "prep" in qc.columns else {}
        passed = dict(qc.select("sample", "passed").iter_rows()) if qc is not None else {}
        low = dict(qc.select("sample", "low_intensity").iter_rows()) if qc is not None else {}
        clock_prep, available = {}, None
        if species_of(spec) == "human":
            cqc = None
            if (reprocessed.series_dir(paths, acc, CLOCK_TRACK) / f"{acc}.qc.csv").exists():
                cqc = reprocessed.read_qc(paths, acc, sesame_dir=CLOCK_TRACK)
                if "prep" in cqc.columns:
                    clock_prep = dict(cqc.select("sample", pl.col("prep").cast(pl.String)).iter_rows())
            clocks = paths.root / CLOCK_SCORES / f"{acc}_clocks.csv"
            betas = paths.root / SERIES / f"{acc}_cpgs.parquet"
            unmasked = paths.root / CLOCK_SCORES / f"{acc}_cpgs.parquet"
            available = set()
            if clocks.exists() and betas.exists():
                available = set(
                    pl.read_csv(clocks, columns=["sample"], schema_overrides={"sample": pl.String})["sample"]
                )
                available &= set(pl.read_parquet_schema(betas))
                if unmasked.exists():
                    available &= set(pl.read_parquet_schema(unmasked))
                    if cqc is not None:
                        available &= set(cqc.filter(pl.col("passed"))["sample"])
        chosen = selection(paths, name, spec, available) if rebuilt else None
        jaffe: set[str] = set()
        if spec.get("exclude_jaffe_donors"):
            j = pl.read_parquet(cohort_layout("jaffe", paths).samples)
            jaffe = set(j.filter((pl.col("group") == "Control") & (pl.col("age") >= 20))["donor"].to_list())
        mine = core.filter(pl.col("series") == name)
        included = set(mine.select("group", "sample_id").unique().iter_rows())
        inverse = {v: k for k, v in ids.items()}
        parts.append(
            mine.with_columns(
                pl.col("sample_id").replace_strict({k: donor_of.get(v) for k, v in inverse.items()}, default=None)
                .alias("donor_key"),
                pl.col("sample_id").replace_strict({k: prep.get(v) for k, v in inverse.items()}, default=None)
                .alias("drift_prep"),
                pl.col("sample_id").replace_strict({k: clock_prep.get(v) for k, v in inverse.items()}, default=None)
                .alias("clock_prep"),
                pl.lit("included").alias("status"),
                pl.lit(None, pl.String).alias("status_reason"),
            )
        )
        for label, filters in spec["groups"].items():
            cand = group_filter(table, filters)
            picked = set(chosen[0][label]["sample"].to_list()) if chosen is not None and label in chosen[0] else set()
            for sample, donor, age, sex in cand.select("sample", "donor", "age", "sex").iter_rows():
                sid = ids.get(sample)
                if (label, sid) in included:
                    continue
                if not rebuilt or sample not in passed:
                    why = "not in the rebuild"
                elif low.get(sample):
                    why = "failed QC: low intensity"
                elif not passed[sample]:
                    why = "failed QC: detection under 90 %"
                elif donor in jaffe:
                    why = "donor of the discovery cortex cohort (left out)"
                elif available is not None and sample not in available:
                    why = "missing from the clock track (its QC, clock scores or clock betas)"
                elif sample in picked:
                    if species_of(spec) != "human":
                        why = "mouse series: per-person scores are computed for human series only"
                    elif len(picked) < MIN_ARRAYS:
                        why = f"group under {MIN_ARRAYS} arrays"
                    else:
                        why = "group not scored (no housekeeping CpGs measured in every array)"
                else:
                    why = "repeat array of a donor (one array per donor)"
                extra.append(
                    {
                        "series": name,
                        "accession": acc,
                        "group": label,
                        "sensitivity": bool(SENSITIVITY.search(name)),
                        "sample_id": sid,
                        "age": float(age),
                        "sex": sex_of(sex),
                        "donor_key": donor_key(acc, donor),
                        "drift_prep": prep.get(sample),
                        "clock_prep": clock_prep.get(sample),
                        "status": "excluded",
                        "status_reason": why,
                    }
                )
    t = pl.concat([*parts, pl.DataFrame(extra, infer_schema_length=None)], how="diagonal_relaxed")
    if precision == "year":
        t = t.with_columns(pl.col("age").floor())
    elif precision == "5y":
        t = t.with_columns(pl.col("age") // 5 * 5)
    t = t.with_columns(
        pl.lit(precision).alias("age_precision"),
        pl.col("sample_id").alias("record_id"),
        pl.col("accession").alias("source_accession"),
        pl.col("sample_id").alias("source_record_id"),
        pl.concat_str(pl.lit("SeSAMe "), pl.col("drift_prep")).alias("drift_track"),
        pl.when(pl.col("clock").is_null())
        .then(pl.lit(None, pl.String))
        .otherwise(pl.concat_str(pl.lit("SeSAMe "), pl.col("clock_prep")))
        .alias("clock_track"),
        licence_expr(pl.col("accession")).alias("licence_access"),
        pl.lit("per array: deposited age and sex with derived scores").alias("table_role"),
        pl.lit("scripts/atlas_export.py (clock_drift's selection; health_states.loo_scores)").alias("producing_module"),
        pl.lit(commit_head()).alias("producing_commit"),
        pl.lit(built).alias("produced_at"),
        pl.lit(commit_head()).alias("export_commit"),
        pl.lit(BUILD_DATE).alias("export_date"),
        row_id(["series", "group", "sample_id", "clock"]).alias("row_id"),
    ).drop("drift_prep", "clock_prep")
    first = ["row_id", "record_id", "donor_key", "series", "accession", "group", "status", "status_reason"]
    t = t.select(*first, *[c for c in t.columns if c not in first]).sort(
        "series", "group", "status", "record_id", "clock", nulls_last=True
    )
    write(t, out, "people")


CARD_FILE = re.compile(r"^## (\S+\.(?:parquet|json))\s*$", re.M)


def check_card(out: Path) -> None:
    """The build fails when a file the card describes is missing."""
    missing = [f for f in CARD_FILE.findall(README) if f != "manifest.json" and not (out / f).exists()]
    if missing:
        raise FileNotFoundError(f"atlas_export: the card names files that are not built: {missing}")


#: inputs pinned in the manifest; a file above HASH_LIMIT is pinned by size and modification time
INPUTS = [
    "config/series/*.toml",
    "config/series/lists/*.csv",
    "docs/forest_series.csv",
    "data/output/series_*/gaps.csv",
    "data/output/series_*/sizes.csv",
    "data/output/clock_atlas/clocks.csv",
    "data/output/clock_drift/correlations.csv",
    "data/output/clock_drift_classes/summary.csv",
    "data/interim/clocks/series_unmasked/*_clocks.csv",
    "data/interim/clocks/series/*_cpgs.parquet",
    "data/interim/clocks/series_unmasked/*_cpgs.parquet",
    "data/interim/sesame/*/*.qc.csv",
    "data/interim/sesame_clock/*/*.qc.csv",
    "data/input/tcga/*/samples.csv",
    "data/input/mitotic/*.csv",
    "data/interim/gene_sets.parquet",
]
HASH_LIMIT = 300_000_000


def pinned_inputs() -> list[dict]:
    from entropy_age.data.annotation import AnnotationLayout

    paths = ProjectPaths.default()
    accs = {spec["accession"] for spec in specs().values()}
    files = {p for g in INPUTS for p in ROOT.glob(g)}
    lay = AnnotationLayout(paths)
    files |= {lay.probes_450k, lay.links_450k}
    for spec in specs().values():
        if "samples_csv" not in spec:
            acc = spec["accession"]
            files.add(paths.input / "geo" / acc / spec.get("matrix", f"{acc}_series_matrix.txt.gz"))
    pinned = []
    for p in sorted(files):
        rel = str(p.resolve().relative_to(ROOT)) if p.resolve().is_relative_to(ROOT) else p.name
        if not p.exists():
            continue
        per_acc = any(seg in rel for seg in ("clocks/", "sesame/", "sesame_clock/"))
        if per_acc and not any(a in p.name for a in accs):
            continue  # inputs of accessions no description file names
        st = p.stat()
        entry = {"path": rel, "bytes": st.st_size}
        if st.st_size <= HASH_LIMIT:
            entry["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
        else:
            entry["pinned_by"] = "size and modification time"
            entry["mtime"] = dt.datetime.fromtimestamp(st.st_mtime, dt.UTC).isoformat(timespec="seconds")
        pinned.append(entry)
    return pinned


def uncommitted() -> list[dict]:
    """Code and configuration that differ from the commit (a build from a working tree says so)."""
    import subprocess

    try:
        r = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain", "--", "scripts", "src", "config"],
            capture_output=True, text=True, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return [{"path": "unknown (no git)"}]
    changed = []
    for line in r.stdout.splitlines():
        rel = line[3:].strip()
        p = ROOT / rel
        sha = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        changed.append({"path": rel, "state": line[:2].strip(), "sha256": sha})
    return changed


def manifest(out: Path) -> None:
    files = []
    for p in sorted(out.iterdir()):
        if p.name == "manifest.json" or p.name.startswith(".") or not p.is_file():
            continue
        entry = {"file": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
        if p.suffix == ".parquet":
            entry["rows"] = pl.scan_parquet(p).select(pl.len()).collect().item()
            entry["columns"] = len(pl.read_parquet_schema(p))
        files.append(entry)
    doc = {
        "name": "entropy-age atlas export",
        "version": "atlas_v1",
        "repository_commit": commit_head(),
        "generated": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "generator": "scripts/atlas_export.py",
        "validator": "scripts/atlas_validate.py",
        "working_tree_changes": uncommitted(),
        "files": files,
        "inputs": pinned_inputs(),
    }
    (out / "manifest.json").write_text(json.dumps(doc, indent=2) + "\n")


README = r"""---
pretty_name: entropy-age atlas v1
license: other
license_name: per-row (column licence_access); a licence for the derived tables is the owner's choice and still open
tags: [dna-methylation, ageing, epigenetic-clocks, illumina-arrays]
configs:
  - {config_name: series, data_files: series.parquet}
  - {config_name: gaps, data_files: gaps.parquet}
  - {config_name: people, data_files: people.parquet}
  - {config_name: clock_atlas, data_files: clock_atlas.parquet}
  - {config_name: clock_drift, data_files: clock_drift.parquet}
  - {config_name: clock_drift_classes, data_files: clock_drift_classes.parquet}
---

# entropy-age atlas export, v1

Tables behind a public explorer of the entropy-age results: how promoter CpGs of housekeeping genes drift with age
against matched CpGs, series by series, and how each person's drift relates to the methylation clocks. Written by
`scripts/atlas_export.py` in the entropy-age repository (`uv run python scripts/atlas_export.py OUT_DIR`); the
commit it ran from is in `manifest.json`, with each file's size, row count and sha256.

Every table holds summary statistics, public sample ids (GEO `GSM` accessions, TCGA sample barcodes) and the age and
sex the depositors published beside them. The script refuses to write a string that looks like a file path, a host
or an account name, and any sample id that is not a GEO or TCGA id.

**Licence and access.** Each row carries `licence_access`: GEO series are public (NCBI places no restrictions on
their use; submitters may hold rights); TCGA arrays and their age and sex come from the open-access tier of the
GDC. Cite the source series of any row you use. The repository's code is AGPL-3.0; the licence of these derived
tables has not been chosen yet.

**Citation.** The entropy-age repository (github.com/longevity-genie/entropy-age), at the commit in `manifest.json`,
and its manuscript (in preparation); and each source series.

**Build and checks.** `uv run python scripts/atlas_export.py OUT_DIR` (parts: `light`, `people` (heavy, Slurm),
`enrich`, `manifest`); the build fails when a file this card describes is missing. `uv run python
scripts/atlas_validate.py OUT_DIR` reads the files back and checks checksums, keys, joins, ages, counts against the
QC tables and licences. `manifest.json` pins every input by sha256 (inputs above 300 MB by size and modification
time) and lists any code that differed from the commit at build time.

**Columns every table carries.**

| column | meaning |
|---|---|
| row_id | stable id of the row: a hash of its natural key (series, group, sample, clock, ...); it never moves while the key holds |
| source_accession | GEO series or TCGA project the row comes from (null for summaries over many) |
| licence_access | the access tier and terms of the source |
| preprocessing_track | how the arrays were processed (SeSAMe prep code, or deposited) |
| table_role | `metadata`, `derived` (computed by a module of the repository) or, for people, per array |
| producing_module | the module that computes the row |
| producing_commit, produced_at | the commit (and its date) that last wrote the module's output table read here; for series and people, the export's own commit and the scores' build time |
| export_commit, export_date | the commit and time of this export |

Terms used throughout:

* **series**: one description file in `config/series` (one GEO or TCGA accession, one array platform), analysed by
  the series runner `entropy_age.experiments.series`. Rebuilt series were reprocessed from the raw IDAT files with
  SeSAMe (`scripts/sesame_preprocess.R`); arrays pass QC when at least 90 % of probes are detected and their mean
  intensity is not far below the series' median (`entropy_age.data.reprocessed.read_qc`).
* **group**: a set of arrays within a series (a tissue, a disease state), one array per donor.
* **gene set**: `housekeeping`, `core_essential` or `non_essential` promoter CpGs (450k annotation, autosomes).
* **matched gap**: the set's mean age slope of a drift measure minus that of CpGs in the same stratum (CpG island
  status x decile of the young level), in units of the SD of the slope over all CpGs. Negative: the set drifts less
  than its look-alikes.
* **B, F, P** (per person): with each person left out in turn, a reference level is fitted on the others; `d` is a
  CpG's deviation from it, `|x - m| / sqrt(m (1 - m))`. `P` is the stratum-matched mean `d` of the housekeeping
  promoter CpGs, `B` that of their matched look-alikes (the person's general drift), and `F = P / B` the shield share
  (below 1: housekeeping promoters deviate less than their look-alikes).

## series.parquet

One row per description file in `config/series`, plus one row per cohort of the paper's forest plot
(`docs/forest_series.csv`, `scripts/forest_series.py`) analysed outside that runner. Produced by
`series_rows` here from `entropy_age.experiments.series` (`specs`, `series_table`, `selection`),
`entropy_age.data.reprocessed.read_qc` and the runner's `data/output/series_<name>/sizes.csv`.

| column | type | meaning |
|---|---|---|
| series | str | description file name (forest rows: the forest's label) |
| source | str | `series runner` (a description file in `config/series`) or `forest plot cohort` (a row of `docs/forest_series.csv` from outside the runner) |
| accession | str | GEO series (GSE) or TCGA project; null for forest rows whose label names none |
| repository | str | `GEO` or `GDC (TCGA)` |
| species | str | human or mouse |
| platform | str | array platform(s) of the rebuilt arrays (450k, EPIC, EPIC v2, mammalian 40k, MM285) or RRBS |
| processing | str | `rebuilt` (SeSAMe from IDATs), `deposited` (authors' betas) or `deposited counts` (RRBS) |
| sensitivity | bool | a rerun of another file with CpGs left out (`_no_field_d10/_d20`); never an independent cohort |
| age_unit | str | years; days (mouse); or the code of an age-group label |
| arrays_with_age | int | arrays in the sample table with a numeric age |
| arrays_rebuilt | int | of those, arrays present in the rebuild |
| arrays_passing_qc | int | of those, arrays passing QC (detection and intensity) |
| arrays_analysed | int | distinct arrays in any group after QC, exclusions and one array per donor |
| groups | list[str] | the file's groups |
| group_tissues | list[str] | tissue word of each group label |
| donors_per_group | list[int] | donors per group: the runner's `sizes.csv` where it ran (arrays with every CpG complete), else the selection's count; forest rows: the forest's n |
| group_age_min, group_age_max | list[float] | age range per group (unit: age_unit) |
| age_min, age_max | float | age range over the analysed arrays |
| n_female, n_male, n_sex_unknown | int | sex over the analysed arrays, as deposited |
| gaps_available | bool | the runner's gap table exists (false: the run is pending or the group had too few arrays) |
| arrays_failed_detection, arrays_low_intensity | int | arrays failing QC, by rule |

## gaps.parquet

Long format: one row per series x group x gene set x metric x CpG band x covariate set. Rows with
`source = series runner` are `data/output/series_<name>/gaps.csv` as written by
`entropy_age.experiments.series.run_series` (gaps from `brain_replication.group_gaps` and
`chromatin_match.matched_gaps`; donor intervals from `series.donor_intervals`). Rows with
`source = forest (module tables)` are the forest's other cohorts (`docs/forest_series.csv`: Hannum blood and Jaffe
cortex from `tissue_contrast`, RRBS cohorts from `jena_rrbs`, `mouse_atlas_rrbs`, `mouse_rrbs_tissues`,
`meer_tissues`), housekeeping and relative-entropy slope only.

| column | type | meaning |
|---|---|---|
| series, accession | str | as in series.parquet |
| source | str | see above |
| group | str | the series group |
| gene_set | str | housekeeping, core_essential or non_essential |
| metric | str | the per-CpG drift measure whose age slope is compared: `kl_slope` (relative entropy from the young level, bits per year), `entropy_slope` (binary entropy), `dispersion_slope` (absolute deviation) |
| band | str | `all` CpGs, or `10-90 %`: CpGs whose young level lies between 10 and 90 % |
| covariates | str | `held` (sex and the file's `hold` columns) or `held + chip` (also chip intercepts) |
| matching | str | the strata of the matched comparison (`island x level`) |
| cpgs | int | CpGs measured in every array of the group, in the band |
| gap_sd | float | matched gap, in SDs of the metric's slope over all CpGs |
| ci_low, ci_high | float | 95 % interval resampling the set's genes, arrays fixed |
| share_removed, share_ci_low, share_ci_high | float | share of the unmatched gap removed by matching, with interval (0 where one matching is run) |
| n_set | int | set CpGs in the comparison |
| donor_ci_low, donor_ci_high | float | 95 % interval resampling donors (animals for mouse), strata fixed; NaN when too few resamples kept |
| donor_resamples | int | donor resamples kept |
| sensitivity | bool | as in series.parquet |
| in_forest | bool | the row the paper's forest plot shows for this cohort |
| interval_kind | str | what the donor interval resamples: donors, animals |

## people.parquet

Long format: one row per series x group x array x clock for included arrays, and one row per array for every array
a group's filters select that the scores leave out (`status = excluded`, clock columns null, `status_reason` says
why: not in the rebuild, failed QC by detection or intensity, a repeat array of a donor, missing from the clock
track, a group under 25 arrays, a mouse series). Rows are never dropped between versions; an array that leaves the
analysis stays with its reason. Human series only, chosen exactly as
`entropy_age.experiments.clock_drift` chooses them (`series.selection`: QC on both rebuilds, one array per donor,
arrays present in the drift betas, the clock-track betas and the clock scores; groups with at least 25 arrays). An
array in two groups of one series (a `controls` subset of `every donor`) appears under each, with scores against
that group's own references. Drift CpGs: autosomal 450k promoter CpGs (gene-set members and a seeded promoter
background, `tissue_contrast._features`) measured in every array of the group, without any CpG a mitotic clock reads
and without the file's `exclude_cpgs`. The drift scores are `entropy_age.experiments.health_states.loo_scores`
(through `per_person.score_against`).

| column | type | meaning |
|---|---|---|
| series, accession, group, sensitivity | | as above |
| record_id, sample_id, source_record_id | str | the array's public id: GEO sample accession (GSM) or TCGA sample barcode |
| donor_key | str | join key of the person: a hash of the accession and the deposited donor label; shared by a donor's arrays across groups, tissues and platform files of one accession |
| status | str | included or excluded |
| status_reason | str | why an array is excluded (null when included) |
| group_arrays | int | arrays in the group |
| drift_cpgs | int | CpGs behind B, F and P in that group |
| age | float | age as deposited (years; TCGA: age at diagnosis from GDC), at `age_precision` |
| age_precision | str | exact, year (floored) or 5y (lower bound of a 5-year bin) |
| drift_track, clock_track | str | SeSAMe prep code of the array's drift rebuild and clock rebuild |
| sex | str | female, male or unknown, as deposited |
| B | float | look-alike drift: matched mean deviation of the housekeeping look-alikes from the leave-one-out reference (unitless) |
| F | float | shield share, P / B |
| P | float | matched mean deviation of housekeeping promoter CpGs (unitless) |
| clock | str | clock name |
| clock_source | str | `biolearn panel` (`data/interim/clocks/series_unmasked/<accession>_clocks.csv`, `scripts/clocks_biolearn.py` on the clock-track rebuild: SeSAMe without the design-flag mask) or `mitotic_clocks` (epiTOC2/3, HypoClock, RepliTali, epiCMIT scored per group by `entropy_age.experiments.mitotic_clocks.score`) |
| generation | str | clock generation (`clock_increment.GENERATION`); `other` for trait and other predictors |
| valid | bool | the clock was trained on this tissue or is pan-tissue (`clock_atlas.valid`); false rows are negative controls |
| output | str | the clock's output (e.g. `Age (Years)`); null for mitotic clocks |
| measured_weight_share | float | share of the clock's absolute weight measured on this array (biolearn panel) |
| prediction | float | the clock's value, in its output's unit; NaN when the array measures too little of the clock |
| acceleration | float | residual of the prediction on age (linear fit within the group over the arrays the clock scored); NaN when fewer than 25 arrays carry a score or the scores do not vary |

## clock_atlas.parquet

`data/output/clock_atlas/clocks.csv` (`entropy_age.experiments.clock_atlas`): per series group and clock.

| column | type | meaning |
|---|---|---|
| series, accession, group | str | as above |
| clock, generation, valid | | as in people.parquet |
| arrays | int | arrays scored |
| age_min, age_max | float | age range |
| r_age | float | Pearson correlation of the prediction with age |
| median_abs_error | float | median absolute prediction error in years (meaningful where the output is years) |

## clock_drift.parquet

`data/output/clock_drift/correlations.csv` (`entropy_age.experiments.clock_drift`): per group and clock, the
Pearson correlation of the clock's acceleration with B and with F over the people of people.parquet.

| column | type | meaning |
|---|---|---|
| series, group, sensitivity | | as above |
| arrays | int | arrays in the group |
| scored | int | arrays with a score for this clock |
| cpgs | int | drift CpGs |
| clock, generation, valid | | as above |
| r_B, r_F | float | correlation of the acceleration with B and with F |

## clock_drift_classes.parquet

`data/output/clock_drift_classes/summary.csv` (`entropy_age.experiments.clock_drift_classes`): per clock, over one
primary group per accession (non-brain), whose CpG class's drift tracks the clock. Classes: `hk` housekeeping
promoters, `op` matched other promoters, `np` matched non-promoter CpGs, each matched on island status x mean
level x SD x probe chemistry and scored per person against leave-one-out references, every clock input CpG left out.
`r_<class>` is the partial correlation of the clock with the class score given age, quality terms and the draw's own
background score, averaged over draws.

| column | meaning |
|---|---|
| clock | clock name |
| mitotic | a mitotic clock |
| accessions | accessions (primary groups) summarised |
| interval_scope | what the bootstrap intervals hold fixed |
| r_<c>_median | median over groups of class c's r |
| r_<c>_pooled_bg_median | the same against one pooled background (a check) |
| <contrast>_median | median over groups of the contrast; contrasts `op_minus_np` (promoters beyond non-promoters) and `hk_minus_op` (housekeeping beyond other promoters) |
| <contrast>_re, _re_lo, _re_hi | DerSimonian-Laird random-effects mean and normal 95 % interval (within-group SEs from the donor bootstrap) |
| <contrast>_re_hk_lo, _re_hk_hi | Hartung-Knapp 95 % interval of that mean |
| <contrast>_tau2 | between-group variance |
| <contrast>_median_lo, _median_hi | percentile interval of the median, resampling accessions |
| <contrast>_below0, _above0 | groups whose interval lies below or above zero |
| <contrast>_pooled_bg_median | median contrast against the pooled background |
| <contrast>_full_groups | groups covered by the full procedure (donor resamples rebuild strata, draws and references) |
| <contrast>_re_full, _re_full_lo, _re_full_hi, _re_full_hk_lo, _re_full_hk_hi | the random-effects mean and both intervals with the full procedure's SEs (epiTOC2 and HypoClock only) |
| <contrast>_full_se_ratio | median ratio of the full procedure's SE to the conditional bootstrap SE |

## manifest.json

Repository commit (read from `.git`), generation time (UTC), and per file: bytes, sha256, and for parquet files rows
and columns.
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--part", choices=["light", "people", "enrich", "manifest", "all"], default="all")
    ap.add_argument(
        "--age-precision",
        choices=["exact", "year", "5y"],
        default="exact",
        help="people.age as deposited, floored to whole years, or to the lower bound of a 5-year bin",
    )
    a = ap.parse_args()
    out = a.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    paths = ProjectPaths.default()
    if a.part in ("light", "all"):
        write(series_rows(paths), out, "series")
        write(gap_rows(), out, "gaps")
        clock_tables(out)
        (out / "README.md").write_text(README)
    if a.part in ("people", "all"):
        people(out)  # heavy: Slurm
    if a.part in ("enrich", "all"):
        enrich(out, a.age_precision)
    if a.part in ("enrich", "manifest", "all"):
        check_card(out)
        manifest(out)


if __name__ == "__main__":
    main()
