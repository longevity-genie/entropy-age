"""Matching variables from another cohort's young donors (planned test P6, the young-SD question).

``level_match`` found that matching on the young between-donor SD as well as the young level halves the blood
housekeeping gap (Hannum). Low young SD may itself be part of the protection: the same donors' young SD then carries
the set's slower drift, and matching on it removes some of the gap it should measure. Here the matching variables
(young mean level, young SD, deciles) come from the young donors of a different 450k blood cohort, so no donor and
no measurement enters both the matching variables and the drift:

* **Drift in Hannum** (GSE40279, whole blood, 656 people aged 19-101; the paper's blood estimate), against
  references from EPIC-Italy's youngest fifth (GSE51032, buffy coat, 34-47 y), from the CHARGE/Kabuki series'
  controls (GSE97362, whole blood, 124 children and adolescents aged 2-19) and from the young adults of GSE73103
  (Voisin 2015, whole blood of healthy young Swedes: the 216 aged 18-30, all weight categories), and from E-Risk
  (GSE105018, whole blood of British twins all aged 18, deposited betas made with wateRmelon's pfilter and dasen:
  another processing than GSE73103's NOOB + BMIQ). E-Risk gives one twin per family, drawn at random once (seeded),
  so every reference row is an independent family; resamples then draw those twins.
* **Drift in EPIC-Italy** (all arrays with intact signal, 34-72 y, covariates sex, later cancer diagnosis and chip),
  against Hannum's youngest fifth (19-52 y), the GSE97362 controls, the GSE73103 young adults and E-Risk.
* ``same donors``: the drift cohort's own youngest fifth, as in the paper, on the same CpGs.

By default every estimate is 450k, deposited betas (Hannum has no raw arrays), on the CpGs measured in all five cohorts
(``tissue_contrast``'s promoter features; at most 5 % missing per cohort, the rest median-imputed). The gaps are
``level_match.gaps`` (pool ``rest``): island x decile (the paper's matching), the finer strata, and the nearest five
on young level alone and on young level and young SD, in SD of the drift measure over all CpGs of the drift cohort
(one divisor per cohort and measure). Below zero means the set drifts less than its matched CpGs. Intervals resample
donors of the drift cohort and of the reference cohort independently (``ENTROPY_AGE_BOOT`` resamples,
``ENTROPY_AGE_WORKERS`` spawned processes); the young fifth and every stratum and neighbour are recomputed in each.

Equal-size references (``subsample``): every reference, the same-donor youngest fifth included, cut to
``SUBSAMPLE_SIZE`` donors drawn with replacement from its reference rows (an m-out-of-n bootstrap, so the share of a
cohort a draw uses does not set its spread), the drift cohort kept whole; as many draws as resamples. The spread over
draws is the sampling variation of each reference's matching variables, and of the gaps matched on them, at one
reference size; ``sd_diff`` is the set-minus-matched young SD, a summary of those SDs, and its spread over draws is
that summary's precision at 100 donors.

**Rebuilt betas** (``ENTROPY_AGE_BETAS=sesame``; output ``level_match_reference_rsesame``): every cohort read from its
raw IDATs through the one SeSAMe pipeline (``scripts/sesame_preprocess.R``), arrays failing its QC left out. Hannum has
no raw arrays, so the drift cohort is EPIC-Italy alone, against its own youngest fifth, the GSE97362 controls, the
GSE73103 young adults and E-Risk; one processing per estimate, and the run refuses to start unless all four series are
rebuilt. The same run writes ``processing``: the housekeeping-minus-matched young SD (island x decile) of the two
young-adult cohorts read both ways, deposited and rebuilt, on the same donors (arrays passing the rebuild's QC) and on
the CpGs every one of the four versions measures, with donor resamples shared by the two versions of a cohort (paired
rebuilt-minus-deposited) and drawn independently between cohorts (GSE73103 minus E-Risk under each processing).

``entropy-age figures level_match_reference`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import gzip
import io
import itertools
import multiprocessing
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.incident_cancer import chips
from entropy_age.experiments.incident_cancer import read_people as epic_italy_people
from entropy_age.experiments.level_match import (
    KEYS,
    MIN_VALID_SHARE,
    _finite_quantile,
    gaps,
    strata,
    stratum_gap,
    young_rows,
)
from entropy_age.experiments.syndromes import CONTROL, matrix_path, read_betas, read_deposited, usable
from entropy_age.experiments.syndromes import read_people as syndrome_people
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    MAX_MISSING,
    SETS,
    Prepared,
    _features,
    _hannum,
    _lifespan_metrics,
    _settings,
)

EPIC_ITALY, CHARGE, YOUNG_ADULTS, ERISK = "GSE51032", "GSE97362", "GSE73103", "GSE105018"
#: E-Risk's betas: a supplementary CSV (CpGs by the study's sample numbers), cached once as parquet parts
ERISK_CSV = "GSE105018_NormalisedData.csv.gz"
ERISK_AGE = 18.0  # every E-Risk participant was sampled at age 18
ERISK_ROWS_PER_PART = 20_000
SUBSAMPLE_SIZE = 100
#: the age band of the young-adult reference (GSE73103 also holds 130 adolescents aged 14-16 and ten people aged 31-34)
ADULT_AGES = (18.0, 30.0)
POOLS = ("rest",)
SAME = "same donors"
#: drift cohort -> references; a reference is a cohort key and whether its youngest fifth (True) or all its
#: donors (False) give the matching variables
PLAN = {
    "hannum": {
        SAME: ("hannum", True),
        "EPIC-Italy youngest fifth": ("epic_italy", True),
        "GSE97362 controls": ("charge", False),
        "GSE73103 young adults": ("young_adults", False),
        "E-Risk one twin per family": ("erisk", False),
    },
    "epic_italy": {
        SAME: ("epic_italy", True),
        "Hannum youngest fifth": ("hannum", True),
        "GSE97362 controls": ("charge", False),
        "GSE73103 young adults": ("young_adults", False),
        "E-Risk one twin per family": ("erisk", False),
    },
}
#: rebuilt betas: Hannum has no raw arrays, so EPIC-Italy is the only drift cohort
SESAME_PLAN = {"epic_italy": {k: v for k, v in PLAN["epic_italy"].items() if v[0] != "hannum"}}
#: the series every cohort of a rebuilt run is read from
REBUILT = {"epic_italy": EPIC_ITALY, "charge": CHARGE, "young_adults": YOUNG_ADULTS, "erisk": ERISK}
CSV_ROWS = 5_000
#: the drift cohort's own young mean and SD, carried as checks: after matching on another cohort's young donors,
#: how far do the set and its matched CpGs still differ in the drift cohort's own young donors?
OWN = ("own_level", "own_sd")


# ----------------------------------------------------------------------------- cohorts
def processing(paths=None) -> str:
    """``deposited`` or ``sesame`` (``ENTROPY_AGE_BETAS``): one processing per estimate. A rebuilt run needs every one
    of its cohorts rebuilt (Hannum, which has no raw arrays, is left out of it); any other value is refused."""
    mode = os.environ.get("ENTROPY_AGE_BETAS", "deposited")
    if mode not in ("deposited", "sesame"):
        raise RuntimeError(f"level_match_reference reads deposited or sesame betas, not {mode!r}")
    if mode == "sesame" and paths is not None:
        missing = [g for g in REBUILT.values() if not reprocessed.available(paths, g)]
        if missing:
            raise RuntimeError(f"level_match_reference: no SeSAMe rebuild of {missing}; a rebuilt run needs all four")
        for gse in REBUILT.values():
            recipe(paths, gse)
    return mode


RECIPE, PLATFORM = "QCDPB", "HM450"


def recipe(paths, gse: str) -> str:
    """The SeSAMe recipe of ``gse``'s rebuild, checked to be the one recipe of a rebuilt run: the default rebuild
    folder, HM450 arrays only, and ``prep`` QCDPB for every array. A QC table without ``prep`` was written before the
    script chose a recipe per platform (and before ``SESAME_PREP`` existed), every array prepared with QCDPB."""
    if reprocessed.series_dir(paths, gse).parent.name != "sesame":
        raise RuntimeError(
            "level_match_reference: a rebuilt run reads the default rebuild (ENTROPY_AGE_SESAME_DIR unset)"
        )
    qc = reprocessed.read_qc(paths, gse)
    platforms = set(qc["platform"].cast(pl.String).to_list()) if "platform" in qc.columns else set()
    if platforms != {PLATFORM}:
        raise RuntimeError(f"level_match_reference: {gse} rebuild platforms {sorted(platforms)}, not {PLATFORM} alone")
    if "prep" not in qc.columns:
        return f"{RECIPE} (fixed recipe of the first script)"
    preps = set(qc["prep"].cast(pl.String).to_list())
    if preps != {RECIPE}:
        raise RuntimeError(f"level_match_reference: {gse} rebuilt with {sorted(preps)}, not {RECIPE} alone")
    return RECIPE


def active_plan(mode: str) -> dict:
    return SESAME_PLAN if mode == "sesame" else PLAN


def _series(paths, gse: str, people: pl.DataFrame, features: list[str]):
    """``people`` with deposited betas (arrays whose raw signal failed left out), their matrix over the CpGs with at
    most ``MAX_MISSING`` missing (median-imputed), the same matrix with its gaps kept (NaN), and those CpGs."""
    matrix = matrix_path(paths, gse)
    keep = usable(matrix, people["sample"].to_list())
    people = people.filter(pl.col("sample").is_in(keep))
    x = read_betas(matrix, features, people["sample"].to_list())
    return people, *impute(x, features)


def impute(x: np.ndarray, features: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Columns missing in more than ``MAX_MISSING`` of the rows dropped; the rest median-imputed (FeatureMatrix's
    rule, for the drift fit) and as measured, gaps kept (for the matching variables)."""
    missing = np.isnan(x).mean(axis=0)
    keep = missing <= MAX_MISSING
    raw = x[:, keep]
    x = np.where(np.isnan(raw), np.nanmedian(raw, axis=0, keepdims=True), raw) if np.isnan(raw).any() else raw
    return x.astype(np.float32), raw.astype(np.float32), [f for f, k in zip(features, keep, strict=True) if k]


def _epic_italy(paths, features) -> Prepared:
    matrix = matrix_path(paths, EPIC_ITALY)
    people = epic_italy_people(matrix).drop_nulls("age")
    chip = chips(matrix)
    people = people.with_columns(chip=pl.col("sample").replace_strict(chip, default=None))
    people, x, raw, kept = _series(paths, EPIC_ITALY, people.sort("age", "sample"), features)
    cov = np.column_stack(
        [
            people["female"].cast(pl.Float64).to_numpy(),
            people["icd"].is_not_null().cast(pl.Float64).to_numpy(),
            dummies(people["chip"].to_numpy()),
        ]
    )
    cov = cov[:, np.ptp(cov, axis=0) > 0]
    donors = people["sample"].to_numpy()
    return Prepared("blood (EPIC-Italy)", x, people["age"].to_numpy(), cov, kept, donors, cov), raw


def _charge_controls(paths, features) -> Prepared:
    people = syndrome_people(matrix_path(paths, CHARGE))
    people = people.filter((pl.col("disease") == CONTROL) & pl.col("age").is_not_null()).sort("age", "sample")
    people, x, raw, kept = _series(paths, CHARGE, people, features)
    sex = dummies(people["sex"].to_numpy())
    label = "blood (GSE97362 controls)"
    return Prepared(label, x, people["age"].to_numpy(), sex, kept, people["sample"].to_numpy(), sex), raw


def young_adult_people(matrix: Path) -> pl.DataFrame:
    """GSE73103's people aged within ``ADULT_AGES`` (inclusive): accession, sex, age, weight category."""
    people = _characteristics(matrix).select(
        "sample",
        pl.col("Sex").alias("sex"),
        pl.col("age").cast(pl.Float64, strict=False),
        pl.col("weight category").alias("weight"),
    )
    lo, hi = ADULT_AGES
    return people.filter(pl.col("age").is_between(lo, hi)).sort("age", "sample")


def _young_adults(paths, features) -> Prepared:
    people = young_adult_people(matrix_path(paths, YOUNG_ADULTS))
    people, x, raw, kept = _series(paths, YOUNG_ADULTS, people, features)
    sex = dummies(people["sex"].to_numpy())
    label = "blood (GSE73103 young adults)"
    return Prepared(label, x, people["age"].to_numpy(), sex, kept, people["sample"].to_numpy(), sex), raw


def erisk_people(matrix: Path) -> pl.DataFrame:
    """E-Risk's arrays: accession, the study's sample number (the betas CSV's column names, read from the titles
    ``whole blood taken from sample N``), family, twin, sex and zygosity."""
    titles: list[str] = []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
                break
    prefix = "whole blood taken from sample "
    if not titles or not all(t.startswith(prefix) for t in titles):
        raise ValueError(f"{matrix.name}: sample titles do not name the study's sample numbers")
    chars = _characteristics(matrix)
    if len(titles) != chars.height:
        raise ValueError(f"{matrix.name}: {len(titles)} titles for {chars.height} samples")
    return chars.select(
        "sample",
        pl.Series("column", [t.removeprefix(prefix) for t in titles]),
        pl.col("familyid").alias("family"),
        pl.col("twinid").alias("twin"),
        pl.col("gender").alias("sex"),
        pl.col("zygosity"),
    )


def one_per_family(people: pl.DataFrame, rng: np.random.Generator, usable: list[str] | None = None) -> pl.DataFrame:
    """One row per ``family``, drawn uniformly at random among its rows (the order of ``people`` does not matter).

    With ``usable`` (the arrays with betas in the processing read) the draw is made over every row first and then
    restricted to the usable ones, so a family keeps the twin drawn without the restriction whenever that array is
    usable, takes its co-twin when it is not, and drops out only when none of its arrays is usable."""
    people = people.sort("family", "sample")
    draw = pl.Series("draw", rng.random(people.height))
    people = people.with_columns(draw)
    if usable is not None:
        people = people.filter(pl.col("sample").is_in(usable))
    return people.filter(pl.col("draw") == pl.col("draw").min().over("family")).drop("draw").sort("sample")


def erisk_parts(paths) -> Path:
    """E-Risk's betas as parquet parts (``cpg`` plus one Float32 column per sample number), written once by streaming
    the gzip CSV in blocks of ``ERISK_ROWS_PER_PART`` CpGs. Each writer fills its own temporary folder, renamed into
    place only when complete and non-empty; a writer that loses the race to another drops its copy."""
    source = paths.input / "geo" / ERISK / ERISK_CSV
    done = paths.interim / "cohorts" / "erisk_betas"
    if done.exists():
        return done
    done.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(tempfile.mkdtemp(prefix="erisk_betas.partial.", dir=done.parent))
    with gzip.open(source, "rb") as handle:
        header = handle.readline()
        names = [c.strip().strip('"') for c in header.decode().rstrip("\r\n").split(",")]
        schema = {names[0]: pl.String} | {n: pl.Float32 for n in names[1:]}
        part = 0
        while True:
            block = list(itertools.islice(handle, ERISK_ROWS_PER_PART))
            if not block:
                break
            frame = pl.read_csv(
                io.BytesIO(header + b"".join(block)), schema=schema, null_values=["", "NA", "NaN"]
            ).rename({names[0]: "cpg"})
            frame.write_parquet(partial / f"part-{part:04d}.parquet")
            part += 1
    if part == 0:
        shutil.rmtree(partial)
        raise ValueError(f"{source.name}: no CpG rows")
    try:
        partial.rename(done)
    except OSError:  # another writer finished first
        shutil.rmtree(partial)
        if not done.exists():
            raise
    return done


def read_erisk_betas(parts: Path, features: list[str], columns: list[str]) -> np.ndarray:
    """``columns`` (sample numbers) by ``features`` (NaN where E-Risk lacks a CpG), scanning only those CpGs."""
    frame = pl.scan_parquet(parts / "*.parquet").filter(pl.col("cpg").is_in(features)).select("cpg", *columns).collect()
    if frame["cpg"].n_unique() != frame.height:
        raise ValueError("E-Risk betas: a CpG appears twice")
    frame = pl.DataFrame({"cpg": features}).join(frame, on="cpg", how="left", maintain_order="left")
    if frame["cpg"].to_list() != list(features):
        raise ValueError("E-Risk betas: CpG order lost in the join")
    return frame.select(columns).to_numpy().T.astype(np.float32)


def _erisk(paths, features, rng: np.random.Generator) -> Prepared:
    """E-Risk, one twin per family among the arrays usable in the processing read (arrays whose raw signal failed are
    left out of the deposited betas too, as ``usable`` does for every series): deposited betas from the supplementary
    CSV, rebuilt betas by accession."""
    matrix = matrix_path(paths, ERISK)
    people = erisk_people(matrix)
    people = one_per_family(people, rng, usable(matrix, people["sample"].to_list()))
    if reprocessed.in_use():
        x = read_betas(matrix, features, people["sample"].to_list())
    else:
        x = read_erisk_betas(erisk_parts(paths), features, people["column"].to_list())
    x, raw, kept = impute(x, features)
    sex = dummies(people["sex"].to_numpy())
    ages = np.full(people.height, ERISK_AGE)
    label = "blood (E-Risk, one twin per family)"
    return Prepared(label, x, ages, sex, kept, people["sample"].to_numpy(), sex), raw


def restrict(p: Prepared, features: list[str], raw: np.ndarray | None = None):
    """``p`` on ``features`` (all of which it holds), in that order; with ``raw`` (``p``'s matrix with gaps kept)
    also that matrix on the same columns."""
    index = {f: j for j, f in enumerate(p.features)}
    cols = np.array([index[f] for f in features])
    out = dataclasses.replace(p, x=np.ascontiguousarray(p.x[:, cols]), features=list(features))
    return out if raw is None else (out, np.ascontiguousarray(raw[:, cols]))


MIN_OBSERVED = 10  # a reference CpG measured in fewer of the reference rows has no matching variables (NaN)


def reference_stats(raw: np.ndarray, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and SD (ddof 1) per CpG over the measured values of ``rows`` (gaps skipped, never filled from other
    donors); NaN where fewer than ``MIN_OBSERVED`` are measured."""
    x = raw[rows].astype(float)
    ok = np.isfinite(x)
    count = ok.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(ok, x, 0.0).sum(axis=0) / count
        sd = np.sqrt(np.where(ok, (x - mean) ** 2, 0.0).sum(axis=0) / (count - 1))
    thin = count < MIN_OBSERVED
    return np.where(thin, np.nan, mean), np.where(thin, np.nan, sd)


# ----------------------------------------------------------------------------- one sample
def reference_rows(p: Prepared, rows: np.ndarray, youngest: bool) -> np.ndarray:
    """The rows of ``rows`` that give the matching variables: the youngest fifth, or all."""
    return rows[young_rows(p.ages[rows])] if youngest else rows


def tables(
    drift: Prepared,
    rows: np.ndarray,
    references: dict[str, tuple[Prepared, np.ndarray, np.ndarray]],
    probes: pl.DataFrame,
    drift_raw: np.ndarray | None = None,
) -> dict[str, pl.DataFrame]:
    """Per reference ``(cohort, its matrix with gaps kept, rows)``: drift measures of ``drift`` over ``rows`` with
    young mean (``baseline``) and young SD over the reference's rows, island status and probe design; also the drift
    cohort's own young mean and SD (``own_level``, ``own_sd``; measured values of ``drift_raw``, default ``drift.x``).
    All cohorts share the CpG order. Only CpGs with finite values in every reference's table are kept, so every
    reference is read on the same CpGs."""
    m = _lifespan_metrics(drift, rows)
    own_level, own_sd = reference_stats(drift.x if drift_raw is None else drift_raw, rows[young_rows(drift.ages[rows])])
    m = m.with_columns(own_level=pl.Series(own_level), own_sd=pl.Series(own_sd))
    annotation = probes.select(pl.col("cpg").alias("feature"), "cpg_island", "probe_type")
    out = {}
    for name, (ref, raw, ref_rows) in references.items():
        if ref.features != drift.features or m["feature"].to_list() != drift.features:
            raise ValueError("level_match_reference: cohorts must share the CpG order")
        mean, sd = reference_stats(raw, ref_rows)
        t = m.with_columns(baseline=pl.Series(mean), young_sd=pl.Series(sd))
        t = t.join(annotation, on="feature", how="left")
        if t["probe_type"].null_count():
            raise ValueError("level_match_reference: a CpG has no probe design")
        out[name] = t.with_columns(pl.col("cpg_island").fill_null(False))
    columns = (*LIFESPAN_METRICS, "baseline", "young_sd", *OWN)
    shared = np.logical_and.reduce([np.all(np.isfinite(t.select(columns).to_numpy()), axis=1) for t in out.values()])
    return {name: t.filter(pl.Series(shared)) for name, t in out.items()}


def sample_rows(
    key: str,
    cohorts: dict,
    probes,
    membership,
    draws: dict[str, np.ndarray],
    size: int | None = None,
    rng: np.random.Generator | None = None,
    plan: dict | None = None,
):
    """Gaps of every reference for drift cohort ``key``; ``cohorts`` maps a key to (Prepared, matrix with gaps kept)
    and ``draws`` gives each cohort's rows (a resample or all). With ``size`` every reference's rows are cut to that
    many, drawn by ``rng`` with replacement. Also the SD of each drift measure over the drift cohort's CpGs (the
    divisor, used from the full sample only). ``plan`` maps drift cohorts to references (default ``PLAN``)."""
    drift, drift_raw = cohorts[key]
    references = {}
    for name, (ref_key, youngest) in (plan or PLAN)[key].items():
        ref, raw = cohorts[ref_key]
        rows = reference_rows(ref, draws[ref_key], youngest)
        if size is not None:
            rows = np.sort(rng.choice(rows, size=size, replace=True))
        references[name] = (ref, raw, rows)
    out, spread = [], []
    for name, t in tables(drift, draws[key], references, probes, drift_raw).items():
        long = _long(gaps(t, membership, POOLS, OWN))
        out.append(long.with_columns(reference=pl.lit(name)))
        if name == SAME:
            for metric in LIFESPAN_METRICS:
                v = t[metric].filter(t[metric].is_finite())
                spread.append({"tissue": drift.label, "metric": metric, "sd": float(v.std()) if v.len() > 1 else None})
    return pl.concat(out).with_columns(tissue=pl.lit(drift.label)), pl.DataFrame(spread)


def _long(rows: list[dict]) -> pl.DataFrame:
    wide = pl.DataFrame(rows)
    keys = ["pool", "scope", "matching", "group", "n", "n_matched", "decile_top", "baseline", "young_sd"]
    keys += ["abs_level_diff", *OWN]
    return wide.unpivot(
        index=keys, on=list(LIFESPAN_METRICS), variable_name="metric", value_name="matched_diff"
    ).rename({"baseline": "level_diff", "young_sd": "sd_diff", "own_level": "own_level_diff", "own_sd": "own_sd_diff"})


# ----------------------------------------------------------------------------- bootstrap
_STATE: dict = {}


def _init(state: dict) -> None:
    _STATE.clear()
    _STATE.update(state)


def _boot(job: tuple[str, int]) -> list[dict]:
    """One resample: donors of every cohort drawn independently (each cohort holds one array per person)."""
    key, seed = job
    cohorts = _STATE["cohorts"]
    rng = np.random.default_rng(seed)
    draws = {k: np.sort(rng.choice(p.x.shape[0], size=p.x.shape[0], replace=True)) for k, (p, _) in cohorts.items()}
    rows, _ = sample_rows(key, cohorts, _STATE["probes"], _STATE["membership"], draws, plan=_STATE.get("plan"))
    return rows.with_columns(boot=pl.lit(seed)).to_dicts()


def _subsample(job: tuple[str, int]) -> list[dict]:
    """One equal-size draw: every cohort whole, every reference cut to ``SUBSAMPLE_SIZE`` of its reference rows."""
    key, seed = job
    cohorts = _STATE["cohorts"]
    full = {k: np.arange(p.x.shape[0]) for k, (p, _) in cohorts.items()}
    rng = np.random.default_rng(seed)
    plan = _STATE.get("plan")
    rows, _ = sample_rows(key, cohorts, _STATE["probes"], _STATE["membership"], full, SUBSAMPLE_SIZE, rng, plan)
    keep = (pl.col("scope") == "all CpGs") & pl.col("matching").is_in(SUBSAMPLE_MATCHINGS)
    return rows.filter(keep).with_columns(draw=pl.lit(seed)).to_dicts()


SUBSAMPLE_MATCHINGS = ("island x decile", "nearest 5 (level)", "nearest 5 (level, SD)")


def _parallel(jobs, workers: int, state: dict, task=_boot) -> list[dict]:
    """Spawned workers (a forked polars pool deadlocks: ``tissue_contrast._run_parallel``)."""
    if workers <= 1:
        _init(state)
        return [r for j in jobs for r in task(j)]
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=_init, initargs=(state,)) as pool:
        return [r for rows in pool.map(task, jobs) for r in rows]


def _with_delta(frame: pl.DataFrame) -> pl.DataFrame:
    """``delta_sd``: the gap minus the same-donor reference's gap under the same matching in the same sample. Only
    for the scope ``all CpGs``: the other scopes are cut by each reference's own young level, so their CpGs differ
    between references."""
    on = [k for k in KEYS if k != "reference"] + (["boot"] if "boot" in frame.columns else [])
    base = frame.filter(pl.col("reference") == SAME).select(*on, pl.col("gap_sd").alias("same_sd"))
    delta = pl.when(pl.col("scope") == "all CpGs").then(pl.col("gap_sd") - pl.col("same_sd"))
    return frame.join(base, on=on, how="left").with_columns(delta_sd=delta.otherwise(None))


def summarise(point: pl.DataFrame, boot: pl.DataFrame, spread: pl.DataFrame) -> pl.DataFrame:
    """Point gaps in SD with 95 % donor-resample intervals of the gap and of its paired change against the
    same-donor reference; intervals only when at least half the resamples are finite."""
    scaled = pl.col("matched_diff") / pl.col("sd")
    on = ["tissue", "metric"]
    point = _with_delta(point.join(spread, on=on, how="left").with_columns(gap_sd=scaled))
    if not boot.height:
        return point.sort(KEYS)
    b = _with_delta(boot.join(spread, on=on, how="left").with_columns(gap_sd=scaled))
    ci = b.group_by(KEYS).agg(
        n_boot=pl.len(),
        n_valid=pl.col("gap_sd").is_finite().sum(),
        n_valid_delta=pl.col("delta_sd").is_finite().sum(),
        gap_lo=_finite_quantile("gap_sd", 0.025),
        gap_hi=_finite_quantile("gap_sd", 0.975),
        delta_lo=_finite_quantile("delta_sd", 0.025),
        delta_hi=_finite_quantile("delta_sd", 0.975),
    )
    gates = {"gap": "n_valid", "delta": "n_valid_delta"}  # each interval gated on its own finite resamples
    ci = ci.with_columns(
        [
            pl.when(pl.col(gates[c.split("_")[0]]) < MIN_VALID_SHARE * pl.col("n_boot"))
            .then(None)
            .otherwise(c)
            .alias(c)
            for c in ("gap_lo", "gap_hi", "delta_lo", "delta_hi")
        ]
    )
    return point.join(ci, on=KEYS, how="left").sort(KEYS)


#: how each cohort's betas were made: the depositors' processing, or the one SeSAMe pipeline from raw IDATs
SESAME = "rebuilt from raw IDATs (SeSAMe prepSesame QCDPB; arrays failing detection or intensity QC left out)"
PROCESSING = {
    "hannum": "deposited (depositors' processing)",
    "epic_italy": "deposited (depositors' processing)",
    "charge": "deposited (depositors' processing)",
    "young_adults": "deposited (NOOB + BMIQ; SNP, cross-reactive and sex-chromosome probes removed)",
    "erisk": "deposited (methylumi import; wateRmelon pfilter and dasen); one twin per family, seeded",
}


def summarise_subsamples(draws: pl.DataFrame, spread: pl.DataFrame) -> pl.DataFrame:
    """Per reference at ``SUBSAMPLE_SIZE`` donors: mean and 2.5-97.5 % range over the finite draws of the gap (SD
    units) and of the set-minus-matched young SD (``sd_diff``, beta units), the SD of each over draws (its precision)
    and the number of finite draws behind each; null where fewer than ``MIN_VALID_SHARE`` of the draws are finite."""
    on = ["tissue", "metric"]
    d = draws.join(spread, on=on, how="left").with_columns(gap_sd=pl.col("matched_diff") / pl.col("sd"))
    aggs = {"n_draws": pl.len(), "reference_size": pl.lit(SUBSAMPLE_SIZE)}
    for name, col in (("gap", "gap_sd"), ("sd_diff", "sd_diff")):
        finite = pl.col(col).filter(pl.col(col).is_finite())
        aggs |= {
            f"{name}_valid": pl.col(col).is_finite().sum(),
            f"{name}_mean": finite.mean(),
            f"{name}_sd_over_draws": finite.std(),
            f"{name}_lo": _finite_quantile(col, 0.025),
            f"{name}_hi": _finite_quantile(col, 0.975),
        }
    out = d.group_by(KEYS).agg(**aggs)
    gated = []
    for name in ("gap", "sd_diff"):
        thin = pl.col(f"{name}_valid") < MIN_VALID_SHARE * pl.col("n_draws")
        gated += [
            pl.when(thin).then(None).otherwise(pl.col(c)).alias(c)
            for c in (f"{name}_mean", f"{name}_sd_over_draws", f"{name}_lo", f"{name}_hi")
        ]
    return out.with_columns(gated).sort(KEYS)


def reference_rule(key: str) -> bool:
    """Whether cohort ``key``, as a reference to another cohort, gives its youngest fifth (True) or all donors."""
    rules = {ref: youngest for refs in PLAN.values() for name, (ref, youngest) in refs.items() if name != SAME}
    return rules[key]


def cohort_summary(
    cohorts: dict,
    cpgs: int,
    own: dict[str, set[str]] | None = None,
    mode: str = "deposited",
    recipes: dict[str, str] | None = None,
) -> list[dict]:
    """One row per cohort: arrays and ages, and the rows that give matching variables when it is another cohort's
    reference (the same ``reference_rows`` call the matching uses). With ``own`` (each cohort's CpGs passing its own
    missingness rule, before the intersection) also that count and how many CpGs the intersection of the other
    cohorts would hold without this one (what this cohort costs the common set)."""
    out = []
    for k, (p, raw) in cohorts.items():
        youngest = reference_rule(k)
        ref = reference_rows(p, np.arange(p.x.shape[0]), youngest)
        out.append(
            {
                "cohort": p.label,
                "key": k,
                "processing": (
                    f"{SESAME}; recipe {(recipes or {}).get(k, RECIPE)}"
                    if mode == "sesame"
                    else PROCESSING.get(k, "deposited")
                ),
                "arrays": p.x.shape[0],
                "age_min": float(p.ages.min()),
                "age_max": float(p.ages.max()),
                "reference_rule": "youngest fifth" if youngest else "all donors",
                "reference_arrays": int(ref.size),
                "reference_age_min": float(p.ages[ref].min()),
                "reference_age_max": float(p.ages[ref].max()),
                "cpgs": cpgs,
                "own_cpgs": len(own[k]) if own else None,
                "common_cpgs_without": len(set.intersection(*(v for j, v in own.items() if j != k))) if own else None,
                "missing_share": float(np.isnan(raw).mean()),
            }
        )
    return out


# ----------------------------------------------------------------------------- deposited against rebuilt
#: the two young-adult references read both ways (``processing`` table of a rebuilt run)
YOUNG_LABELS = {"young_adults": "GSE73103 young adults", "erisk": "E-Risk one twin per family"}
VERSIONS = ("deposited", "sesame")


def deposited_twin(paths, key: str, p: Prepared) -> np.ndarray:
    """The deposited betas of the donors (accessions) and CpGs of the rebuilt cohort ``p``, gaps kept."""
    donors, features = p.donors.tolist(), list(p.features)
    if key == "erisk":
        column = dict(erisk_people(matrix_path(paths, ERISK)).select("sample", "column").iter_rows())
        return read_erisk_betas(erisk_parts(paths), features, [column[d] for d in donors])
    return read_deposited(matrix_path(paths, REBUILT[key]), features, donors).astype(np.float32)


def paired_mask(versions: dict[tuple[str, str], np.ndarray]) -> dict[tuple[str, str], np.ndarray]:
    """Each cohort's versions with a cell kept only where every version of that cohort measures it (finite), so the
    two processings are compared on the same donor cells."""
    out = {}
    for cohort in {c for c, _ in versions}:
        keys = [k for k in versions if k[0] == cohort]
        both = np.logical_and.reduce([np.isfinite(versions[k]) for k in keys])
        out |= {k: np.where(both, versions[k], np.nan).astype(np.float32) for k in keys}
    return out


def shared_columns(versions: dict, features: list[str]) -> np.ndarray:
    """Columns measured (finite) in all but at most ``MAX_MISSING`` of the rows of every version (on ``features``)."""
    return np.logical_and.reduce([(~np.isfinite(x)).mean(axis=0) <= MAX_MISSING for x in versions.values()])


def young_sd_gaps(
    versions: dict[tuple[str, str], np.ndarray],
    rows: dict[str, np.ndarray],
    features: list[str],
    probes: pl.DataFrame,
    membership: pl.DataFrame,
) -> dict[tuple[str, str, str, str], float]:
    """Set-minus-matched young SD (island x young-level decile, the paper's matching; beta units) of every
    (cohort, processing) version over its cohort's ``rows``, per set, on the CpGs every version can describe. Two
    stratifications: ``own`` (deciles of that version's young level, each version rematched) and ``deposited``
    (the strata of the cohort's deposited version applied to both versions, so a rebuilt-minus-deposited change holds
    the matching fixed)."""
    annotation = probes.select(pl.col("cpg").alias("feature"), "cpg_island").with_columns(
        pl.col("cpg_island").fill_null(False)
    )
    frames = {}
    for (cohort, version), x in versions.items():
        mean, sd = reference_stats(x, rows[cohort])
        frames[(cohort, version)] = pl.DataFrame({"feature": features, "baseline": mean, "young_sd": sd})
    finite = np.logical_and.reduce(
        [np.isfinite(f.select("baseline", "young_sd").to_numpy()).all(axis=1) for f in frames.values()]
    )
    frames = {
        k: f.filter(pl.Series(finite))
        .join(annotation, on="feature", how="left", maintain_order="left")
        .with_columns(pl.col("cpg_island").fill_null(False))
        for k, f in frames.items()
    }
    own = {k: strata(f, by_type=False) for k, f in frames.items()}
    out = {}
    for group in SETS:
        members = membership.filter(pl.col("group") == group)["feature"].unique().implode()
        for (cohort, version), f in frames.items():
            in_group = f["feature"].is_in(members).to_numpy()
            fixed = own.get((cohort, "deposited"), own[(cohort, version)])
            for name, codes in (("own", own[(cohort, version)]), ("deposited", fixed)):
                out[(cohort, version, group, name)] = stratum_gap(f, in_group, codes, ("young_sd",))["young_sd"]
    return out


def _contrasts(g: dict[tuple[str, str, str, str], float]) -> dict[tuple[str, str, str, str], float]:
    """Each version's value, rebuilt minus deposited per cohort (paired), GSE73103 minus E-Risk per processing, and
    the change of that cohort difference on rebuilding; under each stratification."""
    out = dict(g)
    for group in SETS:
        for st in ("own", "deposited"):
            for c in YOUNG_LABELS:
                out[(c, "rebuilt minus deposited", group, st)] = (
                    g[(c, "sesame", group, st)] - g[(c, "deposited", group, st)]
                )
            for v in VERSIONS:
                out[("GSE73103 minus E-Risk", v, group, st)] = (
                    g[("young_adults", v, group, st)] - g[("erisk", v, group, st)]
                )
            out[("GSE73103 minus E-Risk", "rebuilt minus deposited", group, st)] = (
                out[("GSE73103 minus E-Risk", "sesame", group, st)]
                - out[("GSE73103 minus E-Risk", "deposited", group, st)]
            )
    return out


def processing_table(
    versions: dict[tuple[str, str], np.ndarray],
    features: list[str],
    probes: pl.DataFrame,
    membership: pl.DataFrame,
    n_boot: int,
    rng: np.random.Generator,
) -> pl.DataFrame:
    """``young_sd_gaps`` and their contrasts with 95 % donor-resample intervals: each cohort's donors drawn with
    replacement, the same draw for its two versions, independently between cohorts. In percentage points (beta x 100;
    below zero = the set varies less between donors than its matched CpGs)."""
    versions = paired_mask(versions)
    keep = shared_columns(versions, features)
    features = [f for f, k in zip(features, keep, strict=True) if k]
    versions = {k: np.ascontiguousarray(x[:, keep]) for k, x in versions.items()}
    n = {c: x.shape[0] for (c, _), x in versions.items()}
    point = _contrasts(young_sd_gaps(versions, {c: np.arange(m) for c, m in n.items()}, features, probes, membership))
    boots = []
    for _ in range(n_boot):
        rows = {c: np.sort(rng.choice(m, size=m, replace=True)) for c, m in n.items()}
        boots.append(_contrasts(young_sd_gaps(versions, rows, features, probes, membership)))
    out = []
    for key, value in point.items():
        cohort, version, group, stratification = key
        b = np.array([d[key] for d in boots], dtype=float)
        ok = np.isfinite(b)
        enough = b.size and ok.sum() >= MIN_VALID_SHARE * b.size
        out.append(
            {
                "cohort": YOUNG_LABELS.get(cohort, cohort),
                "processing": version,
                "group": group,
                "strata": stratification,
                "sd_diff_pp": 100 * value,
                "lo_pp": 100 * float(np.quantile(b[ok], 0.025)) if enough else None,
                "hi_pp": 100 * float(np.quantile(b[ok], 0.975)) if enough else None,
                "n_boot": int(b.size),
                "n_valid": int(ok.sum()),
                "donors": n.get(cohort),
                "cpgs": len(features),
            }
        )
    return pl.DataFrame(out)


def run(ctx: Context) -> list[Path]:
    mode = processing(ctx.paths)
    plan = active_plan(mode)
    out = ctx.out_dir("level_match_reference")
    workers, n_boot = _settings(ctx)
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))  # the background draw of `figures tissues`
    cohorts = {}
    if mode == "deposited":
        hannum = _hannum(ctx.paths, features)  # deposited Hannum betas have no gaps
        hannum = dataclasses.replace(hannum, x=hannum.x.astype(np.float32))
        cohorts["hannum"] = (hannum, hannum.x)
    cohorts |= {
        "epic_italy": _epic_italy(ctx.paths, features),
        "charge": _charge_controls(ctx.paths, features),
        "young_adults": _young_adults(ctx.paths, features),
        "erisk": _erisk(ctx.paths, features, ctx.rng(60_018)),
    }
    if mode == "sesame":  # the two young-adult cohorts read both ways, same donors, before the five-cohort cut
        versions = {}
        for key in YOUNG_LABELS:
            p, raw = cohorts[key]
            idx = {f: j for j, f in enumerate(p.features)}
            full = np.full((raw.shape[0], len(features)), np.nan, dtype=np.float32)
            cols = [j for j, f in enumerate(features) if f in idx]
            full[:, cols] = raw[:, [idx[features[j]] for j in cols]]
            versions[(key, "sesame")] = full
            deposited = dataclasses.replace(p, features=list(features))
            versions[(key, "deposited")] = deposited_twin(ctx.paths, key, deposited)
        write_table(
            processing_table(versions, features, probes, membership, n_boot, ctx.rng(60_020)), out / "processing"
        )
    own = {k: set(p.features) for k, (p, _) in cohorts.items()}
    common = sorted(set.intersection(*own.values()))
    cohorts = {k: restrict(p, common, raw) for k, (p, raw) in cohorts.items()}
    membership = membership.filter(pl.col("feature").is_in(common))
    recipes = {k: recipe(ctx.paths, REBUILT[k]) for k in cohorts} if mode == "sesame" else None
    summary = cohort_summary(cohorts, len(common), own, mode, recipes)
    write_table(pl.DataFrame(summary), out / "cohorts")
    print(f"level_match_reference ({mode}): {summary}", file=sys.stderr, flush=True)
    point, spread = [], []
    full = {k: np.arange(p.x.shape[0]) for k, (p, _) in cohorts.items()}
    for key in plan:
        rows, sd = sample_rows(key, cohorts, probes, membership, full, plan=plan)
        point.append(rows)
        spread.append(sd)
    point, spread = pl.concat(point), pl.concat(spread)
    bad = spread.filter(pl.col("sd").is_null() | ~pl.col("sd").is_finite() | (pl.col("sd") <= 0))
    if bad.height:
        raise ValueError(f"level_match_reference: no spread of a drift measure: {bad.to_dicts()}")
    rng = ctx.rng(60_016)
    jobs = [(key, int(s)) for key in plan for s in rng.integers(0, 2**31 - 1, size=n_boot)]
    state = {"cohorts": cohorts, "probes": probes, "membership": membership, "plan": plan}
    boot = pl.DataFrame(_parallel(jobs, workers, state)) if jobs else pl.DataFrame()
    table = summarise(point, boot, spread).with_columns(processing=pl.lit(mode))
    write_table(table, out / "gaps", csv_max_rows=CSV_ROWS)
    if boot.height:
        write_table(boot, out / "boot")
    sub_jobs = [(key, int(s)) for key in plan for s in ctx.rng(60_019).integers(0, 2**31 - 1, size=n_boot)]
    if sub_jobs:
        draws = pl.DataFrame(_parallel(sub_jobs, workers, state, _subsample))
        write_table(summarise_subsamples(draws, spread), out / "subsample")
    return []
