"""Matching variables from another cohort's young donors (planned test P6, the young-SD question).

``level_match`` found that matching on the young between-donor SD as well as the young level halves the blood
housekeeping gap (Hannum). Low young SD may itself be part of the protection: the same donors' young SD then carries
the set's slower drift, and matching on it removes some of the gap it should measure. Here the matching variables
(young mean level, young SD, deciles) come from the young donors of a different 450k blood cohort, so no donor and
no measurement enters both the matching variables and the drift:

* **Drift in Hannum** (GSE40279, whole blood, 656 people aged 19-101; the paper's blood estimate), against
  references from EPIC-Italy's youngest fifth (GSE51032, buffy coat, 34-47 y), from the CHARGE/Kabuki series'
  controls (GSE97362, whole blood, 124 children and adolescents aged 2-19) and from the young adults of GSE73103
  (Voisin 2015, whole blood of healthy young Swedes: the 216 aged 18-30, all weight categories).
* **Drift in EPIC-Italy** (all arrays with intact signal, 34-72 y, covariates sex, later cancer diagnosis and chip),
  against Hannum's youngest fifth (19-52 y), the GSE97362 controls and the GSE73103 young adults.
* ``same donors``: the drift cohort's own youngest fifth, as in the paper, on the same CpGs.

Every estimate is 450k, deposited betas (Hannum has no raw arrays), on the CpGs measured in all four cohorts
(``tissue_contrast``'s promoter features; at most 5 % missing per cohort, the rest median-imputed). The gaps are
``level_match.gaps`` (pool ``rest``): island x decile (the paper's matching), the finer strata, and the nearest five
on young level alone and on young level and young SD, in SD of the drift measure over all CpGs of the drift cohort
(one divisor per cohort and measure). Below zero means the set drifts less than its matched CpGs. Intervals resample
donors of the drift cohort and of the reference cohort independently (``ENTROPY_AGE_BOOT`` resamples,
``ENTROPY_AGE_WORKERS`` spawned processes); the young fifth and every stratum and neighbour are recomputed in each.

``entropy-age figures level_match_reference`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.incident_cancer import chips
from entropy_age.experiments.incident_cancer import read_people as epic_italy_people
from entropy_age.experiments.level_match import KEYS, MIN_VALID_SHARE, _finite_quantile, gaps, young_rows
from entropy_age.experiments.syndromes import CONTROL, matrix_path, read_betas, usable
from entropy_age.experiments.syndromes import read_people as syndrome_people
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    MAX_MISSING,
    Prepared,
    _features,
    _hannum,
    _lifespan_metrics,
    _settings,
)

EPIC_ITALY, CHARGE, YOUNG_ADULTS = "GSE51032", "GSE97362", "GSE73103"
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
    },
    "epic_italy": {
        SAME: ("epic_italy", True),
        "Hannum youngest fifth": ("hannum", True),
        "GSE97362 controls": ("charge", False),
        "GSE73103 young adults": ("young_adults", False),
    },
}
CSV_ROWS = 5_000
#: the drift cohort's own young mean and SD, carried as checks: after matching on another cohort's young donors,
#: how far do the set and its matched CpGs still differ in the drift cohort's own young donors?
OWN = ("own_level", "own_sd")


# ----------------------------------------------------------------------------- cohorts
def deposited_only() -> None:
    """Hannum exists only as deposited betas, so every cohort here is read deposited: one processing per estimate."""
    if os.environ.get("ENTROPY_AGE_BETAS", "deposited") != "deposited":
        raise RuntimeError("level_match_reference reads deposited betas only (Hannum has no raw arrays)")


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


def sample_rows(key: str, cohorts: dict, probes, membership, draws: dict[str, np.ndarray]):
    """Gaps of every reference for drift cohort ``key``; ``cohorts`` maps a key to (Prepared, matrix with gaps kept)
    and ``draws`` gives each cohort's rows (a resample or all). Also the SD of each drift measure over the drift
    cohort's CpGs (the divisor, used from the full sample only)."""
    drift, drift_raw = cohorts[key]
    references = {}
    for name, (ref_key, youngest) in PLAN[key].items():
        ref, raw = cohorts[ref_key]
        references[name] = (ref, raw, reference_rows(ref, draws[ref_key], youngest))
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
    rows, _ = sample_rows(key, cohorts, _STATE["probes"], _STATE["membership"], draws)
    return rows.with_columns(boot=pl.lit(seed)).to_dicts()


def _parallel(jobs, workers: int, state: dict) -> list[dict]:
    """Spawned workers (a forked polars pool deadlocks: ``tissue_contrast._run_parallel``)."""
    if workers <= 1:
        _init(state)
        return [r for j in jobs for r in _boot(j)]
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=_init, initargs=(state,)) as pool:
        return [r for rows in pool.map(_boot, jobs) for r in rows]


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


#: how each cohort's deposited betas were made (the depositors' processing; none rebuilt here)
PROCESSING = {
    "hannum": "deposited (depositors' processing)",
    "epic_italy": "deposited (depositors' processing)",
    "charge": "deposited (depositors' processing)",
    "young_adults": "deposited (NOOB + BMIQ; SNP, cross-reactive and sex-chromosome probes removed)",
}


def reference_rule(key: str) -> bool:
    """Whether cohort ``key``, as a reference to another cohort, gives its youngest fifth (True) or all donors."""
    rules = {ref: youngest for refs in PLAN.values() for name, (ref, youngest) in refs.items() if name != SAME}
    return rules[key]


def cohort_summary(cohorts: dict, cpgs: int) -> list[dict]:
    """One row per cohort: arrays and ages, and the rows that give matching variables when it is another cohort's
    reference (the same ``reference_rows`` call the matching uses)."""
    out = []
    for k, (p, raw) in cohorts.items():
        youngest = reference_rule(k)
        ref = reference_rows(p, np.arange(p.x.shape[0]), youngest)
        out.append(
            {
                "cohort": p.label,
                "key": k,
                "processing": PROCESSING.get(k, "deposited"),
                "arrays": p.x.shape[0],
                "age_min": float(p.ages.min()),
                "age_max": float(p.ages.max()),
                "reference_rule": "youngest fifth" if youngest else "all donors",
                "reference_arrays": int(ref.size),
                "reference_age_min": float(p.ages[ref].min()),
                "reference_age_max": float(p.ages[ref].max()),
                "cpgs": cpgs,
                "missing_share": float(np.isnan(raw).mean()),
            }
        )
    return out


def run(ctx: Context) -> list[Path]:
    deposited_only()
    out = ctx.out_dir("level_match_reference")
    workers, n_boot = _settings(ctx)
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))  # the background draw of `figures tissues`
    hannum = _hannum(ctx.paths, features)  # deposited Hannum betas have no gaps
    hannum = dataclasses.replace(hannum, x=hannum.x.astype(np.float32))
    cohorts = {
        "hannum": (hannum, hannum.x),
        "epic_italy": _epic_italy(ctx.paths, features),
        "charge": _charge_controls(ctx.paths, features),
        "young_adults": _young_adults(ctx.paths, features),
    }
    common = sorted(set.intersection(*(set(p.features) for p, _ in cohorts.values())))
    cohorts = {k: restrict(p, common, raw) for k, (p, raw) in cohorts.items()}
    membership = membership.filter(pl.col("feature").is_in(common))
    write_table(pl.DataFrame(cohort_summary(cohorts, len(common))), out / "cohorts")
    print(f"level_match_reference: {cohort_summary(cohorts, len(common))}", file=sys.stderr, flush=True)
    point, spread = [], []
    full = {k: np.arange(p.x.shape[0]) for k, (p, _) in cohorts.items()}
    for key in PLAN:
        rows, sd = sample_rows(key, cohorts, probes, membership, full)
        point.append(rows)
        spread.append(sd)
    point, spread = pl.concat(point), pl.concat(spread)
    bad = spread.filter(pl.col("sd").is_null() | ~pl.col("sd").is_finite() | (pl.col("sd") <= 0))
    if bad.height:
        raise ValueError(f"level_match_reference: no spread of a drift measure: {bad.to_dicts()}")
    rng = ctx.rng(60_016)
    jobs = [(key, int(s)) for key in PLAN for s in rng.integers(0, 2**31 - 1, size=n_boot)]
    state = {"cohorts": cohorts, "probes": probes, "membership": membership}
    boot = pl.DataFrame(_parallel(jobs, workers, state)) if jobs else pl.DataFrame()
    table = summarise(point, boot, spread).with_columns(processing=pl.lit("deposited"))
    write_table(table, out / "gaps", csv_max_rows=CSV_ROWS)
    if boot.height:
        write_table(boot, out / "boot")
    return []
