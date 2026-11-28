"""Finer level matching in the lowest decile (planned test P6; claim C1, alternative (a), the floor).

Housekeeping and core-essential promoter CpGs sit near 0 % methylation. The headline matching (island status x
decile of the young level, ``tissue_contrast``) keeps the lowest decile as one stratum, so inside it a set's CpGs can
still start lower than their look-alikes, and a lower start alone gives a smaller age drift on a bounded scale. This
module measures how much of the gap that could be:

* **Residual level after matching.** For every matching, the set's young mean level (and young SD) minus that of its
  matched controls, weighted as the gap is (signed), and the mean absolute residual level per set CpG.
* **Finer strata in the lowest decile.** Island x probe design x decile, the lowest decile cut into 20 or 50 quantile
  bins of the young level (the other deciles as before).
* **Nearest neighbours.** Each set CpG against its ``NN_K`` nearest controls of the same island status and probe
  design, on the young mean level alone, and on young mean level and young SD (each scaled by its SD over the
  scope's CpGs); controls with replacement. Young SD may carry early drift itself, so the level-only row isolates
  what the level does.
* **Scopes.** The whole promoter universe (the gap the paper reports), the lowest decile alone (set and controls
  both inside it), and floors: set and controls both above a young mean of 2 % or 5 % (strata recomputed there).
* **Reference donors.** ``same donors`` uses the youngest fifth of the cohort for the matching variables, as the
  paper does; ``held-out young half`` takes the matching variables (young mean, young SD, deciles, floors) from a
  random half of the young donors and fits the drift on the other donors only, so no donor enters both (the shared
  young donors could otherwise tie the matching variable to the drift through regression to the mean).
* **Control pool.** ``rest``: every CpG outside the set (the paper's estimator; other gene sets included);
  ``non-set``: only CpGs in none of the three sets.

Gaps are in SD of the drift measure over all CpGs of the full sample for that reference split (one divisor for every
row and every resample); below zero means the set drifts less than its matched CpGs. Each resample also gives each
matching's paired change against the scope's coarse matching. Hannum blood (deposited betas; no raw IDATs exist) and
Jaffe cortex (deposited, or the R SeSAMe rebuild under ``ENTROPY_AGE_BETAS=sesame``), 450k each, with
``tissue_contrast``'s features (the background draw of ``figures tissues``), covariates and drift measures.
Intervals resample donors (``ENTROPY_AGE_BOOT`` resamples, ``ENTROPY_AGE_WORKERS`` spawned processes); strata,
neighbours and the held-out split are rebuilt in every resample, copies of a donor kept on one side.

``entropy-age figures level_match`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import multiprocessing
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from scipy.spatial import cKDTree

from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    SETS,
    Prepared,
    _features,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _settings,
)
from entropy_age.stochasticity import _matched_difference, baseline_strata

FINE_BINS = (20, 50)
NN_K = 5
FLOORS = (0.02, 0.05)
YOUNG_QUANTILE = 0.2  # drift_metrics' young reference: the youngest fifth by age
LEVEL_COLUMNS = ("baseline", "young_sd")
REFERENCES = ("same donors", "held-out young half")
POOLS = ("rest", "non-set")
MIN_VALID_SHARE = 0.5  # an interval only when at least half the resamples are finite (the series rule)
CSV_ROWS = 5_000


# ----------------------------------------------------------------------------- per-CpG table
def young_rows(ages: np.ndarray) -> np.ndarray:
    """Positions of the youngest fifth by age (drift_metrics' young reference)."""
    return np.flatnonzero(ages <= np.quantile(ages, YOUNG_QUANTILE))


def young_sd(x: np.ndarray, ages: np.ndarray) -> np.ndarray:
    """SD across the young reference's samples (the same youngest fifth that sets ``baseline``)."""
    return x[young_rows(ages)].std(axis=0, ddof=1)


def held_out_split(ages: np.ndarray, donors: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """(reference rows, drift rows): a random half of the young donors, and every row of the other donors.

    Young is the youngest fifth of the rows passed. Rows of one donor (copies in a resample) stay on one side.
    """
    young = np.unique(donors[young_rows(ages)])
    picked = rng.choice(young, size=young.size // 2, replace=False)
    ref = np.isin(donors, picked)
    return np.flatnonzero(ref), np.flatnonzero(~ref)


def cpg_table(
    p: Prepared, probes: pl.DataFrame, rows: np.ndarray | None = None, reference: np.ndarray | None = None
) -> pl.DataFrame:
    """Drift measures per CpG with young mean (``baseline``), young SD, island status and probe design.

    Drift comes from ``rows``. Without ``reference`` the matching variables come from the youngest fifth of ``rows``
    (as in the paper); with it, from the rows ``reference`` (mean and SD over them), which the caller keeps apart.
    """
    rows = np.arange(p.x.shape[0]) if rows is None else rows
    m = _lifespan_metrics(p, rows)
    if reference is None:
        m = m.with_columns(young_sd=pl.Series(young_sd(p.x[rows], p.ages[rows])))
    else:
        x = p.x[reference].astype(float)
        m = m.with_columns(baseline=pl.Series(x.mean(axis=0)), young_sd=pl.Series(x.std(axis=0, ddof=1)))
    annotation = probes.select(pl.col("cpg").alias("feature"), "cpg_island", "probe_type")
    m = m.join(annotation, on="feature", how="left")
    if m["probe_type"].null_count():
        raise ValueError("level_match: a CpG has no probe design")
    return m.with_columns(pl.col("cpg_island").fill_null(False))


# ----------------------------------------------------------------------------- strata
def lowest_decile(baseline: np.ndarray) -> np.ndarray:
    """True for CpGs in the lowest decile of the young level over the frame passed."""
    return baseline_strata(baseline, 10) == 0


def strata(frame: pl.DataFrame, *, by_type: bool, fine: int = 0) -> np.ndarray:
    """Island [x probe design] x young-level decile; with ``fine`` the lowest decile is cut into that many bins.

    The fine bins are quantiles of the young level inside the lowest decile (pooled over island and design; tied
    edges merge bins).
    """
    level = frame["baseline"].to_numpy()
    decile = baseline_strata(level, 10)
    sub = np.zeros(level.size, dtype=int)
    if fine:
        low = decile == 0
        sub[low] = baseline_strata(level[low], fine)
    parts = [frame["cpg_island"].cast(pl.String).to_numpy(), decile.astype(str), sub.astype(str)]
    if by_type:
        parts.insert(1, frame["probe_type"].to_numpy())
    return _codes(parts)


def within_decile(frame: pl.DataFrame, *, by_type: bool, fine: int) -> np.ndarray:
    """Strata inside one decile: island [x probe design] [x ``fine`` quantile bins of the young level]."""
    parts = [frame["cpg_island"].cast(pl.String).to_numpy()]
    if by_type:
        parts.append(frame["probe_type"].to_numpy())
    if fine:
        parts.append(baseline_strata(frame["baseline"].to_numpy(), fine).astype(str))
    return _codes(parts)


def _codes(parts: list[np.ndarray]) -> np.ndarray:
    if not parts[0].size:
        return np.zeros(0, dtype=int)
    keys = np.array(["|".join(t) for t in zip(*parts, strict=True)])
    return np.unique(keys, return_inverse=True)[1]


def stratum_gap(frame: pl.DataFrame, in_group: np.ndarray, codes: np.ndarray, columns) -> dict:
    """Set minus controls within strata, weighted by set CpGs per stratum, for each of ``columns``.

    Strata that hold only one side are dropped (``compare_groups``' rule); ``n_matched`` counts the set CpGs left.
    ``abs_level_diff`` is the mean over set CpGs of |young level - its stratum's control mean|.
    """
    nan = {c: float("nan") for c in columns} | {"abs_level_diff": float("nan")}
    if not codes.size:
        return {"n_matched": 0} | nan
    n_strata = int(codes.max()) + 1
    n_in = np.bincount(codes, weights=in_group.astype(float), minlength=n_strata)
    n_all = np.bincount(codes, minlength=n_strata)
    mixed = np.flatnonzero((n_in > 0) & (n_in < n_all))
    keep = np.isin(codes, mixed)
    out = {"n_matched": int(in_group[keep].sum())}
    if in_group[keep].sum() < 2 or (~in_group[keep]).sum() < 2:
        return out | nan
    _, recoded = np.unique(codes[keep], return_inverse=True)
    k = int(recoded.max()) + 1
    for c in columns:
        out[c] = _matched_difference(frame[c].to_numpy()[keep], in_group[keep], recoded, k)
    level, member = frame["baseline"].to_numpy()[keep], in_group[keep]
    control_mean = np.bincount(recoded[~member], weights=level[~member], minlength=k) / np.maximum(
        np.bincount(recoded[~member], minlength=k), 1
    )
    out["abs_level_diff"] = float(np.abs(level[member] - control_mean[recoded[member]]).mean())
    return out


def nn_gap(
    frame: pl.DataFrame,
    in_group: np.ndarray,
    columns,
    k: int = NN_K,
    scale: dict | None = None,
    on: tuple[str, ...] = LEVEL_COLUMNS,
) -> dict:
    """Each set CpG minus the mean of its ``k`` nearest controls of the same island status and probe design.

    Distance on the columns ``on``, each divided by ``scale[column]`` (default: its SD over the frame; a column with
    no spread is left out of the distance). Controls are drawn with replacement. Set CpGs in a cell with fewer than
    ``k`` controls are left out. ``abs_level_diff`` is the mean |young level - its neighbours' mean level|.
    """
    scale = scale or {c: float(frame[c].std() or 0.0) for c in on}
    on = tuple(c for c in on if np.isfinite(scale.get(c, np.nan)) and scale[c] > 0)
    nan = {"n_matched": 0, **{c: float("nan") for c in columns}, "abs_level_diff": float("nan")}
    if not on or not frame.height:
        return nan
    z = np.column_stack([frame[c].to_numpy() / scale[c] for c in on])
    values = {c: frame[c].to_numpy() for c in {*columns, "baseline"}}
    cell = _codes([frame["cpg_island"].cast(pl.String).to_numpy(), frame["probe_type"].to_numpy()])
    diffs: dict[str, list] = {c: [] for c in {*columns, "baseline"}}
    for code in np.unique(cell):
        s = np.flatnonzero((cell == code) & in_group)
        r = np.flatnonzero((cell == code) & ~in_group)
        if s.size == 0 or r.size < k:
            continue
        _, idx = cKDTree(z[r]).query(z[s], k=k)
        idx = r[np.asarray(idx).reshape(s.size, k)]  # tree positions -> frame rows
        for c in diffs:
            diffs[c].append(values[c][s] - values[c][idx].mean(axis=1))
    if not diffs["baseline"]:
        return nan
    flat = {c: np.concatenate(v) for c, v in diffs.items()}
    out = {"n_matched": int(flat["baseline"].size)} | {c: float(flat[c].mean()) for c in columns}
    return out | {"abs_level_diff": float(np.abs(flat["baseline"]).mean())}


# ----------------------------------------------------------------------------- one sample
MATCHINGS = {  # scope "all CpGs"; the first is the paper's matching and the reference for paired changes
    "island x decile": {"by_type": False, "fine": 0},
    "island x probe x decile": {"by_type": True, "fine": 0},
    "lowest decile in 20 bins": {"by_type": True, "fine": 20},
    "lowest decile in 50 bins": {"by_type": True, "fine": 50},
}
WITHIN = {  # scope "lowest decile": strata inside the decile
    "island": {"by_type": False, "fine": 0},
    "island x probe": {"by_type": True, "fine": 0},
    "island x probe x 20 bins": {"by_type": True, "fine": 20},
    "island x probe x 50 bins": {"by_type": True, "fine": 50},
}
NN_PLANS = {f"nearest {NN_K} (level)": ("baseline",), f"nearest {NN_K} (level, SD)": LEVEL_COLUMNS}


def gaps(frame: pl.DataFrame, membership: pl.DataFrame, pools=POOLS, extra: tuple[str, ...] = ()) -> list[dict]:
    """Unscaled gaps and residual level differences of every pool, scope, set and matching for one sample's table;
    ``extra`` columns get the same set-minus-matched difference (never used for matching)."""
    columns = (*LIFESPAN_METRICS, *LEVEL_COLUMNS, *extra)
    frame = frame.filter(pl.all_horizontal(pl.col(c).is_finite() for c in columns))
    level = frame["baseline"].to_numpy()
    low = lowest_decile(level)
    top = float(level[low].max()) if low.any() else float("nan")
    any_set = frame["feature"].is_in(membership.filter(pl.col("group").is_in(SETS))["feature"].unique().implode())
    rows = []
    for pool in pools:
        for scope, keep in {
            "all CpGs": np.ones(frame.height, dtype=bool),
            "lowest decile": low,
            **{f"young level > {f:.2f}": level > f for f in FLOORS},
        }.items():
            f_scope = frame.filter(pl.Series(keep))
            for group in SETS:
                members = membership.filter(pl.col("group") == group)["feature"].unique().implode()
                f = f_scope
                if pool == "non-set":  # controls only from CpGs in no set
                    f = f_scope.filter(f_scope["feature"].is_in(members) | ~any_set.filter(pl.Series(keep)))
                in_group = f["feature"].is_in(members).to_numpy()
                if scope == "lowest decile":
                    codes = {name: within_decile(f, **plan) for name, plan in WITHIN.items()}
                elif scope == "all CpGs":
                    codes = {name: strata(f, **plan) for name, plan in MATCHINGS.items()}
                else:
                    codes = {"island x decile": strata(f, **MATCHINGS["island x decile"])}
                head = {"pool": pool, "scope": scope, "group": group, "n": int(in_group.sum()), "decile_top": top}
                for name, c in codes.items():
                    rows.append(head | {"matching": name} | stratum_gap(f, in_group, c, columns))
                scale = {c: float(f[c].std() or 0.0) for c in LEVEL_COLUMNS}
                for name, on in NN_PLANS.items():
                    rows.append(head | {"matching": name} | nn_gap(f, in_group, columns, NN_K, scale, on=on))
    return rows


def long_form(rows: list[dict]) -> pl.DataFrame:
    """One row per (pool, scope, matching, group, measure): the drift gap and the residual young level and SD."""
    wide = pl.DataFrame(rows)
    keys = ["pool", "scope", "matching", "group", "n", "n_matched", "decile_top", "baseline", "young_sd"]
    keys.append("abs_level_diff")
    return wide.unpivot(
        index=keys, on=list(LIFESPAN_METRICS), variable_name="metric", value_name="matched_diff"
    ).rename({"baseline": "level_diff", "young_sd": "sd_diff"})


def sample_rows(
    p: Prepared, probes, membership, rows: np.ndarray, rng: np.random.Generator
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Both reference splits for the rows ``rows`` (a full sample or a resample), and each table's SD per measure."""
    ref, drift = held_out_split(p.ages[rows], p.donors[rows], rng)
    tables = {
        REFERENCES[0]: cpg_table(p, probes, rows),
        REFERENCES[1]: cpg_table(p, probes, rows[drift], reference=rows[ref]),
    }
    out = [long_form(gaps(t, membership)).with_columns(reference=pl.lit(r)) for r, t in tables.items()]
    spread = []
    for r, t in tables.items():
        for m in LIFESPAN_METRICS:
            v = t[m].filter(t[m].is_finite())
            spread.append(
                {"tissue": p.label, "reference": r, "metric": m, "sd": float(v.std()) if v.len() > 1 else None}
            )
    return pl.concat(out).with_columns(tissue=pl.lit(p.label)), pl.DataFrame(spread)


# ----------------------------------------------------------------------------- bootstrap
_STATE: dict = {}


def _init(state: dict) -> None:
    _STATE.clear()
    _STATE.update(state)


def _boot(job: tuple[str, int]) -> list[dict]:
    key, seed = job
    p, probes, membership = _STATE["prepared"][key], _STATE["probes"], _STATE["membership"]
    rng = np.random.default_rng(seed)
    draw = np.sort(rng.choice(p.x.shape[0], size=p.x.shape[0], replace=True))
    rows, _ = sample_rows(p, probes, membership, draw, rng)
    return rows.with_columns(boot=pl.lit(seed)).to_dicts()


def _parallel(jobs, workers: int, state: dict) -> list[dict]:
    """Spawned workers (a forked polars pool deadlocks: ``tissue_contrast._run_parallel``)."""
    if workers <= 1:
        _init(state)
        return [r for j in jobs for r in _boot(j)]
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=_init, initargs=(state,)) as pool:
        return [r for rows in pool.map(_boot, jobs) for r in rows]


KEYS = ["tissue", "reference", "pool", "scope", "matching", "group", "metric"]


def _with_delta(frame: pl.DataFrame) -> pl.DataFrame:
    """``delta_sd``: the gap minus the scope's first (coarsest) matching's gap in the same sample."""
    first = {"all CpGs": "island x decile", "lowest decile": "island"}
    ref = pl.col("scope").replace_strict(first, default="island x decile")
    base = (
        frame.filter(pl.col("matching") == ref)
        .select(*[k for k in KEYS if k != "matching"], *(["boot"] if "boot" in frame.columns else []), "gap_sd")
        .rename({"gap_sd": "coarse_sd"})
    )
    on = [k for k in KEYS if k != "matching"] + (["boot"] if "boot" in frame.columns else [])
    return frame.join(base, on=on, how="left").with_columns(delta_sd=pl.col("gap_sd") - pl.col("coarse_sd"))


def _finite_quantile(col: str, q: float) -> pl.Expr:
    return pl.col(col).filter(pl.col(col).is_finite()).quantile(q)


def summarise(point: pl.DataFrame, boot: pl.DataFrame, spread: pl.DataFrame) -> pl.DataFrame:
    """Point gaps on the all-CpG scale with 95 % donor-resample intervals of the gap, its paired change against the
    coarse matching, and the residual level; intervals only from finite resamples and only when at least half are."""
    scaled = pl.col("matched_diff") / pl.col("sd")
    on = ["tissue", "reference", "metric"]
    point = _with_delta(point.join(spread, on=on, how="left").with_columns(gap_sd=scaled))
    if not boot.height:
        return point.sort(KEYS)
    b = _with_delta(boot.join(spread, on=on, how="left").with_columns(gap_sd=scaled))
    ci = b.group_by(KEYS).agg(
        n_boot=pl.len(),
        n_valid=pl.col("gap_sd").is_finite().sum(),
        gap_lo=_finite_quantile("gap_sd", 0.025),
        gap_hi=_finite_quantile("gap_sd", 0.975),
        delta_lo=_finite_quantile("delta_sd", 0.025),
        delta_hi=_finite_quantile("delta_sd", 0.975),
        level_lo=_finite_quantile("level_diff", 0.025),
        level_hi=_finite_quantile("level_diff", 0.975),
    )
    thin = pl.col("n_valid") < MIN_VALID_SHARE * pl.col("n_boot")
    ci = ci.with_columns(
        [pl.when(thin).then(None).otherwise(pl.col(c)).alias(c) for c in ci.columns if c.endswith(("_lo", "_hi"))]
    )
    return point.join(ci, on=KEYS, how="left").sort(KEYS)


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("level_match")
    workers, n_boot = _settings(ctx)
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))  # the background draw of `figures tissues`
    prepared, point, spread = {}, [], []
    split_rng = ctx.rng(60_007)
    for key, name, loader in (("blood", "hannum", _hannum), ("cortex", "jaffe", _jaffe)):
        p = loader(ctx.paths, features)
        p = dataclasses.replace(p, x=p.x.astype(np.float32))
        prepared[key] = p
        rows, sd = sample_rows(p, probes, membership, np.arange(p.x.shape[0]), split_rng)
        point.append(rows.with_columns(processing=pl.lit(cohort_layout(name, ctx.paths).processing)))
        spread.append(sd)
        print(f"level_match: {p.label}, {p.x.shape[0]} donors", file=sys.stderr, flush=True)
    point, spread = pl.concat(point), pl.concat(spread)
    bad = spread.filter(pl.col("sd").is_null() | ~pl.col("sd").is_finite() | (pl.col("sd") <= 0))
    if bad.height:  # the divisor of every gap: refuse rather than divide by nothing
        raise ValueError(f"level_match: no spread of a drift measure: {bad.to_dicts()}")
    rng = ctx.rng(60_006)
    jobs = [(key, int(s)) for key in prepared for s in rng.integers(0, 2**31 - 1, size=n_boot)]
    state = {"prepared": prepared, "probes": probes, "membership": membership}
    boot = pl.DataFrame(_parallel(jobs, workers, state)) if jobs else pl.DataFrame()
    write_table(summarise(point, boot, spread), out / "gaps", csv_max_rows=CSV_ROWS)
    if boot.height:
        write_table(boot, out / "boot")
    return []
