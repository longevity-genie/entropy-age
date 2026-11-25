"""Clock disagreement as a signal (RESEARCH_PLAN idea 4): does the spread between a person's valid age clocks rise
with drift, and does it carry risk beyond their mean?

Panels. The eligible clocks are those trained to estimate chronological age whose output is in years, listed in
``CHRONOLOGICAL``: Horvath's two, Hannum, Lin, Vidal-Bralo, Weidner, Garagnani, Bocklandt, Ying's CausAge, DamAge and
AdaptAge (elastic nets on age over causal, damaging or adaptive CpGs), PEDBE (buccal) and the cortical clock.
Phenotypic and mortality targets (PhenoAge, HRS PhenoAge, GrimAge; years of a phenotype) and risks (Zhang_10, whose
output is a mortality risk), paces, mitotic and stochastic clocks are left out, so the spread compares estimates of
one quantity. Four panels per group, each reported with its exact clock list:

* ``valid`` (primary): the eligible clocks valid in the group's tissue (``clock_atlas.valid``);
* ``with clocks outside their tissue``: also the eligible blood clocks in solid tissue, so groups where only
  Horvath's two clocks are valid still have a spread;
* ``fixed shared: Horvath 1 and 2``: the same two pan-tissue clocks in every group with no in-sample age screen, so
  groups are compared on one panel (with two clocks the SD is the absolute difference of their z-scores over √2);
* ``with phenotypic-age clocks``: the valid panel plus PhenoAge, HRS PhenoAge and GrimAge where valid, beside the
  earlier panel definition.

In the first, second and fourth panels a clock enters only when its prediction tracks age in the group (r > 0.5, the
guard ``clock_increment`` uses: a single-CpG clock such as Bocklandt runs against age on some arrays); the screen
depends on the group's age range, which the fixed panel avoids. A clock whose residual on age (and on the quality
terms below) has near-zero variance is left out before z-scoring. Per clock, the acceleration is that residual,
z-scored within the group. Per person, the disagreement is the SD across the panel's clocks of these z-scores, and
the mean acceleration their mean.

(a) Drift. Per series group of ``config/series``, the arrays ``series.selection`` chooses: passing QC (detection and
the intensity rule), the description's donor exclusions, present in both the masked betas and the clock scores,
then one array per donor (so a donor whose first array lacks a score keeps another), and one platform per group;
the description's ``exclude_cpgs`` are left out of B and F. Groups need at least 25 arrays. Each person's drift at
housekeeping look-alike promoters (B) and shield share (F) against a reference fitted on everyone else in the group,
on the CpGs measured in every array (``loo_scores``, as ``clock_drift``). Per group and panel: the correlation of
the disagreement with B and with F, the same for the mean acceleration, the partial correlation of the disagreement
with B and F given the mean (and given age and the mean), and the quality sensitivity: accelerations residualised on
age plus array chip, log mean intensity and detection fraction (``reprocessed.read_qc``), the panel's mean measured
weight share per person and, where stored, 11 cell fractions, then the partial correlation of that disagreement with
B and F given age, its mean and the same terms. A missing intensity column or value is reported as ``unknown`` and
the term left out; chip is dropped when it would leave under 10 residual degrees of freedom. Rows from sensitivity
configs (``series.SENSITIVITY``, the tumour-field reruns of the lung groups) are marked ``sensitivity`` and, with the
EPIC-Italy women (a subset of everyone), kept out of the summary across groups (``drift_summary``), which takes
medians over the remaining groups. Groups of one series can share donors (brain regions), so they are separate
groups, not independent samples. EPIC-Italy (GSE51032, rebuilt; B and F from ``incident_cancer_rsesame``,
leave-one-out over everyone) and SJLIFE (GSE314261; B and F against one reference on everyone, as ``sjlife_deaths``)
are added as blood groups.

(b) Risk. EPIC-Italy women on arrays passing QC, with disagreement and mean computed among the women: breast cancer
(C50) diagnosed more than 2 years after the draw (and, as a check, at any time) against women with no cancer
diagnosis; the estimand is the odds of an incident diagnosis in the cohort's case and control women. Unconditional
logistic models of case status on disagreement and mean acceleration together, each per SD in the modelled women,
adjusted for age ("age"), also for mean beta and 11 cell fractions ("covariates"), also for log intensity and
detection ("covariates + quality"); both measures centred within chip ("chip-centred"); chip intercepts on the chips
holding at least one case and one control ("chip effects", retention reported); and conditional models on
age-matched sets, same chip first and same chip only (``cancer_matched.match``; retention and the share of controls
on their case's chip reported). Every fit is checked: one outcome only, a design short of full rank, a coefficient
beyond ±15 on the log-odds scale (separation), a Hessian with condition number above 1e10, or a score not at zero
gives an unavailable estimate (NaN) with the reason in ``fit``. Wald 95 % intervals.

SJLIFE: death among survivors (no control has died), from the clock accelerations ``sjlife_deaths`` stored (Horvath
2013, Hannum and Zhang 2019 in the primary panel, PhenoAge in the phenotypic one; already residualised on age, so the
age screen is not applied), adjusted for age, sex, the four main exposures and diagnosis; all deaths and deaths more
than 2 years after the draw. Survivors whose diagnosis has no death among the eligible are left out of the fit
(their intercept would separate; eligible and modelled counts are both reported). The estimand is the odds of a
recorded death by the data freeze, on GEO's deposited betas (``sjlife_deaths`` streams the processed table; no SeSAMe
rebuild and no array QC), as ``provenance`` says. GEO gives no follow-up or censoring dates, so the > 2 y analysis
cannot exclude living people observed for less than 2 years: they enter as controls, as ``followup`` says in every
SJLIFE row. Mean beta, cell fractions and array quality
are not stored for SJLIFE, so it has no covariate or quality version.

Writes ``clock_disagreement/drift`` (per group and panel), ``drift_summary``, ``per_person``, ``risk``,
``epic_italy_women`` and a two-panel figure. ``entropy-age figures clock_disagreement`` (a Slurm job: leave-one-out).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.cancer_matched import _score, clogit, match
from entropy_age.experiments.clock_atlas import TISSUE_ONLY, valid
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.incident_cancer import _group_means, logit
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import SENSITIVITY, selection, specs
from entropy_age.experiments.sjlife_deaths import EXPOSURES
from entropy_age.experiments.tissue_contrast import _features

SERIES = Path("data/interim/clocks/series")
CLOCK_SCORES = Path("data/interim/clocks/series_unmasked")
EPIC_ITALY = "GSE51032"
SJLIFE_GSE = "GSE314261"
GROUP = "housekeeping"
B, F = f"B {GROUP}", f"F {GROUP}"
MIN_ARRAYS = 25
MIN_CLOCKS = 2
MIN_AGE_R = 0.5
MIN_RESID_DF = 10
#: a clock's residual SD at or below this share of its SD (at least 1 unit) is treated as constant
RESID_TOL = 1e-6
LAG = 2.0
#: fit diagnostics: separation and a singular Hessian give an unavailable estimate
MAX_COEF = 15.0
MAX_CONDITION = 1e10
#: clocks trained on chronological age with output in years (biolearn names)
CHRONOLOGICAL = (
    "Horvathv1",
    "Horvathv2",
    "Hannum",
    "Lin",
    "VidalBralo",
    "Weidner",
    "Garagnani",
    "Bocklandt",
    "YingCausAge",
    "YingDamAge",
    "YingAdaptAge",
    "PEDBE",
    "DNAmClockCortical",
)
#: output in years, target a phenotype or mortality: only in the labelled sensitivity panel
PHENOTYPIC_YEARS = ("PhenoAge", "HRSInCHPhenoAge", "GrimAgeV1", "GrimAgeV2")
FIXED_PANEL = ("Horvathv1", "Horvathv2")
PANELS = ("valid", "with clocks outside their tissue", "fixed shared: Horvath 1 and 2", "with phenotypic-age clocks")
#: SJLIFE's stored accelerations (sjlife_deaths.CLOCKS_USED): chronological (Zhang 2019 is an elastic net on age);
#: PhenoAge only in the phenotypic panel; DunedinPACE (a pace) never
SJLIFE_CLOCKS = ("Horvath 2013", "Hannum", "Zhang 2019")
SJLIFE_PHENOTYPIC = ("PhenoAge",)
SJLIFE_ESTIMAND = "odds of a recorded death by the data freeze"
SJLIFE_FOLLOWUP = (
    "no follow-up or censoring dates in the table: living people enter the > 2 y analysis without evidence of 2 years "
    "of observation"
)
EPIC_ESTIMAND = "odds of incident breast cancer, cases against cancer-free women of the cohort"
EPIC_PROVENANCE = "SeSAMe rebuild of the IDATs, arrays passing detection and intensity QC"
FRACTIONS = [f"fraction_{i}" for i in range(1, 12)]


def in_years(unit: str | None) -> bool:
    return "(years)" in str(unit or "").strip().lower()


def age_clocks(outputs: dict[str, str], group: str, panel: str = PANELS[0]) -> list[str]:
    """The clocks of ``outputs`` (clock -> output unit) in ``panel`` for ``group``, in the order of ``outputs``."""
    years = [c for c, unit in outputs.items() if in_years(unit)]
    if panel == PANELS[2]:
        return list(FIXED_PANEL) if all(c in years for c in FIXED_PANEL) else []
    if panel == PANELS[0]:
        return [c for c in years if c in CHRONOLOGICAL and valid(c, group)]
    if panel == PANELS[1]:
        return [c for c in years if c in CHRONOLOGICAL and (valid(c, group) or c not in TISSUE_ONLY)]
    if panel == PANELS[3]:
        return [c for c in years if (c in CHRONOLOGICAL or c in PHENOTYPIC_YEARS) and valid(c, group)]
    raise ValueError(f"clock_disagreement: unknown panel {panel!r}")


def disagreement(
    preds: np.ndarray,
    age: np.ndarray,
    *,
    min_age_r: float | None = MIN_AGE_R,
    covariates: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per person (rows of ``preds``, people x clocks): the SD across clocks of the z-scored accelerations (residual
    on age and ``covariates``), their mean, and the mask of clocks used (finite, r with age above ``min_age_r``,
    residual variance not near zero)."""
    preds = np.asarray(preds, dtype=float)
    age = np.asarray(age, dtype=float)
    parts = [np.ones(age.size), age]
    if covariates is not None and covariates.size:
        parts.append(covariates)
    design = np.column_stack(parts)
    use = np.zeros(preds.shape[1], dtype=bool)
    for j in range(preds.shape[1]):
        v = preds[:, j]
        if not np.isfinite(v).all() or np.std(v) == 0:
            continue
        if min_age_r is not None and not np.corrcoef(v, age)[0, 1] > min_age_r:
            continue
        resid = v - design @ np.linalg.lstsq(design, v, rcond=None)[0]
        use[j] = np.std(resid) > RESID_TOL * max(1.0, float(np.std(v)))
    if use.sum() < MIN_CLOCKS:
        nan = np.full(age.size, np.nan)
        return nan, nan, use
    coef, *_ = np.linalg.lstsq(design, preds[:, use], rcond=None)
    acc = preds[:, use] - design @ coef
    acc = (acc - acc.mean(axis=0)) / acc.std(axis=0)
    return acc.std(axis=1, ddof=1), acc.mean(axis=1), use


def zs(v: np.ndarray) -> np.ndarray:
    return (v - v.mean()) / v.std(ddof=1)


def partial_r(x: np.ndarray, y: np.ndarray, given: np.ndarray) -> float:
    d = np.column_stack([np.ones(len(x)), given])
    rx = x - d @ np.linalg.lstsq(d, x, rcond=None)[0]
    ry = y - d @ np.linalg.lstsq(d, y, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def quality_terms(
    samples: list[str],
    qc: pl.DataFrame | None,
    *,
    share: np.ndarray | None = None,
    extra: dict[str, np.ndarray] | None = None,
    base_columns: int = 3,
) -> tuple[np.ndarray, str, str]:
    """The quality design for ``samples`` (in that order): chip dummies, log mean intensity and detection fraction
    from ``qc`` (``read_qc``), the measured weight ``share`` and ``extra`` columns (cell fractions). Returns the
    columns, their description and ``known`` or ``unknown`` intensity. Chip is dropped when it would leave fewer than
    ``MIN_RESID_DF`` residual degrees of freedom beside ``base_columns`` (intercept, age, mean)."""
    n = len(samples)
    cols: list[np.ndarray] = []
    names: list[str] = []
    notes: list[str] = []
    intensity = "unknown"
    chip = np.zeros((n, 0))
    if qc is not None:
        q = pl.DataFrame({"sample": samples}).join(qc, on="sample", how="left", validate="1:1")
        if q["sample"].to_list() != samples:
            raise ValueError("clock_disagreement: QC rows out of order")
        if "mean_intensity" in q.columns:
            v = q["mean_intensity"].cast(pl.Float64, strict=False).to_numpy()
            if np.isfinite(v).all() and (v > 0).all():
                cols.append(np.log(v)[:, None])
                names.append("log mean intensity")
                intensity = "known"
        detection = next((c for c in ("frac_dt_species", "frac_dt") if c in q.columns), None)
        if detection is not None:
            v = q[detection].cast(pl.Float64, strict=False).to_numpy()
            if np.isfinite(v).all() and np.std(v) > 0:
                cols.append(v[:, None])
                names.append("detection fraction")
        if "chip" in q.columns and q["chip"].null_count() == 0:
            chip = dummies(q["chip"].cast(pl.String).to_numpy())
    else:
        notes.append("no QC table")
    if share is not None and np.isfinite(share).all() and np.std(share) > 0:
        cols.append(share[:, None])
        names.append("measured weight share")
    for label, v in (extra or {}).items():
        v = np.asarray(v, dtype=float)
        if np.isfinite(v).all() and np.std(v) > 0:
            cols.append(v[:, None])
            names.append(label)
    other = np.hstack(cols) if cols else np.zeros((n, 0))
    if chip.shape[1]:
        full = np.column_stack([np.ones(n), other, chip])
        if n - np.linalg.matrix_rank(full) - (base_columns - 1) >= MIN_RESID_DF:
            other = np.column_stack([other, chip])
            names.insert(0, f"chip ({chip.shape[1] + 1} chips)")
        else:
            notes.append("chip dropped: too few arrays per chip")
    if intensity == "unknown":
        notes.append("intensity unknown")
    return other, ", ".join(names + notes) or "none", intensity


def drift_rows(
    series: str,
    group: str,
    people: pl.DataFrame,
    b,
    f,
    preds,
    clocks,
    age,
    panel,
    min_age_r,
    *,
    quality: tuple[np.ndarray, str, str] | None = None,
    sensitivity: bool = False,
    independent: bool = True,
    provenance: str = "",
) -> tuple:
    """(row of the drift table, per-person frame) for one group and panel, or (None, None) when too few clocks."""
    dis, mean, use = disagreement(preds, age, min_age_r=min_age_r)
    if not np.isfinite(dis).all():
        return None, None
    used = [c for c, u in zip(clocks, use, strict=True) if u]
    b, f = np.asarray(b, dtype=float), np.asarray(f, dtype=float)

    def r(x, y):
        return float(np.corrcoef(x, y)[0, 1]) if np.std(y) > 0 and np.std(x) > 0 else float("nan")

    q_cols, q_terms, intensity = quality if quality is not None else (np.zeros((age.size, 0)), "none", "unknown")
    nan = float("nan")
    quality_r = {"r_dis_B_quality": nan, "partial_r_dis_B_quality": nan, "partial_r_dis_F_quality": nan}
    dis_q = np.full(age.size, np.nan)
    if q_cols.shape[1]:
        dis_q, mean_q, _ = disagreement(preds, age, min_age_r=min_age_r, covariates=q_cols)
        if np.isfinite(dis_q).all():
            held = np.column_stack([age, q_cols])
            quality_r = {
                "r_dis_B_quality": partial_r(dis_q, b, held),
                "partial_r_dis_B_quality": partial_r(dis_q, b, np.column_stack([held, mean_q])),
                "partial_r_dis_F_quality": partial_r(dis_q, f, np.column_stack([held, mean_q])),
            }
    row = {
        "series": series,
        "group": group,
        "panel": panel,
        "sensitivity": sensitivity,
        "independent": independent and not sensitivity,
        "arrays": int(age.size),
        "clocks": len(used),
        "r_dis_B": r(dis, b),
        "r_dis_F": r(dis, f),
        "r_mean_B": r(mean, b),
        "r_mean_F": r(mean, f),
        "partial_r_dis_B": partial_r(dis, b, mean),
        "partial_r_dis_F": partial_r(dis, f, mean),
        "partial_r_dis_B_age": partial_r(dis, b, np.column_stack([age, mean])),
        "partial_r_dis_F_age": partial_r(dis, f, np.column_stack([age, mean])),
        **quality_r,
        "quality_terms": q_terms,
        "intensity": intensity,
        "p_dis_B": float(stats.pearsonr(dis, b).pvalue) if np.std(b) > 0 else nan,
        "p_dis_F": float(stats.pearsonr(dis, f).pvalue) if np.std(f) > 0 else nan,
        "clock_names": ", ".join(used),
        "panel_offered": ", ".join(clocks),
        "provenance": provenance,
    }
    frame = people.select(pl.col("sample").cast(pl.String)).with_columns(
        series=pl.lit(series),
        group=pl.lit(group),
        panel=pl.lit(panel),
        sensitivity=pl.lit(sensitivity),
        age=pl.Series(age),
        B=pl.Series(b),
        F=pl.Series(f),
        disagreement=pl.Series(dis),
        mean_acceleration=pl.Series(mean),
        disagreement_quality=pl.Series(dis_q),
        clocks=pl.lit(len(used)),
    )
    return row, frame


def read_clock_scores(path: Path) -> tuple[pl.DataFrame, dict[str, str]]:
    """The long clock table with string sample ids, checked one row per sample and clock, and clock -> output."""
    long = pl.read_csv(path, schema_overrides={"sample": pl.String})
    if long.select("sample", "clock").is_duplicated().any():
        raise ValueError(f"clock_disagreement: {path} has repeated sample and clock rows")
    outputs = dict(long.select("clock", "output").unique(subset="clock", keep="first", maintain_order=True).iter_rows())
    return long, outputs


def aligned(samples: list[str], long: pl.DataFrame, clocks: list[str], value: str = "prediction") -> np.ndarray:
    """People x clocks of ``value`` in the order of ``samples`` (one-to-one join, order checked)."""
    wide = long.filter(pl.col("clock").is_in(clocks)).pivot(on="clock", index="sample", values=value)
    for c in clocks:
        if c not in wide.columns:
            wide = wide.with_columns(pl.lit(None, dtype=pl.Float64).alias(c))
    frame = pl.DataFrame({"sample": samples}).join(wide, on="sample", how="left", validate="1:1")
    if frame["sample"].to_list() != samples:
        raise ValueError("clock_disagreement: clock rows out of order")
    return frame.select(clocks).to_numpy().astype(float)


def mean_share(samples: list[str], long: pl.DataFrame, clocks: list[str]) -> np.ndarray | None:
    if "measured_weight_share" not in long.columns or not clocks:
        return None
    return np.nanmean(aligned(samples, long, clocks, "measured_weight_share"), axis=1)


def series_groups(ctx: Context) -> tuple[list[dict], list[pl.DataFrame]]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    rows, frames = [], []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv"
        if not (betas_path.exists() and clocks_path.exists()):
            continue
        long, outputs = read_clock_scores(clocks_path)
        available = (set(pl.read_parquet_schema(betas_path)) - {"cpg"}) & set(long["sample"].to_list())
        chosen = selection(ctx.paths, name, spec, available)
        if chosen is None:
            continue
        groups, dropped = chosen
        sensitivity = bool(SENSITIVITY.search(name))
        betas = pl.read_parquet(betas_path).filter(
            pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(sorted(dropped))
        )
        qc = reprocessed.read_qc(ctx.paths, gse).with_columns(pl.col("sample").cast(pl.String))
        provenance = f"SeSAMe rebuild ({gse}), QC passed" + (
            f"; {len(dropped)} listed CpGs left out of B and F" if dropped else ""
        )
        for label, s in groups.items():
            samples = s["sample"].cast(pl.String).to_list()
            if len(samples) < MIN_ARRAYS:
                continue
            panels = {p: age_clocks(outputs, label, p) for p in PANELS}
            if max(len(c) for c in panels.values()) < MIN_CLOCKS:
                continue
            x = betas.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(betas["cpg"].to_list(), ok, strict=True) if k]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                continue
            print(f"clock_disagreement: {name} / {label}: {len(samples)} arrays", file=sys.stderr, flush=True)
            per = loo_scores(x[:, ok], _islands(kept, probes), sets)
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            age = s["age"].to_numpy().astype(float)
            for panel, clocks in panels.items():
                if len(clocks) < MIN_CLOCKS:
                    continue
                row, frame = drift_rows(
                    name,
                    label,
                    s,
                    b,
                    f,
                    aligned(samples, long, clocks),
                    clocks,
                    age,
                    panel,
                    None if panel == PANELS[2] else MIN_AGE_R,
                    quality=quality_terms(samples, qc, share=mean_share(samples, long, clocks)),
                    sensitivity=sensitivity,
                    provenance=provenance,
                )
                if row is not None:
                    rows.append(row)
                    frames.append(frame)
    return rows, frames


def epic_italy(ctx: Context) -> tuple[pl.DataFrame | None, list[dict], list[pl.DataFrame], list[dict]]:
    """EPIC-Italy people with clocks, the drift rows (everyone; women), the per-person frames and the risk rows."""
    src = ctx.paths.output / "incident_cancer_rsesame" / "per_person.parquet"
    clocks_path = ctx.paths.root / CLOCK_SCORES / f"{EPIC_ITALY}_clocks.csv"
    if not (src.exists() and clocks_path.exists() and reprocessed.available(ctx.paths, EPIC_ITALY)):
        print("clock_disagreement: EPIC-Italy needs incident_cancer_rsesame per_person and its clocks", file=sys.stderr)
        return None, [], [], []
    qc = reprocessed.read_qc(ctx.paths, EPIC_ITALY).with_columns(pl.col("sample").cast(pl.String))
    passed = qc.filter(pl.col("passed"))["sample"].to_list()
    long, outputs = read_clock_scores(clocks_path)
    people = (
        pl.read_parquet(src)
        .with_columns(pl.col("sample").cast(pl.String))
        .select("sample", "female", "age", "icd", "years_to_diagnosis", "chip", "mean_beta", *FRACTIONS, B, F)
        .filter(pl.col("sample").is_in(passed) & pl.col("sample").is_in(long["sample"].unique().to_list()))
        .drop_nulls(["age", B, F])
        .sort("sample")
    )
    if people["sample"].is_duplicated().any():
        raise ValueError("clock_disagreement: EPIC-Italy has repeated samples")
    qc_terms = qc.select("sample", *(c for c in ("chip", "mean_intensity", "frac_dt") if c in qc.columns))
    people = people.join(
        qc_terms.select("sample", "mean_intensity", "frac_dt"), on="sample", how="left", validate="1:1"
    )
    rows, frames = [], []
    for label, part, independent in (
        ("blood, everyone", people, True),
        ("blood, women", people.filter(pl.col("female")), False),
    ):
        samples = part["sample"].to_list()
        for panel in PANELS:
            clocks = age_clocks(outputs, "blood", panel)
            if len(clocks) < MIN_CLOCKS:
                continue
            quality = quality_terms(
                samples,
                qc_terms,
                share=mean_share(samples, long, clocks),
                extra={c: part[c].to_numpy() for c in FRACTIONS},
            )
            row, frame = drift_rows(
                "EPIC-Italy (GSE51032)",
                label,
                part,
                part[B].to_numpy(),
                part[F].to_numpy(),
                aligned(samples, long, clocks),
                clocks,
                part["age"].to_numpy().astype(float),
                panel,
                None if panel == PANELS[2] else MIN_AGE_R,
                quality=quality,
                independent=independent,
                provenance=EPIC_PROVENANCE,
            )
            if row is not None:
                rows.append(row)
                frames.append(frame)
    women = people.filter(pl.col("female"))
    clocks = age_clocks(outputs, "blood")
    dis, mean, use = disagreement(aligned(women["sample"].to_list(), long, clocks), women["age"].to_numpy())
    women = women.with_columns(disagreement=pl.Series(dis), mean_acceleration=pl.Series(mean))
    used = [c for c, u in zip(clocks, use, strict=True) if u]
    print(f"clock_disagreement: EPIC-Italy women {women.height}, clocks {', '.join(used)}", file=sys.stderr)
    return women, rows, frames, breast_rows(women, FRACTIONS, used)


def _diagnose(beta, se, info, grad, n) -> tuple[np.ndarray, np.ndarray, str]:
    nan = np.full(beta.size, np.nan)
    if not (np.isfinite(beta).all() and np.isfinite(se).all() and np.isfinite(info).all()):
        return nan, nan, "not finite"
    if np.abs(beta).max() > MAX_COEF:
        return nan, nan, f"separation: a coefficient beyond ±{MAX_COEF:g} on the log-odds scale"
    condition = float(np.linalg.cond(info))
    if not condition < MAX_CONDITION:
        return nan, nan, f"singular Hessian (condition number {condition:.2g})"
    if np.abs(grad).max() > 1e-4 * max(1, n):
        return nan, nan, "did not converge"
    return beta, se, "ok"


def fit_logit(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:
    """``logit`` with diagnostics: coefficients and SEs, or NaN with the reason."""
    nan = np.full(x.shape[1], np.nan)
    if y.min() == y.max():
        return nan, nan, "one outcome only"
    if np.linalg.matrix_rank(x) < x.shape[1]:
        return nan, nan, "design not of full rank"
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        beta, se = logit(y, x)
        p = 1.0 / (1.0 + np.exp(-(x @ beta)))
        info = x.T @ (x * (p * (1 - p))[:, None])
        grad = x.T @ (y - p)
    return _diagnose(beta, se, info, grad, y.size)


def fit_clogit(x: np.ndarray, sets: list[list[int]]) -> tuple[np.ndarray, np.ndarray, str]:
    """``clogit`` with diagnostics: coefficients and SEs, or NaN with the reason."""
    nan = np.full(x.shape[1], np.nan)
    if len(sets) < 2:
        return nan, nan, "fewer than 2 matched sets"
    centred = np.vstack([x[s] - x[s].mean(axis=0) for s in sets])
    if np.linalg.matrix_rank(centred) < x.shape[1]:
        return nan, nan, "design not of full rank within sets"
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        beta, se = clogit(x, sets)
        _, grad, info = _score(x, sets, beta)
    return _diagnose(beta, se, info, grad, len(sets))


def _or_row(cohort, outcome, lag, model, terms, beta, se, status, n_clocks, **extra) -> list[dict]:
    out = []
    for i, term in enumerate(terms):
        ok = bool(np.isfinite(beta[i]) and np.isfinite(se[i]) and se[i] > 0)
        out.append(
            {
                "cohort": cohort,
                "outcome": outcome,
                "lag": lag,
                "model": model,
                "term": term,
                **extra,
                "clocks": n_clocks,
                "odds_ratio": float(np.exp(beta[i])) if ok else float("nan"),
                "low": float(np.exp(beta[i] - 1.96 * se[i])) if ok else float("nan"),
                "high": float(np.exp(beta[i] + 1.96 * se[i])) if ok else float("nan"),
                "p": float(2 * stats.norm.sf(abs(beta[i] / se[i]))) if ok else float("nan"),
                "fit": status,
            }
        )
    return out


def breast_rows(women: pl.DataFrame, fractions: list[str], clocks: list[str]) -> list[dict]:
    breast = pl.col("icd").fill_null("").str.starts_with("C50")
    controls = women.filter(pl.col("icd").is_null())
    rows: list[dict] = []
    terms = ("disagreement", "mean acceleration")
    n_clocks = len(clocks)
    base = {
        "estimand": EPIC_ESTIMAND,
        "provenance": EPIC_PROVENANCE,
        "panel": ", ".join(clocks),
        "followup": "controls: no cancer diagnosis recorded in the cohort data",
    }

    def add(lag, model, used_terms, beta, se, status, **sizes):
        extra = {"chips": None, "controls_on_case_chip": None} | sizes | base
        cohort, outcome = "EPIC-Italy women", "breast cancer"
        rows.extend(_or_row(cohort, outcome, lag, model, used_terms, beta, se, status, n_clocks, **extra))

    def covariates(pool: pl.DataFrame) -> np.ndarray:
        return np.column_stack([zs(pool[c].to_numpy()) for c in ("mean_beta", *fractions)])

    for lag, cond in ((f"> {LAG:g} y", pl.col("years_to_diagnosis") > LAG), ("any", pl.lit(True))):
        cases = women.filter(breast & cond)
        eligible = {"eligible_cases": cases.height, "eligible_controls": controls.height}
        # unconditional models on every eligible woman
        pool = pl.concat([cases, controls], how="vertical")
        y = np.r_[np.ones(cases.height), np.zeros(controls.height)]
        dis, mean = zs(pool["disagreement"].to_numpy()), zs(pool["mean_acceleration"].to_numpy())
        chip = pool["chip"].cast(pl.String).to_numpy()
        dis_c, mean_c = zs(dis - _group_means(dis, chip)), zs(mean - _group_means(mean, chip))
        age = zs(pool["age"].to_numpy())
        cov = covariates(pool)
        quality = np.column_stack(
            [zs(np.log(pool["mean_intensity"].cast(pl.Float64).to_numpy())), zs(pool["frac_dt"].to_numpy())]
        )
        one = np.ones(pool.height)
        sizes = eligible | {"cases": cases.height, "controls": controls.height}
        designs = {
            "age": np.column_stack([one, dis, mean, age]),
            "covariates": np.column_stack([one, dis, mean, age, cov]),
            "covariates + quality": np.column_stack([one, dis, mean, age, cov, quality]),
            "covariates, chip-centred": np.column_stack([one, dis_c, mean_c, age, cov]),
        }
        for model, design in designs.items():
            beta, se, status = fit_logit(y, design)
            add(lag, model, terms, beta[1:], se[1:], status, **sizes)
        beta, se, status = fit_logit(y, np.column_stack([one, dis, age, cov]))
        add(lag, "covariates, disagreement alone", terms[:1], beta[1:], se[1:], status, **sizes)

        # chip intercepts, on the chips holding at least one case and one control (the others carry no
        # within-chip contrast and would separate)
        both = sorted(set(cases["chip"].cast(pl.String)) & set(controls["chip"].cast(pl.String)))
        kc = cases.filter(pl.col("chip").cast(pl.String).is_in(both))
        kn = controls.filter(pl.col("chip").cast(pl.String).is_in(both))
        pool = pl.concat([kc, kn], how="vertical")
        y = np.r_[np.ones(kc.height), np.zeros(kn.height)]
        design = np.column_stack(
            [
                np.ones(pool.height),
                zs(pool["disagreement"].to_numpy()),
                zs(pool["mean_acceleration"].to_numpy()),
                zs(pool["age"].to_numpy()),
                covariates(pool),
                dummies(pool["chip"].cast(pl.String).to_numpy()),
            ]
        )
        beta, se, status = fit_logit(y, design)
        sizes = eligible | {"cases": kc.height, "controls": kn.height, "chips": len(both)}
        add(lag, "covariates, chip effects", terms, beta[1:3], se[1:3], status, **sizes)

        # conditional models on age-matched sets: same chip first, then same chip only
        pool = pl.concat([cases, controls], how="vertical")
        x = np.column_stack(
            [
                zs(pool["disagreement"].to_numpy()),
                zs(pool["mean_acceleration"].to_numpy()),
                zs(pool["age"].to_numpy()),
                covariates(pool),
            ]
        )
        chips = pool["chip"].cast(pl.String).to_list()
        for model, exact in (("matched sets, covariates", False), ("matched sets, same chip only, covariates", True)):
            sets = match(cases, controls, same_chip_only=exact)
            pairs = [chips[s[0]] == chips[j] for s in sets for j in s[1:]]
            beta, se, status = fit_clogit(x, sets)
            sizes = eligible | {
                "cases": len(sets),
                "controls": sum(len(s) - 1 for s in sets),
                "controls_on_case_chip": float(np.mean(pairs)) if pairs else float("nan"),
            }
            add(lag, model, terms, beta[:2], se[:2], status, **sizes)
    return rows


def sjlife(ctx: Context) -> tuple[list[dict], list[pl.DataFrame], list[dict]]:
    src = ctx.paths.output / "sjlife_deaths" / "per_person.parquet"
    if not src.exists():
        print("clock_disagreement: no sjlife_deaths per_person; SJLIFE skipped", file=sys.stderr)
        return [], [], []
    provenance = f"GEO deposited betas ({SJLIFE_GSE} processed table, read by sjlife_deaths); no SeSAMe rebuild, no QC"
    if reprocessed.available(ctx.paths, SJLIFE_GSE):
        provenance += "; a rebuild exists but sjlife_deaths has not been rerun on it"
        print(f"clock_disagreement: {SJLIFE_GSE} has a rebuild; sjlife_deaths still reads the deposited betas")
    people = pl.read_parquet(src).with_columns(pl.col("sample").cast(pl.String)).sort("sample")
    clocks = [c for c in SJLIFE_CLOCKS if c in people.columns]
    pheno = clocks + [c for c in SJLIFE_PHENOTYPIC if c in people.columns]
    people = people.drop_nulls(["age", B, F, *pheno])
    if people["sample"].is_duplicated().any():
        raise ValueError("clock_disagreement: SJLIFE has repeated samples")
    age = people["age"].to_numpy().astype(float)
    rows, frames = [], []
    no_quality = (np.zeros((people.height, 0)), "none: deposited betas, no QC table", "unknown")
    for panel, use_clocks in ((PANELS[0], clocks), (PANELS[3], pheno)):
        row, frame = drift_rows(
            "SJLIFE (GSE314261)",
            "PBMC, everyone",
            people,
            people[B].to_numpy(),
            people[F].to_numpy(),
            people.select(use_clocks).to_numpy().astype(float),
            use_clocks,
            age,
            panel,
            None,
            quality=no_quality,
            provenance=provenance,
        )
        if row is not None:
            rows.append(row | {"followup": SJLIFE_FOLLOWUP})
            frames.append(frame)
    if not rows:
        return [], [], []
    dis, mean, _ = disagreement(people.select(clocks).to_numpy().astype(float), age, min_age_r=None)
    died = people["age_at_death"].is_not_null().to_numpy()
    years = (people["age_at_death"].fill_null(np.inf) - people["age"]).to_numpy()
    survivor = people["survivor"].to_numpy()
    risk = []
    diagnosis = people["diagnosis"].cast(pl.String).fill_null("NA").to_numpy()
    for lag, eligible in (("any", survivor), (f"> {LAG:g} y", survivor & (~died | (years > LAG)))):
        # a diagnosis with no death among the eligible separates (its intercept runs to minus infinity): its
        # people carry no information on the terms and are left out, as the limit of the fit would leave them
        with_death = set(diagnosis[eligible & died].tolist())
        keep = eligible & np.isin(diagnosis, sorted(with_death))
        n_eligible = (int((eligible & died).sum()), int((eligible & ~died).sum()))
        y = died[keep].astype(float)
        diag = dummies(people["diagnosis"].to_numpy()[keep])
        exposures = [people[e].to_numpy()[keep] for e in EXPOSURES if e in people.columns]
        design = np.column_stack(
            [
                np.ones(keep.sum()),
                zs(dis[keep]),
                zs(mean[keep]),
                age[keep],
                people["female"].to_numpy()[keep].astype(float),
                *exposures,
                diag[:, diag.std(axis=0) > 0],
            ]
        )
        beta, se, status = fit_logit(y, design)
        risk += _or_row(
            "SJLIFE survivors",
            "recorded death",
            lag,
            "age, sex, exposures, diagnosis",
            ("disagreement", "mean acceleration"),
            beta[1:],
            se[1:],
            status,
            len(clocks),
            eligible_cases=n_eligible[0],
            eligible_controls=n_eligible[1],
            cases=int(y.sum()),
            controls=int((1 - y).sum()),
            chips=None,
            controls_on_case_chip=None,
            estimand=SJLIFE_ESTIMAND,
            provenance=provenance,
            panel=", ".join(clocks),
            followup=SJLIFE_FOLLOWUP,
        )
    return rows, frames, risk


def summarise(drift: pl.DataFrame) -> pl.DataFrame:
    """Per panel, over the independent groups (no sensitivity configs, no subsets): medians of the correlations, and
    paired medians on the groups that have the quality-adjusted value."""
    with_q = pl.col("partial_r_dis_B_quality").is_not_nan()
    return (
        drift.filter(pl.col("independent"))
        .group_by("panel", maintain_order=True)
        .agg(
            pl.len().alias("groups"),
            pl.col("r_dis_B").median().alias("median_r_dis_B"),
            pl.col("r_dis_F").median().alias("median_r_dis_F"),
            pl.col("partial_r_dis_B").median().alias("median_partial_r_dis_B"),
            pl.col("partial_r_dis_B_age").median().alias("median_partial_r_dis_B_age"),
            with_q.sum().alias("groups_with_quality"),
            pl.col("partial_r_dis_B_age").filter(with_q).median().alias("median_partial_r_dis_B_age_same_groups"),
            pl.col("partial_r_dis_B_quality").filter(with_q).median().alias("median_partial_r_dis_B_quality"),
            pl.col("partial_r_dis_F_quality").filter(with_q).median().alias("median_partial_r_dis_F_quality"),
            (pl.col("r_dis_B") > 0).mean().alias("share_r_dis_B_positive"),
        )
    )


def run(ctx: Context) -> list[Path]:
    women, e_rows, e_frames, risk = epic_italy(ctx)
    s_rows, s_frames, s_risk = sjlife(ctx)
    rows, frames = series_groups(ctx)
    rows, frames, risk = rows + e_rows + s_rows, frames + e_frames + s_frames, risk + s_risk
    if not rows:
        print("clock_disagreement: nothing to score", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_disagreement")
    drift = pl.DataFrame(rows, infer_schema_length=None)
    write_table(drift, out / "drift")
    write_table(summarise(drift), out / "drift_summary")
    write_table(pl.concat(frames, how="vertical"), out / "per_person", csv_max_rows=10_000)
    risk_table = pl.DataFrame(risk, infer_schema_length=None) if risk else pl.DataFrame()
    if risk:
        write_table(risk_table, out / "risk")
    if women is not None:
        write_table(
            women.select("sample", "age", "chip", "icd", "years_to_diagnosis", "disagreement", "mean_acceleration"),
            out / "epic_italy_women",
        )
    return [draw(drift.filter(pl.col("independent")), risk_table, ctx.figure_path("clock_disagreement"))]


def draw(drift: pl.DataFrame, risk: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    ax = axes[0]
    for i, panel in enumerate(PANELS[:2]):
        d = drift.filter(pl.col("panel") == panel)
        for j, col in enumerate(("r_dis_B", "r_dis_F")):
            v = d[col].to_numpy()
            n = d["clocks"].to_numpy()
            xpos = j * 2.6 + i + np.random.default_rng(j * 10 + i).normal(0, 0.07, v.size)
            ax.scatter(xpos[n >= 3], v[n >= 3], s=14, alpha=0.7, color=plots.PALETTE[i])
            ax.scatter(xpos[n < 3], v[n < 3], s=14, alpha=0.7, facecolors="none", edgecolors=plots.PALETTE[i])
            if v.size:
                ax.plot([j * 2.6 + i - 0.3, j * 2.6 + i + 0.3], [np.median(v)] * 2, color=plots.INK, lw=2)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_xticks(
        [0, 1, 2.6, 3.6],
        ["valid\nclocks", "+ outside\ntissue", "valid\nclocks", "+ outside\ntissue"],
        fontsize=7,
    )
    ax.text(0.5, -0.2, "with drift B", transform=ax.get_xaxis_transform(), ha="center", fontsize=8)
    ax.text(3.1, -0.2, "with shield share F", transform=ax.get_xaxis_transform(), ha="center", fontsize=8)
    ax.set_ylabel("r, disagreement with the drift measure\n(dot: one group; open: two clocks)")
    ax.set_title("Disagreement and drift, per tissue group", loc="left", fontsize=9)
    plots.label_panel(ax, "a")
    ax = axes[1]
    if risk.height:
        pick = risk.filter(
            pl.col("model").is_in(["covariates", "age, sex, exposures, diagnosis"]) & pl.col("odds_ratio").is_not_nan()
        )
        labels = []
        for y, r in enumerate(pick.iter_rows(named=True)):
            color = plots.PALETTE[0] if r["term"] == "disagreement" else plots.PALETTE[1]
            ax.errorbar(
                r["odds_ratio"],
                y,
                xerr=[[r["odds_ratio"] - r["low"]], [r["high"] - r["odds_ratio"]]],
                fmt="o",
                ms=4,
                color=color,
            )
            labels.append(f"{r['cohort']}, {r['outcome']} ({r['lag']}): {r['term']}  [{r['cases']}/{r['controls']}]")
        ax.set_yticks(range(len(labels)), labels, fontsize=6.5)
        ax.invert_yaxis()
        ax.axvline(1.0, color=plots.MUTED, lw=0.8)
        ax.set_xscale("log")
        ax.set_xticks([0.7, 1.0, 1.5, 2.0], ["0.7", "1", "1.5", "2"])
        ax.xaxis.set_minor_formatter(plots.plt.NullFormatter())
    ax.set(xlabel="odds ratio per SD (95 % CI), both terms in one model")
    ax.set_title("Disagreement beside the mean acceleration", loc="left", fontsize=9)
    plots.label_panel(ax, "b")
    fig.tight_layout()
    return plots.save(fig, path)
