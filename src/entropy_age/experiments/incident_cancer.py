"""Does the shield predict cancer years ahead? EPIC-Italy (roadmap RM28, second cohort)

The first open cohort with follow-up (docs/LONGITUDINAL_SURVEY.md): EPIC-Italy (GSE51032, Polidoro et al.), 845
people whose blood was drawn at enrolment (450k array), 421 of whom were later diagnosed with cancer (ICD-10 code,
mostly breast C50 and colorectal C18-C20), with the time from blood draw to diagnosis in years.

Per person: ``F``, ``B`` and ``S`` of the gene sets against everyone else (leave-one-out, ``health_states``) and
published clocks' linear predictors, each standardised. Each enters a logistic regression of later cancer with age
and sex (``logit``, Newton steps), reported as an odds ratio per SD:

- all cancers; breast; colorectal;
- cases diagnosed more than 2 years after the draw only, so an undiagnosed tumour present at the draw cannot drive
  the result;
- each also with every person's mean methylation over all array CpGs (``mean_betas``; genome-wide hypomethylation
  before breast cancer was reported in this cohort, van Veldhoven 2015), 11 of 12 blood-cell fractions (EpiDISH
  cent12CT, as in ``stochastic_component``) and the array chip held (``adjusted``).

``entropy-age figures cancer`` (as a Slurm job; the 3 GB series matrix holds the betas).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.shield_clocks import linear_score
from entropy_age.experiments.stochastic_component import EPIDISH, deconvolve
from entropy_age.experiments.syndromes import matrix_path, read_betas
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE51032"
CLOCKS_USED = ("Horvath 2013", "Hannum", "PhenoAge", "DunedinPACE", "Zhang 2019", "StocH", "StocP", "StocZ")
GROUPS = ("core_essential", "housekeeping")
LAG = 2.0  # years between draw and diagnosis for the lagged analysis
OUTCOMES = {"any cancer": None, "breast": ("C50",), "colorectal": ("C18", "C19", "C20")}


def logit(y: np.ndarray, x: np.ndarray, *, steps: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Logistic regression coefficients and SEs by Newton steps; ``x`` carries its own intercept column."""
    beta = np.zeros(x.shape[1])
    for _ in range(steps):
        p = 1.0 / (1.0 + np.exp(-(x @ beta)))
        w = p * (1 - p)
        info = x.T @ (x * w[:, None])
        step = np.linalg.solve(info, x.T @ (y - p))
        beta += step
        if np.abs(step).max() < 1e-10:
            break
    p = 1.0 / (1.0 + np.exp(-(x @ beta)))
    info = x.T @ (x * (p * (1 - p))[:, None])
    return beta, np.sqrt(np.diag(np.linalg.inv(info)))


def read_people(matrix: Path) -> pl.DataFrame:
    """``sample``, ``female``, ``age``, ``icd`` (null for no cancer) and ``years_to_diagnosis``."""
    chars = _characteristics(matrix)
    icd = next(c for c in chars.columns if c.startswith("cancer type"))
    ttd = next(c for c in chars.columns if c.startswith("time to diagnosis"))
    return chars.select(
        "sample",
        (pl.col("gender") == "F").alias("female"),
        pl.col("age").cast(pl.Float64, strict=False),
        pl.col(icd).replace("", None).alias("icd"),
        pl.col(ttd).cast(pl.Float64, strict=False).alias("years_to_diagnosis"),
    )


def mean_betas(matrix: Path, samples: list[str]) -> np.ndarray:
    """Each sample's mean beta over every CpG row of the series matrix's table (missing values skipped)."""
    total, count, columns = None, None, None
    with gzip.open(matrix, "rt") as handle:
        inside = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not inside:
                continue
            cells = line.rstrip("\n").split("\t")
            if columns is None:
                header = [c.strip('"') for c in cells[1:]]
                columns = np.array([header.index(x) for x in samples])
                total, count = np.zeros(len(samples)), np.zeros(len(samples))
                continue
            if not cells[0].strip('"').startswith("cg"):
                continue
            row = np.array([float(c) if c not in ("", "null", "NA", '""') else np.nan for c in cells[1:]])[columns]
            ok = ~np.isnan(row)
            total[ok] += row[ok]
            count[ok] += 1
    return total / np.maximum(count, 1)


def chips(matrix: Path) -> dict[str, str]:
    """GEO accession -> array chip, from titles like ``6969568099_R02C02``."""
    acc = titles = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_geo_accession"):
                acc = cells
            if acc and titles:
                break
    return {a: t.split("_")[0] for a, t in zip(acc, titles, strict=True)}


def outcome_rows(people: pl.DataFrame, measures: dict[str, np.ndarray], extra: np.ndarray | None = None) -> list[dict]:
    """Odds ratios per SD; with ``extra`` (people x covariates) each is also fitted with them held."""
    rows = []
    age = people["age"].to_numpy()
    female = people["female"].to_numpy().astype(float)
    icd = people["icd"].fill_null("").to_numpy()
    ttd = people["years_to_diagnosis"].fill_null(np.inf).to_numpy()
    control = icd == ""
    for outcome, codes in OUTCOMES.items():
        case = ~control if codes is None else np.array([c.startswith(codes) for c in icd])
        for lagged in (False, True):
            keep = control | (case & ((ttd > LAG) if lagged else True))
            y = case[keep].astype(float)
            if y.sum() < 20:
                continue
            for name, v in measures.items():
                z = (v - v.mean()) / v.std(ddof=1)
                for adjusted in (False, True) if extra is not None else (False,):
                    covs = [np.ones(keep.sum()), z[keep], age[keep]]
                    if np.unique(female[keep]).size > 1:
                        covs.append(female[keep])
                    if adjusted:
                        e = extra[keep]
                        covs.append(e[:, e.std(axis=0) > 0])
                    rows.append(_odds_row(outcome, lagged, adjusted, name, y, np.column_stack(covs)))
    return rows


def _odds_row(outcome, lagged, adjusted, name, y, design) -> dict:
    beta, se = logit(y, design)
    return {
        "outcome": outcome,
        "lag_over_2y": lagged,
        "adjusted": adjusted,
        "measure": name,
        "cases": int(y.sum()),
        "controls": int((1 - y).sum()),
        "odds_ratio": float(np.exp(beta[1])),
        "low": float(np.exp(beta[1] - 1.96 * se[1])),
        "high": float(np.exp(beta[1] + 1.96 * se[1])),
        "p": float(2 * stats.norm.sf(abs(beta[1] / se[1]))),
    }


def run(ctx: Context) -> list[Path]:
    matrix = matrix_path(ctx.paths, GSE)
    if not matrix.exists():
        print(f"cancer: {matrix} missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("incident_cancer")
    people = read_people(matrix).drop_nulls("age")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    lay = ClockLayout(ctx.paths)
    clocks = {c: read_clock(lay.file(c)) for c in CLOCKS_USED if lay.file(c).exists()}
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    import rdata

    ref_cpgs = [
        str(c) for c in next(iter(rdata.read_rda(str(ctx.paths.input / EPIDISH)).values())).coords["dim_0"].values
    ]
    wanted = sorted(set(features) | set(clock_cpgs) | set(ref_cpgs))
    x_all = read_betas(matrix, wanted, people["sample"].to_list())
    col = {f: j for j, f in enumerate(wanted)}
    x = x_all[:, [col[f] for f in features]]
    ok = np.isnan(x).mean(axis=0) == 0
    kept = [f for f, k in zip(features, ok, strict=True) if k]
    print(
        f"cancer: {people.height} people, {people['icd'].is_not_null().sum()} cases, {len(kept)} CpGs",
        file=sys.stderr,
        flush=True,
    )
    scores = loo_scores(x[:, ok], _islands(kept, probes), _sets(kept, membership))
    measures = {f"{m} {g}": scores[f"{m}_{g}"] for g in GROUPS for m in ("F", "B", "S") if f"{m}_{g}" in scores}
    xc = x_all[:, [col[f] for f in clock_cpgs]]
    age = people["age"].to_numpy()
    for name, weights in clocks.items():
        score, used = linear_score(xc, clock_cpgs, weights)
        if used and np.isfinite(score).all():
            design = np.column_stack([np.ones(age.size), age])
            coef, *_ = np.linalg.lstsq(design, score, rcond=None)
            measures[name] = score - design @ coef  # acceleration
    extra = _adjustments(ctx, matrix, people, x_all, wanted)
    table = pl.DataFrame(outcome_rows(people, measures, extra))
    write_table(table, out / "effects")
    write_table(people.with_columns(**{k: pl.Series(v) for k, v in measures.items()}), out / "per_person")
    return [_figure(ctx, table)]


def _adjustments(ctx, matrix: Path, people: pl.DataFrame, x_all: np.ndarray, wanted: list[str]) -> np.ndarray:
    """People x [mean beta, 11 cell fractions, chip dummies]."""
    import rdata

    samples = people["sample"].to_list()
    mean = mean_betas(matrix, samples)
    ref = next(iter(rdata.read_rda(str(ctx.paths.input / EPIDISH)).values()))
    ref_cpgs = [str(c) for c in ref.coords["dim_0"].values]
    col = {f: j for j, f in enumerate(wanted)}
    keep = [i for i, c in enumerate(ref_cpgs) if c in col]
    xr = x_all[:, [col[ref_cpgs[i]] for i in keep]]
    xr = np.where(np.isnan(xr), np.nanmean(xr, axis=0), xr)
    fractions = deconvolve(xr, ref.values[keep])
    chip = chips(matrix)
    print(
        f"cancer: mean beta {mean.mean():.3f}; {len(keep)} reference CpGs; neutrophils {fractions[:, -2].mean():.2f}",
        file=sys.stderr,
        flush=True,
    )
    return np.column_stack([mean, fractions[:, :-1], dummies(np.array([chip[x] for x in samples]))])


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    outcomes = table["outcome"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(len(outcomes), ncols=len(outcomes), size=4.0)
    measures = table["measure"].unique(maintain_order=True).to_list()
    for i, outcome in enumerate(outcomes):
        ax = axes[i]
        for y, m in enumerate(measures):
            for j, lagged in enumerate((False, True)):
                r = table.filter(
                    (pl.col("outcome") == outcome)
                    & (pl.col("measure") == m)
                    & (pl.col("lag_over_2y") == lagged)
                    & (pl.col("adjusted") == lagged)
                )
                if r.height:
                    ax.errorbar(
                        r["odds_ratio"][0],
                        y + 0.25 * j,
                        xerr=[[r["odds_ratio"][0] - r["low"][0]], [r["high"][0] - r["odds_ratio"][0]]],
                        fmt="os"[j],
                        ms=3.5,
                        color=plots.PALETTE[j],
                    )
        ax.axvline(1.0, color=plots.MUTED, lw=0.8)
        ax.set_xscale("log")
        ax.set_yticks(range(len(measures)), measures if i == 0 else [""] * len(measures), fontsize=7.5)
        ax.invert_yaxis()
        n = table.filter(pl.col("outcome") == outcome)["cases"].max()
        ax.set(xlabel="odds ratio per SD (95 % CI)", title=f"{outcome} ({n} cases)")
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "EPIC-Italy: circles all cases; squares diagnosed over 2 years later, mean beta, cell mix and chip held",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("incident_cancer"))
