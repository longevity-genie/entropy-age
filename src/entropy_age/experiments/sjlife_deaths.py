"""Do drift and the shield in blood mark later death? St. Jude Lifetime cohort (roadmap RM28, third cohort)

SJLIFE (GSE314261): blood (PBMC) of about 4,700 childhood-cancer survivors and 340 community controls on the EPIC
array, with age, sex, diagnosis, 19 treatment exposures and, for those who have died, the age at death. GEO gives
no censoring date for the living, so this is a logistic model of having died by the data freeze, with age at the
draw held (older survivors have had longer to die), sex, survivor or control, and the main exposures (chest and
brain radiation, anthracyclines, alkylating agents). Each person's first draw is used.

The breast-cancer result in EPIC-Italy (C12, drift B at unprotected promoters) is the reason to look here first:
B, the shield's share F and the signed shift S per person, against one reference fitted on everyone (with
thousands of people a person's own weight in it is negligible), and published clocks' accelerations, each per SD.
Deaths more than 2 years after the draw are tested separately. No control has died, so models use survivors only.
Each is also fitted with mean methylation of the loaded CpGs and 11 blood-cell fractions held and the measure
centred within its array chip (as in ``incident_cancer``).

``entropy-age figures sjlife`` (as a Slurm job; the 28 GB processed table).
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
from entropy_age.experiments.incident_cancer import _group_means, logit
from entropy_age.experiments.per_person import _islands, _sets, score_against
from entropy_age.experiments.shield_clocks import linear_score
from entropy_age.experiments.stochastic_component import EPIDISH, deconvolve
from entropy_age.experiments.syndromes import matrix_path
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import fit_reference

GSE = "GSE314261"
TABLE = "GSE314261_SJLIFE_IlluminaEPIC_proccessed.txt.gz"
CLOCKS_USED = ("Horvath 2013", "Hannum", "PhenoAge", "DunedinPACE", "Zhang 2019")
GROUPS = ("core_essential", "housekeeping")
EXPOSURES = (
    "treatment_potential_chest_5",
    "treatment_potential_brain_5",
    "treatment_anthracyclines_5",
    "treatment_alkylating_classic_5",
)
LAG = 2.0


def read_people(matrix: Path) -> pl.DataFrame:
    """One row per sample: accession, barcode, title, age, sex, group, diagnosis, exposures, age at death."""
    chars = _characteristics(matrix)
    barcodes = titles = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_description") and barcodes is None:
                barcodes = cells
            elif line.startswith("!Sample_title"):
                titles = cells
            if barcodes is not None and titles is not None:
                break
    acc = chars["sample"].to_list()
    death = "age at death" if "age at death" in chars.columns else None
    return chars.with_columns(barcode=pl.Series(barcodes[: len(acc)]), title=pl.Series(titles[: len(acc)])).select(
        "sample",
        "barcode",
        "title",
        pl.col("age").cast(pl.Float64, strict=False),
        (pl.col("Sex") == "female").alias("female"),
        (pl.col("group") == "cancer survivor").alias("survivor"),
        pl.col("diagnostic group").fill_null("control").alias("diagnosis"),
        *(pl.col(e).cast(pl.Float64, strict=False).fill_null(0.0) for e in EXPOSURES if e in chars.columns),
        (pl.col(death).cast(pl.Float64, strict=False) if death else pl.lit(None, dtype=pl.Float64)).alias(
            "age_at_death"
        ),
    )


def first_draws(people: pl.DataFrame) -> pl.DataFrame:
    """Each person's earliest draw: titles of repeated draws share the part before ``_``."""
    return (
        people.with_columns(person=pl.col("title").str.split("_").list.first())
        .sort("age")
        .unique(subset="person", keep="first", maintain_order=True)
    )


def stream_header(table: Path) -> list[str]:
    with gzip.open(table, "rt") as handle:
        return handle.readline().rstrip("\n").split(",")


def stream_rows(table: Path, wanted: set[str], columns: list[str]) -> tuple[list[str], np.ndarray]:
    """The rows of ``wanted`` CpGs for ``columns`` (samples x CpGs), reading the gzip line by line."""
    with gzip.open(table, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
        index = np.array([header.index(c) for c in columns])
        names, rows = [], []
        for line in handle:
            cpg, _, rest = line.partition(",")
            if cpg.strip('"') not in wanted:
                continue
            cells = np.array(rest.rstrip("\n").split(","), dtype=object)
            values = np.array([float(v) if v not in ("", "NA") else np.nan for v in cells[index - 1]])
            names.append(cpg.strip('"'))
            rows.append(values)
    return names, np.array(rows).T


def reference_scores(x: np.ndarray, island: np.ndarray, sets: dict) -> dict[str, np.ndarray]:
    """Everyone's scores against one reference fitted on everyone (for cohorts of thousands)."""
    member = np.any([s for s, _ in sets.values()], axis=0)
    return score_against(x, fit_reference(x, island, ~member), sets)


def death_rows(people: pl.DataFrame, measures: dict[str, np.ndarray], extra: np.ndarray | None = None) -> list[dict]:
    """Odds ratio of having died per SD, among survivors (no control has died); with ``extra`` (people x
    covariates: mean beta, cell fractions) also adjusted, the measure then centred within its array chip."""
    age = people["age"].to_numpy()
    died = people["age_at_death"].is_not_null().to_numpy()
    years = (people["age_at_death"].fill_null(np.inf) - people["age"]).to_numpy()
    base = [age, people["female"].to_numpy().astype(float)]
    exposures = [people[e].to_numpy() for e in EXPOSURES if e in people.columns]
    survivor = people["survivor"].to_numpy()
    chip = people["barcode"].str.split("_").list.first().to_numpy() if "barcode" in people.columns else None
    rows = []
    for lagged in (False, True):
        keep = survivor & (~died | (years > LAG)) if lagged else survivor
        y = died[keep].astype(float)
        if y.sum() < 20:
            continue
        diag = dummies(people["diagnosis"].to_numpy()[keep])
        for name, v in measures.items():
            z = (v - v.mean()) / v.std(ddof=1)
            for adjusted in (False, True) if extra is not None else (False,):
                zz = z
                if adjusted and chip is not None:
                    zz = z - _group_means(z, chip)
                    zz = zz / zz.std(ddof=1)
                cols = [np.ones(keep.sum()), zz[keep], *(b[keep] for b in base), *(e[keep] for e in exposures)]
                if adjusted:
                    e = extra[keep]
                    cols.append(e[:, e.std(axis=0) > 0])
                design = np.column_stack([*cols, diag[:, diag.std(axis=0) > 0]])
                beta, se = logit(y, design)
                rows.append(
                    {
                        "subset": "survivors",
                        "lag_over_2y": lagged,
                        "adjusted": adjusted,
                        "measure": name,
                        "deaths": int(y.sum()),
                        "alive": int((1 - y).sum()),
                        "odds_ratio": float(np.exp(beta[1])),
                        "low": float(np.exp(beta[1] - 1.96 * se[1])),
                        "high": float(np.exp(beta[1] + 1.96 * se[1])),
                        "p": float(2 * stats.norm.sf(abs(beta[1] / se[1]))),
                    }
                )
    return rows


def run(ctx: Context) -> list[Path]:
    matrix = matrix_path(ctx.paths, GSE)
    table = ctx.paths.input / "geo" / GSE / TABLE
    if not (matrix.exists() and table.exists()):
        print(f"sjlife: {matrix.name} or {TABLE} missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("sjlife_deaths")
    people = first_draws(read_people(matrix).drop_nulls("age"))
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    lay = ClockLayout(ctx.paths)
    clocks = {c: read_clock(lay.file(c)) for c in CLOCKS_USED if lay.file(c).exists()}
    clock_cpgs = sorted({c for w in clocks.values() for c in w["cpg"].to_list()})
    import rdata

    reference = next(iter(rdata.read_rda(str(ctx.paths.input / EPIDISH)).values()))
    ref_cpgs = [str(c) for c in reference.coords["dim_0"].values]
    wanted = sorted(set(features) | set(clock_cpgs) | set(ref_cpgs))
    header = stream_header(table)
    people = people.filter(pl.col("barcode").is_in(header))
    cpgs, x_all = stream_rows(table, set(wanted), people["barcode"].to_list())
    print(
        f"sjlife: {people.height} people, {people['age_at_death'].is_not_null().sum()} deaths, {len(cpgs)} CpGs",
        file=sys.stderr,
        flush=True,
    )
    col = {c: j for j, c in enumerate(cpgs)}
    feats = [f for f in features if f in col]
    x = x_all[:, [col[f] for f in feats]]
    ok = np.isnan(x).mean(axis=0) == 0
    kept = [f for f, k in zip(feats, ok, strict=True) if k]
    scores = reference_scores(x[:, ok], _islands(kept, probes), _sets(kept, membership))
    measures = {f"{m} {g}": scores[f"{m}_{g}"] for g in GROUPS for m in ("F", "B", "S") if f"{m}_{g}" in scores}
    age = people["age"].to_numpy()
    for name, weights in clocks.items():
        present = [c for c in clock_cpgs if c in col]
        score, used = linear_score(x_all[:, [col[c] for c in present]], present, weights)
        if used and np.isfinite(score).all():
            design = np.column_stack([np.ones(age.size), age])
            coef, *_ = np.linalg.lstsq(design, score, rcond=None)
            measures[name] = score - design @ coef
    keep_ref = [i for i, c in enumerate(ref_cpgs) if c in col]
    xr = x_all[:, [col[ref_cpgs[i]] for i in keep_ref]]
    xr = np.where(np.isnan(xr), np.nanmean(xr, axis=0), xr)
    fractions = deconvolve(xr, reference.values[keep_ref])
    extra = np.column_stack([np.nanmean(x_all, axis=1), fractions[:, :-1]])
    print(f"sjlife: {len(keep_ref)} reference CpGs; neutrophils {fractions[:, -2].mean():.2f}", file=sys.stderr)
    result = pl.DataFrame(death_rows(people, measures, extra))
    write_table(result, out / "effects")
    write_table(people.with_columns(**{k: pl.Series(v) for k, v in measures.items()}), out / "per_person")
    return [_figure(ctx, result)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.2)
    measures = table["measure"].unique(maintain_order=True).to_list()
    for i, adjusted in enumerate((False, True)):
        ax = axes[i]
        for y, m in enumerate(measures):
            for j, lagged in enumerate((False, True)):
                r = table.filter(
                    (pl.col("adjusted") == adjusted) & (pl.col("measure") == m) & (pl.col("lag_over_2y") == lagged)
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
        n = table["deaths"].max()
        title = "chip, cell mix, mean beta held" if adjusted else "age, sex, treatment, diagnosis held"
        ax.set(xlabel="odds ratio of death per SD (95 % CI)", title=f"{title} ({n} deaths)")
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "SJLIFE: blood drawn before death; circles all deaths, squares over 2 years after the draw",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("sjlife_deaths"))
