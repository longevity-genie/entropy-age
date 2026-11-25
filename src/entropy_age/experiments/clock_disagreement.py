"""Clock disagreement as a signal (RESEARCH_PLAN idea 4): does the spread between a person's valid age clocks rise
with drift, and does it carry risk beyond their mean?

The panel, per group: every scored clock that is valid in the group's tissue (``clock_atlas.valid``) and estimates
age, that is a first- or second-generation clock (``clock_increment.GENERATION``) or any clock whose output is in
years. Pace (DunedinPACE), mitotic and stochastic clocks are left out: they do not estimate age, and a rate or a
division count placed beside age estimates would make the spread measure the difference in kind. A clock enters only
when its prediction tracks age in the group (r > 0.5, the guard ``clock_increment`` uses: a single-CpG clock such as
Bocklandt runs against age on some arrays and would inflate the spread). Per clock, the acceleration is the residual
of its prediction on age within the group, z-scored within the group. Per person, the disagreement is the SD across
the panel's clocks of these z-scores, and the mean acceleration their mean. A second panel adds the adult age clocks
used outside their tissue (blood clocks in solid tissue), so groups where only Horvath's two clocks are valid still
have a spread; it is labelled and never the primary one.

(a) Drift. Per series group of ``config/series`` (one array per donor, arrays passing QC, at least 25 with ages),
each person's drift at housekeeping look-alike promoters (B) and shield share (F) against a reference fitted on
everyone else in the group, on the CpGs measured in every array (``loo_scores``, as ``clock_drift``). Per group and
panel: the correlation of the disagreement with B and with F, the same for the mean acceleration, and the partial
correlation of the disagreement with B and F given the mean. EPIC-Italy (GSE51032, rebuilt; B and F from
``incident_cancer_rsesame``, leave-one-out over everyone) and SJLIFE (GSE314261; B and F against one reference on
everyone, as ``sjlife_deaths``) are added as blood groups.

(b) Risk. EPIC-Italy women on arrays passing QC (``reprocessed.read_qc``: detection and the low-intensity rule;
clock scores exist only for those arrays), with disagreement and mean computed among the women: breast cancer (C50)
diagnosed more than 2 years after the draw (and, as a check, at any time) against women with no cancer diagnosis.
Unconditional logistic models of case status on disagreement and mean acceleration together, each per SD in the
modelled women, adjusted for age (model "age") and also for mean beta and 11 cell fractions (model "covariates", the
set ``clock_increment`` holds); the same with both measures centred within array chip ("chip-centred"); and a
conditional model on age-matched sets, same chip first (``cancer_matched.match``, as ``clock_increment``), with the
covariates. Wald 95 % intervals. SJLIFE: death among survivors (no control has died), from the clock accelerations
``sjlife_deaths`` stored (Horvath 2013, Hannum, PhenoAge and Zhang 2019, linear scores already residualised on age,
so the age guard is not applied there), adjusted for age, sex, the four main exposures and diagnosis; all deaths and
deaths more than 2 years after the draw. Mean beta and cell fractions are not stored for SJLIFE, so that model has no
covariate version.

Writes ``clock_disagreement/drift`` (per group and panel), ``clock_disagreement/per_person``,
``clock_disagreement/risk`` and a two-panel figure. Not registered as a figure; run as a Slurm job (leave-one-out).
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
from entropy_age.experiments.cancer_matched import clogit, match
from entropy_age.experiments.clock_atlas import valid
from entropy_age.experiments.clock_increment import GENERATION
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.incident_cancer import _group_means, logit
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import groups_of, series_table, specs
from entropy_age.experiments.sjlife_deaths import EXPOSURES
from entropy_age.experiments.tissue_contrast import _features

SERIES = Path("data/interim/clocks/series")
CLOCK_SCORES = Path("data/interim/clocks/series_unmasked")
EPIC_ITALY = "GSE51032"
GROUP = "housekeeping"
B, F = f"B {GROUP}", f"F {GROUP}"
MIN_ARRAYS = 25
MIN_CLOCKS = 2
MIN_AGE_R = 0.5
LAG = 2.0
LEFT_OUT = ("3rd: pace of ageing", "mitotic and stochastic")
PANELS = ("valid", "with clocks outside their tissue")
#: SJLIFE's stored accelerations (sjlife_deaths.CLOCKS_USED without DunedinPACE, a pace); all blood or pan-tissue
SJLIFE_CLOCKS = ("Horvath 2013", "Hannum", "PhenoAge", "Zhang 2019")


def age_clocks(outputs: dict[str, str], group: str, *, outside_tissue: bool = False) -> list[str]:
    """The clocks of ``outputs`` (clock -> output unit) that estimate age and are valid in ``group`` (or, with
    ``outside_tissue``, adult age clocks used anywhere); pace, mitotic and stochastic clocks are left out."""
    kept = []
    for clock, unit in outputs.items():
        gen = GENERATION.get(clock, "other")
        if gen in LEFT_OUT:
            continue
        if not (gen.startswith(("1st", "2nd")) or str(unit).strip().lower() == "age (years)"):
            continue
        if valid(clock, group) or (outside_tissue and gen != "other"):
            kept.append(clock)
    return kept


def disagreement(
    preds: np.ndarray, age: np.ndarray, *, min_age_r: float | None = MIN_AGE_R
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per person (rows of ``preds``, people x clocks): the SD across clocks of the z-scored accelerations (residual
    on age), their mean, and the mask of clocks used (finite, not constant, r with age above ``min_age_r``)."""
    preds = np.asarray(preds, dtype=float)
    age = np.asarray(age, dtype=float)
    use = np.zeros(preds.shape[1], dtype=bool)
    for j in range(preds.shape[1]):
        v = preds[:, j]
        if not np.isfinite(v).all() or np.std(v) == 0:
            continue
        use[j] = min_age_r is None or np.corrcoef(v, age)[0, 1] > min_age_r
    design = np.column_stack([np.ones(age.size), age])
    coef, *_ = np.linalg.lstsq(design, preds[:, use], rcond=None)
    acc = preds[:, use] - design @ coef
    acc = (acc - acc.mean(axis=0)) / acc.std(axis=0)
    if acc.shape[1] < MIN_CLOCKS:
        nan = np.full(age.size, np.nan)
        return nan, nan, use
    return acc.std(axis=1, ddof=1), acc.mean(axis=1), use


def zs(v: np.ndarray) -> np.ndarray:
    return (v - v.mean()) / v.std(ddof=1)


def partial_r(x: np.ndarray, y: np.ndarray, given: np.ndarray) -> float:
    d = np.column_stack([np.ones(len(given)), given])
    rx = x - d @ np.linalg.lstsq(d, x, rcond=None)[0]
    ry = y - d @ np.linalg.lstsq(d, y, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def drift_rows(series: str, group: str, people: pl.DataFrame, b, f, preds, clocks, age, panel, min_age_r) -> tuple:
    """(row of the drift table, per-person frame) for one group and panel, or (None, None) when too few clocks."""
    dis, mean, use = disagreement(preds, age, min_age_r=min_age_r)
    if not np.isfinite(dis).all():
        return None, None
    used = [c for c, u in zip(clocks, use, strict=True) if u]

    def r(x, y):
        return float(np.corrcoef(x, y)[0, 1]) if np.std(y) > 0 else float("nan")

    row = {
        "series": series,
        "group": group,
        "panel": panel,
        "arrays": int(age.size),
        "clocks": len(used),
        "r_dis_B": r(dis, b),
        "r_dis_F": r(dis, f),
        "r_mean_B": r(mean, b),
        "r_mean_F": r(mean, f),
        "partial_r_dis_B": partial_r(dis, b, mean),
        "partial_r_dis_F": partial_r(dis, f, mean),
        "p_dis_B": float(stats.pearsonr(dis, b).pvalue),
        "p_dis_F": float(stats.pearsonr(dis, f).pvalue),
        "clock_names": ", ".join(used),
    }
    frame = people.select("sample").with_columns(
        series=pl.lit(series),
        group=pl.lit(group),
        panel=pl.lit(panel),
        age=pl.Series(age),
        B=pl.Series(b),
        F=pl.Series(f),
        disagreement=pl.Series(dis),
        mean_acceleration=pl.Series(mean),
        clocks=pl.lit(len(used)),
    )
    return row, frame


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
        table = series_table(ctx.paths, spec)
        if not (betas_path.exists() and clocks_path.exists() and table is not None):
            continue
        if not reprocessed.available(ctx.paths, gse):
            continue
        long = pl.read_csv(clocks_path)
        outputs = dict(long.select("clock", "output").unique(subset="clock").iter_rows())
        scores = long.pivot(on="clock", index="sample", values="prediction")
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features))
        usable = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))["sample"].to_list()
        for label, s in groups_of(table, spec, usable, set()).items():
            s = s.drop_nulls("age").filter(
                pl.col("sample").is_in(betas.columns) & pl.col("sample").is_in(scores["sample"])
            )
            if s.height < MIN_ARRAYS:
                continue
            panels = {p: age_clocks(outputs, label, outside_tissue=p != PANELS[0]) for p in PANELS}
            if len(panels[PANELS[1]]) < MIN_CLOCKS:
                continue
            samples = s["sample"].to_list()
            sub = betas.select("cpg", *samples)
            x = sub.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                continue
            print(f"clock_disagreement: {name} / {label}: {len(samples)} arrays", file=sys.stderr, flush=True)
            per = loo_scores(x[:, ok], _islands(kept, probes), sets)
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            age = s["age"].to_numpy().astype(float)
            clk = s.select("sample").join(scores, on="sample", how="left")
            for panel, clocks in panels.items():
                if len(clocks) < MIN_CLOCKS:
                    continue
                preds = clk.select(clocks).to_numpy().astype(float)
                row, frame = drift_rows(name, label, s, b, f, preds, clocks, age, panel, MIN_AGE_R)
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
    passed = reprocessed.read_qc(ctx.paths, EPIC_ITALY).filter(pl.col("passed"))["sample"].to_list()
    long = pl.read_csv(clocks_path)
    outputs = dict(long.select("clock", "output").unique(subset="clock").iter_rows())
    clocks = age_clocks(outputs, "blood")
    wide = long.pivot(on="clock", index="sample", values="prediction").select("sample", *clocks)
    fractions = [f"fraction_{i}" for i in range(1, 12)]
    people = (
        pl.read_parquet(src)
        .select("sample", "female", "age", "icd", "years_to_diagnosis", "chip", "mean_beta", *fractions, B, F)
        .filter(pl.col("sample").is_in(passed))
        .join(wide, on="sample", how="inner", validate="1:1")
        .drop_nulls(["age", B, F])
    )
    rows, frames = [], []
    for label, part in (("blood, everyone", people), ("blood, women", people.filter(pl.col("female")))):
        row, frame = drift_rows(
            "EPIC-Italy (GSE51032)",
            label,
            part,
            part[B].to_numpy(),
            part[F].to_numpy(),
            part.select(clocks).to_numpy().astype(float),
            clocks,
            part["age"].to_numpy().astype(float),
            PANELS[0],
            MIN_AGE_R,
        )
        if row is not None:
            rows.append(row)
            frames.append(frame)
    women = people.filter(pl.col("female"))
    dis, mean, use = disagreement(women.select(clocks).to_numpy().astype(float), women["age"].to_numpy())
    women = women.with_columns(disagreement=pl.Series(dis), mean_acceleration=pl.Series(mean))
    used = [c for c, u in zip(clocks, use, strict=True) if u]
    print(f"clock_disagreement: EPIC-Italy women {women.height}, clocks {', '.join(used)}", file=sys.stderr)
    return women, rows, frames, breast_rows(women, fractions, len(used))


def _or_row(cohort, outcome, lag, model, terms, beta, se, cases, controls, n_clocks) -> list[dict]:
    out = []
    for i, term in enumerate(terms):
        out.append(
            {
                "cohort": cohort,
                "outcome": outcome,
                "lag": lag,
                "model": model,
                "term": term,
                "cases": cases,
                "controls": controls,
                "clocks": n_clocks,
                "odds_ratio": float(np.exp(beta[i])),
                "low": float(np.exp(beta[i] - 1.96 * se[i])),
                "high": float(np.exp(beta[i] + 1.96 * se[i])),
                "p": float(2 * stats.norm.sf(abs(beta[i] / se[i]))),
            }
        )
    return out


def breast_rows(women: pl.DataFrame, fractions: list[str], n_clocks: int) -> list[dict]:
    breast = pl.col("icd").fill_null("").str.starts_with("C50")
    controls = women.filter(pl.col("icd").is_null())
    rows = []
    terms = ("disagreement", "mean acceleration")
    for lag, cond in ((f"> {LAG:g} y", pl.col("years_to_diagnosis") > LAG), ("any", pl.lit(True))):
        cases = women.filter(breast & cond)
        pool = pl.concat([cases, controls], how="vertical")
        y = np.r_[np.ones(cases.height), np.zeros(controls.height)]
        dis, mean = zs(pool["disagreement"].to_numpy()), zs(pool["mean_acceleration"].to_numpy())
        chip = pool["chip"].to_numpy()
        dis_c, mean_c = zs(dis - _group_means(dis, chip)), zs(mean - _group_means(mean, chip))
        age = zs(pool["age"].to_numpy())
        cov = np.column_stack([zs(pool["mean_beta"].to_numpy()), *(zs(pool[c].to_numpy()) for c in fractions)])
        one = np.ones(pool.height)
        designs = {
            "age": np.column_stack([one, dis, mean, age]),
            "covariates": np.column_stack([one, dis, mean, age, cov]),
            "covariates, chip-centred": np.column_stack([one, dis_c, mean_c, age, cov]),
        }
        for model, design in designs.items():
            beta, se = logit(y, design)
            rows += _or_row(
                "EPIC-Italy women",
                "breast cancer",
                lag,
                model,
                terms,
                beta[1:],
                se[1:],
                cases.height,
                controls.height,
                n_clocks,
            )
        beta, se = logit(y, np.column_stack([one, dis, age, cov]))
        rows += _or_row(
            "EPIC-Italy women",
            "breast cancer",
            lag,
            "covariates, disagreement alone",
            terms[:1],
            beta[1:],
            se[1:],
            cases.height,
            controls.height,
            n_clocks,
        )
        sets = match(cases, controls)
        beta, se = clogit(np.column_stack([dis, mean, age, cov]), sets)
        rows += _or_row(
            "EPIC-Italy women",
            "breast cancer",
            lag,
            "matched sets, covariates",
            terms,
            beta,
            se,
            len(sets),
            sum(len(s) - 1 for s in sets),
            n_clocks,
        )
    return rows


def sjlife(ctx: Context) -> tuple[list[dict], list[pl.DataFrame], list[dict]]:
    src = ctx.paths.output / "sjlife_deaths" / "per_person.parquet"
    if not src.exists():
        print("clock_disagreement: no sjlife_deaths per_person; SJLIFE skipped", file=sys.stderr)
        return [], [], []
    people = pl.read_parquet(src)
    clocks = [c for c in SJLIFE_CLOCKS if c in people.columns]
    people = people.drop_nulls(["age", B, F, *clocks])
    age = people["age"].to_numpy().astype(float)
    preds = people.select(clocks).to_numpy().astype(float)
    row, frame = drift_rows(
        "SJLIFE (GSE314261)",
        "PBMC, everyone",
        people,
        people[B].to_numpy(),
        people[F].to_numpy(),
        preds,
        clocks,
        age,
        PANELS[0],
        None,
    )
    if row is None:
        return [], [], []
    dis, mean, _ = disagreement(preds, age, min_age_r=None)
    died = people["age_at_death"].is_not_null().to_numpy()
    years = (people["age_at_death"].fill_null(np.inf) - people["age"]).to_numpy()
    survivor = people["survivor"].to_numpy()
    rows = []
    for lag, keep in (("any", survivor), (f"> {LAG:g} y", survivor & (~died | (years > LAG)))):
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
        beta, se = logit(y, design)
        rows += _or_row(
            "SJLIFE survivors",
            "death",
            lag,
            "age, sex, exposures, diagnosis",
            ("disagreement", "mean acceleration"),
            beta[1:],
            se[1:],
            int(y.sum()),
            int((1 - y).sum()),
            len(clocks),
        )
    return [row], [frame], rows


def run(ctx: Context) -> list[Path]:
    women, e_rows, e_frames, risk = epic_italy(ctx)
    s_rows, s_frames, s_risk = sjlife(ctx)
    rows, frames = series_groups(ctx)
    rows, frames, risk = rows + e_rows + s_rows, frames + e_frames + s_frames, risk + s_risk
    if not rows:
        print("clock_disagreement: nothing to score", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_disagreement")
    drift = pl.DataFrame(rows)
    write_table(drift, out / "drift")
    write_table(pl.concat(frames, how="vertical"), out / "per_person", csv_max_rows=10_000)
    risk_table = pl.DataFrame(risk) if risk else pl.DataFrame()
    if risk:
        write_table(risk_table, out / "risk")
    if women is not None:
        write_table(
            women.select("sample", "age", "chip", "icd", "years_to_diagnosis", "disagreement", "mean_acceleration"),
            out / "epic_italy_women",
        )
    return [draw(drift, risk_table, ctx.figure_path("clock_disagreement"))]


def draw(drift: pl.DataFrame, risk: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    ax = axes[0]
    for i, panel in enumerate(PANELS):
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
        pick = risk.filter(pl.col("model").is_in(["covariates", "age, sex, exposures, diagnosis"]))
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
