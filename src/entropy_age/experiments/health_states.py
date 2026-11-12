"""Does the shield track health? Cross-sectional tests on cohorts already on disk (roadmap RM24)

Asked by antonkulaga (2026-11-12): with the data already downloaded, can we estimate whether the shield predicts
health? None of these cohorts follows people forward, so this asks the cross-sectional question: at equal age, do
ill people have a weaker shield (higher ``F``), more drift (``B``) or a shifted set (``S``)?

- Lunnon (GSE59685): Alzheimer's disease against controls and Braak stage 0-VI (continuous), in whole blood and
  four brain regions of the same donors; age and sex as covariates.
- Jaffe (GSE74193): schizophrenia against controls in prefrontal cortex, adults (17 and over), one sample per
  donor; age, sex, race, neuron fraction and plate as covariates.
- GSE41037 (blood, 27k array): schizophrenia against controls; age and sex.
- GTEx (8 tissues): current against never smokers; age and sex.

Every person is scored against all others in the same contrast (leave-one-out reference, as controls are in
``consortium_interventions.arm_scores``), so small control groups (Lunnon blood has 9) do not set the reference
alone; a reference that mixes cases pulls effects toward zero, never away from it. Scored raw and with 5
composition components of the background CpGs removed. A weaker shield in cases may be a consequence of disease or
its treatment, not a predictor.

``entropy-age figures health`` (as a Slurm job).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _islands, _sets, score_against
from entropy_age.experiments.syndromes import matrix_path, read_betas
from entropy_age.experiments.tissue_contrast import _features, _matrix
from entropy_age.protection import fit_reference

COMPONENTS = (0, 5)
MEASURES = ("B", "F", "S")
MIN_GROUP = 8
BRAAK = {"0": 0, "1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6}


@dataclass(frozen=True, slots=True)
class Contrast:
    """One health contrast: ``x`` people by CpGs, ``outcome`` (0/1 or continuous) and ``covariates``."""

    cohort: str
    tissue: str
    outcome: str
    x: np.ndarray
    y: np.ndarray
    covariates: np.ndarray
    features: list[str]
    binary: bool


def loo_scores(x: np.ndarray, island: np.ndarray, sets: dict, *, components: int = 0) -> dict[str, np.ndarray]:
    """Every person's scores against a reference fitted on everyone else (``components`` removed first)."""
    member = np.any([s for s, _ in sets.values()], axis=0)
    rows = []
    for i in range(x.shape[0]):
        others = np.arange(x.shape[0]) != i
        ref = fit_reference(x[others], island, ~member, components=components)
        rows.append({k: float(v[0]) for k, v in score_against(x[i : i + 1], ref, sets).items()})
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def outcome_effects(scores: dict[str, np.ndarray], y: np.ndarray, covariates: np.ndarray) -> list[dict]:
    """Per set and measure: the outcome's OLS coefficient with covariates, its SE and p, standardised effect."""
    design = np.column_stack([np.ones(y.size), y.astype(float), covariates])
    rows = []
    for key, v in scores.items():
        measure, _, group = key.partition("_")
        if measure not in MEASURES or not np.isfinite(v).all():
            continue
        coef, se, p = ols(v, design, 1)
        rows.append(
            {"group": group, "measure": measure, "effect": coef, "se": se, "p": p, "effect_sd": coef / v.std(ddof=1)}
        )
    return rows


def _covariates(frame: pl.DataFrame, numeric: list[str], categorical: list[str]) -> np.ndarray:
    parts = [frame[c].to_numpy().astype(float)[:, None] for c in numeric]
    parts += [dummies(frame[c].to_numpy()) for c in categorical]
    return np.column_stack(parts)


def _lunnon(paths, features) -> list[Contrast]:
    s = pl.read_parquet(cohort_layout("lunnon", paths).samples).drop_nulls("age")
    out = []
    for (tissue,), t in s.group_by("tissue", maintain_order=True):
        t = t.sort("donor")
        m = _matrix("lunnon", paths, t["sample"].to_list(), features)
        ad = t.with_row_index("i").filter(pl.col("ad_status").is_in(["AD", "C"]))
        out.append(
            Contrast(
                "Lunnon",
                tissue,
                "Alzheimer's (AD vs control)",
                m.values[ad["i"].to_numpy()],
                (ad["ad_status"] == "AD").to_numpy(),
                _covariates(ad, ["age"], ["sex"]),
                m.features,
                True,
            )
        )
        br = t.with_row_index("i").filter(pl.col("braak_stage").is_in(list(BRAAK)))
        out.append(
            Contrast(
                "Lunnon",
                tissue,
                "Braak stage (per stage)",
                m.values[br["i"].to_numpy()],
                br["braak_stage"].replace_strict(BRAAK).to_numpy().astype(float),
                _covariates(br, ["age"], ["sex"]),
                m.features,
                False,
            )
        )
    return out


def _jaffe(paths, features) -> list[Contrast]:
    s = (
        pl.read_parquet(cohort_layout("jaffe", paths).samples)
        .filter((pl.col("age") >= 17) & ~pl.col("drop_sample") & pl.col("group").is_in(["Control", "Schizo"]))
        .sort("best_qc", descending=True)
        .unique(subset="donor", keep="first", maintain_order=True)
        .drop_nulls(["age", "neuron_fraction"])
        .sort("sample")
    )
    m = _matrix("jaffe", paths, s["sample"].to_list(), features)
    cov = _covariates(s, ["age", "neuron_fraction"], ["sex", "race", "plate"])
    return [
        Contrast(
            "Jaffe",
            "prefrontal cortex",
            "schizophrenia vs control",
            m.values,
            (s["group"] == "Schizo").to_numpy(),
            cov,
            m.features,
            True,
        )
    ]


def _gse41037(paths, features) -> list[Contrast]:
    from entropy_age.data.cohorts import _characteristics

    matrix = matrix_path(paths, "GSE41037")
    if not matrix.exists():
        return []
    chars = _characteristics(matrix)
    s = chars.select(
        "sample",
        pl.col("gender").alias("sex"),
        pl.col("age").cast(pl.Float64, strict=False),
        pl.col("diseasestatus").str.starts_with("2").alias("case"),
    ).drop_nulls("age")
    x = read_betas(matrix, features, s["sample"].to_list())
    ok = np.isnan(x).mean(axis=0) == 0
    return [
        Contrast(
            "GSE41037",
            "whole blood",
            "schizophrenia vs control",
            x[:, ok],
            s["case"].to_numpy(),
            _covariates(s, ["age"], ["sex"]),
            [f for f, k in zip(features, ok, strict=True) if k],
            True,
        )
    ]


def _gtex(paths, features) -> list[Contrast]:
    s = pl.read_parquet(cohort_layout("gtex", paths).samples).filter(pl.col("smoker").is_in(["Current", "Never"]))
    out = []
    for (tissue,), t in s.group_by("tissue", maintain_order=True):
        t = t.sort("sample")
        m = _matrix("gtex", paths, t["sample"].to_list(), features)
        cov = _covariates(t, ["age"], ["sex"]) if t["sex"].n_unique() > 1 else _covariates(t, ["age"], [])
        out.append(
            Contrast(
                "GTEx",
                tissue,
                "current vs never smoker",
                m.values,
                (t["smoker"] == "Current").to_numpy(),
                cov,
                m.features,
                True,
            )
        )
    return out


LOADERS = (_lunnon, _jaffe, _gse41037, _gtex)


def _enough(c: Contrast) -> bool:
    if c.binary:
        return min(int(c.y.sum()), int((~c.y).sum())) >= MIN_GROUP
    return c.y.size >= 2 * MIN_GROUP and np.unique(c.y).size >= 3


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("health_states")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    rows, scored = [], []
    for loader in LOADERS:
        for c in loader(ctx.paths, features):
            if not _enough(c):
                print(f"health: {c.cohort} {c.tissue} {c.outcome}: too few; left out", file=sys.stderr)
                continue
            island = _islands(c.features, probes)
            sets = _sets(c.features, membership)
            for k in COMPONENTS:
                scores = loo_scores(c.x, island, sets, components=k)
                for r in outcome_effects(scores, c.y, c.covariates):
                    rows.append(
                        {
                            "cohort": c.cohort,
                            "tissue": c.tissue,
                            "outcome": c.outcome,
                            "components": k,
                            "n": int(c.y.size),
                            "cases": int(c.y.sum()) if c.binary else None,
                            **r,
                        }
                    )
                scored.append(
                    pl.DataFrame(
                        {
                            "cohort": c.cohort,
                            "tissue": c.tissue,
                            "outcome": c.outcome,
                            "components": k,
                            "y": c.y.astype(float),
                            **{key: v for key, v in scores.items() if key.split("_")[0] in MEASURES},
                        }
                    )
                )
            print(f"health: {c.cohort} {c.tissue} {c.outcome}: n={c.y.size} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    write_table(pl.concat(scored, how="diagonal"), out / "person_scores")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    """Standardised effects of each health state on F, B, S (housekeeping and core essential), raw and adjusted."""
    contrasts = table.select("cohort", "tissue", "outcome").unique(maintain_order=True)
    fig, axes = plots.panels(3, ncols=3, size=4.6)
    groups = [g for g in ("core_essential", "housekeeping") if g in table["group"].unique()]
    for i, measure in enumerate(MEASURES):
        ax = axes[i]
        labels = []
        for y, (cohort, tissue, outcome) in enumerate(contrasts.iter_rows()):
            for g, group in enumerate(groups):
                for j, k in enumerate(COMPONENTS):
                    r = table.filter(
                        (pl.col("cohort") == cohort)
                        & (pl.col("tissue") == tissue)
                        & (pl.col("outcome") == outcome)
                        & (pl.col("group") == group)
                        & (pl.col("measure") == measure)
                        & (pl.col("components") == k)
                    )
                    if r.is_empty():
                        continue
                    sd = r["effect"][0] / r["effect_sd"][0] if r["effect_sd"][0] else 1.0
                    ax.errorbar(
                        r["effect_sd"][0],
                        y + 0.2 * g + 0.1 * j,
                        xerr=1.96 * r["se"][0] / abs(sd),
                        fmt="os"[j],
                        ms=3.5,
                        color=plots.PALETTE[g],
                        alpha=1.0 if k else 0.45,
                    )
            labels.append((y + 0.15, f"{cohort} {tissue}: {outcome}"))
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels] if i == 0 else [""] * len(labels), fontsize=6.5)
        ax.invert_yaxis()
        ax.set(
            xlabel="Effect in SD of the score (95 % CI)",
            title={"F": "Shield share F (> 0: weaker)", "B": "Drift B", "S": "Signed shift S"}[measure],
        )
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Health states and the shield, cross-sectional (circles raw, squares 5 components)",
        x=0.02,
        y=1.03,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("health_states"))
