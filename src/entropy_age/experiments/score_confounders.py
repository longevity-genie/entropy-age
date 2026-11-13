"""Is the per-person shield score confounded the way outlier burden is? (roadmap RM13 (5), Seeboth 2020)

Seeboth et al. 2020 found a person's burden of methylation outliers driven by smoking, batch and cell mix, and
not by health. Our drift ``B`` (how far a person's look-alike CpGs sit from the young reference) is a burden of
that kind; the shield share ``F`` and the signed shift ``S`` are contrasts within one sample, set against matched
look-alikes, which should cancel what moves all CpGs alike. From the per-person tables already written:

* **SATSA** (``figures within_person``; 1,469 draws of 596 people): smoking (current, former) and the array type
  (450k or EPIC);
* **SJLIFE** (``figures sjlife_deaths``; 4,518 people): the array slide and its row (batch and position, from the
  barcode), the diagnosis (each childhood cancer, or community control), and alkylating agents and anthracyclines
  in the survivors' treatment, within diagnosis. Alkylating chemotherapy shrinks the blood stem-cell pool, a
  bottleneck: selection predicts a stronger shield (lower ``F``) after it.

Each measure is standardised, age and sex held. Per block, its partial R² (the share of what age and sex leave
that the block explains, beside ``chance_r2``, what as many random columns would explain); per term, its
coefficient in SD with an interval from resampling people.

``entropy-age figures score_confounders`` (a few minutes; reads two small tables).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, dummies, write_table

GROUPS = ("core_essential", "housekeeping")
MEASURES = ("B", "F", "S")
N_BOOT = {"paper": 500, "quick": 50}


def _rss(y: np.ndarray, design: np.ndarray) -> float:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    r = y - design @ coef
    return float(r @ r)


def partial_r2(y: np.ndarray, base: np.ndarray, block: np.ndarray) -> float:
    """Share of the residual sum of squares of ``y ~ base`` that adding ``block`` removes."""
    reduced = _rss(y, base)
    return (reduced - _rss(y, np.column_stack([base, block]))) / reduced if reduced > 0 else float("nan")


def coefficients(y: np.ndarray, design: np.ndarray, columns: list[int]) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coef[columns]


def _satsa(path: Path) -> tuple[pl.DataFrame, dict, dict] | None:
    if not path.exists():
        return None
    t = pl.read_parquet(path).filter(pl.col("smoke") != "NA")
    smoke = t["smoke"].to_numpy()
    blocks = {
        "smoking": np.column_stack([smoke == "Current", smoke == "Previous"]).astype(float),
        "array type": t["epic"].to_numpy().astype(float)[:, None],
    }
    terms = {"current smoker": ("smoking", 0), "EPIC array": ("array type", 0)}
    return t, blocks, terms


def _sjlife(path: Path) -> tuple[pl.DataFrame, dict, dict] | None:
    if not path.exists():
        return None
    t = pl.read_parquet(path).drop_nulls(["age", "female"])
    slide = t["barcode"].str.extract(r"^(\d+)_").to_numpy()
    row = t["barcode"].str.extract(r"_R(\d+)C").to_numpy()
    alk = t["treatment_alkylating_classic_5"].fill_null(0).to_numpy().astype(float)
    anth = t["treatment_anthracyclines_5"].fill_null(0).to_numpy().astype(float)
    blocks = {
        "slide (batch)": dummies(slide),
        "row on slide": dummies(row),
        "diagnosis": dummies(t["diagnosis"].fill_null("unknown").to_numpy()),
        "chemotherapy": np.column_stack([alk, anth]),
    }
    terms = {
        "alkylating agents": ("chemotherapy", 0),
        "anthracyclines": ("chemotherapy", 1),
    }
    return t, blocks, terms


def analyse(cohort: str, t: pl.DataFrame, blocks: dict, terms: dict, rng, n_boot: int) -> tuple[list, list]:
    age, female = t["age"].to_numpy().astype(float), t["female"].to_numpy().astype(float)
    base = np.column_stack([np.ones(t.height), age, female])
    full = np.column_stack([base, *blocks.values()])
    offsets, start = {}, base.shape[1]
    for name, block in blocks.items():
        offsets[name] = start
        start += block.shape[1]
    term_cols = [offsets[b] + k for b, k in terms.values()]
    people = t["person"].to_numpy() if "person" in t.columns else np.arange(t.height)
    uniq = np.unique(people)
    rows_of = {p: np.flatnonzero(people == p) for p in uniq}
    r2_rows, term_rows = [], []
    for group in GROUPS:
        for measure in MEASURES:
            col = f"{measure} {group}" if f"{measure} {group}" in t.columns else f"{measure}_{group}"
            if col not in t.columns:
                continue
            y = t[col].to_numpy().astype(float)
            ok = np.isfinite(y)
            y = (y - y[ok].mean()) / y[ok].std()
            yb, fb, bb = y[ok], full[ok], base[ok]
            for name, block in blocks.items():
                r2_rows.append(
                    {
                        "cohort": cohort,
                        "group": group,
                        "measure": measure,
                        "block": name,
                        "partial_r2": partial_r2(yb, bb, block[ok]),
                        "chance_r2": block.shape[1] / int(ok.sum()),
                    }
                )
            est = coefficients(yb, fb, term_cols)
            draws = []
            for _ in range(n_boot):
                pick = np.concatenate([rows_of[p] for p in rng.choice(uniq, uniq.size)])
                pick = pick[ok[pick]]
                draws.append(coefficients(y[pick], full[pick], term_cols))
            draws = np.array(draws)
            for k, term in enumerate(terms):
                term_rows.append(
                    {
                        "cohort": cohort,
                        "group": group,
                        "measure": measure,
                        "term": term,
                        "sd_effect": float(est[k]),
                        "low": float(np.quantile(draws[:, k], 0.025)),
                        "high": float(np.quantile(draws[:, k], 0.975)),
                        "n": int(ok.sum()),
                    }
                )
    return r2_rows, term_rows


def run(ctx: Context) -> list[Path]:
    sources = {
        "SATSA": _satsa(ctx.paths.output / "within_person" / "per_draw.parquet"),
        "SJLIFE": _sjlife(ctx.paths.output / "sjlife_deaths" / "per_person.parquet"),
    }
    if all(v is None for v in sources.values()):
        print("score_confounders: no per-person tables; run within_person and sjlife_deaths first", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("score_confounders")
    rng = ctx.rng(270_000)
    r2, terms = [], []
    for cohort, src in sources.items():
        if src is None:
            continue
        a, b = analyse(cohort, *src, rng, N_BOOT[str(ctx.scale.name)])
        r2 += a
        terms += b
    r2, terms = pl.DataFrame(r2), pl.DataFrame(terms)
    write_table(r2, out / "partial_r2")
    write_table(terms, out / "terms")
    return [_figure(ctx, r2, terms)]


def _figure(ctx: Context, r2: pl.DataFrame, terms: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    ax = axes[0]
    sub = r2.filter(pl.col("group") == "housekeeping")
    keys = sub.select("cohort", "block").unique(maintain_order=True)
    y = np.arange(keys.height)
    for i, measure in enumerate(MEASURES):
        t = keys.join(sub.filter(pl.col("measure") == measure), on=["cohort", "block"], how="left")
        ax.barh(y + 0.27 * i, 100 * t["partial_r2"], height=0.25, color=plots.PALETTE[i], label=measure)
    chance = keys.join(sub.filter(pl.col("measure") == "B"), on=["cohort", "block"], how="left")["chance_r2"]
    ax.scatter(100 * chance, y + 0.27, marker="|", s=300, color="k", zorder=4, label="chance")
    ax.set_yticks(y + 0.27, [f"{c} · {b}" for c, b in keys.iter_rows()])
    ax.invert_yaxis()
    ax.set(xlabel="partial R² beyond age and sex (%)", title="What the score picks up (housekeeping)")
    ax.legend(fontsize=7, title="measure", title_fontsize=7)
    plots.label_panel(ax, "a")
    ax = axes[1]
    sub = terms.filter(pl.col("group") == "housekeeping")
    keys = sub.select("cohort", "term").unique(maintain_order=True)
    y = np.arange(keys.height)
    for i, measure in enumerate(MEASURES):
        t = keys.join(sub.filter(pl.col("measure") == measure), on=["cohort", "term"], how="left")
        ax.errorbar(
            t["sd_effect"],
            y + 0.2 * (i - 1),
            xerr=[t["sd_effect"] - t["low"], t["high"] - t["sd_effect"]],
            fmt="o",
            ms=4,
            lw=1,
            color=plots.PALETTE[i],
        )
    ax.set_yticks(y, [f"{c} · {term}" for c, term in keys.iter_rows()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="effect (SD); age, sex and the other blocks held", title="Smoking, array, chemotherapy")
    plots.label_panel(ax, "b")
    fig.subplots_adjust(wspace=0.9)
    return plots.save(fig, ctx.figure_path("score_confounders"))
