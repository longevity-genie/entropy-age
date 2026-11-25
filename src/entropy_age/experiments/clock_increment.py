"""Does the shield add to a validated clock? EPIC-Italy breast cancer (protocol P4 and change 2).

Women of EPIC-Italy (GSE51032, rebuilt), breast cancer diagnosed more than 2 years after the draw against women
without a cancer diagnosis, matched on age within 3 years as in ``cancer_matched``: primary design same chip first;
sensitivity design same chip only (a set's intercept then absorbs chip). Clocks are scored with biolearn at the
pinned commit (scripts/clocks_biolearn.py) and kept when at least 90 % of their weight is measured in 90 % of arrays.
Per clock, conditional logistic models on identical matched sets, every predictor per SD among all women, adjusted
for age, mean beta and 11 cell fractions:

* M0: clock acceleration (the clock's residual on age among all women);
* M1: M0 + B (drift at housekeeping look-alike promoters);
* M2: M1 + F (the housekeeping shield share, F = P / B);
* M2': M1 + P (the absolute housekeeping gap), a sensitivity check of the ratio.

Reported per clock and design: the clock's odds ratio in M0, B's in M1, F's in M2 and P's in M2' (bootstrap over
matched sets; failed fits are left out and counted), and the likelihood-ratio tests. PhenoAge is the primary clock;
every other kept clock is exploratory. Disagreement: per woman, the SD across an explicit adult-blood panel of
chronological and phenotypic clocks of their accelerations (residual on age, cell fractions, mean beta and chip,
standardised); it enters a model beside the panel mean (exploratory) and its correlations with B and F are reported.
``entropy-age figures clock_increment`` (needs data/interim/clocks/GSE51032_clocks_unmasked.csv).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.cancer_matched import N_BOOT, _score, clogit, match
from entropy_age.experiments.common import Context, dummies, write_table

PRIMARY = "PhenoAge"
MIN_AGE_R = 0.5
LAG = 2.0
#: the clock track (protocol change 4): SeSAMe without the design-flag mask, where 44 clocks pass coverage against 35
#: masked (DunedinPACE, GrimAge, DunedinPoAm38 and DNAmTL among the gain); clocks passing both agree at r >= 0.997
CLOCKS = Path("data/interim/clocks/GSE51032_clocks_unmasked.csv")
B, F, P = "B housekeeping", "F housekeeping", "P housekeeping"
#: adult blood chronological and phenotypic clocks for the disagreement score (gestational, paediatric, placental,
#: cortical and trait predictors are left out)
ADULT_BLOOD = (
    "Horvathv1",
    "Horvathv2",
    "Hannum",
    "Lin",
    "PhenoAge",
    "HRSInCHPhenoAge",
    "GrimAgeV1",
    "GrimAgeV2",
    "YingCausAge",
    "YingDamAge",
    "YingAdaptAge",
    "VidalBralo",
    "Weidner",
    "Garagnani",
    "Bocklandt",
)
#: clock generations (Anton, 2026-11-24): 1st trained on chronological age; 2nd on phenotypic age, mortality or
#: causal/damage components; 3rd on the pace of ageing; plus mitotic and stochastic clocks; any other kept model
#: (traits, gestational, cortical) is "other"
GENERATION = {
    **dict.fromkeys(("Horvathv1", "Horvathv2", "Hannum", "Lin", "VidalBralo", "Weidner", "Garagnani", "Bocklandt"),
                    "1st: chronological"),
    **dict.fromkeys(("PhenoAge", "HRSInCHPhenoAge", "Zhang_10", "GrimAgeV1", "GrimAgeV2", "YingCausAge", "YingDamAge",
                     "YingAdaptAge"), "2nd: phenotypic and mortality"),
    **dict.fromkeys(("DunedinPACE", "DunedinPoAm38"), "3rd: pace of ageing"),
    **dict.fromkeys(("EpiTOC1", "epiTOC2", "epiTOC3", "HypoClock", "RepliTali", "epiCMIT", "StocH", "StocZ", "StocP"), "mitotic and stochastic"),
}
META = ["sample", "female", "age", "icd", "years_to_diagnosis", "chip", "mean_beta", B, F]


def z(v: np.ndarray) -> np.ndarray:
    return (v - np.nanmean(v)) / np.nanstd(v)


def residual(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(y.size), design])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def fit(design: np.ndarray, sets: list[list[int]], steps: int = 50) -> tuple[np.ndarray, float] | None:
    """Coefficients and log-likelihood, or None when the fit is not finite (a failed or rank-deficient fit)."""
    beta, se = clogit(design, sets, steps=steps)
    if not (np.all(np.isfinite(beta)) and np.all(np.isfinite(se))):
        return None
    return beta, float(_score(design, sets, beta)[0])


def models(pool: pl.DataFrame, clock_cols: list[str], sets, base, rng, n_boot: int) -> dict:
    first = np.column_stack([pool[c].to_numpy() for c in clock_cols])
    k = first.shape[1]
    designs = {
        "m0": np.column_stack([first, base]),
        "m1": np.column_stack([first, pool[B].to_numpy(), base]),
        "m2": np.column_stack([first, pool[B].to_numpy(), pool[F].to_numpy(), base]),
        "m2p": np.column_stack([first, pool[B].to_numpy(), pool[P].to_numpy(), base]),
    }
    fits = {name: fit(d, sets) for name, d in designs.items()}
    if any(v is None for v in fits.values()):
        return {}
    draws = {"clock": [], "B": [], "F": [], "P": []}
    for _ in range(n_boot):
        pick = [sets[i] for i in rng.integers(0, len(sets), len(sets))]
        got = {name: fit(d, pick, steps=20) for name, d in designs.items()}
        if any(v is None for v in got.values()):
            continue
        draws["clock"].append(got["m0"][0][0])
        draws["B"].append(got["m1"][0][k])
        draws["F"].append(got["m2"][0][k + 1])
        draws["P"].append(got["m2p"][0][k + 1])

    def ci(v):
        return (float(np.exp(np.quantile(v, 0.025))), float(np.exp(np.quantile(v, 0.975)))) if v else (np.nan, np.nan)

    (b0, l0), (b1, l1), (b2, l2), (b2p, l2p) = fits["m0"], fits["m1"], fits["m2"], fits["m2p"]
    return {
        "clock_or": float(np.exp(b0[0])),
        "clock_low": ci(draws["clock"])[0],
        "clock_high": ci(draws["clock"])[1],
        "second_or": float(np.exp(b0[1])) if k == 2 else float("nan"),
        "B_or_given_clock": float(np.exp(b1[k])),
        "B_low": ci(draws["B"])[0],
        "B_high": ci(draws["B"])[1],
        "lr_p_B": float(stats.chi2.sf(2 * (l1 - l0), 1)),
        "F_or_given_B": float(np.exp(b2[k + 1])),
        "F_low": ci(draws["F"])[0],
        "F_high": ci(draws["F"])[1],
        "lr_p_F": float(stats.chi2.sf(2 * (l2 - l1), 1)),
        "P_or_given_B": float(np.exp(b2p[k + 1])),
        "P_low": ci(draws["P"])[0],
        "P_high": ci(draws["P"])[1],
        "lr_p_P": float(stats.chi2.sf(2 * (l2p - l1), 1)),
        "valid_resamples": len(draws["clock"]),
    }


def run(ctx: Context) -> list[Path]:
    src = ctx.paths.output / "incident_cancer_rsesame" / "per_person.parquet"
    clocks_path = ctx.paths.root / CLOCKS
    if not (src.exists() and clocks_path.exists()):
        print("clock_increment: needs incident_cancer_rsesame per_person and the clock panel", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_increment")
    people = pl.read_parquet(src)
    fractions = [c for c in people.columns if c.startswith("fraction_")]
    people = people.select(*META, *fractions)  # only metadata and B/F: clock columns come from the panel alone
    panel = pl.read_csv(clocks_path).pivot(on="clock", index="sample", values="prediction")
    clocks = [c for c in panel.columns if c != "sample"]
    panel = panel.rename({c: f"clock {c}" for c in clocks})
    women = people.filter(pl.col("female")).join(panel, on="sample", how="inner", validate="1:1")
    women = women.with_columns(pl.Series(P, women[F].to_numpy() * women[B].to_numpy()))  # F = P / B
    finite = np.all(np.isfinite(women.select([f"clock {c}" for c in clocks]).to_numpy()), axis=1)
    women = women.filter(pl.Series(finite))
    age = women["age"].to_numpy()
    tech = np.column_stack(
        [
            age,
            women["mean_beta"].to_numpy(),
            *(women[c].to_numpy() for c in fractions),
            dummies(women["chip"].to_numpy()),
        ]
    )
    for c in clocks:  # model accelerations: residual on age; disagreement accelerations: on the technical design too
        women = women.with_columns(
            pl.Series(f"acc {c}", z(residual(women[f"clock {c}"].to_numpy(), age[:, None]))),
            pl.Series(f"accd {c}", z(residual(women[f"clock {c}"].to_numpy(), tech))),
        )
    # a clock enters the disagreement panel only when it tracks age here (r > 0.5): a single-CpG clock such as
    # Bocklandt runs against age on these arrays and would inflate the spread (figures review 2026-11-24)
    panel_clocks = [
        c for c in ADULT_BLOOD if c in clocks and np.corrcoef(women[f"clock {c}"].to_numpy(), age)[0, 1] > MIN_AGE_R
    ]
    accd = np.column_stack([women[f"accd {c}"].to_numpy() for c in panel_clocks])
    women = women.with_columns(disagreement=pl.Series(z(accd.std(axis=1))), panel_mean=pl.Series(z(accd.mean(axis=1))))
    for name in (B, F, P):
        women = women.with_columns(pl.Series(name, z(women[name].to_numpy())))
    write_table(
        women.select("sample", "age", "chip", B, F, P, "disagreement", "panel_mean", *[f"acc {c}" for c in clocks]),
        out / "per_woman",
    )
    breast = pl.col("icd").fill_null("").str.starts_with("C50")
    cases = women.filter(breast & (pl.col("years_to_diagnosis") > LAG))
    controls = women.filter(pl.col("icd").is_null())
    pool = pl.concat([cases, controls], how="vertical")
    base = np.column_stack(
        [z(pool["age"].to_numpy()), z(pool["mean_beta"].to_numpy()), *(z(pool[c].to_numpy()) for c in fractions)]
    )
    rng = ctx.rng(380_010)
    n_boot = N_BOOT[str(ctx.scale.name)]
    rows = []
    for design, same_chip in (("same chip first", False), ("same chip only", True)):
        sets = match(cases, controls, same_chip_only=same_chip)
        shared = float(np.mean([pool["chip"][s[0]] == pool["chip"][j] for s in sets for j in s[1:]]))
        entries = [(c, [f"acc {c}"]) for c in [PRIMARY, *[k for k in clocks if k != PRIMARY]]]
        entries.append(("panel mean + disagreement", ["panel_mean", "disagreement"]))
        for name, cols in entries:
            result = models(pool, cols, sets, base, rng, n_boot)
            rows.append(
                {
                    "design": design,
                    "clock": name,
                    "generation": GENERATION.get(name, "panel" if name.startswith("panel") else "other"),
                    "role": "primary" if (name == PRIMARY and not same_chip) else "exploratory",
                    "sets": len(sets),
                    "controls_on_case_chip": shared,
                    **result,
                }
            )
        print(f"clock_increment: {design}: {len(sets)} sets, {len(entries)} models", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "nested_models")
    # best clock per generation and design: the strongest clock odds ratio away from 1 (exploratory)
    if "clock_or" in table.columns:
        best = (
            table.filter(pl.col("generation").str.contains(r"^\d"))
            .with_columns(strength=pl.col("clock_or").log().abs())
            .sort("strength", descending=True)
            .group_by("design", "generation", maintain_order=True)
            .first()
            .drop("strength")
            .sort("design", "generation")
        )
        write_table(best, out / "best_per_generation")
    corr = [
        {
            "pair": f"disagreement ~ {m}",
            "r": float(np.corrcoef(women["disagreement"].to_numpy(), women[m].to_numpy())[0, 1]),
        }
        for m in (B, F, P)
    ]
    corr.append({"pair": "disagreement panel: " + ", ".join(panel_clocks), "r": float(len(panel_clocks))})
    write_table(pl.DataFrame(corr), out / "disagreement_correlations")
    return []
