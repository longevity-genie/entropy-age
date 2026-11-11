"""Clocks built from protected and from unprotected CpGs (roadmap RM8, second part).

If a CpG's drift is set by its chromatin, a clock's CpGs decide which tissues it can read. Three candidate
pools of autosomal 450k CpGs, defined by Roadmap chromHMM states in both blood (E062) and prefrontal cortex
(E073), are each subsampled to the same size:

* **active in both** - TSS states (TssA, TssAFlnk) in both tissues, the protected kind;
* **Polycomb in both** - bivalent or Polycomb-repressed states in both, the kind published clocks favour;
* **any** - CpGs drawn at random, a reference.

An elastic net is trained on each pool in one cohort (Hannum blood or Jaffe cortex) and applied to the other
and to Lunnon's same-donor tissues. Reported per pool, training cohort and test tissue: the correlation with
age, the median absolute error, and the mean offset (prediction minus age). In Lunnon the brain-minus-blood
offset of the same donors is free of lab and platform differences between cohorts, so it is the number that
says how tissue-specific a pool's clock is. Each pool is drawn ``N_DRAWS`` times for an interval.

Run as a job: ``entropy-age figures chromatin_clocks``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import ClockConfig, evaluate, train_clock
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.chromatin import ChromatinLayout, cpg_states, read_segmentation
from entropy_age.experiments.clock_cpgs import ACTIVE, POLYCOMB
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import LUNNON_BRAIN, _hannum, _inputs_missing, _jaffe, _lunnon

POOL_SIZE = {"paper": 4000, "quick": 800}
N_DRAWS = {"paper": 5, "quick": 2}
CANDIDATES = 150_000
POOLS = ("active in both", "Polycomb in both", "any")


def pools(features: list[str], blood: np.ndarray, brain: np.ndarray) -> dict[str, np.ndarray]:
    """Column indices of each candidate pool (``blood`` and ``brain`` are the CpGs' chromHMM states)."""
    active = np.isin(blood, ACTIVE) & np.isin(brain, ACTIVE)
    polycomb = np.isin(blood, POLYCOMB) & np.isin(brain, POLYCOMB)
    return {
        "active in both": np.flatnonzero(active),
        "Polycomb in both": np.flatnonzero(polycomb),
        "any": np.arange(len(features)),
    }


def lunnon_offsets(predictions: dict[str, np.ndarray], ages: dict[str, np.ndarray], donors: dict[str, np.ndarray]):
    """Per brain region: mean over shared donors of (region offset - blood offset), and its SE."""
    blood = dict(zip(donors["whole blood"], predictions["whole blood"] - ages["whole blood"], strict=True))
    out = []
    for region in LUNNON_BRAIN:
        diff = np.array(
            [
                p - a - blood[d]
                for p, a, d in zip(predictions[region], ages[region], donors[region], strict=True)
                if d in blood
            ]
        )
        out.append(
            {
                "region": region,
                "donors": diff.size,
                "brain_minus_blood": float(diff.mean()),
                "se": float(diff.std(ddof=1) / np.sqrt(diff.size)),
            }
        )
    return out


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "chromatin_clocks"):
        return []
    plots.apply_style()
    out = ctx.out_dir("chromatin_clocks")
    size, draws = POOL_SIZE[str(ctx.scale.name)], N_DRAWS[str(ctx.scale.name)]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    config = ClockConfig(l1_ratio=(0.1, 0.5, 0.9), n_alphas=30, cv=5, n_jobs=workers)

    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(
        ~pl.col("chr").is_in(["X", "Y"]) & pl.col("pos").is_not_null()
    )
    rng = ctx.rng(150_000)
    candidates = sorted(
        rng.choice(probes["cpg"].to_numpy(), size=min(CANDIDATES, probes.height), replace=False).tolist()
    )
    cohorts = {"blood (Hannum)": _hannum(ctx.paths, candidates), "cortex (Jaffe)": _jaffe(ctx.paths, candidates)}
    lunnon = _lunnon(ctx.paths, candidates)
    shared = set(candidates)
    for p in [*cohorts.values(), *lunnon.values()]:
        shared &= set(p.features)
    features = sorted(shared)
    frame = pl.DataFrame({"cpg": features}).join(probes.select("cpg", "chr", "pos"), on="cpg", how="left")
    layout = ChromatinLayout(ctx.paths)
    state = {
        e: cpg_states(frame, read_segmentation(layout.segmentation(e))).fill_null("none").to_numpy().astype(str)
        for e in ("E062", "E073")
    }
    pool = pools(features, state["E062"], state["E073"])
    sizes = {k: int(v.size) for k, v in pool.items()}
    size = min(size, *sizes.values())
    print(
        f"chromatin_clocks: {len(features)} shared CpGs, pools {sizes}, {size} per clock", file=sys.stderr, flush=True
    )

    def columns(p) -> np.ndarray:
        index = {f: j for j, f in enumerate(p.features)}
        return np.array([index[f] for f in features])

    cols = {k: columns(p) for k, p in [*cohorts.items(), *lunnon.items()]}
    accuracy, offsets = [], []
    for name in POOLS:
        for draw in range(draws):
            picked = np.sort(
                ctx.rng(151_000 + 100 * POOLS.index(name) + draw).choice(pool[name], size=size, replace=False)
            )
            for train_name, train in cohorts.items():
                clock = train_clock(train.x[:, cols[train_name][picked]], train.ages, config=config)
                tests = {k: v for k, v in [*cohorts.items(), *lunnon.items()]}
                predictions, ages, donors = {}, {}, {}
                for test_name, test in tests.items():
                    pred = clock.predict(test.x[:, cols[test_name][picked]])
                    m = evaluate(test.ages, pred)
                    accuracy.append(
                        {
                            "pool": name,
                            "draw": draw,
                            "trained_on": train_name,
                            "tested_on": test_name,
                            "in_sample": test_name == train_name,
                            "cpgs_used": clock.n_nonzero,
                            "pearson_r": m.pearson_r,
                            "median_abs_error": m.median_abs_error,
                            "mean_offset": float(np.mean(pred - test.ages)),
                        }
                    )
                    predictions[test_name], ages[test_name], donors[test_name] = pred, test.ages, test.donors
                lun = {k: predictions[k] for k in lunnon}
                offsets += [
                    {"pool": name, "draw": draw, "trained_on": train_name, **r}
                    for r in lunnon_offsets(lun, {k: ages[k] for k in lunnon}, {k: donors[k] for k in lunnon})
                ]
            print(f"chromatin_clocks: {name}, draw {draw} done", file=sys.stderr, flush=True)
    accuracy, offsets = pl.DataFrame(accuracy), pl.DataFrame(offsets)
    write_table(accuracy, out / "accuracy")
    write_table(offsets, out / "lunnon_offsets")
    write_table(
        pl.DataFrame([{"pool": k, "candidates": v, "per_clock": size} for k, v in sizes.items()]), out / "pools"
    )
    return [_figure(ctx, accuracy, offsets)]


def _figure(ctx: Context, accuracy: pl.DataFrame, offsets: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    tests = accuracy["tested_on"].unique(maintain_order=True).to_list()
    regions = list(LUNNON_BRAIN)
    for i, train in enumerate(("blood (Hannum)", "cortex (Jaffe)")):
        ax = axes[i]
        sel = accuracy.filter(pl.col("trained_on") == train)
        for j, name in enumerate(POOLS):
            agg = (
                sel.filter(pl.col("pool") == name)
                .group_by("tested_on")
                .agg(pl.col("pearson_r").mean(), pl.col("pearson_r").std().alias("sd"))
            )
            agg = {r["tested_on"]: r for r in agg.iter_rows(named=True)}
            xs = np.arange(len(tests)) + (j - 1) * 0.27
            ax.bar(
                xs,
                [agg[t]["pearson_r"] for t in tests],
                yerr=[agg[t]["sd"] or 0 for t in tests],
                width=0.27,
                color=plots.PALETTE[j],
                label=name,
            )
        ax.set_xticks(range(len(tests)), [t.replace(" (", "\n(") for t in tests], fontsize=6, rotation=30, ha="right")
        ax.set(ylabel="r with age", ylim=(0, 1), title=f"Trained on {train}")
        ax.legend(fontsize=6, loc="lower left")
        plots.label_panel(ax, "ab"[i])

        ax = axes[2 + i]
        sel = offsets.filter(pl.col("trained_on") == train)
        for j, name in enumerate(POOLS):
            agg = (
                sel.filter(pl.col("pool") == name)
                .group_by("region")
                .agg(pl.col("brain_minus_blood").mean(), pl.col("brain_minus_blood").std().alias("sd"))
            )
            agg = {r["region"]: r for r in agg.iter_rows(named=True)}
            xs = np.arange(len(regions)) + (j - 1) * 0.27
            ax.bar(
                xs,
                [agg[r]["brain_minus_blood"] for r in regions],
                yerr=[agg[r]["sd"] or 0 for r in regions],
                width=0.27,
                color=plots.PALETTE[j],
                label=name,
            )
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set_xticks(range(len(regions)), [r.replace(" ", "\n", 1) for r in regions], fontsize=6)
        ax.set(ylabel="Brain − blood age, same donors (years)", title=f"Lunnon, clocks trained on {train}")
        ax.legend(fontsize=6)
        plots.label_panel(ax, "cd"[i])
    fig.suptitle(
        "Clocks from protected (active) and unprotected (Polycomb) CpGs, equal pool sizes",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("chromatin_clocks"))
