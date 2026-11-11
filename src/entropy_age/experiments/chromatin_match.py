"""Does chromatin state account for the essential-gene protection? (roadmap RM7)

Every cohort test so far matches gene-set CpGs to background CpGs on CpG-island status and starting level.
If essential genes drift less because their promoters sit in active, open chromatin, matching on the
chromatin state as well should remove most of the gap; if something beyond chromatin protects them (for
instance selection among cells), the gap should survive.

Each 450k promoter CpG takes the state of a Roadmap Epigenomics 15-state chromHMM segmentation
(``entropy_age.data.chromatin``): peripheral-blood mononuclear cells (E062) for blood, dorsolateral
prefrontal cortex (E073) for brain, and blood neutrophils (E030) as a second blood reference. Drift
metrics are those of ``figures tissues``: in Hannum blood and Jaffe cortex the partial age slopes of
entropy and of dispersion; in Lunnon's same-donor tissues the between-person dispersion. For each set and
tissue the matched gap is computed with strata

* island x level (the cohort tests' matching),
* island x level x the tissue's own chromatin state,
* island x level x the other tissue's state (a control: a state that is not the tissue's own should remove
  less),

in SDs of the metric. Its interval comes from resampling the set's genes (the same draws for every
variant, so the share of the gap removed has a paired interval). The interval covers which genes are in
the set, not which people are in the cohort.

Run as a job after ``entropy-age fetch annotation`` (it adds the probe positions) and
``entropy-age fetch chromatin``: ``entropy-age figures chromatin``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.chromatin import EPIGENOMES, STATES, ChromatinLayout, cpg_states, read_segmentation
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import (
    LUNNON_BRAIN,
    SET_LABELS,
    SETS,
    _dispersion_metrics,
    _features,
    _hannum,
    _inputs_missing,
    _jaffe,
    _lifespan_metrics,
    _lunnon,
)
from entropy_age.protection import matched_scores, resampled_gaps
from entropy_age.stochasticity import baseline_strata

#: tissue -> (its own epigenome, the other tissue's epigenome)
OWN_AND_OTHER = {
    "blood (Hannum)": ("E062", "E073"),
    "cortex (Jaffe)": ("E073", "E062"),
    "whole blood": ("E062", "E073"),
    **{region: ("E073", "E062") for region in LUNNON_BRAIN},
}
N_RESAMPLES = {"paper": 200, "quick": 40}
NO_STATE = "none"


def strata_variants(island: np.ndarray, baseline: np.ndarray, states: dict[str, np.ndarray], tissue: str) -> dict:
    """The stratum label of every CpG under each matching: island x level, then with a chromatin state added."""
    base = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(baseline, 10).astype(str))
    own, other = OWN_AND_OTHER[tissue]
    variants = {"island x level": base}
    variants[f"+ own state ({own})"] = np.char.add(np.char.add(base, "|"), states[own])
    variants[f"+ other state ({other})"] = np.char.add(np.char.add(base, "|"), states[other])
    if own == "E062":
        variants["+ own state (E030)"] = np.char.add(np.char.add(base, "|"), states["E030"])
    return variants


def matched_gaps(
    values: np.ndarray, sets: dict, variants: dict[str, np.ndarray], n_resamples: int, seed: int
) -> list[dict]:
    """Per set and matching: the gap in SDs of ``values``, its gene-resampling draws, and the share removed."""
    sd = float(values.std())
    rows = []
    for group, (in_set, genes) in sets.items():
        draws_by_variant = {}
        for variant, strata in variants.items():
            gap = float(matched_scores(values[None, :], in_set, strata)[0][0]) / sd
            draws = resampled_gaps(values[None, :], genes, strata, n_resamples, np.random.default_rng(seed))[:, 0] / sd
            draws_by_variant[variant] = (gap, draws)
        base_gap, base_draws = draws_by_variant["island x level"]
        for variant, (gap, draws) in draws_by_variant.items():
            share = 1 - draws / base_draws
            rows.append(
                {
                    "group": group,
                    "matching": variant,
                    "gap_sd": gap,
                    "ci_low": float(np.quantile(draws, 0.025)),
                    "ci_high": float(np.quantile(draws, 0.975)),
                    "share_removed": 1 - gap / base_gap,
                    "share_ci_low": float(np.quantile(share, 0.025)),
                    "share_ci_high": float(np.quantile(share, 0.975)),
                    "n_set": int(in_set.sum()),
                }
            )
    return rows


def state_composition(states: np.ndarray, sets: dict) -> list[dict]:
    """Share of each set's CpGs, and of the CpGs outside every set, in each chromatin state."""
    member = np.any([s for s, _ in sets.values()], axis=0)
    rows = []
    for group, in_set in [*((g, s) for g, (s, _) in sets.items()), ("background", ~member)]:
        n = int(in_set.sum())
        for state in (*STATES, NO_STATE):
            rows.append({"group": group, "state": state, "share": float((states[in_set] == state).sum() / n)})
    return rows


def _states(features: list[str], probes: pl.DataFrame, paths) -> dict[str, np.ndarray]:
    if "pos" not in probes.columns:
        raise ValueError("the 450k probe table has no positions: rerun `entropy-age fetch annotation`")
    frame = pl.DataFrame({"cpg": features}).join(probes.select("cpg", "chr", "pos"), on="cpg", how="left")
    layout = ChromatinLayout(paths)
    out = {}
    for epigenome in EPIGENOMES:
        state = cpg_states(frame, read_segmentation(layout.segmentation(epigenome)))
        out[epigenome] = state.fill_null(NO_STATE).to_numpy().astype(str)
    return out


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "chromatin"):
        return []
    missing = [str(ChromatinLayout(ctx.paths).segmentation(e)) for e in EPIGENOMES]
    missing = [m for m in missing if not Path(m).exists()]
    if missing:
        print(
            f"chromatin: run `entropy-age fetch chromatin` first ({', '.join(missing)}); nothing was run",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("chromatin_match")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island = dict(probes.select("cpg", "cpg_island").iter_rows())

    tissues = {"blood (Hannum)": _hannum(ctx.paths, features), "cortex (Jaffe)": _jaffe(ctx.paths, features)}
    lunnon = _lunnon(ctx.paths, features)
    rows, composition = [], []
    for label, p, metric_fn, metric_names in [
        *((k, v, _lifespan_metrics, ("entropy_slope", "dispersion_slope")) for k, v in tissues.items()),
        *((k, v, _dispersion_metrics, ("dispersion",)) for k, v in lunnon.items()),
    ]:
        metrics = metric_fn(p)
        feats = metrics["feature"].to_list()
        states = _states(feats, probes, ctx.paths)
        sets = _sets(feats, membership)
        isl = np.array([bool(island.get(f, False)) for f in feats])
        variants = strata_variants(isl, metrics["baseline"].to_numpy(), states, label)
        for metric in metric_names:
            values = metrics[metric].to_numpy()
            finite = np.isfinite(values)
            if not finite.all():
                raise ValueError(f"{label}: {int((~finite).sum())} CpGs with no {metric}")
            rows += [
                {"tissue": label, "metric": metric, **r}
                for r in matched_gaps(values, sets, variants, n_resamples, seed=ctx.seed)
            ]
        own = OWN_AND_OTHER[label][0]
        composition += [{"tissue": label, "epigenome": own, **r} for r in state_composition(states[own], sets)]
        print(f"chromatin: {label} done ({len(feats)} CpGs)", file=sys.stderr, flush=True)

    gaps = pl.DataFrame(rows)
    composition = pl.DataFrame(composition)
    write_table(gaps, out / "matched_gaps")
    write_table(composition, out / "state_composition")
    return [_figure(ctx, gaps, composition)]


def _figure(ctx: Context, gaps: pl.DataFrame, composition: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.4)

    ax = axes[0]
    comp = composition.filter(pl.col("tissue") == "blood (Hannum)")
    groups = [*SETS, "background"]
    bottom = np.zeros(len(groups))
    shown = [s for s in STATES if comp.filter(pl.col("state") == s)["share"].max() >= 0.02]
    cmap = plots.plt.get_cmap("tab20")
    for k, state in enumerate([*shown, "other"]):
        if state == "other":
            share = np.array(
                [1 - comp.filter((pl.col("group") == g) & pl.col("state").is_in(shown))["share"].sum() for g in groups]
            )
        else:
            share = np.array(
                [comp.filter((pl.col("group") == g) & (pl.col("state") == state))["share"].sum() for g in groups]
            )
        ax.bar(range(len(groups)), share, bottom=bottom, color=cmap(k), label=state)
        bottom += share
    ax.set_xticks(range(len(groups)), [SET_LABELS.get(g, g) for g in groups], rotation=20, ha="right")
    ax.set(ylabel="Share of promoter CpGs", title="Chromatin state in blood (E062)", ylim=(0, 1))
    ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
    plots.label_panel(ax, "a")

    for i, metric in enumerate(("dispersion_slope", "dispersion")):
        ax = axes[1 + i]
        sel = gaps.filter((pl.col("metric") == metric) & pl.col("group").is_in(["housekeeping", "core_essential"]))
        tissues = sel["tissue"].unique(maintain_order=True).to_list()
        y, ticks, labels = 0.0, [], []
        for tissue in tissues:
            for group in ("housekeeping", "core_essential"):
                sub = sel.filter((pl.col("tissue") == tissue) & (pl.col("group") == group))
                for j, r in enumerate(sub.iter_rows(named=True)):
                    kind = r["matching"].split(" (")[0]
                    ax.errorbar(
                        r["gap_sd"],
                        y + 0.18 * j,
                        xerr=[[r["gap_sd"] - r["ci_low"]], [r["ci_high"] - r["gap_sd"]]],
                        fmt="o",
                        ms=3,
                        color=plots.PALETTE[["island x level", "+ own state", "+ other state"].index(kind)],
                        label=kind if (tissue, group, j) == (tissues[0], "housekeeping", j) else None,
                    )
                ticks.append(y + 0.25)
                labels.append(f"{SET_LABELS[group]} · {tissue}")
                y += 1
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(ticks, labels, fontsize=6)
        ax.invert_yaxis()
        ax.legend(fontsize=6)
        title = "Lifespan cohorts: dispersion slope" if metric == "dispersion_slope" else "Lunnon, same donors"
        ax.set(xlabel="Set minus matched CpGs (SD), 95 % gene-resampling CI", title=title)
        plots.label_panel(ax, "bc"[i])

    ax = axes[3]
    own = gaps.filter(
        pl.col("matching").str.starts_with("+ own state (E0")
        & ~pl.col("matching").str.contains("E030")
        & (pl.col("metric") != "entropy_slope")
    )
    for j, group in enumerate(("housekeeping", "core_essential")):
        sub = own.filter(pl.col("group") == group)
        xs = np.arange(sub.height) + (j - 0.5) * 0.35
        ax.bar(xs, sub["share_removed"], width=0.35, color=plots.PALETTE[j], label=SET_LABELS[group])
        ax.errorbar(
            xs,
            sub["share_removed"],
            yerr=[sub["share_removed"] - sub["share_ci_low"], sub["share_ci_high"] - sub["share_removed"]],
            fmt="none",
            ecolor="black",
            lw=0.8,
        )
    ticks = own.filter(pl.col("group") == "housekeeping")["tissue"].to_list()
    ax.set_xticks(range(len(ticks)), [t.replace(" (", "\n(") for t in ticks], fontsize=6, rotation=30, ha="right")
    ax.axhline(0.5, color=plots.MUTED, lw=0.8, ls=":")
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set(ylabel="Share of the gap removed by own chromatin", title="How much of the gap is chromatin?")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "d")

    fig.suptitle(
        "Does chromatin state account for the essential-gene protection? (450k promoters, Roadmap chromHMM)",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("chromatin_match"))
