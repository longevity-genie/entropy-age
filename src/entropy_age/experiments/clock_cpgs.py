"""Where published clocks sit on the protection map (roadmap RM8, first part).

The results so far say that a CpG's drift follows its chromatin: CpGs in active promoter chromatin of a
tissue drift little there, CpGs in Polycomb, bivalent or quiet chromatin drift more. Published clocks pick
their CpGs by how well they track age, so they should come from the second kind, and a clock's behaviour in
another tissue should follow how its CpGs' chromatin differs between the tissues. For 16 clocks from biolearn
(``entropy-age fetch clocks``): chronological (Horvath, Hannum, Lin, the cortical clock), phenotypic and pace
(PhenoAge, DunedinPACE, telomere length, Zhang's mortality CpGs), the mitotic clock epiTOC1, the stochastic
clocks StocH, StocP, StocZ and the causal clocks CausAge, DamAge and AdaptAge:

* **chromatin** - the share of each clock's 450k CpGs in active TSS states (TssA, TssAFlnk) and in
  Polycomb or bivalent states, in blood (E062) and prefrontal cortex (E073), against background CpGs;
* **gene sets** - the share in promoters of core-essential, housekeeping and non-essential genes;
* **drift** - the clock CpGs' dispersion and entropy slopes against background CpGs of the same island
  relation and level (Hannum blood, Jaffe cortex), in SDs of the metric, with and without the tissue's
  chromatin added to the matching; the interval resamples the clock's CpGs;
* **tissue offset** - each clock's weighted sum in Lunnon's same donors, brain region minus blood, in SDs of
  the blood values, against the share of its CpGs whose state is active in one tissue and not the other.

Run as a job: ``entropy-age fetch clocks``, then ``entropy-age figures clocks``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.chromatin import ChromatinLayout, cpg_states, read_segmentation
from entropy_age.data.clocks import CLOCKS, ClockLayout, read_clock
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import (
    LUNNON_BRAIN,
    SETS,
    _hannum,
    _inputs_missing,
    _jaffe,
    _lifespan_metrics,
    _lunnon,
)
from entropy_age.protection import contrast_weights, resampled_gaps
from entropy_age.stochasticity import baseline_strata

BACKGROUND = 60_000
ACTIVE = ("1_TssA", "2_TssAFlnk")
POLYCOMB = ("10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk")
ISLAND_CLASS = {"Island": "island", "N_Shore": "shore", "S_Shore": "shore", "N_Shelf": "shelf", "S_Shelf": "shelf"}
N_RESAMPLES = {"paper": 200, "quick": 40}
TISSUE_STATE = {"blood (Hannum)": "E062", "cortex (Jaffe)": "E073"}


def clock_table(paths) -> pl.DataFrame:
    """``clock``, ``kind``, ``cpg``, ``weight`` for every clock CpG."""
    frames = []
    for clock, (_, _, kind) in CLOCKS.items():
        frame = read_clock(ClockLayout(paths).file(clock))
        frames.append(frame.with_columns(clock=pl.lit(clock), kind=pl.lit(kind)))
    return pl.concat(frames).select("clock", "kind", "cpg", "weight")


def chromatin_shares(states: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    """Share of the masked CpGs in active TSS states and in Polycomb or bivalent states."""
    n = max(int(mask.sum()), 1)
    return {
        "active": float(np.isin(states[mask], ACTIVE).sum() / n),
        "polycomb": float(np.isin(states[mask], POLYCOMB).sum() / n),
    }


def clock_gap(
    values: np.ndarray, in_clock: np.ndarray, background: np.ndarray, strata: np.ndarray, n_resamples: int, seed: int
) -> dict[str, float]:
    """Clock CpGs minus matched background CpGs in SDs of ``values``, with a CpG-resampling interval.

    Only clock and background CpGs enter (other clocks' CpGs are left out of both sides).
    """
    keep = in_clock | background
    v, s, c = values[keep], strata[keep], in_clock[keep]
    sd = float(v.std())
    p, _ = contrast_weights(c.astype(float), s)
    draws = resampled_gaps(
        v[None, :], [np.array([j]) for j in np.flatnonzero(c)], s, n_resamples, np.random.default_rng(seed)
    )
    return {
        "gap_sd": float(v @ p) / sd,
        "ci_low": float(np.quantile(draws[:, 0], 0.025)) / sd,
        "ci_high": float(np.quantile(draws[:, 0], 0.975)) / sd,
        "n_matched": int((p[c] > 0).sum()),
    }


def tissue_offsets(scores: dict[str, np.ndarray], donors: dict[str, np.ndarray]) -> list[dict]:
    """Per region: the mean over donors of (region - blood), in SDs of the blood values across donors."""
    blood, blood_donors = scores["whole blood"], donors["whole blood"]
    sd = float(blood.std(ddof=1))
    index = {d: i for i, d in enumerate(blood_donors)}
    out = []
    for region in LUNNON_BRAIN:
        pairs = [(scores[region][i], blood[index[d]]) for i, d in enumerate(donors[region]) if d in index]
        diff = np.array([a - b for a, b in pairs]) / sd
        out.append(
            {
                "region": region,
                "donors": diff.size,
                "offset_sd": float(diff.mean()),
                "offset_se": float(diff.std(ddof=1) / np.sqrt(diff.size)),
            }
        )
    return out


def run(ctx: Context) -> list[Path]:
    if _inputs_missing(ctx, "clocks"):
        return []
    missing = [str(ClockLayout(ctx.paths).file(c)) for c in CLOCKS if not ClockLayout(ctx.paths).file(c).exists()]
    if missing:
        print(
            f"clocks: run `entropy-age fetch clocks` first ({len(missing)} files missing); nothing run", file=sys.stderr
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("clock_cpgs")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    annotation = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(
        ~pl.col("chr").is_in(["X", "Y"]) & pl.col("pos").is_not_null()
    )
    clocks = clock_table(ctx.paths)
    on_array = set(probes["cpg"].to_list())
    coverage = clocks.group_by("clock", "kind", maintain_order=True).agg(
        pl.len().alias("cpgs"), pl.col("cpg").is_in(list(on_array)).sum().alias("on_450k_autosomes")
    )
    write_table(coverage, out / "coverage")
    clock_cpgs = sorted(set(clocks["cpg"].to_list()) & on_array)
    rest = sorted(on_array - set(clock_cpgs))
    rng = ctx.rng(140_000)
    background_cpgs = sorted(rng.choice(rest, size=min(BACKGROUND, len(rest)), replace=False).tolist())
    features = sorted(set(clock_cpgs) | set(background_cpgs))

    frame = pl.DataFrame({"cpg": features}).join(
        probes.select("cpg", "chr", "pos", "island_relation"), on="cpg", how="left"
    )
    layout = ChromatinLayout(ctx.paths)
    states = {
        e: cpg_states(frame, read_segmentation(layout.segmentation(e))).fill_null("none").to_numpy().astype(str)
        for e in ("E062", "E073")
    }
    membership = link_membership(pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.gene_sets))
    index = {f: j for j, f in enumerate(features)}
    is_background = np.isin(features, background_cpgs)
    in_clock = {
        c: np.isin(features, sel["cpg"].to_list()) for (c,), sel in clocks.group_by("clock", maintain_order=True)
    }

    rows = []
    for clock, mask in [*in_clock.items(), ("background", is_background)]:
        row = {"clock": clock, "kind": CLOCKS[clock][2] if clock in CLOCKS else "background", "cpgs": int(mask.sum())}
        for e, tissue in (("E062", "blood"), ("E073", "brain")):
            row |= {f"{k}_{tissue}": v for k, v in chromatin_shares(states[e], mask).items()}
        switch = np.isin(states["E062"], ACTIVE) != np.isin(states["E073"], ACTIVE)
        row["state_switch"] = float(switch[mask].mean()) if mask.any() else float("nan")
        row["active_blood_not_brain"] = float(
            (np.isin(states["E062"], ACTIVE) & ~np.isin(states["E073"], ACTIVE))[mask].mean()
        )
        row["active_brain_not_blood"] = float(
            (np.isin(states["E073"], ACTIVE) & ~np.isin(states["E062"], ACTIVE))[mask].mean()
        )
        for group in SETS:
            members = membership.filter(pl.col("group") == group)["feature"].to_list()
            row[f"share_{group}"] = float(np.isin(np.array(features)[mask], members).mean())
        rows.append(row)
    composition = pl.DataFrame(rows)
    write_table(composition, out / "composition")

    island = (
        frame["island_relation"]
        .fill_null("OpenSea")
        .replace_strict(ISLAND_CLASS, default="open sea")
        .to_numpy()
        .astype(str)
    )
    drift = []
    for label, loader in (("blood (Hannum)", _hannum), ("cortex (Jaffe)", _jaffe)):
        p = loader(ctx.paths, features)
        metrics = _lifespan_metrics(p)
        feats = metrics["feature"].to_list()
        cols = np.array([index[f] for f in feats])
        base = np.char.add(
            np.char.add(island[cols], "|"), baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
        )
        matched = {
            "island x level": base,
            "+ own state": np.char.add(np.char.add(base, "|"), states[TISSUE_STATE[label]][cols]),
        }
        for metric in ("dispersion_slope", "entropy_slope"):
            values = metrics[metric].to_numpy()
            for clock, mask in in_clock.items():
                for matching, strata in matched.items():
                    gap = clock_gap(values, mask[cols], is_background[cols], strata, n_resamples, ctx.seed)
                    drift.append({"tissue": label, "metric": metric, "clock": clock, "matching": matching, **gap})
        print(f"clocks: {label} done ({len(feats)} CpGs)", file=sys.stderr, flush=True)
    drift = pl.DataFrame(drift)
    write_table(drift, out / "drift")

    lunnon = _lunnon(ctx.paths, features)
    weights = {
        c: dict(sel.select("cpg", "weight").iter_rows()) for (c,), sel in clocks.group_by("clock", maintain_order=True)
    }
    offsets = []
    for clock, w in weights.items():
        scores, donors, used = {}, {}, 0
        for tissue, p in lunnon.items():
            cols = [j for j, f in enumerate(p.features) if f in w]
            used = len(cols)
            scores[tissue] = p.x[:, cols] @ np.array([w[p.features[j]] for j in cols])
            donors[tissue] = p.donors
        offsets += [{"clock": clock, "cpgs_used": used, **r} for r in tissue_offsets(scores, donors)]
    offsets = pl.DataFrame(offsets)
    write_table(offsets, out / "tissue_offsets")

    mean_offset = offsets.group_by("clock").agg(pl.col("offset_sd").mean().alias("offset_sd"))
    joined = composition.join(mean_offset, on="clock")
    rho = stats.spearmanr(joined["state_switch"], joined["offset_sd"].abs())
    summary = pl.DataFrame(
        [
            {
                "x": "state_switch",
                "y": "|mean brain - blood offset|",
                "clocks": joined.height,
                "spearman": float(rho.statistic),
                "p": float(rho.pvalue),
            }
        ]
    )
    write_table(summary, out / "offset_vs_chromatin")
    return [_figure(ctx, composition, drift, joined)]


def _figure(ctx: Context, composition: pl.DataFrame, drift: pl.DataFrame, joined: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.4)
    order = [*CLOCKS, "background"]

    ax = axes[0]
    comp = {r["clock"]: r for r in composition.iter_rows(named=True)}
    ys = np.arange(len(order))
    for j, (col, label) in enumerate(
        (
            ("active_blood", "active TSS, blood"),
            ("polycomb_blood", "Polycomb/bivalent, blood"),
            ("active_brain", "active TSS, cortex"),
            ("polycomb_brain", "Polycomb/bivalent, cortex"),
        )
    ):
        ax.barh(ys + (j - 1.5) * 0.2, [comp[c][col] for c in order], height=0.2, color=plots.PALETTE[j], label=label)
    ax.set_yticks(ys, order, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="Share of the clock's CpGs", title="Chromatin of clock CpGs (Roadmap E062, E073)")
    ax.legend(fontsize=6, loc="lower right")
    plots.label_panel(ax, "a")

    for i, tissue in enumerate(TISSUE_STATE):
        ax = axes[1 + i]
        sel = drift.filter((pl.col("tissue") == tissue) & (pl.col("metric") == "dispersion_slope"))
        for j, matching in enumerate(("island x level", "+ own state")):
            sub = {r["clock"]: r for r in sel.filter(pl.col("matching") == matching).iter_rows(named=True)}
            clocks = [c for c in CLOCKS if c in sub]
            y = np.arange(len(clocks)) + (j - 0.5) * 0.3
            g = np.array([sub[c]["gap_sd"] for c in clocks])
            ax.errorbar(
                g,
                y,
                xerr=[
                    g - np.array([sub[c]["ci_low"] for c in clocks]),
                    np.array([sub[c]["ci_high"] for c in clocks]) - g,
                ],
                fmt="os"[j],
                ms=3,
                color=plots.PALETTE[j + 4],
                label=matching,
            )
        ax.set_yticks(range(len(clocks)), clocks, fontsize=7)
        ax.invert_yaxis()
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.legend(fontsize=6)
        ax.set(xlabel="Clock CpGs minus matched CpGs (SD), 95 % CI", title=f"Dispersion slope, {tissue}")
        plots.label_panel(ax, "bc"[i])

    ax = axes[3]
    for k, kind in enumerate(sorted(joined["kind"].unique().to_list())):
        sub = joined.filter(pl.col("kind") == kind)
        ax.scatter(sub["state_switch"], sub["offset_sd"], s=24, color=plots.PALETTE[k % len(plots.PALETTE)], label=kind)
        for r in sub.iter_rows(named=True):
            ax.annotate(
                r["clock"], (r["state_switch"], r["offset_sd"]), fontsize=5, xytext=(2, 2), textcoords="offset points"
            )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set(
        xlabel="Share of clock CpGs whose TSS activity differs, blood vs cortex",
        ylabel="Brain − blood, same donors (SD of blood)",
        title="Lunnon: a clock's tissue offset",
    )
    ax.legend(fontsize=6)
    plots.label_panel(ax, "d")

    fig.suptitle(
        "Published clocks on the protection map (450k, autosomes)", x=0.02, ha="left", fontsize=11, fontweight="bold"
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("clock_cpgs"))
