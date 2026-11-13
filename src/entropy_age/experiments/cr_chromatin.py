"""Does caloric restriction's effect on the blood shield survive matching on chromatin? (roadmap RM12)

In mouse whole blood (GSE80672, ``figures interventions``) CR weakens the shield's share ``F`` by about a
quarter in the one within-batch comparison (B6D2F1), robust to composition. Two readings:

* **chromatin**: CR changes how open essential promoters are, or works through what open chromatin does; then
  matching every set CpG to background CpGs of the same openness takes most of the CR effect away;
* **selection**: CR slows the divisions that purge cells drifting at essential genes; then the effect lies
  beyond chromatin and survives the matching.

Openness is read from ImmGen ATAC-seq (Yoshida et al. 2019, GSE100738, mm10) of six cell types that make up
mouse blood: B cells, naive CD4 and CD8 T cells, NK cells, Ly6C-high monocytes and granulocytes. Each RRBS
promoter CpG gets the mean signal over the six tracks in ``±HALF_WINDOW`` bp. The strata of the
per-mouse scores (``entropy_age.protection``) are then island × young level decile (as before) × openness class
(``ACCESS_BINS`` quantile classes of the signal), and every comparison of ``figures interventions`` is refitted.
Also refitted: the age gap of the standard C57BL/6 mice (entropy and dispersion slopes, set minus matched).

``entropy-age figures cr_chromatin`` as a Slurm job, after ``entropy-age fetch rrbs`` and the six bigWigs in
``data/input/geo/GSE100738``.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.rrbs import RrbsLayout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.interventions import (
    COMPARISONS,
    COMPONENTS,
    CR_WINDOW,
    REFERENCE_MONTHS,
    _select,
    _sets,
    effect,
    membership,
)
from entropy_age.experiments.per_person import score_against
from entropy_age.protection import fit_reference
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics

GSE = "GSE100738"
TRACKS = ("B.Sp", "T.4.Nve.Sp", "T.8.Nve.Sp", "NK.27+11b+.Sp", "Mo.6C+II-.Bl", "GN.Sp")
HALF_WINDOW = 250
ACCESS_BINS = 4


def openness(chrom: np.ndarray, pos: np.ndarray, tracks: list[Path], half: int = HALF_WINDOW) -> np.ndarray:
    """Per CpG, the mean over ``tracks`` of each bigWig's mean signal in ``[pos - half, pos + half]``."""
    import pyBigWig

    total = np.zeros(pos.size)
    for path in tracks:
        bw = pyBigWig.open(str(path))
        sizes = bw.chroms()
        vals = np.zeros(pos.size)
        for c in np.unique(chrom):
            name = c if str(c).startswith("chr") else f"chr{c}"
            if name not in sizes:
                continue
            rows = np.flatnonzero(chrom == c)
            for i in rows:
                lo, hi = max(0, int(pos[i]) - half), min(sizes[name], int(pos[i]) + half)
                v = bw.stats(name, lo, hi, type="mean")[0] if hi > lo else None
                vals[i] = v if v is not None else 0.0
        bw.close()
        total += vals
    return total / len(tracks)


def access_class(signal: np.ndarray, bins: int = ACCESS_BINS) -> np.ndarray:
    """Quantile class (0 = least open) of each CpG's signal; ties at zero all go to class 0."""
    edges = np.quantile(signal, np.linspace(0, 1, bins + 1)[1:-1])
    return np.searchsorted(edges, signal, side="left")


def with_access(strata: np.ndarray, access: np.ndarray) -> np.ndarray:
    return np.char.add(np.char.add(strata.astype(str), "a"), access.astype(str))


def _aging(x, ages, features, island, access, members, rows) -> list[dict]:
    xs, t = x[rows], ages[rows]
    metrics = drift_metrics(xs, t, features, spearman=False)
    young = (t >= REFERENCE_MONTHS[0]) & (t <= REFERENCE_MONTHS[1])
    level = xs[young].mean(axis=0) if young.sum() >= 3 else metrics["baseline"].to_numpy()
    base = np.char.add(island.astype(str), baseline_strata(level, 10).astype(str))
    out = []
    for matching, strata in (("island x level", base), ("+ openness", with_access(base, access))):
        m = metrics.with_columns(stratum=pl.Series(strata))
        for metric in ("entropy_slope", "dispersion_slope"):
            table = compare_groups(m, members, metric=metric, stratify_by="stratum", n_permutations=0)
            sd = float(m[metric].std())
            out += [
                {"matching": matching, "metric": metric, "group": r["group"], "gap_sd": r["matched_diff"] / sd}
                for r in table.to_dicts()
            ]
    return out


def run(ctx: Context) -> list[Path]:
    lay = RrbsLayout(ctx.paths)
    tracks = [ctx.paths.input / "geo" / GSE / f"{GSE}_{t}.bw" for t in TRACKS]
    missing = [p.name for p in tracks if not p.exists()]
    if not lay.betas.exists() or missing:
        print(f"cr_chromatin: RRBS betas or bigWigs missing ({missing}); nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("cr_chromatin")
    samples = pl.read_parquet(lay.samples).drop_nulls("age_months")
    frame = pl.read_parquet(lay.betas)
    m = FeatureMatrix.from_frame(frame, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.1)
    features = m.features
    cpgs = pl.read_parquet(lay.cpgs).filter(pl.col("cpg").is_in(features))
    cpgs = pl.DataFrame({"cpg": features}).join(cpgs, on="cpg", how="left", maintain_order="left")
    isl = cpgs["cpg_island"].fill_null(False).to_numpy()
    members = membership(ctx.paths, features, pl.read_parquet(lay.cpgs))
    sets = _sets(features, members)
    print(f"cr_chromatin: {len(features)} promoter CpGs, {len(tracks)} ATAC tracks", file=sys.stderr, flush=True)
    signal = openness(cpgs["chr"].to_numpy().astype(str), cpgs["pos"].to_numpy(), tracks)
    access = access_class(signal)
    member = np.any([s for s, _ in sets.values()], axis=0)
    write_table(
        pl.DataFrame(
            [
                {"group": g, "access_class": int(c), "share": float((access[s] == c).mean()), "cpgs": int(s.sum())}
                for g, (s, _) in {**sets, "background": (~member, [])}.items()
                for c in range(ACCESS_BINS)
            ]
        ),
        out / "openness_by_set",
    )

    x = m.values
    ages = samples["age_months"].to_numpy()
    wt = np.flatnonzero(
        (
            (samples["strain"] == "C57BL/6") & (samples["diet"] == "standard") & (samples["condition"] == "Wild Type")
        ).to_numpy()
    )
    aging = pl.DataFrame(_aging(x, ages, features, isl, access, members, wt))
    write_table(aging, out / "aging_gaps")

    ref = wt[(ages[wt] >= REFERENCE_MONTHS[0]) & (ages[wt] <= REFERENCE_MONTHS[1])]
    scored = np.setdiff1d(np.arange(x.shape[0]), ref)
    rows = []
    for k in COMPONENTS:
        base = fit_reference(x[ref], isl, ~member, components=k)
        for matching, reference in (
            ("island x level", base),
            ("+ openness", dataclasses.replace(base, strata=with_access(base.strata, access))),
        ):
            scores = pl.DataFrame(
                {"sample": samples["sample"].to_numpy()[scored], **score_against(x[scored], reference, sets)}
            ).join(samples, on="sample")
            for name, (strain, treat, control, adjust) in COMPARISONS.items():
                sel = scores.filter(_select(scores, strain, treat) | _select(scores, strain, control))
                if name == "CR, C57BL/6":
                    sel = sel.filter(pl.col("age_months").is_between(*CR_WINDOW))
                treated = sel.select(_select(sel, strain, treat).alias("t"))["t"].to_numpy()
                if treated.sum() < 3 or (~treated).sum() < 3:
                    continue
                for group in sets:
                    for measure in ("B", "F", "P"):
                        rows.append(
                            {
                                "comparison": name,
                                "components": k,
                                "matching": matching,
                                "group": group,
                                "measure": measure,
                                **effect(sel, f"{measure}_{group}", treated, adjust),
                            }
                        )
    effects = pl.DataFrame(rows)
    write_table(effects, out / "intervention_effects")
    return [_figure(ctx, aging, effects)]


def _figure(ctx: Context, aging: pl.DataFrame, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    groups = ["core_essential", "housekeeping"]
    ax = axes[0]
    sub = effects.filter((pl.col("measure") == "F") & (pl.col("components") == 5) & pl.col("group").is_in(groups)).sort(
        "comparison", "group"
    )
    labels = sub.filter(pl.col("matching") == "island x level").select("comparison", "group")
    y = np.arange(labels.height)
    for i, matching in enumerate(("island x level", "+ openness")):
        t = labels.join(sub.filter(pl.col("matching") == matching), on=["comparison", "group"], how="left")
        ax.errorbar(
            t["effect"],
            y + 0.18 * (i - 0.5),
            xerr=1.96 * t["se"],
            fmt="o",
            color=plots.PALETTE[i],
            ms=4,
            lw=1,
            label=f"matched on {matching}",
        )
    ax.set_yticks(y, [f"{c} · {g.replace('_', ' ')}" for c, g in labels.iter_rows()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="effect on shield share F (5 components removed)", title="Interventions, mouse blood")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
    plots.label_panel(ax, "a")
    ax = axes[1]
    a = aging.filter(pl.col("group").is_in(groups)).sort("metric", "group")
    keys = a.filter(pl.col("matching") == "island x level").select("metric", "group")
    y = np.arange(keys.height)
    for i, matching in enumerate(("island x level", "+ openness")):
        t = keys.join(a.filter(pl.col("matching") == matching), on=["metric", "group"], how="left")
        ax.barh(y + 0.38 * i, t["gap_sd"], height=0.36, color=plots.PALETTE[i])
    ax.set_yticks(y + 0.19, [f"{m.replace('_slope', '')} · {g.replace('_', ' ')}" for m, g in keys.iter_rows()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="age gap, set minus matched (SD)", title="Aging, standard C57BL/6")
    plots.label_panel(ax, "b")
    fig.subplots_adjust(wspace=1.1)
    return plots.save(fig, ctx.figure_path("cr_chromatin"))
