"""Interventions that slow aging: do they change the shield's share? (roadmap RM3)

Mouse whole-blood RRBS of Petkovich et al. 2017 (GSE80672, ``entropy-age fetch rrbs``): about 150 C57BL/6
mice on a standard diet aged 0.7-35 months, caloric restriction (CR) in C57BL/6 (10-27 months) and
B6D2F1 (20-27 months), and two long-lived dwarf models with their littermate controls at 5-6 months
(GHRKO, Snell). Two readings make opposite predictions for an intervention that slows proliferation:

* **chromatin**: total drift may slow, but the shield stays the same share of it;
* **selection**: fewer divisions mean less purging, so the shield's share shrinks (``F`` moves toward 0).

Two parts, both on promoter CpGs (TSS -1500 to +500), each gene set against background CpGs of the same
CpG-island status and young-reference level decile:

1. **Aging in standard C57BL/6** - the age slopes of entropy and dispersion (``drift_metrics``, age in
   months), set minus matched CpGs in SD of the metric, with a bootstrap over mice: does mouse blood
   repeat the human blood shield, with exact ages?
2. **Per-mouse scores** (``entropy_age.protection``) against the young standard C57BL/6 (3-6 months):
   the matched background deviation ``B`` (how far a mouse has drifted) and the relative gap ``F``. Each
   intervention is a term in an OLS within its strain: ``B`` or ``F`` ~ intervention (+ age where the
   ages vary). Comparisons are only within a strain, since strain differences (and their SNPs in RRBS
   reads) shift every CpG. Every score is also computed with five composition components of the young
   reference's background CpGs removed (CR shifts the blood cell mix), and each comparison says whether
   treated and control mice share a sample series (``same_series``, from the titles' prefixes): the
   C57BL/6 CR mice were sequenced as a series of their own, the B6D2F1 CR mice with their controls.

Run as a job: ``entropy-age fetch rrbs`` (downloads and parses about 3 GB), then
``entropy-age figures interventions``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.mouse import MouseAnnotationLayout
from entropy_age.data.rrbs import TRACKS, RrbsLayout, cpg_genes, promoters
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import score_against
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.protection import fit_reference, gene_counts
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics

REFERENCE_MONTHS = (3.0, 6.0)
N_BOOT = {"paper": 200, "quick": 30}
#: per-mouse scores without and with five composition components fitted on the young reference
COMPONENTS = (0, 5)
#: comparison -> (strain, treated rows, control rows, age adjusted)
COMPARISONS = {
    "CR, C57BL/6": ("C57BL/6", {"diet": "Calorie Restricted"}, {"diet": "standard"}, True),
    "CR, B6D2F1": ("B6D2F1", {"diet": "Calorie Restricted"}, {"diet": "standard"}, True),
    "GHRKO": ("(C57BL/6J x BALB/cByJ)/F2", {"condition": "GHRKO"}, {"condition": "GHRKO WT"}, False),
    "Snell dwarf": ("(DW/J x C3H/HEJ)/F2", {"condition": "Snell Dwarf"}, {"condition": "Snell Dwarf Control"}, True),
}
#: age window of the C57BL/6 CR comparison (CR mice are 10-27 months)
CR_WINDOW = (10.0, 28.0)


def membership(paths, features: list[str], cpgs: pl.DataFrame) -> pl.DataFrame:
    lay = RrbsLayout(paths)
    prom = promoters(lay.track(TRACKS[0]), lay.track(TRACKS[1]))
    links = cpg_genes(cpgs.filter(pl.col("cpg").is_in(features)), prom)
    sets = pl.read_parquet(MouseAnnotationLayout(paths).gene_sets).filter(pl.col("gene_set").is_in(SETS))
    return (
        links.join(sets.select(pl.col("gene_set").alias("group"), "entrez"), on="entrez")
        .select("feature", "group", pl.col("entrez").cast(pl.String).alias("gene"))
        .unique(subset=["feature", "group"])
    )


def _sets(features: list[str], members: pl.DataFrame) -> dict:
    out = {}
    for group in SETS:
        rows = members.filter(pl.col("group") == group)
        by_gene = {g: sel["feature"].to_list() for (g,), sel in rows.group_by("gene")}
        in_set, genes = gene_counts(features, by_gene)
        if in_set.sum() >= 2:
            out[group] = (in_set, genes)
    return out


def aging_gaps(x: np.ndarray, ages: np.ndarray, features: list[str], island: np.ndarray, members, rows=None):
    """Set minus matched CpGs (SD) for the entropy and dispersion slopes, strata from the young mice."""
    rows = np.arange(x.shape[0]) if rows is None else rows
    xs, t = x[rows], ages[rows]
    metrics = drift_metrics(xs, t, features, spearman=False)
    young = (t >= REFERENCE_MONTHS[0]) & (t <= REFERENCE_MONTHS[1])
    level = xs[young].mean(axis=0) if young.sum() >= 3 else metrics["baseline"].to_numpy()
    strata = np.char.add(island.astype(str), baseline_strata(level, 10).astype(str))
    metrics = metrics.with_columns(stratum=pl.Series(strata))
    out = []
    for metric in ("entropy_slope", "dispersion_slope"):
        table = compare_groups(metrics, members, metric=metric, stratify_by="stratum", n_permutations=0)
        sd = float(metrics[metric].std())
        out += [{"metric": metric, "group": r["group"], "gap_sd": r["matched_diff"] / sd} for r in table.to_dicts()]
    return out


def effect(frame: pl.DataFrame, column: str, treated: np.ndarray, age_adjusted: bool) -> dict:
    """OLS coefficient of the treatment on ``column`` (with age in months when ``age_adjusted``)."""
    y = frame[column].to_numpy()
    cols = [np.ones(y.size), treated.astype(float)]
    if age_adjusted and np.ptp(frame["age_months"].to_numpy()) > 0:
        cols.append(frame["age_months"].to_numpy())
    design = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    dof = y.size - design.shape[1]
    se = float(np.sqrt(resid @ resid / dof * np.linalg.inv(design.T @ design)[1, 1]))
    control_mean = float(y[~treated].mean())
    return {
        "effect": float(coef[1]),
        "se": se,
        "p": float(2 * stats.t.sf(abs(coef[1] / se), dof)),
        "control_mean": control_mean,
        "relative": float(coef[1] / control_mean) if control_mean else float("nan"),
        "treated": int(treated.sum()),
        "controls": int((~treated).sum()),
    }


def _select(samples: pl.DataFrame, strain: str, spec: dict) -> pl.Expr:
    expr = pl.col("strain") == strain
    for k, v in spec.items():
        expr &= pl.col(k) == v
    return expr


def run(ctx: Context) -> list[Path]:
    lay = RrbsLayout(ctx.paths)
    if not lay.betas.exists():
        print("interventions: run `entropy-age fetch rrbs` first; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("interventions")
    samples = pl.read_parquet(lay.samples).drop_nulls("age_months")
    frame = pl.read_parquet(lay.betas)
    m = FeatureMatrix.from_frame(frame, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.1)
    features = m.features
    cpgs = pl.read_parquet(lay.cpgs)
    island = dict(cpgs.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island[f]) for f in features])
    members = membership(ctx.paths, features, cpgs)
    write_table(
        members.group_by("group").agg(pl.len().alias("cpgs"), pl.col("gene").n_unique().alias("genes")),
        out / "set_sizes",
    )
    x = m.values
    ages = samples["age_months"].to_numpy()
    print(f"interventions: {len(features)} promoter CpGs, {x.shape[0]} blood samples", file=sys.stderr, flush=True)

    # part 1: aging in standard C57BL/6
    wt = np.flatnonzero(
        (
            (samples["strain"] == "C57BL/6") & (samples["diet"] == "standard") & (samples["condition"] == "Wild Type")
        ).to_numpy()
    )
    point = aging_gaps(x, ages, features, isl, members, wt)
    rng = ctx.rng(160_000)
    boot = []
    for b in range(N_BOOT[str(ctx.scale.name)]):
        rows = np.sort(rng.choice(wt, size=wt.size, replace=True))
        boot += [{"boot": b, **r} for r in aging_gaps(x, ages, features, isl, members, rows)]
    boot = pl.DataFrame(boot)
    aging = (
        pl.DataFrame(point)
        .join(
            boot.group_by("metric", "group").agg(
                pl.col("gap_sd").quantile(0.025).alias("ci_low"), pl.col("gap_sd").quantile(0.975).alias("ci_high")
            ),
            on=["metric", "group"],
        )
        .with_columns(mice=pl.lit(wt.size))
    )
    write_table(aging, out / "aging_gaps")

    # part 2: per-mouse scores
    sets = _sets(features, members)
    member = np.any([s for s, _ in sets.values()], axis=0)
    ref = wt[(ages[wt] >= REFERENCE_MONTHS[0]) & (ages[wt] <= REFERENCE_MONTHS[1])]
    scored = np.setdiff1d(np.arange(x.shape[0]), ref)
    samples = samples.with_columns(series=samples["title"].str.extract(r"^(\d*[A-Za-z]+)").str.to_uppercase())
    all_scores, rows = [], []
    for k in COMPONENTS:
        reference = fit_reference(x[ref], isl, ~member, components=k)
        scores = pl.DataFrame(
            {"sample": samples["sample"].to_numpy()[scored], **score_against(x[scored], reference, sets)}
        ).join(samples, on="sample")
        all_scores.append(scores.with_columns(components=pl.lit(k)))
        for name, (strain, treat, control, adjust) in COMPARISONS.items():
            sel = scores.filter(_select(scores, strain, treat) | _select(scores, strain, control))
            if name == "CR, C57BL/6":
                sel = sel.filter(pl.col("age_months").is_between(*CR_WINDOW))
            treated = sel.select(_select(sel, strain, treat).alias("t"))["t"].to_numpy()
            if treated.sum() < 3 or (~treated).sum() < 3:
                print(
                    f"interventions: {name}: {int(treated.sum())} treated, {int((~treated).sum())} controls; skipped",
                    file=sys.stderr,
                )
                continue
            same = set(sel.filter(pl.Series(treated))["series"]) == set(sel.filter(pl.Series(~treated))["series"])
            for group in sets:
                for measure in ("B", "F", "P"):
                    est = effect(sel, f"{measure}_{group}", treated, adjust)
                    rows.append(
                        {
                            "comparison": name,
                            "components": k,
                            "same_series": same,
                            "group": group,
                            "measure": measure,
                            **est,
                        }
                    )
    scores = pl.concat(all_scores)
    write_table(scores, out / "mouse_scores")
    effects = pl.DataFrame(rows)
    write_table(effects, out / "intervention_effects")
    return [_figure(ctx, aging, scores, effects)]


def _figure(ctx: Context, aging: pl.DataFrame, scores: pl.DataFrame, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    ax = axes[0]
    for j, metric in enumerate(("entropy_slope", "dispersion_slope")):
        sub = {r["group"]: r for r in aging.filter(pl.col("metric") == metric).iter_rows(named=True)}
        groups = [g for g in SETS if g in sub]
        xs = np.arange(len(groups)) + (j - 0.5) * 0.3
        g = np.array([sub[k]["gap_sd"] for k in groups])
        ax.errorbar(
            xs,
            g,
            yerr=[g - np.array([sub[k]["ci_low"] for k in groups]), np.array([sub[k]["ci_high"] for k in groups]) - g],
            fmt="os"[j],
            color=plots.PALETTE[j],
            label=metric.replace("_", " "),
        )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set_xticks(range(len(groups)), [SET_LABELS[k] for k in groups])
    ax.set(ylabel="Set minus matched CpGs (SD), 95 % bootstrap CI", title="Aging, standard C57BL/6 blood")
    ax.legend(fontsize=7)
    plots.label_panel(ax, "a")

    for i, measure in enumerate(("B", "F")):
        ax = axes[1 + i]
        b6 = scores.filter((pl.col("strain") == "C57BL/6") & (pl.col("components") == 0))
        for j, diet in enumerate(("standard", "Calorie Restricted")):
            sub = b6.filter(pl.col("diet") == diet)
            ax.scatter(
                sub["age_months"], sub[f"{measure}_housekeeping"], s=10, color=plots.PALETTE[j], label=diet, alpha=0.8
            )
        ax.set(
            xlabel="Age (months)",
            ylabel={"B": "B: matched drift from young", "F": "F: shield as a share of drift"}[measure],
            title=f"C57BL/6, housekeeping: {measure}",
        )
        ax.legend(fontsize=7)
        plots.label_panel(ax, "bc"[i])

    ax = axes[3]
    sel = effects.filter((pl.col("group") == "housekeeping") & pl.col("measure").is_in(["B", "F"]))
    comps = sel["comparison"].unique(maintain_order=True).to_list()
    for j, (measure, k) in enumerate((("B", 0), ("F", 0), ("F", 5))):
        sub = {
            r["comparison"]: r
            for r in sel.filter((pl.col("measure") == measure) & (pl.col("components") == k)).iter_rows(named=True)
        }
        ys = np.arange(len(comps)) + (j - 1) * 0.25
        rel = np.array([sub[c]["effect"] / abs(sub[c]["control_mean"]) for c in comps])
        err = np.array([1.96 * sub[c]["se"] / abs(sub[c]["control_mean"]) for c in comps])
        ax.errorbar(
            rel,
            ys,
            xerr=err,
            fmt="osD"[j],
            color=plots.PALETTE[j + 2],
            label={"B": "drift B", "F": "shield share F"}[measure] + (f", {k} components" if k else ""),
        )
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(len(comps)), comps, fontsize=7)
    ax.invert_yaxis()
    ax.set(
        xlabel="Intervention effect, relative to the controls' |mean| (95 % CI)",
        title="Housekeeping: what interventions change",
    )
    ax.legend(fontsize=7)
    plots.label_panel(ax, "d")
    fig.suptitle(
        "Interventions in mouse blood RRBS (GSE80672): drift and the shield's share",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("interventions"))
