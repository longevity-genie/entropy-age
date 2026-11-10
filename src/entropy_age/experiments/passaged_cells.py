"""The direct selection test: does the protection gap grow with cell divisions? (roadmap RM5)

Endicott et al. (GSE197512) passaged seven primary human cell lines to replicative senescence and
measured each culture on EPIC at many population doublings (PD); in three lines they also arrested
growth with mitomycin C (MMC) next to dividing controls, and in one they withdrew serum, lowered oxygen
or immortalised with TERT. Divisions and time in culture therefore come apart.

**The two readings.** If essential genes are protected because their promoters are maintained better
(chromatin), each division adds less error there than at matched CpGs, but in the same proportion at
every PD: the *fraction* by which set CpGs deviate less than matched CpGs stays constant as deviation
accumulates. If they are protected because cells that damage an essential gene are outcompeted
(selection), the purge accumulates with divisions: the fraction grows with PD. In MMC-arrested cells
there are no divisions, so selection predicts no change at all, whatever the time in culture.

**The per-sample score** (``docs/PER_PERSON_SCORE.md``). In each experiment (one cell line in one
sub-experiment) the reference is the mean of its samples at the earliest day. For every later sample
and promoter CpG, ``d = |x - m| / sqrt(m (1 - m))`` with ``m`` the reference level. A gene set's
matched gap ``P`` is the stratum-weighted mean of ``d`` over set CpGs minus that over all other CpGs
of the same stratum (island status x reference decile), with weights the set's share per stratum,
exactly the contrast ``compare_groups`` takes over features. ``B`` is the matched background
deviation itself and ``F = P / B`` the relative gap: -0.2 means set CpGs have moved 20 % less than
matched CpGs. Each arm's ``P``, ``B`` and ``F`` are regressed on PD (per 10 PD) and on days (per 10 days).

**Caveat.** Loss of methylation in culture is concentrated in late-replicating, CpG-poor domains;
promoter CpGs replicate early, set CpGs perhaps earlier than matched ones. Matching on island status
and level removes part of this, not all; replication timing is not controlled.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS, _features, _matrix
from entropy_age.stochasticity import baseline_strata

ARM_COLUMNS = ("oxygen", "agent", "serum", "vector")
COURSES = ("Baseline profiling", "Reduced oxygen", "TERT immortalization")


def deviations(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """``|x - m| / sqrt(m (1 - m))`` per sample and CpG, ``m`` clipped to [0.01, 0.99]."""
    m = np.clip(reference, 0.01, 0.99)
    return np.abs(x - reference) / np.sqrt(m * (1 - m))


def noise_floor(x: np.ndarray, groups: list[np.ndarray]) -> np.ndarray:
    """Per-CpG variance of one array: the within-group variance of replicate samples, pooled over groups.

    Replicates are parallel cultures measured on the same day; their spread is array noise plus what
    the cultures diverged since they were split, the floor every deviation from a reference sits on.
    """
    usable = [g for g in groups if g.size >= 2]
    if not usable:
        raise ValueError("no group holds two replicates")
    total = sum(x[g].var(axis=0, ddof=1) * (g.size - 1) for g in usable)
    return total / sum(g.size - 1 for g in usable)


def _replicate_groups(unit: pl.DataFrame) -> tuple[list[np.ndarray], float | None]:
    """Rows of the replicate samples (same day, same arm) at the earliest day that has any."""
    groups = unit.group_by("days", "arm").agg(pl.col("row")).filter(pl.col("row").list.len() >= 2)
    if groups.height == 0:
        return [], None
    day = groups["days"].min()
    return [np.array(r) for r in groups.filter(pl.col("days") == day)["row"].to_list()], float(day)


def matched_scores(d: np.ndarray, in_group: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per sample: the set-minus-matched gap ``P`` and the matched background deviation ``B``.

    Strata holding only set CpGs or only other CpGs are dropped, as in ``compare_groups``.
    """
    codes = np.unique(strata, return_inverse=True)[1]
    k = int(codes.max()) + 1
    n_in = np.bincount(codes, weights=in_group.astype(float), minlength=k)
    n_out = np.bincount(codes, weights=(~in_group).astype(float), minlength=k)
    usable = (n_in > 0) & (n_out > 0)
    if not usable.any():
        raise ValueError("no stratum holds both set and other CpGs")
    onehot = np.zeros((codes.size, k))
    onehot[np.arange(codes.size), codes] = 1.0
    mean_in = (d[:, in_group] @ onehot[in_group])[:, usable] / n_in[usable]
    mean_out = (d[:, ~in_group] @ onehot[~in_group])[:, usable] / n_out[usable]
    weights = n_in[usable] / n_in[usable].sum()
    return (mean_in - mean_out) @ weights, mean_out @ weights


def _arm(row: dict) -> str:
    parts = [f"{c}={row[c]}" for c in ARM_COLUMNS if row[c] not in (None, "NA", "Ambient", "Full")]
    return ", ".join(parts) or "standard"


def score_samples(
    samples: pl.DataFrame, x: np.ndarray, features: list[str], probes: pl.DataFrame, membership: pl.DataFrame
) -> pl.DataFrame:
    """One row per non-reference sample: experiment, arm, PD, days and the scores of every set.

    ``P``, ``B`` and ``F`` use absolute deviations. ``Q``, ``C`` and ``G`` are the same contrasts on
    squared deviations with the noise floor removed: ``(x - m)^2 - v (1 + 1/k)``, scaled by
    ``m (1 - m)``, with ``v`` the per-CpG replicate variance (``noise_floor``) and ``k`` the reference
    samples. Array noise is the same for set and matched CpGs, so in ``F`` it pulls the fraction toward
    zero wherever deviations are small; ``G`` is free of that floor. An experiment without replicates
    uses the floor pooled over the others (``noise_source``).
    """
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    islands = np.array([bool(island.get(f, False)) for f in features])
    sets = {g: np.isin(features, membership.filter(pl.col("group") == g)["feature"].to_list()) for g in SETS}
    sets = {g: s for g, s in sets.items() if s.sum() >= 2}
    samples = samples.with_row_index("row").with_columns(
        arm=pl.struct(*ARM_COLUMNS).map_elements(_arm, return_dtype=pl.String)
    )
    units = list(samples.group_by("subexperiment", "donor", maintain_order=True))
    floors: dict[tuple, tuple[np.ndarray, str]] = {}
    for key, unit in units:
        groups, day = _replicate_groups(unit)
        if groups:
            floors[key] = (noise_floor(x, groups), f"own replicates, day {day:g}")
    if not floors:
        raise ValueError("no experiment has replicate samples to estimate the noise floor from")
    pooled = np.mean([v for v, _ in floors.values()], axis=0)
    rows = []
    for key, unit in units:
        subexperiment, donor = key
        first = unit["days"].min()
        ref = unit.filter(pl.col("days") == first)
        later = unit.filter(pl.col("days") > first)
        if later.height == 0:
            continue
        reference = x[ref["row"].to_numpy()].mean(axis=0)
        strata = np.char.add(islands.astype(str), baseline_strata(reference, 10).astype(str))
        d = deviations(x[later["row"].to_numpy()], reference)
        floor, source = floors.get(key, (pooled, "pooled over experiments"))
        m = np.clip(reference, 0.01, 0.99)
        sq = ((x[later["row"].to_numpy()] - reference) ** 2 - floor * (1 + 1 / ref.height)) / (m * (1 - m))
        out = later.select("sample", "arm", "doublings", "days").with_columns(
            noise_source=pl.lit(source),
            subexperiment=pl.lit(subexperiment),
            donor=pl.lit(donor),
            doublings_gained=pl.col("doublings") - float(ref["doublings"].mean()),
            reference_samples=pl.lit(ref.height),
        )
        scores = {}
        for group, in_group in sets.items():
            p, b = matched_scores(d, in_group, strata)
            q, c = matched_scores(sq, in_group, strata)
            scores |= {f"P_{group}": p, f"B_{group}": b, f"F_{group}": p / b}
            scores |= {f"Q_{group}": q, f"C_{group}": c, f"G_{group}": q / c}
        rows.append(out.with_columns(**{k: pl.Series(v) for k, v in scores.items()}))
    return pl.concat(rows) if rows else pl.DataFrame()


def _groups(scores: pl.DataFrame) -> list[str]:
    return [g for g in SETS if f"P_{g}" in scores.columns]


def arm_slopes(scores: pl.DataFrame) -> pl.DataFrame:
    """OLS slope of ``P``, ``B`` and ``F`` on PD gained (per 10 PD) and on days (per 10 days), per arm."""
    rows = []
    for (subexperiment, donor, arm), sel in scores.group_by("subexperiment", "donor", "arm", maintain_order=True):
        for axis, column in (("doublings", "doublings_gained"), ("days", "days")):
            t = sel[column].to_numpy() / 10.0
            if sel.height < 3 or np.ptp(t) == 0:
                continue
            for group in _groups(scores):
                for measure in ("P", "B", "F", "Q", "C", "G"):
                    fit = stats.linregress(t, sel[f"{measure}_{group}"].to_numpy())
                    rows.append(
                        {
                            "subexperiment": subexperiment,
                            "donor": donor,
                            "arm": arm,
                            "axis": axis,
                            "group": group,
                            "measure": measure,
                            "samples": sel.height,
                            "slope_per_10": float(fit.slope),
                            "se": float(fit.stderr),
                            "p": float(fit.pvalue),
                            "start": float(fit.intercept),
                            "range": float(np.ptp(t) * 10),
                        }
                    )
    return pl.DataFrame(rows)


def across_lines(slopes: pl.DataFrame) -> pl.DataFrame:
    """The passaged lines (``Baseline profiling``) as replicates: mean slope, t interval, signs."""
    rows = []
    base = slopes.filter((pl.col("subexperiment") == "Baseline profiling") & (pl.col("axis") == "doublings"))
    for (group, measure), sel in base.group_by("group", "measure", maintain_order=True):
        v = sel["slope_per_10"].to_numpy()
        if v.size < 3:
            continue
        mean, se = float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))
        half = float(stats.t.ppf(0.975, v.size - 1) * se)
        rows.append(
            {
                "group": group,
                "measure": measure,
                "lines": v.size,
                "mean_slope_per_10_pd": mean,
                "ci_low": mean - half,
                "ci_high": mean + half,
                "negative": int((v < 0).sum()),
                "p_t": float(stats.ttest_1samp(v, 0.0).pvalue),
            }
        )
    return pl.DataFrame(rows).sort("measure", "group") if rows else pl.DataFrame()


def drift_coupling(scores: pl.DataFrame) -> pl.DataFrame:
    """Per arm: ``Q ~ 1 + C + C * (PD gained / 10)``, the gap against the noise-free matched drift.

    ``coupling`` is the gap per unit of matched drift at the reference (``-c`` under chromatin
    protection), ``change_per_10_pd`` how that coupling moves with every 10 doublings: zero when the
    protection is a fixed fraction of drift, negative when it grows with divisions, as selection
    predicts. Fitting ``Q`` on ``C`` avoids dividing by drifts near zero, where ``G = Q / C`` explodes.
    """
    rows = []
    for (subexperiment, donor, arm), sel in scores.group_by("subexperiment", "donor", "arm", maintain_order=True):
        if sel.height < 5:
            continue
        t = sel["doublings_gained"].to_numpy() / 10.0
        for group in _groups(scores):
            c, q = sel[f"C_{group}"].to_numpy(), sel[f"Q_{group}"].to_numpy()
            design = np.column_stack([np.ones_like(c), c, c * t])
            coef, *_ = np.linalg.lstsq(design, q, rcond=None)
            dof = sel.height - 3
            residual = q - design @ coef
            cov = np.linalg.pinv(design.T @ design) * float(residual @ residual) / max(dof, 1)
            se = np.sqrt(np.clip(np.diag(cov), 0, None))
            rows.append(
                {
                    "subexperiment": subexperiment,
                    "donor": donor,
                    "arm": arm,
                    "group": group,
                    "samples": sel.height,
                    "pd_range": float(np.ptp(t) * 10),
                    "coupling": float(coef[1]),
                    "coupling_se": float(se[1]),
                    "change_per_10_pd": float(coef[2]),
                    "change_se": float(se[2]),
                    "p_change": float(2 * stats.t.sf(abs(coef[2] / se[2]), dof)) if se[2] > 0 else float("nan"),
                }
            )
    return pl.DataFrame(rows)


def coupling_across_lines(coupling: pl.DataFrame) -> pl.DataFrame:
    """The passaged lines as replicates: mean change of the coupling per 10 PD, t interval, signs."""
    rows = []
    base = coupling.filter(pl.col("subexperiment") == "Baseline profiling")
    for group, sel in base.group_by("group", maintain_order=True):
        v = sel["change_per_10_pd"].to_numpy()
        if v.size < 3:
            continue
        mean, se = float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))
        half = float(stats.t.ppf(0.975, v.size - 1) * se)
        rows.append(
            {
                "group": group[0] if isinstance(group, tuple) else group,
                "lines": v.size,
                "mean_coupling": float(sel["coupling"].mean()),
                "mean_change_per_10_pd": mean,
                "ci_low": mean - half,
                "ci_high": mean + half,
                "negative": int((v < 0).sum()),
                "p_t": float(stats.ttest_1samp(v, 0.0).pvalue),
            }
        )
    return pl.DataFrame(rows).sort("group") if rows else pl.DataFrame()


def arrest_contrast(scores: pl.DataFrame) -> pl.DataFrame:
    """MMC-arrested minus dividing controls at the last day, per line: ``B`` and ``F`` of each set."""
    mmc = scores.filter(pl.col("subexperiment") == "Growth arrest: MMC")
    rows = []
    for donor, sel in mmc.group_by("donor", maintain_order=True):
        last = sel["days"].max()
        at = sel.filter(pl.col("days") == last)
        arrested, control = at.filter(pl.col("arm") == "agent=MMC"), at.filter(pl.col("arm") == "agent=control")
        if arrested.height < 2 or control.height < 2:
            continue
        for group in _groups(scores):
            row = {"donor": donor[0] if isinstance(donor, tuple) else donor, "day": last, "group": group}
            row |= {
                "pd_gained_mmc": float(arrested["doublings_gained"].mean()),
                "pd_gained_control": float(control["doublings_gained"].mean()),
            }
            for measure in ("B", "F", "C", "G"):
                a, c = arrested[f"{measure}_{group}"].to_numpy(), control[f"{measure}_{group}"].to_numpy()
                row |= {
                    f"{measure}_mmc": float(a.mean()),
                    f"{measure}_control": float(c.mean()),
                    f"{measure}_p_welch": float(stats.ttest_ind(a, c, equal_var=False).pvalue),
                }
            rows.append(row)
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    lay = cohort_layout("endicott", ctx.paths)
    needed = [lay.betas, AnnotationLayout(ctx.paths).links_450k]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        print(
            f"passage: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation` and "
            "`entropy-age fetch cohort endicott` first. Nothing was run.",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("passaged_cells")
    features, probes, membership = _features(ctx.paths, ctx.rng(160_000))
    samples = pl.read_parquet(lay.samples).sort("subexperiment", "donor", "days", "sample")
    m = _matrix("endicott", ctx.paths, samples["sample"].to_list(), features)
    membership = membership.filter(pl.col("feature").is_in(m.features))
    scores = score_samples(samples, m.values, m.features, probes, membership)
    write_table(scores, out / "sample_scores")
    slopes = arm_slopes(scores)
    write_table(slopes, out / "arm_slopes")
    lines = across_lines(slopes)
    write_table(lines, out / "across_lines")
    coupling = drift_coupling(scores)
    write_table(coupling, out / "drift_coupling")
    write_table(coupling_across_lines(coupling), out / "coupling_across_lines")
    arrest = arrest_contrast(scores)
    write_table(arrest, out / "mmc_arrest")
    write_table(
        membership.group_by("group").agg(pl.len().alias("cpgs"), pl.col("gene").n_unique().alias("genes")),
        out / "set_sizes",
    )
    return [_figure(ctx, scores, slopes, arrest)]


def _figure(ctx: Context, scores: pl.DataFrame, slopes: pl.DataFrame, arrest: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=5.2)
    base = scores.filter(pl.col("subexperiment") == "Baseline profiling")
    donors = base["donor"].unique(maintain_order=True).to_list()
    colors = {d: plots.PALETTE[i % len(plots.PALETTE)] for i, d in enumerate(donors)}
    for ax, measure, title in (
        (axes[0], "B_housekeeping", "Deviation from the earliest sample (matched CpGs)"),
        (axes[1], "G_housekeeping", "Housekeeping gap as a fraction of it, noise floor removed"),
    ):
        for donor in donors:
            sel = base.filter(pl.col("donor") == donor).sort("doublings_gained")
            ax.plot(
                sel["doublings_gained"],
                sel[measure],
                "o-",
                ms=3,
                lw=0.8,
                color=colors[donor],
                label=donor.replace("_", " "),
            )
        ax.set(xlabel="Population doublings since the reference", title=title)
        if measure.startswith("G"):
            ax.axhline(0, color=plots.MUTED, lw=0.8)
    axes[0].set_ylabel("B (binomial-scaled |Δβ|)")
    axes[1].set_ylabel("G = Q / C (squared, floor removed)")
    axes[0].legend(fontsize=6, loc="upper left")
    plots.label_panel(axes[0], "a")
    plots.label_panel(axes[1], "b")

    ax = axes[2]
    sel = slopes.filter((pl.col("axis") == "doublings") & (pl.col("measure") == "G")).sort("subexperiment", "donor")
    keys = sel.select("subexperiment", "donor", "arm").unique(maintain_order=True).rows()
    for j, group in enumerate(SETS):
        g = sel.filter(pl.col("group") == group)
        for i, key in enumerate(keys):
            r = g.filter((pl.col("subexperiment") == key[0]) & (pl.col("donor") == key[1]) & (pl.col("arm") == key[2]))
            if r.height == 0:
                continue
            y = i + (j - 1) * 0.25
            v, se = r["slope_per_10"][0], r["se"][0]
            ax.plot([v - 1.96 * se, v + 1.96 * se], [y, y], color=plots.PALETTE[j], lw=0.9)
            ax.scatter([v], [y], color=plots.PALETTE[j], s=10, label=SET_LABELS[group] if i == 0 else None)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(
        range(len(keys)),
        [f"{k[1].replace('_', ' ')} · {k[0].split(':')[-1].strip()} · {k[2]}" for k in keys],
        fontsize=5.5,
    )
    ax.invert_yaxis()
    ax.set(xlabel="Slope of G per 10 PD (< 0 = gap grows with divisions)", title="Every arm, on divisions")
    ax.legend(fontsize=6, loc="lower right")
    plots.label_panel(ax, "c")

    ax = axes[3]
    if arrest.height:
        hk = arrest.filter(pl.col("group") == "housekeeping")
        for i, r in enumerate(hk.iter_rows(named=True)):
            ax.scatter([i - 0.12], [r["G_mmc"]], color=plots.ACCENT, s=22, label="MMC arrested" if i == 0 else None)
            ax.scatter([i + 0.12], [r["G_control"]], color=plots.PALETTE[0], s=22, label="dividing" if i == 0 else None)
        ax.set_xticks(
            range(hk.height),
            [d.replace("_", " ") for d in hk["donor"].to_list()],
            rotation=20,
            ha="right",
            fontsize=6.5,
        )
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.legend(fontsize=7)
    ax.set(ylabel="G, housekeeping, last day", title="Arrested against dividing, same line and day")
    plots.label_panel(ax, "d")
    fig.suptitle(
        "Does the protection gap grow with divisions? Primary cells to senescence (GSE197512)",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("passaged_cells"))
