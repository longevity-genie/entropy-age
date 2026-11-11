"""Interventions across tissues of the same mice: which move the shield? (roadmap RM12, RM3)

The Mammalian Methylation Consortium (GSE223748, fetched for ``figures species``) includes mouse intervention
studies: caloric restriction in liver, GHRKO dwarf mice with wild-type littermates in brain, kidney and liver,
partial reprogramming (OSKM with doxycycline) in kidney and skin, and sleep disruption in cerebellum, heart and
liver. For each study and tissue, every animal is scored against the controls of that study and tissue
(``entropy_age.protection``): drift ``B`` and the shield as a share of it, ``F``. Controls are scored against the
other controls (leave one out), so treated and control animals are scored alike. The effect of the intervention
is the OLS coefficient of ``treated`` on ``B`` or ``F`` with age as a covariate where ages vary.

Selection predicts that an intervention which slows proliferation weakens the shield in renewing tissue and
leaves non-dividing tissue (brain) alone; chromatin predicts a change, if any, alike across tissues. The arms
are small (3-59 animals); each effect comes with its standard error and its animals, and groups under
``MIN_ARM`` per arm are left out.

Controls come two ways (column ``controls``). ``pooled``: every control of the study and tissue, which in the
"Genetic perturbation" study pools the wild types of several knockouts. ``same_chip``: only the controls
hybridised on the same array chips as the treated animals (the chip is the first part of the barcode), with the
chip as a covariate; on a chip shared by a knockout and its wild types these are its littermate controls, and a
batch effect cannot separate them.

Run as a job after ``entropy-age fetch mammals``: ``entropy-age figures mmc_interventions``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols, probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.per_person import score_against
from entropy_age.experiments.species_shield import layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.protection import fit_reference

CONTROLS = {"Control", "WT", "1.Control", "Normal"}
MIN_ARM = 5


def same_chip(chips: np.ndarray, treated: np.ndarray, is_control: np.ndarray) -> np.ndarray:
    """Controls on a chip that also holds a treated animal, and the treated animals."""
    shared = np.isin(chips, chips[treated])
    return treated | (is_control & shared)


def chip_design(chips: np.ndarray) -> np.ndarray:
    """Dummy columns for every chip but the first (the intercept carries it); none for a single chip."""
    levels = sorted(set(chips.tolist()))
    if len(levels) == 1:
        return np.empty((chips.size, 0))
    return np.column_stack([chips == c for c in levels[1:]]).astype(float)


def arm_scores(x: np.ndarray, treated: np.ndarray, island: np.ndarray, sets: dict) -> dict[str, np.ndarray]:
    """Scores of every animal: treated against all controls, each control against the other controls."""
    member = np.any([s for s, _ in sets.values()], axis=0)
    out: dict[str, np.ndarray] = {}
    controls = np.flatnonzero(~treated)
    rows = {}
    ref = fit_reference(x[controls], island, ~member)
    for i, sc in zip(np.flatnonzero(treated), _rows(score_against(x[treated], ref, sets)), strict=True):
        rows[i] = sc
    for c in controls:
        others = controls[controls != c]
        sc = score_against(x[c : c + 1], fit_reference(x[others], island, ~member), sets)
        rows[int(c)] = {k: float(v[0]) for k, v in sc.items()}
    for k in rows[controls[0]]:
        out[k] = np.array([rows[i][k] for i in range(x.shape[0])])
    return out


def _rows(scores: dict[str, np.ndarray]) -> list[dict[str, float]]:
    n = len(next(iter(scores.values())))
    return [{k: float(v[i]) for k, v in scores.items()} for i in range(n)]


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print("mmc_interventions: run `entropy-age fetch mammals` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("mmc_interventions")
    samples = read_samples(matrix).filter(
        (pl.col("species") == "Mouse") & pl.col("experiment").is_not_null() & pl.col("barcode").is_not_null()
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    samples = samples.filter(pl.col("barcode").is_in(header))
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([header[0], *samples["barcode"].unique().to_list()])
        .filter(pl.col(header[0]).is_in(probes["cpg"].implode()))
        .collect()
        .rename({header[0]: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    sets = _sets(features, members)
    rows, scored = [], []
    for (experiment, tissue), sel in samples.group_by("experiment", "tissue", maintain_order=True):
        is_control = sel["intervention"].is_in(list(CONTROLS)).to_numpy()
        if is_control.sum() < MIN_ARM:
            continue
        x_all = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x_all).any(axis=0)
        x_all = x_all[:, ok]
        sub_sets = {
            g: (s[ok], [np.searchsorted(np.flatnonzero(ok), gi[ok[gi]]) for gi in genes if ok[gi].any()])
            for g, (s, genes) in sets.items()
        }
        chips_all = sel["barcode"].str.extract(r"^(\d+)_").to_numpy()
        for arm, controls in (
            (a, c) for a in sorted(set(sel["intervention"].to_list()) - CONTROLS) for c in ("pooled", "same_chip")
        ):
            is_arm = (sel["intervention"] == arm).to_numpy()
            keep = is_control | is_arm
            if controls == "same_chip":
                keep = same_chip(chips_all, is_arm, is_control)
            treated = is_arm[keep]
            if treated.sum() < MIN_ARM or (~treated).sum() < MIN_ARM:
                print(
                    f"mmc_interventions: {experiment}, {tissue}, {arm}, {controls}: {int(treated.sum())} treated, "
                    f"{int((~treated).sum())} controls; left out",
                    file=sys.stderr,
                )
                continue
            x = x_all[keep]
            scores = arm_scores(x, treated, isl[ok], sub_sets)
            ages = sel["age"].to_numpy()[keep]
            design = [np.ones(treated.size), treated.astype(float)]
            if np.ptp(ages) > 0:
                design.append(ages)
            design = np.column_stack(design)
            if controls == "same_chip":
                design = np.column_stack([design, chip_design(chips_all[keep])])
            scored.append(
                pl.DataFrame(
                    {
                        "sample": sel["sample"].to_numpy()[keep],
                        "experiment": experiment,
                        "tissue": tissue,
                        "arm": arm,
                        "controls": controls,
                        "treated": treated,
                        "age": ages,
                        **scores,
                    }
                )
            )
            for group in sub_sets:
                for measure in ("B", "F"):
                    y = scores[f"{measure}_{group}"]
                    coef, se, p = ols(y, design, 1)
                    rows.append(
                        {
                            "experiment": experiment,
                            "tissue": tissue,
                            "arm": arm,
                            "controls": controls,
                            "group": group,
                            "measure": measure,
                            "effect": coef,
                            "se": se,
                            "p": p,
                            "control_mean": float(y[~treated].mean()),
                            "treated_n": int(treated.sum()),
                            "control_n": int((~treated).sum()),
                        }
                    )
            print(f"mmc_interventions: {experiment}, {tissue}, {arm}, {controls}: done", file=sys.stderr, flush=True)
    effects = pl.DataFrame(rows)
    write_table(effects, out / "effects")
    write_table(pl.concat(scored, how="diagonal"), out / "animal_scores")
    return [_figure(ctx, effects)]


def _figure(ctx: Context, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=5.0)
    sel = effects.filter(pl.col("group") == "housekeeping")
    labels = [
        f"{r['experiment']} · {r['arm']} · {r['tissue']}" + (" · same chip" if r["controls"] == "same_chip" else "")
        for r in sel.filter(pl.col("measure") == "F").iter_rows(named=True)
    ]
    for i, measure in enumerate(("B", "F")):
        ax = axes[i]
        sub = sel.filter(pl.col("measure") == measure)
        rel = (sub["effect"] / sub["control_mean"].abs()).to_numpy()
        err = (1.96 * sub["se"] / sub["control_mean"].abs()).to_numpy()
        ys = np.arange(sub.height)
        ax.errorbar(rel, ys, xerr=err, fmt="o", ms=4, color=plots.PALETTE[i + 2])
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(ys, labels if i == 0 else [""] * len(labels), fontsize=6)
        ax.invert_yaxis()
        title = {"B": "Drift B (> 0: more drift)", "F": "Shield share F (> 0: weaker shield)"}[measure]
        ax.set(xlabel="Effect relative to the controls' |mean| (95 % CI)", title=title)
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Mouse interventions in the Mammalian Methylation Consortium (housekeeping genes)",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("mmc_interventions"))
