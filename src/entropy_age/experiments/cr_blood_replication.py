"""A second caloric-restriction blood series: does CR weaken the shield again? (roadmap RM16 step 4, RM3)

In GSE80672 (``figures interventions``) CR weakened the shield's share ``F`` in B6D2F1 blood by about a quarter,
in the one comparison whose treated and control mice were sequenced together; CR in liver did not repeat it.
Selection reads this as fewer divisions purging fewer drifting cells; chromatin predicts no change of ``F``.

GSE282661 (Jayne et al. 2025, Nat Aging; ``data/input/geo/GSE282661``) holds C57BL/6J whole blood on the
mammalian array (HorvathMammalMethylChip320), all of one array run: 8 controls, 8 on CR, 8 on CR with a
torpor-like state (TLS+CR), 8 at the start (Time 0), and a separate long-term pair (TLS 7, its controls 6). The
CR is pair feeding at about a third of control intake four days a week for 12 weeks, at 32 °C, in mice of about
5-6 months: short and in young mice, so drift has little time to differ, and the paper saw no clock effect.

Every animal is scored on autosomal promoter probes (``consortium_interventions.arm_scores``): treated animals
against all controls of the comparison, each control against the other controls, with 0 and 3 composition
components of the controls' background CpGs removed. The effect is the OLS difference of treated minus control in
drift ``B``, the shield's share ``F``, the gap ``P`` and the signed shift ``S``, with each array's compression held
for ``S``; its p also from ``N_PERM`` label permutations.

``entropy-age figures cr_blood_replication`` (a small Slurm job).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.consortium_interventions import arm_scores, compression
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.tissue_contrast import SETS

GSE = "GSE282661"
BETAS = f"{GSE}_ET0301_Bvalues.csv.gz"
#: comparison -> (treated label, control label)
COMPARISONS = {
    "CR": ("CR", "Control"),
    "TLS+CR": ("TLS+CR", "Control"),
    "TLS, long-term": ("TLS (long-term)", "Control (long-term)"),
    "controls, 12 weeks vs start": ("Control", "Time 0"),
}
COMPONENTS = (0, 3)
N_PERM = {"paper": 2000, "quick": 200}


def read_groups(matrix: Path) -> pl.DataFrame:
    """``sample`` (the betas' column name, the title before ``__``) and ``group`` (the title after it)."""
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
                break
    parts = [t.split("__", 1) for t in titles]
    return pl.DataFrame({"sample": [p[0] for p in parts], "group": [p[1] if len(p) > 1 else "" for p in parts]})


def permutation_p(y: np.ndarray, treated: np.ndarray, n: int, rng: np.random.Generator) -> float:
    """Two-sided p of the treated-minus-control mean difference from ``n`` label shuffles."""
    obs = y[treated].mean() - y[~treated].mean()
    hits = 0
    for _ in range(n):
        t = rng.permutation(treated)
        hits += abs(y[t].mean() - y[~t].mean()) >= abs(obs) - 1e-12
    return (hits + 1) / (n + 1)


def run(ctx: Context) -> list[Path]:
    base = ctx.paths.input / "geo" / GSE
    if not (base / BETAS).exists():
        print(f"cr_blood_replication: {base / BETAS} missing; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("cr_blood_replication")
    groups = read_groups(base / f"{GSE}_series_matrix.txt.gz")
    lay = MouseAnnotationLayout(ctx.paths)
    probes = pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    header = pl.read_csv(base / BETAS, n_rows=0).columns
    frame = (
        pl.scan_csv(base / BETAS, null_values=["NA", ""])
        .rename({header[0]: "cpg"})
        .filter(pl.col("cpg").is_in(probes["cpg"].implode()))
        .collect()
        .drop_nulls()
    )
    features = frame["cpg"].to_list()
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    sets = _sets(features, members)
    groups = groups.filter(pl.col("sample").is_in(header))
    print(
        f"cr_blood_replication: {len(features)} promoter probes, {groups.height} samples", file=sys.stderr, flush=True
    )
    rng = ctx.rng(240_000)
    rows, per_animal = [], []
    for name, (treat, control) in COMPARISONS.items():
        sel = groups.filter(pl.col("group").is_in([treat, control]))
        treated = (sel["group"] == treat).to_numpy()
        if treated.sum() < 3 or (~treated).sum() < 4:
            continue
        x = frame.select(sel["sample"].to_list()).to_numpy().T.astype(float)
        squeeze = compression(x)
        for k in COMPONENTS:
            scores = arm_scores(x, treated, isl, sets, components=k)
            per_animal.append(
                pl.DataFrame({"sample": sel["sample"], "group": sel["group"], **scores}).with_columns(
                    comparison=pl.lit(name), components=pl.lit(k)
                )
            )
            for group in sets:
                for measure in ("B", "F", "P", "S"):
                    y = scores[f"{measure}_{group}"]
                    cols = [np.ones(y.size), treated.astype(float)]
                    if measure == "S":
                        cols.append(squeeze)
                    eff, se, p = ols(y, np.column_stack(cols), 1)
                    rows.append(
                        {
                            "comparison": name,
                            "components": k,
                            "group": group,
                            "measure": measure,
                            "effect": eff,
                            "se": se,
                            "p": p,
                            "p_perm": permutation_p(y, treated, N_PERM[str(ctx.scale.name)], rng),
                            "control_mean": float(y[~treated].mean()),
                            "treated": int(treated.sum()),
                            "controls": int((~treated).sum()),
                        }
                    )
    effects = pl.DataFrame(rows)
    write_table(effects, out / "effects")
    write_table(pl.concat(per_animal, how="diagonal"), out / "animal_scores")
    return [_figure(ctx, effects)]


def _figure(ctx: Context, effects: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for ax, measure, label in zip(axes, ("F", "B"), ("shield share F", "drift B"), strict=True):
        sub = effects.filter(
            (pl.col("measure") == measure) & (pl.col("components") == 3) & pl.col("group").is_in(list(SETS)[:2])
        ).sort("comparison", "group")
        y = np.arange(sub.height)
        ax.errorbar(sub["effect"], y, xerr=1.96 * sub["se"], fmt="o", color=plots.PALETTE[0], ms=4, lw=1)
        ax.set_yticks(y, [f"{c} · {g.replace('_', ' ')}" for c, g in sub.select("comparison", "group").iter_rows()])
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.set(xlabel=f"treated minus control, {label}", title=f"{label}, mouse blood (GSE282661)")
        ax.locator_params(axis="x", nbins=4)
        if measure == "B":
            ax.set_yticklabels([])
        plots.label_panel(ax, "ab"[measure == "B"])
    fig.subplots_adjust(wspace=0.15)
    return plots.save(fig, ctx.figure_path("cr_blood_replication"))
