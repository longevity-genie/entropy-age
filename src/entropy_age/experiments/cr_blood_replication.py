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
for ``S``.

**Randomisation p (audit 2026-11-18).** Each animal's score depends on which animals are controls (the reference,
the composition components and the leave-one-out control scores are all fitted on them), so shuffling labels on
fixed scores is not a test. ``p_perm`` reassigns treatment within each physical array chip, keeping every chip's
treated count, refits ``arm_scores`` for each reassignment and recomputes the same adjusted effect; all
assignments are enumerated when there are at most ``MAX_EXACT`` (CR against controls: 720), else ``MAX_EXACT``
are drawn, the observed one included. The observed assignment is checked to reproduce the reported effect.
Exchangeability is assumed only within a chip. Chips holding only one arm add no information; a comparison
whose arms share no chip (controls after 12 weeks against Time 0, all of Time 0 on one chip) has no
randomisation p and is flagged ``within_chip_overlap = 0``: its time contrast is confounded with chip.
``effect_chip`` is the OLS effect with chip intercepts where they are identifiable (descriptive).

``entropy-age figures cr_blood_replication`` (a small Slurm job).
"""

from __future__ import annotations

import gzip
import itertools
import math
import re
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
MAX_EXACT = {"paper": 1000, "quick": 50}


def read_groups(matrix: Path) -> pl.DataFrame:
    """``sample`` (the betas' column name, the title before ``__``), ``group`` (the title after it) and ``chip``
    (the physical array, from the deposited IDAT names; null if none)."""
    titles, chips = [], []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_supplementary_file") and not chips and any("_Grn" in c for c in cells):
                chips = [m.group(1) if (m := re.search(r"_(\d+)_R\d+C\d+_Grn", c)) else None for c in cells]
    parts = [t.split("__", 1) for t in titles]
    return pl.DataFrame(
        {
            "sample": [p[0] for p in parts],
            "group": [p[1] if len(p) > 1 else "" for p in parts],
            "chip": chips or [None] * len(parts),
        }
    )


def assignments(treated: np.ndarray, chips: np.ndarray, max_exact: int, rng) -> tuple[list[np.ndarray], int, bool]:
    """Reassignments of treatment that keep each chip's treated count: the observed one first, then all others
    when there are at most ``max_exact`` in total, else ``max_exact - 1`` random ones (repeats possible)."""
    per_chip = []
    for c in np.unique(chips):
        idx = np.flatnonzero(chips == c)
        per_chip.append((idx, int(treated[idx].sum())))
    total = int(np.prod([math.comb(len(idx), k) for idx, k in per_chip]))
    out = [treated.copy()]
    if total <= max_exact:
        options = [[np.array(c, int) for c in itertools.combinations(idx, k)] for idx, k in per_chip]
        for combo in itertools.product(*options):
            t = np.zeros(treated.size, bool)
            for chosen in combo:
                t[chosen] = True
            if not np.array_equal(t, treated):
                out.append(t)
        return out, total, True
    while len(out) < max_exact:
        t = np.zeros(treated.size, bool)
        for idx, k in per_chip:
            t[rng.choice(idx, size=k, replace=False)] = True
        out.append(t)
    return out, total, False


def overlap_chips(treated: np.ndarray, chips: np.ndarray) -> int:
    """Chips holding animals of both arms: the only ones a within-chip reassignment can move."""
    return int(sum(treated[chips == c].any() and (~treated[chips == c]).any() for c in np.unique(chips)))


def effect(y: np.ndarray, treated: np.ndarray, extra: np.ndarray | None = None) -> tuple[float, float, float]:
    cols = [np.ones(y.size), treated.astype(float)]
    if extra is not None:
        cols.append(extra)
    return ols(y, np.column_stack(cols), 1)


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
        chips = sel["chip"].fill_null("unknown").to_numpy()
        overlap = overlap_chips(treated, chips)
        chip_design = np.column_stack([(chips == c).astype(float) for c in np.unique(chips)[1:]])
        for k in COMPONENTS:
            scores = arm_scores(x, treated, isl, sets, components=k)
            per_animal.append(
                pl.DataFrame(
                    {"sample": sel["sample"], "group": sel["group"], "chip": sel["chip"], **scores}
                ).with_columns(comparison=pl.lit(name), components=pl.lit(k))
            )
            keys = [(g, m) for g in sets for m in ("B", "F", "P", "S")]
            observed = {(g, m): effect(scores[f"{m}_{g}"], treated, squeeze if m == "S" else None) for g, m in keys}
            null: dict[tuple[str, str], list[float]] = {key: [] for key in keys}
            draws, total, exact = [treated], 1, True
            if overlap:
                draws, total, exact = assignments(treated, chips, MAX_EXACT[str(ctx.scale.name)], rng)
                for t in draws:
                    sc = scores if t is draws[0] else arm_scores(x, t, isl, sets, components=k)
                    for g, m in keys:
                        null[(g, m)].append(effect(sc[f"{m}_{g}"], t, squeeze if m == "S" else None)[0])
                for key in keys:  # the observed assignment must reproduce the reported effect
                    assert abs(null[key][0] - observed[key][0]) < 1e-9, key
            for g, m in keys:
                eff, se, p = observed[(g, m)]
                y = scores[f"{m}_{g}"]
                extra = [squeeze] if m == "S" else []
                design = np.column_stack([np.ones(y.size), treated.astype(float), *extra, chip_design])
                identifiable = np.linalg.matrix_rank(design) == design.shape[1]
                v = np.array(null[(g, m)])
                rows.append(
                    {
                        "comparison": name,
                        "components": k,
                        "group": g,
                        "measure": m,
                        "effect": eff,
                        "se": se,
                        "p": p,
                        "p_perm": float((np.abs(v) >= abs(eff) - 1e-12).mean()) if overlap else None,
                        "assignments": len(draws) if overlap else 0,
                        "assignments_possible": total if overlap else 0,
                        "exact": exact if overlap else None,
                        "within_chip_overlap": overlap,
                        "chips": int(np.unique(chips).size),
                        "effect_chip": float(ols(y, design, 1)[0]) if identifiable else None,
                        "control_mean": float(y[~treated].mean()),
                        "treated": int(treated.sum()),
                        "controls": int((~treated).sum()),
                    }
                )
            print(
                f"cr_blood_replication: {name}, {k} components: {len(draws)} assignments", file=sys.stderr, flush=True
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
