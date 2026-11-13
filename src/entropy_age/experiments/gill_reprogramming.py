"""Transient reprogramming of human fibroblasts, from IDATs: does a reset spare essential promoters? (roadmap RM17)

Gill et al. 2022 (eLife; GSE165179, EPIC) took dermal fibroblasts of three middle-aged donors (38-53 years) through
maturation-phase transient reprogramming (OSKM for 10-17 days, then recovery). Some cells transiently reprogrammed
(SSEA4-positive); others went through the same treatment and failed to; untreated cells are the negative controls.
The authors report clocks about 30 years younger in the reprogrammed cells.

The series has no young donors, so the shield's restoration cannot be read against youth. What it can say is
whether a large, rejuvenating methylation change respects the shield: chromatin maintenance predicts that
essential promoters move less than their look-alikes under reprogramming, as they do under aging.

Betas come from the IDATs (``scripts/reprocess_idats.py``, SeSAMe-like; run R65a) in
``data/interim/cohorts/gill_betas_sesame.parquet``. Per donor, the negative controls are the reference
(``entropy_age.protection``; each control scored against the donor's other controls), and every sample gets drift
``B``, the gap ``P``, the shield share ``F``, the entropy gap ``E`` and the signed shift ``S`` for each gene set, on
450k promoter CpGs (gene sets and look-alikes of ``figures tissues``). Contrasts, within donor: transiently
reprogrammed against failed (reprogramming itself, the treatment held), and failed against controls (the
treatment without reprogramming); p from permuting labels within donors. Horvath's skin and blood clock (raw
weighted sum, higher = older) checks that the rejuvenation is seen in these betas.

``entropy-age figures gill_reprogramming`` as a Slurm job, after R65a.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _islands, _sets, score_against
from entropy_age.experiments.shield_clocks import linear_score
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import fit_reference

GSE = "GSE165179"
BETAS = "cohorts/gill_betas_sesame.parquet"
ARMS = {"Transiently reprogrammed": "reprogrammed", "Failed to transiently": "failed", "Negative control": "control"}
CONTRASTS = {
    "reprogramming (reprogrammed − failed)": ("reprogrammed", "failed"),
    "treatment (failed − control)": ("failed", "control"),
}
MEASURES = ("B", "P", "F", "E", "S")
N_PERM = {"paper": 5000, "quick": 500}
MAX_MISSING = 0.1


def read_samples(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM), ``donor`` (O1..), ``arm``, ``days`` and ``donor_age`` from the series matrix."""
    rows: dict[str, list[str]] = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            key, values = cells[0], cells[1:]
            if key == "!Sample_geo_accession":
                rows["sample"] = values
            elif key == "!Sample_title":
                rows["donor"] = [v.split("_", 1)[0] for v in values]
            elif key == "!Sample_characteristics_ch1" and values and ":" in values[0]:
                name = values[0].split(":", 1)[0].strip()
                field = {
                    "cell type": "cell",
                    "length of reprogramming phase (days)": "days",
                    "donor age (years)": "donor_age",
                }.get(name)
                if field:
                    rows[field] = [v.split(":", 1)[1].strip() for v in values]
    frame = pl.DataFrame(rows)
    arm = pl.lit(None, dtype=pl.String)
    for prefix, key in ARMS.items():
        arm = pl.when(pl.col("cell").str.starts_with(prefix)).then(pl.lit(key)).otherwise(arm)
    return frame.with_columns(arm=arm, donor_age=pl.col("donor_age").cast(pl.Float64, strict=False)).drop("cell")


def within_donor_effect(y: np.ndarray, treated: np.ndarray, donor: np.ndarray) -> float:
    """Mean over donors having both arms of (treated mean − untreated mean)."""
    diffs = [
        y[(donor == d) & treated].mean() - y[(donor == d) & ~treated].mean()
        for d in np.unique(donor)
        if ((donor == d) & treated).any() and ((donor == d) & ~treated).any()
    ]
    return float(np.mean(diffs)) if diffs else float("nan")


def donor_permutation_p(y: np.ndarray, treated: np.ndarray, donor: np.ndarray, n: int, rng) -> float:
    """Two-sided p of ``within_donor_effect``, labels shuffled within each donor."""
    observed = abs(within_donor_effect(y, treated, donor))
    hits = 0
    for _ in range(n):
        perm = treated.copy()
        for d in np.unique(donor):
            idx = np.flatnonzero(donor == d)
            perm[idx] = rng.permutation(treated[idx])
        hits += abs(within_donor_effect(y, perm, donor)) >= observed - 1e-12
    return (hits + 1) / (n + 1)


def run(ctx: Context) -> list[Path]:
    path = ctx.paths.interim / BETAS
    matrix = ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not path.exists() or not matrix.exists():
        print(f"gill_reprogramming: {path} or the series matrix missing (run R65a); nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("gill_reprogramming")
    samples = read_samples(matrix)
    have = set(pl.read_parquet_schema(path))
    samples = samples.filter(pl.col("sample").is_in(list(have)) & pl.col("arm").is_not_null()).sort(
        "donor", "arm", "sample"
    )
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    names = samples["sample"].to_list()
    frame = pl.scan_parquet(path).select("cpg", *names).filter(pl.col("cpg").is_in(features)).collect()
    m = FeatureMatrix.from_frame(frame, feature_col="cpg", samples=names, max_missing=MAX_MISSING)
    x = np.where(np.isnan(m.values), np.nanmean(m.values, axis=0), m.values)
    feats = m.features
    island = _islands(feats, probes)
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    member = np.any([s > 0 for s, _ in sets.values()], axis=0)
    clock = read_clock(ClockLayout(ctx.paths).file("Horvath skin & blood"))
    cframe = pl.scan_parquet(path).select("cpg", *names).filter(pl.col("cpg").is_in(clock["cpg"].implode())).collect()
    clock_x = cframe.drop("cpg").to_numpy().T.astype(float)
    clock_sum, clock_cpgs = linear_score(clock_x, cframe["cpg"].to_list(), clock)
    arm, donor = samples["arm"].to_numpy(), samples["donor"].to_numpy()
    rows = []
    for d in np.unique(donor):
        ctrl = np.flatnonzero((donor == d) & (arm == "control"))
        if ctrl.size < 2:
            print(f"gill_reprogramming: donor {d} has {ctrl.size} control(s), left out", file=sys.stderr)
            continue
        ref = fit_reference(x[ctrl], island, ~member)
        others = np.flatnonzero((donor == d) & (arm != "control"))
        sc = score_against(x[others], ref, sets)
        for k, i in enumerate(others):
            rows.append({"i": int(i), **{key: float(v[k]) for key, v in sc.items()}})
        for i in ctrl:
            rest = ctrl[ctrl != i]
            if rest.size < 1:
                continue
            sc1 = score_against(x[i : i + 1], fit_reference(x[rest], island, ~member), sets)
            rows.append({"i": int(i), **{key: float(v[0]) for key, v in sc1.items()}})
    scores = pl.DataFrame(rows).sort("i")
    idx = scores["i"].to_numpy()
    scores = scores.drop("i").with_columns(
        sample=pl.Series([names[i] for i in idx]),
        donor=pl.Series(donor[idx]),
        arm=pl.Series(arm[idx]),
        days=pl.Series(samples["days"].to_numpy()[idx]),
        clock_sum=pl.Series(clock_sum[idx]),
    )
    rng = ctx.rng(650_000)
    effects = []
    for contrast, (a, b) in CONTRASTS.items():
        sub = scores.filter(pl.col("arm").is_in([a, b]))
        treated = (sub["arm"] == a).to_numpy()
        dn = sub["donor"].to_numpy()
        cols = [f"{mm}_{g}" for g in sets for mm in MEASURES] + ["clock_sum"]
        for col in cols:
            y = sub[col].to_numpy()
            effects.append(
                {
                    "contrast": contrast,
                    "score": col,
                    "effect": within_donor_effect(y, treated, dn),
                    "p_perm": donor_permutation_p(y, treated, dn, N_PERM[str(ctx.scale.name)], rng),
                    "n": f"{int(treated.sum())} vs {int((~treated).sum())}",
                    "donors": int(len(np.unique(dn))),
                    "cpgs": len(feats),
                    "clock_cpgs": clock_cpgs,
                }
            )
        print(f"gill_reprogramming: {contrast}: {int(treated.sum())} vs {int((~treated).sum())}", file=sys.stderr)
    effects = pl.DataFrame(effects)
    write_table(effects, out / "effects")
    write_table(scores, out / "sample_scores")
    return [_figure(ctx, scores)]


def _figure(ctx: Context, scores: pl.DataFrame) -> Path:
    arms = ["control", "failed", "reprogrammed"]
    panels = [
        ("clock_sum", "Horvath skin & blood (raw sum)"),
        ("B_housekeeping", "drift B, housekeeping"),
        ("F_housekeeping", "shield share F, housekeeping"),
        ("E_housekeeping", "entropy gap E, housekeeping"),
    ]
    fig, axes = plots.panels(len(panels), ncols=4, size=3.2)
    donors = scores["donor"].unique(maintain_order=True).to_list()
    for ax, (col, title), letter in zip(axes, panels, "abcd", strict=True):
        for j, d in enumerate(donors):
            for k, a in enumerate(arms):
                v = scores.filter((pl.col("donor") == d) & (pl.col("arm") == a))[col].to_numpy()
                xs = np.full(v.size, k + 0.2 * (j - 1))
                ax.scatter(xs, v, s=14, color=plots.PALETTE[j], label=d if k == 0 else None)
        ax.set_xticks(range(len(arms)), arms, fontsize=7)
        ax.set_title(title, fontsize=9)
        if col != "clock_sum":
            ax.axhline(0, color=plots.MUTED, lw=0.8)
        if letter == "a":
            ax.legend(fontsize=7, title="donor", title_fontsize=7)
        plots.label_panel(ax, letter)
    fig.suptitle(
        "Transient reprogramming of human fibroblasts (Gill 2022, IDATs)", x=0.02, y=1.05, ha="left", fontsize=11
    )
    fig.subplots_adjust(wspace=0.45)
    return plots.save(fig, ctx.figure_path("gill_reprogramming"))
