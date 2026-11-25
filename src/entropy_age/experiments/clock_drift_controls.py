"""Can a broad array shift make the clock-shield link? Confounder and placebo controls for ``clock_drift``.

``clock_drift`` finds that a mitotic clock's acceleration tracks the per-person shield share F (r about -0.5 outside
the brain; HypoClock, whose value falls with division, about +0.3). F is a ratio: with A the housekeeping promoters'
deviation from the leave-one-out reference and B that of their look-alikes, ``F = (A - B) / B``. A shift of the
whole array (a methylation level change, beta compression toward 0.5, chip, signal or detection quality, cell
mixture) that moves A and B by different amounts moves F, and the same shift can move a clock (compression raises
low-beta promoter clock CpGs and lowers high-beta HypoClock CpGs at once). Codex's second review (item 3) proposed
the check built here.

Per series group, with the arrays, CpGs and clock scores of ``clock_drift`` (``series.selection``, the same
availability rule, every mitotic and EpiTOC1 CpG kept out of the drift CpGs, sensitivity configs flagged):

1. **Drift measures.** ``loo_scores`` per person: ``P_housekeeping`` (the matched gap), ``B_housekeeping`` and
   ``F_housekeeping``; ``A = P + B`` is the housekeeping CpGs' own stratum-weighted deviation.
2. **A third, disjoint probe panel.** A seeded random sample of 20,000 CpGs of the same rebuild, measured in every
   array of the group, outside the drift CpGs, outside every mitotic clock's CpGs and outside the CpGs of every
   biolearn clock with a coefficient file (``clock_cpgs.clock_table``). Per person: its mean beta (``level``),
   ``compression`` = mean ``|beta - 0.5|`` (lower is more compressed toward 0.5) and the SD of its betas (reported).
3. **Covariates.** Age; chip dummies, log mean intensity and detection fraction of the rebuild the drift CpGs come
   from (``reprocessed.read_qc`` through ``clock_disagreement.quality_terms``, chip dropped when it would leave under
   10 residual degrees of freedom); ``batch`` or ``center`` dummies where the series description holds them. The
   series groups carry no cell fractions (only EPIC-Italy has EpiDISH fractions, and it is not a series group), so
   none enter; level and compression stand in for a mixture that moves the whole array.
4. **Same adjustment on both sides.** The clock's acceleration (residual on age, as in ``clock_drift``), A, B and F
   are all residualised on age + covariates + level + compression, then correlated. Reported per group and clock:
   ``r_F_raw`` (``clock_drift``'s value), ``r_F_adj``, ``r_A_given_B`` (acceleration with A at a given B, the
   adjustment set plus B), ``r_level``, ``r_compression``, ``r_sd`` (the acceleration against the panel indices),
   ``r_F_level`` and ``r_F_compression`` (F against them), and n.
5. **Placebo sets.** K = 20 sets per group, each as many CpGs as the housekeeping set in every stratum of island
   status x decile of the group's mean beta (the strata ``fit_reference`` matches on), drawn gene by gene from
   promoter CpGs of genes in none of the three gene sets (the pool ``placebo_sets`` draws from), restricted to the
   drift CpGs of the group so the reference and look-alikes are unchanged. They are scored in the same
   leave-one-out pass as the real set (``loo_drift``: ``loo_scores``' reference, deviations and stratum weights,
   computing only P, B and F; no components, so the reference does not depend on the sets). Per group and
   clock: the median and 95th percentile of ``|r(acc, F_placebo)|``, raw and adjusted, and the share of placebos at
   or beyond the real r in its own direction. ``placebo_fill`` is the share of the housekeeping quota the draw filled.
   A placebo's look-alikes are the CpGs outside it and outside the housekeeping set, the pool the housekeeping set
   is compared with; a first run that left the housekeeping CpGs among a placebo's look-alikes gave placebo r of the
   real r's size and opposite sign (the placebo's B then carries the housekeeping CpGs), which is no null.

Writes ``clock_drift_controls/correlations`` (per group and clock), ``placebo`` (per group, clock and draw),
``summary`` (per mitotic clock, medians over independent non-brain groups: sensitivity configs out, brain groups out,
the largest group per series) and a figure. A Slurm job: the leave-one-out scoring runs once per group, with the
placebo sets inside it.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments import mitotic_clocks
from entropy_age.experiments.clock_cpgs import clock_table
from entropy_age.experiments.clock_disagreement import partial_r, quality_terms
from entropy_age.experiments.clock_drift import CLOCK_SCORES, CLOCK_TRACK, GROUP, MIN_ARRAYS, SERIES, residual
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import SENSITIVITY, selection, specs
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import contrast_weights, deviations, fit_reference
from entropy_age.stochasticity import baseline_strata

PANEL = 20_000
N_PLACEBO = {"paper": 20, "quick": 5}
MITOTIC = ("EpiTOC1", "epiTOC2", "epiTOC3", "HypoClock", "RepliTali", "epiCMIT")
HELD = ("batch", "center")  # technical columns a series description may hold
BRAIN = re.compile(r"brain|cortex|cerebell|hippocamp|amygdala|neuron|gyrus|striatum|caudate|glia|entorhinal", re.I)


# ----------------------------------------------------------------------------- pure helpers
def panel_indices(x: np.ndarray) -> dict[str, np.ndarray]:
    """Per person (rows of ``x``, people x CpGs, no NaN): mean beta, mean ``|beta - 0.5|`` and SD of the betas."""
    x = np.asarray(x, dtype=float)
    return {"level": x.mean(axis=1), "compression": np.abs(x - 0.5).mean(axis=1), "sd": x.std(axis=1)}


def corr(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def matched_placebo(
    strata: np.ndarray, real: np.ndarray, pool_genes: list[np.ndarray], rng: np.random.Generator
) -> tuple[np.ndarray, list[np.ndarray], float]:
    """A placebo set with ``real``'s count of CpGs in every stratum, filled gene by gene from ``pool_genes`` (index
    arrays of candidate CpGs, one per gene) in random order; a gene adds those of its CpGs whose stratum still has
    room. Returns (membership mask, the genes' index arrays as used, share of the quota filled)."""
    levels, codes = np.unique(strata, return_inverse=True)
    quota = np.bincount(codes[real], minlength=levels.size)
    want = int(quota.sum())
    taken = np.zeros(strata.size, dtype=bool)
    genes = []
    for g in rng.permutation(len(pool_genes)):
        if quota.sum() == 0:
            break
        use = []
        for j in pool_genes[g]:
            if not taken[j] and not real[j] and quota[codes[j]] > 0:
                quota[codes[j]] -= 1
                taken[j] = True
                use.append(j)
        if use:
            genes.append(np.array(use, dtype=int))
    return taken, genes, float(taken.sum() / want) if want else float("nan")


def placebo_pool(kept: list[str], member: np.ndarray, links: pl.DataFrame, set_genes: set[str]) -> list[np.ndarray]:
    """Per gene in none of the gene sets: the indices (in ``kept``) of its promoter CpGs outside every set."""
    index = {c: j for j, c in enumerate(kept)}
    rows = links.filter(pl.col("promoter") & pl.col("cpg").is_in(kept) & ~pl.col("symbol").is_in(list(set_genes)))
    inside = set(np.flatnonzero(member).tolist())
    out = []
    for _, g in rows.group_by("symbol", maintain_order=True):
        idx = sorted({index[c] for c in g["cpg"].to_list()} - inside)
        if idx:
            out.append(np.array(idx, dtype=int))
    out.sort(key=lambda a: int(a[0]))  # a fixed order before the seeded shuffle
    return out


def loo_drift(
    x: np.ndarray, island: np.ndarray, sets: dict, drop: dict[str, np.ndarray] | None = None
) -> dict[str, np.ndarray]:
    """``loo_scores``' ``P``, ``B`` and ``F`` of every set (no components), computed once per person for all sets:
    the same leave-one-out reference (``fit_reference``), deviations and stratum weights (``contrast_weights``), with
    the entropy, relative-entropy and signed gaps left out (``F`` is all the placebo sets need). ``drop`` maps a set
    to CpGs left out of its contrast altogether (neither in it nor among its look-alikes)."""
    drop = drop or {}
    member = np.any([s for s, _ in sets.values()], axis=0)
    names = list(sets)
    out = {f"{m}_{g}": np.empty(x.shape[0]) for g in names for m in "PBF"}
    for i in range(x.shape[0]):
        others = np.arange(x.shape[0]) != i
        ref = fit_reference(x[others], island, ~member)
        codes = np.unique(ref.strata, return_inverse=True)[1]
        d = deviations(x[i], ref.level)
        for g in names:
            keep = ~drop[g] if g in drop else np.ones(d.size, dtype=bool)
            p, b = contrast_weights(np.asarray(sets[g][0][keep], dtype=float), codes[keep])
            out[f"P_{g}"][i], out[f"B_{g}"][i] = d[keep] @ p, d[keep] @ b
    for g in names:
        out[f"F_{g}"] = out[f"P_{g}"] / out[f"B_{g}"]
    return out


def adjust_design(age: np.ndarray, quality: np.ndarray, held: np.ndarray, level: np.ndarray, comp: np.ndarray):
    return np.column_stack([age, quality, held, level, comp])


# ----------------------------------------------------------------------------- run
def _panel(path: Path, samples: list[str], exclude: set[str], rng: np.random.Generator) -> np.ndarray:
    """People x PANEL betas of a seeded sample of CpGs measured in every array, outside ``exclude``."""
    ids = pl.scan_parquet(path).select("cpg").collect()["cpg"]
    pool = np.array(sorted(set(ids.to_list()) - exclude))
    order = rng.permutation(pool.size)[: min(pool.size, 3 * PANEL)]
    pick = pool[np.sort(order)]
    frame = pl.scan_parquet(path).filter(pl.col("cpg").is_in(pick.tolist())).select("cpg", *samples).collect()
    frame = frame.drop_nulls().filter(pl.all_horizontal(pl.col(samples).is_not_nan()))
    rank = {c: i for i, c in enumerate(pool[order].tolist())}  # keep the first PANEL complete CpGs in draw order
    frame = frame.with_columns(r=pl.col("cpg").replace_strict(rank, return_dtype=pl.Int64)).sort("r").head(PANEL)
    return frame.select(samples).to_numpy().T.astype(float)


def _held(s: pl.DataFrame) -> np.ndarray:
    parts = [dummies(s[c].cast(pl.String).to_numpy()) for c in HELD if c in s.columns and s[c].null_count() == 0]
    parts = [p for p in parts if p.shape[1]]
    return np.hstack(parts) if parts else np.zeros((s.height, 0))


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # clock_drift's stream: the same drift CpGs
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("cpg").is_in(features))
    links = pl.read_parquet(ann.links_450k)
    set_genes = set(membership["gene"].to_list())
    mitotic = mitotic_clocks.load(ctx.paths.root)
    if mitotic is None:
        raise FileNotFoundError(f"clock_drift_controls: no mitotic coefficients in {mitotic_clocks.DIR}")
    clock_cpgs = mitotic_clocks.cpgs(mitotic)
    other_clocks = set(clock_table(ctx.paths)["cpg"].to_list())
    n_placebo = N_PLACEBO[str(ctx.scale.name)]
    rows, draws = [], []
    for gi, (name, spec) in enumerate(specs().items()):
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv"
        unmasked = ctx.paths.root / CLOCK_SCORES / f"{gse}_cpgs.parquet"
        if not (betas_path.exists() and clocks_path.exists()):
            continue
        long = pl.read_csv(clocks_path, schema_overrides={"sample": pl.String}).select("sample", "clock", "prediction")
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(clock_cpgs))
        mb = None
        available = set(betas.columns) & set(long["sample"].to_list())
        if unmasked.exists():
            mb = pl.scan_parquet(unmasked).filter(pl.col("cpg").is_in(clock_cpgs)).collect()
            available &= set(mb.columns)
            available &= set(
                reprocessed.read_qc(ctx.paths, gse, sesame_dir=CLOCK_TRACK).filter(pl.col("passed"))["sample"]
            )
        chosen = selection(ctx.paths, name, spec, available)
        if chosen is None:
            continue
        groups, dropped = chosen
        if dropped:
            betas = betas.filter(~pl.col("cpg").is_in(dropped))
        base = long.cast({"prediction": pl.Float64}).pivot(on="clock", index="sample", values="prediction")
        qc = reprocessed.read_qc(ctx.paths, gse)  # the rebuild the drift CpGs and the panel come from
        exclude = set(features) | clock_cpgs | other_clocks | set(dropped)
        for li, (label, s) in enumerate(groups.items()):
            s = s.drop_nulls("age")
            if s.height < MIN_ARRAYS:
                continue
            samples = s["sample"].to_list()
            scores = base
            if mb is not None:
                extra, _ = mitotic_clocks.score(mb, samples, mitotic)
                if extra.height:
                    scores = base.join(
                        extra.pivot(on="clock", index="sample", values="prediction"), on="sample", how="left"
                    )
            sub = betas.select("cpg", *samples)
            x = sub.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            x = x[:, ok]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                continue
            island = _islands(kept, probes)
            member = np.any([m for m, _ in sets.values()], axis=0)
            strata = np.char.add(island.astype(str), baseline_strata(x.mean(axis=0), 10).astype(str))
            rng = ctx.rng(130_000 + 100 * gi + li)
            pool = placebo_pool(kept, member, links, set_genes)
            fills = []
            scored_sets = dict(sets)
            for k in range(n_placebo):
                mask, genes, fill = matched_placebo(strata, sets[GROUP][0], pool, rng)
                scored_sets[f"placebo{k}"] = (mask, genes)
                fills.append(fill)
            panel = _panel(betas_path, samples, exclude, rng)
            print(
                f"clock_drift_controls: {name} / {label}: {len(samples)} arrays, {len(kept)} CpGs, panel "
                f"{panel.shape[1]}, placebo fill {np.mean(fills):.3f}",
                file=sys.stderr,
                flush=True,
            )
            # a placebo's look-alikes leave the housekeeping CpGs out, so placebo and housekeeping sets are compared
            # with the same pool (with them in, a placebo's B carries the housekeeping CpGs and mirrors the real F)
            hk = sets[GROUP][0]
            per = loo_drift(x, island, scored_sets, {k: hk for k in scored_sets if k.startswith("placebo")})
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            a = per[f"P_{GROUP}"] + b
            ind = panel_indices(panel)
            age = s["age"].to_numpy().astype(float)
            held = _held(s)
            quality, qdesc, _ = quality_terms(samples, qc, base_columns=6 + held.shape[1])
            clk = s.select("sample").join(scores, on="sample", how="left", validate="1:1")
            if clk["sample"].to_list() != samples:
                raise ValueError(f"clock_drift_controls: {name} / {label}: clock scores out of the arrays' order")
            for clock in [c for c in clk.columns if c != "sample"]:
                v = clk[clock].to_numpy().astype(float)
                use = np.isfinite(v)
                if use.sum() < MIN_ARRAYS or np.std(v[use]) == 0:
                    continue
                acc = residual(v[use], age[use])
                design = adjust_design(age, quality, held, ind["level"], ind["compression"])[use]
                r_raw = corr(acc, f[use])
                r_adj = partial_r(acc, f[use], design)
                rows.append(
                    {
                        "series": name,
                        "group": label,
                        "sensitivity": bool(SENSITIVITY.search(name)),
                        "brain": bool(BRAIN.search(f"{name} {label}")),
                        "n": int(use.sum()),
                        "cpgs": len(kept),
                        "panel_cpgs": int(panel.shape[1]),
                        "covariates": qdesc + (", batch/center" if held.shape[1] else ""),
                        "clock": clock,
                        "mitotic": clock in MITOTIC,
                        "r_F_raw": r_raw,
                        "r_B_raw": corr(acc, b[use]),
                        "r_F_adj": r_adj,
                        "r_B_adj": partial_r(acc, b[use], design),
                        "r_A_given_B": partial_r(acc, a[use], np.column_stack([design, b[use]])),
                        "r_level": corr(acc, ind["level"][use]),
                        "r_compression": corr(acc, ind["compression"][use]),
                        "r_sd": corr(acc, ind["sd"][use]),
                        "r_F_level": corr(f[use], ind["level"][use]),
                        "r_F_compression": corr(f[use], ind["compression"][use]),
                        "placebo_fill": float(np.mean(fills)),
                    }
                )
                for k in range(n_placebo):
                    fp = per[f"F_placebo{k}"][use]
                    draws.append(
                        {
                            "series": name,
                            "group": label,
                            "clock": clock,
                            "draw": k,
                            "r_raw": corr(acc, fp),
                            "r_adj": partial_r(acc, fp, design),
                        }
                    )
    if not rows:
        print("clock_drift_controls: nothing to score", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_drift_controls")
    table = pl.DataFrame(rows)
    placebo = pl.DataFrame(draws)
    table = table.join(placebo_band(table, placebo), on=["series", "group", "clock"], how="left")
    write_table(table, out / "correlations", csv_max_rows=10_000)
    write_table(placebo, out / "placebo", csv_max_rows=10_000)
    summary = summarise(table)
    write_table(summary, out / "summary", csv_max_rows=10_000)
    return [draw(independent(table), ctx.figure_path("clock_drift_controls"))]


def placebo_band(table: pl.DataFrame, placebo: pl.DataFrame) -> pl.DataFrame:
    """Per group and clock: median and 95th percentile of |r| over the placebo draws, raw and adjusted, and the share
    of draws at or beyond the real r in its own direction."""
    joined = placebo.join(
        table.select("series", "group", "clock", "r_F_raw", "r_F_adj"), on=["series", "group", "clock"]
    )
    beyond = lambda real, p: ((pl.col(p) * pl.col(real).sign()) >= pl.col(real).abs()).mean()  # noqa: E731
    return joined.group_by("series", "group", "clock", maintain_order=True).agg(
        placebo_abs_median_raw=pl.col("r_raw").abs().median(),
        placebo_abs_p95_raw=pl.col("r_raw").abs().quantile(0.95, "linear"),
        placebo_abs_median_adj=pl.col("r_adj").abs().median(),
        placebo_abs_p95_adj=pl.col("r_adj").abs().quantile(0.95, "linear"),
        placebo_mean_raw=pl.col("r_raw").mean(),
        placebo_beyond_raw=beyond("r_F_raw", "r_raw"),
        placebo_beyond_adj=beyond("r_F_adj", "r_adj"),
        placebos=pl.len(),
    )


def independent(table: pl.DataFrame) -> pl.DataFrame:
    """Non-brain, non-sensitivity groups, the largest group per series (groups of one series can share donors)."""
    t = table.filter(~pl.col("sensitivity") & ~pl.col("brain"))
    largest = (
        t.group_by("series", "group")
        .agg(pl.col("n").max())
        .sort(["series", "n", "group"], descending=[False, True, False])
        .unique("series", keep="first")
        .select("series", "group")
    )
    return t.join(largest, on=["series", "group"])


def summarise(table: pl.DataFrame) -> pl.DataFrame:
    cols = (
        "r_F_raw",
        "r_F_adj",
        "r_A_given_B",
        "r_level",
        "r_compression",
        "r_F_level",
        "r_F_compression",
        "placebo_abs_median_raw",
        "placebo_abs_p95_raw",
        "placebo_abs_median_adj",
        "placebo_abs_p95_adj",
        "placebo_beyond_raw",
        "placebo_beyond_adj",
    )
    return (
        independent(table)
        .filter(pl.col("mitotic"))
        .group_by("clock", maintain_order=True)
        .agg(pl.len().alias("groups"), *(pl.col(c).median().alias(c) for c in cols))
        .sort("clock")
    )


def draw(table: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    d = table.filter(pl.col("mitotic"))
    clocks = [c for c in MITOTIC if c in d["clock"].to_list()]
    fig, axes = plots.plt.subplots(1, 2, figsize=(12, 4.2), sharey=False)
    ax = axes[0]
    for i, c in enumerate(clocks):
        t = d.filter(pl.col("clock") == c)
        for off, col, color in ((-0.15, "r_F_raw", plots.MUTED), (0.15, "r_F_adj", plots.PALETTE[0])):
            v = t[col].drop_nans().to_numpy()
            jitter = np.random.default_rng(i).normal(0, 0.04, v.size)
            ax.scatter(np.full(v.size, i + off) + jitter, v, s=10, alpha=0.6, color=color)
            if v.size:
                ax.plot([i + off - 0.12, i + off + 0.12], [np.median(v)] * 2, color=plots.INK, lw=2)
        for raw, adj in t.select("r_F_raw", "r_F_adj").iter_rows():
            ax.plot([i - 0.15, i + 0.15], [raw, adj], color=plots.MUTED, lw=0.4, alpha=0.5)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_xticks(range(len(clocks)), clocks, fontsize=7)
    ax.set_ylabel("r, clock acceleration with F\n(grey: age only; blue: age + quality + level + compression)")
    ax.set_title("a  raw and adjusted, independent non-brain groups", loc="left", fontsize=9)
    ax = axes[1]
    for i, c in enumerate(clocks):
        t = d.filter(pl.col("clock") == c).sort("r_F_raw")
        xs = i + np.linspace(-0.3, 0.3, max(t.height, 1))[: t.height]
        ax.vlines(xs, -t["placebo_abs_p95_raw"], t["placebo_abs_p95_raw"], color=plots.MUTED, lw=2, alpha=0.5)
        ax.scatter(xs, t["r_F_raw"], s=12, color=plots.PALETTE[1], zorder=3)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_xticks(range(len(clocks)), clocks, fontsize=7)
    ax.set_ylabel("r with F (dot: housekeeping; bar: ±95th pct of |r|\nover matched placebo sets)")
    ax.set_title("b  housekeeping against matched placebo promoters", loc="left", fontsize=9)
    fig.tight_layout()
    return plots.save(fig, path)
