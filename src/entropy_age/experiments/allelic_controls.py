"""Loci held at half methylation by allele: a built-in control for noise (roadmap RM35 item 1).

At an imprinted region one parental allele is methylated and the other is not, in every cell; in a woman's cells
one X is silenced and its CpG-island promoters methylated, the other active and unmethylated. Bulk methylation sits
near one half there because of how the two alleles are marked, not because cells disagree, and nothing in aging is
known to move it in a direction. So the rise of spread with age at these loci is drift with no directed target:
the noise floor of a CpG that is held by design. Set against level-matched look-alikes, it says whether loci that a
maintenance machinery holds (ZFP57/KAP1 at imprints, XIST at the inactive X) are shielded like essential promoters.

**Imprinted arm** (Hannum blood, Jaffe cortex): CpGs linked to a gene that geneimprint lists as ``Imprinted`` (131
human genes, list read 2026-11-16), in an island or its shores, whose young level (youngest 20 %) is between 0.25
and 0.75. **X arm** (Jaffe cortex, Lunnon blood and frontal cortex; Hannum's deposited matrix has no X probes):
X-linked island CpGs, in women (where silencing puts them near one half) and, as a contrast, in men (one active X,
unmethylated). Each group is compared with every other autosomal CpG within strata of island relation and young
level (20 bins), the same matched gap as the shield's (``tissue_contrast._matched_sd``), with the essential sets
beside it for scale; the background is 60,000 random autosomal CpGs, as in ``tissue_contrast``. A donor bootstrap
gives the intervals.

``entropy-age figures allelic_controls`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, rescale, write_table
from entropy_age.experiments.tissue_contrast import (
    SETS,
    Prepared,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _lunnon,
    _settings,
)
from entropy_age.stochasticity import baseline_strata, compare_groups

IMPRINT_LIST = "geneimprint_human_2026-11-16.csv"
METRICS = ("dispersion_slope", "kl_slope", "abs_slope")
HALF = (0.25, 0.75)
ISLAND_OR_SHORE = ("Island", "N_Shore", "S_Shore")
BINS = 20
N_BOOT_MAX = 100
BACKGROUND = 60_000


def imprinted_genes(path: Path) -> set[str]:
    t = pl.read_csv(path)
    return set(t.filter(pl.col("status") == "Imprinted")["gene"].to_list())


def half_methylated(metrics: pl.DataFrame, candidates: set[str]) -> set[str]:
    """Candidates whose young level sits near one half."""
    m = metrics.filter(pl.col("feature").is_in(list(candidates)) & pl.col("baseline").is_between(*HALF))
    return set(m["feature"].to_list())


def strata(metrics: pl.DataFrame, probes: pl.DataFrame) -> pl.DataFrame:
    m = metrics.join(probes.select(pl.col("cpg").alias("feature"), "island_relation"), on="feature", how="left")
    b = baseline_strata(m["baseline"].to_numpy(), BINS)
    return m.with_columns(stratum=pl.format("{}|{}", pl.col("island_relation").fill_null("OpenSea"), pl.Series(b)))


def gaps(metrics: pl.DataFrame, membership: pl.DataFrame, permutations: int, seed: int) -> pl.DataFrame:
    """Matched gap of each group against the autosomal rest, per metric, in SD of the metric."""
    rows = []
    for metric in METRICS:
        t = compare_groups(
            metrics, membership, metric=metric, stratify_by="stratum", n_permutations=permutations, seed=seed
        )
        sd = float(metrics[metric].std())
        rows += [{"metric": metric, **r, "sd": sd, "matched_diff_sd": r["matched_diff"] / sd} for r in t.to_dicts()]
    return pl.DataFrame(rows)


def _subset(p: Prepared, rows: np.ndarray, label: str) -> Prepared:
    cov = None if p.covariates is None else p.covariates[rows]
    return dataclasses.replace(p, label=label, x=p.x[rows], ages=p.ages[rows], covariates=cov, donors=p.donors[rows])


def _one(p: Prepared, probes, base_membership, groups_fn, permutations, seed, rows=None) -> pl.DataFrame:
    q = p if rows is None else _subset(p, rows, p.label)
    m = strata(_lifespan_metrics(q), probes)
    membership = pl.concat([base_membership.select("feature", "group"), groups_fn(m)], how="vertical")
    keep = set(membership["feature"].to_list()) | set(m.filter(~pl.col("feature").is_in(XLINKED))["feature"])
    m = m.filter(pl.col("feature").is_in(list(keep)))
    return gaps(m, membership, permutations, seed)


XLINKED: list[str] = []  # filled in run(): X CpGs leave the background, so it stays autosomal


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("allelic_controls")
    ann = AnnotationLayout(ctx.paths)
    listing = ann.input / IMPRINT_LIST
    if not listing.exists():
        print(f"allelic_controls: {listing} missing", file=sys.stderr)
        return []
    _, n_boot = _settings(ctx)
    n_boot = min(n_boot, N_BOOT_MAX)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("chr") != "Y")
    links = pl.read_parquet(ann.links_450k)
    sets = link_membership(links, pl.read_parquet(ann.gene_sets)).filter(pl.col("group").is_in(SETS))
    promoters = probes.filter((pl.col("chr") != "X") & pl.col("promoter"))
    sets = sets.filter(pl.col("feature").is_in(promoters["cpg"].implode()))  # the shield's own CpGs, for scale
    imprinted = imprinted_genes(listing)
    island = probes.filter(pl.col("island_relation").is_in(ISLAND_OR_SHORE))
    imp_cands = set(
        links.filter(pl.col("symbol").is_in(list(imprinted)))
        .join(island.filter(pl.col("chr") != "X").select("cpg"), on="cpg")["cpg"]
        .to_list()
    )
    x_cands = set(probes.filter((pl.col("chr") == "X") & (pl.col("island_relation") == "Island"))["cpg"].to_list())
    XLINKED[:] = probes.filter(pl.col("chr") == "X")["cpg"].to_list()
    autosomal = probes.filter(pl.col("chr") != "X")["cpg"].to_numpy()
    background = ctx.rng(34).choice(autosomal, size=min(BACKGROUND, autosomal.size), replace=False)
    features = sorted(set(background.tolist()) | imp_cands | x_cands | set(sets["feature"].to_list()))

    def imprint_groups(m):
        half = half_methylated(m, imp_cands)
        return pl.DataFrame({"feature": sorted(half), "group": "imprinted, half-methylated"})

    def x_groups(m):
        half = half_methylated(m, x_cands)
        low = set(m.filter(pl.col("feature").is_in(list(x_cands)) & (pl.col("baseline") < HALF[0]))["feature"])
        return pl.concat(
            [
                pl.DataFrame({"feature": sorted(half), "group": "X island, half-methylated"}),
                pl.DataFrame({"feature": sorted(low), "group": "X island, unmethylated"}),
            ]
        )

    def both(m):
        return pl.concat([imprint_groups(m), x_groups(m)])

    jobs: list[tuple[str, Prepared, np.ndarray | None, object]] = []
    hannum = _hannum(ctx.paths, features)
    jobs.append(("blood (Hannum)", hannum, None, imprint_groups))
    jaffe = _jaffe(ctx.paths, features)
    sex = _sex("jaffe", ctx, jaffe)
    jobs.append(("cortex (Jaffe)", jaffe, None, imprint_groups))
    for s in ("F", "M"):
        jobs.append((f"cortex (Jaffe), {'women' if s == 'F' else 'men'}", jaffe, np.flatnonzero(sex == s), x_groups))
    lunnon = _lunnon(ctx.paths, features)
    for tissue in ("whole blood", "frontal cortex"):
        p = lunnon[tissue]
        p = dataclasses.replace(p, covariates=p.technical)  # age is the regressor here, not a covariate
        sx = _sex("lunnon", ctx, p)
        for s in ("F", "M"):
            jobs.append((f"{tissue} (Lunnon), {'women' if s == 'F' else 'men'}", p, np.flatnonzero(sx == s), both))

    point, boot = [], []
    rng = ctx.rng(35)
    for label, p, rows, fn in jobs:
        p = dataclasses.replace(p, x=p.x.astype(np.float32))
        idx = np.arange(p.x.shape[0]) if rows is None else rows
        t = _one(p, probes, sets, fn, 1000, 1, idx)
        point.append(t.with_columns(tissue=pl.lit(label), samples=pl.lit(idx.size)))
        for b in range(n_boot):
            draw = rng.choice(idx, size=idx.size, replace=True)
            tb = _one(p, probes, sets, fn, 0, b, draw)
            boot.append(tb.select("metric", "group", "matched_diff").with_columns(tissue=pl.lit(label), b=pl.lit(b)))
        print(f"allelic_controls: {label}: {idx.size} samples", file=sys.stderr, flush=True)
    point = pl.concat(point, how="diagonal")
    boot = pl.concat(boot) if boot else pl.DataFrame()
    if boot.height:
        boot = rescale(boot, point, ["tissue", "metric", "group"])  # the full sample's sd (audit 2026-11-18)
        ci = boot.group_by("tissue", "metric", "group").agg(
            lo=pl.col("matched_diff_sd").quantile(0.025), hi=pl.col("matched_diff_sd").quantile(0.975)
        )
        point = point.join(ci, on=["tissue", "metric", "group"], how="left")
        boot.write_parquet(out / "bootstrap.parquet")
    write_table(point, out / "gaps")
    return [_figure(ctx, point)]


def _sex(cohort: str, ctx: Context, p: Prepared) -> np.ndarray:
    """``F``/``M`` per sample of ``p``, read from the cohort's samples table by donor."""
    s = pl.read_parquet(cohort_layout(cohort, ctx.paths).samples)
    col = "sex"
    first = s.group_by("donor").agg(pl.col(col).first())
    lookup = dict(zip(first["donor"].to_list(), first[col].to_list(), strict=True))
    raw = [str(lookup.get(d, "")).strip().upper()[:1] for d in p.donors]
    return np.array(raw)


def _figure(ctx: Context, point: pl.DataFrame) -> Path:
    t = point.filter(pl.col("metric") == "dispersion_slope")
    groups = ["imprinted, half-methylated", "X island, half-methylated", "X island, unmethylated", *SETS]
    tissues = t["tissue"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(1, ncols=1, size=max(4.5, 0.42 * len(tissues) * 2))
    ax = axes[0]
    y, labels, k = [], [], 0
    for tissue in tissues:
        for i, g in enumerate(groups):
            r = t.filter((pl.col("tissue") == tissue) & (pl.col("group") == g))
            if r.height == 0:
                continue
            v = r["matched_diff_sd"][0]
            lo = r["lo"][0] if "lo" in r.columns else None
            hi = r["hi"][0] if "hi" in r.columns else None
            color = plots.PALETTE[i % len(plots.PALETTE)]
            ax.plot([v], [k], "o", color=color)
            if lo is not None and hi is not None:
                ax.plot([lo, hi], [k, k], color=color, lw=1.2)
            y.append(k)
            labels.append(f"{tissue}: {g}")
            k += 1
        k += 0.6
    ax.set_yticks(y, labels, fontsize=6)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="rise of spread with age, matched gap (SD)", title="Loci held at one half by allele")
    return plots.save(fig, ctx.figure_path("allelic_controls"))
