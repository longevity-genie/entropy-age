"""Fibroblasts in culture: does drift with divisions follow the gene's expression in the same cultures? (RM40 item 4)

Sturm et al. 2022 grew dermal fibroblasts of healthy donors for months, with methylation on the EPIC array
(GSE179847, rebuilt from IDATs) and RNA-seq (GSE179848) taken from the same cultures. In monocytes of living people
the shield follows expression measured in the same cells (R122). Here the same question is asked where time is
population doublings and nothing but division and culture acts on the cells.

Untreated cultures of the healthy lines at ambient oxygen are used. Per promoter CpG, the slopes of entropy and
dispersion on doublings with a line intercept (``drift_metrics``, doublings in the place of age). Expression per
gene is the mean log2(TPM + 1) over the untreated healthy-line RNA samples at 21 % oxygen. Then, as in
``mesa_expression``: set CpGs against matched CpGs with and without the gene's expression decile in the matching,
and, among genes in no set, the slope relative to matched CpGs by expression decile.

Limits: one study; chip is not held (lines and time points are spread over chips unevenly, R94); expression is a
culture-wide mean, not paired sample by sample. ``entropy-age figures sturm_expression`` (needs the SeSAMe rebuild).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, chip_terms, dummies, write_table
from entropy_age.experiments.mesa_expression import BINS, METRICS, NO_PROBE, cpg_expression, expression_bins, gradient
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import Prepared, _features, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

METHYLATION, EXPRESSION = "GSE179847", "GSE179848"
TPM = "GSE179848_txi_gene_abundance.csv.gz"
MIN_PER_LINE = 5


def control_cultures(samples: pl.DataFrame) -> pl.DataFrame:
    """Untreated cultures of healthy lines at ambient oxygen with doublings recorded; lines with few samples out."""
    keep = samples.filter(
        (pl.col("agent") == "control")
        & (pl.col("oxygen") == "Ambient")
        & pl.col("donor").str.starts_with("HC")
        & pl.col("doublings").is_not_null()
    )
    sizes = keep.group_by("donor").len().filter(pl.col("len") >= MIN_PER_LINE)["donor"].to_list()
    return keep.filter(pl.col("donor").is_in(sizes)).sort("donor", "doublings")


def control_rna(chars: pl.DataFrame) -> list[str]:
    """RNA-seq sample numbers (``rnaseq_sampleid``) of untreated healthy-line cultures at 21 % oxygen."""
    keep = chars.filter(
        (pl.col("treatments") == "Control")
        & (pl.col("clinical_condition") == "Normal")
        & (pl.col("percent_oxygen").cast(pl.Float64, strict=False) == 21)
    )
    return keep["rnaseq_sampleid"].to_list()


def gene_means(path: Path, sample_ids: list[str]) -> dict[str, float]:
    """Gene symbol -> mean log2(TPM + 1) over the columns ``<id>_lenscaleTPM`` of ``sample_ids``."""
    with gzip.open(path, "rt") as handle:
        header = [c.strip('"') for c in handle.readline().rstrip("\n").split(",")]
    wanted = [c for c in header[1:] if c.split("_")[0] in set(sample_ids)]
    frame = pl.read_csv(path, columns=[0, *[header.index(c) for c in wanted]], null_values=["NA", ""])
    values = np.log2(frame.select(frame.columns[1:]).to_numpy().astype(float) + 1.0)
    return dict(zip(frame[frame.columns[0]].to_list(), np.nanmean(values, axis=1).tolist(), strict=True))


def chip_overlap(doublings: np.ndarray, lines: np.ndarray, chips: np.ndarray) -> pl.DataFrame:
    """Per chip: cultures, lines and the range of doublings on it; and, on every row, the share of the variance of
    doublings left within line and chip (0 means chip and line fix the doublings: nothing can be estimated)."""
    design = np.column_stack([np.ones(doublings.size), dummies(lines), dummies(chips)])
    coef, *_ = np.linalg.lstsq(design, doublings, rcond=None)
    within = float(((doublings - design @ coef) ** 2).sum() / ((doublings - doublings.mean()) ** 2).sum())
    rows = [
        {
            "chip": c,
            "cultures": int((chips == c).sum()),
            "lines": int(np.unique(lines[chips == c]).size),
            "doublings_min": float(doublings[chips == c].min()),
            "doublings_max": float(doublings[chips == c].max()),
            "variance_left_within_line_and_chip": within,
        }
        for c in sorted(set(chips.tolist()))
    ]
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    rna_dir = paths.input / "geo" / EXPRESSION
    if not (reprocessed.available(paths, METHYLATION) and (rna_dir / TPM).exists()):
        print("sturm_expression: needs the SeSAMe rebuild of GSE179847 and GSE179848's TPM table.", file=sys.stderr)
        return []
    out = ctx.out_dir("sturm_expression")
    scale = str(ctx.scale.name)
    features, probes, membership = _features(paths, ctx.rng(140_000))
    qc = reprocessed.read_qc(paths, METHYLATION).filter(pl.col("passed"))
    passed = set(qc["sample"].to_list())
    cultures = control_cultures(pl.read_parquet(cohort_layout("sturm", paths).samples))
    cultures = cultures.filter(pl.col("sample").is_in(list(passed)))
    betas = reprocessed.by_accession(paths, METHYLATION, {s: s for s in cultures["sample"].to_list()})
    x, feats = complete_matrix(betas, features, cultures["sample"].to_list())
    rna = control_rna(_characteristics(rna_dir / f"{EXPRESSION}_series_matrix.txt.gz"))
    gene_mean = gene_means(rna_dir / TPM, rna)
    print(
        f"sturm_expression: {cultures.height} cultures of {cultures['donor'].n_unique()} lines, {len(feats)} CpGs, "
        f"{len(rna)} RNA samples, {len(gene_mean)} genes",
        file=sys.stderr,
        flush=True,
    )
    lines = dummies(cultures["donor"].to_numpy())
    doublings = cultures["doublings"].to_numpy().astype(float)
    chip_of = dict(qc.select("sample", pl.col("chip").cast(pl.String)).iter_rows())
    chips = np.array([chip_of[s] for s in cultures["sample"].to_list()])
    write_table(chip_overlap(doublings, cultures["donor"].to_numpy(), chips), out / "chip_overlap")
    if chip_terms():  # doublings estimated within chips (RM43 item 5)
        lines = np.column_stack([lines, dummies(chips)])
    donors = cultures["donor"].to_numpy()
    metrics = _lifespan_metrics(Prepared("fibroblasts (Sturm)", x, doublings, lines, feats, donors, lines))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    level, gene = cpg_expression(feats, pl.read_parquet(AnnotationLayout(paths).links_450k), gene_mean)
    bins = expression_bins(level)
    levels = baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
    base = np.char.add(np.char.add(island.astype(str), "|"), levels)
    variants = {"island x level": base, "+ own expression decile": np.char.add(np.char.add(base, "|"), bins)}
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    in_any = np.any([s for s, _ in sets.values()], axis=0)
    gaps, grad = [], []
    for metric in METRICS:
        values = metrics[metric].to_numpy()
        gaps += [{"metric": metric, **r} for r in matched_gaps(values, sets, variants, N_RESAMPLES[scale], ctx.seed)]
        grad += [{"metric": metric, **r} for r in gradient(values, base, bins, gene, ~in_any)]
    write_table(pl.DataFrame(gaps), out / "matched_gaps")
    write_table(pl.DataFrame(grad), out / "gradient")
    where = [
        {"group": g, "bin": b, "cpgs": int(((bins == b) & s).sum()), "share": float(((bins == b) & s).sum() / s.sum())}
        for g, (s, _) in {**sets, "outside the sets": (~in_any, [])}.items()
        for b in [NO_PROBE, *map(str, range(1, BINS + 1))]
    ]
    write_table(pl.DataFrame(where), out / "set_expression")
    sizes = cultures.group_by("donor", maintain_order=True).agg(
        pl.len().alias("cultures"), pl.col("doublings").min().alias("min"), pl.col("doublings").max().alias("max")
    )
    write_table(sizes.with_columns(cpgs=pl.lit(len(feats)), rna_samples=pl.lit(len(rna))), out / "sizes")
    return []
