"""Epidermis of 121 women: does the shield follow the gene's expression in the same tissue of the same people?

Found by the external dataset search (RM40): epidermis of healthy women aged 21 to 79, methylation on the EPIC array
(E-MTAB-8993, 65 women) or the 450k array (E-MTAB-8992, 56 other women), and RNA-seq of all of them (E-MTAB-8995;
gene TPM table of the depositors, batch-corrected by them). The arrays are rebuilt from IDATs with the project's
recipe. The two array arms are never pooled: each is one analysis, reported side by side.

The epidermis renews within weeks, so the shield is expected to be large here. Per arm, as in ``mesa_expression``:
per promoter CpG the age slopes of entropy and dispersion (processing batch held; chip too under
``ENTROPY_AGE_CHIP=1``); set CpGs against matched CpGs, with and without the gene's expression decile in the same
women and the chromatin state of keratinocytes (Roadmap E057) in the matching; among genes in no set, the slope
relative to matched CpGs by expression decile; and, as a check of the join, the across-women correlation of promoter
methylation and expression per gene.

The depositors renamed the IDAT files, so the chip is read from each file's header (``scripts/idat_barcodes.py``,
written to ``chips.csv`` beside the SDRF) and its overlap with age is tabulated.

Limits: one laboratory; women only; cell composition of the sample is not held; expression is the depositors'
table, not requantified. ``entropy-age figures epidermis_expression`` (needs the SeSAMe rebuilds).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.chromatin import ChromatinLayout, cpg_states, read_segmentation
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, chip_terms, dummies, write_table
from entropy_age.experiments.mesa_expression import (
    BINS,
    METRICS,
    NO_PROBE,
    PERMUTATIONS,
    coupling,
    cpg_expression,
    expression_bins,
    gradient,
)
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import Prepared, _features, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

ARMS = {"EPIC": "E-MTAB-8993", "450k": "E-MTAB-8992"}
RNA = "E-MTAB-8995"
TPM = "SIT_gene_TPMs.tsv"
KERATINOCYTES = "E057"
NO_STATE = "none"


def read_women(sdrf: Path) -> pl.DataFrame:
    """``sample`` (the woman's id, which names her IDAT pair and her RNA column), ``age`` and ``batch``."""
    frame = pl.read_csv(sdrf, separator="\t", infer_schema_length=0)
    return (
        frame.select(
            pl.col("Characteristics[individual]").alias("sample"),
            pl.col("Characteristics[age]").cast(pl.Float64, strict=False).alias("age"),
            pl.col("Characteristics[batch]").alias("batch"),
        )
        .unique(subset="sample", keep="first", maintain_order=True)
        .drop_nulls("age")
    )


def read_tpm(path: Path, women: list[str]) -> tuple[list[str], np.ndarray]:
    """Ensembl gene ids (version dropped) and log2(TPM + 1), genes x ``women``; the table's id column is its last.
    The depositors batch-corrected TPM on the linear scale (ComBat), which leaves 6 % of values below zero; they are
    set to zero before the log."""
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    values = frame.select(pl.col(women).cast(pl.Float64, strict=False)).to_numpy()
    genes = [g.split(".")[0] for g in frame["gene_id"].to_list()]
    return genes, np.log2(np.clip(values, 0.0, None) + 1.0)


def to_symbols(genes: list[str], expr: np.ndarray, hgnc: pl.DataFrame) -> tuple[list[str], np.ndarray]:
    """Rows of ``expr`` renamed from Ensembl ids to approved symbols; genes without one, or sharing one, dropped."""
    symbol_of = dict(hgnc.filter(pl.col("ensembl").is_not_null()).select("ensembl", "symbol").iter_rows())
    named = [(i, symbol_of[g]) for i, g in enumerate(genes) if g in symbol_of]
    counts: dict[str, int] = {}
    for _, s in named:
        counts[s] = counts.get(s, 0) + 1
    keep = [(i, s) for i, s in named if counts[s] == 1]
    return [s for _, s in keep], expr[[i for i, _ in keep]]


def age_by_chip(ages: np.ndarray, chips: np.ndarray) -> pl.DataFrame:
    """Per chip: women, their age range, and the share of the variance of age that chip explains (one value)."""
    grand = ((ages - ages.mean()) ** 2).sum()
    within = sum(((ages[chips == c] - ages[chips == c].mean()) ** 2).sum() for c in np.unique(chips))
    r2 = float(1 - within / grand) if grand > 0 else float("nan")
    rows = [
        {
            "chip": str(c),
            "women": int((chips == c).sum()),
            "age_min": float(ages[chips == c].min()),
            "age_max": float(ages[chips == c].max()),
            "r2_age_by_chip": r2,
        }
        for c in np.unique(chips)
    ]
    return pl.DataFrame(rows)


def keratinocyte_state(feats: list[str], probes: pl.DataFrame, paths) -> np.ndarray | None:
    """Chromatin state of each CpG in keratinocytes (Roadmap E057), or ``None`` when the map is not on disk."""
    path = ChromatinLayout(paths).segmentation(KERATINOCYTES)
    if not path.exists():
        return None
    frame = pl.DataFrame({"cpg": feats}).join(probes.select("cpg", "chr", "pos"), on="cpg", how="left")
    return cpg_states(frame, read_segmentation(path)).fill_null(NO_STATE).to_numpy().astype(str)


def one_arm(ctx: Context, arm: str, accession: str, features, probes, membership, links, hgnc) -> dict[str, list]:
    paths, scale = ctx.paths, str(ctx.scale.name)
    folder = paths.input / "arrayexpress" / accession
    qc = reprocessed.read_qc(paths, accession).filter(pl.col("passed"))
    women = read_women(folder / f"{accession}.sdrf.txt").filter(pl.col("sample").is_in(qc["sample"].implode()))
    header = (paths.input / "arrayexpress" / RNA / TPM).open().readline().rstrip("\n").split("\t")
    women = women.filter(pl.col("sample").is_in(header)).sort("sample")
    ids = women["sample"].to_list()
    betas = reprocessed.by_accession(paths, accession, {s: s for s in ids})
    x, feats = complete_matrix(betas, features, ids)
    genes, expr = to_symbols(*read_tpm(paths.input / "arrayexpress" / RNA / TPM, ids), hgnc)
    gene_mean = dict(zip(genes, np.nanmean(expr, axis=1).tolist(), strict=True))
    ages = women["age"].to_numpy()
    cov = dummies(women["batch"].to_numpy())
    chips_file = folder / "chips.csv"
    tables: dict[str, list] = {}
    if chips_file.exists():
        chip_of = dict(pl.read_csv(chips_file, infer_schema_length=0).select("sample", "chip").iter_rows())
        chips = np.array([chip_of[s] for s in ids])
        tables["age_by_chip"] = age_by_chip(ages, chips).to_dicts()
        if chip_terms():
            cov = np.column_stack([cov, dummies(chips)])
    elif chip_terms():
        raise FileNotFoundError(f"{chips_file}: run scripts/idat_barcodes.py first")
    print(
        f"epidermis_expression: {arm}: {len(ids)} women, {len(feats)} CpGs, {len(genes)} genes",
        file=sys.stderr,
        flush=True,
    )
    metrics = _lifespan_metrics(Prepared(f"epidermis ({arm})", x, ages, cov, feats, np.array(ids), cov))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    level, gene = cpg_expression(feats, links, gene_mean)
    bins = expression_bins(level)
    levels = baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
    base = np.char.add(np.char.add(island.astype(str), "|"), levels)
    with_expr = np.char.add(np.char.add(base, "|"), bins)
    variants = {"island x level": base, "+ own expression decile": with_expr}
    state = keratinocyte_state(feats, probes, paths)
    if state is not None:
        variants[f"+ chromatin state ({KERATINOCYTES})"] = np.char.add(np.char.add(base, "|"), state)
        variants["+ expression + state"] = np.char.add(np.char.add(with_expr, "|"), state)
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    in_any = np.any([s for s, _ in sets.values()], axis=0)
    gaps, grad = [], []
    for metric in METRICS:
        values = metrics[metric].to_numpy()
        gaps += [{"metric": metric, **r} for r in matched_gaps(values, sets, variants, N_RESAMPLES[scale], ctx.seed)]
        grad += [{"metric": metric, **r} for r in gradient(values, base, bins, gene, ~in_any)]
    tables["matched_gaps"], tables["gradient"] = gaps, grad
    tables["set_expression"] = [
        {"group": g, "bin": b, "cpgs": int(((bins == b) & s).sum()), "share": float(((bins == b) & s).sum() / s.sum())}
        for g, (s, _) in {**sets, "outside the sets": (~in_any, [])}.items()
        for b in [NO_PROBE, *map(str, range(1, BINS + 1))]
    ]
    row_of = {g: i for i, g in enumerate(genes)}
    by_gene: dict[str, list[int]] = {}
    for j, g in enumerate(gene):
        if g is not None:
            by_gene.setdefault(g, []).append(j)
    used = sorted(by_gene)
    meth = np.vstack([x[:, by_gene[g]].mean(axis=1) for g in used])
    design = np.column_stack([np.ones(len(ids)), ages, cov[:, np.ptp(cov, axis=0) > 1e-12]])
    row = coupling(meth, expr[[row_of[g] for g in used]], design, PERMUTATIONS[scale], ctx.rng(151_000))
    tables["coupling"] = [row]
    tables["sizes"] = [
        {
            "accession": accession,
            "women": len(ids),
            "age_min": float(ages.min()),
            "age_max": float(ages.max()),
            "batches": int(women["batch"].n_unique()),
            "cpgs": len(feats),
            "genes": len(genes),
        }
    ]
    return tables


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    out = ctx.out_dir("epidermis_expression")
    if not (paths.input / "arrayexpress" / RNA / TPM).exists():
        print(f"epidermis_expression: no {TPM} of {RNA}; nothing run", file=sys.stderr)
        return []
    features, probes, membership = _features(paths, ctx.rng(140_000))
    annotation = AnnotationLayout(paths)
    links, hgnc = pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.hgnc)
    collected: dict[str, list] = {}
    for arm, accession in ARMS.items():
        if not reprocessed.available(paths, accession):
            print(f"epidermis_expression: no SeSAMe rebuild for {accession}; arm skipped", file=sys.stderr)
            continue
        for name, rows in one_arm(ctx, arm, accession, features, probes, membership, links, hgnc).items():
            collected.setdefault(name, []).extend({"arm": arm, **r} for r in rows)
    for name, rows in collected.items():
        write_table(pl.DataFrame(rows), out / name)
    return []
