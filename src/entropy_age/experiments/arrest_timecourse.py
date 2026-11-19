"""Drift and the shield in cells that do not divide: fibroblasts held confluent for six months (roadmap RM40).

E-MTAB-15859 (found by the external dataset search) sampled dermal fibroblast cultures of four donors every month for
six months while they were kept confluent and non-proliferating, on the EPIC array. R111 could not compare division
with arrest because its arrested cultures showed no drift; here time passes without division, so two questions can
be asked of each donor: does methylation drift with time at all, and if it does, do essential and housekeeping
promoters drift less than matched CpGs?

Used: the untreated ("Control") fibroblast cultures on EPIC v1 of the three donors whose months sit on one chip
(donor 1: the six arrays of its main chip; donors 3 and 4: all six). Donor 2 is on EPIC v2 with months split over
three chips, so time is confounded with chip there and it is left out. Limits: six arrays per donor and no
replicates; in donors 3 and 4 the month runs in order down the chip, so time is tied to array position (donor 1's
positions are scrambled); the drug arms sit on another array version and are not analysed.

Per donor: the slope of each promoter CpG on month. (1) Drift: the mean squared slope over CpGs against its
distribution under every reordering of the months (720 for six arrays), which array noise alone would give.
(2) Shield: the matched gap of the absolute slope and of the relative-entropy slope (reference: month 0), set minus
CpGs of the same island status and starting-level decile, with gene-resampling intervals.

``entropy-age figures arrest_timecourse`` (needs the SeSAMe rebuild of E-MTAB-15859).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.stochasticity import baseline_strata, relative_entropy

ACCESSION = "E-MTAB-15859"
DONORS = ("Donor #1", "Donor #3", "Donor #4")
MIN_ARRAYS = 5


def read_samples(sdrf: Path) -> pl.DataFrame:
    """``sample`` (chip_position), ``donor``, ``month``, ``chip``, ``position``, ``line``, ``compound`` per array."""
    s = pl.read_csv(sdrf, separator="\t", infer_schema_length=0)
    return (
        s.select(
            pl.col("Characteristics[cell line]").alias("line"),
            pl.col("Characteristics[individual]").alias("donor"),
            pl.col("Array Data File").str.replace(r"_(Grn|Red)\.idat$", "").alias("sample"),
            pl.col("Factor Value[compound]").alias("compound"),
            pl.col("Factor Value[time]").str.extract(r"(\d+)").cast(pl.Int64).alias("month"),
        )
        .unique(subset="sample")
        .with_columns(
            chip=pl.col("sample").str.extract(r"^(\d+)_"), position=pl.col("sample").str.extract(r"_(R\d+C\d+)")
        )
    )


def slopes(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Least-squares slope of every column of ``x`` on ``t``."""
    c = t - t.mean()
    return c @ (x - x.mean(axis=0)) / (c @ c)


def drift_test(x: np.ndarray, t: np.ndarray) -> dict:
    """Mean squared slope over CpGs, and its place among every reordering of ``t`` (exact; n <= 8 arrays)."""
    observed = float((slopes(x, t) ** 2).mean())
    null = np.array([float((slopes(x, np.array(p)) ** 2).mean()) for p in itertools.permutations(t)])
    return {
        "mean_sq_slope": observed,
        "null_median": float(np.median(null)),
        "ratio_to_null": observed / float(np.median(null)),
        "p_reorder": float((null >= observed - 1e-18).mean()),
        "reorderings": int(null.size),
    }


def run(ctx: Context) -> list[Path]:
    sdrf = ctx.paths.input / "arrayexpress" / ACCESSION / f"{ACCESSION}.sdrf.txt"
    if not (reprocessed.available(ctx.paths, ACCESSION) and sdrf.exists()):
        print(f"arrest_timecourse: no SeSAMe rebuild or SDRF for {ACCESSION}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("arrest_timecourse")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, ACCESSION).filter(pl.col("passed") & (pl.col("platform") == "EPIC"))
    samples = read_samples(sdrf).filter(
        (pl.col("line") == "dermal fibroblasts")
        & (pl.col("compound") == "Control")
        & pl.col("donor").is_in(DONORS)
        & pl.col("sample").is_in(qc["sample"].implode())
    )
    betas = reprocessed.by_accession(ctx.paths, ACCESSION, {s: s for s in samples["sample"].to_list()})
    drift_rows, gap_rows, design = [], [], []
    for donor in DONORS:
        d = samples.filter(pl.col("donor") == donor)
        main = d.group_by("chip").len().sort("len", descending=True)["chip"][0]
        d = d.filter(pl.col("chip") == main).sort("month")  # one chip only: time is never a chip contrast
        if d.height < MIN_ARRAYS:
            continue
        x, feats = complete_matrix(betas, features, d["sample"].to_list())
        t = d["month"].to_numpy().astype(float)
        rows = np.array([int(p[1:3]) for p in d["position"].to_list()], dtype=float)
        design.append(
            {
                "donor": donor,
                "arrays": d.height,
                "months": ",".join(str(int(m)) for m in t),
                "chip": main,
                "cpgs": len(feats),
                "r_month_row": float(np.corrcoef(t, rows)[0, 1]),
            }
        )
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        base = x[t == t.min()].mean(axis=0)
        sets = _sets(feats, membership)
        strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(base, 10).astype(str))
        drift_rows.append({"donor": donor, **drift_test(x, t)})
        slope = slopes(x, t)
        kl_slope = slopes(relative_entropy(np.clip(x, 0.0, 1.0), base), t)
        for metric, values in (("abs_slope", np.abs(slope)), ("kl_slope", kl_slope), ("signed_slope", slope)):
            gap_rows += [
                {"donor": donor, "metric": metric, **r}
                for r in matched_gaps(values, sets, {"island x level": strata}, n_resamples, seed=ctx.seed)
            ]
        print(f"arrest_timecourse: {donor}: {d.height} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True)
    write_table(pl.DataFrame(design), out / "design")
    write_table(pl.DataFrame(drift_rows), out / "drift")
    write_table(pl.DataFrame(gap_rows), out / "gaps")
    return []
