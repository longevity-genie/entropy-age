"""One fibroblast line with and without telomerase: is the shield per doubling tied to approaching senescence?

GSE268145 (found by the external dataset search, RM40) passaged BJ foreskin fibroblasts (one male line) along two
trajectories, A and B, as mortal cells and with hTERT, on the EPIC array; the series is rebuilt from IDATs. Mortal
cells divide towards replicative senescence; hTERT cells go on dividing without it. If the shield (essential
promoters drifting less than matched CpGs) comes with division itself, it should appear per doubling in both; if it
comes with the approach to senescence, hTERT cells should lack it.

Three groups, each analysed alone: hTERT cells; mortal cells hybridised on the same chips as the hTERT cells; mortal
cells of an earlier batch on other chips. The two deeply senescent arrays (doubling 104) are left out. Per promoter
CpG, slopes of relative entropy, entropy and dispersion on doublings with a trajectory intercept (and chip
intercepts under ``ENTROPY_AGE_CHIP=1``); then set CpGs against CpGs matched on island status and starting level.
The hTERT and same-chip mortal groups share chips, so their slopes are comparable; the overlap of doublings with
chip is tabulated per group.

Limits: one cell line; no arrested control; within a trajectory about half of the spread of doublings lies between
chips. ``entropy-age figures bj_telomerase`` (needs the SeSAMe rebuild of GSE268145).
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, chip_terms, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.sturm_expression import chip_overlap
from entropy_age.experiments.tissue_contrast import Prepared, _features, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

GSE = "GSE268145"
METRICS = ("kl_slope", "entropy_slope", "dispersion_slope")
TITLE = re.compile(r"Cell_BJF_(hTERT_)?Trajectory_(A|B|DeepSen)_PDL_([0-9.]+)")
CHIP = re.compile(r"_(\d{9,})_(R\d\dC\d\d)_Grn")
HTERT, MORTAL_SAME, MORTAL_EARLIER = "hTERT", "mortal, chips shared with hTERT", "mortal, earlier batch"
MIN_ARRAYS = 12


def read_cultures(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM), ``telomerase``, ``trajectory``, ``doublings`` and ``chip`` from the series matrix."""
    fields: dict[str, list[str]] = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            key, _, rest = line.rstrip("\n").partition("\t")
            if key in ("!Sample_title", "!Sample_geo_accession") or (
                key == "!Sample_supplementary_file" and "_Grn" in rest
            ):
                fields[key] = [v.strip('"') for v in rest.split("\t")]
    rows = []
    for gsm, title, file in zip(
        fields["!Sample_geo_accession"], fields["!Sample_title"], fields["!Sample_supplementary_file"], strict=True
    ):
        t, c = TITLE.search(title), CHIP.search(file)
        if t is None or c is None:
            continue
        rows.append(
            {
                "sample": gsm,
                "telomerase": t.group(1) is not None,
                "trajectory": t.group(2),
                "doublings": float(t.group(3)),
                "chip": c.group(1),
            }
        )
    return pl.DataFrame(rows)


def groups(cultures: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """The three groups: hTERT, mortal on chips that also carry hTERT arrays, mortal on other chips; deeply
    senescent arrays (trajectory ``DeepSen``) left out."""
    kept = cultures.filter(pl.col("trajectory").is_in(["A", "B"]))
    htert = kept.filter(pl.col("telomerase"))
    shared = set(htert["chip"].to_list())
    mortal = kept.filter(~pl.col("telomerase"))
    return {
        HTERT: htert,
        MORTAL_SAME: mortal.filter(pl.col("chip").is_in(list(shared))),
        MORTAL_EARLIER: mortal.filter(~pl.col("chip").is_in(list(shared))),
    }


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    matrix = paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(paths, GSE) and matrix.exists()):
        print(f"bj_telomerase: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("bj_telomerase")
    scale = str(ctx.scale.name)
    features, probes, membership = _features(paths, ctx.rng(140_000))
    passed = reprocessed.read_qc(paths, GSE).filter(pl.col("passed"))["sample"].to_list()
    cultures = read_cultures(matrix).filter(pl.col("sample").is_in(passed))
    betas = reprocessed.by_accession(paths, GSE, {s: s for s in cultures["sample"].to_list()})
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    gaps, overlap, sizes, drift = [], [], [], []
    for name, group in groups(cultures).items():
        group = group.sort("trajectory", "doublings")
        if group.height < MIN_ARRAYS:
            print(f"bj_telomerase: {name}: {group.height} arrays, skipped", file=sys.stderr)
            continue
        x, feats = complete_matrix(betas, features, group["sample"].to_list())
        doublings = group["doublings"].to_numpy()
        trajectory, chips = group["trajectory"].to_numpy(), group["chip"].to_numpy()
        overlap += [{"cells": name, **r} for r in chip_overlap(doublings, trajectory, chips).to_dicts()]
        cov = dummies(trajectory)
        if chip_terms():
            cov = np.column_stack([cov, dummies(chips)])
        metrics = _lifespan_metrics(Prepared(name, x, doublings, cov, feats, trajectory, cov))
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        levels = baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
        variants = {"island x level": np.char.add(np.char.add(island.astype(str), "|"), levels)}
        sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
        for metric in METRICS:
            values = metrics[metric].to_numpy()
            gaps += [
                {"cells": name, "metric": metric, **r}
                for r in matched_gaps(values, sets, variants, N_RESAMPLES[scale], ctx.seed)
            ]
            drift.append({"cells": name, "metric": metric, "mean": float(values.mean()), "sd": float(values.std())})
        sizes.append(
            {
                "cells": name,
                "arrays": group.height,
                "chips": int(np.unique(chips).size),
                "doublings_min": float(doublings.min()),
                "doublings_max": float(doublings.max()),
                "cpgs": len(feats),
            }
        )
        print(f"bj_telomerase: {name}: {group.height} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True)
    for name, rows in (("matched_gaps", gaps), ("chip_overlap", overlap), ("sizes", sizes), ("drift", drift)):
        if rows:
            write_table(pl.DataFrame(rows), out / name)
    return []
