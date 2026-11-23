"""Neurons against oligodendrocytes of the same brains by whole-genome bisulfite sequencing (GSE107729).

Mendizabal et al. 2019 (Genome Biol) sequenced NeuN+ and OLIG2+ nuclei sorted from prefrontal cortex (BA46) of
adults with and without schizophrenia; per CpG they deposit methylated and total reads (hg19). This is a second,
platform-independent test of the glia-against-neuron contrast that the arrays gave (R110: one glial fraction more
shielded than paired neurons; R115: same sign, not resolved).

Each library is reduced to the 450k probes' CpG positions (GRCh37 MAPINFO; offset checked against one library: 301,041
probes at offset 0, none at +-1), covered at least ``MIN_COVERAGE`` times, so the matching, gene sets and strata are
those of every array analysis. Donors with both fractions are kept; CpGs measured in every library of both fractions
are used; age slopes hold sex and diagnosis. As in ``nuclei_paired``: the gap per cell type and the oligodendrocyte
minus neuron difference, with a donor bootstrap that draws donors once and takes both of their libraries, scaled by
the full sample's SD. A second analysis uses the controls alone.

``entropy-age figures wgbs_nuclei`` (unpacks ``GSE107729_RAW.tar`` once; reads ``samples_brief.txt``, the series'
sample records in GEO's text form).
"""

from __future__ import annotations

import os
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.nuclei_paired import N_BOOT, gaps, paired_difference
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE107729"
MIN_COVERAGE = 10
MIN_DONORS = 15
MANIFEST = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"


def read_samples(path: Path) -> pl.DataFrame:
    """``sample``, ``cell_type``, ``donor``, ``age``, ``sex``, ``diagnosis`` from GEO's text records."""
    rows, current = [], {}
    for line in path.read_text().splitlines():
        if line.startswith("^SAMPLE"):
            if current:
                rows.append(current)
            current = {"sample": line.split("=", 1)[1].strip()}
        elif line.startswith("!Sample_characteristics_ch1"):
            key, _, value = line.split("=", 1)[1].strip().partition(":")
            current[key.strip()] = value.strip()
    if current:
        rows.append(current)
    frame = pl.DataFrame(rows)
    return frame.select(
        "sample",
        pl.col("cell type").alias("cell_type"),
        pl.col("subject id").alias("donor"),
        pl.col("age").cast(pl.Float64, strict=False),
        pl.col("Sex").alias("sex"),
        pl.col("disease state").alias("diagnosis"),
    ).drop_nulls("age")


def probe_positions(paths) -> pl.DataFrame:
    """``cpg``, ``chromosome`` (``chr1`` ...) and ``CpG_position`` (GRCh37) of every 450k cg probe."""
    m = pl.read_csv(
        AnnotationLayout(paths).input / MANIFEST,
        skip_rows=7,
        columns=["IlmnID", "CHR", "MAPINFO"],
        infer_schema_length=0,
        truncate_ragged_lines=True,
    )
    return m.filter(pl.col("IlmnID").str.starts_with("cg") & pl.col("MAPINFO").is_not_null()).select(
        pl.col("IlmnID").alias("cpg"),
        ("chr" + pl.col("CHR")).alias("chromosome"),
        pl.col("MAPINFO").cast(pl.Int64).alias("CpG_position"),
    )


def library_betas(path: Path, probes: pl.DataFrame) -> pl.DataFrame:
    """``cpg`` and ``beta`` at the probe positions covered ``MIN_COVERAGE`` times in one library."""
    reads = pl.read_csv(path, separator=" ", schema_overrides={"chromosome": pl.String, "CpG_position": pl.Int64})
    return (
        reads.filter(pl.col("Total_counts") >= MIN_COVERAGE)
        .join(probes, on=["chromosome", "CpG_position"])
        .select("cpg", (pl.col("M_counts") / pl.col("Total_counts")).alias("beta"))
    )


def matrix(ctx: Context, samples: pl.DataFrame) -> pl.DataFrame:
    """``cpg`` x sample betas, cached in ``data/interim/wgbs_gse107729``."""
    cache = ctx.paths.interim / "wgbs_gse107729" / "betas.parquet"
    if cache.exists():
        return pl.read_parquet(cache)
    folder = ctx.paths.input / "geo" / GSE
    raw = folder / "raw"
    with tarfile.open(folder / f"{GSE}_RAW.tar") as archive:
        archive.extractall(raw, filter="data")
    probes = probe_positions(ctx.paths)
    columns = []
    for path in sorted(raw.glob("GSM*_CpG_WGBS.txt.gz")):
        gsm = re.match(r"(GSM\d+)_", path.name).group(1)  # type: ignore[union-attr]
        if gsm in set(samples["sample"].to_list()):
            columns.append(library_betas(path, probes).rename({"beta": gsm}))
            print(f"wgbs_nuclei: {gsm}: {columns[-1].height} probe CpGs", file=sys.stderr, flush=True)
    frame = columns[0]
    for c in columns[1:]:
        frame = frame.join(c, on="cpg", how="full", coalesce=True)
    cache.parent.mkdir(parents=True, exist_ok=True)
    frame.write_parquet(cache)
    return frame


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE
    if not (folder / f"{GSE}_RAW.tar").exists() or not (folder / "samples_brief.txt").exists():
        print(f"wgbs_nuclei: {GSE} archive or samples_brief.txt missing; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("wgbs_nuclei")
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", N_BOOT[str(ctx.scale.name)]))
    samples = read_samples(folder / "samples_brief.txt")
    betas = matrix(ctx, samples)
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    samples = samples.filter(pl.col("sample").is_in(betas.columns))
    rng = ctx.rng(312_000)
    rows, sizes = [], []
    for label, subset in (
        ("every donor (diagnosis held)", samples),
        ("controls", samples.filter(pl.col("diagnosis") == "Control")),
    ):
        per_type = {t: subset.filter(pl.col("cell_type") == t).unique("donor", keep="first") for t in ("Olig2", "NeuN")}
        both = sorted(set(per_type["Olig2"]["donor"].to_list()) & set(per_type["NeuN"]["donor"].to_list()))
        if len(both) < MIN_DONORS:
            print(f"wgbs_nuclei: {label}: {len(both)} donors with both, skipped", file=sys.stderr)
            continue
        o, n = (per_type[t].filter(pl.col("donor").is_in(both)).sort("donor") for t in ("Olig2", "NeuN"))
        cols = [*o["sample"].to_list(), *n["sample"].to_list()]
        complete = betas.select("cpg", *cols).drop_nulls()
        feats = [f for f in complete["cpg"].to_list() if f in set(features)]
        complete = complete.filter(pl.col("cpg").is_in(feats))
        feats = complete["cpg"].to_list()
        x = complete.select(cols).to_numpy().T
        xo, xn = x[: o.height], x[o.height :]
        ages = o["age"].to_numpy()
        cov = np.column_stack([dummies(o["sex"].to_numpy()), dummies(o["diagnosis"].to_numpy())])
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        sets = _sets(feats, membership)
        point = (gaps(xo, ages, cov, feats, island, sets), gaps(xn, ages, cov, feats, island, sets))
        boots = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(both), len(both))
            if np.ptp(ages[idx]) < 20:
                continue
            boots.append(
                (
                    gaps(xo[idx], ages[idx], cov[idx], feats, island, sets),
                    gaps(xn[idx], ages[idx], cov[idx], feats, island, sets),
                )
            )
        for r in paired_difference(point[0], point[1], boots):
            renamed = {k.replace("glia", "oligodendrocytes"): v for k, v in r.items()}
            rows.append({"donors_used": label, **renamed})
        sizes.append(
            {
                "donors_used": label,
                "donors": len(both),
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
            }
        )
        print(f"wgbs_nuclei: {label}: {len(both)} donors, {len(feats)} CpGs", file=sys.stderr, flush=True)
    if rows:
        write_table(pl.DataFrame(sizes), out / "sizes")
        write_table(pl.DataFrame(rows), out / "paired")
    return []
