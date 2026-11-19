"""A renewing and a post-mitotic tissue of one mouse cohort by RRBS: colon and hippocampus (roadmap RM40).

The mouse-array turnover ordering was withdrawn because age sat on different chips (R92). The external dataset search
found the one public pair that can carry the question by sequencing: C57BL/6J/Ukj males aged 3 to 28 months,
colon (GSE233734, 82 libraries in the main set) and hippocampus (GSE304754), 16-18 mice per age, RRBS on a NovaSeq.
For the colon series every one of its five flowcells holds every age bin (R2 of age by flowcell 0.02;
``data/input/geo/audit/GSE233734``); the hippocampus series was audited by the search's track only.

Both are processed per-CpG tables: colon as methylation rates with a coverage matrix, hippocampus as methylated and
unmethylated read counts per mouse. As for the other RRBS series, CpGs in promoters (``data/rrbs.py``) covered at
least 10 times are used, kept when present in 90 % of a tissue's mice. Per tissue: the matched gap of each gene set
(age slope of relative entropy, entropy and dispersion; strata of island status x young-level decile, young = the
3-month mice), with a mouse bootstrap. The contrast is colon minus hippocampus, from independent resamples: the
records do not link a colon library to a hippocampus library of the same mouse (the library names look alike, but
that is unverified), so pairing is not assumed.

``entropy-age figures jena_rrbs`` as a Slurm job (reads about 1.4 GB of compressed tables).
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.rrbs import MIN_COVERAGE, MIN_PRESENT, TRACKS, RrbsLayout, in_intervals, islands, promoters
from entropy_age.experiments.common import Context, scaled_intervals, write_table
from entropy_age.experiments.interventions import membership
from entropy_age.experiments.mouse_rrbs_tissues import tissue_gaps
from entropy_age.stochasticity import baseline_strata

COLON, HIPPOCAMPUS = "GSE233734", "GSE304754"
YOUNG_MONTHS = 3.0
N_BOOT = {"paper": 200, "quick": 20}
CHROMOSOMES = [str(k) for k in range(1, 20)]


def _promoter_rows(chrom: pl.Series, pos0: pl.Series, prom: pl.DataFrame) -> np.ndarray:
    c = chrom.str.strip_prefix("chr").to_numpy()
    return np.isin(c, CHROMOSOMES) & in_intervals(c, pos0.to_numpy(), prom)


def read_colon(rates: Path, coverage: Path, prom: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Promoter CpGs of the colon set as ``cpg`` (``chr:pos``, 1-based) plus one beta column per library (null below
    ``MIN_COVERAGE``), and ``sample``, ``age_months`` read from the library names (``15M.N3``)."""
    opts = {"separator": "\t", "null_values": ["NA"], "infer_schema_length": 10_000}
    r = pl.read_csv(rates, **opts)
    c = pl.read_csv(coverage, **opts)
    if r.columns != c.columns or r.height != c.height:
        raise ValueError("colon rates and coverage tables do not line up")
    keep = _promoter_rows(r["chrom"], r["start"], prom)
    r, c = r.filter(pl.Series(keep)), c.filter(pl.Series(keep))
    names = r.columns[3:]
    cpg = (r["chrom"].str.strip_prefix("chr") + ":" + (r["start"] + 1).cast(pl.String)).alias("cpg")
    joined = pl.concat(
        [
            r.select(pl.col(n).cast(pl.Float64) for n in names),
            c.select(pl.col(n).cast(pl.Float64).alias(f"cov::{n}") for n in names),
        ],
        how="horizontal",
    )
    betas = joined.select(pl.when(pl.col(f"cov::{n}") >= MIN_COVERAGE).then(pl.col(n)).alias(n) for n in names)
    betas = betas.insert_column(0, cpg)
    samples = pl.DataFrame(
        {"sample": names, "age_months": [float(re.match(r"(\d+)M", n).group(1)) for n in names]}  # type: ignore[union-attr]
    )
    return betas, samples


def read_hippocampus(counts: Path, matrix: Path, prom: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The same for hippocampus, from per-mouse ``<id>_5mC`` and ``<id>_C`` read counts; ages from the GEO record,
    whose description line carries the mouse id (``3M_4, J14208``)."""
    t = pl.read_csv(counts, separator="\t", infer_schema_length=10_000, truncate_ragged_lines=True)
    chrom, pos = t["Chromosome"].cast(pl.String), t["Position"].cast(pl.Int64)
    keep = _promoter_rows(chrom, pos - 1, prom)
    t, chrom, pos = t.filter(pl.Series(keep)), chrom.filter(pl.Series(keep)), pos.filter(pl.Series(keep))
    ids = [c.removesuffix("_5mC") for c in t.columns if c.endswith("_5mC")]
    counts = t.select(pl.col(f"{i}_5mC", f"{i}_C").cast(pl.Float64) for i in ids)
    betas = counts.select(
        pl.when((pl.col(f"{i}_5mC") + pl.col(f"{i}_C")) >= MIN_COVERAGE)
        .then(pl.col(f"{i}_5mC") / (pl.col(f"{i}_5mC") + pl.col(f"{i}_C")))
        .alias(i)
        for i in ids
    )
    betas = betas.insert_column(0, (chrom.str.strip_prefix("chr") + ":" + pos.cast(pl.String)).alias("cpg"))
    ages, descriptions = [], []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1") and cells and cells[0].startswith("age:"):
                ages = [float(re.search(r"(\d+)", c).group(1)) for c in cells]  # type: ignore[union-attr]
            elif line.startswith("!Sample_description") and cells and re.search(r"J\d+", cells[0]):
                descriptions = [re.search(r"J\d+", c).group(0) for c in cells]  # type: ignore[union-attr]
    age_of = dict(zip(descriptions, ages, strict=True))
    samples = pl.DataFrame({"sample": ids, "age_months": [age_of.get(i) for i in ids]}).drop_nulls("age_months")
    return betas.select("cpg", *samples["sample"].to_list()), samples


def present_enough(betas: pl.DataFrame) -> pl.DataFrame:
    """CpGs measured in at least ``MIN_PRESENT`` of the libraries."""
    n = len(betas.columns) - 1
    present = pl.sum_horizontal(pl.col(c).is_not_null().cast(pl.Int64) for c in betas.columns[1:])
    return betas.filter(present >= MIN_PRESENT * n)


def run(ctx: Context) -> list[Path]:
    geo = ctx.paths.input / "geo"
    files = {
        "colon": (geo / COLON / f"{COLON}_colon_CpG_methylationrates.bedGraph.gz",
                  geo / COLON / f"{COLON}_colon_CpG_coverages.bedGraph.gz"),
        "hippocampus": (geo / HIPPOCAMPUS / f"{HIPPOCAMPUS}_hippocampus_combined.cov.gz",
                        geo / HIPPOCAMPUS / f"{HIPPOCAMPUS}_series_matrix.txt.gz"),
    }  # fmt: skip
    rr = RrbsLayout(ctx.paths)
    missing = [str(p) for pair in files.values() for p in pair if not p.exists()]
    missing += [str(rr.track(t)) for t in TRACKS if not rr.track(t).exists()]
    if missing:
        print(f"jena_rrbs: inputs missing ({', '.join(missing)}); nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("jena_rrbs")
    prom = promoters(rr.track(TRACKS[0]), rr.track(TRACKS[1]))
    isl = islands(rr.track(TRACKS[2]))
    readers = {"colon": read_colon, "hippocampus": read_hippocampus}
    turnover = {"colon": "renewing", "hippocampus": "post-mitotic"}
    rng = ctx.rng(330_000)
    gaps, boots, sizes, measured = [], {}, [], {}
    for tissue, (a, b) in files.items():
        betas, samples = readers[tissue](a, b, prom)
        betas = present_enough(betas)
        samples = samples.sort("age_months", "sample")
        m = FeatureMatrix.from_frame(betas, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.1)
        chrom = pl.Series([f.split(":")[0] for f in m.features])
        pos = pl.Series([int(f.split(":")[1]) for f in m.features])
        cpgs = pl.DataFrame({"cpg": m.features, "chr": chrom, "pos": pos}).with_columns(
            cpg_island=pl.Series(in_intervals(chrom.to_numpy(), pos.to_numpy() - 1, isl))
        )
        members = membership(ctx.paths, m.features, cpgs)
        measured[tissue] = set(m.features)
        x, ages = m.values, samples["age_months"].to_numpy()
        young = ages <= YOUNG_MONTHS
        island = cpgs["cpg_island"].to_numpy()
        strata = np.char.add(island.astype(str), baseline_strata(x[young].mean(axis=0), 10).astype(str))
        cov = np.zeros((x.shape[0], 0))
        sizes.append(
            {
                "tissue": tissue,
                "mice": x.shape[0],
                "young": int(young.sum()),
                "ages": ",".join(f"{a:g}x{int((ages == a).sum())}" for a in np.unique(ages)),
                "cpgs": len(m.features),
                "island_share": float(island.mean()),
            }
        )
        print(f"jena_rrbs: {tissue}: {x.shape[0]} mice, {len(m.features)} promoter CpGs", file=sys.stderr, flush=True)
        point = pl.DataFrame(tissue_gaps(x, ages, cov, m.features, strata, members))
        boot = []
        for k in range(N_BOOT[str(ctx.scale.name)]):
            rows = np.sort(rng.choice(x.shape[0], size=x.shape[0], replace=True))
            if (ages[rows] <= YOUNG_MONTHS).sum() < 3:
                continue
            boot += [{"boot": k, **r} for r in tissue_gaps(x[rows], ages[rows], cov[rows], m.features, strata, members)]
        boot = pl.DataFrame(boot)
        gaps.append(
            scaled_intervals(point, boot, ["metric", "group"]).with_columns(
                tissue=pl.lit(tissue), turnover=pl.lit(turnover[tissue])
            )
        )
        boots[tissue] = boot.join(
            point.select("metric", "group", pl.col("sd").alias("sd_full")), on=["metric", "group"]
        )
    gaps = pl.concat(gaps)
    # a coordinate check: both tissues are RRBS of one strain, so most promoter CpGs should be shared; an off-by-one
    # in either table would show as a larger overlap after shifting hippocampus positions by one
    shifted = {f"{c}:{int(q) + s}" for c, q in (f.split(":") for f in measured["hippocampus"]) for s in (-1, 1)}
    print(
        f"jena_rrbs: CpGs shared {len(measured['colon'] & measured['hippocampus'])}, "
        f"shared after a one-base shift {len(measured['colon'] & shifted)}",
        file=sys.stderr,
        flush=True,
    )
    contrast = []
    for (metric, group), c in boots["colon"].group_by("metric", "group"):
        h = boots["hippocampus"].filter((pl.col("metric") == metric) & (pl.col("group") == group))
        dc = (c["matched_diff"] / c["sd_full"]).to_numpy()
        dh = (h["matched_diff"] / h["sd_full"]).to_numpy()
        n = min(dc.size, dh.size)
        diff = dc[:n] - dh[:n]  # independent resamples of the two tissues
        point = gaps.filter((pl.col("metric") == metric) & (pl.col("group") == group))
        est = dict(zip(point["tissue"], point["gap_sd"], strict=True))
        tail = min((diff <= 0).mean(), (diff >= 0).mean())
        contrast.append(
            {
                "metric": metric,
                "group": group,
                "colon": est["colon"],
                "hippocampus": est["hippocampus"],
                "difference": est["colon"] - est["hippocampus"],
                "ci_low": float(np.quantile(diff, 0.025)),
                "ci_high": float(np.quantile(diff, 0.975)),
                "p_bootstrap": float(min(1.0, max(2 * tail, 1 / n))),
            }
        )
    write_table(pl.DataFrame(sizes), out / "sizes")
    write_table(gaps, out / "gaps")
    write_table(pl.DataFrame(contrast).sort("metric", "group"), out / "contrast")
    return []
