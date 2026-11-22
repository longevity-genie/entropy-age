"""Many tissues of the same mice at three ages, by RRBS on lanes that mix the ages (roadmap RM40, GSE225166).

GSE225166 (found by the external dataset search) sequenced up to 14 tissues and sorted cell types of 43 mice aged 3
to 24 months. Each sequencing lane holds libraries of several ages, read from the file names
(``<GSM>_lane<k>_<barcode>_<mouse>_<age>_<sex>_<tissue>_...``), so lane can be held and cannot make an age effect.
Liver, heart, cortex and lung have 10 to 12 mice at each of 3, 15 and 24 months; the other tissues have few mice per
age and are reported only when they reach ``MIN_MICE``.

As for the other RRBS series: CpGs in promoters (``data/rrbs.py``) covered at least 10 times, kept when present in
90 % of a tissue's libraries. Per tissue: the matched gap of each gene set on the age slopes of relative entropy,
entropy and dispersion (sex and lane held; strata of island status x young-level decile, young = the youngest age
of the tissue), with a mouse bootstrap.

``entropy-age figures mouse_atlas_rrbs`` as a Slurm job (unpacks ``GSE225166_RAW.tar`` once).
"""

from __future__ import annotations

import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.rrbs import MIN_COVERAGE, TRACKS, RrbsLayout, in_intervals, islands, promoters
from entropy_age.experiments.common import Context, dummies, scaled_intervals, write_table
from entropy_age.experiments.interventions import membership
from entropy_age.experiments.jena_rrbs import CHROMOSOMES, present_enough
from entropy_age.experiments.mouse_rrbs_tissues import tissue_gaps

GSE = "GSE225166"
NAME = re.compile(
    r"(GSM\d+)_lane(\d+)_[ACGTN]+_(M[0-9A-Za-z-]+)_(\d+)_([MF])_(.+?)(?:_plate\d+_pool\d+)?_i(?:PCRtag|TAG|Tag)"
)
MIN_MICE = 12  # libraries a tissue needs, with at least MIN_PER_END at its youngest and at its oldest age
MIN_PER_END = 3
N_BOOT = {"paper": 200, "quick": 10}


def parse_names(names: list[str]) -> pl.DataFrame:
    """``sample`` (GSM), ``file``, ``lane``, ``mouse``, ``age_months``, ``sex`` and ``tissue`` from the file names;
    names that do not follow the pattern are left out."""
    rows = []
    for name in names:
        found = NAME.search(Path(name).name)
        if found is None:
            continue
        gsm, lane, mouse, age, sex, tissue = found.groups()
        rows.append(
            {
                "sample": gsm,
                "file": name,
                "lane": f"lane{lane}",
                "mouse": mouse,
                "age_months": float(age),
                "sex": sex,
                "tissue": tissue,
            }
        )
    return pl.DataFrame(rows)


def read_library(path: Path, prom: pl.DataFrame) -> pl.DataFrame:
    """``cpg`` (``chr:pos``) and ``beta`` of one Bismark coverage file, promoter CpGs covered ``MIN_COVERAGE`` times."""
    cov = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chr", "pos", "end", "pct", "meth", "unmeth"],
        schema_overrides={"chr": pl.String, "pct": pl.Float64},
    ).filter((pl.col("meth") + pl.col("unmeth")) >= MIN_COVERAGE)
    chrom = cov["chr"].str.strip_prefix("chr")
    c = chrom.to_numpy()
    keep = np.isin(c, CHROMOSOMES) & in_intervals(c, cov["pos"].to_numpy() - 1, prom)
    return cov.filter(pl.Series(keep)).select(
        (pl.col("chr").str.strip_prefix("chr") + ":" + pl.col("pos").cast(pl.String)).alias("cpg"),
        (pl.col("meth") / (pl.col("meth") + pl.col("unmeth"))).alias("beta"),
    )


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE
    tar = folder / f"{GSE}_RAW.tar"
    rr = RrbsLayout(ctx.paths)
    if not tar.exists() or any(not rr.track(t).exists() for t in TRACKS):
        print(f"mouse_atlas_rrbs: {tar} or the mm10 tracks are missing; nothing run", file=sys.stderr)
        return []
    raw = folder / f"{GSE}_RAW"
    if not raw.exists():
        with tarfile.open(tar) as archive:
            archive.extractall(raw, filter="data")
    out = ctx.out_dir("mouse_atlas_rrbs")
    prom = promoters(rr.track(TRACKS[0]), rr.track(TRACKS[1]))
    isl = islands(rr.track(TRACKS[2]))
    samples_all = parse_names([p.name for p in sorted(raw.glob("*.cov.gz"))])
    rng = ctx.rng(340_000)
    gaps, sizes = [], []
    for (tissue,), samples in samples_all.group_by("tissue", maintain_order=True):
        samples = samples.unique(subset="mouse", keep="first", maintain_order=True).sort("age_months", "sample")
        ages = samples["age_months"].to_numpy()
        ends_ok = (ages == ages.min()).sum() >= MIN_PER_END and (ages == ages.max()).sum() >= MIN_PER_END
        if samples.height < MIN_MICE or ages.min() == ages.max() or not ends_ok:
            sizes.append({"tissue": tissue, "mice": samples.height, "used": False})
            continue
        long = pl.concat(
            [
                read_library(raw / f, prom).with_columns(sample=pl.lit(s))
                for s, f in zip(samples["sample"].to_list(), samples["file"].to_list(), strict=True)
            ]
        )
        betas = present_enough(long.pivot(on="sample", index="cpg", values="beta"))
        m = FeatureMatrix.from_frame(betas, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=0.1)
        chrom = pl.Series([f.split(":")[0] for f in m.features])
        pos = pl.Series([int(f.split(":")[1]) for f in m.features])
        cpgs = pl.DataFrame({"cpg": m.features, "chr": chrom, "pos": pos}).with_columns(
            cpg_island=pl.Series(in_intervals(chrom.to_numpy(), pos.to_numpy() - 1, isl))
        )
        members = membership(ctx.paths, m.features, cpgs)
        x = m.values
        young = ages == ages.min()
        island = cpgs["cpg_island"].to_numpy()
        strata = np.char.add(island.astype(str), baseline_strata_of(x[young]))
        cov = np.column_stack([dummies(samples["sex"].to_numpy()), dummies(samples["lane"].to_numpy())])
        lanes = samples.group_by("lane").agg(pl.col("age_months").n_unique().alias("ages"))
        sizes.append(
            {
                "tissue": tissue,
                "mice": x.shape[0],
                "used": True,
                "ages": ",".join(f"{a:g}x{int((ages == a).sum())}" for a in np.unique(ages)),
                "cpgs": len(m.features),
                "lanes": lanes.height,
                "lanes_with_2_or_more_ages": int((lanes["ages"] >= 2).sum()),
            }
        )
        print(
            f"mouse_atlas_rrbs: {tissue}: {x.shape[0]} mice, {len(m.features)} promoter CpGs",
            file=sys.stderr,
            flush=True,
        )
        point = pl.DataFrame(tissue_gaps(x, ages, cov, m.features, strata, members))
        boot = []
        for k in range(N_BOOT[str(ctx.scale.name)]):
            rows = np.sort(rng.choice(x.shape[0], size=x.shape[0], replace=True))
            if (ages[rows] == ages.min()).sum() < 2 or np.unique(ages[rows]).size < 2:
                continue
            boot += [{"boot": k, **r} for r in tissue_gaps(x[rows], ages[rows], cov[rows], m.features, strata, members)]
        gaps.append(
            scaled_intervals(point, pl.DataFrame(boot), ["metric", "group"]).with_columns(tissue=pl.lit(tissue))
        )
    write_table(pl.DataFrame(sizes, infer_schema_length=None), out / "sizes")
    if gaps:
        write_table(pl.concat(gaps), out / "gaps")
    return []


def baseline_strata_of(young: np.ndarray) -> np.ndarray:
    """Decile of the young level per CpG, as strings (CpGs missing in every young library fall in one class)."""
    from entropy_age.stochasticity import baseline_strata

    level = np.nanmean(young, axis=0)
    return baseline_strata(np.nan_to_num(level, nan=0.0), 10).astype(str)
