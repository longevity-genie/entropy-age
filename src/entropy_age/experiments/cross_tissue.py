"""Four tissues of the same adults on shared chips: does the shield differ between tissues within a person?

GSE280465 (Apsley et al. 2025; found by the external dataset search, RM40) took buccal swabs, saliva, dried blood
spots and blood mononuclear cells from the same adults (19 to 70 years) and ran them on the EPIC v2 array with the
tissues of a person spread over the same chips (17 of 27 chips carry all four). It is rebuilt from IDATs. Unlike
the blood-and-brain series of Lunnon, where slide goes with tissue, a tissue contrast is identifiable here.

Per pair of tissues, the donors who have both are kept and the CpGs measured in every array of both are used. Two
measures of the gap (set CpGs minus CpGs matched on island status and that tissue's own level decile):

* ``dispersion``: between-person spread at each CpG with age and sex regressed out (and chip under
  ``ENTROPY_AGE_CHIP=1``), as in the Lunnon analysis; needs no age slope;
* the age slopes of relative entropy, entropy and absolute level (``nuclei_paired.gaps``), sex held.

The statistic is the gap in the first tissue minus the gap in the second, with an interval from a donor bootstrap
that draws donors once and takes both of each donor's arrays, each resample scaled by the full sample's SD.

Limits: all four tissues renew; saliva and buccal swabs are mixtures of epithelial cells and leukocytes whose
proportions differ between people and are not held; about 40 donors per pair, so age slopes are imprecise.
``entropy-age figures cross_tissue`` (needs the SeSAMe rebuild of GSE280465).
"""

from __future__ import annotations

import gzip
import itertools
import os
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.common import Context, chip_terms, dummies, write_table
from entropy_age.experiments.nuclei_paired import gaps as slope_gaps
from entropy_age.experiments.nuclei_paired import paired_difference
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, dispersion_levels

GSE = "GSE280465"
TISSUES = ("Buccal", "Saliva", "DBS", "PBMC")
N_BOOT = {"paper": 400, "quick": 30}
MIN_DONORS = 25


def read_people(matrix: Path) -> pl.DataFrame:
    """``sample``, ``donor`` (from the title ``<donor>_<tissue>``), ``tissue``, ``age`` and ``sex``."""
    with gzip.open(matrix, "rt") as handle:
        titles = next(
            [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            for line in handle
            if line.startswith("!Sample_title")
        )
    chars = _characteristics(matrix)
    chars = chars.rename({c: c.lower() for c in chars.columns})
    return chars.select(
        "sample",
        pl.Series("donor", [t.rsplit("_", 1)[0] for t in titles]),
        "tissue",
        pl.col("age").cast(pl.Float64, strict=False),
        pl.col("sex").fill_null("NA"),
    ).drop_nulls("age")


def dispersion_gaps(x: np.ndarray, cov: np.ndarray, feats: list[str], island: np.ndarray, sets: dict) -> dict:
    """``("dispersion", group)`` -> (unscaled matched difference, SD of dispersion over CpGs) for one tissue."""
    cov = cov[:, np.ptp(cov, axis=0) > 1e-12]
    table = dispersion_levels(x, feats, covariates=cov if cov.shape[1] else None)
    values = table["dispersion"].to_numpy()
    strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(x.mean(axis=0), 10).astype(str))
    return {
        ("dispersion", group): (float(matched_scores(values[None, :], in_set, strata)[0][0]), float(values.std()))
        for group, (in_set, _) in sets.items()
    }


def both_measures(x, ages, sex, technical, feats, island, sets) -> dict:
    """The dispersion gap (age, sex and ``technical`` held) and the age-slope gaps (sex held) of one tissue."""
    out = dispersion_gaps(x, np.column_stack([ages, sex, technical]), feats, island, sets)
    out.update(slope_gaps(x, ages, sex, feats, island, sets))
    return out


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    matrix = paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(paths, GSE) and matrix.exists()):
        print(f"cross_tissue: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("cross_tissue")
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", N_BOOT[str(ctx.scale.name)]))
    features, _, membership = _features(paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(paths, GSE).filter(pl.col("passed"))
    people = read_people(matrix).join(qc.select("sample", "chip"), on="sample")
    betas = reprocessed.by_accession(paths, GSE, {s: s for s in people["sample"].to_list()})
    people = people.filter(pl.col("sample").is_in(betas.columns))
    rng = ctx.rng(320_000)
    rows, sizes = [], []
    for a, b in itertools.combinations(TISSUES, 2):
        pa = people.filter(pl.col("tissue") == a).unique(subset="donor", keep="first", maintain_order=True)
        pb = people.filter(pl.col("tissue") == b).unique(subset="donor", keep="first", maintain_order=True)
        both = sorted(set(pa["donor"].to_list()) & set(pb["donor"].to_list()))
        if len(both) < MIN_DONORS:
            print(f"cross_tissue: {a} and {b}: {len(both)} donors with both, skipped", file=sys.stderr)
            continue
        pa, pb = (p.filter(pl.col("donor").is_in(both)).sort("donor") for p in (pa, pb))
        x_all, feats = complete_matrix(betas, features, [*pa["sample"].to_list(), *pb["sample"].to_list()])
        xa, xb = x_all[: pa.height], x_all[pa.height :]
        ages = pa["age"].to_numpy()
        sex = dummies(pa["sex"].to_numpy())
        empty = np.empty((len(both), 0))
        ta = dummies(pa["chip"].to_numpy()) if chip_terms() else empty
        tb = dummies(pb["chip"].to_numpy()) if chip_terms() else empty
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        sets = _sets(feats, membership)
        point = (
            both_measures(xa, ages, sex, ta, feats, island, sets),
            both_measures(xb, ages, sex, tb, feats, island, sets),
        )
        boots = []
        for _ in range(n_boot):
            i = rng.integers(0, len(both), len(both))
            boots.append(
                (
                    both_measures(xa[i], ages[i], sex[i], ta[i], feats, island, sets),
                    both_measures(xb[i], ages[i], sex[i], tb[i], feats, island, sets),
                )
            )
        for r in paired_difference(point[0], point[1], boots):
            renamed = {k.replace("glia", "first").replace("neurons", "second"): v for k, v in r.items()}
            rows.append({"first": a, "second": b, **renamed})
        shared_chip = float((pa["chip"].to_numpy() == pb["chip"].to_numpy()).mean())
        sizes.append(
            {
                "first": a,
                "second": b,
                "donors": len(both),
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "share_of_donors_with_both_on_one_chip": shared_chip,
                "resamples": n_boot,
            }
        )
        print(f"cross_tissue: {a} and {b}: {len(both)} donors, {len(feats)} CpGs", file=sys.stderr, flush=True)
    if rows:
        write_table(pl.DataFrame(sizes), out / "sizes")
        write_table(pl.DataFrame(rows), out / "paired")
    return []
