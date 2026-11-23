"""Maintenance or selection? Drift inside single-cell clones (protocol change 1, P5).

Flunkert et al. 2018 (GSE112812 strain 1, GSE112873 strain 2): human fetal fibroblasts; single cells were cloned and
each clone grown about 20 doublings without other clones beside it; the starting culture was arrayed twice. Both
series are rebuilt from IDATs with the shared recipe. Primary: untreated (control) clones; irradiated clones are
exploratory.

Per clone and promoter CpG, the relative entropy of the clone's level from the starting culture's mean level; per
CpG its mean over the clones of a strain. The housekeeping gap is the set's mean minus matched CpGs (island status x
starting-culture level decile), in SD over CpGs; the interval resamples clones within strain, strata and SD fixed
from the full sample. The same statistic for one starting-culture array against the other is the technical floor:
array noise alone, which a clone effect has to exceed.

Prediction (protocol): replication-coupled maintenance keeps the gap inside a lineage, so clones show a negative gap
beyond the technical floor; selection among lineages cannot act inside one clone, so clones show none.
``entropy-age figures fibroblast_clones``.
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
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, relative_entropy

SERIES = {"strain 1": "GSE112812", "strain 2": "GSE112873"}
GROUP = re.compile(r"fetal fibroblasts_(control clone|irradiated clone|starting culture)_")
N_BOOT = {"paper": 400, "quick": 30}


def read_samples(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM) and ``kind`` (control clone, irradiated clone, starting culture) from the titles."""
    fields: dict[str, list[str]] = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            key, _, rest = line.rstrip("\n").partition("\t")
            if key in ("!Sample_title", "!Sample_geo_accession"):
                fields[key] = [v.strip('"') for v in rest.split("\t")]
    rows = []
    for gsm, title in zip(fields["!Sample_geo_accession"], fields["!Sample_title"], strict=True):
        m = GROUP.search(title)
        if m:
            rows.append({"sample": gsm, "kind": m.group(1)})
    return pl.DataFrame(rows)


def gap(values: np.ndarray, sets: dict, strata: np.ndarray, sd: float) -> dict[str, float]:
    return {
        group: float(matched_scores(values[None, :], in_set, strata)[0][0]) / sd for group, (in_set, _) in sets.items()
    }


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    out = ctx.out_dir("fibroblast_clones")
    features, probes, membership = _features(paths, ctx.rng(160_000))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    rng = ctx.rng(161_000)
    rows, sizes = [], []
    for strain, gse in SERIES.items():
        matrix = paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
        if not (reprocessed.available(paths, gse) and matrix.exists()):
            print(f"fibroblast_clones: {gse} not rebuilt; skipped", file=sys.stderr)
            continue
        passed = reprocessed.read_qc(paths, gse).filter(pl.col("passed"))["sample"].to_list()
        samples = read_samples(matrix).filter(pl.col("sample").is_in(passed))
        betas = reprocessed.by_accession(paths, gse, {s: s for s in samples["sample"].to_list()})
        samples = samples.filter(pl.col("sample").is_in(betas.columns))
        parents = samples.filter(pl.col("kind") == "starting culture")["sample"].to_list()
        if len(parents) < 2:
            print(f"fibroblast_clones: {strain}: {len(parents)} starting-culture arrays; skipped", file=sys.stderr)
            continue
        x, feats = complete_matrix(betas, features, samples["sample"].to_list())
        index = {s: i for i, s in enumerate(samples["sample"].to_list())}
        parent = x[[index[p] for p in parents]].mean(axis=0)
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(parent, 10).astype(str))
        sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
        # technical floor: one starting-culture array against the other
        floor = relative_entropy(x[index[parents[0]]], x[index[parents[1]]])
        floor_sd = float(floor.std())
        for group, value in gap(floor, sets, strata, floor_sd).items():
            rows.append(
                {
                    "strain": strain,
                    "clones": "starting culture, array vs array",
                    "group": group,
                    "gap_sd": value,
                    "ci_low": float("nan"),
                    "ci_high": float("nan"),
                    "n": 2,
                }
            )
        for kind in ("control clone", "irradiated clone"):
            clones = [index[s] for s in samples.filter(pl.col("kind") == kind)["sample"].to_list()]
            if len(clones) < 5:
                continue
            per_clone = np.vstack([relative_entropy(x[i], parent) for i in clones])
            values = per_clone.mean(axis=0)
            sd = float(values.std())
            point = gap(values, sets, strata, sd)
            draws: dict[str, list[float]] = {g: [] for g in sets}
            for _ in range(N_BOOT[str(ctx.scale.name)]):
                pick = rng.integers(0, len(clones), len(clones))
                for g, v in gap(per_clone[pick].mean(axis=0), sets, strata, sd).items():
                    draws[g].append(v)
            for g, v in point.items():
                rows.append(
                    {
                        "strain": strain,
                        "clones": kind,
                        "group": g,
                        "gap_sd": v,
                        "ci_low": float(np.quantile(draws[g], 0.025)),
                        "ci_high": float(np.quantile(draws[g], 0.975)),
                        "n": len(clones),
                    }
                )
            sizes.append(
                {
                    "strain": strain,
                    "clones": kind,
                    "n": len(clones),
                    "cpgs": len(feats),
                    "mean_distance": float(values.mean()),
                    "floor_mean_distance": float(floor.mean()),
                }
            )
        print(f"fibroblast_clones: {strain}: {len(feats)} CpGs, {samples.height} arrays", file=sys.stderr, flush=True)
    if rows:
        write_table(pl.DataFrame(rows), out / "gaps")
        write_table(pl.DataFrame(sizes), out / "sizes")
    return []
