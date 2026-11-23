"""Exercise training and unloading in human skeletal muscle, the same people before and after (roadmap RM40).

Two public series of vastus lateralis biopsies on Illumina arrays, rebuilt from IDATs:

* GSE171140 (Gene SMART; Voisin et al.): young adults before (``PRE``) and after 4, 8 and 12 weeks of high-intensity
  interval training (``4WP``, ``8WP``, ``12WP``), and at the start of a control period without training (``CON``)
  taken weeks before ``PRE``, so ``CON`` against ``PRE`` is elapsed time without training.
* GSE114763 (Seaborne et al. 2018): eight men at baseline, after one bout, after 7 weeks of loading, 7 weeks of
  unloading and 7 weeks of reloading.

Per CpG, methylation is regressed on condition (the baseline as reference) with person and chip intercepts. For each
condition the squared coefficient minus its squared standard error estimates the squared change from baseline
without array noise (``B2`` over CpGs matched to the set, ``P2`` the set minus them, ``F2 = P2 / B2``), as in
``senescence_donors``. A negative ``F2`` means essential promoters changed less than look-alike CpGs. Strata are island
status x decile of the baseline level. Each condition is also refitted with each person left out in turn.

Limits: young adults only (no age contrast); training moves few CpGs, so ``B2`` may be near zero and ``F2`` then
means nothing (read ``B2`` first); conditions are not on balanced chips in every person.
``entropy-age figures muscle_training`` (needs the SeSAMe rebuilds).
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.senescence_donors import arm_effects, squared_share
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import LEVEL_CLIP
from entropy_age.stochasticity import baseline_strata

SERIES = {
    "GSE171140": {
        "pattern": r"Muscle_([^_]+)_(GS\d+)_",
        "reference": "PRE",
        "conditions": {
            "4WP": "4 weeks of training",
            "8WP": "8 weeks",
            "12WP": "12 weeks",
            "CON": "no training (control period)",
        },
    },
    "GSE114763": {
        "pattern": r"Muscle_(.+?)_P(?:articipant_)?(\d+)",
        "reference": "Baseline",
        "conditions": {
            "Acute": "one bout",
            "7wk_Loading": "7 weeks of loading",
            "7wk_Unloading": "7 weeks of unloading",
            "7wk_Re-oading": "7 weeks of reloading",
        },
    },
}
MIN_PEOPLE = 6


def read_titles(matrix: Path) -> tuple[list[str], list[str]]:
    rows = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith(("!Sample_title", "!Sample_geo_accession")):
                rows[line.split("\t", 1)[0]] = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
    return rows["!Sample_geo_accession"], rows["!Sample_title"]


def design(samples: list[str], titles: list[str], pattern: str) -> pl.DataFrame:
    """``sample``, ``condition`` and ``person`` read from the titles; titles that do not match raise."""
    rx = re.compile(pattern)
    rows = []
    for s, t in zip(samples, titles, strict=True):
        found = rx.search(t)
        if found is None:
            raise ValueError(f"muscle_training: title {t!r} does not match {pattern!r}")
        rows.append({"sample": s, "condition": found.group(1), "person": found.group(2).lstrip("0")})
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("muscle_training")
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    shares, sizes = [], []
    for gse, spec in SERIES.items():
        matrix = ctx.paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
        if not (reprocessed.available(ctx.paths, gse) and matrix.exists()):
            print(f"muscle_training: no rebuild of {gse}; skipped", file=sys.stderr)
            continue
        qc = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))
        table = design(*read_titles(matrix), spec["pattern"]).join(qc.select("sample", "chip"), on="sample")
        keep = [spec["reference"], *spec["conditions"]]
        table = (
            table.filter(pl.col("condition").is_in(keep)).unique(["person", "condition"], keep="first").sort("sample")
        )
        with_reference = set(table.filter(pl.col("condition") == spec["reference"])["person"].to_list())
        table = table.filter(pl.col("person").is_in(list(with_reference)))
        betas = reprocessed.by_accession(ctx.paths, gse, {s: s for s in table["sample"].to_list()})
        table = table.filter(pl.col("sample").is_in(betas.columns))
        x, feats = complete_matrix(betas, features, table["sample"].to_list())
        arm, person, chip = (table[c].to_numpy() for c in ("condition", "person", "chip"))
        base = np.clip(x[arm == spec["reference"]].mean(axis=0), LEVEL_CLIP, 1 - LEVEL_CLIP)
        scale = base * (1 - base)
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(base, 10).astype(str))
        sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
        arms = list(spec["conditions"])
        effects = arm_effects(x, arm, person, chip.astype(str), arms=arms, reference=spec["reference"])
        for a, (b, se) in effects.items():
            people = int(np.unique(person[arm == a]).size)
            if people < MIN_PEOPLE:
                continue
            for group, (in_set, _) in sets.items():
                shares.append({"series": gse, "condition": spec["conditions"][a], "people": people, "group": group,
                               "left_out": "none", **squared_share(b, se, scale, in_set, strata)})  # fmt: skip
        for p in np.unique(person):
            rest = person != p
            try:
                loo = arm_effects(
                    x[rest], arm[rest], person[rest], chip[rest].astype(str), arms=arms, reference=spec["reference"]
                )
            except ValueError:
                continue
            for a, (b, se) in loo.items():
                if np.unique(person[rest & (arm == a)]).size < MIN_PEOPLE:
                    continue
                for group, (in_set, _) in sets.items():
                    shares.append({"series": gse, "condition": spec["conditions"][a], "people": None, "group": group,
                                   "left_out": str(p), **squared_share(b, se, scale, in_set, strata)})  # fmt: skip
        sizes.append({"series": gse, "arrays": table.height, "people": int(np.unique(person).size), "cpgs": len(feats),
                      "chips": int(np.unique(chip).size)})  # fmt: skip
        print(
            f"muscle_training: {gse}: {table.height} arrays, {np.unique(person).size} people",
            file=sys.stderr,
            flush=True,
        )
    if shares:
        write_table(pl.DataFrame(shares), out / "squared_share")
        write_table(pl.DataFrame(sizes), out / "sizes")
    return []
