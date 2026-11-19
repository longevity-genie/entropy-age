"""A second look at the breast-cancer lead: drift at unprotected promoters in GSE283951 (roadmap RM40, claim C12).

In EPIC-Italy, drift at matched, unprotected promoter CpGs (``B``) was associated with breast cancer diagnosed years
later, but the adjusted association did not survive the rebuild from raw IDATs (R96). GSE283951, found by the
external dataset search, is an independent series on the EPIC array with raw IDATs: blood of 34 Polish women with
breast cancer and 56 healthy women, cases and controls sharing chips. It is rebuilt here with the project's recipe.

Per woman: ``B`` (matched CpGs' deviation), ``F`` (shield share) and ``S`` (signed shift) of the housekeeping and
core-essential sets, scored against the controls (each control against the other controls), as in
``consortium_interventions.arm_scores``. Each score is centred within its chip and standardised, then entered in a
logistic regression of case status with age: odds ratio per SD. A second model adds each woman's mean methylation
over the promoter CpGs analysed, also chip-centred. Cell fractions are not held (no reference fit in this module).

Limits from the record: 34 cases; the time from blood draw to diagnosis is not deposited per sample, so this is a
case-control contrast and cannot say the signal precedes disease; 21 of the controls are arrays shared with a sister
series whose pre-diagnostic status is unverified.

``entropy-age figures breast_replication`` (needs the SeSAMe rebuild of GSE283951).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.incident_cancer import _group_means, _odds_row
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE283951"
CASE = "Polish women with BC"
GROUPS = ("housekeeping", "core_essential")
MEASURES = ("B", "F", "S")


def standardised_within_chip(v: np.ndarray, chips: np.ndarray) -> np.ndarray:
    """``v`` minus its chip mean, in SDs of the centred values."""
    centred = v - _group_means(v, chips)
    return centred / centred.std()


def run(ctx: Context) -> list[Path]:
    if not reprocessed.available(ctx.paths, GSE):
        print(f"breast_replication: no SeSAMe rebuild for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("breast_replication")
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, GSE).filter(pl.col("passed"))
    chars = _characteristics(ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz")
    people = (
        qc.select("sample", "chip")
        .join(chars.select("sample", "subject status", pl.col("age").cast(pl.Float64, strict=False)), on="sample")
        .drop_nulls("age")
        .sort("sample")
    )
    betas = reprocessed.by_accession(ctx.paths, GSE, {s: s for s in people["sample"].to_list()})
    x, feats = complete_matrix(betas, features, people["sample"].to_list())
    case = (people["subject status"] == CASE).to_numpy()
    chips, age = people["chip"].to_numpy(), people["age"].to_numpy()
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    sets = {g: s for g, s in _sets(feats, membership).items() if g in GROUPS}
    scores = arm_scores(x, case, island, sets)
    y = case.astype(float)
    age_c = (age - age.mean()) / age.std()
    mean_meth = standardised_within_chip(x.mean(axis=1), chips)
    rows = []
    for group in GROUPS:
        for measure in MEASURES:
            z = standardised_within_chip(scores[f"{measure}_{group}"], chips)
            for adjusted, extra in (("age", []), ("age + mean methylation", [mean_meth])):
                design = np.column_stack([np.ones_like(z), z, age_c, *extra])
                rows.append({**_odds_row("breast (case-control)", False, adjusted, f"{measure} {group}", y, design)})
    table = pl.DataFrame(rows).drop("lag_over_2y")
    write_table(table, out / "odds")
    both = sum(1 for c in np.unique(chips) if case[chips == c].any() and (~case[chips == c]).any())
    write_table(
        pl.DataFrame(
            {
                "cases": [int(case.sum())],
                "controls": [int((~case).sum())],
                "chips": [int(np.unique(chips).size)],
                "chips_with_both": [both],
                "cpgs": [len(feats)],
                "age_cases": [float(age[case].mean())],
                "age_controls": [float(age[~case].mean())],
            }
        ),
        out / "sizes",
    )
    write_table(
        people.with_columns(**{k: pl.Series(v) for k, v in scores.items() if k.split("_", 1)[0] in MEASURES}),
        out / "per_person",
    )
    print(
        f"breast_replication: {int(case.sum())} cases, {int((~case).sum())} controls, {len(feats)} CpGs",
        file=sys.stderr,
    )
    return []
