"""Senescence by division and by damage in the same donors' cells, with chip held (roadmap RM40, GSE227160).

GSE227160 (Kwiatkowska et al. 2023; found by the external dataset search) took skin fibroblasts of two donors and
mesenchymal stromal cells of four, and ran each donor's cells on the EPIC array as a young culture, a
middle-passage culture, two cultures in replicative senescence and two each made senescent early by irradiation and
by doxorubicin. Every donor has one array on each of the 8 chips, with the treatments spread over the chips, so
donor and chip can both be held. The series is rebuilt from IDATs.

If promoter drift and the shield come with division, cultures that divided to senescence should have drifted from
the young state and carry the gap, and cultures arrested early by damage should show little drift.

Per CpG one linear model, methylation on donor, chip and treatment (young as the reference). For each treatment the
coefficient ``b`` is that treatment's change from young with donor and chip held, and ``se`` its standard error.
Two readings, both scaled by the young level ``m`` as ``1 / (m (1 - m))``:

* ``abs_change``: ``|b| / sqrt(m (1 - m))`` per CpG, then set CpGs against CpGs matched on island status and young
  level (gap in SD, interval from resampling the set's genes);
* the share of squared drift, noise removed: ``b^2 - se^2`` estimates the squared change without the part that
  array noise adds, ``B2`` is its mean over matched CpGs, ``P2`` the set's mean minus the matched CpGs', and
  ``F2 = P2 / B2`` (−0.5 means the set carries half the squared drift of its look-alikes). ``F2`` is also given
  with each donor left out in turn.

Limits: six donors of two cell types pooled with a donor intercept; one young culture per donor; doublings are not
deposited; the model assumes chip acts alike on every donor.
``entropy-age figures senescence_donors`` (needs the SeSAMe rebuild of GSE227160).
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
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import LEVEL_CLIP, matched_scores
from entropy_age.stochasticity import baseline_strata

GSE = "GSE227160"
YOUNG = "Young"
ARMS = {
    "Middle-passaged": "middle passage",
    "Replicative sensescent": "replicative senescence",  # the depositors' spelling
    "Irradiation-induced prematurely senescent": "irradiation",
    "Doxorubicin-induced prematurely senescent": "doxorubicin",
}


def arm_effects(
    x: np.ndarray, arm: np.ndarray, donor: np.ndarray, chip: np.ndarray, arms=None, reference: str = YOUNG
) -> dict[str, tuple]:
    """Per treatment, (coefficient, standard error) per CpG of methylation on treatment (young as reference) with
    donor and chip intercepts held. Treatments absent from ``arm`` are left out."""
    arms = list(ARMS) if arms is None else list(arms)
    unknown = set(arm.tolist()) - set(arms) - {reference}
    if unknown:  # a label outside the known arms would silently join the young reference
        raise ValueError(f"senescence_donors: unknown treatment labels {sorted(unknown)}")
    present = [a for a in arms if (arm == a).any()]
    indicators = np.column_stack([(arm == a).astype(float) for a in present])
    nuisance = np.column_stack([np.ones(arm.size), dummies(donor), dummies(chip)])
    design = np.column_stack([indicators, nuisance])
    design = design[:, np.concatenate([np.ones(len(present), bool), np.ptp(nuisance, axis=0) > 1e-12])]
    design = np.column_stack([design, np.ones(arm.size)])
    rank = np.linalg.matrix_rank(design)
    if rank < design.shape[1]:  # a treatment fixed by donor and chip is unidentifiable, not to be adjusted
        raise ValueError(f"senescence_donors: design of {design.shape[1]} columns has rank {rank}")
    pinv = np.linalg.pinv(design)
    coef = pinv @ x
    resid = x - design @ coef
    sigma2 = (resid**2).sum(axis=0) / (arm.size - rank)
    unscaled = np.einsum("ij,ij->i", pinv, pinv)
    return {a: (coef[i], np.sqrt(sigma2 * unscaled[i])) for i, a in enumerate(present)}


def squared_share(b: np.ndarray, se: np.ndarray, scale: np.ndarray, in_set: np.ndarray, strata: np.ndarray) -> dict:
    """``B2`` (matched CpGs' mean squared change, noise removed), ``P2`` (set minus matched) and ``F2 = P2 / B2``."""
    p2, b2 = matched_scores(((b**2 - se**2) / scale)[None, :], in_set, strata)
    return {"B2": float(b2[0]), "P2": float(p2[0]), "F2": float(p2[0] / b2[0])}


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    matrix = paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(paths, GSE) and matrix.exists()):
        print(f"senescence_donors: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("senescence_donors")
    features, _, membership = _features(paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(paths, GSE).filter(pl.col("passed"))
    chars = _characteristics(matrix).select("sample", pl.col("treatment").alias("arm"), "donor")
    samples = chars.join(qc.select("sample", pl.col("chip").cast(pl.String)), on="sample").sort("donor", "arm")
    betas = reprocessed.by_accession(paths, GSE, {s: s for s in samples["sample"].to_list()})
    x, feats = complete_matrix(betas, features, samples["sample"].to_list())
    arm, donor, chip = (samples[c].to_numpy() for c in ("arm", "donor", "chip"))
    young = np.clip(x[arm == YOUNG].mean(axis=0), LEVEL_CLIP, 1 - LEVEL_CLIP)
    scale = young * (1 - young)
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(young, 10).astype(str))
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    gaps, shares = [], []
    for a, (b, se) in arm_effects(x, arm, donor, chip).items():
        gaps += [
            {"arm": ARMS[a], **r}
            for r in matched_gaps(
                np.abs(b) / np.sqrt(scale), sets, {"island x level": strata}, N_RESAMPLES[str(ctx.scale.name)], ctx.seed
            )
        ]
        for group, (in_set, _) in sets.items():
            shares.append(
                {"arm": ARMS[a], "group": group, "left_out": "none", **squared_share(b, se, scale, in_set, strata)}
            )
    for d in np.unique(donor):
        keep = donor != d
        for a, (b, se) in arm_effects(x[keep], arm[keep], donor[keep], chip[keep]).items():
            for group, (in_set, _) in sets.items():
                shares.append(
                    {"arm": ARMS[a], "group": group, "left_out": d, **squared_share(b, se, scale, in_set, strata)}
                )
    layout = samples.group_by("chip", "arm").len().pivot(on="arm", index="chip", values="len").fill_null(0)
    write_table(pl.DataFrame(gaps), out / "matched_gaps")
    write_table(pl.DataFrame(shares), out / "squared_share")
    write_table(layout.sort("chip"), out / "chip_by_arm")
    write_table(
        pl.DataFrame(
            {
                "cultures": [samples.height],
                "donors": [int(np.unique(donor).size)],
                "chips": [int(np.unique(chip).size)],
                "cpgs": [len(feats)],
            }
        ),
        out / "sizes",
    )
    print(f"senescence_donors: {samples.height} cultures, {len(feats)} CpGs", file=sys.stderr, flush=True)
    return []
