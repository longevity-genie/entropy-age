"""Division against arrest in culture: does the shield depend on how cells became senescent? (roadmap RM40)

GSE197723 holds human fibroblasts on the EPIC array in a balanced layout: each of 12 chips carries two untreated
cultures, two driven to replicative senescence (many divisions), two arrested by X-rays and two by oncogenic RAS
(arrest with few divisions). Every comparison is therefore within chip. The series was found by the external
dataset search (``docs/DATASET_SEARCH_2026-11-18.md``) and is rebuilt here from its IDATs.

If the shield is built by selection among dividing cells, or by any division-coupled process, cultures that reached
senescence by dividing should show a larger shield (essential promoters deviating less than matched CpGs) than
cultures arrested without dividing, at equal distance from the untreated state.

Per chip, the reference is the mean of its two untreated cultures. Each treated culture is scored against it as in
``per_person``: ``d = |x - m| / sqrt(m (1 - m))`` per CpG; ``B`` the matched CpGs' mean deviation, ``P`` the set's
minus the matched CpGs', ``F = P / B``, ``K`` the same matched gap of relative entropy, ``S`` of the signed
deviation; strata are island status x decile of the reference level over all untreated cultures. Because ``F`` is a
share of total deviation, array noise dilutes it where real drift is small; ``B2``, ``P2`` and ``F2`` repeat the
scores on squared deviations with a noise floor subtracted, estimated per chip from its two untreated replicates.
Scores are averaged
over a chip's two replicates, giving one value per chip and arm. Arms are compared by the mean of the within-chip
differences, with an exact sign-flip test over the 12 chips (4,096 assignments) and a t interval.

Limits read from the record: donor or cell-line identity per chip is not deposited, nor are the doublings of the
senescent cultures; "chip" is the unit, and whether chips are different lines is unknown.

``entropy-age figures senescence_modes`` (needs the SeSAMe rebuild of GSE197723).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import LEVEL_CLIP, deviations, matched_scores, signed_deviations
from entropy_age.stochasticity import baseline_strata, relative_entropy

GSE = "GSE197723"
CONTROL = "None"
ARMS = ("Replication", "X-ray", "RAS")
CONTRASTS = (("Replication", "X-ray"), ("Replication", "RAS"), ("X-ray", "RAS"))
MEASURES = ("B", "P", "F", "K", "S")


def noise_floor(controls: np.ndarray, reference: np.ndarray) -> np.ndarray | None:
    """Per CpG, the squared scaled deviation a treated culture would show from array noise alone, estimated from a
    chip's two untreated replicates: ``(u1 - u2)^2`` has variance ``2 s^2`` and a treated culture against their mean
    ``1.5 s^2``, so the floor is ``0.75 (u1 - u2)^2 / (m (1 - m))``. None unless the chip has exactly two."""
    if controls.shape[0] != 2:
        return None
    m = np.clip(reference, LEVEL_CLIP, 1 - LEVEL_CLIP)
    return 0.75 * (controls[0] - controls[1]) ** 2 / (m * (1 - m))


def sample_scores(
    x: np.ndarray, reference: np.ndarray, strata: np.ndarray, sets: dict, floor: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Scores of each row of ``x`` against ``reference``: ``<measure>_<set>`` -> one value per row. With ``floor``
    (see :func:`noise_floor`) also ``B2``, ``P2``, ``F2``: the same on squared deviations with the noise floor
    subtracted, so that a share of drift is not diluted by array noise where drift is small."""
    d = deviations(x, reference)
    signed = signed_deviations(x, reference)
    kl = relative_entropy(np.clip(x, 0.0, 1.0), reference)
    out = {}
    for group, (in_set, _) in sets.items():
        p, b = matched_scores(d, in_set, strata)
        out[f"B_{group}"], out[f"P_{group}"], out[f"F_{group}"] = b, p, p / b
        out[f"K_{group}"] = matched_scores(kl, in_set, strata)[0]
        out[f"S_{group}"] = matched_scores(signed, in_set, strata)[0]
        if floor is not None:
            p2, b2 = matched_scores(d**2, in_set, strata)
            fp, fb = matched_scores(floor[None, :], in_set, strata)
            out[f"B2_{group}"], out[f"P2_{group}"] = b2 - fb[0], p2 - fp[0]
            out[f"F2_{group}"] = (p2 - fp[0]) / (b2 - fb[0])
    return out


def sign_flip(diff: np.ndarray) -> float:
    """Exact two-sided sign-flip p of the mean of paired differences (all 2^n assignments, n <= 16)."""
    diff = diff[np.isfinite(diff)]
    if diff.size < 4 or diff.size > 16:
        return float("nan")
    signs = np.array(list(itertools.product((1.0, -1.0), repeat=diff.size)))
    null = np.abs(signs @ diff) / diff.size
    return float((null >= abs(diff.mean()) - 1e-15).mean())


def paired(values: pl.DataFrame, a: str, b: str, column: str) -> dict:
    """Arm ``a`` minus arm ``b`` within chips: mean difference, t interval, exact sign-flip p."""
    wide = values.filter(pl.col("arm").is_in([a, b])).pivot(on="arm", index="chip", values=column).drop_nulls()
    diff = (wide[a] - wide[b]).to_numpy()
    n = diff.size
    se = diff.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    half = stats.t.ppf(0.975, n - 1) * se if n > 1 else np.nan
    return {
        "contrast": f"{a} - {b}",
        "chips": n,
        "mean_a": float(wide[a].mean()),
        "mean_b": float(wide[b].mean()),
        "difference": float(diff.mean()),
        "ci_low": float(diff.mean() - half),
        "ci_high": float(diff.mean() + half),
        "p_sign_flip": sign_flip(diff),
        "chips_negative": int((diff < 0).sum()),
    }


def run(ctx: Context) -> list[Path]:
    if not reprocessed.available(ctx.paths, GSE):
        print(f"senescence_modes: no SeSAMe rebuild for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("senescence_modes")
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, GSE).filter(pl.col("passed"))
    chars = _characteristics(ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz")
    samples = (
        qc.select("sample", "chip")
        .join(chars.select("sample", "cell type", pl.col("condition").alias("arm")), on="sample")
        .filter((pl.col("cell type") == "Fibroblast") & pl.col("arm").is_in([CONTROL, *ARMS]))
    )
    complete = (
        samples.group_by("chip")
        .agg(*((pl.col("arm") == a).sum().alias(a) for a in (CONTROL, *ARMS)))
        .filter(pl.all_horizontal(pl.col(a) >= 1 for a in (CONTROL, *ARMS)))
    )
    samples = samples.filter(pl.col("chip").is_in(complete["chip"].implode())).sort("chip", "arm", "sample")
    betas = reprocessed.by_accession(ctx.paths, GSE, {s: s for s in samples["sample"].to_list()})
    x, feats = complete_matrix(betas, features, samples["sample"].to_list())
    arm, chip = samples["arm"].to_numpy(), samples["chip"].to_numpy()
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    sets = _sets(feats, membership)
    overall = x[arm == CONTROL].mean(axis=0)
    strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(overall, 10).astype(str))
    rows = []
    for c in np.unique(chip):
        on = chip == c
        reference = x[on & (arm == CONTROL)].mean(axis=0)
        treated = on & (arm != CONTROL)
        floor = noise_floor(x[on & (arm == CONTROL)], reference)
        scores = sample_scores(x[treated], reference, strata, sets, floor)
        for i, (s, a) in enumerate(zip(samples["sample"].to_numpy()[treated], arm[treated], strict=True)):
            rows.append({"sample": s, "chip": c, "arm": a, **{k: float(v[i]) for k, v in scores.items()}})
    per_sample = pl.DataFrame(rows)
    columns = [c for c in per_sample.columns if c not in ("sample", "chip", "arm")]
    per_chip = per_sample.group_by("chip", "arm").agg(pl.col(columns).mean())
    effects = []
    for column in columns:
        measure, group = column.split("_", 1)
        for a, b in CONTRASTS:
            effects.append({"measure": measure, "group": group, **paired(per_chip, a, b, column)})
        for a in ARMS:  # each arm against zero: is there a shield at all in this arm?
            v = per_chip.filter(pl.col("arm") == a)[column].to_numpy()
            effects.append(
                {
                    "measure": measure,
                    "group": group,
                    "contrast": a,
                    "chips": v.size,
                    "mean_a": float(v.mean()),
                    "mean_b": 0.0,
                    "difference": float(v.mean()),
                    "ci_low": float(v.mean() - stats.t.ppf(0.975, v.size - 1) * v.std(ddof=1) / np.sqrt(v.size)),
                    "ci_high": float(v.mean() + stats.t.ppf(0.975, v.size - 1) * v.std(ddof=1) / np.sqrt(v.size)),
                    "p_sign_flip": sign_flip(v),
                    "chips_negative": int((v < 0).sum()),
                }
            )
    write_table(per_sample, out / "per_sample")
    write_table(per_chip.sort("chip", "arm"), out / "per_chip")
    write_table(pl.DataFrame(effects), out / "effects")
    write_table(
        pl.DataFrame({"chips": [int(np.unique(chip).size)], "cultures": [int(samples.height)], "cpgs": [len(feats)]}),
        out / "sizes",
    )
    print(
        f"senescence_modes: {np.unique(chip).size} chips, {samples.height} cultures, {len(feats)} CpGs",
        file=sys.stderr,
    )
    return []
