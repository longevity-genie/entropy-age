"""Housekeeping promoters, other promoters and non-promoter CpGs: which class's drift tracks a clock?

Codex's fifth review of ``clock_drift_controls`` ended with the analysis it would run next: a three-way matched
comparison that tests the promoter-class reading directly and asks whether the small housekeeping margin survives
without the ratio F's shared denominator. Built here, on the donors and primary groups of ``clock_drift_controls``
(its ``group_registry``: arrays passing the clock track's QC, read from the same physical array in both rebuilds,
one primary group per accession fixed before any prediction is read). Registered as the figure
``clock_drift_classes``.

Per primary group:

1. **Universe.** The drift CpGs of ``clock_drift`` (autosomal 450k promoter CpGs: the gene-set members and a seeded
   promoter background) and a seeded sample of ``NON_PROMOTER`` eligible non-promoter CpGs (autosomal 450k CpGs
   without the promoter flag and without a promoter link in the 450k annotation), each measured in every array of the
   group. Every prediction input of every clock is left out (``clock_drift_controls.excluded_cpgs``: the mitotic
   clocks' CpGs, the local clock files, every data file a scored gallery clock reads, DunedinPACE's normalisation
   probes included), and the CpGs the series description drops. EPIC-only probes carry no 450k annotation and stay
   out.
2. **Three target classes**, all matched on island status x decile of the group's mean beta x decile of the SD of
   beta across its arrays x probe chemistry (``clock_drift_controls.match_strata``) and balanced to one CpG per gene
   or cluster: *housekeeping promoters* (one CpG of each housekeeping gene, at random), *matched non-housekeeping
   promoters* (promoter CpGs of genes in none of the gene sets, filled gene by gene to the housekeeping draw's
   stratum counts, one CpG per gene) and *matched non-promoter CpGs* (filled cluster by cluster, one CpG per
   cluster). Both promoter classes keep one rule for CpGs annotated to several genes (review 6, finding 3): a CpG
   enters only while every gene its promoter links name holds no CpG of the class yet
   (``clock_drift_controls.promoter_links``); ``sizes`` reports the largest count per linked gene. A cluster is a run
   of eligible non-promoter CpGs on one chromosome with gaps under ``CLUSTER_GAP`` bp, built on the whole eligible
   annotation before anything is left out or sampled (``nonpromoter_universe``; review 6, finding 4), so promoter
   CpGs never bridge two clusters and the labels hold through filtering and matching. The three draws are then
   thinned to their smallest count in every stratum (the common support), so each class holds the same number of
   CpGs, genes or clusters, in the same strata. ``K`` draws per group (``N_DRAWS``).
3. **Backgrounds.** Each draw has its own background: every universe CpG outside every gene set's members and
   outside that draw's three classes. With each person left out in turn, the reference level is the other donors'
   mean (as ``loo_scores``) and ``d = |x - m| / sqrt(m (1 - m))``. Each class's score is its stratum-weighted mean of
   ``d`` and the background's ``B`` the same with identical stratum weights (``shared_weights``: within a scoring
   stratum, island x decile of the reference level, every class and the background get the pooled class count's
   share), so the classes differ only in their CpGs.
4. **Model.** For each clock and draw, the partial correlation of the clock's value with a class's score given age,
   the quality terms of ``clock_drift_controls`` (chip, log mean intensity, detection), batch or centre where
   recorded, the disjoint panel's level and compression, and that draw's own background score B (``draw_rs``; age in
   the design makes this the acceleration's association, and the clock enters with its own scale). The reported r
   per class (``r_hk``, ``r_op`` other promoters, ``r_np`` non-promoters) is the mean over draws of these within-draw
   r, so every estimate rests on a background disjoint from its targets (review 6, finding 2). A check: the scores
   averaged over draws against one pooled background disjoint from the union of every draw's targets
   (``class_draws``), its stratum support verified (``support_gap``, ``balance_no_background``), give ``_pooled_bg``.
5. **Contrasts.** (a) ``r_op - r_np``: promoters beyond non-promoters; (b) ``r_hk - r_op``: housekeeping beyond other
   promoters. For a clock whose acceleration tracks lower deviation (epiTOC, r below 0) a negative contrast means the
   first class tracks it more strongly; for HypoClock (r above 0) a positive one does. Uncertainty: ``N_BOOT``
   resamples of donors within the group (each class's draw-mean r recomputed on the resampled rows; percentile 95 %
   interval and bootstrap SE). These intervals are **conditional on the constructed scores** (``interval_scope``):
   the strata, references, draws and scores stay as built on the whole cohort. The draws' spread is reported beside
   them (median, range and sign change of the per-draw contrasts). **Full procedure** (review 6, finding 5): for
   ``FULL_CLOCKS`` (epiTOC2 and HypoClock) in every non-brain primary group, ``FULL_BOOT`` donor resamples rebuild the
   matching strata, ``FULL_DRAWS`` class draws, the donor-aware leave-one-out references (every copy of a donor left
   out of its own reference) and the scores within each resample (``donor_resample``), written as
   ``donor_resampling`` and as the ``_full_`` columns.
6. **Study level.** One primary group per accession, so per clock over the non-brain primary groups: the median
   contrast, a DerSimonian-Laird random-effects mean with its normal 95 % interval and the Hartung-Knapp small-study
   interval (t on k - 1 df; ``hartung_knapp``), both with the conditional bootstrap SEs as within-study errors and,
   for the full-procedure clocks, again with the donor-resampling SEs; a percentile interval of the median from
   resampling accessions; and the count of accessions whose interval lies above or below zero.

Writes ``clock_drift_classes/associations`` (per group and clock), ``draws`` (per draw), ``sizes`` (per group: class
sizes, background sizes and support, per-gene caps, balance), ``donor_resampling`` (per group and full-procedure
clock), ``summary`` (per clock) and one figure (``clock_drift_classes``).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments import mitotic_clocks
from entropy_age.experiments.clock_disagreement import quality_terms
from entropy_age.experiments.clock_drift import GROUP, MIN_ARRAYS, residual
from entropy_age.experiments.clock_drift_controls import (
    MIN_RESID_DF,
    RESID_TOL,
    _held,
    _panel,
    adjust_design,
    corr,
    excluded_cpgs,
    flag_primary,
    gene_balanced,
    group_clock_scores,
    group_registry,
    match_strata,
    matched_placebo,
    per_gene_counts,
    placebo_pool,
    promoter_links,
    scaled_basis,
    series_scores,
)
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import deviations
from entropy_age.stochasticity import baseline_strata

CLASSES = ("hk", "op", "np")  # housekeeping promoters, other promoters, non-promoter CpGs
LABELS = {"hk": "housekeeping promoters", "op": "other promoters", "np": "non-promoter CpGs"}
MITOTIC = ("EpiTOC1", "epiTOC2", "epiTOC3", "HypoClock", "epiCMIT", "RepliTali")
CLOCKS = (*MITOTIC, "Hannum", "PhenoAge", "Horvathv1", "DunedinPACE")
N_DRAWS = {"paper": 20, "quick": 3}
N_BOOT = {"paper": 1000, "quick": 100}
NON_PROMOTER = 60_000
CLUSTER_GAP = 1000  # bp: non-promoter CpGs closer than this on one chromosome form one cluster
PER_UNIT = 1  # CpGs per gene or cluster
FULL_CLOCKS = ("epiTOC2", "HypoClock")  # the donor resampling of the full procedure (review 6, finding 5)
FULL_BOOT = {"paper": 100, "quick": 3}
FULL_DRAWS = {"paper": 4, "quick": 2}
INTERVAL_SCOPE = "conditional on the constructed scores (donors resampled, scores and strata fixed)"


# ----------------------------------------------------------------------------- pure helpers
def clusters(chrom: np.ndarray, pos: np.ndarray, gap: int = CLUSTER_GAP) -> np.ndarray:
    """A cluster label per CpG: CpGs on one chromosome whose neighbour (in position order) lies under ``gap`` bp away
    share a label."""
    chrom = np.asarray(chrom).astype(str)
    pos = np.asarray(pos, dtype=np.int64)
    order = np.lexsort((pos, chrom))
    c, p = chrom[order], pos[order]
    new = np.ones(order.size, dtype=bool)
    new[1:] = (c[1:] != c[:-1]) | (p[1:] - p[:-1] >= gap)
    out = np.empty(order.size, dtype=np.int64)
    out[order] = np.cumsum(new) - 1
    return out


def units(labels: np.ndarray, members: np.ndarray) -> list[np.ndarray]:
    """Index arrays of ``members`` (indices) grouped by ``labels`` of those members, in a fixed order."""
    members = np.asarray(members, dtype=int)
    if members.size == 0:
        return []
    lab = np.asarray(labels)[members]
    order = np.argsort(lab, kind="stable")
    lab, mem = lab[order], members[order]
    cuts = np.flatnonzero(lab[1:] != lab[:-1]) + 1
    return [np.sort(g) for g in np.split(mem, cuts)]


def thin(strata: np.ndarray, mask: np.ndarray, want: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """``mask`` reduced at random to ``want[k]`` CpGs in stratum code ``k`` (``strata`` integer codes)."""
    out = np.zeros(mask.size, dtype=bool)
    for k in np.flatnonzero(want):
        cand = np.flatnonzero(mask & (strata == k))
        out[rng.choice(cand, int(want[k]), replace=False)] = True
    return out


def class_draw(
    strata: np.ndarray,
    hk_genes: list[np.ndarray],
    op_genes: list[np.ndarray],
    np_units: list[np.ndarray],
    rng: np.random.Generator,
    per_unit: int = PER_UNIT,
    linked: list[np.ndarray] | None = None,
) -> dict[str, np.ndarray]:
    """One draw of the three classes (masks ``hk``, ``op``, ``np``): ``per_unit`` CpGs per housekeeping gene, the
    other promoters and the non-promoter CpGs filled to that draw's stratum counts with ``per_unit`` per gene or
    cluster, then all three thinned to their smallest count in every stratum (identical counts, the common
    support). With ``linked`` (per CpG, the genes its promoter links name; ``clock_drift_controls.promoter_links``)
    both promoter classes keep the cap over every linked gene (review 6, finding 3)."""
    codes = np.unique(strata, return_inverse=True)[1]
    k = int(codes.max()) + 1 if codes.size else 0
    hk = gene_balanced(hk_genes, codes.size, rng, per_unit, linked)
    op, _, _ = matched_placebo(codes, hk, op_genes, rng, per_unit, linked)
    np_, _, _ = matched_placebo(codes, hk, np_units, rng, per_unit)
    counts = [np.bincount(codes[m], minlength=k) for m in (hk, op, np_)]
    want = np.min(counts, axis=0)
    return {c: thin(codes, m, want, rng) for c, m in zip(CLASSES, (hk, op, np_), strict=True)}


def shared_weights(
    sets: list[np.ndarray], background: np.ndarray, codes: np.ndarray
) -> tuple[list[np.ndarray], np.ndarray, dict[str, float]]:
    """Positive weights scoring disjoint ``sets`` and a ``background`` with one set of stratum weights: within stratum
    ``k`` (non-negative integer ``codes``) each set's CpGs get ``w_k / n_set_k`` and the background's ``w_k / n_bg_k``,
    ``w_k`` the sets' pooled count in ``k`` over the strata holding every set and the background (summing to 1).
    Returns the per-set weights, the background's and ``balance``: the largest total-variation distance between two
    sets' stratum counts (``tvd``) and the share of set CpGs in strata left out (``unusable``)."""
    sets = [np.asarray(s, dtype=bool) for s in sets]
    background = np.asarray(background, dtype=bool)
    union = np.zeros(background.size, dtype=bool)
    for s in sets:
        if (union & s).any():
            raise ValueError("shared_weights: the sets must be disjoint")
        union |= s
    if (union & background).any():
        raise ValueError("shared_weights: the background must be disjoint from the sets")
    k = int(codes.max()) + 1 if codes.size else 0
    n = [np.bincount(codes[s], minlength=k).astype(float) for s in sets]
    n_bg = np.bincount(codes[background], minlength=k).astype(float)
    usable = np.all([c > 0 for c in n], axis=0) & (n_bg > 0)
    if not usable.any():
        raise ValueError("shared_weights: no stratum holds every set and the background")
    w = np.where(usable, np.sum(n, axis=0), 0.0)
    w = w / w.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        per = [np.where(usable, w / c, 0.0) for c in n]
        per_bg = np.where(usable, w / n_bg, 0.0)
    weights = [np.where(s, p[codes], 0.0) for s, p in zip(sets, per, strict=True)]
    bw = np.where(background, per_bg[codes], 0.0)
    shares = [c / max(c.sum(), 1) for c in n]
    tvd = max((0.5 * np.abs(a - b).sum() for i, a in enumerate(shares) for b in shares[i + 1 :]), default=0.0)
    total = sum(c.sum() for c in n)
    unusable = 1 - sum(c[usable].sum() for c in n) / total if total else float("nan")
    return weights, bw, {"tvd": float(tvd), "unusable": float(unusable)}


def loo_codes(level: np.ndarray, island: np.ndarray) -> np.ndarray:
    """Scoring strata of a reference level: island x decile of the level (``fit_reference``'s partition), as codes."""
    return np.asarray(island, dtype=int) * 11 + baseline_strata(level, 10)


def class_draws(
    strata: np.ndarray,
    hk_genes: list[np.ndarray],
    op_genes: list[np.ndarray],
    np_units: list[np.ndarray],
    member: np.ndarray,
    rng: np.random.Generator,
    n_draws: int,
    linked: list[np.ndarray] | None = None,
) -> tuple[list[tuple[list[np.ndarray], np.ndarray]], np.ndarray]:
    """``n_draws`` class draws (``class_draw``), each with its own background (every CpG outside the gene sets'
    members and outside that draw's three classes), and the pooled background: every CpG outside the members and
    outside the union of all draws' classes, disjoint from every target of every draw (review 6, finding 2). Returns
    ``(draws, pooled)``, ``draws`` as ``(sets, background)``."""
    member = np.asarray(member, dtype=bool)
    draws, used = [], member.copy()
    for _ in range(n_draws):
        m = class_draw(strata, hk_genes, op_genes, np_units, rng, linked=linked)
        sets = [m[c] for c in CLASSES]
        targets = np.any(sets, axis=0)
        used |= targets
        draws.append((sets, ~(member | targets)))
    return draws, ~used


def support_gap(strata: np.ndarray, draws: list[tuple[list[np.ndarray], np.ndarray]], background: np.ndarray) -> int:
    """Matching strata that hold a target CpG of some draw and no CpG of ``background`` (0: shared support)."""
    codes = np.unique(strata, return_inverse=True)[1]
    k = int(codes.max()) + 1 if codes.size else 0
    targets = np.any([np.any(sets, axis=0) for sets, _ in draws], axis=0)
    held = np.bincount(codes[targets], minlength=k) > 0
    return int((held & (np.bincount(codes[np.asarray(background, dtype=bool)], minlength=k) == 0)).sum())


N_CODES = 2 * 11  # loo_codes: island x (10 deciles + the non-finite stratum)


def class_scores(
    x: np.ndarray,
    island: np.ndarray,
    draws: list[tuple[list[np.ndarray], np.ndarray]],
    counts: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per draw and person, each set's score and the background's, against the leave-one-out reference (the other
    people's mean), with ``shared_weights``' stratum weights. ``draws`` holds ``(sets, background)``. ``counts`` (a
    donor resample's multiplicity per row of ``x``) makes the leave-out donor-aware: a donor's reference is the
    count-weighted mean of the other donors, every copy of the donor left out; donors with count 0 are not scored
    (NaN). Returns ``A`` (draws x people x sets), ``B`` (draws x people) and the per-draw mean balance (draws x 3: tvd,
    unusable, ``no_background``: the share of set CpGs in strata without a background CpG)."""
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    c = np.ones(n) if counts is None else np.asarray(counts, dtype=float)
    n_sets = len(draws[0][0]) if draws else 0
    a = np.full((len(draws), n, n_sets), np.nan)
    b = np.full((len(draws), n), np.nan)
    bal = np.zeros((len(draws), 3))
    keys: dict[bytes, int] = {}
    bg_of, bgs, set_idx = [], [], []
    for sets, background in draws:
        sets = [np.asarray(m, dtype=bool) for m in sets]
        background = np.asarray(background, dtype=bool)
        union = np.zeros(background.size, dtype=bool)
        for m in sets:
            if (union & m).any():
                raise ValueError("class_scores: the sets must be disjoint")
            union |= m
        if (union & background).any():
            raise ValueError("class_scores: the background must be disjoint from the sets")
        key = np.packbits(background).tobytes()
        if key not in keys:
            keys[key] = len(bgs)
            bgs.append(np.flatnonzero(background))
        bg_of.append(keys[key])
        set_idx.append([np.flatnonzero(m) for m in sets])
    total = c @ x
    size = c.sum()
    scored = np.flatnonzero(c > 0)
    for i in scored:
        level = (total - c[i] * x[i]) / (size - c[i])
        d = deviations(x[i], level)
        codes = loo_codes(level, island)
        bg_stats = [
            (np.bincount(codes[ix], minlength=N_CODES), np.bincount(codes[ix], weights=d[ix], minlength=N_CODES))
            for ix in bgs
        ]
        for j, sets in enumerate(set_idx):
            nb, sb = bg_stats[bg_of[j]]
            ns = np.array([np.bincount(codes[ix], minlength=N_CODES) for ix in sets], dtype=float)
            ss = np.array([np.bincount(codes[ix], weights=d[ix], minlength=N_CODES) for ix in sets])
            usable = (ns > 0).all(axis=0) & (nb > 0)
            if not usable.any():
                raise ValueError("class_scores: no stratum holds every set and the background")
            w = np.where(usable, ns.sum(axis=0), 0.0)
            w = w / w.sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                a[j, i] = (w * np.where(usable, ss / ns, 0.0)).sum(axis=1)
                b[j, i] = (w * np.where(usable, sb / nb, 0.0)).sum()
            shares = ns / np.maximum(ns.sum(axis=1, keepdims=True), 1)
            tvd = max(
                (0.5 * np.abs(shares[p] - shares[q]).sum() for p in range(len(sets)) for q in range(p + 1, len(sets))),
                default=0.0,
            )
            tot = ns.sum()
            unusable = 1 - ns[:, usable].sum() / tot if tot else np.nan
            no_bg = ns[:, nb == 0].sum() / tot if tot else np.nan
            bal[j] += np.array([tvd, unusable, no_bg]) / scored.size
    return a, b, bal


def partial_rs(
    v: np.ndarray, targets: np.ndarray, given: np.ndarray, min_df: int = MIN_RESID_DF
) -> tuple[np.ndarray, str]:
    """Partial correlations of ``v`` with each column of ``targets`` given an intercept and ``given``, with one
    column-scaled design for rank and fit (``scaled_basis``) and the guards of ``partial_r_checked``: NaN for every
    column when ``v`` is explained or the residual df fall under ``min_df``, NaN for a column the design explains."""
    v = np.asarray(v, dtype=float)
    t = np.asarray(targets, dtype=float).reshape(v.size, -1)
    given = np.asarray(given, dtype=float).reshape(v.size, -1)
    nan = np.full(t.shape[1], np.nan)
    if not (np.isfinite(v).all() and np.isfinite(t).all() and np.isfinite(given).all()):
        return nan, "non-finite value"
    basis, rank = scaled_basis(np.column_stack([np.ones(v.size), given]))
    if v.size - rank < min_df:
        return nan, f"residual df {v.size - rank} < {min_df}"
    rv = v - basis @ (basis.T @ v)
    size = np.linalg.norm(v - v.mean())
    if size == 0 or np.linalg.norm(rv) <= RESID_TOL * size:
        return nan, "clock explained by the design"
    rt = t - basis @ (basis.T @ t)
    sizes = np.linalg.norm(t - t.mean(axis=0), axis=0)
    out = np.array([corr(rv, rt[:, j]) for j in range(t.shape[1])])
    out[np.linalg.norm(rt, axis=0) <= RESID_TOL * np.where(sizes > 0, sizes, np.inf)] = np.nan
    return out, ""


def contrasts(r: np.ndarray) -> np.ndarray:
    """``(r_op - r_np, r_hk - r_op)`` from ``r`` ordered as ``CLASSES`` (last axis)."""
    r = np.asarray(r, dtype=float)
    return np.stack([r[..., 1] - r[..., 2], r[..., 0] - r[..., 1]], axis=-1)


def bootstrap_rs(
    v: np.ndarray, targets: np.ndarray, given: np.ndarray, rng: np.random.Generator, n_boot: int
) -> np.ndarray:
    """``partial_rs`` over ``n_boot`` resamples of the rows (donors) with replacement: n_boot x targets."""
    n = v.size
    out = np.empty((n_boot, np.asarray(targets).reshape(n, -1).shape[1]))
    for k in range(n_boot):
        idx = rng.integers(0, n, n)
        out[k] = partial_rs(v[idx], targets[idx], given[idx])[0]
    return out


def random_effects(est: np.ndarray, se: np.ndarray) -> tuple[float, float, float]:
    """DerSimonian-Laird pooled mean, its SE and tau^2 over studies with finite estimates and positive SEs."""
    est, se = np.asarray(est, dtype=float), np.asarray(se, dtype=float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    if est.size == 0:
        return float("nan"), float("nan"), float("nan")
    w = 1 / se**2
    fixed = (w * est).sum() / w.sum()
    q = (w * (est - fixed) ** 2).sum()
    c = w.sum() - (w**2).sum() / w.sum()
    tau2 = max(0.0, (q - (est.size - 1)) / c) if c > 0 else 0.0
    w = 1 / (se**2 + tau2)
    return float((w * est).sum() / w.sum()), float(np.sqrt(1 / w.sum())), float(tau2)


def draw_rs(v: np.ndarray, a: np.ndarray, b: np.ndarray, fixed: np.ndarray, min_df: int = MIN_RESID_DF) -> np.ndarray:
    """Per draw ``j``, ``partial_rs(v, a[j], [fixed, b[j]])`` (draws x targets): each draw's classes against its own
    background score, so every draw keeps its disjoint background (review 6, finding 2). One scaled basis of the
    intercept and ``fixed`` is extended by each draw's background column."""
    v = np.asarray(v, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    fixed = np.asarray(fixed, dtype=float).reshape(v.size, -1)
    out = np.full(a.shape[::2], np.nan)
    if not (np.isfinite(v).all() and np.isfinite(fixed).all()):
        return out
    q0, r0 = scaled_basis(np.column_stack([np.ones(v.size), fixed]))
    size = np.linalg.norm(v - v.mean())
    tol = max(v.size, r0 + 1) * np.finfo(float).eps * 10
    for j in range(a.shape[0]):
        t, bj = a[j], b[j]
        if not (np.isfinite(t).all() and np.isfinite(bj).all()):
            continue
        rb = bj - q0 @ (q0.T @ bj)
        nb = np.linalg.norm(bj)
        q = np.column_stack([q0, rb / np.linalg.norm(rb)]) if nb > 0 and np.linalg.norm(rb) > tol * nb else q0
        if v.size - q.shape[1] < min_df:
            continue
        rv = v - q @ (q.T @ v)
        if size == 0 or np.linalg.norm(rv) <= RESID_TOL * size:
            continue
        rt = t - q @ (q.T @ t)
        sizes = np.linalg.norm(t - t.mean(axis=0), axis=0)
        r = np.array([corr(rv, rt[:, k]) for k in range(t.shape[1])])
        r[np.linalg.norm(rt, axis=0) <= RESID_TOL * np.where(sizes > 0, sizes, np.inf)] = np.nan
        out[j] = r
    return out


def bootstrap_draw_rs(
    v: np.ndarray, a: np.ndarray, b: np.ndarray, fixed: np.ndarray, rng: np.random.Generator, n_boot: int
) -> np.ndarray:
    """The draw-mean of ``draw_rs`` over ``n_boot`` resamples of the scored rows (donors) with replacement: n_boot x
    targets. The scores stay as built on the full cohort, so the spread is conditional on them."""
    n = v.size
    out = np.empty((n_boot, a.shape[2]))
    for k in range(n_boot):
        idx = rng.integers(0, n, n)
        with np.errstate(invalid="ignore"):
            out[k] = np.nanmean(draw_rs(v[idx], a[:, idx], b[:, idx], fixed[idx]), axis=0)
    return out


def hartung_knapp(est: np.ndarray, se: np.ndarray, level: float = 0.95) -> tuple[float, float, float]:
    """The DerSimonian-Laird mean with the Hartung-Knapp small-study interval: SE from the weighted residual spread
    around the mean, ``sum w (y - mu)^2 / ((k - 1) sum w)`` with ``w = 1 / (se^2 + tau^2)``, and a t quantile on
    ``k - 1`` df. Returns ``(mu, lo, hi)``; NaN under two studies."""
    from scipy.stats import t as student

    est, se = np.asarray(est, dtype=float), np.asarray(se, dtype=float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    if est.size < 2:
        return float("nan"), float("nan"), float("nan")
    mu, _, tau2 = random_effects(est, se)
    w = 1 / (se**2 + tau2)
    se_hk = float(np.sqrt((w * (est - mu) ** 2).sum() / ((est.size - 1) * w.sum())))
    q = float(student.ppf(0.5 + level / 2, est.size - 1))
    return mu, mu - q * se_hk, mu + q * se_hk


def donor_resample(
    x: np.ndarray,
    island: np.ndarray,
    chemistry: np.ndarray,
    units_: tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]],
    member: np.ndarray,
    linked: list[np.ndarray] | None,
    clocks: dict[str, np.ndarray],
    fixed: np.ndarray,
    rng: np.random.Generator,
    n_boot: int,
    n_draws: int,
) -> dict[str, np.ndarray]:
    """The full procedure within donor resamples (review 6, finding 5): per resample of the donors with replacement,
    the matching strata are rebuilt from the resampled cohort, ``n_draws`` class draws are made on them, every donor
    is scored against the donor-aware leave-one-out reference of the resample (``class_scores`` with ``counts``),
    and each clock's draw-mean partial r (``draw_rs``) is fitted on the resampled rows. Returns per clock an
    ``n_boot`` x classes array."""
    hk_genes, op_genes, np_units = units_
    n = x.shape[0]
    out = {k: np.full((n_boot, len(CLASSES)), np.nan) for k in clocks}
    for k in range(n_boot):
        idx = rng.integers(0, n, n)
        counts = np.bincount(idx, minlength=n)
        strata = match_strata(x[idx], island, chemistry)
        draws, _ = class_draws(strata, hk_genes, op_genes, np_units, member, rng, n_draws, linked)
        a, b, _ = class_scores(x, island, draws, counts)
        for clock, v in clocks.items():
            use = np.isfinite(v[idx])
            rows = idx[use]
            with np.errstate(invalid="ignore"):
                out[clock][k] = np.nanmean(draw_rs(v[rows], a[:, rows], b[:, rows], fixed[rows]), axis=0)
    return out


# ----------------------------------------------------------------------------- run
def nonpromoter_universe(
    probes: pl.DataFrame,
    linked: set[str],
    exclude: set[str],
    rng: np.random.Generator,
    size: int = NON_PROMOTER,
    gap: int = CLUSTER_GAP,
) -> tuple[list[str], dict[str, int]]:
    """The eligible non-promoter annotation (autosomal 450k CpGs without the promoter flag and without a promoter link
    in ``linked``) clustered on its own (``clusters``) before anything is left out or sampled, so a promoter CpG never
    bridges two non-promoter CpGs and leaving a CpG out never splits a cluster (review 6, finding 4); then a seeded
    sample of ``size`` of its CpGs outside ``exclude``. Returns the sample and the cluster label of every eligible
    CpG."""
    elig = probes.filter(
        ~pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]) & ~pl.col("cpg").is_in(list(linked))
    ).sort("cpg")
    cluster_of = dict(zip(elig["cpg"].to_list(), clusters(elig["chr"], elig["pos"], gap).tolist(), strict=True))
    pool = sorted(set(cluster_of) - exclude)
    pick = sorted(rng.choice(pool, size=min(size, len(pool)), replace=False).tolist())
    return pick, cluster_of


def _universe(ctx: Context, features: list[str], exclude: set[str]) -> tuple[list[str], pl.DataFrame, dict[str, int]]:
    """The drift CpGs and a seeded sample of eligible non-promoter CpGs (``nonpromoter_universe``), less ``exclude``;
    their annotation (chromosome, position, island, probe type, promoter flag); and the non-promoter cluster labels."""
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(~pl.col("chr").is_in(["X", "Y"]))
    linked = set(pl.read_parquet(ann.links_450k).filter(pl.col("promoter"))["cpg"].to_list())
    pick, cluster_of = nonpromoter_universe(probes, linked, exclude, ctx.rng(150_000))
    cpgs = sorted((set(features) - exclude) | set(pick))
    return cpgs, probes.filter(pl.col("cpg").is_in(cpgs)), cluster_of


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # clock_drift's drift CpGs
    links = pl.read_parquet(AnnotationLayout(ctx.paths).links_450k)
    set_genes = set(membership["gene"].to_list())
    mitotic = mitotic_clocks.load(ctx.paths.root)
    clock_cpgs, other_clocks = excluded_cpgs(ctx, "clock_drift_classes")
    universe, probes, cluster_of = _universe(ctx, features, clock_cpgs | other_clocks)
    chemistry = dict(probes.select("cpg", "probe_type").iter_rows())
    entries, registry = group_registry(ctx, "clock_drift_classes")
    registry = flag_primary(registry, registry)
    primary = set(registry.filter(pl.col("primary")).select("series", "group").iter_rows())
    brain_of = {(s, g): b for s, g, b in registry.select("series", "group", "brain").iter_rows()}
    scale = str(ctx.scale.name)
    n_draws, n_boot = N_DRAWS[scale], N_BOOT[scale]
    rows, draw_rows, size_rows, full_rows = [], [], [], []
    for e in entries:
        inp, label, s, samples = e.series, e.label, e.s, e.samples
        if (inp.name, label) not in primary:
            continue
        sc = series_scores(inp, clock_cpgs)
        if sc.base is None and sc.mb is None:
            print(f"clock_drift_classes: {inp.name} / {label}: {sc.note}", file=sys.stderr)
            continue
        want = [c for c in universe if c not in inp.dropped]
        frame = pl.scan_parquet(inp.betas_path).filter(pl.col("cpg").is_in(want)).select("cpg", *samples).collect()
        x = frame.select(samples).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        kept = [c for c, k in zip(frame["cpg"].to_list(), ok, strict=True) if k]
        x = x[:, ok]
        sets = _sets(kept, membership)
        if GROUP not in sets:
            print(f"clock_drift_classes: {inp.name} / {label}: no housekeeping set", file=sys.stderr)
            continue
        member = np.any([m for m, _ in sets.values()], axis=0)
        hk_genes = sets[GROUP][1]
        op_genes = placebo_pool(kept, member, links, set_genes)
        linked, gene_cpgs = promoter_links(kept, links)
        is_np = np.array([c in cluster_of for c in kept]) & ~member
        np_units = units(np.array([cluster_of.get(c, -1) for c in kept]), np.flatnonzero(is_np))
        island = _islands(kept, probes)
        chem = np.array([chemistry.get(c, "NA") for c in kept])
        strata = match_strata(x, island, chem)
        rng = ctx.rng(160_000 + 100 * inp.gi + e.li)
        draws, pooled = class_draws(strata, hk_genes, op_genes, np_units, member, rng, n_draws, linked)
        gap = support_gap(strata, draws, pooled)
        if gap:
            print(f"clock_drift_classes: {inp.name} / {label}: {gap} strata without pooled background", file=sys.stderr)
        # one pass: every draw against its own background, and every draw against the pooled background
        a_all, b_all, bal_all = class_scores(x, island, [*draws, *[(sets, pooled) for sets, _ in draws]])
        a, b, bal = a_all[:n_draws], b_all[:n_draws], bal_all[:n_draws]
        a_pool, b_pool, bal_pool = a_all[n_draws:].mean(axis=0), b_all[n_draws:].mean(axis=0), bal_all[n_draws:]
        size_rows.append(
            {
                "accession": inp.gse,
                "series": inp.name,
                "group": label,
                "donors": len(samples),
                "universe": len(kept),
                "hk_genes": len(hk_genes),
                "op_genes_pool": len(op_genes),
                "np_clusters_pool": len(np_units),
                "per_class_cpgs": float(np.mean([sets[0].sum() for sets, _ in draws])),
                "background_cpgs": float(np.mean([bg.sum() for _, bg in draws])),
                "pooled_background_cpgs": int(pooled.sum()),
                "pooled_support_gap_strata": gap,
                "hk_max_per_linked_gene": int(
                    max(per_gene_counts(sets[0], gene_cpgs).max(initial=0) for sets, _ in draws)
                ),
                "op_max_per_linked_gene": int(
                    max(per_gene_counts(sets[1], gene_cpgs).max(initial=0) for sets, _ in draws)
                ),
                "np_max_per_cluster": int(max(per_gene_counts(sets[2], np_units).max(initial=0) for sets, _ in draws)),
                "balance_tvd": float(bal[:, 0].mean()),
                "balance_unusable": float(bal[:, 1].mean()),
                "balance_no_background": float(bal[:, 2].mean()),
                "pooled_balance_unusable": float(bal_pool[:, 1].mean()),
                "pooled_balance_no_background": float(bal_pool[:, 2].mean()),
            }
        )
        exclude = set(universe) | clock_cpgs | other_clocks | inp.dropped
        panel, _ = _panel(inp.betas_path, samples, exclude, ctx.rng(170_000 + 100 * inp.gi + e.li))
        level, comp = panel.mean(axis=1), np.abs(panel - 0.5).mean(axis=1)
        age = s["age"].to_numpy().astype(float)
        held = _held(s)
        quality, qdesc, _ = quality_terms(samples, inp.qc, base_columns=6 + held.shape[1])
        clk = group_clock_scores(samples, sc, mitotic)
        fixed_all = adjust_design(age, quality, held, level, comp)
        print(
            f"clock_drift_classes: {inp.name} / {label}: {len(samples)} arrays, {len(kept)} CpGs, "
            f"{size_rows[-1]['per_class_cpgs']:.0f} CpGs per class, pooled background {int(pooled.sum())}",
            file=sys.stderr,
            flush=True,
        )
        boot_rng = ctx.rng(180_000 + 100 * inp.gi + e.li)
        values = {}
        for clock in [c for c in CLOCKS if c in clk.columns]:
            v = clk[clock].cast(pl.Float64).to_numpy()
            use = np.isfinite(v)
            if use.sum() < MIN_ARRAYS or np.std(v[use]) == 0:
                continue
            values[clock] = v
            fixed = fixed_all[use]
            per_draw = draw_rs(v[use], a[:, use], b[:, use], fixed)
            with np.errstate(invalid="ignore"):
                r = np.nanmean(per_draw, axis=0)
            note = "" if np.isfinite(per_draw).any() else "undefined in every draw"
            con = contrasts(r)
            dcon = contrasts(per_draw)
            boot = bootstrap_draw_rs(v[use], a[:, use], b[:, use], fixed, boot_rng, n_boot)
            bcon = contrasts(boot)
            r_pool, _ = partial_rs(v[use], a_pool[use], np.column_stack([fixed, b_pool[use]]))
            con_pool = contrasts(r_pool)
            for j in range(n_draws):
                draw_rows.append(
                    {
                        "series": inp.name,
                        "group": label,
                        "clock": clock,
                        "draw": j,
                        **{f"r_{c}": float(per_draw[j, i]) for i, c in enumerate(CLASSES)},
                        "op_minus_np": float(dcon[j, 0]),
                        "hk_minus_op": float(dcon[j, 1]),
                    }
                )
            row = {
                "accession": inp.gse,
                "series": inp.name,
                "group": label,
                "brain": brain_of[(inp.name, label)],
                "clock": clock,
                "mitotic": clock in MITOTIC,
                "n": int(use.sum()),
                "covariates": qdesc + (", batch/center" if held.shape[1] else "") + ", panel level/compression, B",
                "export_check": sc.note,
                "export_identity": sc.identity,
                "interval_scope": INTERVAL_SCOPE,
                "note": note,
            }
            for i, c in enumerate(CLASSES):
                row[f"r_{c}"] = float(r[i])
                row[f"r_{c}_lo"], row[f"r_{c}_hi"] = (float(q) for q in np.nanpercentile(boot[:, i], [2.5, 97.5]))
                row[f"r_{c}_pooled_bg"] = float(r_pool[i])
            for i, name in enumerate(("op_minus_np", "hk_minus_op")):
                row[name] = float(con[i])
                row[f"{name}_lo"], row[f"{name}_hi"] = (float(q) for q in np.nanpercentile(bcon[:, i], [2.5, 97.5]))
                row[f"{name}_se"] = float(np.nanstd(bcon[:, i], ddof=1))
                row[f"{name}_draw_median"] = float(np.nanmedian(dcon[:, i]))
                row[f"{name}_draw_min"] = float(np.nanmin(dcon[:, i]))
                row[f"{name}_draw_max"] = float(np.nanmax(dcon[:, i]))
                row[f"{name}_draw_sign_changes"] = bool(np.nanmin(dcon[:, i]) < 0 < np.nanmax(dcon[:, i]))
                row[f"{name}_pooled_bg"] = float(con_pool[i])
            row["r_B"] = corr(residual(v[use], age[use]), b.mean(axis=0)[use])
            rows.append(row)
        full_clocks = {k: v for k, v in values.items() if k in FULL_CLOCKS}
        if full_clocks and not brain_of[(inp.name, label)] and FULL_BOOT[scale]:
            started = time.perf_counter()
            res = donor_resample(
                x,
                island,
                chem,
                (hk_genes, op_genes, np_units),
                member,
                linked,
                full_clocks,
                fixed_all,
                ctx.rng(200_000 + 100 * inp.gi + e.li),
                FULL_BOOT[scale],
                FULL_DRAWS[scale],
            )
            for clock, rr in res.items():
                fc = contrasts(rr)
                full_rows.append(
                    {
                        "accession": inp.gse,
                        "series": inp.name,
                        "group": label,
                        "clock": clock,
                        "resamples": FULL_BOOT[scale],
                        "draws_per_resample": FULL_DRAWS[scale],
                        **{f"r_{c}_full_se": float(np.nanstd(rr[:, i], ddof=1)) for i, c in enumerate(CLASSES)},
                        **{
                            f"{name}_full_{stat}": float(val)
                            for i, name in enumerate(("op_minus_np", "hk_minus_op"))
                            for stat, val in zip(
                                ("se", "lo", "hi"),
                                (np.nanstd(fc[:, i], ddof=1), *np.nanpercentile(fc[:, i], [2.5, 97.5])),
                                strict=True,
                            )
                        },
                    }
                )
            print(
                f"clock_drift_classes: {inp.name} / {label}: donor resampling done in "
                f"{time.perf_counter() - started:.0f} s",
                file=sys.stderr,
                flush=True,
            )
    if not rows:
        print("clock_drift_classes: nothing to score", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_drift_classes")
    table = pl.DataFrame(rows)
    full = pl.DataFrame(full_rows) if full_rows else None
    if full is not None:
        table = table.join(
            full.drop("accession", "resamples", "draws_per_resample"), on=["series", "group", "clock"], how="left"
        )
        write_table(full, out / "donor_resampling", csv_max_rows=10_000)
    write_table(table, out / "associations", csv_max_rows=10_000)
    write_table(pl.DataFrame(draw_rows), out / "draws", csv_max_rows=0)
    write_table(pl.DataFrame(size_rows), out / "sizes", csv_max_rows=10_000)
    summary = summarise(table, ctx.rng(190_000), n_boot)
    write_table(summary, out / "summary", csv_max_rows=10_000)
    lacking = sorted({r["series"] for r in rows if not r["export_identity"]})
    print(f"clock_drift_classes: predictions without export identity: {', '.join(lacking) or 'none'}", file=sys.stderr)
    return [draw(table.filter(~pl.col("brain")), summary, ctx.figure_path("clock_drift_classes"))]


def summarise(table: pl.DataFrame, rng: np.random.Generator, n_boot: int = 1000) -> pl.DataFrame:
    """Per clock over non-brain primary groups (one per accession): medians of each class's r and of both contrasts,
    the random-effects mean of each contrast with its normal 95 % interval and the Hartung-Knapp interval (``_hk_``),
    a percentile interval of the median contrast over resampled accessions, and how many accessions' bootstrap
    intervals lie below or above zero. These intervals are conditional on the constructed scores
    (``interval_scope``). Where the donor resampling of the full procedure ran (``_full_se``), the same
    random-effects mean and both intervals with its SEs (``_re_full``), over the groups it covered."""
    rows = []
    t = table.filter(~pl.col("brain"))
    for clock in [c for c in CLOCKS if c in t["clock"].to_list()]:
        d = t.filter(pl.col("clock") == clock)
        row: dict = {
            "clock": clock,
            "mitotic": clock in MITOTIC,
            "accessions": d["accession"].n_unique(),
            "interval_scope": INTERVAL_SCOPE,
        }
        for c in CLASSES:
            row[f"r_{c}_median"] = float(np.nanmedian(d[f"r_{c}"].to_numpy()))
            if f"r_{c}_pooled_bg" in d.columns:
                row[f"r_{c}_pooled_bg_median"] = float(np.nanmedian(d[f"r_{c}_pooled_bg"].to_numpy()))
        for name in ("op_minus_np", "hk_minus_op"):
            est, se = d[name].to_numpy(), d[f"{name}_se"].to_numpy()
            row[f"{name}_median"] = float(np.nanmedian(est))
            mu, mse, tau2 = random_effects(est, se)
            row[f"{name}_re"], row[f"{name}_re_lo"], row[f"{name}_re_hi"] = mu, mu - 1.96 * mse, mu + 1.96 * mse
            _, row[f"{name}_re_hk_lo"], row[f"{name}_re_hk_hi"] = hartung_knapp(est, se)
            row[f"{name}_tau2"] = tau2
            ok = est[np.isfinite(est)]
            meds = [np.median(rng.choice(ok, ok.size)) for _ in range(n_boot)] if ok.size else [np.nan]
            row[f"{name}_median_lo"], row[f"{name}_median_hi"] = (float(q) for q in np.percentile(meds, [2.5, 97.5]))
            row[f"{name}_below0"] = int((d[f"{name}_hi"] < 0).sum())
            row[f"{name}_above0"] = int((d[f"{name}_lo"] > 0).sum())
            if f"{name}_pooled_bg" in d.columns:
                row[f"{name}_pooled_bg_median"] = float(np.nanmedian(d[f"{name}_pooled_bg"].to_numpy()))
            f = d.filter(pl.col(f"{name}_full_se").is_not_null()) if f"{name}_full_se" in d.columns else d.clear()
            if f.height:
                fe, fs = f[name].to_numpy(), f[f"{name}_full_se"].to_numpy()
                mu, mse, _ = random_effects(fe, fs)
                row[f"{name}_full_groups"] = f.height
                row[f"{name}_re_full"], row[f"{name}_re_full_lo"], row[f"{name}_re_full_hi"] = (
                    mu,
                    mu - 1.96 * mse,
                    mu + 1.96 * mse,
                )
                _, row[f"{name}_re_full_hk_lo"], row[f"{name}_re_full_hk_hi"] = hartung_knapp(fe, fs)
                row[f"{name}_full_se_ratio"] = float(np.nanmedian(fs / f[f"{name}_se"].to_numpy()))
        rows.append(row)
    return pl.DataFrame(rows)


def draw(table: pl.DataFrame, summary: pl.DataFrame, path: Path) -> Path:
    """Per clock: each class's r per primary group (a), and the two contrasts per group with the random-effects mean
    and its interval (b, c)."""
    plots.apply_style()
    clocks = [c for c in CLOCKS if c in summary["clock"].to_list()]
    fig, axes = plots.plt.subplots(1, 3, figsize=(15, 0.45 * len(clocks) + 1.8), sharey=True)
    y = np.arange(len(clocks))[::-1]
    ax = axes[0]
    for i, c in enumerate(CLASSES):
        for yy, clock in zip(y, clocks, strict=True):
            v = table.filter(pl.col("clock") == clock)[f"r_{c}"].drop_nans().to_numpy()
            ax.scatter(v, np.full(v.size, yy + (i - 1) * 0.22), s=8, alpha=0.5, color=plots.PALETTE[i])
            if v.size:
                ax.plot([np.median(v)] * 2, [yy + (i - 1) * 0.22 - 0.1, yy + (i - 1) * 0.22 + 0.1], color=plots.INK)
        ax.scatter([], [], color=plots.PALETTE[i], label=LABELS[c])
    ax.axvline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_yticks(y, clocks, fontsize=8)
    ax.set_xlabel("partial r, clock with class deviation | B and covariates")
    ax.set_title("a  per class, primary non-brain groups", loc="left", fontsize=9)
    ax.legend(fontsize=7, loc="lower left")
    for ax, name, title in (
        (axes[1], "op_minus_np", "b  other promoters minus non-promoters"),
        (axes[2], "hk_minus_op", "c  housekeeping minus other promoters"),
    ):
        for yy, clock in zip(y, clocks, strict=True):
            d = table.filter(pl.col("clock") == clock)
            ax.scatter(d[name], np.full(d.height, yy + 0.15), s=8, alpha=0.5, color=plots.MUTED)
            srow = summary.filter(pl.col("clock") == clock).row(0, named=True)
            ax.errorbar(
                srow[f"{name}_re"],
                yy - 0.15,
                xerr=[[srow[f"{name}_re"] - srow[f"{name}_re_lo"]], [srow[f"{name}_re_hi"] - srow[f"{name}_re"]]],
                fmt="o",
                color=plots.PALETTE[0],
                ms=4,
            )
        ax.axvline(0, color=plots.MUTED, lw=0.8, ls="--")
        ax.set_xlabel("difference in partial r (grey: groups; blue: random-effects mean, 95 % CI)")
        ax.set_title(title, loc="left", fontsize=9)
    fig.tight_layout()
    return plots.save(fig, path)
