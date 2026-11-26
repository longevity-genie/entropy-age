"""Can a broad array shift make the clock-shield link? Confounder and placebo controls for ``clock_drift``.

``clock_drift`` finds that a mitotic clock's acceleration tracks the per-person shield share F (r about -0.5 outside
the brain; HypoClock, whose value falls with division, about +0.3). F is a ratio: with A the housekeeping promoters'
deviation from the leave-one-out reference and B that of their look-alikes, ``F = (A - B) / B``. A shift of the
whole array (a methylation level change, beta compression toward 0.5, chip, signal or detection quality, cell
mixture) that moves A and B by different amounts moves F, and the same shift can move a clock (compression raises
low-beta promoter clock CpGs and lowers high-beta HypoClock CpGs at once). Codex's second review (item 3) proposed
the check built here; its fourth review (findings 1, 3, 5, 6 and 7) shaped the placebo comparison, the panel, the
primary groups, the guards on partial correlations and the array checks; its fifth review (findings 3 to 7) added the
gene-balanced placebo pairs, the exclusion of preprocessing probes, the group registry, one scaled design for fit and
rank, and the export check.

Per series group, with the arrays, CpGs and clock scores of ``clock_drift`` (``series.selection``, the same
availability rule, every mitotic and EpiTOC1 CpG kept out of the drift CpGs, sensitivity configs flagged):

0. **Arrays.** The clock track's QC (``sesame_clock``) is required for every series: a series without it is skipped
   with a logged reason, and only arrays passing it enter. For every sample id in both rebuilds, the physical array
   (chip and position from each track's ``qc.csv``, the file prefix where either is missing) must be the same; a
   sample whose two exports read different arrays is dropped, logged and counted (``array_mismatch_dropped``).
   The clock scores are then checked against the current rebuild (``export_check``): where the beta export carries
   a manifest of physical arrays (``scripts/export_clock_cpgs.py`` writes ``<GSE>_cpgs.manifest.csv`` since review
   5), a sample whose exported array is not the one the current clock-track QC keeps is dropped
   (``export_dropped``); without one, the export is flagged stale (``export_stale``, ``export_check``) when it is
   older than the clock-track QC, the scores are older than the export, the scores and the export hold different
   sample ids, or a passing array is missing from the export.
1. **Drift measures.** ``loo_scores`` per person: ``P_housekeeping`` (the matched gap), ``B_housekeeping`` and
   ``F_housekeeping``; ``A = P + B`` is the housekeeping CpGs' own stratum-weighted deviation.
2. **A third, disjoint probe panel.** CpGs of the same rebuild drawn in a seeded order, kept when measured in every
   array of the group, until 20,000 are found or the pool runs out. The pool leaves out the drift CpGs, every
   mitotic clock's CpGs, every prediction input of every clock in the exported biolearn gallery (each clock named in
   the clock scores, at the commit ``scripts/clocks_biolearn.py`` pins: its coefficient file, a PC clock's centre and
   rotation files, and the data files its preprocessing reads, which for DunedinPACE is the 20,000-probe
   gold-standard list of its quantile normalisation; ``gallery_inputs``) and the CpGs of the local clock files
   (``clock_cpgs.clock_table``). The probe ids per group are saved (``panel_probes``) and checked against the
   exclusion after writing (``panel_overlap``; the run stops on any overlap). Per person:
   the panel's mean beta (``level``), ``compression`` = mean ``|beta - 0.5|`` (lower is more compressed toward 0.5)
   and the SD of its betas (reported).
3. **Covariates.** Age; chip dummies, log mean intensity and detection fraction of the rebuild the drift CpGs come
   from (``reprocessed.read_qc`` through ``clock_disagreement.quality_terms``, chip dropped when it would leave under
   10 residual degrees of freedom); ``batch`` or ``center`` dummies where the series description holds them. The
   series groups carry no cell fractions (only EPIC-Italy has EpiDISH fractions, and it is not a series group), so
   none enter; level and compression stand in for a mixture that moves the whole array.
4. **Same adjustment on both sides.** The clock's values, A, B and F are all residualised on age + covariates +
   level + compression (age sits in the joint design, so the clock enters with its own scale), then correlated; the
   raw correlations use the acceleration (residual on age, as in ``clock_drift``), undefined when age explains the
   clock (residual norm at most 1e-8 of the clock's centred norm, ``age_note``). Reported per group and clock:
   ``r_F_raw`` (``clock_drift``'s value), ``r_F_adj``, ``r_A_given_B`` (acceleration with A at a given B, the
   adjustment set plus B), ``r_level``, ``r_compression``, ``r_sd`` (the acceleration against the panel indices),
   ``r_F_level`` and ``r_F_compression`` (F against them), and n. Each partial correlation checks, on that clock's
   own rows, the rank of the column-scaled design and the residual degrees of freedom (at least 10), fits with the
   same scaled design and tolerance (``scaled_basis``), and calls a vector explained by the design when its residual
   norm is at most 1e-8 of its centred norm; such an estimate is
   NaN and its reason is written (``r_F_adj_note``, ``r_A_given_B_note``, ``adj_note`` per placebo draw).
5. **Placebo sets, each scored beside a matched housekeeping subset.** K = 20 draws per group. Matching strata:
   island status x decile of the group's mean beta x decile of the SD of beta across the group's arrays x probe
   chemistry (Infinium I or II, ``probe_type`` of the 450k annotation). A placebo takes the housekeeping set's count
   in every stratum as its quota and is drawn gene by gene from promoter CpGs of genes in none of the three gene
   sets (the pool ``placebo_sets`` draws from), restricted to the drift CpGs of the group; a stratum the pool cannot
   fill stays short (``placebo_fill``). For each draw a housekeeping subset is drawn at random with the placebo's
   attained count in every stratum. The pair is scored in the leave-one-out pass (``loo_drift``: ``loo_scores``'
   reference and deviations) against one background, the CpGs outside every housekeeping CpG and outside the
   placebo, with one set of stratum weights (the pair's pooled count in each scoring stratum that holds both sets and
   background), so B is shared and the two F differ only in their sets. Per draw: both r with acceleration, raw and
   adjusted, and the pair's balance under the scoring reference (total-variation distance of their scoring-stratum
   counts, share of their CpGs in strata one of them lacks), averaged over the left-out people. Per group and clock:
   the housekeeping r on the matched subsets (median over draws, ``hk_matched_raw`` and ``hk_matched_adj``) beside
   the placebo distribution (median and 95th percentile of ``|r|``), the share of draws whose placebo r reaches its
   paired housekeeping r in that r's direction, and whether ``|hk_matched_adj|`` exceeds the placebo 95th
   percentile. Per-stratum quotas and mean fills are written per group (``placebo_strata``).
   **Gene-balanced pairs** (``scheme`` ``gene``, columns ending ``_gb``): the same K draws with both targets sampled
   by gene (``gene_balanced_pair``): one CpG of each housekeeping gene at random, a placebo filled to that draw's
   stratum counts with at most one CpG per pool gene, and the housekeeping draw thinned to the placebo's attained
   counts. Both sets then hold one CpG per gene and the same number of genes, so equal probe counts mean equal numbers
   of promoters. Reported for both schemes: the housekeeping margin on this common support (``margin_raw``,
   ``margin_adj``: median over draws of ``sign(hk r) * (hk r - placebo r)``), the genes per set and the fill.
   Expression and chromatin state are not matched (no per-group data here).
6. **Primary groups.** A registry of every group (``group_registry``, written as ``groups``) is built from arrays
   that pass the clock track's QC, read the same physical array in both rebuilds and are held by the drift betas,
   before any prediction file is opened. One primary group per accession (study): the non-sensitivity eligible group
   (at least 25 arrays with an age) with the most donors (a pooled group wins a tie). It stays primary when a clock
   has no estimate there; every primary group and clock without a row is listed with its reason (``missing``).
   Every row stays in ``correlations`` with a ``primary`` flag. ``summary`` gives, per mitotic clock, medians over
   primary non-brain groups with the count of groups exceeding the placebo band (both schemes) and of positive
   gene-balanced margins, and medians over accessions of each accession's median over its non-sensitivity non-brain
   groups.

Writes ``clock_drift_controls/correlations``, ``placebo`` (per group, clock, scheme and draw), ``placebo_strata``,
``panel_probes``, ``groups``, ``missing``, ``summary`` and a figure. A Slurm job: the leave-one-out scoring runs
once per group, with the placebo pairs inside it.
"""

from __future__ import annotations

import os
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments import mitotic_clocks
from entropy_age.experiments.clock_cpgs import clock_table
from entropy_age.experiments.clock_disagreement import quality_terms
from entropy_age.experiments.clock_drift import CLOCK_SCORES, CLOCK_TRACK, GROUP, MIN_ARRAYS, SERIES, residual
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import SENSITIVITY, selection, specs
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import contrast_weights, deviations, fit_reference
from entropy_age.stochasticity import baseline_strata

PANEL = 20_000
GENE_CAP = 1  # CpGs per gene in the gene-balanced placebo pairs (review 5, finding 3)
N_PLACEBO = {"paper": 20, "quick": 5}
MITOTIC = ("EpiTOC1", "epiTOC2", "epiTOC3", "HypoClock", "RepliTali", "epiCMIT")
HELD = ("batch", "center")  # technical columns a series description may hold
BRAIN = re.compile(r"brain|cortex|cerebell|hippocamp|amygdala|neuron|gyrus|striatum|caudate|glia|entorhinal", re.I)
POOLED = r"(?i)\ball\b|every|pooled"  # a group pooling a series' strata wins a tie on donors
MIN_RESID_DF = 10  # as clock_disagreement's chip guard
RESID_TOL = 1e-8  # a residual at most this share of the centred norm counts as explained by the design
BIOLEARN_SCRIPT = Path("scripts/clocks_biolearn.py")
CPG_ID = re.compile(r"c[gh][\d.]")
UNREAD_PREPROCESS = "<preprocess not read>"


# ----------------------------------------------------------------------------- pure helpers
def panel_indices(x: np.ndarray) -> dict[str, np.ndarray]:
    """Per person (rows of ``x``, people x CpGs, no NaN): mean beta, mean ``|beta - 0.5|`` and SD of the betas."""
    x = np.asarray(x, dtype=float)
    return {"level": x.mean(axis=1), "compression": np.abs(x - 0.5).mean(axis=1), "sd": x.std(axis=1)}


def corr(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def partial_r_checked(
    x: np.ndarray, y: np.ndarray, given: np.ndarray, names: tuple[str, str] = ("x", "y"), min_df: int = MIN_RESID_DF
) -> tuple[float, str]:
    """Correlation of ``x`` and ``y`` after both are regressed on an intercept and ``given``, with the reason when it
    is undefined: a non-finite value, fewer than ``min_df`` residual degrees of freedom (rows minus the rank of the
    column-scaled design), or a vector the design explains (residual norm at most ``RESID_TOL`` of its centred
    norm). Returns ``(r, "")`` or ``(nan, reason)``."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    given = np.asarray(given, dtype=float).reshape(len(x), -1)
    if not (np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(given).all()):
        return float("nan"), "non-finite value"
    basis, rank = scaled_basis(np.column_stack([np.ones(len(x)), given]))
    df = len(x) - rank
    if df < min_df:
        return float("nan"), f"residual df {df} < {min_df} (n {len(x)}, design rank {rank})"
    res = []
    for v, name in zip((x, y), names, strict=True):
        size = float(np.linalg.norm(v - v.mean()))
        if size == 0:
            return float("nan"), f"{name} constant"
        r = v - basis @ (basis.T @ v)
        if np.linalg.norm(r) <= RESID_TOL * size:
            return float("nan"), f"{name} explained by the design"
        res.append(r)
    return corr(res[0], res[1]), ""


def scaled_basis(d: np.ndarray) -> tuple[np.ndarray, int]:
    """An orthonormal basis of the column space of ``d`` after each column is scaled to unit norm, and its rank: the
    singular values of the scaled design above ``matrix_rank``'s default tolerance. Rank check and fit share it, so a
    tiny column the rank counts is also used in the fit (review 5, finding 6)."""
    d = np.asarray(d, dtype=float)
    scale = np.linalg.norm(d, axis=0)
    u, sv, _ = np.linalg.svd(d / np.where(scale > 0, scale, 1.0), full_matrices=False)
    tol = sv.max() * max(d.shape) * np.finfo(float).eps if sv.size else 0.0
    rank = int((sv > tol).sum())
    return u[:, :rank], rank


def age_residual_checked(v: np.ndarray, age: np.ndarray) -> tuple[np.ndarray, str]:
    """The clock's acceleration (``residual`` on age) and, when age explains the clock (the residual's norm at most
    ``RESID_TOL`` of the clock's centred norm), a reason; the acceleration is then all NaN."""
    v = np.asarray(v, dtype=float)
    acc = residual(v, np.asarray(age, dtype=float))
    size = float(np.linalg.norm(v - v.mean()))
    if size == 0 or np.linalg.norm(acc) <= RESID_TOL * size:
        return np.full(v.size, np.nan), "clock explained by age"
    return acc, ""


def clock_stats(
    v: np.ndarray, age: np.ndarray, design: np.ndarray, f: np.ndarray, a: np.ndarray, b: np.ndarray, ind: dict
) -> dict[str, float | str]:
    """One clock's correlations in one group (rows already restricted to the arrays the clock scored). Raw: the
    acceleration (``age_residual_checked``) with F, B and the panel indices. Adjusted: the original clock values and F
    (or A given B) regressed on the joint design, which holds age, so the clock's own scale decides whether the design
    explains it (review 5, finding 6)."""
    acc, note_age = age_residual_checked(v, age)
    r_adj, note_adj = partial_r_checked(v, f, design, ("clock", "F"))
    r_ab, note_ab = partial_r_checked(v, a, np.column_stack([design, b]), ("clock", "A"))
    ok = not note_age
    return {
        "age_note": note_age,
        "r_F_raw": corr(acc, f) if ok else float("nan"),
        "r_B_raw": corr(acc, b) if ok else float("nan"),
        "r_F_adj": r_adj,
        "r_F_adj_note": note_adj,
        "r_B_adj": partial_r(v, b, design),
        "r_A_given_B": r_ab,
        "r_A_given_B_note": note_ab,
        "r_level": corr(acc, ind["level"]) if ok else float("nan"),
        "r_compression": corr(acc, ind["compression"]) if ok else float("nan"),
        "r_sd": corr(acc, ind["sd"]) if ok else float("nan"),
        "r_F_level": corr(f, ind["level"]),
        "r_F_compression": corr(f, ind["compression"]),
    }


def partial_r(x: np.ndarray, y: np.ndarray, given: np.ndarray) -> float:
    return partial_r_checked(x, y, given)[0]


def matched_placebo(
    strata: np.ndarray,
    real: np.ndarray,
    pool_genes: list[np.ndarray],
    rng: np.random.Generator,
    per_gene: int | None = None,
) -> tuple[np.ndarray, list[np.ndarray], float]:
    """A placebo set with ``real``'s count of CpGs in every stratum, filled gene by gene from ``pool_genes`` (index
    arrays of candidate CpGs, one per gene) in random order; a gene adds those of its CpGs whose stratum still has
    room, at most ``per_gene`` of them (taken in random order) when given. Returns (membership mask, the genes' index
    arrays as used, share of the quota filled)."""
    levels, codes = np.unique(strata, return_inverse=True)
    quota = np.bincount(codes[real], minlength=levels.size)
    want = int(quota.sum())
    taken = np.zeros(strata.size, dtype=bool)
    genes = []
    for g in rng.permutation(len(pool_genes)):
        if quota.sum() == 0:
            break
        cand = pool_genes[g] if per_gene is None else rng.permutation(pool_genes[g])
        use = []
        for j in cand:
            if per_gene is not None and len(use) >= per_gene:
                break
            if not taken[j] and not real[j] and quota[codes[j]] > 0:
                quota[codes[j]] -= 1
                taken[j] = True
                use.append(j)
        if use:
            genes.append(np.array(use, dtype=int))
    return taken, genes, float(taken.sum() / want) if want else float("nan")


def gene_balanced(genes: list[np.ndarray], n: int, rng: np.random.Generator, per_gene: int = 1) -> np.ndarray:
    """A mask over ``n`` CpGs with at most ``per_gene`` CpGs of every gene (index arrays), drawn at random; a CpG
    shared by two genes is taken once."""
    out = np.zeros(n, dtype=bool)
    for g in genes:
        free = np.setdiff1d(np.asarray(g, dtype=int), np.flatnonzero(out))
        if free.size:
            out[rng.choice(free, min(per_gene, free.size), replace=False)] = True
    return out


def per_gene_counts(mask: np.ndarray, genes: list[np.ndarray]) -> np.ndarray:
    """How many CpGs of ``mask`` each gene holds (genes holding none left out)."""
    counts = np.array([int(mask[np.asarray(g, dtype=int)].sum()) for g in genes])
    return counts[counts > 0]


def gene_balanced_pair(
    strata: np.ndarray,
    hk_genes: list[np.ndarray],
    pool_genes: list[np.ndarray],
    rng: np.random.Generator,
    per_gene: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Both targets sampled by gene (review 5, finding 3): a housekeeping draw with at most ``per_gene`` CpGs per
    housekeeping gene (``gene_balanced``), a placebo filled to its stratum counts with at most ``per_gene`` CpGs per
    pool gene, and the housekeeping draw thinned to the placebo's attained counts (the common support,
    ``matched_subset``). With ``per_gene`` 1 both sets hold one CpG per gene and as many genes as CpGs. The background
    leaves out every housekeeping CpG and the placebo. Returns ``(subset, placebo, background, fill)``."""
    n = strata.size
    hk_all = np.zeros(n, dtype=bool)
    for g in hk_genes:
        hk_all[np.asarray(g, dtype=int)] = True
    hk = gene_balanced(hk_genes, n, rng, per_gene)
    placebo, _, fill = matched_placebo(strata, hk, pool_genes, rng, per_gene)
    return matched_subset(strata, hk, placebo, rng), placebo, ~(hk_all | placebo), fill


def matched_subset(strata: np.ndarray, real: np.ndarray, placebo: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """CpGs of ``real`` drawn at random, as many in every stratum as ``placebo`` holds there (the attained quotas)."""
    levels, codes = np.unique(strata, return_inverse=True)
    want = np.bincount(codes[placebo], minlength=levels.size)
    idx = np.flatnonzero(real)
    idx = idx[np.argsort(codes[idx], kind="stable")]
    bounds = np.searchsorted(codes[idx], np.arange(levels.size + 1))
    out = np.zeros(strata.size, dtype=bool)
    for s in np.flatnonzero(want):
        cand = idx[bounds[s] : bounds[s + 1]]
        if cand.size < want[s]:
            raise ValueError(f"stratum {levels[s]}: {want[s]} wanted, {cand.size} in the set")
        out[rng.choice(cand, int(want[s]), replace=False)] = True
    return out


def paired_weights(
    a: np.ndarray, b: np.ndarray, background: np.ndarray, strata: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Weights scoring two disjoint sets against one background with one set of stratum weights: within stratum ``k``
    each set's CpGs get ``w_k / n_set_k`` and the background ``-w_k / n_bg_k``, with ``w_k`` the pair's pooled count in
    ``k`` over the strata holding both sets and background. ``strata`` may be labels or non-negative integer codes.
    Returns ``(pa, pb, bw, balance)``: ``P_a = d @ pa``, ``P_b = d @ pb``, the shared ``B = d @ bw``, and ``balance``
    with the total-variation distance of the two sets' stratum counts (``tvd``) and the share of their CpGs in strata
    left out (``unusable``)."""
    a, b, background = (np.asarray(v, dtype=bool) for v in (a, b, background))
    if (a & b).any() or (background & (a | b)).any():
        raise ValueError("paired_weights: the two sets and the background must be disjoint")
    strata = np.asarray(strata)
    if np.issubdtype(strata.dtype, np.integer) and (strata.size == 0 or strata.min() >= 0):
        codes = strata
    else:
        codes = np.unique(strata, return_inverse=True)[1]
    k = int(codes.max()) + 1 if codes.size else 0
    n_a = np.bincount(codes[a], minlength=k).astype(float)
    n_b = np.bincount(codes[b], minlength=k).astype(float)
    n_bg = np.bincount(codes[background], minlength=k).astype(float)
    usable = (n_a > 0) & (n_b > 0) & (n_bg > 0)
    if not usable.any():
        raise ValueError("paired_weights: no stratum holds both sets and background")
    w = np.where(usable, n_a + n_b, 0.0)
    w = w / w.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        per_a = np.where(usable, w / n_a, 0.0)
        per_b = np.where(usable, w / n_b, 0.0)
        per_bg = np.where(usable, w / n_bg, 0.0)
    bw = np.where(background, per_bg[codes], 0.0)
    pa = np.where(a, per_a[codes], 0.0) - bw
    pb = np.where(b, per_b[codes], 0.0) - bw
    total = n_a.sum() + n_b.sum()
    balance = {
        "tvd": float(0.5 * np.abs(n_a / max(n_a.sum(), 1) - n_b / max(n_b.sum(), 1)).sum()),
        "unusable": float(1 - (n_a[usable].sum() + n_b[usable].sum()) / total) if total else float("nan"),
    }
    return pa, pb, bw, balance


def placebo_pair(
    strata: np.ndarray, hk: np.ndarray, pool_genes: list[np.ndarray], rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """One placebo draw (``matched_placebo``), the housekeeping subset matched to its attained quotas
    (``matched_subset``) and the background both are scored against: every CpG outside the whole housekeeping set and
    outside the placebo. Returns ``(subset, placebo, background, fill)``."""
    placebo, _, fill = matched_placebo(strata, hk, pool_genes, rng)
    return matched_subset(strata, hk, placebo, rng), placebo, ~(hk | placebo), fill


def placebo_pool(kept: list[str], member: np.ndarray, links: pl.DataFrame, set_genes: set[str]) -> list[np.ndarray]:
    """Per gene in none of the gene sets: the indices (in ``kept``) of its promoter CpGs outside every set."""
    index = {c: j for j, c in enumerate(kept)}
    rows = links.filter(pl.col("promoter") & pl.col("cpg").is_in(kept) & ~pl.col("symbol").is_in(list(set_genes)))
    inside = set(np.flatnonzero(member).tolist())
    out = []
    for _, g in rows.group_by("symbol", maintain_order=True):
        idx = sorted({index[c] for c in g["cpg"].to_list()} - inside)
        if idx:
            out.append(np.array(idx, dtype=int))
    out.sort(key=lambda a: int(a[0]))  # a fixed order before the seeded shuffle
    return out


def loo_drift(
    x: np.ndarray,
    island: np.ndarray,
    sets: dict,
    drop: dict[str, np.ndarray] | None = None,
    pairs: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] | None = None,
) -> dict[str, np.ndarray]:
    """``loo_scores``' ``P``, ``B`` and ``F`` of every set (no components), computed once per person for all sets:
    the same leave-one-out reference (``fit_reference``), deviations and stratum weights (``contrast_weights``), with
    the entropy, relative-entropy and signed gaps left out. ``drop`` maps a set to CpGs left out of its contrast
    altogether (neither in it nor among its look-alikes). ``pairs`` maps a name to ``(a, b, background)``, two sets
    scored against one background with shared stratum weights (``paired_weights``): ``P_<name>_a``, ``P_<name>_b``,
    the shared ``B_<name>``, ``F_<name>_a``, ``F_<name>_b`` and per person ``tvd_<name>`` and ``unusable_<name>``."""
    drop = drop or {}
    pairs = pairs or {}
    member = np.any([s for s, _ in sets.values()], axis=0)
    names = list(sets)
    out = {f"{m}_{g}": np.empty(x.shape[0]) for g in names for m in "PBF"}
    for g in pairs:
        for key in (f"P_{g}_a", f"P_{g}_b", f"B_{g}", f"tvd_{g}", f"unusable_{g}"):
            out[key] = np.empty(x.shape[0])
    for i in range(x.shape[0]):
        others = np.arange(x.shape[0]) != i
        ref = fit_reference(x[others], island, ~member)
        codes = np.unique(ref.strata, return_inverse=True)[1]
        d = deviations(x[i], ref.level)
        for g in names:
            keep = ~drop[g] if g in drop else np.ones(d.size, dtype=bool)
            p, b = contrast_weights(np.asarray(sets[g][0][keep], dtype=float), codes[keep])
            out[f"P_{g}"][i], out[f"B_{g}"][i] = d[keep] @ p, d[keep] @ b
        for g, (a, b_set, background) in pairs.items():
            pa, pb, bw, bal = paired_weights(a, b_set, background, codes)
            out[f"P_{g}_a"][i], out[f"P_{g}_b"][i], out[f"B_{g}"][i] = d @ pa, d @ pb, d @ bw
            out[f"tvd_{g}"][i], out[f"unusable_{g}"][i] = bal["tvd"], bal["unusable"]
    for g in names:
        out[f"F_{g}"] = out[f"P_{g}"] / out[f"B_{g}"]
    for g in pairs:
        out[f"F_{g}_a"] = out[f"P_{g}_a"] / out[f"B_{g}"]
        out[f"F_{g}_b"] = out[f"P_{g}_b"] / out[f"B_{g}"]
    return out


def adjust_design(age: np.ndarray, quality: np.ndarray, held: np.ndarray, level: np.ndarray, comp: np.ndarray):
    return np.column_stack([age, quality, held, level, comp])


def array_ids(qc: pl.DataFrame) -> pl.DataFrame:
    """``sample`` and its physical ``array``: ``<chip>_<position>`` from a QC table, the file prefix where either is
    missing."""
    cols = set(qc.columns)
    none = pl.lit(None, dtype=pl.String)
    chip = pl.col("chip").cast(pl.String) if "chip" in cols else none
    pos = pl.col("position").cast(pl.String) if "position" in cols else none
    prefix = pl.col("prefix").cast(pl.String) if "prefix" in cols else none
    both = chip.is_not_null() & pos.is_not_null() & (chip != "NA") & (pos != "NA")
    array = pl.when(both).then(pl.concat_str([chip, pos], separator="_")).otherwise(prefix)
    return qc.select(pl.col("sample").cast(pl.String), array=array)


def matched_arrays(series_qc: pl.DataFrame, clock_qc: pl.DataFrame) -> tuple[set[str], list[str]]:
    """Sample ids in both QC tables whose physical arrays agree, and a reason per shared sample id whose arrays
    differ or are unknown in either table."""
    j = array_ids(series_qc).join(array_ids(clock_qc), on="sample", how="inner", suffix="_clock")
    agree = pl.col("array").is_not_null() & (pl.col("array") == pl.col("array_clock")).fill_null(False)
    same = j.filter(agree)
    other = j.filter(~agree)
    reasons = [
        f"{s}: drift array {a}, clock array {c}"
        for s, a, c in other.select("sample", "array", "array_clock").iter_rows()
    ]
    return set(same["sample"].to_list()), reasons


def manifest_path(export: Path) -> Path:
    """The export manifest beside a beta export (``<GSE>_cpgs.parquet`` -> ``<GSE>_cpgs.manifest.csv``):
    ``scripts/export_clock_cpgs.py`` writes each exported sample's physical array there."""
    return export.with_name(export.name.removesuffix(".parquet") + ".manifest.csv")


def export_check(
    *,
    qc_mtime: float,
    export_mtime: float | None,
    clocks_mtime: float,
    export_samples: set[str] | None,
    clock_samples: set[str],
    passing: set[str],
    manifest: pl.DataFrame | None = None,
    current: pl.DataFrame | None = None,
) -> tuple[set[str], str, bool]:
    """Is the clock-score export the current rebuild's (review 5, finding 7)? With a ``manifest`` (``sample`` and
    ``array``, written at export) and the ``current`` clock-track QC: the samples whose exported array differs from the
    current one are returned to be dropped, and the export counts as current. Without one: the export is flagged stale
    when it is older than the clock-track QC, when the clock scores are older than the export, when the scores and
    the export hold different sample ids, or when a passing array of the current QC is missing from the export.
    Returns ``(drop, note, stale)``."""
    if manifest is not None and current is not None:
        j = array_ids(current).join(
            manifest.select(pl.col("sample").cast(pl.String), "array"), on="sample", how="inner", suffix="_export"
        )
        bad = j.filter(~(pl.col("array") == pl.col("array_export")).fill_null(False))["sample"].to_list()
        unlisted = clock_samples - set(manifest["sample"].cast(pl.String).to_list())
        drop = set(bad) | (unlisted & passing)
        note = f"manifest: {len(bad)} arrays differ, {len(unlisted & passing)} scored samples not in the manifest"
        return drop, note, False
    flags = []
    if export_mtime is not None and export_mtime < qc_mtime:
        flags.append("export older than the clock-track QC")
    if export_mtime is not None and clocks_mtime < export_mtime:
        flags.append("clock scores older than the export")
    if export_mtime is None and clocks_mtime < qc_mtime:
        flags.append("clock scores older than the clock-track QC")
    if export_samples is not None and export_samples != clock_samples:
        flags.append(f"scores and export differ in {len(export_samples ^ clock_samples)} sample ids")
    absent = passing - (export_samples if export_samples is not None else clock_samples)
    if absent:
        flags.append(f"{len(absent)} passing arrays absent from the export")
    if flags:
        return set(), "stale: " + "; ".join(flags), True
    return set(), "no manifest; times and ids agree", False


def primary_groups(groups: pl.DataFrame) -> pl.DataFrame:
    """One group per accession: among its non-sensitivity groups the one with the most eligible donors (counted
    before any clock is scored); a pooled group wins a tie, then series and group name. ``groups`` carries
    ``accession``, ``series``, ``group``, ``donors`` and ``sensitivity``."""
    if "eligible" in groups.columns:
        groups = groups.filter(pl.col("eligible"))
    return (
        groups.filter(~pl.col("sensitivity"))
        .select("accession", "series", "group", "donors")
        .unique()
        .with_columns(pooled=pl.col("group").str.contains(POOLED))
        .sort(["accession", "donors", "pooled", "series", "group"], descending=[False, True, True, False, False])
        .unique("accession", keep="first", maintain_order=True)
        .select("series", "group")
    )


def seeded_panel(
    pool: np.ndarray, rng: np.random.Generator, complete: Callable[[list[str]], set[str]], size: int = PANEL
) -> list[str]:
    """The first ``size`` CpGs of ``pool``, in a seeded order, that ``complete`` accepts (measured in every array),
    read in batches until ``size`` are found or the pool runs out."""
    order = pool[rng.permutation(pool.size)]
    picked: list[str] = []
    start = 0
    while len(picked) < size and start < order.size:
        stop = min(order.size, start + max(int(1.5 * (size - len(picked))), 1000))
        batch = order[start:stop].tolist()
        ok = complete(batch)
        picked += [c for c in batch if c in ok]
        start = stop
    return picked[:size]


# ----------------------------------------------------------------------------- clock CpGs kept out of the panel
def biolearn_package(root: Path) -> Path:
    """The ``biolearn`` package at the commit ``scripts/clocks_biolearn.py`` pins, from uv's git cache (or
    ``ENTROPY_AGE_BIOLEARN_SRC``, a checkout's ``biolearn`` folder)."""
    if os.environ.get("ENTROPY_AGE_BIOLEARN_SRC"):
        return Path(os.environ["ENTROPY_AGE_BIOLEARN_SRC"])
    pin = re.search(r"biolearn@([0-9a-f]{7,40})", (root / BIOLEARN_SCRIPT).read_text())
    if pin is None:
        raise ValueError(f"{BIOLEARN_SCRIPT}: no pinned biolearn commit")
    caches = [Path(p) for p in (os.environ.get("UV_CACHE_DIR"), Path.home() / ".cache" / "uv") if p]
    for cache in caches:
        hits = sorted((cache / "git-v1" / "checkouts").glob(f"*/{pin.group(1)[:7]}*/biolearn"))
        if hits:
            return hits[0]
    raise FileNotFoundError(
        f"biolearn {pin.group(1)[:7]} is not in uv's git cache; `uv run --script {BIOLEARN_SCRIPT}` fetches it"
    )


def _function_data_files(pkg: Path, name: str) -> list[str]:
    """The ``get_data_file("...")`` literals that biolearn's function ``name`` reads, in its own body and in the bodies
    of the functions of its module it calls (followed transitively)."""
    import ast

    for src in sorted(pkg.glob("*.py")):
        tree = ast.parse(src.read_text())
        defs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
        if name not in defs:
            continue
        files: list[str] = []
        seen, todo = set(), [name]
        while todo:
            fn = todo.pop()
            if fn in seen:
                continue
            seen.add(fn)
            for node in ast.walk(defs[fn]):
                if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
                    continue
                if node.func.id == "get_data_file" and node.args and isinstance(node.args[0], ast.Constant):
                    files.append(str(node.args[0].value))
                elif node.func.id in defs:
                    todo.append(node.func.id)
        return sorted(set(files))
    raise ValueError(f"clock_drift_controls: biolearn function {name} not found in {pkg}")


def gallery_inputs(pkg: Path) -> dict[str, list[str]]:
    """Clock name to every data file its prediction reads (review 5, finding 4): the coefficient ``file``, a
    ``center`` or ``rotation`` file (PC clocks), and the data files of its ``preprocess`` function (DunedinPACE's
    quantile normalisation reads a 20,000-probe gold-standard list beside its coefficients). Imputation fills only the
    clock's own CpGs and adds no input."""
    model_py = (pkg / "model.py").read_text()
    keys = list(re.finditer(r'^    "([^"]+)": \{', model_py, re.M))
    out = {}
    for m, nxt in zip(keys, [*keys[1:], None], strict=True):
        body = model_py[m.end() : nxt.start() if nxt else len(model_py)]
        files = [f.group(2) for f in re.finditer(r'"(file|center|rotation)":\s*"([^"]+)"', body)]
        named = re.findall(r'"preprocess":\s*([A-Za-z_]\w*)\s*[,}\n]', body)
        for fn in named:
            files += _function_data_files(pkg, fn)
        if body.count('"preprocess"') > len(named):  # a lambda or factory: its inputs cannot be read here
            files.append(UNREAD_PREPROCESS)
        if files:
            out[m.group(1)] = sorted(set(files))
    return out


def gallery_cpgs(root: Path, scored: set[str]) -> set[str]:
    """Every CpG id in every data file a scored gallery clock's prediction reads (``gallery_inputs``: coefficients and
    preprocessing probe lists; any column, ids starting ``cg`` or ``ch``)."""
    pkg = biolearn_package(root)
    inputs = gallery_inputs(pkg)
    missing = sorted(scored - set(inputs))
    if missing:
        raise ValueError(f"clock_drift_controls: scored clocks not in biolearn's gallery: {missing}")
    out: set[str] = set()
    for f in sorted({f for c in scored for f in inputs[c]}):
        if f == UNREAD_PREPROCESS:
            raise ValueError("clock_drift_controls: a scored clock's preprocessing inputs could not be read")
        path = pkg / "data" / f
        if f.startswith(("http://", "https://")) or not path.exists():
            raise FileNotFoundError(f"clock_drift_controls: no local prediction input {f} in {pkg / 'data'}")
        frame = pl.read_csv(path, infer_schema_length=0)
        for col in frame.columns:
            out |= {v.strip('"') for v in frame[col].drop_nulls().to_list() if CPG_ID.match(v.strip('"'))}
    return out


def panel_overlap(panel: pl.DataFrame, exclude: set[str]) -> int:
    """Panel probe-group entries (``panel_probes`` rows) whose CpG is a prediction input or otherwise excluded."""
    return int(panel.filter(pl.col("cpg").is_in(list(exclude))).height) if panel.height else 0


def scored_clocks(root: Path) -> set[str]:
    """The biolearn clocks named in any series' clock scores."""
    names: set[str] = set()
    for f in sorted((root / CLOCK_SCORES).glob("*_clocks.csv")):
        try:
            frame = pl.read_csv(f, columns=["clock"], infer_schema_length=0)
        except pl.exceptions.NoDataError:  # a series whose export scored nothing
            continue
        names |= set(frame["clock"].to_list())
    return names


# ----------------------------------------------------------------------------- run
def _panel(path: Path, samples: list[str], exclude: set[str], rng: np.random.Generator) -> tuple[np.ndarray, list[str]]:
    """People x CpGs betas of the seeded panel (``seeded_panel``) outside ``exclude``, and its probe ids in order."""
    ids = pl.scan_parquet(path).select("cpg").collect()["cpg"]
    pool = np.array(sorted(set(ids.to_list()) - exclude))

    def complete(batch: list[str]) -> set[str]:
        f = pl.scan_parquet(path).filter(pl.col("cpg").is_in(batch)).select("cpg", *samples).collect()
        f = f.drop_nulls().filter(pl.all_horizontal(pl.col(samples).is_not_nan()))
        return set(f["cpg"].to_list())

    picked = seeded_panel(pool, rng, complete)
    frame = pl.scan_parquet(path).filter(pl.col("cpg").is_in(picked)).select("cpg", *samples).collect()
    rank = {c: i for i, c in enumerate(picked)}
    frame = frame.with_columns(r=pl.col("cpg").replace_strict(rank, return_dtype=pl.Int64)).sort("r")
    return frame.select(samples).to_numpy().T.astype(float), picked


def _held(s: pl.DataFrame) -> np.ndarray:
    parts = [dummies(s[c].cast(pl.String).to_numpy()) for c in HELD if c in s.columns and s[c].null_count() == 0]
    parts = [p for p in parts if p.shape[1]]
    return np.hstack(parts) if parts else np.zeros((s.height, 0))


def match_strata(x: np.ndarray, island: np.ndarray, chemistry: np.ndarray) -> np.ndarray:
    """Island x decile of the mean beta x decile of the SD of beta across the arrays x probe chemistry, per CpG."""
    mean = baseline_strata(x.mean(axis=0), 10)
    sd = baseline_strata(x.std(axis=0), 10)
    return np.array([f"{int(i)}|m{m}|s{v}|{c}" for i, m, v, c in zip(island, mean, sd, chemistry, strict=True)])


@dataclass
class SeriesInput:
    """One series' inputs that do not depend on clock output: the QC of both rebuilds and the paths."""

    gi: int
    name: str
    spec: dict
    gse: str
    qc: pl.DataFrame
    clock_qc: pl.DataFrame
    clock_qc_path: Path
    mismatched: list[str]
    betas_path: Path
    clocks_path: Path
    unmasked: Path
    dropped: set[str]


@dataclass
class GroupEntry:
    """One group of the registry: its arrays (``s``, ages present) and whether it is eligible (enough arrays)."""

    li: int
    label: str
    s: pl.DataFrame
    series: SeriesInput

    @property
    def samples(self) -> list[str]:
        return self.s["sample"].to_list()


def group_registry(ctx: Context, log: str = "clock_drift_controls") -> tuple[list[GroupEntry], pl.DataFrame]:
    """Every group, built from arrays that pass the clock track's QC and whose two rebuilds read the same physical
    array (``matched_arrays``) and that the drift betas hold, before any clock prediction is read (review 5, finding
    5). Returns the entries and the registry table (``accession``, ``series``, ``group``, ``donors``, ``eligible``:
    at least ``MIN_ARRAYS`` arrays with an age, ``sensitivity``, ``brain``, ``array_mismatch_dropped``)."""
    entries, rows = [], []
    for gi, (name, spec) in enumerate(specs().items()):
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        if not (betas_path.exists() and reprocessed.available(ctx.paths, gse)):
            continue
        clock_qc_path = reprocessed.series_dir(ctx.paths, gse, CLOCK_TRACK) / f"{gse}.qc.csv"
        if not clock_qc_path.exists():
            print(f"{log}: {name}: skipped, no clock-track QC {clock_qc_path}", file=sys.stderr)
            continue
        qc = reprocessed.read_qc(ctx.paths, gse)  # the rebuild the drift CpGs and the panel come from
        clock_qc = reprocessed.read_qc(ctx.paths, gse, sesame_dir=CLOCK_TRACK).filter(pl.col("passed"))
        same_array, mismatched = matched_arrays(qc, clock_qc)
        if mismatched:
            print(
                f"{log}: {name}: {len(mismatched)} samples dropped, their clock export read another "
                f"physical array: {'; '.join(mismatched[:5])}",
                file=sys.stderr,
            )
        measured = set(pl.read_parquet_schema(betas_path).names()) - {"cpg"}
        chosen = selection(ctx.paths, name, spec, measured & same_array)
        if chosen is None:
            continue
        groups, dropped = chosen
        inp = SeriesInput(
            gi,
            name,
            spec,
            gse,
            qc,
            clock_qc,
            clock_qc_path,
            mismatched,
            betas_path,
            ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv",
            ctx.paths.root / CLOCK_SCORES / f"{gse}_cpgs.parquet",
            set(dropped),
        )
        for li, (label, s) in enumerate(groups.items()):
            s = s.drop_nulls("age")
            entries.append(GroupEntry(li, label, s, inp))
            rows.append(
                {
                    "accession": gse,
                    "series": name,
                    "group": label,
                    "donors": s.height,
                    "eligible": s.height >= MIN_ARRAYS,
                    "sensitivity": bool(SENSITIVITY.search(name)),
                    "brain": bool(BRAIN.search(f"{name} {label}")),
                    "array_mismatch_dropped": len(mismatched),
                }
            )
    return entries, pl.DataFrame(rows)


@dataclass
class SeriesScores:
    """A series' clock scores (wide, one row per sample, ``None`` when it has none), the unmasked export's mitotic
    CpGs, and the export check (``export_check``)."""

    base: pl.DataFrame | None
    mb: pl.DataFrame | None
    note: str
    stale: bool
    dropped: set[str]


def series_scores(inp: SeriesInput, clock_cpgs: set[str]) -> SeriesScores:
    if not inp.clocks_path.exists():
        return SeriesScores(None, None, "no clock scores", False, set())
    long = pl.read_csv(inp.clocks_path, schema_overrides={"sample": pl.String}).select("sample", "clock", "prediction")
    export_samples, export_mtime, mb = None, None, None
    if inp.unmasked.exists():
        export_samples = set(pl.read_parquet_schema(inp.unmasked).names()) - {"cpg"}
        export_mtime = inp.unmasked.stat().st_mtime
        mb = pl.scan_parquet(inp.unmasked).filter(pl.col("cpg").is_in(list(clock_cpgs))).collect()
    man = manifest_path(inp.unmasked)
    drop, note, stale = export_check(
        qc_mtime=inp.clock_qc_path.stat().st_mtime,
        export_mtime=export_mtime,
        clocks_mtime=inp.clocks_path.stat().st_mtime,
        export_samples=export_samples,
        clock_samples=set(long["sample"].to_list()),
        passing=set(inp.clock_qc["sample"].cast(pl.String).to_list()),
        manifest=pl.read_csv(man, schema_overrides={"sample": pl.String, "array": pl.String}) if man.exists() else None,
        current=inp.clock_qc,
    )
    if stale or drop:
        print(f"clock scores {inp.name}: {note}", file=sys.stderr)
    long = long.filter(~pl.col("sample").is_in(list(drop)))
    if mb is not None and drop:
        mb = mb.drop([c for c in drop if c in mb.columns])
    base = long.cast({"prediction": pl.Float64}).pivot(on="clock", index="sample", values="prediction")
    return SeriesScores(base, mb, note, stale, drop)


def group_clock_scores(samples: list[str], sc: SeriesScores, mitotic: dict) -> pl.DataFrame:
    """One row per array of ``samples`` (in order) with every clock's score, null where it has none."""
    frame = pl.DataFrame({"sample": samples})
    if sc.base is not None:
        frame = frame.join(sc.base, on="sample", how="left", validate="1:1")
    if sc.mb is not None:
        in_mb = [v for v in samples if v in sc.mb.columns]
        if in_mb:
            extra, _ = mitotic_clocks.score(sc.mb, in_mb, mitotic)
            if extra.height:
                frame = frame.join(
                    extra.pivot(on="clock", index="sample", values="prediction"), on="sample", how="left"
                )
    if frame["sample"].to_list() != samples:
        raise ValueError("clock scores out of the arrays' order")
    return frame


def excluded_cpgs(ctx: Context, log: str = "clock_drift_controls") -> tuple[set[str], set[str]]:
    """The mitotic clocks' CpGs, and every other prediction input: the local clock files and every data file a scored
    gallery clock reads (``gallery_cpgs``)."""
    mitotic = mitotic_clocks.load(ctx.paths.root)
    if mitotic is None:
        raise FileNotFoundError(f"{log}: no mitotic coefficients in {mitotic_clocks.DIR}")
    gallery = scored_clocks(ctx.paths.root)
    other = set(clock_table(ctx.paths)["cpg"].to_list()) | gallery_cpgs(ctx.paths.root, gallery)
    print(
        f"{log}: leaving out {len(other)} prediction-input CpGs of {len(gallery)} gallery clocks (coefficients and "
        f"preprocessing probe lists) and the local clock files, and {len(mitotic_clocks.cpgs(mitotic))} mitotic CpGs",
        file=sys.stderr,
        flush=True,
    )
    return mitotic_clocks.cpgs(mitotic), other


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # clock_drift's stream: the same drift CpGs
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("cpg").is_in(features))
    chemistry = dict(probes.select("cpg", "probe_type").iter_rows())
    links = pl.read_parquet(ann.links_450k)
    set_genes = set(membership["gene"].to_list())
    mitotic = mitotic_clocks.load(ctx.paths.root)
    clock_cpgs, other_clocks = excluded_cpgs(ctx)
    out = ctx.out_dir("clock_drift_controls")
    old_panel = out / "panel_probes.parquet"
    if old_panel.exists():  # how far the previous run's saved panel fell inside the full exclusion
        print(
            f"clock_drift_controls: the saved panel holds {panel_overlap(pl.read_parquet(old_panel), other_clocks)} "
            "probe-group entries that are prediction inputs",
            file=sys.stderr,
        )
    entries, registry = group_registry(ctx)
    n_placebo = N_PLACEBO[str(ctx.scale.name)]
    rows, draws, strata_rows, panel_rows = [], [], [], []
    reasons: dict[tuple[str, str, str], str] = {}
    cache: dict[str, SeriesScores] = {}
    betas_cache: tuple[str, pl.DataFrame] | None = None
    for e in entries:
        inp, label, s, samples = e.series, e.label, e.s, e.samples
        name, gse = inp.name, inp.gse
        if s.height < MIN_ARRAYS:
            reasons[(name, label, "*")] = f"{s.height} arrays with an age < {MIN_ARRAYS}"
            continue
        if gse not in cache:
            cache = {gse: series_scores(inp, clock_cpgs)}  # one series at a time in memory
        sc = cache[gse]
        if sc.base is None and sc.mb is None:
            reasons[(name, label, "*")] = sc.note
            continue
        if betas_cache is None or betas_cache[0] != gse:
            betas = pl.read_parquet(inp.betas_path).filter(
                pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(list(clock_cpgs | inp.dropped))
            )
            betas_cache = (gse, betas)
        betas = betas_cache[1]
        exclude = set(features) | clock_cpgs | other_clocks | inp.dropped
        sub = betas.select("cpg", *samples)
        x = sub.select(samples).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
        x = x[:, ok]
        sets = _sets(kept, membership)
        if GROUP not in sets:
            reasons[(name, label, "*")] = "no housekeeping set among the measured CpGs"
            continue
        island = _islands(kept, probes)
        member = np.any([m for m, _ in sets.values()], axis=0)
        hk, hk_genes = sets[GROUP]
        match = match_strata(x, island, np.array([chemistry.get(c, "NA") for c in kept]))
        rng = ctx.rng(130_000 + 100 * inp.gi + e.li)
        pool = placebo_pool(kept, member, links, set_genes)
        fills, pairs, filled, gb = [], {}, [], []
        for k in range(n_placebo):
            sub_k, mask, background, fill = placebo_pair(match, hk, pool, rng)
            pairs[f"placebo{k}"] = (sub_k, mask, background)
            fills.append(fill)
            filled.append(mask)
        gb_rng = ctx.rng(140_000 + 100 * inp.gi + e.li)
        for k in range(n_placebo):
            sub_k, mask, background, fill = gene_balanced_pair(match, hk_genes, pool, gb_rng, GENE_CAP)
            pairs[f"genebal{k}"] = (sub_k, mask, background)
            gb.append(
                {
                    "fill": fill,
                    "cpgs": int(mask.sum()),
                    "hk_genes": int(per_gene_counts(sub_k, hk_genes).size),
                    "placebo_genes": int(per_gene_counts(mask, pool).size),
                    "hk_max_per_gene": int(per_gene_counts(sub_k, hk_genes).max(initial=0)),
                    "placebo_max_per_gene": int(per_gene_counts(mask, pool).max(initial=0)),
                }
            )
        levels, codes = np.unique(match, return_inverse=True)
        quota = np.bincount(codes[hk], minlength=levels.size)
        got = np.mean([np.bincount(codes[m], minlength=levels.size) for m in filled], axis=0)
        strata_rows += [
            {"series": name, "group": label, "stratum": str(st), "quota": int(q), "filled_mean": float(g)}
            for st, q, g in zip(levels, quota, got, strict=True)
            if q
        ]
        panel, picked = _panel(inp.betas_path, samples, exclude, rng)
        panel_rows += [{"series": name, "group": label, "rank": r, "cpg": c} for r, c in enumerate(picked)]
        print(
            f"clock_drift_controls: {name} / {label}: {len(samples)} arrays, {len(kept)} CpGs, panel "
            f"{panel.shape[1]}, placebo fill {np.mean(fills):.3f}, gene-balanced fill "
            f"{np.mean([g['fill'] for g in gb]):.3f} ({np.mean([g['cpgs'] for g in gb]):.0f} genes)",
            file=sys.stderr,
            flush=True,
        )
        per = loo_drift(x, island, {GROUP: sets[GROUP]}, pairs=pairs)
        b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
        a = per[f"P_{GROUP}"] + b
        ind = panel_indices(panel)
        age = s["age"].to_numpy().astype(float)
        held = _held(s)
        quality, qdesc, _ = quality_terms(samples, inp.qc, base_columns=6 + held.shape[1])
        clk = group_clock_scores(samples, sc, mitotic)
        for clock in set(MITOTIC) - set(clk.columns):
            reasons[(name, label, clock)] = "not scored here (coverage rule, or absent from the clock scores)"
        for clock in [c for c in clk.columns if c != "sample"]:
            v = clk[clock].cast(pl.Float64).to_numpy()
            use = np.isfinite(v)
            if use.sum() < MIN_ARRAYS:
                reasons[(name, label, clock)] = f"{int(use.sum())} arrays scored < {MIN_ARRAYS}"
                continue
            if np.std(v[use]) == 0:
                reasons[(name, label, clock)] = "constant clock"
                continue
            design = adjust_design(age, quality, held, ind["level"], ind["compression"])[use]
            stats = clock_stats(v[use], age[use], design, f[use], a[use], b[use], {k: w[use] for k, w in ind.items()})
            rows.append(
                {
                    "accession": gse,
                    "series": name,
                    "group": label,
                    "sensitivity": bool(SENSITIVITY.search(name)),
                    "brain": bool(BRAIN.search(f"{name} {label}")),
                    "donors": len(samples),
                    "n": int(use.sum()),
                    "cpgs": len(kept),
                    "panel_cpgs": int(panel.shape[1]),
                    "array_mismatch_dropped": len(inp.mismatched),
                    "export_check": sc.note,
                    "export_stale": sc.stale,
                    "export_dropped": len(sc.dropped),
                    "covariates": qdesc + (", batch/center" if held.shape[1] else ""),
                    "clock": clock,
                    "mitotic": clock in MITOTIC,
                    **stats,
                    "placebo_fill": float(np.mean(fills)),
                    "placebo_fill_gb": float(np.mean([g["fill"] for g in gb])),
                    "genes_gb": float(np.mean([g["cpgs"] for g in gb])),
                }
            )
            acc, _ = age_residual_checked(v[use], age[use])
            for scheme, prefix, meta in (("cpg", "placebo", None), ("gene", "genebal", gb)):
                for k in range(n_placebo):
                    g = f"{prefix}{k}"
                    fh, fp = per[f"F_{g}_a"][use], per[f"F_{g}_b"][use]
                    rh, nh = partial_r_checked(v[use], fh, design, ("clock", "matched housekeeping F"))
                    rp, npl = partial_r_checked(v[use], fp, design, ("clock", "placebo F"))
                    draws.append(
                        {
                            "series": name,
                            "group": label,
                            "clock": clock,
                            "scheme": scheme,
                            "draw": k,
                            "fill": fills[k] if meta is None else meta[k]["fill"],
                            "hk_genes": None if meta is None else meta[k]["hk_genes"],
                            "placebo_genes": None if meta is None else meta[k]["placebo_genes"],
                            "hk_max_per_gene": None if meta is None else meta[k]["hk_max_per_gene"],
                            "placebo_max_per_gene": None if meta is None else meta[k]["placebo_max_per_gene"],
                            "balance_tvd": float(per[f"tvd_{g}"].mean()),
                            "balance_unusable": float(per[f"unusable_{g}"].mean()),
                            "hk_r_raw": corr(acc, fh),
                            "hk_r_adj": rh,
                            "r_raw": corr(acc, fp),
                            "r_adj": rp,
                            "adj_note": "; ".join(n for n in (nh, npl) if n),
                        }
                    )
    if not rows:
        print("clock_drift_controls: nothing to score", file=sys.stderr)
        return []
    placebo = pl.DataFrame(draws, infer_schema_length=None)
    band = placebo_band(placebo.filter(pl.col("scheme") == "cpg")).join(
        placebo_band(placebo.filter(pl.col("scheme") == "gene"), suffix="_gb"),
        on=["series", "group", "clock"],
        how="left",
    )
    table = pl.DataFrame(rows).join(band, on=["series", "group", "clock"], how="left")
    table = flag_primary(table, registry)
    registry = flag_primary(registry, registry)
    clocks = sorted(set(table["clock"].to_list()) | set(MITOTIC))
    missing = missing_estimates(registry, table, clocks, reasons)
    panel_table = pl.DataFrame(panel_rows)
    write_table(table, out / "correlations", csv_max_rows=10_000)
    write_table(placebo, out / "placebo", csv_max_rows=0)
    write_table(pl.DataFrame(strata_rows), out / "placebo_strata", csv_max_rows=10_000)
    write_table(panel_table, out / "panel_probes", csv_max_rows=0)
    write_table(registry, out / "groups", csv_max_rows=10_000)
    write_table(missing, out / "missing", csv_max_rows=10_000)
    write_table(summarise(table), out / "summary", csv_max_rows=10_000)
    saved = pl.read_parquet(out / "panel_probes.parquet")
    overlap = panel_overlap(saved, other_clocks | clock_cpgs | set(features))
    print(f"clock_drift_controls: saved panel entries inside the exclusion: {overlap}", file=sys.stderr)
    if overlap:
        raise ValueError(f"clock_drift_controls: {overlap} saved panel entries are prediction inputs or drift CpGs")
    return [draw(table.filter(pl.col("primary") & ~pl.col("brain")), ctx.figure_path("clock_drift_controls"))]


def flag_primary(table: pl.DataFrame, registry: pl.DataFrame | None = None) -> pl.DataFrame:
    """``table`` with a ``primary`` flag. The primary groups come from ``registry`` (the groups and their eligible
    donors, built before any prediction is read; review 5, finding 5), so a group that emits no clock row keeps its
    place and a smaller group never stands in for it. Without a registry, from ``table`` itself."""
    primary = primary_groups(table if registry is None else registry).with_columns(primary=pl.lit(True))
    return table.join(primary, on=["series", "group"], how="left").with_columns(pl.col("primary").fill_null(False))


def missing_estimates(
    registry: pl.DataFrame, table: pl.DataFrame, clocks: list[str], reasons: dict[tuple[str, str, str], str]
) -> pl.DataFrame:
    """Per primary group of ``registry`` and clock of ``clocks`` without a row in ``table``: the reason (from
    ``reasons``, keyed by series, group and clock, or by series, group and ``*`` for the whole group)."""
    have = set(table.select("series", "group", "clock").iter_rows()) if table.height else set()
    rows = []
    prim = registry.join(primary_groups(registry), on=["series", "group"], how="semi")
    for acc, series, group, brain in prim.select("accession", "series", "group", "brain").iter_rows():
        for clock in clocks:
            if (series, group, clock) in have:
                continue
            why = reasons.get((series, group, clock)) or reasons.get((series, group, "*")) or "no estimate"
            rows.append(
                {"accession": acc, "series": series, "group": group, "brain": brain, "clock": clock, "reason": why}
            )
    schema = {"accession": pl.String, "series": pl.String, "group": pl.String, "brain": pl.Boolean}
    return pl.DataFrame(rows, schema={**schema, "clock": pl.String, "reason": pl.String})


def placebo_band(placebo: pl.DataFrame, suffix: str = "") -> pl.DataFrame:
    """Per group and clock: the housekeeping r on the matched subsets (median over draws) and the placebo
    distribution of |r| (median, 95th percentile), raw and adjusted; the share of draws whose placebo r reaches its
    paired housekeeping r in that r's direction (draws with both defined); the housekeeping margin on the common
    support (median over draws of ``sign(hk r) * (hk r - placebo r)``, positive when the matched housekeeping subset
    tracks the clock more strongly than its paired placebo); whether the matched housekeeping |r| exceeds the placebo
    95th percentile; mean fill and balance. ``suffix`` is appended to every statistic's name."""

    def both(hk: str, p: str) -> pl.Expr:
        return pl.col(hk).is_not_nan() & pl.col(p).is_not_nan()

    def beyond(hk: str, p: str) -> pl.Expr:
        hit = (pl.col(p) * pl.col(hk).sign()) >= pl.col(hk).abs()
        return hit.filter(both(hk, p)).mean()

    def margin(hk: str, p: str) -> pl.Expr:
        return (pl.col(hk).sign() * (pl.col(hk) - pl.col(p))).filter(both(hk, p)).median()

    band = placebo.group_by("series", "group", "clock", maintain_order=True).agg(
        hk_matched_raw=pl.col("hk_r_raw").drop_nans().median(),
        hk_matched_adj=pl.col("hk_r_adj").drop_nans().median(),
        placebo_abs_median_raw=pl.col("r_raw").drop_nans().abs().median(),
        placebo_abs_p95_raw=pl.col("r_raw").drop_nans().abs().quantile(0.95, "linear"),
        placebo_abs_median_adj=pl.col("r_adj").drop_nans().abs().median(),
        placebo_abs_p95_adj=pl.col("r_adj").drop_nans().abs().quantile(0.95, "linear"),
        placebo_mean_raw=pl.col("r_raw").drop_nans().mean(),
        placebo_beyond_raw=beyond("hk_r_raw", "r_raw"),
        placebo_beyond_adj=beyond("hk_r_adj", "r_adj"),
        margin_raw=margin("hk_r_raw", "r_raw"),
        margin_adj=margin("hk_r_adj", "r_adj"),
        placebo_balance_tvd=pl.col("balance_tvd").mean(),
        placebo_balance_unusable=pl.col("balance_unusable").mean(),
        placebos=pl.len(),
        placebos_undefined_adj=(pl.col("adj_note") != "").sum(),
    )
    band = band.with_columns(
        exceeds_p95_raw=pl.col("hk_matched_raw").abs() > pl.col("placebo_abs_p95_raw"),
        exceeds_p95_adj=pl.col("hk_matched_adj").abs() > pl.col("placebo_abs_p95_adj"),
    )
    keys = ("series", "group", "clock")
    return band.rename({c: f"{c}{suffix}" for c in band.columns if c not in keys}) if suffix else band


SUMMARY_COLS = (
    "r_F_raw",
    "r_F_adj",
    "r_A_given_B",
    "hk_matched_raw",
    "hk_matched_adj",
    "r_level",
    "r_compression",
    "r_F_level",
    "r_F_compression",
    "placebo_abs_median_raw",
    "placebo_abs_p95_raw",
    "placebo_abs_median_adj",
    "placebo_abs_p95_adj",
    "placebo_beyond_raw",
    "placebo_beyond_adj",
    "placebo_fill",
    "placebo_balance_tvd",
    "margin_raw",
    "margin_adj",
    "hk_matched_adj_gb",
    "placebo_abs_p95_adj_gb",
    "margin_raw_gb",
    "margin_adj_gb",
    "placebo_fill_gb",
    "genes_gb",
)


def summarise(table: pl.DataFrame) -> pl.DataFrame:
    """Per mitotic clock: medians over primary non-brain groups (``level`` = ``primary groups``, with the counts of
    groups whose matched housekeeping |r| exceeds the placebo 95th percentile and of undefined adjusted r) and
    medians over accessions of each accession's median over its non-sensitivity non-brain groups (``level`` =
    ``accession medians``)."""
    cols = [c for c in SUMMARY_COLS if c in table.columns]
    t = table.filter(pl.col("mitotic") & ~pl.col("brain"))
    prim = (
        t.filter(pl.col("primary"))
        .group_by("clock", maintain_order=True)
        .agg(
            pl.len().alias("groups"),
            pl.col("r_F_adj").is_nan().sum().alias("undefined_adj"),
            pl.col("exceeds_p95_raw").sum().alias("exceeding_raw"),
            pl.col("exceeds_p95_adj").sum().alias("exceeding_adj"),
            *(
                [pl.col("exceeds_p95_adj_gb").sum().alias("exceeding_adj_gb")]
                if "exceeds_p95_adj_gb" in table.columns
                else []
            ),
            *(
                [(pl.col("margin_adj_gb") > 0).sum().alias("margin_gb_positive")]
                if "margin_adj_gb" in table.columns
                else []
            ),
            *(pl.col(c).drop_nans().median().alias(c) for c in cols),
        )
        .with_columns(level=pl.lit("primary groups"))
    )
    acc = (
        t.filter(~pl.col("sensitivity"))
        .group_by("clock", "accession", maintain_order=True)
        .agg(*(pl.col(c).drop_nans().median().alias(c) for c in cols))
        .group_by("clock", maintain_order=True)
        .agg(pl.len().alias("groups"), *(pl.col(c).drop_nans().median().alias(c) for c in cols))
        .with_columns(level=pl.lit("accession medians"))
    )
    out = pl.concat([prim, acc], how="diagonal_relaxed")
    return out.select("level", *[c for c in out.columns if c != "level"]).sort(
        ["level", "clock"], descending=[True, False]
    )


def draw(table: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    d = table.filter(pl.col("mitotic"))
    clocks = [c for c in MITOTIC if c in d["clock"].to_list()]
    fig, axes = plots.plt.subplots(1, 2, figsize=(12, 4.2), sharey=False)
    ax = axes[0]
    for i, c in enumerate(clocks):
        t = d.filter(pl.col("clock") == c)
        for off, col, color in ((-0.15, "r_F_raw", plots.MUTED), (0.15, "r_F_adj", plots.PALETTE[0])):
            v = t[col].drop_nans().drop_nulls().to_numpy()
            jitter = np.random.default_rng(i).normal(0, 0.04, v.size)
            ax.scatter(np.full(v.size, i + off) + jitter, v, s=10, alpha=0.6, color=color)
            if v.size:
                ax.plot([i + off - 0.12, i + off + 0.12], [np.median(v)] * 2, color=plots.INK, lw=2)
        for raw, adj in t.select("r_F_raw", "r_F_adj").iter_rows():
            ax.plot([i - 0.15, i + 0.15], [raw, adj], color=plots.MUTED, lw=0.4, alpha=0.5)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_xticks(range(len(clocks)), clocks, fontsize=7)
    ax.set_ylabel("r, clock acceleration with F\n(grey: age only; blue: age + quality + level + compression)")
    ax.set_title("a  raw and adjusted, primary non-brain groups", loc="left", fontsize=9)
    ax = axes[1]
    for i, c in enumerate(clocks):
        t = d.filter(pl.col("clock") == c).sort("hk_matched_adj")
        xs = i + np.linspace(-0.3, 0.3, max(t.height, 1))[: t.height]
        ax.vlines(xs, -t["placebo_abs_p95_adj"], t["placebo_abs_p95_adj"], color=plots.MUTED, lw=2, alpha=0.5)
        ax.scatter(xs, t["hk_matched_adj"], s=12, color=plots.PALETTE[1], zorder=3)
    ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
    ax.set_xticks(range(len(clocks)), clocks, fontsize=7)
    ax.set_ylabel("adjusted r with F (dot: matched housekeeping subset;\nbar: ±95th pct of |r| over placebo sets)")
    ax.set_title("b  matched housekeeping subsets against placebo promoters", loc="left", fontsize=9)
    fig.tight_layout()
    return plots.save(fig, path)
