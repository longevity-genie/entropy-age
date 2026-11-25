"""Can a broad array shift make the clock-shield link? Confounder and placebo controls for ``clock_drift``.

``clock_drift`` finds that a mitotic clock's acceleration tracks the per-person shield share F (r about -0.5 outside
the brain; HypoClock, whose value falls with division, about +0.3). F is a ratio: with A the housekeeping promoters'
deviation from the leave-one-out reference and B that of their look-alikes, ``F = (A - B) / B``. A shift of the
whole array (a methylation level change, beta compression toward 0.5, chip, signal or detection quality, cell
mixture) that moves A and B by different amounts moves F, and the same shift can move a clock (compression raises
low-beta promoter clock CpGs and lowers high-beta HypoClock CpGs at once). Codex's second review (item 3) proposed
the check built here; its fourth review (findings 1, 3, 5, 6 and 7) shaped the placebo comparison, the panel, the
primary groups, the guards on partial correlations and the array checks.

Per series group, with the arrays, CpGs and clock scores of ``clock_drift`` (``series.selection``, the same
availability rule, every mitotic and EpiTOC1 CpG kept out of the drift CpGs, sensitivity configs flagged):

0. **Arrays.** The clock track's QC (``sesame_clock``) is required for every series: a series without it is skipped
   with a logged reason, and only arrays passing it enter. For every sample id in both rebuilds, the physical array
   (chip and position from each track's ``qc.csv``, the file prefix where either is missing) must be the same; a
   sample whose two exports read different arrays is dropped, logged and counted (``array_mismatch_dropped``).
1. **Drift measures.** ``loo_scores`` per person: ``P_housekeeping`` (the matched gap), ``B_housekeeping`` and
   ``F_housekeeping``; ``A = P + B`` is the housekeeping CpGs' own stratum-weighted deviation.
2. **A third, disjoint probe panel.** CpGs of the same rebuild drawn in a seeded order, kept when measured in every
   array of the group, until 20,000 are found or the pool runs out. The pool leaves out the drift CpGs, every
   mitotic clock's CpGs, every CpG of every clock in the exported biolearn gallery (each clock named in the clock
   scores, read from its coefficient file at the commit ``scripts/clocks_biolearn.py`` pins) and the CpGs of the
   local clock files (``clock_cpgs.clock_table``). The probe ids per group are saved (``panel_probes``). Per person:
   the panel's mean beta (``level``), ``compression`` = mean ``|beta - 0.5|`` (lower is more compressed toward 0.5)
   and the SD of its betas (reported).
3. **Covariates.** Age; chip dummies, log mean intensity and detection fraction of the rebuild the drift CpGs come
   from (``reprocessed.read_qc`` through ``clock_disagreement.quality_terms``, chip dropped when it would leave under
   10 residual degrees of freedom); ``batch`` or ``center`` dummies where the series description holds them. The
   series groups carry no cell fractions (only EPIC-Italy has EpiDISH fractions, and it is not a series group), so
   none enter; level and compression stand in for a mixture that moves the whole array.
4. **Same adjustment on both sides.** The clock's acceleration (residual on age, as in ``clock_drift``), A, B and F
   are all residualised on age + covariates + level + compression, then correlated. Reported per group and clock:
   ``r_F_raw`` (``clock_drift``'s value), ``r_F_adj``, ``r_A_given_B`` (acceleration with A at a given B, the
   adjustment set plus B), ``r_level``, ``r_compression``, ``r_sd`` (the acceleration against the panel indices),
   ``r_F_level`` and ``r_F_compression`` (F against them), and n. Each partial correlation checks, on that clock's
   own rows, the rank of the column-scaled design and the residual degrees of freedom (at least 10), and calls a
   vector explained by the design when its residual norm is at most 1e-8 of its centred norm; such an estimate is
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
6. **Primary groups.** One per accession (study): the non-sensitivity group with the most eligible donors (arrays
   after ``series.selection``, counted before any clock is scored; a pooled group wins a tie). Every row stays in
   ``correlations`` with a ``primary`` flag. ``summary`` gives, per mitotic clock, medians over primary non-brain
   groups with the count of groups exceeding the placebo band, and medians over accessions of each accession's
   median over its non-sensitivity non-brain groups.

Writes ``clock_drift_controls/correlations``, ``placebo`` (per group, clock and draw), ``placebo_strata``,
``panel_probes``, ``summary`` and a figure. A Slurm job: the leave-one-out scoring runs once per group, with the
placebo pairs inside it.
"""

from __future__ import annotations

import os
import re
import sys
from collections.abc import Callable
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
N_PLACEBO = {"paper": 20, "quick": 5}
MITOTIC = ("EpiTOC1", "epiTOC2", "epiTOC3", "HypoClock", "RepliTali", "epiCMIT")
HELD = ("batch", "center")  # technical columns a series description may hold
BRAIN = re.compile(r"brain|cortex|cerebell|hippocamp|amygdala|neuron|gyrus|striatum|caudate|glia|entorhinal", re.I)
POOLED = r"(?i)\ball\b|every|pooled"  # a group pooling a series' strata wins a tie on donors
MIN_RESID_DF = 10  # as clock_disagreement's chip guard
RESID_TOL = 1e-8  # a residual at most this share of the centred norm counts as explained by the design
BIOLEARN_SCRIPT = Path("scripts/clocks_biolearn.py")
CPG_ID = re.compile(r"c[gh][\d.]")


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
    d = np.column_stack([np.ones(len(x)), given])
    scale = np.linalg.norm(d, axis=0)
    rank = int(np.linalg.matrix_rank(d / np.where(scale > 0, scale, 1.0)))
    df = len(x) - rank
    if df < min_df:
        return float("nan"), f"residual df {df} < {min_df} (n {len(x)}, design rank {rank})"
    res = []
    for v, name in zip((x, y), names, strict=True):
        size = float(np.linalg.norm(v - v.mean()))
        if size == 0:
            return float("nan"), f"{name} constant"
        r = v - d @ np.linalg.lstsq(d, v, rcond=None)[0]
        if np.linalg.norm(r) <= RESID_TOL * size:
            return float("nan"), f"{name} explained by the design"
        res.append(r)
    return corr(res[0], res[1]), ""


def partial_r(x: np.ndarray, y: np.ndarray, given: np.ndarray) -> float:
    return partial_r_checked(x, y, given)[0]


def matched_placebo(
    strata: np.ndarray, real: np.ndarray, pool_genes: list[np.ndarray], rng: np.random.Generator
) -> tuple[np.ndarray, list[np.ndarray], float]:
    """A placebo set with ``real``'s count of CpGs in every stratum, filled gene by gene from ``pool_genes`` (index
    arrays of candidate CpGs, one per gene) in random order; a gene adds those of its CpGs whose stratum still has
    room. Returns (membership mask, the genes' index arrays as used, share of the quota filled)."""
    levels, codes = np.unique(strata, return_inverse=True)
    quota = np.bincount(codes[real], minlength=levels.size)
    want = int(quota.sum())
    taken = np.zeros(strata.size, dtype=bool)
    genes = []
    for g in rng.permutation(len(pool_genes)):
        if quota.sum() == 0:
            break
        use = []
        for j in pool_genes[g]:
            if not taken[j] and not real[j] and quota[codes[j]] > 0:
                quota[codes[j]] -= 1
                taken[j] = True
                use.append(j)
        if use:
            genes.append(np.array(use, dtype=int))
    return taken, genes, float(taken.sum() / want) if want else float("nan")


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


def primary_groups(groups: pl.DataFrame) -> pl.DataFrame:
    """One group per accession: among its non-sensitivity groups the one with the most eligible donors (counted
    before any clock is scored); a pooled group wins a tie, then series and group name. ``groups`` carries
    ``accession``, ``series``, ``group``, ``donors`` and ``sensitivity``."""
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


def gallery_files(model_py: str) -> dict[str, str]:
    """Clock name to coefficient file, from biolearn's ``model_definitions`` source (each top-level key and the first
    ``"file"`` of its entry)."""
    keys = list(re.finditer(r'^    "([^"]+)": \{', model_py, re.M))
    out = {}
    for m, nxt in zip(keys, [*keys[1:], None], strict=True):
        body = model_py[m.end() : nxt.start() if nxt else len(model_py)]
        f = re.search(r'"file":\s*"([^"]+)"', body)
        if f:
            out[m.group(1)] = f.group(1)
    return out


def gallery_cpgs(root: Path, scored: set[str]) -> set[str]:
    """Every CpG id in the coefficient file of each scored gallery clock (any column, ids starting ``cg`` or ``ch``)."""
    pkg = biolearn_package(root)
    files = gallery_files((pkg / "model.py").read_text())
    missing = sorted(scored - set(files))
    if missing:
        raise ValueError(f"clock_drift_controls: scored clocks not in biolearn's gallery: {missing}")
    out: set[str] = set()
    for f in sorted({files[c] for c in scored}):
        path = pkg / "data" / f
        if not path.exists():
            raise FileNotFoundError(f"clock_drift_controls: no local coefficient file {path}")
        frame = pl.read_csv(path, infer_schema_length=0)
        for col in frame.columns:
            out |= {v for v in frame[col].drop_nulls().to_list() if CPG_ID.match(v)}
    return out


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


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # clock_drift's stream: the same drift CpGs
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("cpg").is_in(features))
    chemistry = dict(probes.select("cpg", "probe_type").iter_rows())
    links = pl.read_parquet(ann.links_450k)
    set_genes = set(membership["gene"].to_list())
    mitotic = mitotic_clocks.load(ctx.paths.root)
    if mitotic is None:
        raise FileNotFoundError(f"clock_drift_controls: no mitotic coefficients in {mitotic_clocks.DIR}")
    clock_cpgs = mitotic_clocks.cpgs(mitotic)
    gallery = scored_clocks(ctx.paths.root)
    other_clocks = set(clock_table(ctx.paths)["cpg"].to_list()) | gallery_cpgs(ctx.paths.root, gallery)
    print(
        f"clock_drift_controls: the panel leaves out {len(other_clocks)} CpGs of {len(gallery)} gallery clocks and "
        f"the local clock files, and {len(clock_cpgs)} mitotic CpGs",
        file=sys.stderr,
        flush=True,
    )
    n_placebo = N_PLACEBO[str(ctx.scale.name)]
    rows, draws, strata_rows, panel_rows = [], [], [], []
    for gi, (name, spec) in enumerate(specs().items()):
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv"
        unmasked = ctx.paths.root / CLOCK_SCORES / f"{gse}_cpgs.parquet"
        if not (betas_path.exists() and clocks_path.exists() and reprocessed.available(ctx.paths, gse)):
            continue
        clock_qc_path = reprocessed.series_dir(ctx.paths, gse, CLOCK_TRACK) / f"{gse}.qc.csv"
        if not clock_qc_path.exists():
            print(f"clock_drift_controls: {name}: skipped, no clock-track QC {clock_qc_path}", file=sys.stderr)
            continue
        qc = reprocessed.read_qc(ctx.paths, gse)  # the rebuild the drift CpGs and the panel come from
        clock_qc = reprocessed.read_qc(ctx.paths, gse, sesame_dir=CLOCK_TRACK).filter(pl.col("passed"))
        same_array, mismatched = matched_arrays(qc, clock_qc)
        if mismatched:
            print(
                f"clock_drift_controls: {name}: {len(mismatched)} samples dropped, their clock export read another "
                f"physical array: {'; '.join(mismatched[:5])}",
                file=sys.stderr,
            )
        long = pl.read_csv(clocks_path, schema_overrides={"sample": pl.String}).select("sample", "clock", "prediction")
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(clock_cpgs))
        mb = None
        available = set(betas.columns) & set(long["sample"].to_list()) & same_array
        if unmasked.exists():
            mb = pl.scan_parquet(unmasked).filter(pl.col("cpg").is_in(clock_cpgs)).collect()
            available &= set(mb.columns)
        chosen = selection(ctx.paths, name, spec, available)
        if chosen is None:
            continue
        groups, dropped = chosen
        if dropped:
            betas = betas.filter(~pl.col("cpg").is_in(dropped))
        base = long.cast({"prediction": pl.Float64}).pivot(on="clock", index="sample", values="prediction")
        exclude = set(features) | clock_cpgs | other_clocks | set(dropped)
        for li, (label, s) in enumerate(groups.items()):
            s = s.drop_nulls("age")
            if s.height < MIN_ARRAYS:
                continue
            samples = s["sample"].to_list()
            scores = base
            if mb is not None:
                extra, _ = mitotic_clocks.score(mb, samples, mitotic)
                if extra.height:
                    scores = base.join(
                        extra.pivot(on="clock", index="sample", values="prediction"), on="sample", how="left"
                    )
            sub = betas.select("cpg", *samples)
            x = sub.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            x = x[:, ok]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                continue
            island = _islands(kept, probes)
            member = np.any([m for m, _ in sets.values()], axis=0)
            hk = sets[GROUP][0]
            match = match_strata(x, island, np.array([chemistry.get(c, "NA") for c in kept]))
            rng = ctx.rng(130_000 + 100 * gi + li)
            pool = placebo_pool(kept, member, links, set_genes)
            fills, pairs, filled = [], {}, []
            for k in range(n_placebo):
                sub, mask, background, fill = placebo_pair(match, hk, pool, rng)
                pairs[f"placebo{k}"] = (sub, mask, background)
                fills.append(fill)
                filled.append(mask)
            levels, codes = np.unique(match, return_inverse=True)
            quota = np.bincount(codes[hk], minlength=levels.size)
            got = np.mean([np.bincount(codes[m], minlength=levels.size) for m in filled], axis=0)
            strata_rows += [
                {"series": name, "group": label, "stratum": str(st), "quota": int(q), "filled_mean": float(g)}
                for st, q, g in zip(levels, quota, got, strict=True)
                if q
            ]
            panel, picked = _panel(betas_path, samples, exclude, rng)
            panel_rows += [{"series": name, "group": label, "rank": r, "cpg": c} for r, c in enumerate(picked)]
            print(
                f"clock_drift_controls: {name} / {label}: {len(samples)} arrays, {len(kept)} CpGs, panel "
                f"{panel.shape[1]}, placebo fill {np.mean(fills):.3f}",
                file=sys.stderr,
                flush=True,
            )
            per = loo_drift(x, island, {GROUP: sets[GROUP]}, pairs=pairs)
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            a = per[f"P_{GROUP}"] + b
            ind = panel_indices(panel)
            age = s["age"].to_numpy().astype(float)
            held = _held(s)
            quality, qdesc, _ = quality_terms(samples, qc, base_columns=6 + held.shape[1])
            clk = s.select("sample").join(scores, on="sample", how="left", validate="1:1")
            if clk["sample"].to_list() != samples:
                raise ValueError(f"clock_drift_controls: {name} / {label}: clock scores out of the arrays' order")
            for clock in [c for c in clk.columns if c != "sample"]:
                v = clk[clock].to_numpy().astype(float)
                use = np.isfinite(v)
                if use.sum() < MIN_ARRAYS or np.std(v[use]) == 0:
                    continue
                acc = residual(v[use], age[use])
                design = adjust_design(age, quality, held, ind["level"], ind["compression"])[use]
                r_adj, note_adj = partial_r_checked(acc, f[use], design, ("acceleration", "F"))
                r_ab, note_ab = partial_r_checked(acc, a[use], np.column_stack([design, b[use]]), ("acceleration", "A"))
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
                        "array_mismatch_dropped": len(mismatched),
                        "covariates": qdesc + (", batch/center" if held.shape[1] else ""),
                        "clock": clock,
                        "mitotic": clock in MITOTIC,
                        "r_F_raw": corr(acc, f[use]),
                        "r_B_raw": corr(acc, b[use]),
                        "r_F_adj": r_adj,
                        "r_F_adj_note": note_adj,
                        "r_B_adj": partial_r(acc, b[use], design),
                        "r_A_given_B": r_ab,
                        "r_A_given_B_note": note_ab,
                        "r_level": corr(acc, ind["level"][use]),
                        "r_compression": corr(acc, ind["compression"][use]),
                        "r_sd": corr(acc, ind["sd"][use]),
                        "r_F_level": corr(f[use], ind["level"][use]),
                        "r_F_compression": corr(f[use], ind["compression"][use]),
                        "placebo_fill": float(np.mean(fills)),
                    }
                )
                for k in range(n_placebo):
                    g = f"placebo{k}"
                    fh, fp = per[f"F_{g}_a"][use], per[f"F_{g}_b"][use]
                    rh, nh = partial_r_checked(acc, fh, design, ("acceleration", "matched housekeeping F"))
                    rp, npl = partial_r_checked(acc, fp, design, ("acceleration", "placebo F"))
                    draws.append(
                        {
                            "series": name,
                            "group": label,
                            "clock": clock,
                            "draw": k,
                            "fill": fills[k],
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
    out = ctx.out_dir("clock_drift_controls")
    placebo = pl.DataFrame(draws)
    table = pl.DataFrame(rows).join(placebo_band(placebo), on=["series", "group", "clock"], how="left")
    table = flag_primary(table)
    write_table(table, out / "correlations", csv_max_rows=10_000)
    write_table(placebo, out / "placebo", csv_max_rows=10_000)
    write_table(pl.DataFrame(strata_rows), out / "placebo_strata", csv_max_rows=10_000)
    write_table(pl.DataFrame(panel_rows), out / "panel_probes", csv_max_rows=0)
    write_table(summarise(table), out / "summary", csv_max_rows=10_000)
    return [draw(table.filter(pl.col("primary") & ~pl.col("brain")), ctx.figure_path("clock_drift_controls"))]


def flag_primary(table: pl.DataFrame) -> pl.DataFrame:
    primary = primary_groups(table).with_columns(primary=pl.lit(True))
    return table.join(primary, on=["series", "group"], how="left").with_columns(pl.col("primary").fill_null(False))


def placebo_band(placebo: pl.DataFrame) -> pl.DataFrame:
    """Per group and clock: the housekeeping r on the matched subsets (median over draws) and the placebo
    distribution of |r| (median, 95th percentile), raw and adjusted; the share of draws whose placebo r reaches its
    paired housekeeping r in that r's direction (draws with both defined); whether the matched housekeeping |r|
    exceeds the placebo 95th percentile; mean fill and balance."""

    def beyond(hk: str, p: str) -> pl.Expr:
        hit = (pl.col(p) * pl.col(hk).sign()) >= pl.col(hk).abs()
        return hit.filter(pl.col(hk).is_not_nan() & pl.col(p).is_not_nan()).mean()

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
        placebo_balance_tvd=pl.col("balance_tvd").mean(),
        placebo_balance_unusable=pl.col("balance_unusable").mean(),
        placebos=pl.len(),
        placebos_undefined_adj=(pl.col("adj_note") != "").sum(),
    )
    return band.with_columns(
        exceeds_p95_raw=pl.col("hk_matched_raw").abs() > pl.col("placebo_abs_p95_raw"),
        exceeds_p95_adj=pl.col("hk_matched_adj").abs() > pl.col("placebo_abs_p95_adj"),
    )


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
