"""Mitotic loss at solo-WCGW CpGs in cultures with counted doublings, and the shield there (planned test P10).

Zhou et al. (2018) showed that CpGs with A or T on both sides and no other CpG within 35 bp (solo-WCGW) lose
methylation with each cell division, most in late-replicating, partially methylated domains. Two questions are asked
of every culture series already rebuilt from IDATs:

1. **Does the loss follow doublings, and not time?** Loss is the drop of methylation from the culture's reference
   (the mean of its arrays at the earliest time point), in percentage points, averaged over all autosomal
   solo-WCGW CpGs measured in every array of the culture (and over the 678 PMD solo-WCGWs of HypoClock). Its slope
   on doublings, per 10 doublings, is taken per culture; higher means more loss per division. Arrested arms
   (mitomycin C in Endicott GSE197512, contact inhibition in Sturm GSE179847, six months confluent in E-MTAB-15859)
   are the control: there days pass without doublings. Where arms mix, loss is regressed jointly on doublings and
   days with one intercept per culture.
2. **Is the loss smaller next to housekeeping promoters (the shield)?** Housekeeping promoter-flanking solo-WCGWs
   (450k-annotated CpGs in a housekeeping gene's TSS1500, TSS200, 5'UTR or first exon) against solo-WCGWs matched on
   island relation (island, shore, shelf, open sea) and starting-level decile, from two pools: every other
   solo-WCGW, and promoter solo-WCGWs of other genes. The gap is housekeeping minus matched loss (negative: the
   housekeeping CpGs lose less); its slope per 10 doublings below zero means the shield widens with division, and
   the ratio of the two loss slopes below one means they lose less per division.

Every arm of a culture is its own unit with its own reference (its earliest arrays) and intercept; arrested arms
are paired with the dividing arm of the same culture at the same oxygen, never with a regrowth arm. Intervals: the t
interval over cultures where a series has three or more, and a bootstrap over lines (the biological intervals); a
bootstrap over housekeeping genes for the gap and the ratio, conditional on the cultures. The ratio is given only
where the matched loss slope is positive. Every pooled slope is also given with chip intercepts; a slope whose
exposure the intercepts leave with under 2 % of its variance is not estimated.

Replication timing (housekeeping promoters replicate early, and solo-WCGW loss is a late-replication signal) is
taken up in :func:`timing_shield`: Repli-seq timing of BJ and IMR-90 fibroblasts joins the strata, and the gap is
read within early-, mid- and late-replicating tertiles and within Zhou et al.'s common HMD and PMD bins, with
leave-one-line-out fits and housekeeping genes resampled as blocks. :func:`context_loss` and the WCGW-with-neighbour
and SCGS designs ask whether the division effect and the gap are particular to solo-WCGWs.

``entropy-age figures solo_wcgw_culture`` (needs the SeSAMe rebuilds of the series below).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy import sparse, stats

from entropy_age.data import reprocessed
from entropy_age.data.annotation import PMD_BINS, REPLISEQ_450K, AnnotationLayout, link_membership, solo_wcgw
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments import arrest_timecourse, astrocyte_passage, bj_telomerase
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.passaged_cells import ARM_COLUMNS, _arm
from entropy_age.stochasticity import baseline_strata

ENDICOTT, STURM = "GSE197512", "GSE179847"
#: Sturm arms kept: untreated and contact-inhibited cultures (drugs, 5-azacytidine above all, change methylation)
STURM_AGENTS = ("control", "Contact_Inhibition", "Contact_Inhibition_Regrowth")
ARRESTED = ("agent=MMC", "Contact_Inhibition")  # an arm whose label carries one of these is arrested
DAYS_PER_MONTH = 30.0  # the protocol sampled every 30 days
BANDS = {"all": 0.0, "start >= 0.6": 0.6}
N_BOOT = {"paper": 2000, "quick": 200}
MIN_ARRAYS = 3


def island_class(relation: str) -> str:
    """Island, Shore, Shelf or OpenSea from the manifest's relation (N_/S_ dropped)."""
    return relation.removeprefix("N_").removeprefix("S_")


def matched_weights(count_in: np.ndarray, pool: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Weights ``(s, m)``: ``x @ s`` is the set's mean, ``x @ m`` the pool's mean re-weighted to the set's strata.

    ``count_in`` is how often each CpG counts in the set (a gene resample can count one twice, or not at all); ``pool``
    marks the CpGs a set CpG may be matched to (never a set CpG, drawn or not). Strata without set or without pool
    CpGs are dropped from both sides.
    """
    count_in = np.asarray(count_in, dtype=float)
    pool = np.asarray(pool, dtype=bool) & ~(count_in > 0)
    _, codes = np.unique(strata, return_inverse=True)
    k = codes.max() + 1 if codes.size else 0
    n_in = np.bincount(codes, weights=count_in, minlength=k)
    n_pool = np.bincount(codes, weights=pool.astype(float), minlength=k)
    usable = (n_in > 0) & (n_pool > 0)
    if not usable.any():
        raise ValueError("no stratum holds both set and pool CpGs")
    share = np.where(usable, n_in, 0.0) / n_in[usable].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(usable[codes], count_in * (share / n_in)[codes], 0.0)
        m = np.where(pool & usable[codes], (share / n_pool)[codes], 0.0)
    return np.nan_to_num(s), np.nan_to_num(m)


MIN_LEFT = 0.02  # share of an exposure's variance that must be left after the nuisance terms, else no slope


def left_share(t: np.ndarray, nuisance: np.ndarray) -> float:
    """Share of the variance of ``t`` left after regressing it on ``nuisance`` (intercepts): 0 means unidentified."""
    coef, *_ = np.linalg.lstsq(nuisance, t, rcond=None)
    total = float(((t - t.mean()) ** 2).sum())
    return float(((t - nuisance @ coef) ** 2).sum()) / total if total > 0 else 0.0


def _fit(y: np.ndarray, exposures: list[np.ndarray], nuisance: np.ndarray) -> list[float]:
    """OLS coefficients of ``exposures`` with ``nuisance`` columns; NaN for an exposure the nuisance (and the other
    exposures) leave less than ``MIN_LEFT`` of its variance (a rank-deficient fit gives no slope)."""
    ok = []
    for i, e in enumerate(exposures):
        others = np.column_stack([nuisance, *[x for j, x in enumerate(exposures) if j != i]])
        ok.append(left_share(e, others) >= MIN_LEFT)
    coef, *_ = np.linalg.lstsq(np.column_stack([*exposures, nuisance]), y, rcond=None)
    return [float(c) if good else np.nan for c, good in zip(coef[: len(exposures)], ok, strict=True)]


def within_slope(y: np.ndarray, t: np.ndarray, units: np.ndarray, chips: np.ndarray | None = None) -> float:
    """OLS slope of ``y`` on ``t`` with one intercept per unit (and per chip when ``chips`` is given); NaN when unit
    and chip leave ``t`` (almost) no variance."""
    cols = [np.ones(t.size), dummies(units)]
    if chips is not None:
        cols.append(dummies(chips))
    return _fit(y, [t], np.column_stack(cols))[0]


def joint_rates(y: np.ndarray, doublings: np.ndarray, days: np.ndarray, units: np.ndarray) -> tuple[float, float]:
    """Slopes of ``y`` on doublings and on days fitted together, one intercept per unit (NaN where unidentified)."""
    pd_rate, day_rate = _fit(y, [doublings, days], np.column_stack([np.ones(y.size), dummies(units)]))
    return pd_rate, day_rate


def t_interval(values: np.ndarray) -> dict:
    """Mean, 95 % t interval and two-sided p over replicate slopes (NaN with fewer than three)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return {
            "n": int(v.size),
            "mean": float(v.mean()) if v.size else np.nan,
            "low": np.nan,
            "high": np.nan,
            "p": np.nan,
        }
    se = v.std(ddof=1) / np.sqrt(v.size)
    half = float(stats.t.ppf(0.975, v.size - 1) * se)
    return {
        "n": int(v.size),
        "mean": float(v.mean()),
        "low": float(v.mean() - half),
        "high": float(v.mean() + half),
        "p": float(stats.ttest_1samp(v, 0.0).pvalue),
    }


# --------------------------------------------------------------------------------------------- the culture series


def _endicott(paths) -> pl.DataFrame:
    s = pl.read_parquet(cohort_layout("endicott", paths).samples)
    s = s.with_columns(arm=pl.struct(*ARM_COLUMNS).map_elements(_arm, return_dtype=pl.String))
    return s.select(
        pl.lit("Endicott fibroblasts and others (GSE197512)").alias("series"),
        pl.lit(ENDICOTT).alias("gse"),
        "sample",
        pl.concat_str("subexperiment", pl.lit(" | "), "donor", pl.lit(" | "), "arm").alias("unit"),
        pl.concat_str("subexperiment", pl.lit(" | "), "donor").alias("culture"),
        "subexperiment",
        "donor",
        "arm",
        "doublings",
        "days",
    )


def _sturm(paths) -> pl.DataFrame:
    s = pl.read_parquet(cohort_layout("sturm", paths).samples).filter(
        pl.col("agent").is_in(STURM_AGENTS) & (pl.col("subexperiment") != "SURF1 mutation")  # genetic disease lines
    )
    s = s.with_columns(arm=pl.struct(*ARM_COLUMNS).map_elements(_arm, return_dtype=pl.String))
    return s.select(
        pl.lit("Sturm dermal fibroblasts (GSE179847)").alias("series"),
        pl.lit(STURM).alias("gse"),
        "sample",
        pl.concat_str("subexperiment", pl.lit(" | "), "donor", pl.lit(" | "), "arm").alias("unit"),
        pl.concat_str("subexperiment", pl.lit(" | "), "donor").alias("culture"),
        "subexperiment",
        "donor",
        "arm",
        "doublings",
        "days",
    )


def _bj(paths) -> pl.DataFrame:
    matrix = paths.input / "geo" / bj_telomerase.GSE / f"{bj_telomerase.GSE}_series_matrix.txt.gz"
    if not matrix.exists():
        return pl.DataFrame()
    frames = []
    for name, g in bj_telomerase.groups(bj_telomerase.read_cultures(matrix)).items():
        frames.append(
            g.select(
                pl.lit("BJ fibroblasts, mortal and hTERT (GSE268145)").alias("series"),
                pl.lit(bj_telomerase.GSE).alias("gse"),
                "sample",
                pl.concat_str(pl.lit(name + " | trajectory "), "trajectory").alias("unit"),
                pl.concat_str(pl.lit(name + " | trajectory "), "trajectory").alias("culture"),
                pl.lit(name).alias("subexperiment"),
                pl.concat_str(pl.lit("BJ "), "trajectory").alias("donor"),
                pl.lit("hTERT" if name == bj_telomerase.HTERT else "standard").alias("arm"),
                "doublings",
                pl.lit(None, dtype=pl.Float64).alias("days"),
            )
        )
    return pl.concat(frames)


def _astrocytes(paths, gse: str, select, label: str) -> pl.DataFrame:
    matrix = paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
    if not matrix.exists():
        return pl.DataFrame()
    chars = _characteristics(matrix)
    c = select(chars)
    if "days" in chars.columns:  # days in culture, as recorded, by GSM
        c = c.join(chars.select("sample", pl.col("days").cast(pl.Float64, strict=False)), on="sample", how="left")
    if "days" not in c.columns or c["days"].null_count():
        c = c.with_columns(days=pl.lit(None, dtype=pl.Float64))
    return c.select(
        pl.lit(label).alias("series"),
        pl.lit(gse).alias("gse"),
        "sample",
        pl.col("donor").alias("unit"),
        pl.col("donor").alias("culture"),
        pl.lit(label).alias("subexperiment"),
        "donor",
        pl.lit("standard").alias("arm"),
        "doublings",
        "days",
    )


def _arrest(paths) -> pl.DataFrame:
    acc = arrest_timecourse.ACCESSION
    sdrf = paths.input / "arrayexpress" / acc / f"{acc}.sdrf.txt"
    if not (sdrf.exists() and reprocessed.available(paths, acc)):
        return pl.DataFrame()
    qc = reprocessed.read_qc(paths, acc).filter(pl.col("passed") & (pl.col("platform") == "EPIC"))
    s = arrest_timecourse.read_samples(sdrf).filter(
        (pl.col("line") == "dermal fibroblasts")
        & (pl.col("compound") == "Control")
        & pl.col("donor").is_in(arrest_timecourse.DONORS)
        & pl.col("sample").is_in(qc["sample"].implode())
    )
    frames = []
    for donor in arrest_timecourse.DONORS:  # each donor's main chip only, as in arrest_timecourse
        d = s.filter(pl.col("donor") == donor)
        if d.height == 0:
            continue
        main = d.group_by("chip").len().sort("len", descending=True)["chip"][0]
        frames.append(d.filter(pl.col("chip") == main))
    s = pl.concat(frames)
    return s.select(
        pl.lit("Confluent fibroblasts, six months (E-MTAB-15859)").alias("series"),
        pl.lit(acc).alias("gse"),
        "sample",
        pl.col("donor").alias("unit"),
        pl.col("donor").alias("culture"),
        pl.lit("arrested").alias("subexperiment"),
        "donor",
        pl.lit("confluent, arrested").alias("arm"),
        pl.lit(0.0).alias("doublings"),
        (pl.col("month").cast(pl.Float64) * DAYS_PER_MONTH).alias("days"),
    )


def cultures(paths) -> pl.DataFrame:
    """Every culture array of every rebuilt series: series, gse, sample, unit, arm, doublings, days, chip,
    ``arrested``; arrays failing SeSAMe QC left out."""
    frames = []
    loaders = (
        (ENDICOTT, lambda: _endicott(paths)),
        (STURM, lambda: _sturm(paths)),
        (bj_telomerase.GSE, lambda: _bj(paths)),
        (
            astrocyte_passage.GSE,
            lambda: _astrocytes(
                paths, astrocyte_passage.GSE, astrocyte_passage.bulk_cultures, "fetal astrocytes (GSE202554)"
            ),
        ),
        (
            astrocyte_passage.GSE_HTERT,
            lambda: _astrocytes(
                paths, astrocyte_passage.GSE_HTERT, astrocyte_passage.htert_cultures, "hTERT astrocytes (GSE226079)"
            ),
        ),
        (arrest_timecourse.ACCESSION, lambda: _arrest(paths)),
    )
    for gse, load in loaders:
        if not reprocessed.available(paths, gse):
            print(f"solo_wcgw_culture: no SeSAMe rebuild of {gse}; left out", file=sys.stderr)
            continue
        frame = load()
        if frame.height == 0:
            print(f"solo_wcgw_culture: no sample table for {gse}; left out", file=sys.stderr)
            continue
        qc = reprocessed.read_qc(paths, gse).filter(pl.col("passed")).select("sample", "chip")
        frame = frame.with_columns(pl.col("doublings").cast(pl.Float64), pl.col("days").cast(pl.Float64))
        joined = frame.join(qc, on="sample", how="inner")
        print(f"solo_wcgw_culture: {gse}: {frame.height} arrays selected, {joined.height} pass QC", file=sys.stderr)
        frames.append(joined)
    out = pl.concat(frames, how="vertical_relaxed")
    arrested = pl.lit(False)
    for label in ARRESTED:
        arrested = arrested | (
            pl.col("arm").str.contains(label, literal=True) & ~pl.col("arm").str.contains("Regrowth")
        )
    return out.with_columns(arrested=arrested | (pl.col("arm") == "confluent, arrested"))


# ---------------------------------------------------------------------------------------------------- the CpG sets


def cpg_sets(paths) -> tuple[pl.DataFrame, list[str]]:
    """Per autosomal solo-WCGW CpG: ``island`` class, ``hk`` (housekeeping promoter), ``hk_gene``, ``promoter``,
    ``non_essential`` (promoter of a non-expressed gene) and ``pmd`` (a HypoClock CpG); and the HypoClock list."""
    annotation = AnnotationLayout(paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(pl.col("solo_wcgw") & ~pl.col("chr").is_in(["X", "Y"]))
    membership = link_membership(pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.gene_sets))
    hk = membership.filter(pl.col("group") == "housekeeping").select(pl.col("feature").alias("cpg"), pl.col("gene"))
    hk = hk.unique(subset="cpg", keep="first", maintain_order=True)
    ne = set(membership.filter(pl.col("group") == "non_essential")["feature"].to_list())
    hypo_path = paths.input / "mitotic" / "hypoclock.csv"
    hypo = pl.read_csv(hypo_path)["cpg"].to_list() if hypo_path.exists() else []
    table = (
        probes.select("cpg", "island_relation", "promoter")
        .join(hk.rename({"gene": "hk_gene"}), on="cpg", how="left")
        .with_columns(
            island=pl.col("island_relation").map_elements(island_class, return_dtype=pl.String),
            hk=pl.col("hk_gene").is_not_null(),
            non_essential=pl.col("cpg").is_in(sorted(ne)),
            pmd=pl.col("cpg").is_in(hypo),
        )
        .sort("cpg")
    )
    return table, hypo


# ------------------------------------------------------------------------------------------------------- scoring


def unit_loss(x: np.ndarray, reference_rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(reference, loss)``: the reference level per CpG and the loss of every array, ``reference - x``, in points."""
    reference = x[reference_rows].mean(axis=0)
    return reference, 100.0 * (reference - x)


def score_unit(
    loss: np.ndarray, reference: np.ndarray, table: pl.DataFrame, band_floor: float, counts: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Per array: mean loss over all solo-WCGWs and PMD ones; housekeeping, matched and gap under both pools; and the
    non-essential contrast. ``table`` is aligned with the columns of ``loss``; ``counts`` the housekeeping CpGs'
    multiplicities (a gene resample), default one each."""
    keep = reference >= band_floor
    hk = table["hk"].to_numpy() & keep
    count = (hk.astype(float) if counts is None else counts * keep).astype(float)
    strata = np.char.add(
        np.char.add(table["island"].to_numpy().astype(str), "|"), baseline_strata(reference, 10).astype(str)
    )
    promoter = table["promoter"].to_numpy()
    eligible = keep & ~table["hk"].to_numpy()
    out = {"loss_all": loss[:, keep].mean(axis=1)}
    pmd = table["pmd"].to_numpy() & keep
    out["loss_pmd"] = loss[:, pmd].mean(axis=1) if pmd.sum() >= 20 else np.full(loss.shape[0], np.nan)
    for pool_name, pool in (("other", eligible), ("promoter", eligible & promoter)):
        s, m = matched_weights(count, pool, strata)
        out[f"loss_hk_{pool_name}"] = loss @ s
        out[f"matched_{pool_name}"] = loss @ m
        out[f"gap_{pool_name}"] = loss @ (s - m)
    ne = table["non_essential"].to_numpy() & keep & ~table["hk"].to_numpy()
    if ne.sum() >= 20:
        s, m = matched_weights(ne.astype(float), eligible & ~ne, strata)
        out["gap_non_essential"] = loss @ (s - m)
    else:
        out["gap_non_essential"] = np.full(loss.shape[0], np.nan)
    return out


def reference_rows(unit: pl.DataFrame) -> np.ndarray:
    """Rows of the unit's earliest arrays: earliest day where days are recorded, else fewest doublings."""
    t = unit["days"] if unit["days"].null_count() == 0 else unit["doublings"]
    return np.flatnonzero((t == t.min()).to_numpy())


def with_presplit(unit: pl.DataFrame, culture: pl.DataFrame) -> pl.DataFrame:
    """The unit, preceded by its culture's pre-split arrays when the culture's earliest arrays all belong to the
    ``standard`` arm and this unit is another arm (Endicott's MMC course: day-8 cultures before the split); a unit
    whose culture has no such common start keeps its own earliest arrays as reference."""
    if unit["arm"][0] == "standard" or culture.height == unit.height:
        return unit
    first = culture.filter(pl.int_range(pl.len()).is_in(reference_rows(culture).tolist()))
    if (first["arm"] == "standard").all():
        return pl.concat([first, unit], how="vertical_relaxed")
    return unit


def prepare(paths, samples: pl.DataFrame, table: pl.DataFrame) -> list[dict]:
    """Per unit: the later arrays' metadata, the loss matrix, the reference and the aligned CpG table."""
    units = []
    for gse, frame in samples.group_by("gse", maintain_order=True):
        gse = gse[0]
        betas = reprocessed.by_accession(paths, gse, {s: s for s in frame["sample"].to_list()})
        features = table["cpg"].to_list()
        frame = frame.sort("unit", "days", "doublings")
        for (unit_name,), unit in frame.group_by("unit", maintain_order=True):
            unit = with_presplit(unit, frame.filter(pl.col("culture") == unit["culture"][0]))
            if unit.height < MIN_ARRAYS:
                continue
            x, feats = complete_matrix(betas, features, unit["sample"].to_list())
            ref_rows = reference_rows(unit)
            reference, loss = unit_loss(x, ref_rows)
            later = np.setdiff1d(np.arange(unit.height), ref_rows)
            if later.size < 2:
                continue
            ref = unit.filter(pl.int_range(pl.len()).is_in(ref_rows.tolist()))
            has_days = unit["days"].null_count() == 0
            meta = unit.filter(pl.int_range(pl.len()).is_in(later.tolist())).with_columns(
                doublings_gained=pl.col("doublings") - float(ref["doublings"].mean()),
                days_gained=(pl.col("days") - float(ref["days"].mean()))
                if has_days
                else pl.lit(None, dtype=pl.Float64),
                reference_arrays=pl.lit(int(ref_rows.size)),
                reference_samples=pl.lit(",".join(ref["sample"].to_list())),
                cpgs=pl.lit(len(feats)),
            )
            aligned = pl.DataFrame({"cpg": feats}).join(table, on="cpg", how="left")
            meta = meta.with_columns(
                hk_cpgs=pl.lit(int(aligned["hk"].sum())), pmd_cpgs=pl.lit(int(aligned["pmd"].sum()))
            )
            units.append({"meta": meta, "loss": loss[later], "reference": reference, "table": aligned})
            print(
                f"solo_wcgw_culture: {unit_name}: {later.size} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True
            )
    return units


def scores(units: list[dict], band_floor: float, gene_counts: dict[str, int] | None = None) -> pl.DataFrame:
    frames = []
    for u in units:
        counts = None
        if gene_counts is not None:
            genes = u["table"]["hk_gene"].to_list()
            counts = np.array([gene_counts.get(g, 0) if g is not None else 0 for g in genes], dtype=float)
        s = score_unit(u["loss"], u["reference"], u["table"], band_floor, counts)
        frames.append(u["meta"].with_columns(**{k: pl.Series(v) for k, v in s.items()}))
    return pl.concat(frames, how="vertical_relaxed")


METRICS = (
    "loss_all",
    "loss_pmd",
    "loss_hk_other",
    "matched_other",
    "gap_other",
    "loss_hk_promoter",
    "matched_promoter",
    "gap_promoter",
    "gap_non_essential",
)


# ------------------------------------------------------------------------------------------------------- the fits

#: per series, the arrays read for "loss per doubling": dividing cultures of the standard (or hTERT) course
DIVIDING = {
    "Endicott, 7 primary lines": lambda s: (
        (pl.col("gse") == ENDICOTT) & (pl.col("subexperiment") == "Baseline profiling")
    ),
    "Endicott, TERT line": lambda s: (pl.col("gse") == ENDICOTT) & (pl.col("arm") == "vector=TERT"),
    "Sturm, 5 lines": lambda s: (pl.col("gse") == STURM) & (pl.col("subexperiment") == "Baseline profiling"),
    "BJ mortal, EPICv2 (chips shared with hTERT)": lambda s: (
        (pl.col("gse") == bj_telomerase.GSE) & (pl.col("subexperiment") == bj_telomerase.MORTAL_SAME)
    ),
    "BJ mortal, EPIC (earlier batch)": lambda s: (
        (pl.col("gse") == bj_telomerase.GSE) & (pl.col("subexperiment") == bj_telomerase.MORTAL_EARLIER)
    ),
    "BJ hTERT": lambda s: (pl.col("gse") == bj_telomerase.GSE) & (pl.col("arm") == "hTERT"),
    "astrocytes, 2 donors": lambda s: pl.col("gse") == astrocyte_passage.GSE,
    "hTERT astrocytes, 3 lines": lambda s: pl.col("gse") == astrocyte_passage.GSE_HTERT,
}


def unit_slopes(s: pl.DataFrame) -> pl.DataFrame:
    """Per unit and arm: slope of every metric on doublings (per 10) and on days (per 10)."""
    rows = []
    for (series, unit, arm), g in s.group_by("series", "unit", "arm", maintain_order=True):
        pd_span = float(np.ptp(g["doublings"].to_numpy()))
        for axis, col in (("doublings", "doublings_gained"), ("days", "days_gained")):
            t = g[col].to_numpy().astype(float) / 10.0
            if g.height < 3 or not np.isfinite(t).all() or np.ptp(t) == 0:
                continue
            for metric in METRICS:
                y = g[metric].to_numpy().astype(float)
                if not np.isfinite(y).all():
                    continue
                fit = stats.linregress(t, y)
                rows.append(
                    {
                        "series": series,
                        "unit": unit,
                        "culture": g["culture"][0],
                        "arm": arm,
                        "doublings_in_arm": pd_span,
                        "arrested": bool(g["arrested"][0]),
                        "axis": axis,
                        "metric": metric,
                        "arrays": g.height,
                        "span": float(np.ptp(t) * 10),
                        "slope_per_10": float(fit.slope),
                        "se": float(fit.stderr),
                    }
                )
    return pl.DataFrame(rows)


def cluster_boot(fn, s: pl.DataFrame, cluster: str, n: int, rng) -> np.ndarray:
    """``fn`` over resamples of whole clusters (lines), relabelled so a line drawn twice counts as two units."""
    names = s[cluster].unique(maintain_order=True).to_list()
    parts = {c: s.filter(pl.col(cluster) == c) for c in names}
    out = []
    for _ in range(n):
        draw = rng.choice(len(names), size=len(names), replace=True)
        frame = pl.concat(
            [
                parts[names[i]].with_columns(pl.concat_str("unit", pl.lit(f"#{j}")).alias("unit"))
                for j, i in enumerate(draw)
            ]
        )
        try:
            out.append(fn(frame))
        except (np.linalg.LinAlgError, ValueError):
            continue
    return np.array(out)


def series_fits(s: pl.DataFrame, slopes: pl.DataFrame, n_boot: int, rng) -> pl.DataFrame:
    """Per dividing selection and metric: pooled within-unit slope per 10 doublings (and with chip intercepts), the t
    interval over units and the bootstrap over lines."""
    rows = []
    for label, select in DIVIDING.items():
        g = s.filter(select(s) & ~pl.col("arrested"))
        if g.height < 5:
            continue
        units = g["unit"].to_numpy()
        t = g["doublings_gained"].to_numpy() / 10.0
        lines = g["donor"].n_unique()
        for metric in METRICS:
            y = g[metric].to_numpy().astype(float)
            if not np.isfinite(y).all():
                continue
            per_unit = slopes.filter(
                (pl.col("metric") == metric)
                & (pl.col("axis") == "doublings")
                & pl.col("unit").is_in(g["unit"].unique().implode())
                & pl.col("arm").is_in(g["arm"].unique().implode())
            )["slope_per_10"].to_numpy()
            ti = t_interval(per_unit)
            row = {
                "selection": label,
                "metric": metric,
                "arrays": g.height,
                "units": int(np.unique(units).size),
                "lines": lines,
                "doublings_span": float(np.ptp(g["doublings_gained"].to_numpy())),
                "slope_per_10_pd": within_slope(y, t, units),
                "slope_chip_held": within_slope(y, t, units, g["chip"].to_numpy()),
                "units_mean": ti["mean"],
                "units_low": ti["low"],
                "units_high": ti["high"],
                "units_p": ti["p"],
                "units_positive": int((per_unit > 0).sum()),
            }
            if lines >= 4:
                b = cluster_boot(
                    lambda f, m=metric: within_slope(
                        f[m].to_numpy().astype(float), f["doublings_gained"].to_numpy() / 10.0, f["unit"].to_numpy()
                    ),
                    g,
                    "donor",
                    n_boot,
                    rng,
                )
                row |= {"lines_low": float(np.quantile(b, 0.025)), "lines_high": float(np.quantile(b, 0.975))}
            rows.append(row)
    return pl.DataFrame(rows)


def joint_fits(s: pl.DataFrame, n_boot: int, rng) -> pl.DataFrame:
    """Within each series with arms that part doublings from days (Endicott, Sturm): every metric on doublings (per 10)
    and days (per 10) together, one intercept per culture; bootstrap over lines."""
    rows = []
    for gse in (ENDICOTT, STURM):
        g = s.filter(
            (pl.col("gse") == gse)
            & pl.col("days_gained").is_not_null()
            & ~pl.col("arm").str.starts_with("vector=")  # immortalised and vector-control arms: their own rows
        )
        if g.height < 10:
            continue

        def fit(f, metric):
            return joint_rates(
                f[metric].to_numpy().astype(float),
                f["doublings_gained"].to_numpy() / 10.0,
                f["days_gained"].to_numpy() / 10.0,
                f["unit"].to_numpy(),
            )

        for metric in METRICS:
            if not np.isfinite(g[metric].to_numpy().astype(float)).all():
                continue
            pd_rate, day_rate = fit(g, metric)
            b = cluster_boot(lambda f, m=metric: fit(f, m), g, "donor", n_boot, rng)
            rows.append(
                {
                    "series": g["series"][0],
                    "metric": metric,
                    "arrays": g.height,
                    "units": g["unit"].n_unique(),
                    "lines": g["donor"].n_unique(),
                    "arrested_arrays": int(g["arrested"].sum()),
                    "pd_left_share": left_share(
                        g["doublings_gained"].to_numpy(),
                        np.column_stack(
                            [np.ones(g.height), dummies(g["unit"].to_numpy()), g["days_gained"].to_numpy()]
                        ),
                    ),
                    "per_10_pd": pd_rate,
                    "pd_low": float(np.nanquantile(b[:, 0], 0.025)),
                    "pd_high": float(np.nanquantile(b[:, 0], 0.975)),
                    "per_10_days": day_rate,
                    "days_low": float(np.nanquantile(b[:, 1], 0.025)),
                    "days_high": float(np.nanquantile(b[:, 1], 0.975)),
                }
            )
    return pl.DataFrame(rows)


def _oxygen(arm: str) -> str:
    """The oxygen part of an arm label (``ambient`` when none is named)."""
    return next((p for p in arm.split(", ") if p.startswith("oxygen=")), "ambient")


def arrest_contrast(slopes: pl.DataFrame) -> pl.DataFrame:
    """Loss per 10 days in each arrested arm beside the dividing arm of the same culture at the same oxygen (regrowth
    arms never count as the dividing sister); the confluent series has none."""
    days = slopes.filter((pl.col("axis") == "days") & pl.col("metric").is_in(["loss_all", "loss_pmd", "gap_other"]))
    days = days.with_columns(oxygen=pl.col("arm").map_elements(_oxygen, return_dtype=pl.String))
    arrested = days.filter(pl.col("arrested"))
    rows = []
    for r in arrested.iter_rows(named=True):
        sister = days.filter(
            (pl.col("culture") == r["culture"])
            & ~pl.col("arrested")
            & (pl.col("metric") == r["metric"])
            & ~pl.col("arm").str.contains("Regrowth")
            & (pl.col("oxygen") == r["oxygen"])
        )
        rows.append(
            {
                **{k: r[k] for k in ("series", "culture", "arm", "metric", "arrays", "span", "doublings_in_arm")},
                "arrested_per_10_days": r["slope_per_10"],
                "arrested_se": r["se"],
                "dividing_arm": ", ".join(sister["arm"].to_list()) or None,
                "dividing_per_10_days": float(sister["slope_per_10"].mean()) if sister.height else None,
            }
        )
    return pl.DataFrame(rows)


def gene_bootstrap(units: list[dict], table: pl.DataFrame, n: int, rng) -> pl.DataFrame:
    """Per dividing selection: housekeeping-minus-matched gap per 10 doublings and the ratio of the loss slopes, with
    intervals over resampled housekeeping genes (band ``all``)."""
    genes = sorted(set(table.filter(pl.col("hk"))["hk_gene"].to_list()))
    base = scores(units, 0.0)
    draws = []
    for _ in range(n):
        picked = rng.choice(len(genes), size=len(genes), replace=True)
        counts: dict[str, int] = {}
        for i in picked:
            counts[genes[i]] = counts.get(genes[i], 0) + 1
        draws.append(scores(units, 0.0, counts))
    rows = []
    for label, select in DIVIDING.items():
        for pool in ("other", "promoter"):

            def stat(f, select=select, pool=pool):
                g = f.filter(select(f) & ~pl.col("arrested"))
                if g.height < 5:
                    return None
                t, u = g["doublings_gained"].to_numpy() / 10.0, g["unit"].to_numpy()
                hk = within_slope(g[f"loss_hk_{pool}"].to_numpy(), t, u)
                matched = within_slope(g[f"matched_{pool}"].to_numpy(), t, u)
                return hk - matched, hk / matched if matched > 0 else np.nan, hk, matched

            point = stat(base)
            if point is None:
                continue
            boot = np.array([stat(d) for d in draws], dtype=float)
            rows.append(
                {
                    "selection": label,
                    "pool": pool,
                    "hk_loss_per_10_pd": point[2],
                    "matched_loss_per_10_pd": point[3],
                    "gap_per_10_pd": point[0],
                    "gap_low": float(np.nanquantile(boot[:, 0], 0.025)),
                    "gap_high": float(np.nanquantile(boot[:, 0], 0.975)),
                    "ratio": point[1],
                    "ratio_low": float(np.nanquantile(boot[:, 1], 0.025)),
                    "ratio_high": float(np.nanquantile(boot[:, 1], 0.975)),
                    "ratio_valid_draws": int(np.isfinite(boot[:, 1]).sum()),
                    "genes": len(genes),
                    "resamples": n,
                }
            )
    return pl.DataFrame(rows)


# ------------------------------------------------------------- replication timing, and a comparator outside solo-WCGW

TIMING_HALF_WINDOW = 25_000  # bp each side of the CpG over which the Repli-seq wave signal is averaged
TIMING_STRATA = 5  # quintiles of timing, edges from the autosomal solo-WCGWs, in the matching strata
CONTEXTS = ("solo_wcgw", "wcgw", "scgs")
MIN_SET = 20  # fewer housekeeping CpGs than this in a design (median over units): reported, flagged ``small``


def flank_context(sequence: str, flank: int = 35) -> str | None:
    """The probed CpG's context from the manifest's ``Forward_Sequence``: ``solo_wcgw`` (A or T on both sides, no
    other CpG within ``flank`` bp), ``wcgw`` (A or T on both sides with another CpG within ``flank`` bp), ``scgs`` (C or
    G on both sides) or ``other``; None when the bracket or the context is missing."""
    i = sequence.find("[CG]")
    if i < flank + 1 or len(sequence) < i + 4 + flank + 1:
        return None
    left, right = sequence[i - 1].upper(), sequence[i + 4].upper()
    if left in "AT" and right in "AT":
        return "solo_wcgw" if solo_wcgw(sequence, flank) else "wcgw"
    if left in "CG" and right in "CG":
        return "scgs"
    return "other"


def timing_classes(timing: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """``late``, ``mid`` or ``early`` by the two tertile ``edges`` of the wave signal (higher = earlier); NaN timing
    is ``unknown``."""
    out = np.where(timing > edges[1], "early", np.where(timing > edges[0], "mid", "late")).astype(object)
    out[~np.isfinite(timing)] = "unknown"
    return out.astype(str)


def fixed_strata(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Bin index of every value by fixed inner ``edges``; NaN gets ``len(edges) + 1``."""
    out = np.searchsorted(edges, values, side="right")
    return np.where(np.isfinite(values), out, len(edges) + 1)


def replication_annotation(paths, probes: pl.DataFrame) -> pl.DataFrame:
    """``cpg``, ``timing`` (mean Repli-seq wave signal of BJ and IMR-90 fibroblasts over the CpG's 50 kb; higher =
    earlier) and ``domain`` (``commonPMD``, ``commonHMD`` or ``neither``, Zhou et al. 2018's 100 kb bins) for every
    probe of ``probes`` (``cpg``, ``chr``, ``pos`` on GRCh37)."""
    import pyBigWig

    root = paths.input / "annotation"
    waves = []
    for name in REPLISEQ_450K.values():
        bw = pyBigWig.open(str(root / name))
        try:
            chroms = bw.chroms()
            values = []
            for c, x in zip(probes["chr"].to_list(), probes["pos"].to_list(), strict=True):
                chrom = f"chr{c}"
                if chrom not in chroms or x is None:
                    values.append(np.nan)
                    continue
                at = int(x) - 1  # MAPINFO is 1-based; bigWig intervals are 0-based
                lo, hi = max(at - TIMING_HALF_WINDOW, 0), min(at + TIMING_HALF_WINDOW + 1, chroms[chrom])
                v = bw.stats(chrom, lo, hi, type="mean")[0]
                values.append(np.nan if v is None else v)
            waves.append(np.array(values, dtype=float))
        finally:
            bw.close()
    stacked = np.vstack(waves)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # a CpG neither track covers: NaN
        timing = np.nanmean(stacked, axis=0)
    bins = pl.read_csv(
        root / PMD_BINS,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "depth", "call", "domain"],
        schema_overrides={"start": pl.Int64},
    ).select("chrom", pl.col("start").alias("bin"), "domain")
    out = probes.select("cpg", "chr", "pos").with_columns(
        timing=pl.Series(timing),
        chrom=pl.lit("chr") + pl.col("chr"),
        bin=((pl.col("pos") - 1) // 100_000) * 100_000,
    )
    return out.join(bins, on=["chrom", "bin"], how="left").select(
        "cpg", "timing", pl.col("domain").fill_null("unknown")
    )


def cpg_universe(paths) -> pl.DataFrame:
    """Every autosomal 450k CpG of the three contexts (solo-WCGW, WCGW with a CpG neighbour, SCGS) with the columns
    of :func:`cpg_sets` plus ``context``, ``timing``, ``timing_class``, ``timing_stratum`` and ``domain``. Timing
    edges (tertiles for the class, quintiles for the strata) come from the autosomal solo-WCGWs and are applied to
    every context."""
    from entropy_age.data.annotation import MANIFEST_450K, _assay_block

    annotation = AnnotationLayout(paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(~pl.col("chr").is_in(["X", "Y"]))
    sequences = _assay_block(annotation.input / MANIFEST_450K).select(
        pl.col("Name").alias("cpg"), pl.col("Forward_Sequence").fill_null("").alias("sequence")
    )
    probes = probes.join(sequences, on="cpg", how="inner").with_columns(
        context=pl.col("sequence").map_elements(flank_context, return_dtype=pl.String)
    )
    probes = probes.filter(pl.col("context").is_in(CONTEXTS))
    membership = link_membership(pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.gene_sets))
    hk = membership.filter(pl.col("group") == "housekeeping").select(pl.col("feature").alias("cpg"), pl.col("gene"))
    hk = hk.unique(subset="cpg", keep="first", maintain_order=True)
    ne = membership.filter(pl.col("group") == "non_essential")["feature"].unique()
    hypo_path = paths.input / "mitotic" / "hypoclock.csv"
    hypo = pl.read_csv(hypo_path)["cpg"] if hypo_path.exists() else pl.Series("cpg", [], dtype=pl.String)
    table = (
        probes.select("cpg", "chr", "pos", "island_relation", "promoter", "context")
        .join(hk.rename({"gene": "hk_gene"}), on="cpg", how="left")
        .with_columns(
            island=pl.col("island_relation").str.strip_prefix("N_").str.strip_prefix("S_"),
            hk=pl.col("hk_gene").is_not_null(),
            non_essential=pl.col("cpg").is_in(ne.implode()),
            pmd=pl.col("cpg").is_in(hypo.implode()),
        )
    )
    table = table.join(replication_annotation(paths, table), on="cpg", how="left")
    solo = table.filter(pl.col("context") == "solo_wcgw")["timing"].to_numpy()
    tertiles = np.nanquantile(solo, [1 / 3, 2 / 3])
    inner = np.nanquantile(solo, np.linspace(0, 1, TIMING_STRATA + 1)[1:-1])
    timing = table["timing"].to_numpy().astype(float)
    return (
        table.with_columns(
            timing_class=pl.Series(timing_classes(timing, tertiles)),
            timing_stratum=pl.Series(fixed_strata(timing, inner)),
        )
        .drop("chr", "pos")
        .sort("cpg")
    )


def subset_units(units: list[dict], context: str) -> list[dict]:
    """The units restricted to one context's CpGs, with the CpG counts in their metadata recomputed."""
    out = []
    for u in units:
        keep = (u["table"]["context"] == context).to_numpy()
        table = u["table"].filter(pl.Series(keep))
        meta = u["meta"].with_columns(
            cpgs=pl.lit(table.height), hk_cpgs=pl.lit(int(table["hk"].sum())), pmd_cpgs=pl.lit(int(table["pmd"].sum()))
        )
        out.append(
            {
                "meta": meta,
                "loss": np.ascontiguousarray(u["loss"][:, keep]),
                "reference": u["reference"][keep],
                "table": table,
            }
        )
    return out


#: design -> (context, restriction of both set and pool, timing in the strata). The first repeats the main gap.
TIMING_DESIGNS: dict[str, tuple[str, str | None, bool]] = {
    "island x start": ("solo_wcgw", None, False),
    "island x start x timing": ("solo_wcgw", None, True),
    "early-replicating": ("solo_wcgw", "early", True),
    "mid-replicating": ("solo_wcgw", "mid", True),
    "late-replicating": ("solo_wcgw", "late", True),
    "common HMD": ("solo_wcgw", "commonHMD", True),
    "outside common HMD": ("solo_wcgw", "notHMD", True),
    "common PMD": ("solo_wcgw", "commonPMD", True),
    "WCGW with a CpG neighbour": ("wcgw", None, True),
    "WCGW with a CpG neighbour, early": ("wcgw", "early", True),
    "SCGS": ("scgs", None, True),
    "SCGS, early": ("scgs", "early", True),
}


def restriction(table: pl.DataFrame, rule: str | None) -> np.ndarray:
    """CpGs a design keeps: all, a timing class, a domain, or outside common HMDs (``notHMD``)."""
    if rule is None:
        return np.ones(table.height, dtype=bool)
    if rule in ("early", "mid", "late"):
        return (table["timing_class"] == rule).to_numpy()
    if rule == "notHMD":
        return table["domain"].is_in(["commonPMD", "neither"]).to_numpy()
    return (table["domain"] == rule).to_numpy()


def design_parts(u: dict, design: tuple[str, str | None, bool], genes: dict[str, int]) -> dict | None:
    """What a gene draw needs to score one unit under one design: the housekeeping CpGs' losses, strata and gene
    indices, and per pool the pool's mean loss in every stratum (arrays x strata) and its size."""
    context, rule, timed = design
    table = u["table"]
    keep = (table["context"] == context).to_numpy() & restriction(table, rule)
    if keep.sum() < 2 * MIN_SET:
        return None
    reference = u["reference"]
    level = np.full(table.height, -1)
    level[keep] = baseline_strata(reference[keep], 10)
    parts = [table["island"].to_numpy().astype(str), level.astype(str)]
    if timed:
        parts.append(table["timing_stratum"].to_numpy().astype(str))
    labels = parts[0]
    for p in parts[1:]:
        labels = np.char.add(np.char.add(labels, "|"), p)
    _, codes = np.unique(labels, return_inverse=True)
    k = int(codes.max()) + 1
    hk = keep & table["hk"].to_numpy()
    if not hk.any():
        return None
    gene_names = table["hk_gene"].to_list()
    out = {
        "hk_loss": u["loss"][:, hk],
        "hk_codes": codes[hk],
        "hk_gene": np.array([genes[gene_names[i]] for i in np.flatnonzero(hk)]),
        "hk_timing": table["timing"].to_numpy()[hk],
        "hk_start": reference[hk],
        "hk_pmd": (table["domain"] == "commonPMD").to_numpy()[hk],
        "pools": {},
    }
    eligible = keep & ~table["hk"].to_numpy()
    for pool_name, pool in (("other", eligible), ("promoter", eligible & table["promoter"].to_numpy())):
        n = np.bincount(codes[pool], minlength=k).astype(float)
        if not (n > 0).any():
            continue
        size = int(pool.sum())
        indicator = sparse.csr_matrix((np.ones(size), (np.arange(size), codes[pool])), shape=(size, k))
        sums = np.asarray((indicator.T @ u["loss"][:, pool].T).T)
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = np.where(n > 0, sums / n, 0.0)
        idx = np.flatnonzero(pool)

        def stratum_mean(values, idx=idx, k=k):
            v = values[idx]
            ok = np.isfinite(v)
            with np.errstate(divide="ignore", invalid="ignore"):
                num = np.bincount(codes[idx][ok], weights=v[ok], minlength=k)
                return num / np.bincount(codes[idx][ok], minlength=k)  # NaN where no annotated pool CpG

        out["pools"][pool_name] = {
            "mean": mean,
            "n": n,
            "timing": stratum_mean(table["timing"].to_numpy().astype(float)),
            "start": stratum_mean(reference),
            "pmd": stratum_mean((table["domain"] == "commonPMD").to_numpy().astype(float)),
        }
    return out


def score_parts(parts: dict, pool: str, gene_counts: np.ndarray) -> tuple[np.ndarray, np.ndarray, int] | None:
    """Per array: housekeeping loss and matched loss under ``gene_counts`` (multiplicity of every housekeeping gene),
    and the weighted count of housekeeping CpGs used. Same weights as :func:`matched_weights`."""
    p = parts["pools"].get(pool)
    if p is None:
        return None
    c = gene_counts[parts["hk_gene"]].astype(float) * (p["n"][parts["hk_codes"]] > 0)
    total = c.sum()
    if total <= 0:
        return None
    w = c / total
    share = np.bincount(parts["hk_codes"], weights=w, minlength=p["n"].size)
    return parts["hk_loss"] @ w, p["mean"] @ share, int((c > 0).sum())


def balance(parts: dict, pool: str) -> dict | None:
    """Housekeeping against matched means of timing, starting level and common-PMD share (all genes once)."""
    p = parts["pools"].get(pool)
    if p is None:
        return None
    use = p["n"][parts["hk_codes"]] > 0
    if not use.any():
        return None
    w = use / use.sum()
    share = np.bincount(parts["hk_codes"], weights=w, minlength=p["n"].size)

    def mean(values, weights):  # weighted mean over annotated entries, weights renormalised
        ok = np.isfinite(values) & (weights > 0)
        return float(values[ok] @ weights[ok] / weights[ok].sum()) if ok.any() else np.nan

    return {
        "hk_timing": mean(parts["hk_timing"], w),
        "matched_timing": mean(p["timing"], share),
        "hk_start": mean(parts["hk_start"], w),
        "matched_start": mean(p["start"], share),
        "hk_pmd_share": mean(parts["hk_pmd"].astype(float), w),
        "matched_pmd_share": mean(p["pmd"], share),
    }


class SlopeFrame:
    """Within-unit slopes on doublings for one dividing selection, by sums per line so that line resamples and
    leave-one-line-out fits are sums of per-line parts (a line drawn twice counts as two units, as in
    :func:`cluster_boot`)."""

    def __init__(self, meta: pl.DataFrame):
        self.meta = meta
        t = meta["doublings_gained"].to_numpy().astype(float) / 10.0
        units = meta["unit"].to_numpy()
        centred = t.copy()
        for u in np.unique(units):
            sel = units == u
            centred[sel] = t[sel] - t[sel].mean()
        self.centred = centred
        self.lines = meta["donor"].unique(maintain_order=True).to_list()
        donor = meta["donor"].to_numpy()
        self.line_index = np.array([self.lines.index(d) for d in donor])
        self.den = np.bincount(self.line_index, weights=centred**2, minlength=len(self.lines))
        self.sum_t = np.bincount(self.line_index, weights=t, minlength=len(self.lines))
        self.sum_t2 = np.bincount(self.line_index, weights=t**2, minlength=len(self.lines))
        self.count = np.bincount(self.line_index, minlength=len(self.lines)).astype(float)

    def _guard(self, num, den, sum_t, sum_t2, count):
        """``num / den``, NaN where unit intercepts leave under ``MIN_LEFT`` of the exposure's variance (the rule
        of :func:`_fit`)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            total = sum_t2 - sum_t**2 / count
            ok = (total > 0) & (den / total >= MIN_LEFT)
            return np.where(ok, num / den, np.nan)

    def per_line_num(self, y: np.ndarray) -> np.ndarray:
        return np.bincount(self.line_index, weights=self.centred * y, minlength=len(self.lines))

    def slope(self, y: np.ndarray) -> float:
        return float(
            self._guard(
                self.per_line_num(y).sum(), self.den.sum(), self.sum_t.sum(), self.sum_t2.sum(), self.count.sum()
            )
        )

    def line_slopes(self, y: np.ndarray) -> np.ndarray:
        return self._guard(self.per_line_num(y), self.den, self.sum_t, self.sum_t2, self.count)

    def leave_one_out(self, y: np.ndarray) -> np.ndarray:
        num = self.per_line_num(y)
        return self._guard(
            num.sum() - num,
            self.den.sum() - self.den,
            self.sum_t.sum() - self.sum_t,
            self.sum_t2.sum() - self.sum_t2,
            self.count.sum() - self.count,
        )

    def line_boot(self, y: np.ndarray, draws: np.ndarray) -> np.ndarray:
        """Slopes over line resamples ``draws`` (resamples x lines of line indices)."""
        num = self.per_line_num(y)
        return self._guard(
            num[draws].sum(axis=1),
            self.den[draws].sum(axis=1),
            self.sum_t[draws].sum(axis=1),
            self.sum_t2[draws].sum(axis=1),
            self.count[draws].sum(axis=1),
        )


def _interval(values: np.ndarray) -> tuple[float, float]:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))) if v.size >= 20 else (np.nan, np.nan)


def timing_shield(units: list[dict], n_genes: int, n_lines: int, rng) -> dict[str, pl.DataFrame]:
    """The housekeeping gap per 10 doublings under every design of :data:`TIMING_DESIGNS`, per dividing selection:
    housekeeping and matched loss slopes, the gap, its interval over resampled housekeeping genes (blocks: all CpGs
    of a gene go together; cultures held), over resampled lines (four or more), leave-one-line-out estimates and the
    gap with chip intercepts. Also the balance of timing, start and PMD share after matching."""
    genes_all = sorted({g for u in units for g in u["table"].filter(pl.col("hk"))["hk_gene"].to_list()})
    genes = {g: i for i, g in enumerate(genes_all)}
    parts = {name: [design_parts(u, d, genes) for u in units] for name, d in TIMING_DESIGNS.items()}
    meta = pl.concat([u["meta"].with_columns(row_unit=pl.lit(i)) for i, u in enumerate(units)], how="vertical_relaxed")
    ones = np.ones(len(genes_all))
    n = len(genes_all)
    gene_draws = [np.bincount(rng.choice(n, n), minlength=n) for _ in range(n_genes)]

    def per_array(name, pool, counts):
        hk, matched, used = [], [], []
        for u, p in zip(units, parts[name], strict=True):
            r = score_parts(p, pool, counts) if p is not None else None
            if r is None:
                nan = np.full(u["loss"].shape[0], np.nan)
                hk.append(nan), matched.append(nan), used.append(0)
                continue
            hk.append(r[0]), matched.append(r[1]), used.append(r[2])
        return np.concatenate(hk), np.concatenate(matched), used

    rows, loo_rows, balance_rows = [], [], []
    selections = {}
    for label, select in DIVIDING.items():
        mask = (meta.select((select(meta) & ~pl.col("arrested")).alias("m"))["m"]).to_numpy()
        if mask.sum() >= 5:
            selections[label] = (mask, SlopeFrame(meta.filter(pl.Series(mask))))
    for name, (context, rule, timed) in TIMING_DESIGNS.items():
        for pool in ("other", "promoter"):
            hk, matched, used = per_array(name, pool, ones)
            if not np.isfinite(hk).any():
                continue
            boots = [per_array(name, pool, c)[:2] for c in gene_draws]
            for label, (mask, frame) in selections.items():
                scored = mask & np.isfinite(hk) & np.isfinite(matched)
                dropped = int(np.unique(meta["row_unit"].to_numpy()[mask & ~scored]).size)
                if scored.sum() < 5:
                    rows.append({"selection": label, "design": name, "pool": pool, "status": "no scored units"})
                    continue
                if dropped:
                    mask, frame = scored, SlopeFrame(meta.filter(pl.Series(scored)))
                g_hk, g_m = hk[mask], matched[mask]
                sel_units = np.unique(meta["row_unit"].to_numpy()[mask])
                hk_used = float(np.median([used[i] for i in sel_units]))
                gap = g_hk - g_m
                gb = np.array([frame.slope(b[0][mask] - b[1][mask]) for b in boots])
                hb = np.array([frame.slope(b[0][mask]) for b in boots])
                mb = np.array([frame.slope(b[1][mask]) for b in boots])
                loo = frame.leave_one_out(gap)
                per_line = frame.line_slopes(gap)
                row = {
                    "selection": label,
                    "design": name,
                    "context": context,
                    "restriction": rule or "all",
                    "timing_matched": timed,
                    "pool": pool,
                    "status": "ok",
                    "units_unscored": dropped,
                    "hk_cpgs": hk_used,
                    "small": hk_used < MIN_SET,
                    "arrays": int(mask.sum()),
                    "lines": len(frame.lines),
                    "hk_loss_per_10_pd": frame.slope(g_hk),
                    "hk_low_genes": _interval(hb)[0],
                    "hk_high_genes": _interval(hb)[1],
                    "matched_loss_per_10_pd": frame.slope(g_m),
                    "matched_low_genes": _interval(mb)[0],
                    "matched_high_genes": _interval(mb)[1],
                    "gap_per_10_pd": frame.slope(gap),
                    "gap_low_genes": _interval(gb)[0],
                    "gap_high_genes": _interval(gb)[1],
                    "gene_draws_valid": int(np.isfinite(gb).sum()),
                    "gap_chip_held": within_slope(
                        gap,
                        frame.meta["doublings_gained"].to_numpy() / 10.0,
                        frame.meta["unit"].to_numpy(),
                        frame.meta["chip"].to_numpy(),
                    ),
                    "lines_negative": int((per_line[np.isfinite(per_line)] < 0).sum()),
                    "loo_min": float(np.nanmin(loo)) if len(frame.lines) >= 3 else np.nan,
                    "loo_max": float(np.nanmax(loo)) if len(frame.lines) >= 3 else np.nan,
                }
                if len(frame.lines) >= 4:
                    draws = rng.choice(len(frame.lines), size=(n_lines, len(frame.lines)))
                    row["gap_low_lines"], row["gap_high_lines"] = _interval(frame.line_boot(gap, draws))
                    row["hk_low_lines"], row["hk_high_lines"] = _interval(frame.line_boot(g_hk, draws))
                    row["matched_low_lines"], row["matched_high_lines"] = _interval(frame.line_boot(g_m, draws))
                rows.append(row)
                if len(frame.lines) >= 3:
                    loo_hk, loo_m = frame.leave_one_out(g_hk), frame.leave_one_out(g_m)
                    for i, line in enumerate(frame.lines):
                        loo_rows.append(
                            {
                                "selection": label,
                                "design": name,
                                "pool": pool,
                                "left_out": line,
                                "gap_per_10_pd": float(loo[i]),
                                "hk_loss_per_10_pd": float(loo_hk[i]),
                                "matched_loss_per_10_pd": float(loo_m[i]),
                                "line_gap_per_10_pd": float(per_line[i]),
                            }
                        )
                bal = [balance(parts[name][i], pool) for i in sel_units if parts[name][i] is not None]
                bal = [b for b in bal if b is not None]
                if bal:
                    balance_rows.append(
                        {"selection": label, "design": name, "pool": pool, "units": len(bal)}
                        | {k: float(np.mean([b[k] for b in bal])) for k in bal[0]}
                    )
    return {
        "timing_shield": pl.DataFrame(rows),
        "timing_loo": pl.DataFrame(loo_rows),
        "timing_balance": pl.DataFrame(balance_rows),
    }


START_EDGES = np.linspace(0.1, 0.9, 9)  # fixed starting-level bins for comparing contexts


def context_loss(units: list[dict], n_lines: int, rng) -> pl.DataFrame:
    """Is the loss per division particular to solo-WCGWs? Per dividing selection and context (WCGW with a CpG
    neighbour, SCGS): loss per 10 doublings of the context's non-housekeeping CpGs re-weighted to the solo-WCGWs'
    strata (island relation x fixed starting-level bin x timing quintile), beside the solo-WCGWs' own loss over the
    strata both share; the difference (context minus solo-WCGW) below zero means the context loses less."""
    per_unit = {c: ([], []) for c in CONTEXTS[1:]}
    for u in units:
        table = u["table"]
        strata = np.char.add(
            np.char.add(table["island"].to_numpy().astype(str), "|"),
            np.char.add(
                fixed_strata(u["reference"], START_EDGES).astype(str),
                np.char.add("|", table["timing_stratum"].to_numpy().astype(str)),
            ),
        )
        context = table["context"].to_numpy()
        base = ~table["hk"].to_numpy()
        solo = (base & (context == "solo_wcgw")).astype(float)
        for c in CONTEXTS[1:]:
            try:
                s, m = matched_weights(solo, base & (context == c), strata)
            except ValueError:  # no shared stratum: this unit gives no comparison
                nan = np.full(u["loss"].shape[0], np.nan)
                per_unit[c][0].append(nan), per_unit[c][1].append(nan)
                continue
            per_unit[c][0].append(u["loss"] @ s)
            per_unit[c][1].append(u["loss"] @ m)
    meta = pl.concat([u["meta"] for u in units], how="vertical_relaxed")
    rows = []
    for label, select in DIVIDING.items():
        mask = meta.select((select(meta) & ~pl.col("arrested")).alias("m"))["m"].to_numpy()
        if mask.sum() < 5:
            continue
        for c in CONTEXTS[1:]:
            all_solo, all_ctx = np.concatenate(per_unit[c][0]), np.concatenate(per_unit[c][1])
            scored = mask & np.isfinite(all_solo) & np.isfinite(all_ctx)
            if scored.sum() < 5:
                rows.append({"selection": label, "context": c, "status": "no scored units"})
                continue
            frame = SlopeFrame(meta.filter(pl.Series(scored)))
            draws = rng.choice(len(frame.lines), size=(n_lines, len(frame.lines))) if len(frame.lines) >= 4 else None
            solo_y, ctx_y = all_solo[scored], all_ctx[scored]
            diff = ctx_y - solo_y
            loo = frame.leave_one_out(diff)
            row = {
                "selection": label,
                "context": c,
                "status": "ok",
                "arrays_unscored": int((mask & ~scored).sum()),
                "lines": len(frame.lines),
                "solo_wcgw_per_10_pd": frame.slope(solo_y),
                "context_per_10_pd": frame.slope(ctx_y),
                "difference_per_10_pd": frame.slope(diff),
                "difference_chip_held": within_slope(
                    diff,
                    frame.meta["doublings_gained"].to_numpy() / 10.0,
                    frame.meta["unit"].to_numpy(),
                    frame.meta["chip"].to_numpy(),
                ),
                "loo_min": float(np.nanmin(loo)) if len(frame.lines) >= 3 else np.nan,
                "loo_max": float(np.nanmax(loo)) if len(frame.lines) >= 3 else np.nan,
            }
            if draws is not None:
                row["solo_low"], row["solo_high"] = _interval(frame.line_boot(solo_y, draws))
                row["context_low"], row["context_high"] = _interval(frame.line_boot(ctx_y, draws))
                row["difference_low"], row["difference_high"] = _interval(frame.line_boot(diff, draws))
            rows.append(row)
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    out = ctx.out_dir("solo_wcgw_culture")
    n_boot = N_BOOT[str(ctx.scale.name)]
    table, hypo = cpg_sets(paths)
    universe = cpg_universe(paths)
    samples = cultures(paths)
    write_table(
        samples.group_by("series", "unit", "arm", "arrested", maintain_order=True).agg(
            pl.len().alias("arrays"),
            pl.col("doublings").min().alias("doublings_min"),
            pl.col("doublings").max().alias("doublings_max"),
            pl.col("days").max().alias("days_max"),
            pl.col("chip").n_unique().alias("chips"),
        ),
        out / "design",
    )
    every = prepare(paths, samples, universe)
    units = subset_units(every, "solo_wcgw")
    solo = universe.filter(pl.col("context") == "solo_wcgw")
    if solo["cpg"].to_list() != table["cpg"].to_list() or solo["hk"].to_list() != table["hk"].to_list():
        raise ValueError("the solo-WCGW CpGs of the universe differ from cpg_sets")
    write_table(
        pl.DataFrame(
            {
                "set": [
                    "solo-WCGW",
                    "housekeeping promoter",
                    "promoter, other genes",
                    "non-essential promoter",
                    "HypoClock PMD",
                ],
                "cpgs": [
                    table.height,
                    int(table["hk"].sum()),
                    int((table["promoter"] & ~table["hk"]).sum()),
                    int(table["non_essential"].sum()),
                    int(table["pmd"].sum()),
                ],
                "genes": [None, table.filter("hk")["hk_gene"].n_unique(), None, None, None],
            }
        ),
        out / "set_sizes",
    )
    fits, slope_tables = [], []
    for band, floor in BANDS.items():
        s = scores(units, floor)
        if band == "all":
            write_table(s, out / "array_scores")
        sl = unit_slopes(s)
        slope_tables.append(sl.with_columns(band=pl.lit(band)))
        fits.append(series_fits(s, sl, n_boot, ctx.rng(310_000)).with_columns(band=pl.lit(band)))
        if band == "all":
            write_table(joint_fits(s, n_boot, ctx.rng(310_001)), out / "joint_fits")
            write_table(arrest_contrast(sl), out / "arrest")
    write_table(pl.concat(slope_tables), out / "unit_slopes")
    write_table(pl.concat(fits), out / "series_fits")
    write_table(gene_bootstrap(units, table, n_boot // 2, ctx.rng(310_002)), out / "shield")
    for name, frame in timing_shield(every, n_boot // 2, n_boot, ctx.rng(310_003)).items():
        write_table(frame, out / name)
    write_table(context_loss(every, n_boot, ctx.rng(310_004)), out / "context_loss")
    write_table(
        universe.group_by(keys := ["context", "timing_class", "domain", "hk"]).len().sort(keys),
        out / "timing_set_sizes",
    )
    print(f"solo_wcgw_culture: {len(units)} cultures; HypoClock CpGs listed {len(hypo)}", file=sys.stderr)
    return []
