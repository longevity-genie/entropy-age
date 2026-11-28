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

Not controlled: replication timing (housekeeping promoters replicate early, and solo-WCGW loss is a
late-replication signal), which is part of what the shield could be here.

``entropy-age figures solo_wcgw_culture`` (needs the SeSAMe rebuilds of the series below).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout, link_membership
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


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    out = ctx.out_dir("solo_wcgw_culture")
    n_boot = N_BOOT[str(ctx.scale.name)]
    table, hypo = cpg_sets(paths)
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
    units = prepare(paths, samples, table)
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
    print(f"solo_wcgw_culture: {len(units)} cultures; HypoClock CpGs listed {len(hypo)}", file=sys.stderr)
    return []
