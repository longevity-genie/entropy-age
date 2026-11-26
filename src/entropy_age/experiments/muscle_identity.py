"""Does ageing skeletal muscle lose its cell identity in its methylome? Drift toward adipocytes, at identity CpGs.

References (all rebuilt from IDATs with the project's SeSAMe recipe, arrays passing its QC):

* **myofibre**: pooled type I and type II fibres of GSE304045 at baseline (EPIC v2), the mean of the two fibre-type
  means; whole muscle of the same people is kept as a bulk check;
* **myoblast**: primary myoblasts of GSE221540 (EPIC, older adults);
* **adipocyte**: isolated subcutaneous adipocytes of the lean controls of GSE222595 (450k and EPIC parts, the mean of
  the two platform means);
* **same-body tissues**: GSE48472 (450k), muscle and subcutaneous fat of the same autopsy bodies.

CpGs are the autosomal CpGs of the 450k annotation. Two kinds of identity CpGs:

(a) **identity genes**: CpGs linked to curated muscle identity genes (``IDENTITY_GENES``: myosins, actin, creatine
    kinase, desmin, the myogenic factors, MEF2C and classic MEF2 targets, sarcomere and excitation-contraction genes)
    at their promoter (TSS1500, TSS200, 5'UTR, first exon), or anywhere in the gene when the CpG lies in an enhancer
    state (6_EnhG, 7_Enh) of either Roadmap skeletal-muscle epigenome (E107, E108);
(b) **muscle-vs-adipocyte CpGs** (data-driven): ``|myofibre - adipocyte| > 0.3`` in the sorted references, confirmed
    by the same sign and at least 0.1 in the paired muscle-minus-fat difference of GSE48472.

The **adipocyte axis** is ``w = adipocyte - myofibre`` over the (b) CpGs; a person's identity distance is
``t = sum w (x - myofibre) / sum w^2`` (0 at the myofibre reference, 1 at the adipocyte reference). The same is done
toward the myoblast reference (``|myofibre - myoblast| > 0.3``), along the type I to type II fibre axis (``> 0.2``; a
fibre-type shift with age would move every axis), and toward adipocytes over the (a) CpGs. Because CpGs regress
toward intermediate levels with age in any tissue, each axis also gets a **matched null**: every axis CpG replaced by
background CpGs of the same island status and cohort mean-level stratum (20 bins), the weights kept, averaged over
``N_NULL_DRAWS`` draws; ``t - t_null`` is the shift beyond generic drift.

Per cohort (never pooled): ``t`` regressed on age with the cohort's held covariates (sex where it varies, the
Hertfordshire cohort label, which is also the extraction method), and again with chip intercepts where chips are known
and age stays identifiable; intervals from resampling donors. Drift at identity CpGs: the per-CpG age slopes of
relative entropy, dispersion and entropy (``stochasticity.drift_metrics``) compared with background CpGs matched on
island status and young-level decile (``protection.matched_scores``, as ``series.run_series``), scaled by the SD of the
metric, with donor-resampling intervals (``series.donor_intervals``). In the obese cohort (GSE61454 muscle) ``t`` is
regressed on BMI with age and sex held. GSE61454 muscle passes the 90 % detection cut in only 13 of 60 arrays, so it
is also run, labelled as a sensitivity analysis, at 80 % detection (only CpGs measured in every array are used).
GTEx muscle (GSE213478) enters with its deposited BMIQ betas and decade-midpoint ages, labelled as such.

Limits: no cohort spans young to old with IDATs except 20 post-mortem donors; Hertfordshire is 68-83 only; references
come from other labs and platforms, so ``t`` is comparable within a cohort only; bulk muscle holds fibroblasts,
endothelium and blood, so its ``t`` is above 0 at any age; no sorted FAP reference exists.
Run as ``muscle_identity.run(ctx)``; outputs in ``data/output/muscle_identity``.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.chromatin import URL, ChromatinLayout, cpg_states, read_segmentation
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.data.download import download
from entropy_age.experiments.brain_replication import BANDS, METRICS, MIN_SAMPLES
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.series import (
    N_DONOR_BOOT,
    age_identifiable,
    donor_intervals,
    groups_of,
    selection,
    series_table,
    specs,
    titles,
)
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, drift_metrics

IDENTITY_GENES = (
    # myosins, actin, creatine kinase, desmin, myoglobin
    "MYH1", "MYH2", "MYH4", "MYH7", "MYL1", "MYL2", "MYLPF", "ACTA1", "CKM", "DES", "MB",
    # myogenic regulatory factors and MEF2C
    "MYOD1", "MYOG", "MYF5", "MYF6", "MEF2C",
    # sarcomere and MEF2 targets
    "TTN", "NEB", "TNNT1", "TNNT3", "TNNI1", "TNNI2", "TNNC1", "TNNC2", "ACTN2", "ACTN3", "MYBPC1", "MYBPC2",
    "MYOZ1", "MYOT", "LMOD2", "SMPX", "CAV3", "DMD",
    # excitation-contraction coupling and muscle metabolism
    "RYR1", "ATP2A1", "CASQ1", "CHRNA1", "PYGM",
)  # fmt: skip
PROMOTER_REGIONS = ("TSS1500", "TSS200", "5'UTR", "1stExon")
ENHANCER_STATES = ("6_EnhG", "7_Enh")
MUSCLE_EPIGENOMES = ("E107", "E108")  # Roadmap skeletal muscle, male and female
DIFF = 0.3  # muscle vs adipocyte (and vs myoblast) difference in the references
CONFIRM = 0.1  # the paired GSE48472 muscle-minus-fat difference must agree in sign and reach this
FIBRE_TYPE_DIFF = 0.2
BACKGROUND = 40_000
N_NULL_DRAWS = 20
LEVEL_BINS = 20
N_SLOPE_BOOT = {"paper": 2000, "quick": 100}
RELAXED_DETECTION = 0.80

#: (label, config name or None, group or spec) of every muscle cohort analysed per person; never pooled
GENE_SMART = {
    "accession": "GSE171140",
    "age": "age",
    "sex": "Sex",
    "donor": "title:_(GS\\d+)_",
    "derive": {"timepoint": "title:^Muscle_([^_]+)_"},
    "groups": {"skeletal muscle (Gene SMART, baseline)": {"timepoint": "PRE"}},
}
COHORTS = (
    ("muscle_hertfordshire_epic", "skeletal muscle (Hertfordshire, EPIC)"),
    ("muscle_hertfordshire_450k", "skeletal muscle (Hertfordshire, 450k)"),
    ("muscle_hertfordshire_450k", "cultured myoblasts (Hertfordshire, 450k)"),
    ("muscle_postmortem_korean", "skeletal muscle (post-mortem, GSE244995)"),
    ("obese_four_tissues", "skeletal_muscle"),
    ("gene_smart", "skeletal muscle (Gene SMART, baseline)"),
)
HOLD_OF = {"muscle_hertfordshire_epic": ["cohort"]}


# ---------------------------------------------------------------- pure helpers


def axis_projection(x: np.ndarray, origin: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Per row of ``x`` (samples x CpGs): its position on the line from ``origin`` (0) to ``target`` (1),
    ``sum w (x - origin) / sum w^2`` with ``w = target - origin``, over the CpGs finite in the row and both ends."""
    x = np.atleast_2d(np.asarray(x, dtype=float))
    w = np.asarray(target, dtype=float) - np.asarray(origin, dtype=float)
    d = x - origin
    ok = np.isfinite(d) & np.isfinite(w)
    num = np.where(ok, d * w, 0.0).sum(axis=1)
    den = np.where(ok, w * w, 0.0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def differential(
    a: np.ndarray, b: np.ndarray, threshold: float, confirm: np.ndarray | None = None, confirm_min: float = 0.0
) -> np.ndarray:
    """CpGs where ``|a - b| > threshold``; with ``confirm`` (another estimate of ``a - b``) also of the same sign and
    at least ``confirm_min`` in it. Missing values never pass."""
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    keep = np.isfinite(diff) & (np.abs(diff) > threshold)
    if confirm is not None:
        c = np.asarray(confirm, dtype=float)
        keep &= np.isfinite(c) & (np.sign(c) == np.sign(diff)) & (np.abs(c) >= confirm_min)
    return keep


def identity_gene_cpgs(links: pl.DataFrame, enhancer_cpgs: set[str], genes=IDENTITY_GENES) -> pl.DataFrame:
    """``cpg``, ``symbol``, ``kind`` (promoter or enhancer) of CpGs linked to ``genes`` at their promoter, or anywhere
    in the gene when the CpG is in ``enhancer_cpgs``; a CpG at a promoter counts as promoter."""
    linked = links.filter(pl.col("symbol").is_in(list(genes)))
    promoter = linked.filter(pl.col("region").is_in(list(PROMOTER_REGIONS))).with_columns(kind=pl.lit("promoter"))
    enhancer = linked.filter(pl.col("cpg").is_in(list(enhancer_cpgs))).with_columns(kind=pl.lit("enhancer"))
    both = pl.concat([promoter.select("cpg", "symbol", "kind"), enhancer.select("cpg", "symbol", "kind")])
    return both.sort("kind", descending=True).unique(subset="cpg", keep="first").sort("cpg")


def ols_slope(y: np.ndarray, predictor: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Coefficient of ``predictor`` and its SE in ``y ~ 1 + predictor + cov``."""
    design = (
        np.column_stack([np.ones(y.size), predictor, cov])
        if cov.shape[1]
        else np.column_stack([np.ones(y.size), predictor])
    )
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    df = y.size - np.linalg.matrix_rank(design)
    if df <= 0:
        return float(coef[1]), float("nan")
    sigma2 = float(resid @ resid) / df
    xtx = np.linalg.pinv(design.T @ design)
    return float(coef[1]), float(math.sqrt(max(sigma2 * xtx[1, 1], 0.0)))


def slope_interval(y, predictor, cov, n: int, rng: np.random.Generator) -> dict:
    """OLS slope of ``predictor`` with ``cov`` held, its SE and normal p, and the 2.5/97.5 percentiles over ``n``
    resamples of donors (rows); resamples where the predictor is not identifiable beside ``cov`` are skipped, and with
    fewer than half usable the interval is NaN."""
    y, predictor = np.asarray(y, float), np.asarray(predictor, float)
    cov = np.asarray(cov, float).reshape(y.size, -1)
    cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
    if not age_identifiable(predictor, cov, min_df=5):
        return {
            "slope": float("nan"),
            "se": float("nan"),
            "p": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "resamples": 0,
        }
    slope, se = ols_slope(y, predictor, cov)
    draws = []
    for _ in range(n):
        rows = rng.integers(0, y.size, y.size)
        c = cov[rows]
        c = c[:, np.ptp(c, axis=0) > 1e-12] if c.shape[1] else c
        if not age_identifiable(predictor[rows], c, distinct=np.unique(rows).size, min_df=5):
            continue
        draws.append(ols_slope(y[rows], predictor[rows], c)[0])
    ok = len(draws) >= n / 2
    p = math.erfc(abs(slope / se) / math.sqrt(2)) if se and np.isfinite(se) and se > 0 else float("nan")
    return {
        "slope": slope,
        "se": se,
        "p": p,
        "ci_low": float(np.quantile(draws, 0.025)) if ok else float("nan"),
        "ci_high": float(np.quantile(draws, 0.975)) if ok else float("nan"),
        "resamples": len(draws),
    }


def strata_of(island: np.ndarray, level: np.ndarray, bins: int) -> np.ndarray:
    """Island status x quantile bin of ``level``, as one label per CpG."""
    return np.char.add(np.char.add(np.asarray(island).astype(str), "|"), baseline_strata(level, bins).astype(str))


def matched_partners(strata: np.ndarray, in_axis: np.ndarray, background: np.ndarray, draws: int, rng) -> np.ndarray:
    """``(draws, n_axis)`` indices of background CpGs drawn from each axis CpG's stratum (with replacement); -1 where
    the stratum holds no background CpG."""
    axis_idx = np.flatnonzero(in_axis)
    pool: dict[str, np.ndarray] = {}
    for s in np.unique(strata[axis_idx]):
        pool[s] = np.flatnonzero(background & (strata == s))
    out = np.full((draws, axis_idx.size), -1, dtype=int)
    for j, i in enumerate(axis_idx):
        cand = pool[strata[i]]
        if cand.size:
            out[:, j] = rng.choice(cand, size=draws, replace=True)
    return out


def null_projection(
    x: np.ndarray, origin: np.ndarray, target: np.ndarray, partners: np.ndarray, axis_idx
) -> np.ndarray:
    """Per row: the axis projection with each axis CpG's deviation ``x - origin`` taken at a matched background CpG
    (``partners``, indices into the columns of ``x``; its own origin level there), the weights ``target - origin`` of
    the axis CpG kept; the mean over draws. Pairs without a partner are left out."""
    w = (target - origin)[axis_idx]
    out = []
    for draw in partners:
        ok = draw >= 0
        d = x[:, draw[ok]] - origin[draw[ok]]
        out.append((d * w[ok]).sum(axis=1) / (w[ok] ** 2).sum())
    return np.mean(out, axis=0)


def band_gaps(metrics: pl.DataFrame, island: np.ndarray, sets: dict[str, np.ndarray]) -> list[dict]:
    """Per band, metric and set: the set-minus-matched gap (island x young-level decile strata within the band),
    divided by the metric's SD over the band, as ``series.donor_intervals`` scales its resamples."""
    baseline = metrics["baseline"].to_numpy()
    rows = []
    for band, (lo, hi) in BANDS.items():
        keep = (baseline >= lo) & (baseline <= hi)
        strata = strata_of(island[keep], baseline[keep], 10)
        for metric in METRICS:
            values = metrics[metric].to_numpy()[keep]
            sd = float(np.std(values))
            for name, in_set in sets.items():
                n_in = int(in_set[keep].sum())
                if n_in < 2:
                    continue
                gap = float(matched_scores(values[None, :], in_set[keep].astype(float), strata)[0][0]) / sd
                rows.append({"band": band, "metric": metric, "set": name, "cpgs_in_set": n_in, "gap_sd": gap})
    return rows


# ---------------------------------------------------------------- data


def _universe(paths) -> pl.DataFrame:
    probes = pl.read_parquet(AnnotationLayout(paths).probes_450k)
    return probes.filter(~pl.col("chr").is_in(["X", "Y"])).select("cpg", "chr", "pos", "cpg_island").sort("cpg")


def _enhancer_cpgs(paths, probes: pl.DataFrame) -> set[str]:
    lay = ChromatinLayout(paths)
    out: set[str] = set()
    for e in MUSCLE_EPIGENOMES:
        path = lay.segmentation(e)
        if not path.exists():
            download(URL.format(epigenome=e), path, timeout=600)
        states = cpg_states(probes, read_segmentation(path))
        out |= set(probes.filter(states.is_in(list(ENHANCER_STATES)).fill_null(False))["cpg"].to_list())
    return out


class Betas:
    """Rebuilt betas of each series, all arrays, restricted to the universe CpGs; read once."""

    def __init__(self, paths, cpgs: list[str]):
        self.paths, self.cpgs, self.frames = paths, cpgs, {}

    def frame(self, gse: str) -> pl.DataFrame:
        if gse not in self.frames:
            f = reprocessed.read_betas(self.paths, gse, passed_only=False).rename({"probe": "cpg"})
            self.frames[gse] = f.filter(pl.col("cpg").is_in(self.cpgs))
        return self.frames[gse]

    def matrix(self, gse: str, samples: list[str]) -> np.ndarray:
        return reprocessed.matrix(self.frame(gse), self.cpgs, samples)


def _profile(x: np.ndarray, min_share: float = 0.5) -> np.ndarray:
    measured = np.isfinite(x).mean(axis=0) >= min_share
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m = np.nanmean(x, axis=0)
    return np.where(measured, m, np.nan)


def _passed(paths, gse: str) -> set[str]:
    return set(reprocessed.read_qc(paths, gse).filter(pl.col("passed"))["sample"].to_list())


def references(paths, betas: Betas) -> tuple[dict[str, np.ndarray], list[dict]]:
    """Reference profiles over the universe CpGs and a row per reference saying what went in."""
    geo = paths.input / "geo"
    out, rows = {}, []
    # myofibres and whole muscle, GSE304045 baseline
    m = geo / "GSE304045" / "GSE304045_series_matrix.txt.gz"
    chars = _characteristics(m)
    ok = _passed(paths, "GSE304045")
    base = chars.filter((pl.col("timepoint") == "Pre") & (pl.col("replicate") == "No") & pl.col("sample").is_in(ok))
    types = {}
    for t in ("T1", "T2", "WM"):
        s = base.filter(pl.col("cell type") == t)["sample"].to_list()
        types[t] = _profile(betas.matrix("GSE304045", s))
        rows.append({"reference": f"GSE304045 {t} baseline", "arrays": len(s)})
    out["type I fibre"], out["type II fibre"] = types["T1"], types["T2"]
    out["myofibre"] = np.nanmean(np.vstack([types["T1"], types["T2"]]), axis=0)
    out["whole muscle (GSE304045)"] = types["WM"]
    # myoblasts
    s = sorted(_passed(paths, "GSE221540"))
    out["myoblast"] = _profile(betas.matrix("GSE221540", s))
    rows.append({"reference": "GSE221540 myoblasts", "arrays": len(s)})
    # adipocytes: lean controls, subcutaneous, mean of platform means
    ok = _passed(paths, "GSE222595")
    per_platform = []
    for gpl in ("GPL13534", "GPL21145"):
        mm = geo / "GSE222595" / f"GSE222595-{gpl}_series_matrix.txt.gz"
        frame = _characteristics(mm).with_columns(
            title=pl.Series(titles(mm)), source=pl.Series(titles(mm, "!Sample_source_name_ch1"))
        )
        s = frame.filter(
            pl.col("title").str.contains("control") & (pl.col("source") == "subcutaneous") & pl.col("sample").is_in(ok)
        )["sample"].to_list()
        per_platform.append(_profile(betas.matrix("GSE222595", s)))
        rows.append({"reference": f"GSE222595 {gpl} subcutaneous adipocytes, controls", "arrays": len(s)})
    out["adipocyte"] = np.nanmean(np.vstack(per_platform), axis=0)
    # GSE48472 same bodies
    mm = geo / "GSE48472" / "GSE48472_series_matrix.txt.gz"
    ok = _passed(paths, "GSE48472")
    frame = pl.DataFrame({"sample": _characteristics(mm)["sample"], "title": titles(mm)}).filter(
        pl.col("sample").is_in(ok)
    )
    frame = frame.with_columns(
        tissue=pl.col("title").str.replace(r"_(IT|PT)\d+$", ""), donor=pl.col("title").str.extract(r"_((?:IT|PT)\d+)$")
    )
    muscle = frame.filter(pl.col("tissue") == "Muscle")
    fat = frame.filter(pl.col("tissue") == "SC_Fat")
    pairs = muscle.join(fat, on="donor", suffix="_fat")
    xm = betas.matrix("GSE48472", pairs["sample"].to_list())
    xf = betas.matrix("GSE48472", pairs["sample_fat"].to_list())
    out["muscle minus fat (GSE48472)"] = _profile(xm - xf)
    out["muscle (GSE48472)"], out["subcutaneous fat (GSE48472)"] = _profile(xm), _profile(xf)
    om = frame.filter(pl.col("tissue") == "Omentum")["sample"].to_list()
    out["omentum (GSE48472)"] = _profile(betas.matrix("GSE48472", om))
    rows.append({"reference": "GSE48472 muscle-fat pairs (same bodies)", "arrays": 2 * pairs.height})
    return out, rows


def _cohort(paths, name: str, group: str, relaxed: bool = False) -> tuple[pl.DataFrame, str] | None:
    spec = GENE_SMART if name == "gene_smart" else specs()[name]
    gse = spec["accession"]
    if relaxed:
        table = series_table(paths, spec)
        qc = pl.read_csv(reprocessed.series_dir(paths, gse) / f"{gse}.qc.csv", infer_schema_length=None)
        usable = qc.filter(pl.col("frac_dt").cast(pl.Float64) >= RELAXED_DETECTION)["sample"].unique().to_list()
        groups = groups_of(table, spec, usable, set())
    else:
        chosen = selection(paths, name, spec)
        if chosen is None:
            return None
        groups = chosen[0]
    return groups[group].sort("age", "sample"), gse


def _gtex(paths, cpgs: list[str]) -> tuple[pl.DataFrame, np.ndarray] | None:
    lay = cohort_layout("gtex", paths)
    if not lay.betas.exists():
        return None
    s = (
        pl.read_parquet(lay.samples)
        .filter(pl.col("tissue") == "skeletal muscle")
        .drop_nulls("age")
        .unique(subset="donor", keep="first", maintain_order=True)
        .sort("age", "sample")
    )
    frame = pl.scan_parquet(lay.betas).select("cpg", *s["sample"].to_list()).filter(pl.col("cpg").is_in(cpgs))
    return s, reprocessed.matrix(frame.collect(), cpgs, s["sample"].to_list())


# ---------------------------------------------------------------- run


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    needed = ["GSE304045", "GSE221540", "GSE222595", "GSE48472"]
    missing = [g for g in needed if not reprocessed.available(paths, g)]
    if missing:
        print(f"muscle_identity: no SeSAMe rebuild for references {missing}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("muscle_identity")
    scale = str(ctx.scale.name)
    probes = _universe(paths)
    cpgs = probes["cpg"].to_list()
    island = probes["cpg_island"].fill_null(False).to_numpy()
    betas = Betas(paths, cpgs)
    ref, ref_rows = references(paths, betas)
    print(f"muscle_identity: references read ({len(ref_rows)} groups)", file=sys.stderr, flush=True)

    # identity CpGs
    links = pl.read_parquet(AnnotationLayout(paths).links_450k)
    gene_table = identity_gene_cpgs(links.filter(pl.col("cpg").is_in(cpgs)), _enhancer_cpgs(paths, probes))
    gene_set = np.isin(np.array(cpgs), gene_table["cpg"].to_list())
    raw_adipo = differential(ref["myofibre"], ref["adipocyte"], DIFF)
    adipo = differential(ref["myofibre"], ref["adipocyte"], DIFF, ref["muscle minus fat (GSE48472)"], CONFIRM)
    myoblast = differential(ref["myofibre"], ref["myoblast"], DIFF)
    fibre_type = differential(ref["type II fibre"], ref["type I fibre"], FIBRE_TYPE_DIFF)
    up = adipo & (ref["adipocyte"] > ref["myofibre"])
    set_rows = [
        {
            "set": "identity genes (promoter + muscle enhancer)",
            "cpgs": int(gene_set.sum()),
            "promoter": int((gene_table["kind"] == "promoter").sum()),
            "enhancer": int((gene_table["kind"] == "enhancer").sum()),
            "genes": gene_table["symbol"].n_unique(),
        },
        {"set": "myofibre vs adipocyte > 0.3 (sorted references)", "cpgs": int(raw_adipo.sum())},
        {
            "set": "muscle-vs-adipocyte CpGs (confirmed in GSE48472)",
            "cpgs": int(adipo.sum()),
            "higher_in_adipocyte": int(up.sum()),
            "lower_in_adipocyte": int((adipo & ~up).sum()),
            "overlap_identity_genes": int((adipo & gene_set).sum()),
        },
        {"set": "myofibre vs myoblast > 0.3", "cpgs": int(myoblast.sum())},
        {"set": "type I vs type II fibre > 0.2", "cpgs": int(fibre_type.sum())},
    ]
    write_table(pl.DataFrame(set_rows, infer_schema_length=None), out / "cpg_sets")
    write_table(pl.DataFrame(ref_rows), out / "references")
    write_table(
        pl.DataFrame(
            {
                "cpg": cpgs,
                "identity_gene": gene_set,
                "adipocyte_axis": adipo,
                "myoblast_axis": myoblast,
                "fibre_type_axis": fibre_type,
                "myofibre": ref["myofibre"],
                "adipocyte": ref["adipocyte"],
                "myoblast": ref["myoblast"],
            }
        ).filter(
            pl.col("identity_gene") | pl.col("adipocyte_axis") | pl.col("myoblast_axis") | pl.col("fibre_type_axis")
        ),
        out / "axis_cpgs",
    )

    axes = {
        "adipocyte axis": (adipo, ref["myofibre"], ref["adipocyte"]),
        "adipocyte axis, identity genes": (gene_set, ref["myofibre"], ref["adipocyte"]),
        "myoblast axis": (myoblast, ref["myofibre"], ref["myoblast"]),
        "fibre type axis (I to II)": (fibre_type, ref["type I fibre"], ref["type II fibre"]),
    }
    # where the references and same-body tissues sit on each axis (sanity)
    pos_rows = []
    for rname in (
        "whole muscle (GSE304045)",
        "muscle (GSE48472)",
        "subcutaneous fat (GSE48472)",
        "omentum (GSE48472)",
        "myoblast",
        "adipocyte",
        "myofibre",
    ):
        for aname, (mask, o, t) in axes.items():
            pos_rows.append(
                {
                    "profile": rname,
                    "axis": aname,
                    "position": float(axis_projection(ref[rname][mask][None], o[mask], t[mask])[0]),
                }
            )
    write_table(pl.DataFrame(pos_rows), out / "reference_positions")

    rng_bg = ctx.rng(610_000)
    background_pool = np.flatnonzero(~(gene_set | adipo | myoblast | fibre_type))
    background = np.zeros(len(cpgs), dtype=bool)
    background[rng_bg.choice(background_pool, size=min(BACKGROUND, background_pool.size), replace=False)] = True

    cohorts: list[tuple[str, pl.DataFrame, np.ndarray, np.ndarray, np.ndarray, str]] = []
    for name, group in COHORTS:
        got = _cohort(paths, name, group)
        if got is None:
            print(f"muscle_identity: {name}: no rebuild or sample table; skipped", file=sys.stderr)
            continue
        s, gse = got
        x = betas.matrix(gse, s["sample"].to_list())
        qc = reprocessed.read_qc(paths, gse)
        chip_of = dict(qc.select("sample", "chip").iter_rows())
        chips = np.array([str(chip_of.get(v) or "NA") for v in s["sample"].to_list()])
        held = [dummies(s["sex"].fill_null("NA").to_numpy())]
        held += [dummies(s[c].fill_null("NA").to_numpy()) for c in HOLD_OF.get(name, [])]
        cohorts.append((group, s, x, np.column_stack(held), chips, "SeSAMe rebuild, 90 % detection"))
        if name == "obese_four_tissues":
            got = _cohort(paths, name, group, relaxed=True)
            s2, _ = got
            x2 = betas.matrix(gse, s2["sample"].to_list())
            chips2 = np.array([str(chip_of.get(v) or "NA") for v in s2["sample"].to_list()])
            held2 = dummies(s2["sex"].fill_null("NA").to_numpy())
            cohorts.append(
                (
                    f"{group} (80 % detection, sensitivity)",
                    s2,
                    x2,
                    held2,
                    chips2,
                    "SeSAMe rebuild, 80 % detection (sensitivity)",
                )
            )
    gtex = _gtex(paths, cpgs)
    if gtex is not None:
        s, x = gtex
        cohorts.append(
            (
                "skeletal muscle (GTEx, deposited BMIQ)",
                s,
                x,
                dummies(s["sex"].fill_null("NA").to_numpy()),
                s["plate"].fill_null("NA").to_numpy().astype(str),
                "deposited BMIQ, decade midpoints",
            )
        )

    obese_bmi = {}
    mm = paths.input / "geo" / "GSE61454" / "GSE61454_series_matrix.txt.gz"
    if mm.exists():
        chars = _characteristics(mm)
        obese_bmi = dict(chars.select("sample", pl.col("bmi").cast(pl.Float64, strict=False)).iter_rows())

    size_rows, person_rows, slope_rows, gap_rows, bmi_rows = [], [], [], [], []
    n_boot = N_SLOPE_BOOT[scale]
    for ci, (label, s, x, held, chips, source) in enumerate(cohorts):
        ages = s["age"].to_numpy().astype(float)
        complete = np.isfinite(x).all(axis=0)
        level = np.nanmean(x, axis=0)
        strata = strata_of(island, np.where(complete, level, np.nan), LEVEL_BINS)
        size_rows.append(
            {
                "cohort": label,
                "source": source,
                "donors": s.height,
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs_complete": int(complete.sum()),
                "chips_known": int(np.unique(chips[chips != "NA"]).size),
            }
        )
        designs = {"held": held}
        chip_cov = dummies(chips)
        if (chips != "NA").all() and chip_cov.shape[1]:
            designs["held + chip"] = np.column_stack([held, chip_cov])
        values: dict[str, np.ndarray] = {}
        for ai, (aname, (mask, o, t)) in enumerate(axes.items()):
            use = mask & complete & np.isfinite(o) & np.isfinite(t)
            tt = axis_projection(x[:, use], o[use], t[use])
            idx = np.flatnonzero(use)
            partners = matched_partners(
                strata, use, background & complete & np.isfinite(o), N_NULL_DRAWS, ctx.rng(620_000 + 37 * ci + ai)
            )
            tn = null_projection(x, o, t, partners, idx)
            values[aname] = tt
            values[f"{aname}, beyond matched"] = tt - tn
            for v_name, v in ((aname, tt), (f"{aname}, beyond matched", tt - tn)):
                for design, cov in designs.items():
                    r = slope_interval(v, ages, cov, n_boot, ctx.rng(630_000 + 101 * ci + 7 * ai + len(design)))
                    slope_rows.append(
                        {
                            "cohort": label,
                            "measure": v_name,
                            "covariates": design,
                            "axis_cpgs": int(use.sum()),
                            "mean_position": float(np.mean(v)),
                            "sd_position": float(np.std(v, ddof=1)),
                            **{
                                f"{k}_per_decade" if k in ("slope", "se", "ci_low", "ci_high") else k: (
                                    val * 10 if k in ("slope", "se", "ci_low", "ci_high") else val
                                )
                                for k, val in r.items()
                            },
                        }
                    )
        for i, row in enumerate(s.select("sample", "donor", "age").iter_rows(named=True)):
            person_rows.append({"cohort": label, **row, **{k: float(v[i]) for k, v in values.items()}})
        # BMI in the obese cohort
        if label.startswith("skeletal_muscle"):
            bmi = np.array([obese_bmi.get(v, np.nan) for v in s["sample"].to_list()], dtype=float)
            ok = np.isfinite(bmi)
            for k, v in values.items():
                cov = np.column_stack([ages[ok], held[ok]])
                r = slope_interval(v[ok], bmi[ok], cov, n_boot, ctx.rng(640_000 + ci))
                bmi_rows.append(
                    {
                        "cohort": label,
                        "measure": k,
                        "donors": int(ok.sum()),
                        "bmi_min": float(bmi[ok].min()),
                        "bmi_max": float(bmi[ok].max()),
                        **{
                            f"{kk}_per_5_bmi" if kk in ("slope", "se", "ci_low", "ci_high") else kk: (
                                val * 5 if kk in ("slope", "se", "ci_low", "ci_high") else val
                            )
                            for kk, val in r.items()
                        },
                    }
                )
        # drift at identity CpGs vs matched CpGs
        if s.height < MIN_SAMPLES:
            print(f"muscle_identity: {label}: {s.height} donors, below {MIN_SAMPLES}; no drift gaps", file=sys.stderr)
            continue
        feat = complete & (gene_set | adipo | background)
        xf = x[:, feat]
        isl = island[feat]
        sets = {
            "identity genes": gene_set[feat],
            "muscle-vs-adipocyte CpGs": adipo[feat],
            "adipocyte-higher CpGs": up[feat],
            "adipocyte-lower CpGs": (adipo & ~up)[feat],
        }
        names = [c for c, k in zip(cpgs, feat, strict=True) if k]
        for k, (design, cov) in enumerate(designs.items()):
            cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
            if not age_identifiable(ages, cov):
                print(f"muscle_identity: {label}: {design}: age not identifiable; no drift gaps", file=sys.stderr)
                continue
            full = drift_metrics(xf, ages, names, covariates=cov if cov.shape[1] else None, spearman=False)
            by_donor = donor_intervals(
                xf,
                ages,
                cov,
                names,
                isl,
                {n: (v.astype(float), []) for n, v in sets.items()},
                N_DONOR_BOOT[scale],
                ctx.rng(650_000 + 10 * ci + k),
            )
            for r in band_gaps(full, isl, sets):
                lo, hi, kept = by_donor.get((r["band"], r["metric"], r["set"]), (float("nan"), float("nan"), 0))
                gap_rows.append(
                    {
                        "cohort": label,
                        "covariates": design,
                        **r,
                        "donor_ci_low": lo,
                        "donor_ci_high": hi,
                        "donor_resamples": kept,
                    }
                )
        print(f"muscle_identity: {label}: {s.height} donors done", file=sys.stderr, flush=True)

    write_table(pl.DataFrame(size_rows), out / "sizes")
    write_table(pl.DataFrame(person_rows), out / "positions")
    write_table(pl.DataFrame(slope_rows), out / "age_slopes")
    if gap_rows:
        write_table(pl.DataFrame(gap_rows), out / "drift_gaps")
    if bmi_rows:
        write_table(pl.DataFrame(bmi_rows), out / "bmi")
    return [
        _figure(ctx, pl.DataFrame(person_rows), pl.DataFrame(slope_rows), pl.DataFrame(gap_rows) if gap_rows else None)
    ]


def _figure(ctx: Context, persons: pl.DataFrame, slopes: pl.DataFrame, gaps: pl.DataFrame | None) -> Path:
    plots.apply_style()
    labels = persons["cohort"].unique(maintain_order=True).to_list()
    fig, axs = plots.panels(len(labels) + 1, ncols=4, size=3.0)
    for ax, label in zip(axs, labels, strict=False):
        p = persons.filter(pl.col("cohort") == label)
        ax.scatter(p["age"], p["adipocyte axis"], s=12, color=plots.POINT, alpha=0.8)
        r = slopes.filter(
            (pl.col("cohort") == label) & (pl.col("measure") == "adipocyte axis") & (pl.col("covariates") == "held")
        ).row(0, named=True)
        a = p["age"].to_numpy()
        if np.isfinite(r["slope_per_decade"]):
            grid = np.linspace(a.min(), a.max(), 2)
            ax.plot(
                grid, p["adipocyte axis"].mean() + r["slope_per_decade"] / 10 * (grid - a.mean()), color=plots.ACCENT
            )
        ax.set_title(
            label.replace("skeletal muscle ", "muscle ").replace("skeletal_muscle", "obese muscle (GSE61452)"),
            fontsize=8,
        )
        ax.set_xlabel("Age (years)")
        ax.set_ylabel("Position, myofibre (0) to adipocyte (1)")
        ax.text(
            0.02,
            0.98,
            f"{r['slope_per_decade']:+.3f}/decade\n[{r['ci_low_per_decade']:+.3f}, {r['ci_high_per_decade']:+.3f}]",
            transform=ax.transAxes,
            va="top",
            fontsize=7,
        )
    ax = axs[len(labels)]
    if gaps is not None and gaps.height:
        g = gaps.filter(
            (pl.col("band") == "all")
            & (pl.col("metric") == "kl_slope")
            & (pl.col("covariates") == "held")
            & pl.col("set").is_in(["identity genes", "muscle-vs-adipocyte CpGs"])
        )
        for i, row in enumerate(g.iter_rows(named=True)):
            color = plots.PALETTE[0] if row["set"] == "identity genes" else plots.PALETTE[1]
            ax.errorbar(
                row["gap_sd"],
                i,
                xerr=[[row["gap_sd"] - row["donor_ci_low"]], [row["donor_ci_high"] - row["gap_sd"]]],
                fmt="o",
                color=color,
                ms=4,
            )
        ax.set_yticks(range(g.height))
        ax.set_yticklabels([f"{c[:22]} | {s[:10]}" for c, s in zip(g["cohort"], g["set"], strict=True)], fontsize=6)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_xlabel("Relative-entropy slope gap vs matched (SD)")
        ax.set_title("Drift at identity CpGs", fontsize=8)
    return plots.save(fig, ctx.figure_path("muscle_identity"))
