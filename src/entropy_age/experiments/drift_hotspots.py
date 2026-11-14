"""Stochastic or specific: where does the spread that rises with age sit? (roadmap RM34 item 4)

The omics re-analysis practice: "stochastic" drift is a claim that needs a second test to fail. Take the loci whose
spread rises most with age, then ask whether they are enriched in particular region classes against a background
matched for mean level. Our own tests start from gene sets; this one starts from the loci.

On every autosomal CpG of Hannum blood and Jaffe's adult cortex (``figures tissues``' samples and covariates), each
CpG gets its age slope of dispersion and of relative entropy from the young state (``drift_metrics``). The top 5 %
on each are the hotspots. For every region class (island relation, promoter, solo-WCGW, the tissue's own 15-state
chromatin segmentation, the essential gene sets) the observed number of hotspots is set against the number expected
if hotspots were drawn at random within young-level deciles (the mean-matched background): the ratio, and a z from
the within-decile binomial variance. Broad drift would put every class near 1; a programme would make a few classes
stand far out.

``entropy-age figures drift_hotspots`` as a Slurm job.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.chromatin import ChromatinLayout, cpg_states, read_segmentation
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import SETS, _hannum, _jaffe, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

TOP = 0.05
METRICS = ("dispersion_slope", "kl_slope")
OWN_STATE = {"blood (Hannum)": "E062", "cortex (Jaffe)": "E073"}


def enrichment(hot: np.ndarray, member: np.ndarray, strata: np.ndarray) -> tuple[float, float, int, float]:
    """Observed hotspots in a class against those expected within strata: ratio, z, observed, expected."""
    obs, exp, var = 0.0, 0.0, 0.0
    for s in np.unique(strata):
        k = strata == s
        n_hot, share = hot[k].sum(), member[k].mean()
        obs += (hot & member & k).sum()
        exp += n_hot * share
        var += n_hot * share * (1 - share)
    ratio = obs / exp if exp > 0 else float("nan")
    z = (obs - exp) / np.sqrt(var) if var > 0 else float("nan")
    return float(ratio), float(z), int(obs), float(exp)


def classes(features: list[str], probes: pl.DataFrame, membership: pl.DataFrame, states: np.ndarray) -> dict:
    """Region class name -> boolean mask over ``features``."""
    f = pl.DataFrame({"cpg": features}).join(probes, on="cpg", how="left")
    out = {}
    for rel in f["island_relation"].drop_nulls().unique().sort().to_list():
        out[f"island relation: {rel}"] = (f["island_relation"] == rel).fill_null(False).to_numpy()
    out["promoter"] = f["promoter"].fill_null(False).to_numpy()
    if "solo_wcgw" in f.columns:
        out["solo-WCGW"] = f["solo_wcgw"].fill_null(False).to_numpy()
    for st in np.unique(states):
        out[f"chromatin: {st}"] = states == st
    for g in SETS:
        s = set(membership.filter(pl.col("group") == g)["feature"].to_list())
        out[f"set: {g}"] = np.array([c in s for c in features])
    return out


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("drift_hotspots")
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(~pl.col("chr").is_in(["X", "Y"]))
    membership = link_membership(pl.read_parquet(ann.links_450k), pl.read_parquet(ann.gene_sets))
    membership = membership.filter(pl.col("group").is_in(SETS))
    features = sorted(probes["cpg"].to_list())
    layout = ChromatinLayout(ctx.paths)
    rows = []
    for loader in (_hannum, _jaffe):
        p = loader(ctx.paths, features)
        p = dataclasses.replace(p, x=p.x.astype(np.float32))
        metrics = _lifespan_metrics(p)
        frame = pl.DataFrame({"cpg": p.features}).join(probes.select("cpg", "chr", "pos"), on="cpg", how="left")
        seg = read_segmentation(layout.segmentation(OWN_STATE[p.label]))
        states = cpg_states(frame, seg).fill_null("none").to_numpy().astype(str)
        masks = classes(p.features, probes, membership, states)
        strata = baseline_strata(metrics["baseline"].to_numpy(), 10)
        for metric in METRICS:
            v = metrics[metric].to_numpy()
            ok = np.isfinite(v)
            hot = np.zeros(v.size, dtype=bool)
            hot[ok] = v[ok] >= np.quantile(v[ok], 1 - TOP)
            for name, mask in masks.items():
                if mask[ok].sum() < 50:
                    continue
                ratio, z, obs, exp = enrichment(hot[ok], mask[ok], strata[ok])
                rows.append(
                    {
                        "tissue": p.label,
                        "metric": metric,
                        "class": name,
                        "cpgs": int(mask[ok].sum()),
                        "hotspots": obs,
                        "expected": exp,
                        "ratio": ratio,
                        "z": z,
                    }
                )
        print(f"drift_hotspots: {p.label}: {len(p.features)} CpGs", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "enrichment")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    t = table.filter(pl.col("metric") == "dispersion_slope")
    tissues = t["tissue"].unique(maintain_order=True).to_list()
    names = t["class"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(len(tissues), ncols=len(tissues), size=max(4.4, 0.17 * len(names)))
    for ax, tissue, letter in zip(axes, tissues, "ab", strict=False):
        d = pl.DataFrame({"class": names}).join(t.filter(pl.col("tissue") == tissue), on="class", how="left")
        y = np.arange(len(names))
        r = d["ratio"].fill_null(np.nan).to_numpy()
        colors = [plots.PALETTE[0] if v >= 1 else plots.PALETTE[1] for v in r]
        ax.barh(y, np.log2(r), color=colors)
        ax.set_yticks(y, names if letter == "a" else [""] * len(names), fontsize=6)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.set(xlabel="log2 observed / expected (level-matched)", title=f"{tissue}: top 5 % rising spread")
        plots.label_panel(ax, letter)
    fig.suptitle("Where the spread that rises with age sits", x=0.02, y=1.02, ha="left", fontsize=11)
    return plots.save(fig, ctx.figure_path("drift_hotspots"))
