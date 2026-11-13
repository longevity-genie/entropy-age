"""One cycle of in-vivo partial reprogramming, RRBS: does a reset restore the shield? (roadmap RM17)

Chondronasiou et al. 2022 (GSE156557, Aging Cell) gave 1-year-old female OSKM mice one doxycycline cycle and
two weeks off, against untreated mice of the same age and young controls (13 weeks), about 5 per arm, in
pancreas, liver and spleen (RRBS, Bismark coverage files). Chromatin maintenance predicts the shield restored
along with the clock (a reset of the epigenome); selection predicts nothing (no cells removed). R32 found no
consistent change on the mammalian array with cyclic OSKM (Browder 2022); this is the RRBS replication.

Per tissue: promoter CpGs (RefSeq TSS −1500 to +500, mm10) covered ``MIN_COVERAGE`` times in every sample of the
tissue; each old mouse scored against the young controls (``entropy_age.protection``; each young control against
the other young ones), drift ``B``, the shield share ``F`` and the gap ``P``; the effects of age (old controls
against young) and of OSKM (old OSKM against old controls), with a label-permutation p for 5 against 5.

``entropy-age figures reprogramming_rrbs`` as a Slurm job, after GSE156557's RAW.tar is unpacked into
``data/input/geo/GSE156557/cov``.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.rrbs import TRACKS, RrbsLayout, in_intervals, islands, promoters
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.cr_blood_replication import permutation_p
from entropy_age.experiments.interventions import _sets, membership
from entropy_age.experiments.per_person import score_against
from entropy_age.protection import fit_reference

GSE = "GSE156557"
MIN_COVERAGE = 10
ARMS = {"young": "Young control", "old": "Old control", "oskm": "Old OSKM"}
N_PERM = {"paper": 5000, "quick": 500}


def read_arms(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM), ``tissue`` and ``arm`` (young, old, oskm) from the series matrix's cell-type field."""
    rows = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            if cells[0] == "!Sample_geo_accession":
                rows["sample"] = cells[1:]
            elif cells[0] == "!Sample_source_name_ch1":
                rows["tissue"] = cells[1:]
            elif cells[0] == "!Sample_characteristics_ch1" and cells[1].startswith("cell type:"):
                rows["cell"] = [c.split(":", 1)[1].strip() for c in cells[1:]]
    frame = pl.DataFrame(rows)
    arm = pl.lit(None, dtype=pl.String)
    for key, prefix in ARMS.items():
        arm = pl.when(pl.col("cell").str.starts_with(prefix)).then(pl.lit(key)).otherwise(arm)
    return frame.with_columns(arm=arm).drop("cell").drop_nulls("arm")


def read_cov(path: Path) -> pl.DataFrame:
    """``cpg`` (``chr:pos``, 1-based, no ``chr`` prefix), ``beta`` and ``coverage`` of a Bismark coverage file."""
    frame = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chr", "start", "end", "percent", "m", "u"],
        schema_overrides={"chr": pl.String},
    )
    return frame.select(
        pl.format("{}:{}", pl.col("chr").str.strip_prefix("chr"), pl.col("start")).alias("cpg"),
        (pl.col("percent") / 100).cast(pl.Float32).alias("beta"),
        (pl.col("m") + pl.col("u")).cast(pl.Int32).alias("coverage"),
    )


def effect(y: np.ndarray, treated: np.ndarray) -> tuple[float, float]:
    return float(y[treated].mean() - y[~treated].mean()), float(y[~treated].mean())


def run(ctx: Context) -> list[Path]:
    base = ctx.paths.input / "geo" / GSE
    covdir = base / "cov"
    if not covdir.exists():
        print(f"reprogramming_rrbs: {covdir} missing; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("reprogramming_rrbs")
    lay = RrbsLayout(ctx.paths)
    prom = promoters(lay.track(TRACKS[0]), lay.track(TRACKS[1]))
    isl_tracks = islands(lay.track(TRACKS[2]))
    arms = read_arms(base / f"{GSE}_series_matrix.txt.gz")
    files = {re.match(r"(GSM\d+)", f.name).group(1): f for f in covdir.glob("GSM*cov*") if re.match(r"GSM\d+", f.name)}
    arms = arms.filter(pl.col("sample").is_in(list(files)))
    rng = ctx.rng(280_000)
    rows, animals = [], []
    for (tissue,), sel in arms.group_by("tissue", maintain_order=True):
        long = []
        for s in sel["sample"].to_list():
            t = read_cov(files[s]).filter(pl.col("coverage") >= MIN_COVERAGE)
            chrom = t["cpg"].str.split(":").list.get(0).to_numpy()
            pos = t["cpg"].str.split(":").list.get(1).cast(pl.Int64).to_numpy() - 1
            long.append(
                t.filter(pl.Series(in_intervals(chrom, pos, prom))).select("cpg", "beta").with_columns(sample=pl.lit(s))
            )
        long = pl.concat(long)
        kept = long.group_by("cpg").len().filter(pl.col("len") == sel.height).select("cpg")
        wide = long.join(kept, on="cpg", how="semi").pivot(on="sample", index="cpg", values="beta").sort("cpg")
        features = wide["cpg"].to_list()
        chrom = wide["cpg"].str.split(":").list.get(0)
        pos = wide["cpg"].str.split(":").list.get(1).cast(pl.Int64)
        cpgs = pl.DataFrame({"cpg": wide["cpg"], "chr": chrom, "pos": pos})
        isl = in_intervals(chrom.to_numpy(), pos.to_numpy() - 1, isl_tracks)
        sets = _sets(features, membership(ctx.paths, features, cpgs))
        member = np.any([s for s, _ in sets.values()], axis=0)
        order = sel["sample"].to_list()
        x = wide.select(order).to_numpy().T.astype(float)
        arm = sel["arm"].to_numpy()
        young = np.flatnonzero(arm == "young")
        scores = {}
        ref = fit_reference(x[young], isl, ~member)
        others = np.flatnonzero(arm != "young")
        for i, sc in zip(others, _rows(score_against(x[others], ref, sets)), strict=True):
            scores[int(i)] = sc
        for i in young:
            rest = young[young != i]
            sc = score_against(x[i : i + 1], fit_reference(x[rest], isl, ~member), sets)
            scores[int(i)] = {k: float(v[0]) for k, v in sc.items()}
        table = pl.DataFrame([{"sample": order[i], "arm": arm[i], **scores[i]} for i in range(len(order))])
        animals.append(table.with_columns(tissue=pl.lit(tissue)))
        print(f"reprogramming_rrbs: {tissue}: {len(features)} CpGs, {sel.height} mice", file=sys.stderr, flush=True)
        for contrast, (a, b) in {"age (old − young)": ("old", "young"), "OSKM (OSKM − old)": ("oskm", "old")}.items():
            sub = table.filter(pl.col("arm").is_in([a, b]))
            treated = (sub["arm"] == a).to_numpy()
            for group in sets:
                for measure in ("B", "F", "P"):
                    y = sub[f"{measure}_{group}"].to_numpy()
                    eff, ctrl = effect(y, treated)
                    rows.append(
                        {
                            "tissue": tissue,
                            "contrast": contrast,
                            "group": group,
                            "measure": measure,
                            "effect": eff,
                            "comparison_mean": ctrl,
                            "p_perm": permutation_p(y, treated, N_PERM[str(ctx.scale.name)], rng),
                            "n": f"{int(treated.sum())} vs {int((~treated).sum())}",
                            "cpgs": len(features),
                        }
                    )
    effects = pl.DataFrame(rows)
    write_table(effects, out / "effects")
    write_table(pl.concat(animals, how="diagonal"), out / "animal_scores")
    return [_figure(ctx, pl.concat(animals, how="diagonal"))]


def _rows(scores: dict[str, np.ndarray]) -> list[dict[str, float]]:
    n = len(next(iter(scores.values())))
    return [{k: float(v[i]) for k, v in scores.items()} for i in range(n)]


def _figure(ctx: Context, animals: pl.DataFrame) -> Path:
    tissues = animals["tissue"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for ax, measure, letter, label in zip(axes, ("F", "B"), "ab", ("shield share F", "drift B"), strict=True):
        col = f"{measure}_housekeeping"
        for j, arm in enumerate(ARMS):
            for k, tissue in enumerate(tissues):
                v = animals.filter((pl.col("tissue") == tissue) & (pl.col("arm") == arm))[col].to_numpy()
                xs = np.full(v.size, k + 0.25 * (j - 1))
                ax.scatter(xs, v, s=14, color=plots.PALETTE[j], label=ARMS[arm] if k == 0 else None)
                if v.size:
                    ax.plot([k + 0.25 * (j - 1) - 0.08, k + 0.25 * (j - 1) + 0.08], [v.mean()] * 2, color="k", lw=1.2)
        ax.set_xticks(range(len(tissues)), tissues)
        ax.set(ylabel=f"housekeeping {label}", title=label)
        if letter == "a":
            ax.legend(fontsize=7)
        plots.label_panel(ax, letter)
    fig.suptitle("One OSKM cycle in 1-year-old mice (GSE156557, RRBS)", x=0.02, y=1.03, ha="left", fontsize=11)
    fig.subplots_adjust(wspace=0.35)
    return plots.save(fig, ctx.figure_path("reprogramming_rrbs"))
