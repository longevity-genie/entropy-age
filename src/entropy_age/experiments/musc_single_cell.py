"""Drift between single cells: muscle stem cells of young and old mice, methylation and RNA from the same cell.

Every other analysis here reads drift between people or cultures. GSE121437 (Hernando-Herraez et al. 2019) sequenced
single muscle stem cells of 2-month and 24-month mice for methylation (scBS-seq, GSE121436; 275 cells of three young
and three old mice) and, from the same cells, RNA (GSE121364; counts deposited for a part of them). Muscle is the
tissue in which the human array data show almost no shield, and its stem cells divide rarely.

Per cell and promoter, methylated and total reads are pooled over the promoter's CpGs. Per promoter and age group,
the between-cell variance of its methylation beyond what read sampling gives is estimated by the method of moments,
``tau2 = (sum n_i (p_i - p)^2 - (C - 1) p (1 - p)) / (N - sum n_i^2 / N)``, and scaled by ``m (1 - m)`` of the young
level. Three readings per promoter: ``heterogeneity_young`` (how much cells of young mice differ), ``divergence``
(old minus young: what age adds between cells) and ``abs_change`` (the shift of the mean, scaled). Each is compared
between a gene set's promoters and promoters matched on island status and young-level decile, and again with the
decile of the gene's mean expression in these cells added to the strata. Intervals resample cells within each mouse;
the table also gives each reading with one mouse left out.

Limits: three mice per age, about 1.8 reads per CpG per cell, so only promoters with enough reads in enough cells
are used; scBS-seq reads 5-hydroxymethylcytosine as methylation; expression is the mean over the cells with counts.
``entropy-age figures musc_single_cell`` (reads the unpacked ``GSE121437_RAW`` folder).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.mouse import MouseAnnotationLayout
from entropy_age.data.rrbs import TRACKS, RrbsLayout, _columns, cpg_genes, in_intervals, islands, promoters
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.protection import LEVEL_CLIP, matched_scores
from entropy_age.stochasticity import baseline_strata

GSE = "GSE121437"
RNA_TABLE = "GSE121364_Table_raw_counts_RNA_QSC.txt.gz"
CELL = re.compile(r"_([OY]\d+)_TA_Hi_(\d+)")
MIN_READS = 4  # reads a cell needs in a promoter to count for it
MIN_CELLS = 20  # cells per age group a promoter needs
N_BOOT = {"paper": 200, "quick": 10}
READINGS = ("heterogeneity_young", "divergence", "abs_change")


def between_cell_variance(k: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per promoter (column): pooled level, between-cell variance beyond read sampling, and cells used. ``k`` and
    ``n`` are cells x promoters of methylated and total reads; cells with fewer than ``MIN_READS`` reads are left out
    of a promoter."""
    use = n >= MIN_READS
    n = np.where(use, n, 0.0)
    k = np.where(use, k, 0.0)
    total = n.sum(axis=0)
    cells = use.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k.sum(axis=0) / total
        rate = np.where(use, k / np.where(n > 0, n, 1.0), 0.0)
        weighted = (n * (rate - p) ** 2).sum(axis=0)
        tau2 = (weighted - (cells - 1) * p * (1 - p)) / (total - (n**2).sum(axis=0) / total)
    return p, tau2, cells


def readings(k: np.ndarray, n: np.ndarray, old: np.ndarray) -> dict[str, np.ndarray]:
    """The three per-promoter readings and the young level; NaN where a promoter lacks ``MIN_CELLS`` in a group."""
    p_y, t_y, c_y = between_cell_variance(k[~old], n[~old])
    p_o, t_o, c_o = between_cell_variance(k[old], n[old])
    ok = (c_y >= MIN_CELLS) & (c_o >= MIN_CELLS)
    m = np.clip(p_y, LEVEL_CLIP, 1 - LEVEL_CLIP)
    scale = m * (1 - m)
    out = {
        "heterogeneity_young": t_y / scale,
        "divergence": (t_o - t_y) / scale,
        "abs_change": np.abs(p_o - p_y) / np.sqrt(scale),
    }
    return {name: np.where(ok, v, np.nan) for name, v in out.items()} | {"young": np.where(ok, p_y, np.nan)}


def gaps(values: dict[str, np.ndarray], keep: np.ndarray, sets: dict[str, np.ndarray], strata: np.ndarray) -> dict:
    """``(reading, set)`` -> set minus matched promoters, in SD of the reading over the kept promoters."""
    out = {}
    for name in READINGS:
        v = values[name][keep]
        sd = float(np.nanstd(v))
        for group, in_set in sets.items():
            out[(name, group)] = float(matched_scores(v[None, :], in_set[keep], strata[keep])[0][0]) / sd
    return out


def read_cells(folder: Path, prom: pl.DataFrame) -> tuple[pl.DataFrame, list[str]]:
    """Promoter CpGs of every cell (``cell``, ``chr``, ``pos``, ``meth``, ``total``) and the cell names."""
    frames, names = [], []
    for path in sorted(folder.glob("*.cov.gz")):
        found = CELL.search(path.name)
        if found is None:
            continue
        name = f"{found.group(1)}_TA_Hi_{found.group(2)}"
        cov = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=["chr", "pos", "end", "pct", "meth", "unmeth"],
            schema_overrides={"chr": pl.String},
        )
        chrom = cov["chr"].str.strip_prefix("chr")
        inside = in_intervals(chrom.to_numpy(), cov["pos"].to_numpy() - 1, prom)
        frames.append(
            cov.filter(pl.Series(inside)).select(
                pl.lit(len(names), dtype=pl.Int32).alias("cell"),
                pl.col("chr").str.strip_prefix("chr"),
                "pos",
                pl.col("meth").cast(pl.Float64),
                (pl.col("meth") + pl.col("unmeth")).cast(pl.Float64).alias("total"),
            )
        )
        names.append(name)
    return pl.concat(frames), names


def expression_deciles(table: Path, symbols: dict[str, int]) -> dict[int, int]:
    """Entrez id -> decile (0-9) of the gene's mean log counts per million over the cells of the count table."""
    counts = pl.read_csv(table, separator=" ", infer_schema_length=0)
    cells = [c for c in counts.columns if CELL.search(c)]
    x = counts.select(pl.col(cells).cast(pl.Float64)).to_numpy()
    cpm = np.log2(x / x.sum(axis=0, keepdims=True) * 1e6 + 1.0).mean(axis=1)
    frame = (
        pl.DataFrame({"symbol": counts["geneName"], "level": cpm})
        .with_columns(entrez=pl.col("symbol").replace_strict(symbols, default=None, return_dtype=pl.Int64))
        .drop_nulls("entrez")
        .group_by("entrez")
        .agg(pl.col("level").max())
    )
    deciles = baseline_strata(frame["level"].to_numpy(), 10)
    return dict(zip(frame["entrez"].to_list(), [int(d) for d in deciles], strict=True))


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE / f"{GSE}_RAW"
    rr = RrbsLayout(ctx.paths)
    if not folder.exists() or any(not rr.track(t).exists() for t in TRACKS):
        print(f"musc_single_cell: {folder} or the mm10 tracks are missing; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("musc_single_cell")
    prom = promoters(rr.track(TRACKS[0]), rr.track(TRACKS[1]))
    isl = islands(rr.track(TRACKS[2]))
    cpgs, names = read_cells(folder, prom)
    sites = cpgs.select("chr", "pos").unique().with_columns(cpg=pl.format("{}:{}", "chr", "pos"))
    links = cpg_genes(sites, prom).rename({"feature": "cpg"})
    on_island = dict(
        zip(
            sites["cpg"].to_list(),
            in_intervals(sites["chr"].to_numpy(), sites["pos"].to_numpy() - 1, isl).tolist(),
            strict=True,
        )
    )
    per_gene = (
        cpgs.join(sites, on=["chr", "pos"])
        .join(links, on="cpg")
        .with_columns(island=pl.col("cpg").replace_strict(on_island, return_dtype=pl.Boolean))
        .group_by("cell", "entrez")
        .agg(pl.col("meth").sum(), pl.col("total").sum(), pl.col("island").mean())
    )
    genes = np.sort(per_gene["entrez"].unique().to_numpy())
    column = {int(g): j for j, g in enumerate(genes)}
    k = np.zeros((len(names), genes.size))
    n = np.zeros_like(k)
    rows = per_gene["cell"].to_numpy()
    cols = per_gene["entrez"].replace_strict(column, return_dtype=pl.Int64).to_numpy()
    k[rows, cols], n[rows, cols] = per_gene["meth"].to_numpy(), per_gene["total"].to_numpy()
    island_share = per_gene.group_by("entrez").agg(pl.col("island").mean()).sort("entrez")["island"].to_numpy()
    mouse = np.array([name.split("_")[0] for name in names])
    old = np.array([m.startswith("O") for m in mouse])
    gene_sets = pl.read_parquet(MouseAnnotationLayout(ctx.paths).gene_sets).filter(pl.col("gene_set").is_in(SETS))
    sets = {g: np.isin(genes, gene_sets.filter(pl.col("gene_set") == g)["entrez"].to_numpy()) for g in SETS}
    link = _columns(rr.track(TRACKS[1]), [2, 6], ["symbol", "entrez"], truncate_ragged_lines=True).with_columns(
        pl.col("entrez").cast(pl.Int64, strict=False)
    )
    symbols = dict(link.drop_nulls().unique(subset="symbol").iter_rows())
    decile_of = expression_deciles(ctx.paths.input / "geo" / GSE / RNA_TABLE, symbols)
    expressed = np.array([decile_of.get(int(g), -1) for g in genes])

    point = readings(k, n, old)
    keep = np.isfinite(point["divergence"]) & np.isfinite(point["heterogeneity_young"])
    level = np.full(genes.size, "", dtype=object)
    level[keep] = baseline_strata(point["young"][keep], 10).astype(str)
    base = np.char.add(np.char.add((island_share > 0.5).astype(str), "|"), level.astype(str))
    variants = {
        "island x level": base,
        "+ expression decile": np.char.add(np.char.add(base, "|"), expressed.astype(str)),
    }
    rng = ctx.rng(160_000)
    index = [np.flatnonzero(mouse == m) for m in np.unique(mouse)]
    table = []
    for matching, strata in variants.items():
        full = gaps(point, keep, sets, strata)
        draws: dict[tuple[str, str], list[float]] = {key: [] for key in full}
        for _ in range(N_BOOT[str(ctx.scale.name)]):
            pick = np.concatenate([rng.choice(i, i.size, replace=True) for i in index])
            boot = readings(k[pick], n[pick], old[pick])
            ok = keep & np.isfinite(boot["divergence"]) & np.isfinite(boot["heterogeneity_young"])
            for key, value in gaps(boot, ok, sets, strata).items():
                draws[key].append(value)
        left_out = {}
        for m in np.unique(mouse):
            rest = mouse != m
            loo = readings(k[rest], n[rest], old[rest])
            ok = keep & np.isfinite(loo["divergence"]) & np.isfinite(loo["heterogeneity_young"])
            left_out[m] = gaps(loo, ok, sets, strata)
        for (name, group), value in full.items():
            d = np.array(draws[(name, group)])
            table.append(
                {
                    "matching": matching,
                    "reading": name,
                    "group": group,
                    "gap_sd": value,
                    "ci_low": float(np.nanquantile(d, 0.025)),
                    "ci_high": float(np.nanquantile(d, 0.975)),
                    "promoters_in_set": int((sets[group] & keep).sum()),
                    **{f"without_{m}": left_out[m][(name, group)] for m in np.unique(mouse)},
                }
            )
    write_table(pl.DataFrame(table), out / "gaps")
    write_table(
        pl.DataFrame(
            {
                "cells_young": [int((~old).sum())],
                "cells_old": [int(old.sum())],
                "mice": [int(np.unique(mouse).size)],
                "promoters_with_reads": [int(genes.size)],
                "promoters_used": [int(keep.sum())],
                "used_with_expression": [int((keep & (expressed >= 0)).sum())],
                "median_divergence": [float(np.nanmedian(point["divergence"][keep]))],
                "median_heterogeneity_young": [float(np.nanmedian(point["heterogeneity_young"][keep]))],
            }
        ),
        out / "sizes",
    )
    print(f"musc_single_cell: {len(names)} cells, {int(keep.sum())} promoters used", file=sys.stderr, flush=True)
    return []
