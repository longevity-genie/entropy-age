"""After a stem-cell bottleneck, does the shield grow? Transplant donors and their recipients (roadmap RM12)

A haematopoietic stem-cell transplant passes a donor's blood through a bottleneck: a few stem cells rebuild a
recipient's blood, and clones expand under strong competition. Selection among cells predicts that a
recipient's blood, made of the donor's cells, has drifted from the donor's own blood less at essential promoters
than elsewhere, by more than ordinary aging divergence would give. Chromatin predicts the same share either way.

Holland et al. 2024 (Aging Cell; Mendeley Data j5krtjfj6y; ``data/input/mendeley/j5krtjfj6y``) joined the
Stölzel 2017 and Søraas 2019 cohorts: 450k/EPIC betas (ENmix) in six batch files, sample ``<pair>_<D|R>_<batch>_<n>``,
with the donor's and recipient's ages at transplant and the years after it. Donors were sampled at transplant,
recipients 0-23 years after, several times each. Every comparison stays inside one batch:

* **donor → recipient**: each recipient draw against its own donor's sample (the bottleneck, and the years);
* **recipient → recipient**: consecutive draws of one recipient (years, no new bottleneck);
* **donor → donor**: a donor's repeated samples, the measurement floor.

For each comparison, every CpG's squared difference is taken; the set's matched gap ``P`` against look-alikes of
the same island status and level decile, the look-alikes' mean ``B``, and ``F = P / B``, the share of divergence
the shield removes. ``F`` is regressed on the kind of comparison, log(1 + years apart) and batch; the interval
comes from resampling transplant pairs.

``entropy-age figures hsct_bottleneck`` as a Slurm job (reads about 3.8 GB of CSV).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata

FOLDER = ("mendeley", "j5krtjfj6y")
BATCHES = tuple(f"batch{i}.csv" for i in range(1, 7))
GROUPS = ("core_essential", "housekeeping", "non_essential")
KINDS = ("donor → recipient", "recipient → recipient", "donor → donor")
N_BOOT = {"paper": 1000, "quick": 100}


def parse_ids(ids: list[str]) -> pl.DataFrame:
    """``pair``, ``role`` (D or R), ``batch`` and ``n`` from ``<pair>_<role>_<batch>_<n>``."""
    parts = [i.split("_") for i in ids]
    return pl.DataFrame(
        {
            "sample": ids,
            "pair": [p[0] for p in parts],
            "role": [p[1] for p in parts],
            "batch": [p[2] for p in parts],
            "n": [int(p[3]) for p in parts],
        }
    )


def comparisons(samples: pl.DataFrame) -> pl.DataFrame:
    """Rows ``first``, ``second`` (row indices into ``samples``), ``kind``, ``years`` apart, ``pair``, ``batch``.

    ``samples`` needs ``pair``, ``role``, ``batch`` and ``years`` (after transplant; donors 0).
    """
    s = samples.with_row_index("row")
    rows = []
    for (pair, batch), g in s.group_by("pair", "batch", maintain_order=True):
        donors = g.filter(pl.col("role") == "D").sort("row")
        recips = g.filter(pl.col("role") == "R").sort("years", "row")
        for d in donors.iter_rows(named=True):
            for r in recips.iter_rows(named=True):
                rows.append((d["row"], r["row"], KINDS[0], r["years"], pair, batch))
        r = recips.rows(named=True)
        for a, b in zip(r, r[1:], strict=False):
            rows.append((a["row"], b["row"], KINDS[1], b["years"] - a["years"], pair, batch))
        d = donors.rows(named=True)
        for a, b in zip(d, d[1:], strict=False):
            rows.append((a["row"], b["row"], KINDS[2], 0.0, pair, batch))
    return pl.DataFrame(rows, schema=["first", "second", "kind", "years", "pair", "batch"], orient="row")


def divergence_scores(x: np.ndarray, first: np.ndarray, second: np.ndarray, sets: dict, strata: np.ndarray) -> dict:
    """Per comparison and set: ``P`` (set minus matched squared difference), ``B`` (matched mean) and ``F = P / B``."""
    d = (x[second] - x[first]) ** 2
    out = {}
    for group, (in_set, _) in sets.items():
        p, b = matched_scores(d, in_set, strata)
        out[f"P_{group}"], out[f"B_{group}"] = p, b
        out[f"F_{group}"] = np.where(b > 0, p / np.where(b > 0, b, 1.0), np.nan)
    return out


def kind_effect(table: pl.DataFrame, y: str) -> tuple[float, np.ndarray]:
    """OLS coefficient of donor → recipient against recipient → recipient on ``y``, with log(1 + years) and batch."""
    t = table.filter(pl.col("kind").is_in(KINDS[:2]) & pl.col(y).is_finite())
    yy = t[y].to_numpy()
    design = np.column_stack(
        [
            np.ones(t.height),
            (t["kind"] == KINDS[0]).to_numpy().astype(float),
            np.log1p(np.clip(t["years"].to_numpy(), 0, None)),
            dummies(t["batch"].to_numpy()),
        ]
    )
    coef, *_ = np.linalg.lstsq(design, yy, rcond=None)
    return float(coef[1]), coef


def read_batches(folder: Path, wanted: set[str]) -> tuple[list[str], list[str], np.ndarray]:
    """CpGs measured in every batch (no missing value), sample ids, and the (samples, CpGs) betas."""
    frames = []
    for name in BATCHES:
        path = folder / name
        header = pl.read_csv(path, n_rows=0).columns
        f = (
            pl.scan_csv(path, null_values=["NA", ""])
            .filter(pl.col(header[0]).is_in(list(wanted)))
            .collect()
            .rename({header[0]: "cpg"})
        )
        frames.append(f.drop_nulls())
        print(f"hsct_bottleneck: {name}: {f.height} CpGs, {len(header) - 1} samples", file=sys.stderr, flush=True)
    joined = frames[0]
    for f in frames[1:]:
        joined = joined.join(f, on="cpg", how="inner")
    joined = joined.sort("cpg")
    samples = [c for c in joined.columns if c != "cpg"]
    return joined["cpg"].to_list(), samples, joined.select(samples).to_numpy().T.astype(float)


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input.joinpath(*FOLDER)
    if not all((folder / b).exists() for b in BATCHES) or not (folder / "comRand.csv").exists():
        print(f"hsct_bottleneck: {folder} incomplete; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("hsct_bottleneck")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    cpgs, ids, x = read_batches(folder, set(features))
    meta = pl.read_csv(folder / "comRand.csv").rename({"randID": "sample", "YrsAfterReceivingTX": "years"})
    samples = parse_ids(ids).join(meta, on="sample", how="left").with_columns(pl.col("years").fill_null(0.0))
    sets = {g: v for g, v in _sets(cpgs, membership).items() if g in GROUPS}
    strata = np.char.add(_islands(cpgs, probes).astype(str), baseline_strata(x.mean(axis=0), 10).astype(str))
    comp = comparisons(samples)
    print(
        f"hsct_bottleneck: {len(cpgs)} CpGs, {len(ids)} samples; "
        + ", ".join(f"{k}: {(comp['kind'] == k).sum()}" for k in KINDS),
        file=sys.stderr,
        flush=True,
    )
    scores = divergence_scores(x, comp["first"].to_numpy(), comp["second"].to_numpy(), sets, strata)
    table = comp.with_columns(**{k: pl.Series(v) for k, v in scores.items()})
    write_table(table, out / "comparisons")

    rng = ctx.rng(250_000)
    pairs = table["pair"].unique().to_list()
    rows = []
    for group in sets:
        for measure in ("F", "B"):
            y = f"{measure}_{group}"
            est, _ = kind_effect(table, y)
            draws = []
            for _ in range(N_BOOT[str(ctx.scale.name)]):
                pick = rng.choice(pairs, size=len(pairs), replace=True)
                boot = pl.concat([table.filter(pl.col("pair") == p) for p in pick])
                if boot.filter(pl.col("kind") == KINDS[0]).height < 3:
                    continue
                draws.append(kind_effect(boot, y)[0])
            draws = np.array(draws)
            means = {k: float(table.filter(pl.col("kind") == k)[y].mean() or np.nan) for k in KINDS}
            rows.append(
                {
                    "group": group,
                    "measure": measure,
                    "bottleneck_effect": est,
                    "low": float(np.quantile(draws, 0.025)),
                    "high": float(np.quantile(draws, 0.975)),
                    "p_boot": float(2 * min((draws <= 0).mean(), (draws >= 0).mean())),
                    **{f"mean {k}": v for k, v in means.items()},
                    **{f"n {k}": int((table["kind"] == k).sum()) for k in KINDS},
                }
            )
    effects = pl.DataFrame(rows)
    write_table(effects, out / "effects")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for ax, group, letter in zip(axes, ("housekeeping", "core_essential"), "ab", strict=True):
        for i, kind in enumerate(KINDS):
            t = table.filter((pl.col("kind") == kind) & pl.col(f"F_{group}").is_finite())
            if t.height == 0:
                continue
            ax.scatter(
                np.log1p(np.clip(t["years"].to_numpy(), 0, None)),
                t[f"F_{group}"],
                s=10,
                alpha=0.6,
                color=plots.PALETTE[i],
                label=f"{kind} ({t.height})",
            )
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set(
            xlabel="log(1 + years apart)",
            ylabel="shield share of divergence F",
            title=group.replace("_", " "),
        )
        if letter == "a":
            ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(1.1, -0.15), ncol=3, frameon=False)
        plots.label_panel(ax, letter)
    fig.suptitle(
        "Divergence after a stem-cell transplant (Holland 2024)",
        x=0.02,
        y=1.03,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.subplots_adjust(wspace=0.35)
    return plots.save(fig, ctx.figure_path("hsct_bottleneck"))
