"""Natural experiments in people: does a constitutional variant in a chromatin writer move the shield? (roadmap RM14)

Tet3 heterozygous mice lose the shield in brain against their littermates (run R15). People carry such
variants too, from birth, in every cell. Two public 450k blood series of the Weksberg lab cover four writers:

- GSE74432 (Choufani 2015): Sotos syndrome and other NSD1 variants (H3K36 methyltransferase), Weaver syndrome
  (EZH2, the Polycomb H3K27 methyltransferase), and controls. NSD1 loss is known to speed up epigenetic clocks
  (Martin-Herranz 2019). GEO gives no ages for this series, so its effects carry no age covariate.
- GSE97362 (Butcher 2017): Kabuki syndrome and other KMT2D variants (H3K4 methyltransferase), CHARGE syndrome
  and other CHD7 variants (a chromatin remodeller), and controls, with ages.

Every person is scored against the series' controls (``entropy_age.protection``; controls leave themselves
out, as in ``consortium_interventions.arm_scores``): ``P``, drift ``B``, the shield as a share of drift ``F``,
and the signed shift ``S`` (above 0: set CpGs gained methylation). The effect of a variant is the OLS
coefficient of carrying it, with sex and, where known, age as covariates. Blood cell mix differs between
syndromes, so the scores are computed twice: without composition removal and with the 5 leading components of
the controls' background CpGs removed (``COMPONENTS``).

Chromatin maintenance predicts a weaker shield where the writer that marks active promoters (KMT2D, H3K4) or
keeps their neighbours (NSD1) is halved; selection predicts nothing specific to the writer.

Run as a job after ``entropy-age fetch syndromes``: ``entropy-age figures syndromes``.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.per_person import _islands
from entropy_age.experiments.tissue_contrast import _features

SERIES = {
    "GSE74432": {"Sotos": "NSD1", "NSD1 variant": "NSD1", "Weaver": "EZH2"},
    "GSE97362": {"Kabuki": "KMT2D", "KMT2D variant": "KMT2D", "CHARGE": "CHD7", "CHD7 variant": "CHD7"},
}
CONTROL = "Control"
COMPONENTS = (0, 5)
MIN_ARM = 5


def matrix_path(paths, gse: str) -> Path:
    return paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"


def read_people(matrix: Path) -> pl.DataFrame:
    """One row per sample: accession, disease state, sex, tissue, age (null where GEO gives none)."""
    chars = _characteristics(matrix)
    accessions = []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                accessions = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
                break
    age = next((c for c in chars.columns if c.startswith("age")), None)
    return pl.DataFrame({"sample": accessions}).join(
        chars.select(
            "sample",
            pl.col("disease state").alias("disease"),
            pl.col("gender").alias("sex"),
            pl.col("tissue"),
            (pl.col(age).cast(pl.Float64, strict=False) if age else pl.lit(None, dtype=pl.Float64)).alias("age"),
        ),
        on="sample",
    )


def with_known_age(people: pl.DataFrame) -> pl.DataFrame:
    """Where a series gives ages, only the people with one, so age can be a covariate; else everyone."""
    if people["age"].is_null().all():
        return people
    dropped = people.filter(pl.col("age").is_null())
    if dropped.height:
        print(f"syndromes: {dropped.height} without an age left out: {dropped['sample'].to_list()}", file=sys.stderr)
    return people.filter(pl.col("age").is_not_null())


def read_betas(matrix: Path, features: list[str], samples: list[str]) -> np.ndarray:
    """The series matrix's beta table, samples by ``features`` (NaN where a CpG is missing)."""
    index = {f: i for i, f in enumerate(features)}
    x = np.full((len(samples), len(features)), np.nan)
    columns = None
    with gzip.open(matrix, "rt") as handle:
        inside = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not inside:
                continue
            cells = line.rstrip("\n").split("\t")
            if columns is None:
                header = [c.strip('"') for c in cells[1:]]
                columns = np.array([header.index(s) for s in samples])
                continue
            j = index.get(cells[0].strip('"'))
            if j is not None:
                row = np.array([float(c) if c not in ("", "null", "NA", '""') else np.nan for c in cells[1:]])
                x[:, j] = row[columns]
    if columns is None:
        raise ValueError(f"{matrix}: no beta table")
    return x


def effects(scores: dict[str, np.ndarray], carrier: np.ndarray, covariates: np.ndarray) -> list[dict]:
    """Per set and measure: the carrier coefficient, its SE and p, and the controls' mean."""
    design = np.column_stack([np.ones(carrier.size), carrier.astype(float), covariates])
    rows = []
    for key, y in scores.items():
        measure, _, group = key.partition("_")
        if measure not in ("B", "F", "S"):
            continue
        coef, se, p = ols(y, design, 1)
        rows.append(
            {
                "group": group,
                "measure": measure,
                "effect": coef,
                "se": se,
                "p": p,
                "control_mean": float(y[~carrier].mean()),
            }
        )
    return rows


def run(ctx: Context) -> list[Path]:
    missing = [g for g in SERIES if not matrix_path(ctx.paths, g).exists()]
    if missing:
        print(
            f"syndromes: run `entropy-age fetch syndromes` first ({', '.join(missing)}). Nothing run.", file=sys.stderr
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("syndromes")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    island = _islands(features, probes)
    rows, scored = [], []
    for gse, arms in SERIES.items():
        people = with_known_age(read_people(matrix_path(ctx.paths, gse)).filter(pl.col("tissue") == "whole blood"))
        x_all = read_betas(matrix_path(ctx.paths, gse), features, people["sample"].to_list())
        ok = ~np.isnan(x_all).any(axis=0)
        x_all = x_all[:, ok]
        kept = [f for f, k in zip(features, ok, strict=True) if k]
        sets = _sets(kept, membership)
        print(f"syndromes: {gse}: {people.height} blood samples, {len(kept)} CpGs", file=sys.stderr, flush=True)
        is_control = (people["disease"] == CONTROL).to_numpy()
        for arm, gene in arms.items():
            is_arm = (people["disease"] == arm).to_numpy()
            if is_arm.sum() < MIN_ARM:
                print(f"syndromes: {gse}, {arm}: {int(is_arm.sum())} people; left out", file=sys.stderr)
                continue
            keep = is_control | is_arm
            sel = people.filter(pl.Series(keep))
            carrier = is_arm[keep]
            ages = sel["age"].to_numpy()
            cov = [dummies(sel["sex"].to_numpy())]
            if not np.isnan(ages.astype(float)).any():
                cov.append(ages[:, None].astype(float))
            covariates = np.column_stack(cov)
            for k in COMPONENTS:
                scores = arm_scores(x_all[keep], carrier, island[ok], sets, components=k)
                scored.append(
                    pl.DataFrame(
                        {
                            "sample": sel["sample"],
                            "series": gse,
                            "arm": arm,
                            "gene": gene,
                            "components": k,
                            "carrier": carrier,
                            "age": sel["age"],
                            **scores,
                        }
                    )
                )
                for r in effects(scores, carrier, covariates):
                    rows.append(
                        {
                            "series": gse,
                            "arm": arm,
                            "gene": gene,
                            "components": k,
                            "age_covariate": covariates.shape[1] > dummies(sel["sex"].to_numpy()).shape[1],
                            "carriers": int(carrier.sum()),
                            "controls": int((~carrier).sum()),
                            **r,
                        }
                    )
            print(f"syndromes: {gse}, {arm}: done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    write_table(pl.concat(scored, how="diagonal"), out / "person_scores")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    groups = [g for g in ("core_essential", "housekeeping", "non_essential") if g in table["group"].unique()]
    fig, axes = plots.panels(3, ncols=3, size=4.2)
    for i, measure in enumerate(("F", "B", "S")):
        ax = axes[i]
        for j, k in enumerate(COMPONENTS):
            sub = table.filter((pl.col("measure") == measure) & (pl.col("components") == k))
            labels, ys = [], []
            for r in sub.sort("series", "gene", "arm", "group").iter_rows(named=True):
                if r["group"] not in groups:
                    continue
                y = len(labels)
                labels.append(f"{r['arm']} ({r['gene']}) · {r['group'].replace('_', ' ')}")
                ys.append(y)
                ax.errorbar(
                    r["effect"],
                    y + 0.2 * j,
                    xerr=1.96 * r["se"],
                    fmt="o" if k == 0 else "s",
                    ms=3.5,
                    color=plots.PALETTE[groups.index(r["group"])],
                    alpha=1.0 if k else 0.5,
                )
            if j == 0:
                ax.set_yticks(ys, labels if i == 0 else [""] * len(labels), fontsize=6)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        title = {
            "F": "Shield share F (> 0: weaker shield)",
            "B": "Drift B (> 0: more drift)",
            "S": "Signed shift S (> 0: set CpGs gained methylation)",
        }[measure]
        ax.set(xlabel="Carriers minus controls (95 % CI)", title=title)
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Constitutional chromatin-writer variants in blood: circles raw, squares with 5 composition components",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("syndromes"))
