"""Selection among blood cells in vivo: does clonal haematopoiesis move the shield? (roadmap RM16)

Clonal haematopoiesis (CH) is selection among blood stem cells in a living person: a clone with a driver mutation
(DNMT3A, TET2) outgrows the rest. If part of blood's shield is selection, CH carriers should differ. Basrai et al.
2026 (GSE315366) measured EPIC methylation in blood of 64 patients in cardiogenic shock, with CH called per person
(DNMT3A and TET2; no allele fractions). The patients are acutely ill, so this is a comparison within one sick
cohort, not against healthy people.

Arms (``ARMS``) are scored against the patients without CH exactly as ``syndromes`` scores carriers: controls leave
themselves out, people outside the controls' ages are left out, age and sex are covariates, and the scores are
computed raw and with 5 composition components.

``entropy-age fetch chip`` then ``entropy-age figures chip``.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.per_person import _islands
from entropy_age.experiments.syndromes import COMPONENTS, CONTROL, effects, with_known_age
from entropy_age.experiments.tissue_contrast import SET_LABELS, _features

GSE = "GSE315366"
BETAS = f"{GSE}_cardiogenic_shock_beta_values.txt.gz"
ARMS = {"any CH": "ch", "DNMT3A": "dnmt3a", "TET2": "tet2"}
MIN_ARM = 5


def layout(paths) -> tuple[Path, Path]:
    base = paths.input / "geo" / GSE
    return base / f"{GSE}_series_matrix.txt.gz", base / BETAS


def read_patients(matrix: Path) -> pl.DataFrame:
    """``sample`` (the title, as the beta columns name it), ``age``, ``sex``, flags ``ch``, ``dnmt3a``, ``tet2`` and
    ``disease`` (``Control`` without CH, else ``CH``)."""
    chars = _characteristics(matrix)
    titles = accessions = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_geo_accession"):
                accessions = cells
    flag = {"clonal hematopoiesis": "ch", "dnmt3a mutation": "dnmt3a", "tet2 mutation": "tet2"}
    out = pl.DataFrame({"accession": accessions, "sample": titles}).join(
        chars.select(
            pl.col("sample").alias("accession"),
            pl.col("age").cast(pl.Float64, strict=False),
            pl.col("Sex").alias("sex"),
            *(pl.col(k).str.to_uppercase().eq("TRUE").alias(v) for k, v in flag.items()),
        ),
        on="accession",
    )
    return out.with_columns(disease=pl.when(pl.col("ch")).then(pl.lit("CH")).otherwise(pl.lit(CONTROL)))


def read_beta_table(path: Path, features: list[str], samples: list[str]) -> np.ndarray:
    """A tab-separated ``ID_REF`` x sample beta table, samples by ``features`` (NaN where absent)."""
    index = {f: i for i, f in enumerate(features)}
    x = np.full((len(samples), len(features)), np.nan)
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = np.array([header.index(s) - 1 for s in samples])  # values start after ID_REF
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            j = index.get(cells[0])
            if j is not None:
                row = np.array([float(c) if c not in ("", "NA", "NaN") else np.nan for c in cells[1:]])
                x[:, j] = row[columns]
    return x


def run(ctx: Context) -> list[Path]:
    matrix, betas = layout(ctx.paths)
    if not (matrix.exists() and betas.exists()):
        print("chip: run `entropy-age fetch chip` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("chip_methylation")
    features, probes, membership = _features(ctx.paths, ctx.rng(150_000))
    island = _islands(features, probes)
    people = with_known_age(read_patients(matrix))
    x_all = read_beta_table(betas, features, people["sample"].to_list())
    ok = ~np.isnan(x_all).any(axis=0)
    x_all = x_all[:, ok]
    kept = [f for f, k in zip(features, ok, strict=True) if k]
    sets = _sets(kept, membership)
    print(f"chip: {people.height} patients, {len(kept)} CpGs", file=sys.stderr, flush=True)
    is_control = (people["disease"] == CONTROL).to_numpy()
    rows = []
    for arm, column in ARMS.items():
        is_arm = people[column].to_numpy()
        if is_arm.sum() < MIN_ARM:
            print(f"chip: {arm}: {int(is_arm.sum())} carriers; left out", file=sys.stderr)
            continue
        keep = is_control | is_arm
        sel = people.filter(pl.Series(keep))
        carrier = is_arm[keep]
        covariates = np.column_stack([dummies(sel["sex"].to_numpy()), sel["age"].to_numpy()])
        for k in COMPONENTS:
            scores = arm_scores(x_all[keep], carrier, island[ok], sets, components=k)
            for r in effects(scores, carrier, covariates):
                rows.append(
                    {
                        "arm": arm,
                        "components": k,
                        "carriers": int(carrier.sum()),
                        "controls": int((~carrier).sum()),
                        **r,
                    }
                )
        print(f"chip: {arm}: done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=4.0)
    groups = [g for g in ("core_essential", "housekeeping", "non_essential") if g in table["group"].unique()]
    for i, measure in enumerate(("F", "B", "S")):
        ax = axes[i]
        labels, y = [], 0
        for arm in ARMS:
            for group in groups:
                for j, k in enumerate(COMPONENTS):
                    r = table.filter(
                        (pl.col("arm") == arm)
                        & (pl.col("group") == group)
                        & (pl.col("measure") == measure)
                        & (pl.col("components") == k)
                    )
                    if r.is_empty():
                        continue
                    ax.errorbar(
                        r["effect"][0],
                        y + 0.25 * j,
                        xerr=1.96 * r["se"][0],
                        fmt="os"[j],
                        ms=4,
                        color=plots.PALETTE[groups.index(group)],
                        alpha=1.0 if k else 0.5,
                    )
                labels.append((y + 0.12, f"{arm} · {SET_LABELS[group]}"))
                y += 1
            y += 0.5
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels] if i == 0 else [""] * len(labels), fontsize=7)
        ax.invert_yaxis()
        ax.set(
            xlabel="Carriers minus non-carriers (95 % CI)",
            title={"F": "Shield share F (> 0: weaker)", "B": "Drift B", "S": "Signed shift S"}[measure],
        )
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Clonal haematopoiesis in cardiogenic-shock blood (GSE315366): circles raw, squares 5 components",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("chip_methylation"))
