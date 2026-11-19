"""Tet3 against Tet1 on the CpGs their two arrays share (roadmap RM43 item 8).

Tet3 heterozygous cortex and striatum lost part of the shield (R15: consortium data on the 40k mammalian array,
deposited betas), while Tet1 knockout and catalytic-dead cortex gained it (R99: Foong et al. 2024 on the 285k mouse
array, rebuilt from IDATs). The two results were scored on different CpGs. Here each is scored twice: on every
promoter CpG its own array measures in all its animals, and on the promoter CpGs both arrays measure (matched by the
probe's cg identifier, the design suffix of the mouse array removed). A different sign on the shared CpGs is then a
difference between the models and their series, not between the CpGs each array happens to carry.

The arrays turn out to share almost nothing: 150 CpG positions in the mouse genome, 42 of them at promoters (R120),
so the shared-CpG score rests on a handful of CpGs and is reported only to show that. The usable version is by gene:
each array's own CpGs, with every gene set cut down to the genes that have promoter CpGs on both arrays
(``shared genes``), so both models are scored on the same genes.

What stays different whatever the CpGs: the array chemistry and processing (deposited against rebuilt), the
laboratory, heterozygous against homozygous loss, and age. ``entropy-age figures tet_shared`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments import tet1_cortex
from entropy_age.experiments.bxd_liver import ols, probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.consortium_interventions import CONTROLS, arm_scores, chip_design, same_chip
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.species_shield import layout as consortium_layout
from entropy_age.experiments.species_shield import read_samples
from entropy_age.experiments.tissue_contrast import SETS

TET3_ARM = "Tet3KO.Het"
TET3_TISSUES = ("CerebralCortex", "Striatum")
MEASURES = ("B", "F", "S")
SUFFIX = r"_[A-Z]{2}\d+$"


def base_ids(wide: pl.DataFrame) -> pl.DataFrame:
    """``wide`` (``cpg`` plus sample columns) with the design suffix removed from ``cpg``; the first probe of a CpG
    measured more than once is kept."""
    return wide.with_columns(pl.col("cpg").str.replace(SUFFIX, "")).unique(
        subset="cpg", keep="first", maintain_order=True
    )


def complete(wide: pl.DataFrame, samples: list[str]) -> tuple[np.ndarray, list[str]]:
    """Samples x CpGs matrix of the CpGs measured in every sample, and their names."""
    x = wide.select(samples).to_numpy().T.astype(float)
    ok = ~np.isnan(x).any(axis=0)
    return x[:, ok], [f for f, k in zip(wide["cpg"].to_list(), ok, strict=True) if k]


def effects(
    x: np.ndarray,
    features: list[str],
    keep: list[str],
    treated: np.ndarray,
    design: np.ndarray,
    probes,
    members,
    genes: set[str] | None = None,
) -> list[dict]:
    """Arm effect on B, F and S of every gene set, scored on the CpGs of ``keep`` only; with ``genes``, the sets
    keep only those genes' CpGs and the other set CpGs leave the analysis (they are not matched controls)."""
    wanted = set(keep)
    if genes is not None:
        wanted -= set(members.filter(~pl.col("gene").is_in(list(genes)))["feature"].to_list())
    cols = np.array([f in wanted for f in features])
    feats = [f for f, c in zip(features, cols, strict=True) if c]
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    sets = _sets(feats, members.filter(pl.col("feature").is_in(feats)))
    scores = arm_scores(x[:, cols], treated, island, sets)
    rows = []
    for group, (in_set, _) in sets.items():
        for measure in MEASURES:
            y = scores[f"{measure}_{group}"]
            coef, se, p = ols(y, design, 1)
            rows.append(
                {
                    "group": group,
                    "measure": measure,
                    "effect": coef,
                    "se": se,
                    "p": p,
                    "control_mean": float(y[~treated].mean()),
                    "treated_n": int(treated.sum()),
                    "control_n": int((~treated).sum()),
                    "cpgs": len(feats),
                    "set_cpgs": int(in_set.sum()),
                }
            )
    return rows


def run(ctx: Context) -> list[Path]:
    lay = MouseAnnotationLayout(ctx.paths)
    matrix, betas_path = consortium_layout(ctx.paths)
    tet1_path = tet1_cortex.layout(ctx.paths)
    if not (lay.probes.exists() and betas_path.exists() and tet1_path.exists()):
        print("tet_shared: needs `fetch mammals`, `fetch tet1` and the mouse annotation. Nothing run.", file=sys.stderr)
        return []
    out = ctx.out_dir("tet_shared")
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    members = mouse_membership(probes, pl.read_parquet(lay.gene_sets)).filter(pl.col("group").is_in(SETS))
    promoter = probes["cpg"].implode()

    # Tet1: the 285k mouse array
    if reprocessed.in_use():
        names = tet1_cortex.bisulfite_animals(tet1_path.parent / f"{tet1_cortex.GSE}_series_matrix.txt.gz")
        wide1 = reprocessed.by_accession(ctx.paths, tet1_cortex.GSE, names)
    else:
        wide1 = tet1_cortex.wide_bisulfite(
            pl.read_csv(tet1_path, null_values=["NA"], schema_overrides={"Betas": pl.Float64})
        )
    wide1 = base_ids(wide1).filter(pl.col("cpg").is_in(promoter))
    animals = [c for c in wide1.columns if c != "cpg"]
    x1, feats1 = complete(wide1, animals)
    geno = np.array([tet1_cortex.genotype(a) for a in animals])

    # Tet3: the 40k mammalian array, same-chip controls
    samples = read_samples(matrix).filter(
        (pl.col("species") == "Mouse")
        & pl.col("barcode").is_not_null()
        & pl.col("tissue").is_in(TET3_TISSUES)
        & pl.col("intervention").is_in([*CONTROLS, TET3_ARM])
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    samples = samples.filter(pl.col("barcode").is_in(header))
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([header[0], *samples["barcode"].unique().to_list()])
        .filter(pl.col(header[0]).is_in(promoter))
        .collect()
        .rename({header[0]: "cpg"})
    )
    tet3 = {}
    for tissue in TET3_TISSUES:
        sel = samples.filter(pl.col("tissue") == tissue)
        is_arm = (sel["intervention"] == TET3_ARM).to_numpy()
        chips = sel["barcode"].str.extract(r"^(\d+)_").to_numpy()
        keep = same_chip(chips, is_arm, sel["intervention"].is_in(list(CONTROLS)).to_numpy())
        sel, treated, chips = sel.filter(pl.Series(keep)), is_arm[keep], chips[keep]
        x, feats = complete(frame, sel["barcode"].to_list())
        ages = sel["age"].to_numpy()
        design = [np.ones(treated.size), treated.astype(float)]
        if np.ptp(ages) > 0:
            design.append(ages)
        tet3[tissue] = (x, feats, treated, np.column_stack([*design, chip_design(chips)]))

    shared = set(feats1)
    for _, feats, _, _ in tet3.values():
        shared &= set(feats)
    shared = sorted(shared)
    genes = set(members.filter(pl.col("feature").is_in(feats1))["gene"].to_list())
    for _, feats, _, _ in tet3.values():
        genes &= set(members.filter(pl.col("feature").is_in(feats))["gene"].to_list())
    print(
        f"tet_shared: Tet1 {len(feats1)} CpGs, Tet3 {[len(v[1]) for v in tet3.values()]}, "
        f"shared {len(shared)}, genes {len(genes)}",
        file=sys.stderr,
        flush=True,
    )
    rows = []
    for arm in tet1_cortex.ARMS:
        keep = (geno == tet1_cortex.CONTROL) | (geno == arm)
        treated = geno[keep] == arm
        design = np.column_stack([np.ones(treated.size), treated.astype(float)])
        for label, cpgs in (("own array", feats1), ("shared", shared)):
            for r in effects(x1[keep], feats1, cpgs, treated, design, probes, members):
                rows.append({"model": f"Tet1 {arm}", "tissue": "cortex", "cpg_set": label, **r})
        for r in effects(x1[keep], feats1, feats1, treated, design, probes, members, genes):
            rows.append({"model": f"Tet1 {arm}", "tissue": "cortex", "cpg_set": "shared genes", **r})
    for tissue, (x, feats, treated, design) in tet3.items():
        for label, cpgs in (("own array", feats), ("shared", shared)):
            for r in effects(x, feats, cpgs, treated, design, probes, members):
                rows.append({"model": "Tet3 het", "tissue": tissue, "cpg_set": label, **r})
        for r in effects(x, feats, feats, treated, design, probes, members, genes):
            rows.append({"model": "Tet3 het", "tissue": tissue, "cpg_set": "shared genes", **r})
    write_table(pl.DataFrame(rows), out / "effects")
    return []
