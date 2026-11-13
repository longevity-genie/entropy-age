"""Does the shield hold within a person over time? SATSA repeated draws (roadmap RM28, repeated measures)

Every earlier test compared different people of different ages. SATSA (Swedish Adoption/Twin Study of Aging,
BioStudies S-BSST1206) drew blood from the same 596 people up to 8 times over 20 years (ages 48-99; 392 with two or
more draws), on the 450k and later the EPIC array, preprocessed together (noob, functional normalisation, cell
counts and batch regressed out, CpGs differing by over 0.15 between chips removed). Two caveats follow: cell mix is
already removed (good for the shield, which is not about composition), and ComBat may shrink spread.

**CpG level.** For each pair of consecutive draws of one person, the change ``delta = x_later - x_earlier`` per
CpG. Its mean part is the shared, directed change: ``delta = b * dt + c * dchip`` over all pairs (``dchip`` -1, 0
or +1 when the array changed). Its remainder ``r`` is what this person did differently. A random walk spreads as
``E[r^2] = 2 s^2 + q * dt``: the intercept ``2 s^2`` is the measurement floor (with its own term when the chip
changed), the slope ``q`` the **divergence per year**, the stochastic drift itself, free of measurement noise and of
differences between birth cohorts. Per CpG: ``|b|`` (directed), ``q`` (divergence), the floor. Essential and
housekeeping CpGs are compared with look-alikes matched on island status and level decile (``tissue_contrast``).

A positive control runs the usual cross-sectional comparison (``tissue_contrast``'s lifespan metrics) on each
person's first draw: if the shield is absent there too, the preprocessing removed it, not the time axis.

**Person level.** B, F and S per draw against one reference fitted on all draws; the within-person slope of each
on age (person means removed, chip held), with a person-level bootstrap. Does a person's own shield share change
as they age?

``entropy-age figures within_person`` (as a Slurm job; five 1.35 GB tables).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.sjlife_deaths import reference_scores
from entropy_age.experiments.tissue_contrast import LIFESPAN_METRICS, _features, _matched_sd, _with_strata
from entropy_age.stochasticity import drift_metrics

STUDY = "S-BSST1206"
TABLES = tuple(f"SATSA_DNA_meth_{i}.txt" for i in range(1, 6))
METRICS = ("abs_directed", "divergence", "floor")
GROUPS = ("core_essential", "housekeeping")
BOOTSTRAPS = 1000


def read_people(pheno: Path) -> pl.DataFrame:
    """One row per draw: ``sample``, ``person``, ``age``, ``female``, ``epic``, ``smoke``, ``wave``."""
    p = pl.read_csv(pheno, separator="\t", infer_schema_length=0)
    return p.select(
        pl.col("SID").alias("sample"),
        pl.col("IID").alias("person"),
        pl.col("AGE").cast(pl.Float64).alias("age"),
        (pl.col("SEX") == "Woman").alias("female"),
        (pl.col("CHIP") == "EPIC").alias("epic"),
        pl.col("SMOKE").alias("smoke"),
        pl.col("STUDY").alias("wave"),
    ).sort("person", "age")


def read_rows(tables: list[Path], wanted: set[str], samples: list[str]) -> tuple[list[str], np.ndarray]:
    """The ``wanted`` CpGs for ``samples`` (samples x CpGs) from R-style tables (header one cell short)."""
    names, rows = [], []
    for table in tables:
        with open(table) as handle:
            header = [h.strip('"') for h in handle.readline().rstrip("\n").split("\t")]
            if header and header[0] in ("", "cpg", "CpG", "ID"):
                header = header[1:]
            index = np.array([header.index(s) for s in samples]) + 1
            for line in handle:
                cpg, _, _ = line.partition("\t")
                cpg = cpg.strip('"')
                if cpg not in wanted:
                    continue
                cells = np.array(line.rstrip("\n").split("\t"), dtype=object)[index]
                rows.append(np.array([float(v) if v not in ("", "NA") else np.nan for v in cells]))
                names.append(cpg)
    return names, np.array(rows).T


def consecutive_pairs(person: np.ndarray, age: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Indices (earlier, later) of consecutive draws of each person; rows sorted by person then age."""
    order = np.lexsort((age, person))
    p, o = person[order], order
    same = p[1:] == p[:-1]
    return o[:-1][same], o[1:][same]


def pair_metrics(x: np.ndarray, age: np.ndarray, epic: np.ndarray, first: np.ndarray, second: np.ndarray) -> dict:
    """Per CpG: directed slope ``b`` per year, divergence ``q`` per year and the measurement floor (see module)."""
    dx = x[second] - x[first]
    dt = age[second] - age[first]
    dchip = epic[second].astype(float) - epic[first].astype(float)
    mean_design = np.column_stack([dt, dchip]) if np.any(dchip) else dt[:, None]
    coef, *_ = np.linalg.lstsq(mean_design, dx, rcond=None)
    r2 = (dx - mean_design @ coef) ** 2
    changed = (dchip != 0).astype(float)
    var_design = (
        np.column_stack([np.ones(dt.size), dt, changed]) if changed.any() else np.column_stack([np.ones(dt.size), dt])
    )
    vcoef, *_ = np.linalg.lstsq(var_design, r2, rcond=None)
    return {"directed": coef[0], "divergence": vcoef[1], "floor": vcoef[0]}


def within_slope(v: np.ndarray, age: np.ndarray, person: np.ndarray, epic: np.ndarray) -> float:
    """Slope of ``v`` on age with each person's mean removed (age, chip and ``v`` demeaned within person)."""
    keys, inv = np.unique(person, return_inverse=True)

    def demean(a):
        sums = np.bincount(inv, weights=a, minlength=keys.size)
        counts = np.bincount(inv, minlength=keys.size)
        return a - (sums / counts)[inv]

    cols = [demean(age)]
    chip = demean(epic.astype(float))
    if np.any(np.abs(chip) > 1e-12):
        cols.append(chip)
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), demean(v), rcond=None)
    return float(coef[0])


def person_slopes(people: pl.DataFrame, measures: dict[str, np.ndarray], rng: np.random.Generator) -> list[dict]:
    """Within-person and between-person (first draws) age slopes of each measure, per SD of the measure."""
    person = people["person"].to_numpy()
    age = people["age"].to_numpy()
    epic = people["epic"].to_numpy()
    repeated = people.group_by("person").len().filter(pl.col("len") >= 2)["person"]
    rep = np.isin(person, repeated.to_numpy())
    first = pl.DataFrame({"p": person, "a": age}).with_columns(m=pl.col("a").min().over("p"))["m"].to_numpy() == age
    ids = np.unique(person[rep])
    groups = {k: np.flatnonzero(person == k) for k in ids}
    rows = []
    for name, v in measures.items():
        z = (v - v.mean()) / v.std(ddof=1)
        est = within_slope(z[rep], age[rep], person[rep], epic[rep])
        boots = []
        for _ in range(BOOTSTRAPS):
            pick = rng.choice(ids, ids.size, replace=True)
            idx = np.concatenate([groups[k] for k in pick])
            tag = np.repeat(np.arange(pick.size), [groups[k].size for k in pick])
            boots.append(within_slope(z[idx], age[idx], tag, epic[idx]))
        lo, hi = np.quantile(boots, [0.025, 0.975])
        f = first
        design = np.column_stack([np.ones(f.sum()), age[f], epic[f].astype(float)])
        between, *_ = np.linalg.lstsq(design, z[f], rcond=None)
        rows.append(
            {
                "measure": name,
                "within_per_decade": 10 * est,
                "low": 10 * lo,
                "high": 10 * hi,
                "between_per_decade": 10 * float(between[1]),
                "people": int(ids.size),
                "draws": int(rep.sum()),
            }
        )
    return rows


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "biostudies" / STUDY
    tables = [folder / t for t in TABLES]
    if not (folder / "SATSA_pheno.txt").exists() or not all(t.exists() for t in tables):
        print(f"within_person: {STUDY} tables missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    people = read_people(folder / "SATSA_pheno.txt")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    cpgs, x = read_rows(tables, set(features), people["sample"].to_list())
    return analyse(ctx, people, cpgs, x, probes, membership, "within_person", "SATSA")


def read_sdrf(sdrf: Path) -> pl.DataFrame:
    """One row per array of ArrayExpress E-MTAB-7309: ``sample`` (sentrix_position), ``person``, ``age``, ``female``,
    ``epic`` (all 450k)."""
    t = pl.read_csv(sdrf, separator="\t", infer_schema_length=0)
    t = t.filter(pl.col("Label") == "Cy3")
    return (
        t.select(
            pl.col("Assay Name").str.replace(r"_Grn$", "").alias("sample"),
            pl.col("Characteristics[individual]").alias("person"),
            pl.col("Characteristics[age]").cast(pl.Float64).alias("age"),
            (pl.col("Characteristics[sex]") == "female").alias("female"),
        )
        .with_columns(epic=pl.lit(False))
        .unique(subset="sample", keep="first")
        .sort("person", "age")
    )


IDAT_BETAS = "cohorts/satsa_betas_sesame.parquet"
SDRF = "arrayexpress/E-MTAB-7309/E-MTAB-7309.sdrf.txt"


def run_idats(ctx: Context) -> list[Path]:
    """The same analysis on SATSA's 450k waves reprocessed from their IDATs (E-MTAB-7309; run R66a)."""
    path, sdrf = ctx.paths.interim / IDAT_BETAS, ctx.paths.input / SDRF
    if not path.exists() or not sdrf.exists():
        print(f"satsa_idats: {path} or {sdrf} missing (run R66a). Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    have = set(pl.read_parquet_schema(path))
    people = read_sdrf(sdrf).filter(pl.col("sample").is_in(list(have)))
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    frame = (
        pl.scan_parquet(path).select("cpg", *people["sample"].to_list()).filter(pl.col("cpg").is_in(features)).collect()
    )
    x = frame.drop("cpg").to_numpy().T.astype(float)
    return analyse(ctx, people, frame["cpg"].to_list(), x, probes, membership, "satsa_idats", "SATSA from IDATs")


def analyse(ctx, people, cpgs, x, probes, membership, name: str, label: str) -> list[Path]:
    """Pair metrics, gaps, the cross-sectional control and per-person slopes; tables under ``name``."""
    out = ctx.out_dir(name)
    ok = ~np.isnan(x).any(axis=0)
    cpgs = [c for c, k in zip(cpgs, ok, strict=True) if k]
    x = x[:, ok]
    print(f"{name}: {people.height} draws, {people['person'].n_unique()} people, {len(cpgs)} CpGs", file=sys.stderr)
    age, epic, person = people["age"].to_numpy(), people["epic"].to_numpy(), people["person"].to_numpy()
    first, second = consecutive_pairs(person, age)
    m = pair_metrics(x, age, epic, first, second)
    youngest = pl.DataFrame({"p": person, "a": age}).with_columns(m=pl.col("a").min().over("p"))["m"].to_numpy() == age
    metrics = pl.DataFrame(
        {
            "feature": cpgs,
            "baseline": x[youngest].mean(axis=0),
            "abs_directed": np.abs(m["directed"]),
            "divergence": m["divergence"],
            "floor": m["floor"],
        }
    )
    write_table(metrics, out / "cpgs")
    table = _matched_sd(_with_strata(metrics, probes), membership, METRICS, 1000, ctx.seed)
    table = table.with_columns(
        pairs=pl.lit(int(first.size)), median_years=pl.lit(float(np.median(age[second] - age[first])))
    )
    write_table(table, out / "gaps")
    # positive control: the usual cross-sectional comparison on each person's first draw, same preprocessed data
    chip0 = epic[youngest].astype(float)
    cov0 = chip0[:, None] if chip0.std() > 0 else None
    cross = drift_metrics(x[youngest], age[youngest], cpgs, covariates=cov0, spearman=False)
    cross = _matched_sd(_with_strata(cross, probes), membership, LIFESPAN_METRICS, 1000, ctx.seed)
    write_table(cross, out / "cross_sectional")
    print(
        f"{name}: {first.size} pairs; divergence > 0 for {np.mean(m['divergence'] > 0):.2f} of CpGs",
        file=sys.stderr,
    )
    scores = reference_scores(x, _islands(cpgs, probes), _sets(cpgs, membership))
    measures = {f"{k} {g}": scores[f"{k}_{g}"] for g in GROUPS for k in ("B", "F", "S") if f"{k}_{g}" in scores}
    slopes = pl.DataFrame(person_slopes(people, measures, ctx.rng(7)))
    write_table(slopes, out / "person_slopes")
    write_table(people.with_columns(**{k: pl.Series(v) for k, v in measures.items()}), out / "per_draw")
    return [_figure(ctx, table, slopes, name, label)]


def _figure(ctx: Context, table: pl.DataFrame, slopes: pl.DataFrame, name: str, label: str) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.8)
    fig.subplots_adjust(wspace=0.75)
    ax = axes[0]
    groups = [g for g in ("core_essential", "housekeeping", "non_essential") if g in table["group"].unique()]
    for i, group in enumerate(groups):
        t = table.filter(pl.col("group") == group)
        ax.barh(
            np.arange(t.height) + 0.27 * i,
            t["matched_diff_sd"],
            height=0.25,
            color=plots.PALETTE[i],
            label=group.replace("_", " "),
        )
    names = table.filter(pl.col("group") == groups[0])["metric"].to_list()
    ax.set_yticks(np.arange(len(names)) + 0.27, [n.replace("_", " ") for n in names])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.legend(fontsize=7, loc="upper right")
    ax.locator_params(axis="x", nbins=4)
    ax.set(xlabel="set minus look-alikes (SD)", title="Change between a person's draws")
    plots.label_panel(ax, "a")
    ax = axes[1]
    y = np.arange(slopes.height)
    ax.errorbar(
        slopes["within_per_decade"],
        y,
        xerr=[slopes["within_per_decade"] - slopes["low"], slopes["high"] - slopes["within_per_decade"]],
        fmt="o",
        color=plots.PALETTE[0],
        label="within person",
    )
    ax.plot(slopes["between_per_decade"], y + 0.2, "s", color=plots.PALETTE[1], label="between people")
    ax.set_yticks(y, [m.replace("_", " ") for m in slopes["measure"].to_list()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.set(xlabel="SD per decade of age", title="Scores over a person's own aging")
    plots.label_panel(ax, "b")
    fig.suptitle(
        f"The shield within a person over time, {label}", x=0.02, y=1.04, ha="left", fontsize=11, fontweight="bold"
    )
    return plots.save(fig, ctx.figure_path(name))
