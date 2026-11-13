"""Why do human bone marrow and lymph node reverse the shield? (roadmap RM16, step 2)

In the mammalian array (``figures species``, R12) human bone marrow (20 donors) and lymph node (27) are the only
human tissues whose housekeeping entropy gap is positive (+0.18, +0.20 SD); spleen, liver, heart and blood of the
same organ-donor series are negative. The gap's interval there resamples genes, not donors. Four checks, on every
human tissue of that series:

1. **Donors resampled**: the gap's interval from drawing donors with replacement, not genes.
2. **Chip held**: the slopes with each array chip as a fixed effect (marrow sits on two chips).
3. **Same donors**: donors are matched across tissues by age and sex; each tissue is refitted on the donors it
   shares with bone marrow, so a tissue effect is told apart from a donor-set effect.
4. **Small-sample null**: from the larger tissues, ``N_DRAWS`` draws of 20 donors within marrow's age range;
   how often a gap as high as marrow's appears by chance.

And one about what ages in marrow: across all promoter CpGs, the correlation of each tissue's per-CpG methylation
slope with age against the adipose-minus-blood level difference. Marrow turns fatty with age; a slope pattern
that follows the adipose direction says the age signal there is a change of cell mix.

``entropy-age figures marrow_anomaly`` as a Slurm job (reads the human columns of GSE223748).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import probes_without_suffix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.species_shield import layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, drift_metrics

MIN_DONORS = 15
MIN_SPAN = 15.0  # years of age a tissue's donors span
N_BOOT = {"paper": 300, "quick": 40}
N_DRAWS = {"paper": 400, "quick": 40}
ODD = ("BoneMarrow", "LymphNode")
METRICS = ("entropy_slope", "dispersion_slope")


def gaps(x: np.ndarray, ages: np.ndarray, island: np.ndarray, sets: dict, cov: np.ndarray | None = None) -> dict:
    """``{(metric, group): gap in SD}``: set minus look-alikes of the same island status and young level decile."""
    metrics = drift_metrics(x, ages, None, covariates=cov, spearman=False)
    young = ages <= np.quantile(ages, 0.25)
    strata = np.char.add(island.astype(str), baseline_strata(x[young].mean(axis=0), 10).astype(str))
    out = {}
    for metric in METRICS:
        v = metrics[metric].to_numpy()
        for group, (in_set, _) in sets.items():
            sd = float(np.nanstd(v))
            gap = float(matched_scores(v[None, :], in_set, strata)[0][0])
            out[(metric, group)] = gap / sd if sd > 0 else np.nan
    return out


def donor_key(ages: np.ndarray, female: np.ndarray) -> np.ndarray:
    """A donor's identity across tissues of one series: age to 0.1 year and sex."""
    return np.array([f"{a:.1f}|{f}" for a, f in zip(ages, female, strict=True)])


def mix_correlation(slope: np.ndarray, direction: np.ndarray) -> float:
    ok = np.isfinite(slope) & np.isfinite(direction)
    return float(np.corrcoef(slope[ok], direction[ok])[0, 1])


def _human(ctx: Context):
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        return None
    samples = read_samples(matrix).filter(
        (pl.col("species") == "Human")
        & pl.col("age").is_not_null()
        & pl.col("barcode").is_not_null()
        & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
    )
    sex = _characteristics(matrix).select("sample", pl.col("female").cast(pl.String).fill_null("?"))
    samples = samples.join(sex, on="sample", how="left")
    header = pl.read_csv(betas_path, n_rows=0).columns
    samples = samples.filter(pl.col("barcode").is_in(header))
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([header[0], *samples["barcode"].unique().to_list()])
        .filter(pl.col(header[0]).is_in(probes["cpg"].implode()))
        .collect()
        .rename({header[0]: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    return samples, frame, features, isl, members


def run(ctx: Context) -> list[Path]:
    loaded = _human(ctx)
    if loaded is None:
        print("marrow_anomaly: GSE223748 betas missing; run `entropy-age fetch mammals`. Nothing run.", file=sys.stderr)
        return []
    samples, frame, features, isl, members = loaded
    plots.apply_style()
    out = ctx.out_dir("marrow_anomaly")
    rng = ctx.rng(230_000)
    scale = str(ctx.scale.name)
    tissues = {}
    for (tissue,), sel in samples.group_by("tissue", maintain_order=True):
        if sel.height < MIN_DONORS or float(sel["age"].max() - sel["age"].min()) < MIN_SPAN:
            continue
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        tissues[tissue] = (sel, x)
    # CpGs measured in every marrow and lymph-node sample; a tissue missing over a tenth of them is left out
    ok = np.all([~np.isnan(tissues[t][1]).any(axis=0) for t in ODD if t in tissues], axis=0)
    for t in [t for t in tissues if t not in ODD]:
        full = ~np.isnan(tissues[t][1]).any(axis=0)
        if (full & ok).sum() < 0.9 * ok.sum():
            print(f"marrow_anomaly: {t}: {int((ok & ~full).sum())} CpGs missing; left out", file=sys.stderr)
            del tissues[t]
        else:
            ok &= full
    cols = np.flatnonzero(ok)
    sets = {
        g: (s[cols], [np.searchsorted(cols, gi[np.isin(gi, cols)]) for gi in genes if np.isin(gi, cols).any()])
        for g, (s, genes) in _sets(features, members).items()
    }
    isl = isl[cols]
    tissues = {t: (sel, x[:, cols]) for t, (sel, x) in tissues.items()}
    print(f"marrow_anomaly: {len(tissues)} human tissues, {cols.size} CpGs", file=sys.stderr, flush=True)

    marrow_sel = tissues["BoneMarrow"][0]
    marrow_donors = set(donor_key(marrow_sel["age"].to_numpy(), marrow_sel["female"].to_numpy()))
    lo, hi = float(marrow_sel["age"].min()), float(marrow_sel["age"].max())
    direction = None
    if "Adipose" in tissues and "Blood" in tissues:
        direction = tissues["Adipose"][1].mean(axis=0) - tissues["Blood"][1].mean(axis=0)

    rows, null_rows = [], []
    for tissue, (sel, x) in tissues.items():
        age = sel["age"].to_numpy()
        chips = sel["barcode"].str.extract(r"^(\d+)_").to_numpy()
        base = gaps(x, age, isl, sets)
        chip = gaps(x, age, isl, sets, dummies(chips) if np.unique(chips).size > 1 else None)
        boot = []
        for _ in range(N_BOOT[scale]):
            idx = rng.integers(0, age.size, age.size)
            if np.unique(age[idx]).size < 5:
                continue
            boot.append(gaps(x[idx], age[idx], isl, sets))
        keys = donor_key(age, sel["female"].to_numpy())
        shared = np.isin(keys, list(marrow_donors))
        same = gaps(x[shared], age[shared], isl, sets) if shared.sum() >= 12 else {}
        slope = drift_metrics(x, age, None, spearman=False)
        mix = mix_correlation(slope["slope"].to_numpy(), direction) if direction is not None else np.nan
        for key, value in base.items():
            draws = np.array([b[key] for b in boot])
            rows.append(
                {
                    "tissue": tissue,
                    "donors": sel.height,
                    "chips": int(np.unique(chips).size),
                    "metric": key[0],
                    "group": key[1],
                    "gap_sd": value,
                    "donor_low": float(np.nanquantile(draws, 0.025)),
                    "donor_high": float(np.nanquantile(draws, 0.975)),
                    "chip_held": chip[key],
                    "shared_donors": int(shared.sum()),
                    "same_donor_gap": same.get(key, np.nan),
                    "slope_vs_adipose": mix,
                }
            )
        in_range = np.flatnonzero((age >= lo) & (age <= hi))
        if tissue in ODD or in_range.size < 30:
            continue
        for _ in range(N_DRAWS[scale]):
            idx = rng.choice(in_range, size=marrow_sel.height, replace=False)
            g = gaps(x[idx], age[idx], isl, sets)
            null_rows += [{"tissue": tissue, "metric": k[0], "group": k[1], "gap_sd": v} for k, v in g.items()]
        print(f"marrow_anomaly: {tissue} done", file=sys.stderr, flush=True)

    table = pl.DataFrame(rows)
    null = pl.DataFrame(null_rows)
    observed = table.filter(pl.col("tissue").is_in(ODD)).select("tissue", "metric", "group", observed="gap_sd")
    share = (
        null.join(observed, on=["metric", "group"])
        .group_by("tissue_right", "metric", "group")
        .agg(
            draws=pl.len(),
            share_as_high=(pl.col("gap_sd") >= pl.col("observed")).mean(),
            null_mean=pl.col("gap_sd").mean(),
            null_sd=pl.col("gap_sd").std(),
        )
        .rename({"tissue_right": "odd_tissue"})
        .sort("odd_tissue", "metric", "group")
    )
    write_table(table, out / "checks")
    write_table(null, out / "null_draws")
    write_table(share, out / "null_share")
    return [_figure(ctx, table, null)]


def _figure(ctx: Context, table: pl.DataFrame, null: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    hk = table.filter((pl.col("metric") == "entropy_slope") & (pl.col("group") == "housekeeping")).sort("gap_sd")
    y = np.arange(hk.height)
    ax = axes[0]
    colors = [plots.PALETTE[1] if t in ODD else plots.PALETTE[0] for t in hk["tissue"]]
    ax.errorbar(
        hk["gap_sd"],
        y,
        xerr=[hk["gap_sd"] - hk["donor_low"], hk["donor_high"] - hk["gap_sd"]],
        fmt="none",
        ecolor=plots.MUTED,
        lw=1,
    )
    ax.scatter(hk["gap_sd"], y, c=colors, zorder=3, label="all donors")
    ax.scatter(hk["chip_held"], y, marker="|", s=80, color="k", zorder=4, label="chip held")
    ax.scatter(hk["same_donor_gap"], y, marker="x", color=plots.PALETTE[2], zorder=4, label="marrow's donors")
    ax.set_yticks(y, [f"{t} ({n})" for t, n in zip(hk["tissue"], hk["donors"], strict=True)])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set(xlabel="housekeeping entropy gap (SD), donor interval", title="Human tissues, one donor series")
    ax.legend(fontsize=7, loc="lower right")
    plots.label_panel(ax, "a")
    ax = axes[1]
    sub = null.filter((pl.col("metric") == "entropy_slope") & (pl.col("group") == "housekeeping"))
    ax.hist(sub["gap_sd"], bins=40, color=plots.PALETTE[0], alpha=0.8)
    for t in ODD:
        v = hk.filter(pl.col("tissue") == t)["gap_sd"]
        if v.len():
            ax.axvline(float(v[0]), color=plots.PALETTE[1], lw=1.5)
            ax.text(float(v[0]), ax.get_ylim()[1] * 0.9, t, rotation=90, ha="right", va="top", fontsize=7)
    ax.set(xlabel="gap in 20-donor draws of other tissues (SD)", ylabel="draws", title="Small-sample null")
    plots.label_panel(ax, "b")
    fig.subplots_adjust(wspace=0.6)
    return plots.save(fig, ctx.figure_path("marrow_anomaly"))
