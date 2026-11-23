"""Protocol P3: one table and one forest plot of every series' housekeeping gap, with provenance.

Rows: the relative-entropy slope gap of housekeeping promoter CpGs against matched CpGs, all CpGs, with the interval
that resamples donors, animals or cultures. Shared-program series use the chip-held fit where it is identifiable,
else the fit without chip (column ``fit``). One row per independent cohort and tissue. No pooled estimate and no
moderator claim (protocol P3). ``uv run python scripts/forest_series.py`` writes docs/forest_series.csv and
docs/figures/forest_series.png.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "output"
SERIES_INFO = {  # name -> (species, platform, processing)
    "brain_regions": ("human", "EPIC", "rebuilt"),
    "neurons": ("human", "EPIC", "rebuilt"),
    "colon_mucosa": ("human", "450k", "rebuilt"),
    "lung_450k": ("human", "450k", "rebuilt"),
    "lung_epic": ("human", "EPIC", "rebuilt"),
    "lung_healthy": ("human", "EPIC", "rebuilt"),
    "breast_komen": ("human", "450k", "rebuilt"),
    "breast_cancer_free": ("human", "450k", "rebuilt"),
    "kidney_normal": ("human", "450k", "rebuilt"),
    "myometrium_epic": ("human", "EPIC", "rebuilt"),
    "myometrium_epicv2": ("human", "EPIC v2", "rebuilt"),
    "multi_tissue_mammal40": ("human", "mammalian 40k", "rebuilt"),
    "mouse_hearts": ("mouse", "MM285", "rebuilt"),
    "mouse_liver_field": ("mouse", "MM285", "rebuilt"),
    "mouse_bxd_liver_control": ("mouse", "mammalian 40k", "rebuilt"),
}
SKIP_GROUPS = {
    "neurons, controls",
    "myometrium, fibroid-free",
    "liver, field",
    "hearts, saline",
    "prefrontal cortex, controls",
    "amygdala, controls",
    "hippocampus, controls",
}  # subsets of a row


def runner_rows() -> list[dict]:
    rows = []
    for name, (species, platform, processing) in SERIES_INFO.items():
        path = OUT / f"series_{name}" / "gaps.csv"
        if not path.exists():
            continue
        g = pl.read_csv(path).filter(
            (pl.col("metric") == "kl_slope") & (pl.col("group") == "housekeeping") & (pl.col("band") == "all")
        )
        sizes = pl.read_csv(OUT / f"series_{name}" / "sizes.csv")
        spec = tomllib.loads((ROOT / "config" / "series" / f"{name}.toml").read_text())
        for cohort in g["cohort"].unique(maintain_order=True).to_list():
            if cohort in SKIP_GROUPS:
                continue
            sub = g.filter(pl.col("cohort") == cohort)
            chip = sub.filter((pl.col("covariates") == "held + chip") & pl.col("donor_ci_low").is_not_nan())
            r = (chip if chip.height else sub.filter(pl.col("covariates") == "held")).row(0, named=True)
            size = sizes.filter(pl.col("group") == cohort).row(0, named=True)
            label = "lung next to adenocarcinoma" if name.startswith("lung_") and name != "lung_healthy" else cohort
            rows.append(
                {
                    "series": f"{label.replace('_', ' ')} ({spec['accession']})",
                    "species": species,
                    "platform": platform,
                    "processing": processing,
                    "n": size["donors"],
                    "age_min": size["age_min"],
                    "age_max": size["age_max"],
                    "fit": r["covariates"],
                    "gap": r["gap_sd"],
                    "ci_low": r["donor_ci_low"],
                    "ci_high": r["donor_ci_high"],
                    "interval": "donors or animals",
                }
            )
    return rows


def headline_rows() -> list[dict]:
    """Hannum blood and Jaffe cortex: donor bootstrap scaled by the full sample's SD."""
    boot = pl.read_parquet(OUT / "tissue_contrast" / "bootstrap_lifespan.parquet").filter(
        (pl.col("metric") == "kl_slope") & (pl.col("group") == "housekeeping")
    )
    point = pl.read_parquet(OUT / "tissue_contrast" / "point_estimates.parquet").filter(
        (pl.col("metric") == "kl_slope") & (pl.col("group") == "housekeeping") & (pl.col("part") == "lifespan")
    )
    info = {"blood": ("whole blood (Hannum, GSE40279)", 656), "cortex": ("prefrontal cortex (Jaffe, GSE74193)", 213)}
    rows = []
    for r in point.iter_rows(named=True):
        v = boot.filter(pl.col("cohort") == r["tissue"])["matched_diff"].to_numpy() / r["sd"]
        label, n = info[r["tissue"]]
        rows.append(
            {
                "series": label,
                "species": "human",
                "platform": "450k",
                "processing": "deposited" if r["tissue"] == "blood" else "rebuilt",
                "n": n,
                "age_min": None,
                "age_max": None,
                "fit": "held",
                "gap": r["matched_diff_sd"],
                "ci_low": float(np.quantile(v, 0.025)),
                "ci_high": float(np.quantile(v, 0.975)),
                "interval": "donors",
            }
        )
    return rows


def module_rows() -> list[dict]:
    """Sequencing series whose modules resample animals or donors."""
    rows = []
    for folder, table, label, species, platform in (
        ("jena_rrbs", "gaps.csv", "mouse {tissue} (Jena RRBS)", "mouse", "RRBS"),
        ("mouse_atlas_rrbs", "gaps.csv", "mouse {tissue} (GSE225166)", "mouse", "RRBS"),
        ("mouse_rrbs_tissues", "aging_gaps.csv", "mouse {tissue} (Thompson RRBS)", "mouse", "RRBS"),
        ("meer_tissues", "aging_gaps.csv", "mouse {tissue} (Meer RRBS)", "mouse", "RRBS"),
    ):
        path = OUT / folder / table
        if not path.exists():
            continue
        g = pl.read_csv(path).filter((pl.col("metric") == "kl_slope") & (pl.col("group") == "housekeeping"))
        for r in g.iter_rows(named=True):
            rows.append(
                {
                    "series": label.format(tissue=r["tissue"]),
                    "species": species,
                    "platform": platform,
                    "processing": "deposited counts",
                    "n": r.get("mice"),
                    "age_min": None,
                    "age_max": None,
                    "fit": "held",
                    "gap": r["gap_sd"],
                    "ci_low": r["ci_low"],
                    "ci_high": r["ci_high"],
                    "interval": "animals",
                }
            )
    return rows


def main() -> None:
    table = pl.DataFrame([*headline_rows(), *runner_rows(), *module_rows()], infer_schema_length=None).sort(
        "species", "gap", descending=[True, False]
    )
    table.write_csv(ROOT / "docs" / "forest_series.csv")
    plots.apply_style()
    fig, ax = plots.plt.subplots(figsize=(7.2, 0.24 * table.height + 1.2))
    y = np.arange(table.height)
    resolved = (table["ci_high"] < 0) | (table["ci_low"] > 0)
    colours = [plots.PALETTE[0] if s == "human" else plots.PALETTE[1] for s in table["species"]]
    for i, r in enumerate(table.iter_rows(named=True)):
        lo, hi = r["ci_low"], r["ci_high"]
        if lo is not None and np.isfinite(lo):
            ax.plot([lo, hi], [i, i], color=colours[i], lw=1.2)
        ax.scatter([r["gap"]], [i], color=colours[i], s=22 if resolved[i] else 10, zorder=3)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(y, [f"{s} · {p}" for s, p in zip(table["series"], table["platform"], strict=True)], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel(
        "Housekeeping minus matched CpGs (SD of the relative-entropy slope), 95 % interval over donors/animals"
    )
    ax.set_title("Every series, one measure (blue human, orange mouse; large dots: interval excludes 0)", fontsize=9)
    fig.savefig(ROOT / "docs" / "figures" / "forest_series.png")
    print(f"{table.height} rows; {int(resolved.sum())} resolved, {int((table['gap'] < 0).sum())} negative")


if __name__ == "__main__":
    main()
