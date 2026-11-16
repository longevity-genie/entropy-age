"""The RNA turnover test in a third rodent cohort: Tabula Muris Senis bulk RNA-seq (GSE132040; roadmap RM4).

Schaum et al. 2020 profiled 17 organs of C57BL/6JN mice at ten ages from 1 to 27 months, about 54 animals per
organ, by bulk RNA-seq (counts on GENCODE vM19). The measure, gene mapping, gene sets, matching, adjustments and
resampling are those of ``rna_rodents`` (``analyse_species``); only the reader and the turnover classes differ.

**Turnover classes:** renewing - bone marrow, spleen, small intestine, skin, white blood cells; slow - liver,
kidney, pancreas, lung, bone and the four fat depots (BAT, gonadal, mesenteric, subcutaneous); heart & muscle -
heart, limb muscle; brain - whole brain. Animals younger than 3 months are left out, so that every slope is one
of adult aging, not of growth.

``entropy-age figures rna_tms`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.rna_gtex import contrast_panel
from entropy_age.experiments.rna_rodents import analyse_species
from entropy_age.experiments.tissue_contrast import _settings

GSE = "GSE132040"
COUNTS = "GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz"
METADATA = "GSE132040_MACA_Bulk_metadata.csv.gz"
MIN_AGE = 3.0
TURNOVER: dict[str, str] = (
    dict.fromkeys(("Marrow", "Spleen", "Small_Intestine", "Skin", "WBC"), "renewing")
    | dict.fromkeys(("Liver", "Kidney", "Pancreas", "Lung", "Bone", "BAT", "GAT", "MAT", "SCAT"), "slow")
    | dict.fromkeys(("Heart", "Limb_Muscle"), "heart & muscle")
    | {"Brain": "brain"}
)


def read_tms(metadata: Path, counts: Path) -> tuple[pl.DataFrame, list[str], np.ndarray]:
    """Samples (``sample``, ``tissue``, ``sex``, ``age`` in months), gene symbols and a samples x genes count matrix."""
    meta = pl.read_csv(metadata, infer_schema_length=0)
    samples = (
        meta.select(
            pl.col("Sample name").alias("sample"),
            pl.col("source name").str.replace(r"_\d+$", "").alias("tissue"),
            pl.col("characteristics: sex").str.to_lowercase().alias("sex"),
            pl.col("characteristics: age").cast(pl.Float64, strict=False).alias("age"),
        )
        .drop_nulls(["age", "tissue"])
        .filter(pl.col("tissue").is_in(list(TURNOVER)) & (pl.col("age") >= MIN_AGE))
    )
    frame = pl.read_csv(counts, infer_schema_length=0)
    genes = frame[frame.columns[0]].to_list()
    by_sample = {c.split(".")[0]: c for c in frame.columns[1:]}
    samples = samples.filter(pl.col("sample").is_in(list(by_sample)))
    x = frame.select([by_sample[s] for s in samples["sample"].to_list()]).cast(pl.Float64).to_numpy().T
    return samples, genes, x


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE
    annotation = AnnotationLayout(ctx.paths)
    needed = [folder / COUNTS, folder / METADATA, annotation.gene_sets, annotation.hgnc]
    if any(not p.exists() for p in needed):
        print(f"rna_tms: inputs missing ({', '.join(str(p) for p in needed if not p.exists())})", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("rna_tms")
    workers, n_boot = _settings(ctx)
    hgnc, gene_sets = pl.read_parquet(annotation.hgnc), pl.read_parquet(annotation.gene_sets)
    samples, genes, x = read_tms(folder / METADATA, folder / COUNTS)
    tissues, genes_row, point, boot, contrasts, intervals = analyse_species(
        "mouse (TMS)", samples, genes, x, hgnc, gene_sets, n_boot, workers, turnover=TURNOVER
    )
    write_table(tissues, out / "tissues")
    write_table(pl.DataFrame([genes_row], infer_schema_length=None), out / "genes")
    write_table(point, out / "point_estimates")
    boot.write_parquet(out / "bootstrap.parquet")
    write_table(contrasts, out / "contrasts")
    write_table(intervals, out / "tissue_intervals")
    return [_figure(ctx, contrasts)]


def _figure(ctx: Context, contrasts: pl.DataFrame) -> Path:
    fig, axes = plots.panels(1, ncols=1, size=5.0)
    contrast_panel(axes[0], contrasts.filter(pl.col("adjustment") == "covariates"), "Tabula Muris Senis bulk")
    return plots.save(fig, ctx.figure_path("rna_tms"))
