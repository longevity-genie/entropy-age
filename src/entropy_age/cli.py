"""Command line interface: ``uv run entropy-age --help``."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

# keep BLAS from saturating every core (the figures are many small fits, not one big one)
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "2")

import polars as pl  # noqa: E402
import typer  # noqa: E402

from entropy_age.data.geo import fetch_series_matrix, prepare_series  # noqa: E402
from entropy_age.data.matrix import FeatureMatrix  # noqa: E402
from entropy_age.data.upstream import UpstreamLayout, fetch_upstream, prepare_upstream  # noqa: E402
from entropy_age.paths import ProjectPaths  # noqa: E402

app = typer.Typer(no_args_is_help=True, help="Aging clocks from accumulating stochastic variation.")
fetch_app = typer.Typer(no_args_is_help=True, help="Download and prepare input data.")
app.add_typer(fetch_app, name="fetch")


@fetch_app.command("upstream")
def fetch_upstream_cmd(
    from_dir: Annotated[
        Path | None, typer.Option(help="Copy from a local clone's Data/ folder instead of GitHub")
    ] = None,
    overwrite: Annotated[bool, typer.Option(help="Re-download files that already exist")] = False,
) -> None:
    """Clock coefficients, sample sheets and the dolphin methylation example from the upstream repo."""
    paths = ProjectPaths.default().ensure()
    files = fetch_upstream(paths, from_dir=from_dir, overwrite=overwrite)
    layout = prepare_upstream(paths)
    typer.echo(f"{len(files)} input files in {layout.input}\nparquet in {layout.interim}")


@fetch_app.command("geo")
def fetch_geo_cmd(
    accession: Annotated[str, typer.Argument(help="GEO series, e.g. GSE41037")] = "GSE41037",
    overwrite: bool = False,
) -> None:
    """A GEO series matrix (beta values + sample characteristics), saved as parquet."""
    paths = ProjectPaths.default().ensure()
    raw = fetch_series_matrix(accession, paths, overwrite=overwrite)
    samples, betas = prepare_series(accession, paths)
    typer.echo(f"downloaded {raw}\nsamples -> {samples}\nbetas   -> {betas}")


@fetch_app.command("annotation")
def fetch_annotation_cmd(overwrite: bool = False) -> None:
    """Illumina 27k probe manifest (GEO GPL8490) and curated gene sets (essential, non-essential, housekeeping)."""
    from entropy_age.data.annotation import fetch_annotation, prepare_annotation

    paths = ProjectPaths.default().ensure()
    files = fetch_annotation(paths, overwrite=overwrite)
    layout = prepare_annotation(paths)
    typer.echo(
        f"{len(files)} input files in {layout.input}\nprobes    -> {layout.probes}\ngene sets -> {layout.gene_sets}"
    )


@app.command()
def figures(
    names: Annotated[
        list[str] | None, typer.Argument(help="fig1 fig3 fig4 fig5 extension atlas genesets (default: all)")
    ] = None,
    scale: Annotated[str, typer.Option(help="quick (laptop) or paper (publication sizes)")] = "quick",
    seed: int = 42,
) -> None:
    """Regenerate figures into docs/figures and tables into output/."""
    from entropy_age.experiments import FIGURES, Context, Scale

    ctx = Context(ProjectPaths.default().ensure(), Scale.of(scale), seed)
    for name in names or list(FIGURES):
        if name not in FIGURES:
            raise typer.BadParameter(f"unknown figure {name!r}; choose from {', '.join(FIGURES)}")
        typer.echo(f"[{name}] running at scale={scale} seed={seed}")
        for path in FIGURES[name](ctx):
            typer.echo(f"[{name}] wrote {path}")


@app.command()
def stochasticity(
    betas: Annotated[Path, typer.Argument(help="Parquet matrix frame: feature column + one column per sample")],
    samples: Annotated[Path, typer.Argument(help="Parquet/CSV with sample and age columns")],
    out: Annotated[Path, typer.Option(help="Where to write per-feature metrics (parquet)")] = Path(
        "stochasticity.parquet"
    ),
    sample_col: str = "sample",
    age_col: str = "age",
    feature_col: Annotated[str | None, typer.Option(help="Defaults to the first column of BETAS")] = None,
    groups: Annotated[Path | None, typer.Option(help="Parquet/CSV with feature and group columns")] = None,
    metric: Annotated[str, typer.Option(help="Metric compared between groups")] = "entropy_slope",
) -> None:
    """Per-feature drift, entropy, dispersion and relaxation-rate estimates; optional gene-set contrasts."""
    from entropy_age.stochasticity import compare_groups, drift_metrics, fit_relaxation

    def read(path: Path) -> pl.DataFrame:
        return pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(path)

    table = read(samples).select(pl.col(sample_col).cast(pl.String), pl.col(age_col).cast(pl.Float64))
    frame = pl.read_parquet(betas)
    table = table.filter(pl.col(sample_col).is_in(frame.columns)).drop_nulls()
    matrix = FeatureMatrix.from_frame(frame, feature_col=feature_col, samples=table[sample_col].to_list())
    ages = table[age_col].to_numpy()
    metrics = fit_relaxation(matrix.values, ages, matrix.features).join(
        drift_metrics(matrix.values, ages, matrix.features), on="feature"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    metrics.write_parquet(out)
    typer.echo(f"{metrics.height} features x {len(matrix.samples)} samples -> {out}")
    if groups is not None:
        contrasts = compare_groups(metrics, read(groups).select("feature", "group"), metric=metric)
        contrasts_path = out.with_name(out.stem + f"_groups_{metric}.parquet")
        contrasts.write_parquet(contrasts_path)
        typer.echo(contrasts)


@app.command()
def info() -> None:
    """Show where data lives and what is already prepared."""
    paths = ProjectPaths.default()
    layout = UpstreamLayout(paths)
    checks = {
        "upstream clocks": layout.clock("horvath"),
        "dolphin methylation": layout.dolphin_betas,
        "GSE41037 sample sheet": layout.gse41037_samples,
        "GSE41037 beta values": paths.interim / "geo" / "gse41037_betas.parquet",
        "27k probe annotation + gene sets": paths.interim / "annotation" / "gene_sets.parquet",
    }
    typer.echo(f"root    {paths.root}\ndata    {paths.data}\nfigures {paths.figures}")
    for label, path in checks.items():
        typer.echo(f"{'ok ' if path.exists() else '-- '} {label}: {path}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
