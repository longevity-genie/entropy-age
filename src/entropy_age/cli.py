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


@fetch_app.command("chromatin")
def fetch_chromatin_cmd(overwrite: bool = False) -> None:
    """Roadmap Epigenomics chromHMM 15-state segmentations of blood and brain epigenomes (~10 MB each)."""
    from entropy_age.data.chromatin import ChromatinLayout, fetch_chromatin

    paths = ProjectPaths.default().ensure()
    files = fetch_chromatin(paths, overwrite=overwrite)
    typer.echo("\n".join(str(f) for f in files) + f"\n-> {ChromatinLayout(paths).segmentation('E062').parent}")


@fetch_app.command("clocks")
def fetch_clocks_cmd(overwrite: bool = False) -> None:
    """Published clocks' CpGs and weights (biolearn's data at a pinned commit, each file hash-checked)."""
    from entropy_age.data.clocks import fetch_clocks

    paths = ProjectPaths.default().ensure()
    files = fetch_clocks(paths, overwrite=overwrite)
    typer.echo(f"{len(files)} clocks in {files[0].parent}")


@fetch_app.command("rrbs")
def fetch_rrbs_cmd(overwrite: bool = False) -> None:
    """Mouse blood RRBS of GSE80672 (about 3 GB) and mm10 promoter and CpG-island tracks, parsed to promoter CpGs."""
    from entropy_age.data.rrbs import fetch_rrbs, prepare_rrbs

    paths = ProjectPaths.default().ensure()
    files = fetch_rrbs(paths, overwrite=overwrite)
    layout = prepare_rrbs(paths)
    typer.echo(f"{len(files)} files\nsamples -> {layout.samples}\nbetas   -> {layout.betas}")


@fetch_app.command("mouse-annotation")
def fetch_mouse_annotation_cmd(overwrite: bool = False) -> None:
    """Mammalian-array probe annotation (GEO GPL32543, ~200 MB) and MGI mouse-human orthologs."""
    from entropy_age.data.mouse import fetch_mouse_annotation, prepare_mouse_annotation

    paths = ProjectPaths.default().ensure()
    files = fetch_mouse_annotation(paths, overwrite=overwrite)
    layout = prepare_mouse_annotation(paths)
    typer.echo(
        f"{len(files)} input files in {layout.input}\nprobes    -> {layout.probes}\ngene sets -> {layout.gene_sets}"
    )


@fetch_app.command("rodent-rna")
def fetch_rodent_rna_cmd(overwrite: bool = False) -> None:
    """Mouse and rat aging RNA atlases (GSE314817, GSE314970): series matrices and raw counts, ~80 MB each."""
    from entropy_age.data.download import download
    from entropy_age.data.geo import series_matrix_url
    from entropy_age.experiments.rna_rodents import SPECIES

    paths = ProjectPaths.default().ensure()
    for gse, counts in SPECIES.values():
        folder = paths.input / "geo" / gse
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/{counts}"
        for target, source in (
            (folder / f"{gse}_series_matrix.txt.gz", series_matrix_url(gse)),
            (folder / counts, url),
        ):
            download(source, target, overwrite=overwrite, timeout=3600)
            typer.echo(f"{gse}: {target}")


@fetch_app.command("cohort")
def fetch_cohort_cmd(
    names: Annotated[list[str], typer.Argument(help="hannum jaffe lunnon huang gtex fibres endicott sturm")],
    overwrite: bool = False,
) -> None:
    """450k methylation cohorts for the cell-type test (GB-scale downloads; run on a cluster)."""
    from entropy_age.data.cohorts import fetch_cohort, prepare_cohort

    paths = ProjectPaths.default().ensure()
    for name in names:
        files = fetch_cohort(name, paths, overwrite=overwrite)
        lay = prepare_cohort(name, paths)
        typer.echo(f"[{name}] {len(files)} input files; samples -> {lay.samples}; betas -> {lay.betas}")


@app.command()
def figures(
    names: Annotated[
        list[str] | None, typer.Argument(help="figure names, e.g. fig1 tissues gtex rna_gtex (default: all)")
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
