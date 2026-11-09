"""GEO series-matrix download and parsing (beta values plus per-sample metadata)."""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

_SERIES_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz"
_KEY_VALUE = re.compile(r"^\s*([^:]+?)\s*:\s*(.*)$")


def series_matrix_url(gse: str) -> str:
    gse = gse.upper()
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError(f"not a GEO series accession: {gse!r}")
    return _SERIES_URL.format(prefix=gse[:-3], gse=gse)


def fetch_series_matrix(gse: str, paths: ProjectPaths, *, overwrite: bool = False) -> Path:
    gse = gse.upper()
    target = paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
    return download(series_matrix_url(gse), target, overwrite=overwrite, timeout=600)


@dataclass(frozen=True, slots=True)
class SeriesMatrix:
    samples: pl.DataFrame
    """One row per GSM: title, source and every characteristics line."""
    values: pl.DataFrame
    """Wide matrix frame: ``feature`` column (probe ID) plus one Float32 column per GSM."""


def _split(line: str) -> list[str]:
    return [cell.strip().strip('"') for cell in line.rstrip("\n").split("\t")]


def _clean_key(raw: str) -> str:
    """``age (in years)`` -> ``age``: the text before any parenthesis, lower case, spaces as underscores."""
    return raw.split("(")[0].strip().lower().replace(" ", "_") or raw.strip().lower().replace(" ", "_")


def _characteristic_columns(rows: list[list[str]]) -> dict[str, list[str | None]]:
    """``!Sample_characteristics_ch1`` rows; ``key: value`` cells become a ``key`` column."""
    columns: dict[str, list[str | None]] = {}
    for i, row in enumerate(rows, start=1):
        matches = [_KEY_VALUE.match(cell) for cell in row]
        keys = {_clean_key(m.group(1)) for m in matches if m}
        if len(keys) == 1 and all(m or not cell for m, cell in zip(matches, row, strict=True)):
            name = keys.pop()
            values = [m.group(2) if m else None for m in matches]
        else:
            name = f"characteristic_{i}"
            values = [cell or None for cell in row]
        while name in columns:
            name = f"{name}_{i}"
        columns[name] = values
    return columns


def read_series_matrix(path: Path) -> SeriesMatrix:
    opener = gzip.open if str(path).endswith(".gz") else open
    meta: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    table_lines: list[str] = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if in_table:
                if line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue
                table_lines.append(line)
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!Sample_"):
                key, *cells = _split(line)
                if key == "!Sample_characteristics_ch1":
                    characteristics.append(cells)
                else:
                    meta.setdefault(key.removeprefix("!Sample_"), cells)

    sample_ids = meta["geo_accession"]
    samples = pl.DataFrame(
        {
            "sample": sample_ids,
            **{
                k: v
                for k, v in meta.items()
                if k in ("title", "source_name_ch1", "platform_id") and len(v) == len(sample_ids)
            },
            **_characteristic_columns(characteristics),
        }
    )
    values = pl.read_csv(
        io.StringIO("".join(table_lines)),
        separator="\t",
        quote_char='"',
        null_values=["", "NA", "null", "NaN"],
        infer_schema_length=0,
    )
    values = values.rename({values.columns[0]: "feature"}).with_columns(
        pl.exclude("feature").cast(pl.Float32, strict=False)
    )
    return SeriesMatrix(samples=samples, values=values)


def prepare_series(gse: str, paths: ProjectPaths) -> tuple[Path, Path]:
    """Parse a downloaded series matrix into ``{gse}_samples.parquet`` and ``{gse}_betas.parquet``."""
    gse = gse.upper()
    source = paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
    parsed = read_series_matrix(source)
    out_dir = paths.interim / "geo"
    out_dir.mkdir(parents=True, exist_ok=True)
    samples_path = out_dir / f"{gse.lower()}_samples.parquet"
    betas_path = out_dir / f"{gse.lower()}_betas.parquet"
    parsed.samples.write_parquet(samples_path)
    parsed.values.rename({"feature": "cpg"}).write_parquet(betas_path)
    return samples_path, betas_path


def read_supplementary_betas(
    path: Path,
    *,
    suffix: str | None = None,
    id_column: int = 0,
    preamble: int = 0,
    after_header: int = 0,
) -> pl.DataFrame:
    """A probe-by-sample beta table from a GEO supplementary CSV/TSV (gzip allowed).

    Returns a matrix frame: ``cpg`` plus one Float32 column per sample. With ``suffix`` only
    columns ending in it are kept and the suffix is stripped (GSE74193 interleaves
    ``SampleN_Beta`` with ``SampleN_detectP``). ``preamble`` lines are skipped before the header
    and ``after_header`` rows after it (GSE59685 has four comment lines, a barcode header and two
    rows of GEO accessions). The separator is sniffed from the header.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for _ in range(preamble):
            handle.readline()
        header = handle.readline()
    separator = "\t" if header.count("\t") > header.count(",") else ","
    names = [cell.strip().strip('"') for cell in header.rstrip("\n").split(separator)]
    id_name = names[id_column]
    keep = [n for n in names if n != id_name and (suffix is None or n.endswith(suffix))]
    frame = pl.read_csv(
        path,
        separator=separator,
        skip_rows=preamble,
        skip_rows_after_header=after_header,
        columns=[id_name, *keep],
        schema_overrides={n: pl.Float32 for n in keep} | {id_name: pl.String},
        null_values=["", "NA", "NaN", "null"],
    )
    renamed = {id_name: "cpg"} | ({n: n.removesuffix(suffix) for n in keep} if suffix else {})
    return frame.rename(renamed)
