"""Probe-to-gene annotation for the Illumina 27k array and curated gene sets.

* Illumina HumanMethylation27 manifest (GEO GPL8490 supplementary file): gene,
  Entrez ID, RefSeq accession, distance to TSS and CpG-island flag per probe.
* Hart et al. core-essential (CEGv2, 684 genes) and non-essential (NEGv1, 927
  genes) reference sets from CRISPR screens, pinned to a commit of hart-lab/bagel (MIT).
* Eisenberg & Levanon (2013) human housekeeping genes (3,804 genes).

Hart sets are matched to probes by Entrez ID and the housekeeping set by RefSeq
accession (falling back to the symbol), because the 2007 manifest uses old symbols.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

MANIFEST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8490/suppl/GPL8490_HumanMethylation27_270596_v.1.2.csv.gz"
)
HART_COMMIT = "53388adbb4fb0931e5c9dda135502be19e4555f0"
HART_URL = f"https://raw.githubusercontent.com/hart-lab/bagel/{HART_COMMIT}/{{name}}"
HOUSEKEEPING_URL = "https://www.tau.ac.il/~elieis/HKG/HK_genes.txt"

SOURCES: dict[str, str] = {
    "GPL8490_HumanMethylation27_270596_v.1.2.csv.gz": MANIFEST_URL,
    "CEGv2.txt": HART_URL.format(name="CEGv2.txt"),
    "NEGv1.txt": HART_URL.format(name="NEGv1.txt"),
    "HK_genes.txt": HOUSEKEEPING_URL,
}

#: gene-set name -> (file, description)
GENE_SETS: dict[str, tuple[str, str]] = {
    "core_essential": ("CEGv2.txt", "Hart et al. 2017 core-essential genes (CRISPR screens)"),
    "non_essential": ("NEGv1.txt", "Hart et al. 2014 non-essential genes (not expressed in cell lines)"),
    "housekeeping": ("HK_genes.txt", "Eisenberg & Levanon 2013 housekeeping genes"),
}


@dataclass(frozen=True, slots=True)
class AnnotationLayout:
    paths: ProjectPaths

    @property
    def input(self) -> Path:
        return self.paths.input / "annotation"

    @property
    def interim(self) -> Path:
        return self.paths.interim / "annotation"

    @property
    def probes(self) -> Path:
        return self.interim / "hm27_probes.parquet"

    @property
    def gene_sets(self) -> Path:
        return self.interim / "gene_sets.parquet"


def fetch_annotation(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    layout = AnnotationLayout(paths)
    return [download(url, layout.input / name, overwrite=overwrite) for name, url in SOURCES.items()]


def read_manifest(path: Path) -> pl.DataFrame:
    """The ``[Assay]`` block of an Illumina manifest as one row per probe."""
    with gzip.open(path, "rt", encoding="latin-1") if str(path).endswith(".gz") else open(path) as handle:
        lines = handle.read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("[Assay]")) + 1
    end = next((i for i, line in enumerate(lines) if line.startswith("[Controls]")), len(lines))
    table = pl.read_csv(io.StringIO("\n".join(lines[start:end])), infer_schema_length=0, truncate_ragged_lines=True)
    return table.select(
        pl.col("Name").alias("cpg"),
        pl.col("Chr").alias("chr"),
        pl.col("Symbol").str.strip_chars().alias("symbol"),
        pl.col("Gene_ID").str.strip_prefix("GeneID:").cast(pl.Int64, strict=False).alias("entrez"),
        pl.col("Accession").str.replace(r"\.\d+$", "").alias("refseq"),
        pl.col("Distance_to_TSS").cast(pl.Int64, strict=False).alias("distance_to_tss"),
        (pl.col("CPG_ISLAND").str.to_uppercase() == "TRUE").alias("cpg_island"),
    )


def read_gene_sets(input_dir: Path) -> pl.DataFrame:
    """``gene_set, symbol, entrez, refseq`` for every curated set."""
    frames = []
    for name, (file, _) in GENE_SETS.items():
        path = input_dir / file
        if file == "HK_genes.txt":
            rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
            frame = pl.DataFrame(
                {"symbol": [r[0] for r in rows], "refseq": [r[1] if len(r) > 1 else None for r in rows]},
                schema={"symbol": pl.String, "refseq": pl.String},
            ).with_columns(entrez=pl.lit(None, dtype=pl.Int64))
        else:
            frame = pl.read_csv(path, separator="\t", infer_schema_length=0).select(
                pl.col("GENE").alias("symbol"),
                pl.col("ENTREZ_ID").cast(pl.Int64, strict=False).alias("entrez"),
                pl.lit(None, dtype=pl.String).alias("refseq"),
            )
        frames.append(frame.with_columns(gene_set=pl.lit(name)).select("gene_set", "symbol", "entrez", "refseq"))
    return pl.concat(frames)


def prepare_annotation(paths: ProjectPaths) -> AnnotationLayout:
    layout = AnnotationLayout(paths)
    layout.interim.mkdir(parents=True, exist_ok=True)
    read_manifest(layout.input / "GPL8490_HumanMethylation27_270596_v.1.2.csv.gz").write_parquet(layout.probes)
    read_gene_sets(layout.input).write_parquet(layout.gene_sets)
    return layout


def probe_membership(probes: pl.DataFrame, gene_sets: pl.DataFrame) -> pl.DataFrame:
    """``feature`` (CpG), ``group`` (gene set) and ``gene``: Entrez first, then RefSeq, then symbol."""
    keys = (("entrez", "entrez"), ("refseq", "refseq"), ("symbol", "symbol"))
    matched = []
    for set_key, probe_key in keys:
        genes = gene_sets.filter(pl.col(set_key).is_not_null()).select("gene_set", pl.col(set_key).alias("key"))
        hits = (
            probes.filter(pl.col(probe_key).is_not_null())
            .select("cpg", "symbol", pl.col(probe_key).alias("key"))
            .join(genes, on="key")
        )
        matched.append(
            hits.select(
                pl.col("cpg").alias("feature"), pl.col("gene_set").alias("group"), pl.col("symbol").alias("gene")
            )
        )
    return pl.concat(matched).unique(subset=["feature", "group"], keep="first", maintain_order=True)
