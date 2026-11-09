"""Mouse annotation for the mammalian methylation array, and the human gene sets carried to mouse.

* GEO platform GPL32543 (HorvathMammalMethylChip320: the mm285k mouse array plus the 40k
  mammalian array), whose table maps every probe to mm10 with ChIPseeker-style annotation:
  gene symbol, Entrez ID, a promoter flag and a CpG-island flag. GEO serves the table from its
  query page (about 200 MB); the FTP copy is a 1.3 GB family file with every sample in it.
* MGI's mouse-human homology report (``HOM_MouseHumanSequence.rpt``), which groups orthologous
  genes under one class key. MGI updates it in place, so a fresh download can differ slightly.

The human essential, non-essential and housekeeping sets are matched to the report by Entrez ID
(then symbol) and carried to every mouse gene of the same class.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

PLATFORM = "GPL32543"
PLATFORM_FILE = f"{PLATFORM}_table.txt"
PLATFORM_URL = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={PLATFORM}&targ=self&form=text&view=data"
ORTHOLOG_FILE = "HOM_MouseHumanSequence.rpt"
ORTHOLOG_URL = f"https://www.informatics.jax.org/downloads/reports/{ORTHOLOG_FILE}"


@dataclass(frozen=True, slots=True)
class MouseAnnotationLayout:
    paths: ProjectPaths

    @property
    def input(self) -> Path:
        return self.paths.input / "annotation" / "mouse"

    @property
    def interim(self) -> Path:
        return self.paths.interim / "annotation" / "mouse"

    @property
    def probes(self) -> Path:
        return self.interim / "mammal320_probes.parquet"

    @property
    def gene_sets(self) -> Path:
        return self.interim / "gene_sets.parquet"


def fetch_mouse_annotation(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    lay = MouseAnnotationLayout(paths)
    return [
        download(PLATFORM_URL, lay.input / PLATFORM_FILE, overwrite=overwrite, timeout=1800),
        download(ORTHOLOG_URL, lay.input / ORTHOLOG_FILE, overwrite=overwrite, timeout=600),
    ]


def read_platform_table(path: Path) -> pl.DataFrame:
    """``cpg, symbol, entrez, chr, promoter, cpg_island`` for every probe that maps to mm10."""
    with open(path) as handle:
        skip = next(i for i, line in enumerate(handle) if line.startswith("ID\t"))
    frame = pl.read_csv(
        path,
        separator="\t",
        skip_rows=skip,
        quote_char='"',
        infer_schema_length=0,
        columns=["ID", "SYMBOL", "ENTREZID", "mm10.Chr", "Promoter", "CpG.island"],
        truncate_ragged_lines=True,
    )
    null = ("", "NA")
    return (
        frame.filter(~pl.col("ID").str.starts_with("!"))
        .select(
            pl.col("ID").alias("cpg"),
            pl.when(pl.col("SYMBOL").is_in(null)).then(None).otherwise(pl.col("SYMBOL")).alias("symbol"),
            pl.col("ENTREZID").cast(pl.Int64, strict=False).alias("entrez"),
            pl.when(pl.col("mm10.Chr").is_in(null)).then(None).otherwise(pl.col("mm10.Chr")).alias("chr"),
            (pl.col("Promoter").str.to_uppercase() == "TRUE").fill_null(False).alias("promoter"),
            (pl.col("CpG.island").str.to_lowercase() == "yes").fill_null(False).alias("cpg_island"),
        )
        .filter(pl.col("chr").is_not_null())
    )


def read_orthologs(path: Path) -> pl.DataFrame:
    """``class_key, organism (human|mouse), symbol, entrez`` from MGI's homology report."""
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
    organism = pl.col("NCBI Taxon ID").replace_strict({"9606": "human", "10090": "mouse"}, default=None)
    return frame.select(
        pl.col("DB Class Key").alias("class_key"),
        organism.alias("organism"),
        pl.col("Symbol").alias("symbol"),
        pl.col("EntrezGene ID").cast(pl.Int64, strict=False).alias("entrez"),
    ).filter(pl.col("organism").is_not_null())


def mouse_gene_sets(human: pl.DataFrame, orthologs: pl.DataFrame) -> pl.DataFrame:
    """``gene_set, symbol, entrez`` of mouse genes orthologous to each human set member.

    A human gene is found in the report by Entrez ID, else by symbol; every mouse gene of its
    class is then a member. Human genes without a mouse ortholog drop out.
    """
    people = orthologs.filter(pl.col("organism") == "human")
    mice = orthologs.filter(pl.col("organism") == "mouse").select("class_key", "symbol", "entrez")
    by_entrez = human.filter(pl.col("entrez").is_not_null()).join(people.select("class_key", "entrez"), on="entrez")
    by_symbol = human.join(people.select("class_key", "symbol"), on="symbol")
    classes = pl.concat([by_entrez.select("gene_set", "class_key"), by_symbol.select("gene_set", "class_key")])
    return (
        classes.unique()
        .join(mice, on="class_key")
        .select("gene_set", "symbol", "entrez")
        .unique()
        .sort("gene_set", "symbol")
    )


def mouse_membership(probes: pl.DataFrame, gene_sets: pl.DataFrame, *, promoter_only: bool = True) -> pl.DataFrame:
    """``feature`` (probe), ``group`` (gene set) and ``gene``, matched by mouse Entrez ID."""
    base = probes.filter(pl.col("promoter")) if promoter_only else probes
    return (
        base.filter(pl.col("entrez").is_not_null())
        .join(gene_sets.select("gene_set", "entrez"), on="entrez")
        .select(pl.col("cpg").alias("feature"), pl.col("gene_set").alias("group"), pl.col("symbol").alias("gene"))
        .unique(subset=["feature", "group"], keep="first", maintain_order=True)
    )


def prepare_mouse_annotation(paths: ProjectPaths) -> MouseAnnotationLayout:
    lay = MouseAnnotationLayout(paths)
    human = AnnotationLayout(paths).gene_sets
    if not human.exists():
        raise FileNotFoundError(f"{human} is missing: run `entropy-age fetch annotation` first")
    lay.interim.mkdir(parents=True, exist_ok=True)
    read_platform_table(lay.input / PLATFORM_FILE).write_parquet(lay.probes)
    orthologs = read_orthologs(lay.input / ORTHOLOG_FILE)
    mouse_gene_sets(pl.read_parquet(human), orthologs).write_parquet(lay.gene_sets)
    return lay
