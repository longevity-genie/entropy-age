"""Probe-to-gene annotation for the Illumina 27k array and curated gene sets.

* Illumina HumanMethylation27 manifest (GEO GPL8490 supplementary file): gene,
  Entrez ID, RefSeq accession, distance to TSS and CpG-island flag per probe.
* Hart et al. core-essential (CEGv2, 684 genes) and non-essential (NEGv1, 927
  genes) reference sets from CRISPR screens, pinned to a commit of hart-lab/bagel (MIT).
* Eisenberg & Levanon (2013) human housekeeping genes (3,804 genes).
* HGNC's complete gene set, a pinned monthly release: approved, previous and alias symbols with
  Entrez and Ensembl IDs, used to put the gene sets on the Ensembl IDs of RNA-seq counts
  (``gene_set_ensembl``) and to tell protein-coding genes from the rest.

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
MANIFEST_450K = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
MANIFEST_450K_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/{MANIFEST_450K}"
#: UCSC RefGene regions counted as a promoter (the 27k array's own design targets these)
PROMOTER_REGIONS = ("TSS1500", "TSS200", "5'UTR", "1stExon")
HGNC_FILE = "hgnc_complete_set_2026-10-06.txt"
HGNC_URL = f"https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/{HGNC_FILE}"

SOURCES: dict[str, str] = {
    "GPL8490_HumanMethylation27_270596_v.1.2.csv.gz": MANIFEST_URL,
    "CEGv2.txt": HART_URL.format(name="CEGv2.txt"),
    "NEGv1.txt": HART_URL.format(name="NEGv1.txt"),
    "HK_genes.txt": HOUSEKEEPING_URL,
    MANIFEST_450K: MANIFEST_450K_URL,
    HGNC_FILE: HGNC_URL,
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

    @property
    def probes_450k(self) -> Path:
        return self.interim / "hm450_probes.parquet"

    @property
    def links_450k(self) -> Path:
        return self.interim / "hm450_gene_links.parquet"

    @property
    def hgnc(self) -> Path:
        return self.interim / "hgnc_genes.parquet"


def fetch_annotation(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    layout = AnnotationLayout(paths)
    return [download(url, layout.input / name, overwrite=overwrite) for name, url in SOURCES.items()]


def _assay_block(path: Path) -> pl.DataFrame:
    """The ``[Assay]`` block of an Illumina manifest, every column as text."""
    with gzip.open(path, "rt", encoding="latin-1") if str(path).endswith(".gz") else open(path) as handle:
        lines = handle.read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("[Assay]")) + 1
    end = next((i for i, line in enumerate(lines) if line.startswith("[Controls]")), len(lines))
    return pl.read_csv(io.StringIO("\n".join(lines[start:end])), infer_schema_length=0, truncate_ragged_lines=True)


def read_manifest(path: Path) -> pl.DataFrame:
    """The 27k manifest as one row per probe."""
    table = _assay_block(path)
    return table.select(
        pl.col("Name").alias("cpg"),
        pl.col("Chr").alias("chr"),
        pl.col("Symbol").str.strip_chars().alias("symbol"),
        pl.col("Gene_ID").str.strip_prefix("GeneID:").cast(pl.Int64, strict=False).alias("entrez"),
        pl.col("Accession").str.replace(r"\.\d+$", "").alias("refseq"),
        pl.col("Distance_to_TSS").cast(pl.Int64, strict=False).alias("distance_to_tss"),
        (pl.col("CPG_ISLAND").str.to_uppercase() == "TRUE").alias("cpg_island"),
    )


def solo_wcgw(sequence: str, flank: int = 35) -> bool:
    """Whether the probed CpG is a solo-WCGW: A or T on both sides, and no other CpG within ``flank`` bp.

    ``sequence`` is the manifest's ``Forward_Sequence``, the probed CpG written ``[CG]`` inside its genomic context.
    Solo-WCGW CpGs in late-replicating, CpG-poor domains lose methylation with cell divisions (Zhou et al. 2018).
    A sequence without the bracket, or too short a context, is not counted as one.
    """
    i = sequence.find("[CG]")
    if i < flank + 1 or len(sequence) < i + 4 + flank + 1:
        return False
    left, right = sequence[:i].upper(), sequence[i + 4 :].upper()
    if left[-1] not in "AT" or right[0] not in "AT":
        return False
    return "CG" not in left[-flank:] and "CG" not in right[:flank]


def read_manifest_450k(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The 450k manifest as ``(probes, links)``.

    ``probes`` has one row per CpG: chromosome, position (GRCh37, the manifest's ``MAPINFO``), CpG-island
    relation and whether the CpG sits in a
    promoter region of any gene. ``links`` has one row per CpG and gene: a 450k probe often lies in
    several genes, and its gene, RefSeq accession and region lists are parallel ``;``-separated
    fields. The manifest carries no Entrez IDs, so gene sets are matched by RefSeq, then symbol.
    """
    table = _assay_block(path).select(
        pl.col("Name").alias("cpg"),
        pl.col("CHR").alias("chr"),
        pl.col("MAPINFO").cast(pl.Int64, strict=False).alias("pos"),
        pl.col("Forward_Sequence").fill_null("").alias("sequence"),
        pl.col("Relation_to_UCSC_CpG_Island").fill_null("OpenSea").alias("island_relation"),
        pl.col("UCSC_RefGene_Name").fill_null("").str.split(";").alias("symbol"),
        pl.col("UCSC_RefGene_Accession").fill_null("").str.split(";").alias("refseq"),
        pl.col("UCSC_RefGene_Group").fill_null("").str.split(";").alias("region"),
    )
    links = (
        table.filter(pl.col("symbol").list.len() == pl.col("region").list.len())
        .select("cpg", "symbol", "refseq", "region")
        .explode("symbol", "refseq", "region")
        .filter(pl.col("symbol") != "")
        .with_columns(promoter=pl.col("region").is_in(list(PROMOTER_REGIONS)))
        .unique(maintain_order=True)
    )
    promoter_cpgs = links.filter(pl.col("promoter"))["cpg"].unique()
    probes = table.select(
        "cpg",
        "chr",
        "pos",
        "island_relation",
        (pl.col("island_relation") == "Island").alias("cpg_island"),
        pl.col("sequence").map_elements(solo_wcgw, return_dtype=pl.Boolean).alias("solo_wcgw"),
        pl.col("cpg").is_in(promoter_cpgs.implode()).alias("promoter"),
    )
    return probes, links


def link_membership(links: pl.DataFrame, gene_sets: pl.DataFrame, *, promoter_only: bool = True) -> pl.DataFrame:
    """``feature`` (CpG), ``group`` and ``gene`` from a CpG-to-gene link table: RefSeq first, then symbol."""
    base = links.filter(pl.col("promoter")) if promoter_only else links
    matched = []
    for key in ("refseq", "symbol"):
        genes = gene_sets.filter(pl.col(key).is_not_null()).select("gene_set", pl.col(key).alias("key"))
        hits = base.select("cpg", "symbol", pl.col(key).alias("key")).join(genes, on="key")
        matched.append(
            hits.select(
                pl.col("cpg").alias("feature"), pl.col("gene_set").alias("group"), pl.col("symbol").alias("gene")
            )
        )
    return pl.concat(matched).unique(subset=["feature", "group"], keep="first", maintain_order=True)


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
    probes, links = read_manifest_450k(layout.input / MANIFEST_450K)
    probes.write_parquet(layout.probes_450k)
    links.write_parquet(layout.links_450k)
    read_hgnc(layout.input / HGNC_FILE).write_parquet(layout.hgnc)
    return layout


def read_hgnc(path: Path) -> pl.DataFrame:
    """``symbol, previous, aliases, entrez, ensembl, protein_coding`` per HGNC gene (lists split on ``|``)."""
    raw = pl.read_csv(path, separator="\t", infer_schema_length=0)
    return raw.select(
        "symbol",
        pl.col("prev_symbol").str.split("|").fill_null([]).alias("previous"),
        pl.col("alias_symbol").str.split("|").fill_null([]).alias("aliases"),
        pl.col("entrez_id").cast(pl.Int64, strict=False).alias("entrez"),
        pl.col("ensembl_gene_id").alias("ensembl"),
        (pl.col("locus_group") == "protein-coding gene").alias("protein_coding"),
    )


def gene_set_ensembl(gene_sets: pl.DataFrame, hgnc: pl.DataFrame) -> pl.DataFrame:
    """``feature`` (Ensembl gene), ``group``, ``gene`` (approved symbol) and ``matched_by`` for every set gene.

    A gene is matched by its Entrez ID first, then by its symbol as HGNC approves it today, then as a
    previous symbol, then as an alias. A previous symbol or an alias counts only when it names exactly
    one HGNC gene: the old lists use pre-2019 names (``AARS`` is now ``AARS1``), and an ambiguous
    alias would put a gene into a set it was never in. Genes without an Ensembl ID are left out.
    """
    ids = hgnc.filter(pl.col("ensembl").is_not_null())
    sets = gene_sets.with_row_index("row")

    def unique_key(column: str) -> pl.DataFrame:
        long = ids.select("symbol", "ensembl", pl.col(column).alias("key")).explode("key").drop_nulls("key")
        counts = long.group_by("key").agg(pl.len().alias("k"))
        return long.join(counts.filter(pl.col("k") == 1), on="key").select("key", "symbol", "ensembl")

    lookups = [
        ("entrez", "entrez", ids.select(pl.col("entrez").alias("key"), "symbol", "ensembl").drop_nulls("key")),
        ("symbol", "symbol", ids.select(pl.col("symbol").alias("key"), "symbol", "ensembl")),
        ("previous", "symbol", unique_key("previous")),
        ("alias", "symbol", unique_key("aliases")),
    ]
    found, rest = [], sets
    for how, column, table in lookups:
        hits = (
            rest.filter(pl.col(column).is_not_null())
            .join(table, left_on=column, right_on="key", how="inner", suffix="_hgnc")
            .unique(subset="row", keep="first", maintain_order=True)
        )
        found.append(
            hits.select(
                "row",
                pl.col("ensembl").alias("feature"),
                pl.col("gene_set").alias("group"),
                pl.col("symbol_hgnc").alias("gene"),
                pl.lit(how).alias("matched_by"),
            )
        )
        rest = rest.filter(~pl.col("row").is_in(hits["row"].implode()))
    out = pl.concat(found).sort("row").drop("row")
    return out.unique(subset=["feature", "group"], keep="first", maintain_order=True)


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


def symbol_ensembl(symbols: list[str], hgnc: pl.DataFrame) -> pl.DataFrame:
    """``symbol`` (as given, e.g. a mouse or rat symbol), ``ensembl`` (human) and ``matched_by``.

    The upper-case symbol is looked up as an HGNC approved symbol, then as a previous symbol, then as
    an alias; a previous symbol or alias counts only when it names one HGNC gene. Two given symbols that
    reach the same human gene are both dropped, since neither can stand for it alone. Symbols that
    match nothing are left out.
    """
    ids = hgnc.filter(pl.col("ensembl").is_not_null())

    def unique_key(column: str) -> pl.DataFrame:
        long = ids.select("ensembl", pl.col(column).alias("key")).explode("key").drop_nulls("key")
        counts = long.group_by("key").agg(pl.len().alias("k"))
        return long.join(counts.filter(pl.col("k") == 1), on="key").select("key", "ensembl")

    lookups = [
        ("symbol", ids.select(pl.col("symbol").alias("key"), "ensembl")),
        ("previous", unique_key("previous")),
        ("alias", unique_key("aliases")),
    ]
    rest = (
        pl.DataFrame({"symbol": symbols})
        .unique(maintain_order=True)
        .with_columns(key=pl.col("symbol").str.to_uppercase())
    )
    found = []
    for how, table in lookups:
        hits = rest.join(table, on="key", how="inner").unique(subset="symbol", keep="first", maintain_order=True)
        found.append(hits.select("symbol", "ensembl", pl.lit(how).alias("matched_by")))
        rest = rest.filter(~pl.col("symbol").is_in(hits["symbol"].implode()))
    out = pl.concat(found)
    shared = out.group_by("ensembl").agg(pl.len().alias("k")).filter(pl.col("k") > 1)
    return out.filter(~pl.col("ensembl").is_in(shared["ensembl"].implode()))
