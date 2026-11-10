"""GTEx v10 RNA-seq counts from an external, already parsed copy (read-only).

entropy-age does not download or parse GTEx's 19,788 RNA-seq samples itself. It reads a folder in the
layout the rna-age project writes, named by the setting ``ENTROPY_AGE_GTEX_RNA``:

* ``samples.parquet`` - one row per sample: ``sample_id``, ``donor_id``, ``tissue``, ``sex``, the decade
  bin of age as ``age_min``/``age_max``, and ``exclude_reason`` (cell lines are excluded);
* ``genes.parquet`` - ``gene_id`` (Ensembl without version) in the row order of every count part;
* ``counts/part-NNNNN.parquet`` - gene x sample blocks of raw counts, with ``counts/index.json``
  mapping each sample to its part.

Technical covariates come from GTEx's public annotation files, in a folder named by
``ENTROPY_AGE_GTEX_ATTRIBUTES``: ``GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt`` (RNA
integrity ``SMRIN``, ischemic time ``SMTSISCH``) and ``GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt``
(Hardy death scale ``DTHHRDY``). Either setting missing, or a file missing under it, is a refusal
that names what is absent; nothing is substituted.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

COUNTS_ENV = "ENTROPY_AGE_GTEX_RNA"
ATTRIBUTES_ENV = "ENTROPY_AGE_GTEX_ATTRIBUTES"
SAMPLE_ATTRIBUTES = "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
SUBJECT_PHENOTYPES = "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"


class GtexRnaMissing(RuntimeError):
    """The external GTEx RNA folder or its annotation files are not where the settings say."""


@dataclass(frozen=True, slots=True)
class GtexRnaSource:
    counts: Path  # the parsed folder: samples.parquet, genes.parquet, counts/
    attributes: Path  # the folder with GTEx's two annotation tables

    @classmethod
    def from_environment(cls) -> GtexRnaSource:
        missing = [name for name in (COUNTS_ENV, ATTRIBUTES_ENV) if not os.environ.get(name)]
        if missing:
            raise GtexRnaMissing(
                f"{' and '.join(missing)} not set: point {COUNTS_ENV} at a parsed GTEx v10 RNA folder "
                f"(samples.parquet, genes.parquet, counts/) and {ATTRIBUTES_ENV} at the folder holding "
                f"{SAMPLE_ATTRIBUTES} and {SUBJECT_PHENOTYPES}"
            )
        source = cls(Path(os.environ[COUNTS_ENV]), Path(os.environ[ATTRIBUTES_ENV]))
        source.check()
        return source

    def needed(self) -> list[Path]:
        return [
            self.counts / "samples.parquet",
            self.counts / "genes.parquet",
            self.counts / "counts" / "index.json",
            self.attributes / SAMPLE_ATTRIBUTES,
            self.attributes / SUBJECT_PHENOTYPES,
        ]

    def check(self) -> None:
        absent = [str(p) for p in self.needed() if not p.exists()]
        if absent:
            raise GtexRnaMissing(f"GTEx RNA inputs missing: {', '.join(absent)}")


def read_samples(source: GtexRnaSource) -> pl.DataFrame:
    """Every non-excluded sample with its donor, tissue, age bin midpoint and technical covariates.

    ``age`` is the middle of the decade bin (20-29 -> 25), as in the methylation GTEx test.
    ``rin``, ``ischemic_minutes`` and ``hardy`` come from GTEx's annotation tables; a sample they do
    not list keeps nulls there, and the caller decides what to do with them.
    """
    samples = pl.read_parquet(source.counts / "samples.parquet").filter(pl.col("exclude_reason").is_null())
    attributes = pl.read_csv(
        source.attributes / SAMPLE_ATTRIBUTES,
        separator="\t",
        columns=["SAMPID", "SMRIN", "SMTSISCH"],
        schema_overrides={"SAMPID": pl.String, "SMRIN": pl.Float64, "SMTSISCH": pl.Float64},
        quote_char=None,
        null_values=["", "NA"],
    )
    subjects = pl.read_csv(
        source.attributes / SUBJECT_PHENOTYPES,
        separator="\t",
        columns=["SUBJID", "DTHHRDY"],
        schema_overrides={"SUBJID": pl.String, "DTHHRDY": pl.String},
        quote_char=None,
        null_values=["", "NA"],
    )
    return (
        samples.join(attributes, left_on="sample_id", right_on="SAMPID", how="left")
        .join(subjects, left_on="donor_id", right_on="SUBJID", how="left")
        .select(
            pl.col("sample_id").alias("sample"),
            pl.col("donor_id").alias("donor"),
            "tissue",
            "sex",
            ((pl.col("age_min") + pl.col("age_max") + 1) / 2).alias("age"),
            pl.col("SMRIN").alias("rin"),
            pl.col("SMTSISCH").alias("ischemic_minutes"),
            pl.col("DTHHRDY").alias("hardy"),
        )
        .sort("tissue", "donor", "sample")
    )


def read_genes(source: GtexRnaSource) -> list[str]:
    return pl.read_parquet(source.counts / "genes.parquet", columns=["gene_id"])["gene_id"].to_list()


def read_counts(source: GtexRnaSource, sample_ids: list[str]) -> np.ndarray:
    """Raw counts, samples x genes in ``read_genes`` order, for ``sample_ids`` in the order given."""
    index: dict[str, str] = json.loads((source.counts / "counts" / "index.json").read_text())
    unknown = [s for s in sample_ids if s not in index]
    if unknown:
        raise KeyError(f"{len(unknown)} samples not in the counts index, e.g. {unknown[:3]}")
    by_part: dict[str, list[str]] = {}
    for sample in sample_ids:
        by_part.setdefault(index[sample], []).append(sample)
    columns: dict[str, np.ndarray] = {}
    for part, ids in sorted(by_part.items()):
        frame = pl.read_parquet(source.counts / "counts" / part, columns=ids)
        columns.update({c: frame[c].to_numpy() for c in ids})
    return np.stack([columns[s] for s in sample_ids]) if sample_ids else np.zeros((0, len(read_genes(source))))
