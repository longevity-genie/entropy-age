"""Published epigenetic clocks' CpGs and weights, from biolearn's data at a pinned commit.

biolearn (bio-learn/biolearn) keeps one CSV per clock: ``CpGmarker,CoefficientTraining``, an ``intercept`` row
where the clock has one. Each file is checked against its git blob SHA-1 at ``BIOLEARN_COMMIT``, the same pin
as the just-biomarkers repository uses.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

BIOLEARN_COMMIT = "0d714f5a0c0ad5743b6315f3edab1f86112d394f"
URL = f"https://raw.githubusercontent.com/bio-learn/biolearn/{BIOLEARN_COMMIT}/biolearn/data/{{name}}"
#: clock -> (file, git blob SHA-1 at the pin, kind)
CLOCKS: dict[str, tuple[str, str, str]] = {
    "Horvath 2013": ("Horvath1.csv", "2c889c852087374b6bd0c6cdd19545beb226e029", "chronological"),
    "Horvath skin & blood": ("Horvath2.csv", "fea1667ce4481c4969f50318d72e4b8048b50993", "chronological"),
    "Hannum": ("Hannum.csv", "1a26380822523748d1092d8a9f1ccae4f4494e68", "chronological"),
    "Lin": ("Lin.csv", "f575affa18a0f54dca7bd2792f6387c9c124273b", "chronological"),
    "Zhang 2019": ("Zhang2019.csv", "3f68b5fd8aead7a8308e6a235057e7b6ab77cfaa", "chronological"),
    "Zhang 10": ("Zhang_10.csv", "7ac2526d41e609500fcdb3b85eb4ed45bd8e01f2", "mortality"),
    "PhenoAge": ("PhenoAge.csv", "2362abc9c736996056a95ff5cb3f02f08276ec1f", "phenotypic"),
    "DunedinPACE": ("DunedinPACE.csv", "b65bd0933816d7f074b35907216c235630ce028b", "pace"),
    "Cortical": ("DNAmClockCortical.csv", "fb63f3da68eaa47355f117ee4ebc6c86fe29fdb6", "chronological"),
    "DNAm telomere length": ("DNAmTL.csv", "c8c45081b7c5bdf9f9409fc72216d6e1e258fcfe", "phenotypic"),
    "epiTOC1": ("EpiTOC1.csv", "dbf53392a7180a2bbd1eaeefc7f00338189853bf", "mitotic"),
    "StocH": ("StocH.csv", "e025e52cd84ade2113fa617f7ad4a54390e403a5", "stochastic"),
    "StocP": ("StocP.csv", "3c3f565bbc27675dea83ddb568b9d06e8ee60e83", "stochastic"),
    "StocZ": ("StocZ.csv", "1e37c7940776ce1e2b6131eeb524e3b4b883eb9a", "stochastic"),
    "CausAge": ("YingCausAge.csv", "27fcc77f54110d51c129edea578710160355f6b8", "causal"),
    "DamAge": ("YingDamAge.csv", "42b0fa0086381212ff6651a7bbee1f50ea1164c4", "causal"),
    "AdaptAge": ("YingAdaptAge.csv", "87652423f0c8faa6e1f5450d9b98a2c2ab36375e", "causal"),
}


@dataclass(frozen=True, slots=True)
class ClockLayout:
    paths: ProjectPaths

    def file(self, clock: str) -> Path:
        return self.paths.input / "clocks" / CLOCKS[clock][0]


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def fetch_clocks(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    """Download every clock in ``CLOCKS``; a file whose hash is not the pinned one is removed and refused."""
    out = []
    for clock, (name, sha, _) in CLOCKS.items():
        path = download(URL.format(name=name), ClockLayout(paths).file(clock), overwrite=overwrite)
        got = git_blob_sha1(path.read_bytes())
        if got != sha:
            path.unlink()
            raise ValueError(f"{name}: blob {got} is not the pinned {sha} (biolearn {BIOLEARN_COMMIT[:12]})")
        out.append(path)
    return out


def read_clock(path: Path) -> pl.DataFrame:
    """``cpg`` and ``weight`` of a clock's CpGs; the intercept row and any non-CpG row are dropped."""
    text = path.read_bytes().decode("utf-8-sig")
    frame = pl.read_csv(text.encode(), schema_overrides={"CpGmarker": pl.String, "CoefficientTraining": pl.Float64})
    cpgs = frame.filter(pl.col("CpGmarker").str.contains(r"^c[gh]\d+$"))
    if cpgs.height == 0:
        raise ValueError(f"{path.name}: no CpG rows")
    return cpgs.select(pl.col("CpGmarker").alias("cpg"), pl.col("CoefficientTraining").alias("weight"))
