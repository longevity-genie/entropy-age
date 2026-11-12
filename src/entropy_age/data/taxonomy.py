"""Species lineages from NCBI Taxonomy, and the shared-ancestry covariance they give (roadmap RM10b).

A species' lineage is the list of its NCBI Taxonomy ancestors, root first. Two species share the leading part of
their lineages down to their most recent common ancestor; the length of that shared part, over the longest
lineage, is a covariance in the spirit of Brownian motion on a tree whose branch lengths are all one. It is a
taxonomy, not a dated phylogeny, so it orders relatedness without timing it; and because every branch counts one,
two sibling species share at most half their path below their deepest shared node, so a strong clade effect is
captured only in part (the correction is conservative toward independence).
"""

from __future__ import annotations

import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.download import USER_AGENT

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAUSE = 0.4  # NCBI asks for at most three requests a second without a key


def _get(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def binomial(name: str) -> str:
    """Genus and species of a name: subspecies and trailing words dropped (``Ceratotherium simum simum``)."""
    return " ".join(name.split()[:2])


def lookup(name: str) -> tuple[str, list[str]]:
    """The taxon id of ``name`` and its lineage of ancestor ids, root first, ending with the taxon itself."""
    term = urllib.parse.quote(f"{binomial(name)}[Scientific Name]")
    ids = ET.fromstring(_get(f"{EUTILS}/esearch.fcgi?db=taxonomy&term={term}")).findall("./IdList/Id")
    if len(ids) != 1:
        raise ValueError(f"taxonomy: {name!r} matched {len(ids)} taxa")
    taxid = ids[0].text
    time.sleep(PAUSE)
    taxon = ET.fromstring(_get(f"{EUTILS}/efetch.fcgi?db=taxonomy&id={taxid}")).find("./Taxon")
    lineage = [t.findtext("TaxId") for t in taxon.findall("./LineageEx/Taxon")]
    return taxid, [*lineage, taxid]


def fetch_lineages(names: list[str], cache: Path) -> pl.DataFrame:
    """``organism``, ``taxid`` and ``lineage`` (ids joined by ``;``) for every name, cached in ``cache`` (CSV)."""
    known = (
        pl.read_csv(cache)
        if cache.exists()
        else pl.DataFrame(schema={"organism": pl.String, "taxid": pl.String, "lineage": pl.String})
    )
    have = set(known["organism"].to_list())
    rows = []
    for name in sorted(set(names) - have):
        taxid, lineage = lookup(name)
        rows.append({"organism": name, "taxid": taxid, "lineage": ";".join(lineage)})
        time.sleep(PAUSE)
    out = pl.concat([known, pl.DataFrame(rows, schema=known.schema)]) if rows else known
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.sort("organism").write_csv(cache)
    return out.filter(pl.col("organism").is_in(names))


def shared_depth(a: list[str], b: list[str]) -> int:
    """How many leading ancestors two lineages share."""
    n = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            break
        n += 1
    return n


def lineage_covariance(lineages: list[list[str]]) -> np.ndarray:
    """Shared depth over the longest lineage, so a taxon has 1 with itself when its lineage is the longest.

    Unequal lineage lengths (NCBI ranks some groups more finely) leave the diagonal below one for shallow
    lineages; the diagonal is set to one so every taxon has the same variance.
    """
    n = len(lineages)
    longest = max(len(x) for x in lineages)
    c = np.array([[shared_depth(lineages[i], lineages[j]) for j in range(n)] for i in range(n)], dtype=float)
    c /= longest
    np.fill_diagonal(c, 1.0)
    return c
