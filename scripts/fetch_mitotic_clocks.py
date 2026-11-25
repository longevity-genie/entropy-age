# /// script
# dependencies = ["pyreadr", "pandas"]
# ///
"""Fetch the coefficients of three mitotic clocks from Teschendorff's EpiMitClocks R package
(github.com/aet21/EpiMitClocks, data/) and write them as CSV for ``experiments/mitotic_clocks.py`` (RESEARCH_PLAN
idea 3: is the shield a mitotic clock?).

``uv run --script scripts/fetch_mitotic_clocks.py data/input/mitotic`` writes:

- ``epitoc2.csv``: cpg, delta, beta0 (epiTOC2, Teschendorff 2020 Genome Med; ``estETOC3.m``)
- ``replitali.csv``: term, coefficient (RepliTali, Endicott et al. 2022 Nat Commun; first row the intercept)
- ``epicmit.csv``: cpg, class (epiCMIT hyper and hypo CpGs, Duran-Ferrer et al. 2020 Nat Cancer)

The repository states no licence; the coefficients are those published with the papers.
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import pyreadr

API = "https://api.github.com/repos/aet21/EpiMitClocks/contents/data"
FILES = ("estETOC3.rda", "Replitali.rda", "EpiCMITcpgs.rda")


def main(out: str) -> None:
    dst = Path(out)
    dst.mkdir(parents=True, exist_ok=True)
    listing = {e["name"]: e["download_url"] for e in json.load(urllib.request.urlopen(API, timeout=120))}
    for name in FILES:
        (dst / name).write_bytes(urllib.request.urlopen(listing[name], timeout=300).read())
    est = pyreadr.read_r(dst / "estETOC3.rda")["estETOC3.m"]
    est.rename_axis("cpg").reset_index().to_csv(dst / "epitoc2.csv", index=False)
    r = pyreadr.read_r(dst / "Replitali.rda")
    terms = r["replitali.cpg.v"].iloc[:, 0].tolist()
    coe = r["replitali.coe"].iloc[:, 0].tolist()
    with open(dst / "replitali.csv", "w") as f:
        f.write("term,coefficient\n")
        f.writelines(f"{t},{c!r}\n" for t, c in zip(terms, coe, strict=True))
    cm = pyreadr.read_r(dst / "EpiCMITcpgs.rda")["epiCMIT.df"]
    cm = cm.rename(columns={"Name": "cpg", "epiCMIT.class": "class"})[["cpg", "class"]]
    cm["class"] = cm["class"].str.replace("epiCMIT.", "", regex=False)
    cm.to_csv(dst / "epicmit.csv", index=False)
    for name in FILES:
        (dst / name).unlink()
    print(f"epiTOC2 {len(est)} CpGs, RepliTali {len(terms) - 1} CpGs, epiCMIT {len(cm)} CpGs -> {dst}")


if __name__ == "__main__":
    main(sys.argv[1])
