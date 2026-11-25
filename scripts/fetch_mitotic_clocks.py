# /// script
# dependencies = ["rdata", "pandas"]
# ///
"""Fetch the coefficients of the mitotic clocks in Teschendorff's EpiMitClocks R package
(github.com/aet21/EpiMitClocks, data/, pinned to commit ``COMMIT``) and write them as CSV for
``experiments/mitotic_clocks.py`` (RESEARCH_PLAN idea 3: is the shield a mitotic clock?).

``uv run --script scripts/fetch_mitotic_clocks.py data/input/mitotic`` writes:

- ``epitoc2.csv``: cpg, delta, beta0; epiTOC2 (Teschendorff 2020 Genome Med), 163 CpGs, ``dataETOC3.l$epiTOC2``,
  the table the package's ``epiTOC2()`` uses
- ``epitoc3.csv``: the same for the package's later epiTOC3 (170 CpGs, ``dataETOC3.l$epiTOC3``)
- ``epitoc1_cpgs.csv``: cpg; EpiTOC1's 385 CpGs (``dataETOC3.l$epiTOC``), used only to keep clock CpGs out of the
  drift measures
- ``hypoclock.csv``: cpg; the 678 solo-WCGW CpGs in partially methylated domains (Zhou et al. 2018 Nat Genet),
  ``dataETOC3.l$cm-solo-WCGW``
- ``replitali.csv``: term, coefficient (RepliTali, Endicott et al. 2022 Nat Commun; first row the intercept)
- ``epicmit.csv``: cpg, class (epiCMIT hyper and hypo CpGs, Duran-Ferrer et al. 2020 Nat Cancer)

The repository states no licence; the coefficients are those published with the papers.
"""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import pandas as pd
import rdata

COMMIT = "2c236cd651a111b018eab353bcdf7137bdcc1449"
RAW = f"https://raw.githubusercontent.com/aet21/EpiMitClocks/{COMMIT}/data"
FILES = ("dataETOC3.rda", "Replitali.rda", "EpiCMITcpgs.rda")


def _table(arr) -> pd.DataFrame:
    """An R matrix read by rdata (xarray, rows named by CpG) as cpg, <columns>."""
    frame = pd.DataFrame(arr.values, index=arr.coords["dim_0"].values, columns=arr.coords["dim_1"].values)
    return frame.rename_axis("cpg").reset_index()


def main(out: str) -> None:
    dst = Path(out)
    dst.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        (dst / name).write_bytes(urllib.request.urlopen(f"{RAW}/{name}", timeout=300).read())
    lists = rdata.read_rda(dst / "dataETOC3.rda")["dataETOC3.l"]
    epitoc2, epitoc3 = _table(lists["epiTOC2"]), _table(lists["epiTOC3"])
    epitoc2.to_csv(dst / "epitoc2.csv", index=False)
    epitoc3.to_csv(dst / "epitoc3.csv", index=False)
    pd.DataFrame({"cpg": list(lists["epiTOC"])}).to_csv(dst / "epitoc1_cpgs.csv", index=False)
    pd.DataFrame({"cpg": list(lists["cm-solo-WCGW"])}).to_csv(dst / "hypoclock.csv", index=False)
    r = rdata.read_rda(dst / "Replitali.rda")
    terms, coe = list(r["replitali.cpg.v"]), [float(v) for v in r["replitali.coe"]]
    pd.DataFrame({"term": terms, "coefficient": coe}).to_csv(dst / "replitali.csv", index=False)
    cm = rdata.read_rda(dst / "EpiCMITcpgs.rda")["epiCMIT.df"]
    cm = cm.rename(columns={"Name": "cpg", "epiCMIT.class": "class"})[["cpg", "class"]]
    cm["class"] = cm["class"].astype(str).str.replace("epiCMIT.", "", regex=False)
    cm.to_csv(dst / "epicmit.csv", index=False)
    for name in FILES:
        (dst / name).unlink()
    (dst / "SOURCE").write_text(f"github.com/aet21/EpiMitClocks@{COMMIT}\n")
    print(
        f"epiTOC2 {len(epitoc2)}, epiTOC3 {len(epitoc3)}, EpiTOC1 {len(lists['epiTOC'])}, "
        f"HypoClock {len(lists['cm-solo-WCGW'])}, RepliTali {len(terms) - 1}, epiCMIT {len(cm)} CpGs -> {dst}"
    )


if __name__ == "__main__":
    main(sys.argv[1])
