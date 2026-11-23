"""CpGs whose methylation differs between TCGA lung adenocarcinomas and normal lung (open-access GDC level-3 betas).

Mean beta per CpG in solid-tissue normals and primary tumours; a CpG is tumour-changed when the means differ by at
least the cut (0.1 and 0.2 written). Used to leave the tumour field out of the tumour-adjacent lung series.
``uv run python scripts/tcga_tumour_cpgs.py`` (needs data/input/tcga_luad from logs/entropy_fetch_tcga_luad.py).
"""

from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "input" / "tcga_luad"
OUT = ROOT / "config" / "series" / "lists"


def mean_beta(kind: str) -> pl.DataFrame:
    frames = [
        pl.read_csv(f, separator="\t", has_header=False, new_columns=["cpg", f.stem], null_values=["NA"])
        for f in sorted((DATA / kind).glob("*.txt"))
    ]
    wide = frames[0]
    for f in frames[1:]:
        wide = wide.join(f, on="cpg", how="full", coalesce=True)
    values = wide.drop("cpg")
    print(f"{kind}: {values.width} arrays")
    return wide.select("cpg", pl.mean_horizontal(values.columns).alias(kind))


normal, tumour = mean_beta("normal"), mean_beta("tumour")
both = normal.join(tumour, on="cpg").drop_nulls(["normal", "tumour"]).with_columns(delta=pl.col("tumour") - pl.col("normal"))
OUT.mkdir(parents=True, exist_ok=True)
for cut in (0.1, 0.2):
    sel = both.filter(pl.col("delta").abs() >= cut).select("cpg", "normal", "tumour", "delta").sort("cpg")
    sel.select("cpg").write_csv(OUT / f"luad_tumour_cpgs_d{int(cut * 100)}.csv")  # ids only: keeps the list small
    print(f"|delta| >= {cut}: {sel.height} CpGs ({(sel['delta'] > 0).sum()} gain, {(sel['delta'] < 0).sum()} loss) of {both.height}")
