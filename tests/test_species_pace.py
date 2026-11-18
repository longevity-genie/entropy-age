import polars as pl

from entropy_age.experiments.species_pace import pair_agreement, within_species


def _paces():
    rows = []
    for i in range(6):
        for t, d, n in (("Blood", 3.0, 1.0), ("Liver", 2.0, 2.0), ("Skin", 1.0, 3.0)):
            rows.append({"species": f"s{i}", "tissue": t, "directed_year": d + 0.01 * i, "noise_year": n})
    return pl.DataFrame(rows)


def test_within_species_finds_opposite_orders():
    w = within_species(_paces(), "year")
    assert w.height == 6 and (w["rho"] == -1.0).all()


def test_pairs_count_species_ordered_alike():
    p = pair_agreement(_paces(), "year").filter((pl.col("tissue_a") == "Blood") & (pl.col("tissue_b") == "Liver"))
    r = p.row(0, named=True)
    assert r["species"] == 6 and r["directed_share_a_faster"] == 1.0 and r["noise_share_a_faster"] == 0.0
