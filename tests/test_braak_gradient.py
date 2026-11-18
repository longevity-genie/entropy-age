import polars as pl

from entropy_age.experiments.braak_gradient import ORDER, gradients


def test_gradient_is_the_within_donor_slope_on_region_order():
    rows = []
    for d, step in (("a", 1.0), ("b", -1.0), ("c", 0.3)):
        for region, r in ORDER.items():
            rows.append({"donor": d, "region": region, "components": 0, "B_x": step * r + 0.5 * (d == "c")})
    g = gradients(pl.DataFrame(rows), ["B_x"]).sort("donor")
    # opposite trends give opposite slopes
    s = g["B_x"].to_list()
    assert s[0] > 0 > s[1]
