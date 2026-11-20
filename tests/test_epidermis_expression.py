import struct

import numpy as np
import polars as pl

from entropy_age.experiments.epidermis_expression import age_by_chip, read_tpm, read_women, to_symbols


def test_read_women_takes_one_row_per_woman_with_age_and_batch(tmp_path):
    path = tmp_path / "x.sdrf.txt"
    path.write_text(
        "Source Name\tCharacteristics[individual]\tCharacteristics[age]\tUnit[time unit]\tCharacteristics[batch]\t"
        "Array Data File\n"
        "a\ta\t36\tyear\t2\ta_Grn.idat\n"
        "a\ta\t36\tyear\t2\ta_Red.idat\n"
        "b\tb\t\tyear\t3\tb_Grn.idat\n"
        "c\tc\t70\tyear\t3\tc_Grn.idat\n"
    )
    women = read_women(path)
    assert women["sample"].to_list() == ["a", "c"]
    assert women["age"].to_list() == [36.0, 70.0]
    assert women["batch"].to_list() == ["2", "3"]


def test_read_tpm_logs_values_and_orders_columns_as_asked(tmp_path):
    path = tmp_path / "tpm.tsv"
    path.write_text("w1\tw2\tgene_id\n1\t3\tENSG01.4\n-5\t7\tENSG02\n")
    genes, expr = read_tpm(path, ["w2", "w1"])
    assert genes == ["ENSG01", "ENSG02"]
    assert np.allclose(expr, [[2.0, 1.0], [3.0, 0.0]])


def test_to_symbols_drops_unnamed_and_ambiguous_genes():
    hgnc = pl.DataFrame({"ensembl": ["E1", "E2", "E3", None], "symbol": ["A", "B", "B", "C"]})
    genes, expr = to_symbols(["E1", "E2", "E3", "E9"], np.arange(8.0).reshape(4, 2), hgnc)
    assert genes == ["A"]
    assert expr.tolist() == [[0.0, 1.0]]


def test_age_by_chip_is_one_when_chips_split_the_ages_and_zero_when_they_do_not():
    ages = np.array([20.0, 20.0, 60.0, 60.0])
    split = age_by_chip(ages, np.array(["a", "a", "b", "b"]))
    assert split["r2_age_by_chip"].to_list() == [1.0, 1.0]
    mixed = age_by_chip(ages, np.array(["a", "b", "a", "b"]))
    assert mixed["r2_age_by_chip"].to_list() == [0.0, 0.0]
    assert mixed["women"].to_list() == [2, 2]


def test_idat_header_reads_chip_and_position(tmp_path):
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location("idat_barcodes", Path("scripts/idat_barcodes.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def text(s: str) -> bytes:
        return bytes([len(s)]) + s.encode()

    body = text("202176300072") + text("BeadChip 8x5") + text("R01C01")
    start = 4 + 8 + 4 + 3 * 10
    offsets = [start, start + 13, start + 26]
    head = b"IDAT" + struct.pack("<q", 3) + struct.pack("<i", 3)
    head += b"".join(struct.pack("<Hq", code, off) for code, off in zip((402, 403, 404), offsets, strict=True))
    path = tmp_path / "w_Grn.idat"
    path.write_bytes(head + body)
    assert module.table(tmp_path) == [
        {"sample": "w", "chip": "202176300072", "chip_type": "BeadChip 8x5", "position": "R01C01"}
    ]
