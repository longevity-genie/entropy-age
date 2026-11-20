"""Sentrix chip and position of every IDAT in a folder, read from the file header.

Some depositors rename IDAT files (``1620_9005_Grn.idat``), so the chip is not in the name; the binary header still
carries it (fields 402 barcode, 404 position, 403 chip type). Writes ``sample, chip, position, chip_type`` as CSV.

    uv run python scripts/idat_barcodes.py IDAT_DIR OUT.csv
"""

from __future__ import annotations

import csv
import struct
import sys
from pathlib import Path

FIELDS = {402: "chip", 403: "chip_type", 404: "position"}


def _string(handle) -> str:
    """A .NET length-prefixed string: 7-bit encoded length, then the bytes."""
    length, shift = 0, 0
    while True:
        byte = handle.read(1)[0]
        length |= (byte & 0x7F) << shift
        if not byte & 0x80:
            break
        shift += 7
    return handle.read(length).decode("ascii", errors="replace")


def header(path: Path) -> dict[str, str]:
    """``chip``, ``position`` and ``chip_type`` of one IDAT (version 3) file."""
    with path.open("rb") as handle:
        if handle.read(4) != b"IDAT":
            raise ValueError(f"{path.name}: not an IDAT file")
        (version,) = struct.unpack("<q", handle.read(8))
        if version != 3:
            raise ValueError(f"{path.name}: IDAT version {version}, only 3 is read")
        (n_fields,) = struct.unpack("<i", handle.read(4))
        offsets = dict(struct.unpack("<Hq", handle.read(10)) for _ in range(n_fields))
        out = {}
        for code, name in FIELDS.items():
            if code in offsets:
                handle.seek(offsets[code])
                out[name] = _string(handle)
    return out


def table(folder: Path) -> list[dict[str, str]]:
    rows = []
    for path in sorted(folder.glob("*_Grn.idat")):
        rows.append({"sample": path.name.removesuffix("_Grn.idat"), **header(path)})
    return rows


if __name__ == "__main__":
    rows = table(Path(sys.argv[1]))
    with open(sys.argv[2], "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["sample", "chip", "position", "chip_type"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} arrays on {len({r.get('chip') for r in rows})} chips", file=sys.stderr)
