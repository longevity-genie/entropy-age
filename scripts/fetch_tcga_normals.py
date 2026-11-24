"""Fetch TCGA normal-adjacent 450k arrays (open IDATs) with age and sex, for the field-cancerisation test
(RESEARCH_PLAN, 2026-11-24: does the shield break down in tissue next to a tumour?).

``uv run python scripts/fetch_tcga_normals.py OUT_DIR TCGA-KIRC [TCGA-BRCA ...]`` writes, per project,
``OUT_DIR/<project>/idat/<file name>`` (the GDC's open "Masked Intensities" IDATs of "Solid Tissue Normal" samples on
the 450k) and ``OUT_DIR/<project>/samples.csv``: one row per array pair with the sample barcode, the Sentrix basename
the IDAT files share, age at diagnosis in years, sex and the case id. Downloads go through the GDC data endpoint in
batches (tar.gz), skip files already present, and check each file's md5 against the GDC's.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

API = "https://api.gdc.cancer.gov"
FIELDS = ",".join(
    [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.demographic.sex_at_birth",
        "cases.diagnoses.age_at_diagnosis",
        "cases.demographic.days_to_birth",
    ]
)
BATCH = 40


def post(endpoint: str, payload: dict) -> bytes:
    req = urllib.request.Request(
        f"{API}/{endpoint}", data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.read()


def files(project: str) -> list[dict]:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Solid Tissue Normal"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Intensities"]}},
            {"op": "in", "content": {"field": "platform", "value": ["Illumina Human Methylation 450"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    out = json.loads(post("files", {"filters": filters, "fields": FIELDS, "size": 5000}))
    return out["data"]["hits"]


def age_years(case: dict) -> float | None:
    for d in case.get("diagnoses") or []:
        if d.get("age_at_diagnosis") is not None:
            return round(d["age_at_diagnosis"] / 365.25, 2)
    dtb = (case.get("demographic") or {}).get("days_to_birth")
    return round(-dtb / 365.25, 2) if dtb is not None else None


def main(out_dir: str, projects: list[str]) -> None:
    for project in projects:
        hits = files(project)
        base = Path(out_dir) / project
        idat = base / "idat"
        idat.mkdir(parents=True, exist_ok=True)
        rows = {}
        for h in hits:
            case = h["cases"][0]
            sample = next(s["submitter_id"] for s in case["samples"] if s["sample_type"] == "Solid Tissue Normal")
            basename = h["file_name"].rsplit("_", 1)[0]
            rows[basename] = {
                "sample": sample,
                "basename": basename,
                "case_id": case["case_id"],
                "age": age_years(case),
                "sex": (case.get("demographic") or {}).get("sex_at_birth"),
            }
        todo = [h for h in hits if not (idat / h["file_name"]).exists()]
        print(f"{project}: {len(hits)} IDATs ({len(rows)} arrays), {len(todo)} to fetch", flush=True)
        for i in range(0, len(todo), BATCH):
            chunk = todo[i : i + BATCH]
            blob = post("data", {"ids": [h["file_id"] for h in chunk]})
            md5 = {h["file_name"]: h["md5sum"] for h in chunk}
            if len(chunk) == 1:
                members = [(chunk[0]["file_name"], blob)]
            else:
                with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
                    members = [(Path(m.name).name, tar.extractfile(m).read()) for m in tar.getmembers() if m.isfile()]
            for name, data in members:
                if name not in md5:
                    continue  # the archive's MANIFEST
                if hashlib.md5(data).hexdigest() != md5[name]:
                    raise ValueError(f"{project}: md5 mismatch for {name}")
                (idat / name).write_bytes(data)
            print(f"{project}: {min(i + BATCH, len(todo))}/{len(todo)}", flush=True)
        with open(base / "samples.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["sample", "basename", "case_id", "age", "sex"])
            w.writeheader()
            w.writerows(sorted(rows.values(), key=lambda r: r["sample"]))
        ages = [r["age"] for r in rows.values() if r["age"] is not None]
        print(f"{project}: samples.csv written, {len(rows)} arrays, ages {min(ages):.0f}-{max(ages):.0f}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
