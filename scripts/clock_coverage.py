# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = [
#   "biolearn @ git+https://github.com/bio-learn/biolearn@0d714f5a0c0ad5743b6315f3edab1f86112d394f",
#   "torch", "seaborn", "matplotlib", "requests", "pyyaml", "appdirs", "scikit-learn", "pyarrow",
# ]
# ///
"""Clock coverage on a fixed set of arrays, with the rule of scripts/clocks_biolearn.py and no prediction.

``uv run --script scripts/clock_coverage.py MANIFEST OUT``. MANIFEST is a CSV with ``series, track, parquet,
samples`` (``samples``: the array ids joined by ``;``), written by ``clock_preprocessing.write_manifest`` with the
arrays both tracks exported that pass the current QC, so the two tracks' coverage is judged on the same arrays. OUT
gets one row per series, track and clock: the arrays, the share of them with at least 90 % of the clock's absolute
weight measured, the median measured share, ``kept`` or ``coverage``, and what the result stands for: the parquet
read, its fingerprint (size and modification time, read before the betas) and a SHA-256 of the sorted array ids, which
clock_preprocessing checks against the current arrays and exports before it reuses the result. The weights and the
measured-row rule are
those of clocks_biolearn.py (imported from beside this file); clocks it does not score (other model types) are not
listed, their status stays the one the clock track's panel gives.
"""

import hashlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clocks_biolearn import HYBRID_ROW_SHARE, MIN_ARRAY_SHARE, MIN_WEIGHT_SHARE, grimage_weights, site_weights  # noqa: E402


def clock_weights() -> dict[str, tuple[pd.Series, bool]]:
    from biolearn.model import model_definitions
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()
    out = {}
    for name, entry in model_definitions.items():
        if entry.get("species") != "Human":
            continue
        grim = entry["model"]["type"] == "GrimageModel"
        if not grim and "Methylation" not in entry["model"]["type"]:
            continue
        try:
            model = gallery.get(name)
            weight = grimage_weights(getattr(model, "clock", model)) if grim else site_weights(model)
        except Exception:  # noqa: BLE001 - clocks_biolearn lists these as not built
            continue
        if weight.sum() > 0:
            out[name] = (weight, grim)
    return out


def samples_digest(samples) -> str:
    """SHA-256 of the sorted array ids joined by ``;`` (clock_preprocessing.samples_digest)."""
    return hashlib.sha256(";".join(sorted(str(s) for s in samples)).encode()).hexdigest()


def fingerprint(path: str) -> str:
    """Size and modification time of a betas export (clock_preprocessing.fingerprint)."""
    st = Path(path).stat()
    return f"{st.st_size}:{st.st_mtime_ns}"


def main(manifest: str, dst: str) -> None:
    jobs = pd.read_csv(manifest, dtype=str)
    weights = clock_weights()
    needed = sorted({c for w, _ in weights.values() for c in w.index})
    rows = []
    for job in jobs.itertuples():
        samples = [s for s in str(job.samples).split(";") if s]
        stamp = fingerprint(job.parquet)
        digest = samples_digest(samples)
        betas = pd.read_parquet(job.parquet, columns=["cpg", *samples], filters=[("cpg", "in", needed)])
        betas = betas.set_index("cpg")
        for name, (weight, grim) in weights.items():
            present = betas.reindex(weight.index).notna()
            kept_rows = present.mean(axis=1) >= (0.0 if grim else HYBRID_ROW_SHARE)
            measured = present.loc[kept_rows].mul(weight[kept_rows], axis=0).sum(axis=0) / weight.sum()
            share = float((measured >= MIN_WEIGHT_SHARE).mean()) if len(samples) else float("nan")
            rows.append(
                {
                    "series": job.series,
                    "track": job.track,
                    "clock": name,
                    "arrays": len(samples),
                    "arrays_covered": share,
                    "median_weight_share": float(measured.median()) if len(samples) else float("nan"),
                    "status": "kept" if share >= MIN_ARRAY_SHARE else "coverage",
                    "parquet": job.parquet,
                    "parquet_fingerprint": stamp,
                    "samples_sha256": digest,
                }
            )
        if fingerprint(job.parquet) != stamp:  # rewritten while read: the result stands for neither version
            for row in rows:
                if row["series"] == job.series and row["track"] == job.track:
                    row["parquet_fingerprint"] = "changed while read"
        print(f"{job.series} {job.track}: {len(samples)} arrays", flush=True)
    pd.DataFrame(rows).to_csv(dst, index=False)


if __name__ == "__main__":
    main(*sys.argv[1:3])
