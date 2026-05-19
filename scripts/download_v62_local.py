"""
Download AADR v62.0 1240K public files from Harvard Dataverse to data/input_data/.

Used when the local AADR copy was deleted and v66 (on R2) is in TGENO format
that the existing readers don't support.

Streams directly to disk with progress reporting.
"""
from __future__ import annotations

import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "data" / "input_data"
DEST.mkdir(parents=True, exist_ok=True)

DATAVERSE = "https://dataverse.harvard.edu/api/access/datafile"

# File IDs from Dataverse dataset version 9.1 (AADR v62.0)
FILES = [
    (10537414, "v62.0_1240k_public.ind",   "0.00 GB"),
    (10537413, "v62.0_1240k_public.anno",  "0.01 GB"),
    (10537415, "v62.0_1240k_public.snp",   "0.08 GB"),
    (10537126, "v62.0_1240k_public.geno",  "5.44 GB"),
]


def download(file_id: int, filename: str, size_str: str) -> None:
    dest = DEST / filename
    if dest.exists():
        size_mb = dest.stat().st_size / 1e6
        if size_mb > 10:  # roughly correct size
            print(f"  [SKIP]  {filename} already exists ({size_mb:.1f} MB)")
            return

    url = f"{DATAVERSE}/{file_id}"
    print(f"  downloading {filename}  ({size_str})  → {dest}")

    with requests.get(url, stream=True, timeout=300, allow_redirects=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        written = 0
        last_pct = -1
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                fh.write(chunk)
                written += len(chunk)
                if total:
                    pct = int(written / total * 100)
                    if pct // 10 != last_pct // 10:
                        last_pct = pct
                        print(f"    {pct:3d}%  ({written/1e9:5.2f} / {total/1e9:5.2f} GB)")
    print(f"  done: {filename}")


def main() -> None:
    print(f"Destination: {DEST}")
    print()
    for fid, name, size in FILES:
        download(fid, name, size)
    print()
    print("All v62 files downloaded.")


if __name__ == "__main__":
    main()
