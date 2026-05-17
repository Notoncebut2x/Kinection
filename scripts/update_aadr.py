"""
AADR dataset update checker and R2 uploader.

Uses the Harvard Dataverse API to detect the latest AADR release, then
streams the 1240K public dataset files directly into R2 without writing
them to local disk. Writes a manifest at dataset/current_version.json
which r2_client.py reads at runtime — restart the daemon to pick up
the new version.

Usage:
  python scripts/update_aadr.py              # check and upload if newer
  python scripts/update_aadr.py --check      # print version info only
  python scripts/update_aadr.py --force      # re-upload even if current

Required environment (copy .env.example → .env):
  R2_ENDPOINT_URL, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY, R2_BUCKET
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from pathlib import Path

import requests
from boto3.s3.transfer import TransferConfig
from dotenv import load_dotenv

load_dotenv()

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import r2_client

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Harvard Dataverse — official AADR repository
DATAVERSE_API = "https://dataverse.harvard.edu/api"
DATASET_DOI   = "doi:10.7910/DVN/FFIDCW"

# Match only the 1240K public files (not 2M, HO, compatibility, or MT variants)
AADR_FILE_RE = re.compile(
    r"^v\d+\.1240K\.aadr\.PUB\.(geno|ind|snp|anno)$",
    re.IGNORECASE,
)

# Manifest key in R2 — tracks current version and file paths.
# r2_client.py reads this at runtime, so no source patching is needed.
VERSION_MANIFEST_KEY = "dataset/current_version.json"


# ---------------------------------------------------------------------------
# Dataverse API helpers
# ---------------------------------------------------------------------------

def _dataverse_get(path: str, **params) -> dict:
    # Build query string manually — Dataverse requires unencoded colons and slashes in the DOI
    qs = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"{DATAVERSE_API}/{path}?{qs}" if qs else f"{DATAVERSE_API}/{path}"
    r = requests.get(url, timeout=30, allow_redirects=True)
    r.raise_for_status()
    return r.json()


def get_latest_files() -> list[dict]:
    """
    Return file metadata for the 1240K PUB files in the latest Dataverse version.
    Each dict has keys: filename, file_id, size_bytes.
    """
    data = _dataverse_get(
        "datasets/:persistentId/versions/:latest/files",
        persistentId=DATASET_DOI,
    )
    results = []
    for entry in data["data"]:
        df = entry["dataFile"]
        label = entry.get("label", df.get("filename", ""))
        if AADR_FILE_RE.match(label):
            results.append({
                "filename": label,
                "file_id": df["id"],
                "size_bytes": df.get("filesize", 0),
            })
    return results


def extract_version(files: list[dict]) -> str:
    """Extract AADR version string (e.g. 'v66') from filenames."""
    for f in files:
        m = re.match(r"(v\d+)", f["filename"], re.IGNORECASE)
        if m:
            return m.group(1).lower()
    raise ValueError("Could not determine AADR version from Dataverse filenames")


def dataverse_download_url(file_id: int) -> str:
    return f"{DATAVERSE_API}/access/datafile/{file_id}"


# ---------------------------------------------------------------------------
# R2 version manifest
# ---------------------------------------------------------------------------

def get_r2_version() -> str | None:
    """Return the AADR version string currently in R2, or None."""
    try:
        data = json.loads(r2_client.get_object_text(VERSION_MANIFEST_KEY))
        return data.get("version")
    except Exception:
        return None


def write_manifest(version: str, files: list[dict]) -> None:
    manifest = {
        "version": version,
        "uploaded_at": int(time.time()),
        "files": {
            Path(f["filename"]).suffix.lstrip("."): f"dataset/{version}/{f['filename']}"
            for f in files
        },
    }
    r2_client.put_object(
        VERSION_MANIFEST_KEY,
        json.dumps(manifest, indent=2),
        "application/json",
    )
    log.info("manifest written → %s", VERSION_MANIFEST_KEY)


# ---------------------------------------------------------------------------
# Streaming upload
# ---------------------------------------------------------------------------

class _ProgressReader:
    def __init__(self, raw, total: int, label: str):
        self._raw = raw
        self._total = total
        self._read = 0
        self._label = label
        self._last_pct = -1

    def read(self, amt: int = -1) -> bytes:
        chunk = self._raw.read(amt)
        self._read += len(chunk)
        if self._total:
            pct = int(self._read / self._total * 100)
            if pct // 10 != self._last_pct // 10:
                self._last_pct = pct
                log.info("  %s  %d%%  (%.2f / %.2f GB)",
                         self._label, pct,
                         self._read / 1e9, self._total / 1e9)
        return chunk

    def __getattr__(self, name):
        return getattr(self._raw, name)


def stream_to_r2(file_id: int, filename: str, version: str, size_bytes: int) -> None:
    """Stream a Dataverse file directly to R2 via multipart upload."""
    client = r2_client.get_r2_client()
    url = dataverse_download_url(file_id)
    key = f"dataset/{version}/{filename}"

    log.info("uploading %s  →  r2://%s  (%.2f GB)", filename, key, size_bytes / 1e9)

    config = TransferConfig(
        multipart_threshold=128 * 1024 * 1024,
        multipart_chunksize=128 * 1024 * 1024,
        max_concurrency=4,
        use_threads=True,
    )

    with requests.get(url, stream=True, timeout=120, allow_redirects=True) as resp:
        resp.raise_for_status()
        wrapped = _ProgressReader(resp.raw, size_bytes, filename)
        client.upload_fileobj(wrapped, r2_client.R2_BUCKET, key, Config=config)

    log.info("done: %s", key)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check for AADR updates on Harvard Dataverse and upload to R2"
    )
    parser.add_argument("--check", action="store_true",
                        help="Check version status only — do not download")
    parser.add_argument("--force", action="store_true",
                        help="Re-upload even if R2 already has the latest version")
    args = parser.parse_args()

    current = get_r2_version()
    log.info("R2 current version  : %s", current or "(none)")

    log.info("Querying Harvard Dataverse for latest AADR release...")
    files = get_latest_files()
    if not files:
        log.error("No 1240K PUB files found on Dataverse — check the API or dataset DOI.")
        sys.exit(1)

    target = extract_version(files)
    log.info("Latest on Dataverse : %s", target)
    log.info("Files to upload:")
    for f in sorted(files, key=lambda x: x["size_bytes"]):
        log.info("  %-45s  %.2f GB", f["filename"], f["size_bytes"] / 1e9)

    if args.check:
        if current == target:
            log.info("R2 is up to date.")
        else:
            log.info("Update available: %s → %s", current, target)
        return

    if current == target and not args.force:
        log.info("R2 is already at %s. Use --force to re-upload.", target)
        return

    log.info("Starting upload for %s ...", target)
    # Upload smallest files first (ind, snp, anno) then geno last
    for f in sorted(files, key=lambda x: x["size_bytes"]):
        stream_to_r2(f["file_id"], f["filename"], target, f["size_bytes"])

    write_manifest(target, files)
    log.info("All done. R2 is now at %s. Restart the daemon to pick up the new version.", target)


if __name__ == "__main__":
    main()
