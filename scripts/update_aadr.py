"""
AADR dataset update checker and R2 uploader.

Checks what version is currently stored in R2, probes the Reich Lab server
for newer releases, and streams any new files directly into R2 without
writing them to local disk.

On a successful upload it also patches scripts/utils/r2_client.py so the
analysis pipeline automatically uses the new version.

Usage:
  python scripts/update_aadr.py              # check and upload if newer
  python scripts/update_aadr.py --check      # print version info only
  python scripts/update_aadr.py --force      # re-upload even if current
  python scripts/update_aadr.py --version v63.0   # target a specific version

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

# Reich Lab public release directory
REICH_BASE = "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases"

# Manifest written to R2 to track current version
VERSION_MANIFEST_KEY = "dataset/current_version.json"

# Dataset file extensions to download (in order — smallest first)
DATASET_EXTS = ["ind", "snp", "anno", "geno"]

# Path to r2_client.py for automatic constant patching
R2_CLIENT_PATH = Path(__file__).parent / "utils" / "r2_client.py"


# ---------------------------------------------------------------------------
# URL / key helpers
# ---------------------------------------------------------------------------

def _parse_version(version: str) -> tuple[str, str]:
    """'v62.0' → ('62', '0')"""
    m = re.match(r"v?(\d+)\.(\d+)$", version)
    if not m:
        raise ValueError(f"Invalid version string: {version!r} — expected format like v62.0")
    return m.group(1), m.group(2)


def download_url(version: str, ext: str) -> str:
    major, minor = _parse_version(version)
    filename = f"v{major}.{minor}_1240k_public.{ext}"
    return f"{REICH_BASE}/V{major}/V{major}.{minor}/SHARE/public.dir/{filename}"


def r2_key(version: str, ext: str) -> str:
    major, minor = _parse_version(version)
    return f"dataset/v{major}.{minor}/v{major}.{minor}_1240k_public.{ext}"


# ---------------------------------------------------------------------------
# Version detection
# ---------------------------------------------------------------------------

def get_r2_version() -> str | None:
    """Return the version string currently stored in R2, or None."""
    try:
        data = json.loads(r2_client.get_object_text(VERSION_MANIFEST_KEY))
        return data.get("version")
    except Exception:
        return None


def probe_version(version: str) -> bool:
    """Return True if the .ind file for this version exists on the Reich Lab server."""
    url = download_url(version, "ind")
    try:
        r = requests.head(url, timeout=15, allow_redirects=True)
        return r.status_code == 200
    except Exception:
        return False


def detect_latest_version(known_current: str = "v62.0") -> str:
    """
    Probe the Reich Lab server for versions newer than known_current.
    Increments minor version then major version until the probe fails.
    Returns the highest confirmed version found.
    """
    major, minor = _parse_version(known_current)
    imajor, iminor = int(major), int(minor)

    latest = known_current
    for m in range(imajor, imajor + 5):
        for n in range(0 if m > imajor else iminor + 1, 10):
            candidate = f"v{m}.{n}"
            log.info("  probing %s ...", candidate)
            if probe_version(candidate):
                latest = candidate
                log.info("  found %s", candidate)
            else:
                break  # no more minor versions for this major
    return latest


# ---------------------------------------------------------------------------
# Download → R2 streaming
# ---------------------------------------------------------------------------

class _ProgressReader:
    """Wraps a urllib3 response body to log download progress."""

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
                log.info("  %s  %d%%  (%.1f / %.1f GB)",
                         self._label, pct,
                         self._read / 1e9, self._total / 1e9)
        return chunk

    def __getattr__(self, name):
        return getattr(self._raw, name)


def stream_to_r2(url: str, key: str) -> None:
    """Stream a URL directly to R2 via multipart upload — no local disk required."""
    client = r2_client.get_r2_client()
    filename = url.split("/")[-1]

    head = requests.head(url, timeout=15, allow_redirects=True)
    head.raise_for_status()
    total = int(head.headers.get("content-length", 0))
    log.info("uploading %s  →  r2://%s  (%.2f GB)", filename, key, total / 1e9)

    config = TransferConfig(
        multipart_threshold=128 * 1024 * 1024,
        multipart_chunksize=128 * 1024 * 1024,
        max_concurrency=4,
        use_threads=True,
    )

    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        wrapped = _ProgressReader(resp.raw, total, filename)
        client.upload_fileobj(wrapped, r2_client.R2_BUCKET, key, Config=config)

    log.info("done: %s", key)


# ---------------------------------------------------------------------------
# Version manifest + r2_client.py patch
# ---------------------------------------------------------------------------

def write_manifest(version: str) -> None:
    manifest = {
        "version": version,
        "uploaded_at": int(time.time()),
        "files": {ext: r2_key(version, ext) for ext in DATASET_EXTS},
    }
    r2_client.put_object(
        VERSION_MANIFEST_KEY,
        json.dumps(manifest, indent=2),
        "application/json",
    )
    log.info("manifest written: %s", VERSION_MANIFEST_KEY)


def patch_r2_client(version: str) -> None:
    """Update the DATASET_PREFIX and file key constants in r2_client.py."""
    major, minor = _parse_version(version)
    ver = f"v{major}.{minor}"

    source = R2_CLIENT_PATH.read_text()

    # Replace DATASET_PREFIX
    source = re.sub(
        r"(DATASET_PREFIX\s*=\s*)['\"]dataset/v[\d.]+['\"]",
        f"\\1'dataset/{ver}'",
        source,
    )
    # Replace version string in GENO/IND/ANNO/SNP_KEY references
    source = re.sub(
        r"(v\d+\.\d+)_1240k_public",
        f"{ver}_1240k_public",
        source,
    )

    R2_CLIENT_PATH.write_text(source)
    log.info("patched r2_client.py → %s", ver)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check for AADR updates and upload new versions to R2"
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Check version status only — do not download",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-upload even if R2 already has the target version",
    )
    parser.add_argument(
        "--version",
        help="Upload a specific version (e.g. v63.0) instead of auto-detecting",
    )
    args = parser.parse_args()

    current = get_r2_version()
    log.info("R2 current version : %s", current or "(none)")

    if args.version:
        target = args.version
        log.info("Target version     : %s (manual)", target)
    else:
        log.info("Probing Reich Lab for latest version...")
        target = detect_latest_version(current or "v62.0")
        log.info("Latest available   : %s", target)

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
    for ext in DATASET_EXTS:
        url = download_url(target, ext)
        key = r2_key(target, ext)
        stream_to_r2(url, key)

    write_manifest(target)
    patch_r2_client(target)
    log.info("All done. R2 and r2_client.py are now at %s.", target)


if __name__ == "__main__":
    main()
