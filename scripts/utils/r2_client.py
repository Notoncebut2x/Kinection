"""
Cloudflare R2 client for the Kinection analysis pipeline.

R2 is S3-compatible, so we use boto3 pointed at the R2 endpoint.

Required environment variables:
  R2_ENDPOINT_URL        https://<account_id>.r2.cloudflarestorage.com
  R2_ACCESS_KEY_ID       R2 API token access key ID
  R2_SECRET_ACCESS_KEY   R2 API token secret
  R2_BUCKET              bucket name (default: kinection)

Set USE_R2=1 in any analysis script to activate R2 mode.

The current AADR version and file paths are loaded from
`dataset/current_version.json` in R2 (written by scripts/update_aadr.py).
Accessing `r2_client.GENO_KEY` (etc.) reads that manifest on first use and
caches it for the process lifetime — restart the daemon to pick up a new
AADR release.
"""
from __future__ import annotations

import json
import os
import tempfile
from functools import lru_cache
from pathlib import Path

import boto3
from botocore.config import Config

# ---------------------------------------------------------------------------
# Bucket and fixed paths
# ---------------------------------------------------------------------------

R2_BUCKET = os.environ.get('R2_BUCKET', 'kinection')

# Manifest written by update_aadr.py — tracks the currently uploaded AADR version
VERSION_MANIFEST_KEY = 'dataset/current_version.json'

# Haplogroup marker databases — uploaded once, also cached in Workers KV
YDNA_MARKERS_KEY  = 'markers/ydna_markers.json'
MTDNA_MARKERS_KEY = 'markers/mtdna_markers.json'

# Modern individual DNA files are NEVER stored in R2 — local filesystem only.


# Fallback paths used if the manifest is missing (fresh R2, pre-update_aadr).
# Set to v66 to match the current canonical cloud version; the actual key
# resolution always prefers the manifest at dataset/current_version.json.
_DEFAULT_VERSION = 'v66'
_DEFAULT_KEYS = {
    'version': _DEFAULT_VERSION,
    'prefix':  f'dataset/{_DEFAULT_VERSION}',
    'geno':    f'dataset/{_DEFAULT_VERSION}/{_DEFAULT_VERSION}.1240K.aadr.PUB.geno',
    'ind':     f'dataset/{_DEFAULT_VERSION}/{_DEFAULT_VERSION}.1240K.aadr.PUB.ind',
    'snp':     f'dataset/{_DEFAULT_VERSION}/{_DEFAULT_VERSION}.1240K.aadr.PUB.snp',
    'anno':    f'dataset/{_DEFAULT_VERSION}/{_DEFAULT_VERSION}.1240K.aadr.PUB.anno',
}


def output_key(job_id: str, filename: str) -> str:
    """R2 key for a per-job analysis output file."""
    return f'outputs/{job_id}/{filename}'


# ---------------------------------------------------------------------------
# Client factory
# ---------------------------------------------------------------------------

def get_r2_client():
    """Return a boto3 S3 client pointed at Cloudflare R2."""
    return boto3.client(
        's3',
        endpoint_url=os.environ['R2_ENDPOINT_URL'],
        aws_access_key_id=os.environ['R2_ACCESS_KEY_ID'],
        aws_secret_access_key=os.environ['R2_SECRET_ACCESS_KEY'],
        config=Config(signature_version='s3v4'),
        region_name='auto',
    )


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def get_object_bytes(key: str, byte_range: tuple[int, int] | None = None) -> bytes:
    """Download an R2 object (or a byte range of it) into memory."""
    kwargs: dict = {'Bucket': R2_BUCKET, 'Key': key}
    if byte_range is not None:
        kwargs['Range'] = f'bytes={byte_range[0]}-{byte_range[1]}'
    return get_r2_client().get_object(**kwargs)['Body'].read()


def get_object_text(key: str) -> str:
    """Download an R2 text object as a UTF-8 string."""
    return get_object_bytes(key).decode('utf-8', errors='replace')


def put_object(key: str, data: bytes | str,
               content_type: str = 'application/octet-stream') -> None:
    """Upload bytes or a string to R2."""
    if isinstance(data, str):
        data = data.encode('utf-8')
    get_r2_client().put_object(
        Bucket=R2_BUCKET, Key=key, Body=data, ContentType=content_type,
    )


def upload_file(local_path: Path | str, key: str) -> None:
    """Upload a local file to R2."""
    get_r2_client().upload_file(str(local_path), R2_BUCKET, key)


def download_to_temp(key: str, suffix: str = '') -> Path:
    """
    Download an R2 object to a temporary local file and return its path.
    Caller must delete the file when done (use in a try/finally block).
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.close()
    get_r2_client().download_file(R2_BUCKET, key, tmp.name)
    return Path(tmp.name)


def download_to_named_temp(key: str, filename: str) -> Path:
    """
    Download an R2 object into a fresh temp directory under a chosen, generic
    filename (e.g. 'modern_individual.txt') — so the on-disk working copy is
    never named after any specific person. Returns the file path; caller should
    delete the file and its parent directory when done.
    """
    d = Path(tempfile.mkdtemp(prefix='kinection_'))
    dest = d / filename
    get_r2_client().download_file(R2_BUCKET, key, str(dest))
    return dest


# ---------------------------------------------------------------------------
# Dynamic dataset version resolution
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def get_dataset_keys() -> dict:
    """
    Return current dataset key paths read from the R2 version manifest.
    Falls back to defaults if the manifest is missing or unreadable.

    Cached for the process lifetime — restart the daemon to pick up a
    newly uploaded AADR version.
    """
    try:
        manifest = json.loads(get_object_text(VERSION_MANIFEST_KEY))
        version = manifest['version']
        files = manifest['files']
        return {
            'version': version,
            'prefix':  f'dataset/{version}',
            'geno':    files['geno'],
            'ind':     files['ind'],
            'snp':     files['snp'],
            'anno':    files['anno'],
        }
    except Exception:
        return dict(_DEFAULT_KEYS)


# PEP 562: lazy module-level attributes resolve on first access.
# Existing code using `r2_client.GENO_KEY` etc. keeps working unchanged.
_DYNAMIC = {
    'GENO_KEY':       lambda: get_dataset_keys()['geno'],
    'IND_KEY':        lambda: get_dataset_keys()['ind'],
    'SNP_KEY':        lambda: get_dataset_keys()['snp'],
    'ANNO_KEY':       lambda: get_dataset_keys()['anno'],
    'DATASET_PREFIX': lambda: get_dataset_keys()['prefix'],
}


def __getattr__(name: str):
    if name in _DYNAMIC:
        return _DYNAMIC[name]()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
