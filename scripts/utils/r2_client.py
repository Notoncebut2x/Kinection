"""
Cloudflare R2 client for the Kinection analysis pipeline.

R2 is S3-compatible, so we use boto3 pointed at the R2 endpoint.

Required environment variables:
  R2_ENDPOINT_URL        https://<account_id>.r2.cloudflarestorage.com
  R2_ACCESS_KEY_ID       R2 API token access key ID
  R2_SECRET_ACCESS_KEY   R2 API token secret
  R2_BUCKET              bucket name (default: kinection)

Set USE_R2=1 in any analysis script to activate R2 mode.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import boto3
from botocore.config import Config

# ---------------------------------------------------------------------------
# Bucket and object key layout
# ---------------------------------------------------------------------------

R2_BUCKET = os.environ.get('R2_BUCKET', 'kinection')

# Reference dataset — uploaded once by the operator
DATASET_PREFIX = 'dataset/v62.0'
GENO_KEY = f'{DATASET_PREFIX}/v62.0_1240k_public.geno'
IND_KEY  = f'{DATASET_PREFIX}/v62.0_1240k_public.ind'
ANNO_KEY = f'{DATASET_PREFIX}/v62.0_1240k_public.anno'
SNP_KEY  = f'{DATASET_PREFIX}/v62.0_1240k_public.snp'

# Haplogroup marker databases — uploaded once, also cached in Workers KV
YDNA_MARKERS_KEY  = 'markers/ydna_markers.json'
MTDNA_MARKERS_KEY = 'markers/mtdna_markers.json'


def upload_key(job_id: str) -> str:
    """R2 key for a user's raw DNA upload."""
    return f'uploads/{job_id}/raw.txt'


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
