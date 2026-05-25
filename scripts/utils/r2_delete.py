"""
R2 deletion helper for raw modern DNA uploads.

Two guarantees the Step 5.1.1 lifecycle requires:
  1. After a successful analysis, the raw upload at uploads/<job_id>/raw.txt
     is DeleteObject'd from R2.
  2. Deletion is *verified* — HeadObject must return 404 before we mark the
     job as deleted. A 200 means the object is still recoverable.

When verification fails, retry with exponential backoff. After N retries,
raise — the worker / daemon escalates to an on-call alert. We never claim
deletion success without HeadObject confirmation.

See ADR 0015 (D1 schema — deletion_receipts table).
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass

from botocore.exceptions import ClientError

from . import r2_client

log = logging.getLogger(__name__)

DEFAULT_RETRIES = 4
DEFAULT_BACKOFF_BASE = 1.0  # seconds; doubles each retry


@dataclass(frozen=True)
class DeletionReceipt:
    """Proof of an R2 deletion, suitable for persisting to D1."""
    r2_key:         str
    deleted_at:     int     # unix seconds
    verified:       bool    # True iff HeadObject returned 404
    attempts:       int
    reason:         str     # 'post_analysis' | 'user_request' | 'reaper' | 'failed_analysis'
    requestor:      str     # user_id | 'system' | 'reaper'


def delete_and_verify(
    key: str,
    *,
    reason: str = "post_analysis",
    requestor: str = "system",
    retries: int = DEFAULT_RETRIES,
    backoff_base: float = DEFAULT_BACKOFF_BASE,
) -> DeletionReceipt:
    """
    Delete an R2 object and verify it is gone via HeadObject.

    Returns a DeletionReceipt with verified=True on success.

    Raises RuntimeError after `retries` consecutive verification failures so
    the caller (worker / daemon) escalates rather than silently moving on.
    """
    client = r2_client.get_r2_client()
    bucket = r2_client.R2_BUCKET

    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            client.delete_object(Bucket=bucket, Key=key)
        except ClientError as e:
            last_exc = e
            log.warning("DeleteObject %s failed on attempt %d: %s", key, attempt, e)
            _sleep_backoff(attempt, backoff_base)
            continue

        # Verify
        try:
            client.head_object(Bucket=bucket, Key=key)
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "")
            if code in ("404", "NoSuchKey", "NotFound"):
                log.info("DELETE verified: %s (attempt %d, reason=%s)", key, attempt, reason)
                return DeletionReceipt(
                    r2_key=key,
                    deleted_at=int(time.time()),
                    verified=True,
                    attempts=attempt,
                    reason=reason,
                    requestor=requestor,
                )
            last_exc = e
            log.warning("HeadObject %s after delete returned unexpected error %s",
                        key, code)
        else:
            log.warning("HeadObject %s still found after DeleteObject (attempt %d)",
                        key, attempt)
        _sleep_backoff(attempt, backoff_base)

    raise RuntimeError(
        f"Failed to delete-and-verify {key} after {retries} attempts: {last_exc}"
    )


def _sleep_backoff(attempt: int, base: float) -> None:
    time.sleep(base * (2 ** (attempt - 1)))


def list_orphan_upload_keys(max_age_hours: float = 24.0) -> list[str]:
    """
    Used by the daily reaper. Returns R2 keys under uploads/ older than
    max_age_hours that should not exist (analysis is long complete or
    abandoned). Does NOT cross-reference D1 — that join happens in the
    reaper script which has the D1 binding.

    The daemon side just enumerates R2; the Worker side joins against D1
    rows (uploads.deleted_at IS NULL) to find true orphans.
    """
    client = r2_client.get_r2_client()
    bucket = r2_client.R2_BUCKET
    cutoff = time.time() - max_age_hours * 3600

    orphans: list[str] = []
    paginator = client.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix="uploads/"):
        for obj in page.get("Contents", []):
            last_modified = obj["LastModified"].timestamp()
            if last_modified < cutoff:
                orphans.append(obj["Key"])
    return orphans
