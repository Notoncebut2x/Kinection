"""
Daily reaper — sweeps orphaned modern DNA uploads from R2.

The Step 5.1.1 lifecycle relies on inline deletion after analysis. The
reaper is the backstop: any object under `uploads/` older than
REAPER_MAX_AGE_HOURS that's still in R2 means the inline deletion failed
(crashed job, network blip, code bug). The reaper deletes it and emits a
receipt with reason='reaper'.

Run on a schedule (cron / Cloudflare Cron Trigger):

    python scripts/reaper.py

Set REAPER_DRY_RUN=1 to log what would be deleted without actually
deleting. Recommended for the first few real runs.

Future: when the worker has D1 access, the reaper should join against
the `uploads` and `deletion_receipts` tables before deleting (avoid
deleting an upload whose row says it's still actively being processed).
For v1 the age-based heuristic is acceptable because no job should
legitimately have a raw upload sitting in R2 for >24h.
"""
from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from dotenv import load_dotenv  # noqa: E402
load_dotenv(ROOT / ".env")

from utils import r2_delete  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("reaper")


REAPER_MAX_AGE_HOURS = float(os.environ.get("REAPER_MAX_AGE_HOURS", "24"))
REAPER_DRY_RUN       = os.environ.get("REAPER_DRY_RUN", "").lower() in ("1", "true", "yes")


def main() -> int:
    log.info("Reaper starting (max_age_hours=%.1f, dry_run=%s)",
             REAPER_MAX_AGE_HOURS, REAPER_DRY_RUN)

    orphans = r2_delete.list_orphan_upload_keys(REAPER_MAX_AGE_HOURS)
    log.info("Found %d candidate orphan upload object(s)", len(orphans))

    if not orphans:
        log.info("Nothing to reap.")
        return 0

    failures: list[str] = []
    for key in orphans:
        if REAPER_DRY_RUN:
            log.info("[DRY RUN] would delete %s", key)
            continue
        try:
            receipt = r2_delete.delete_and_verify(
                key, reason="reaper", requestor="reaper",
            )
            log.info("Reaped %s (attempts=%d, verified=%s)",
                     key, receipt.attempts, receipt.verified)
        except RuntimeError as e:
            log.error("REAPER FAILED to delete %s: %s", key, e)
            failures.append(key)

    if failures:
        log.error("Reaper finished with %d FAILURE(s) — escalate to on-call.",
                  len(failures))
        return 1

    log.info("Reaper finished cleanly at %s.", time.strftime("%Y-%m-%dT%H:%M:%S"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
