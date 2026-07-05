"""
Kinection local job daemon.

Polls the Cloudflare Worker API for queued jobs, runs the analysis
pipeline locally (where the modern DNA file lives), and uploads results
to R2.

Required environment variables (copy .env.example → .env and fill in):
  WORKER_API_URL   https://kinection-api.<your-subdomain>.workers.dev
  COMPUTE_API_KEY  same secret set via `wrangler secret put COMPUTE_API_KEY`
  USE_R2=1
  JOB_ID           (set automatically per job — do not set manually)
  R2_ENDPOINT_URL, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY, R2_BUCKET

Usage:
  python scripts/daemon.py
"""
from __future__ import annotations

import os
import sys
import time
import subprocess
import logging
import requests
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

WORKER_API_URL = os.environ["WORKER_API_URL"].rstrip("/")
COMPUTE_API_KEY = os.environ["COMPUTE_API_KEY"]
POLL_INTERVAL = int(os.environ.get("POLL_INTERVAL", "5"))

AUTH_HEADERS = {"Authorization": f"Bearer {COMPUTE_API_KEY}"}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
log = logging.getLogger(__name__)

# Defence-in-depth: scrub any raw genotype data from log records (Step 5.1.1).
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))
from utils.log_redact import RedactGenotypesFilter  # noqa: E402
from utils import r2_delete  # noqa: E402
from utils import r2_client  # noqa: E402
for h in logging.getLogger().handlers:
    h.addFilter(RedactGenotypesFilter())


def upload_key_for(job_id: str) -> str:
    """The R2 key where the modern raw DNA file for this job is stored.
    Matches the convention used by the Worker's presigned-PUT flow."""
    return f"uploads/{job_id}/raw.txt"


def claim_job(job_id: str) -> bool:
    """Mark the job as processing. Returns False if the update fails."""
    r = requests.patch(
        f"{WORKER_API_URL}/jobs/{job_id}/status",
        json={"status": "processing"},
        headers=AUTH_HEADERS,
        timeout=10,
    )
    return r.ok


def complete_job(job_id: str, error: str | None = None) -> None:
    status = "failed" if error else "complete"
    payload: dict = {"status": status}
    if error:
        payload["error"] = error
    requests.patch(
        f"{WORKER_API_URL}/jobs/{job_id}/status",
        json=payload,
        headers=AUTH_HEADERS,
        timeout=10,
    )


def run_analysis(job_id: str) -> None:
    env = {**os.environ, "JOB_ID": job_id, "USE_R2": "1"}

    # Download THIS job's uploaded raw DNA file from R2 and point the pipeline
    # at it via MODERN_DNA. The steps have no default file, so this is the only
    # way they get a modern individual. The on-disk copy is named generically
    # ('modern_individual.txt'), never after any person. Format (AncestryDNA vs
    # 23andMe) is auto-detected downstream.
    raw_key = upload_key_for(job_id)
    raw_local = r2_client.download_to_named_temp(raw_key, "modern_individual.txt")
    env["MODERN_DNA"] = str(raw_local)
    env["OUTPUT_LABEL"] = job_id[:8]     # per-job output dirs, avoid collisions
    log.info("[%s] fetched uploaded DNA file from R2 (%s)", job_id, raw_key)

    steps = [
        "step1_parse_harmonise.py",
        "step2_haplogroup.py",
        "step3_similarity_pca.py",
        "step1_4_tmrca.py",
        "step1_5_admixture.py",
        "step1_6_synthesis.py",
    ]
    try:
        for step in steps:
            log.info("[%s] running %s", job_id, step)
            result = subprocess.run(
                [sys.executable, str(SCRIPTS_DIR / step)],
                env=env,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"{step} exited {result.returncode}:\n{result.stderr[-2000:]}"
                )
    finally:
        # Never leave the plaintext raw DNA on local disk after the run.
        try:
            raw_local.unlink()
            raw_local.parent.rmdir()
        except Exception:
            pass


def post_deletion_receipt(job_id: str, receipt: r2_delete.DeletionReceipt) -> None:
    """POST the receipt back to the Worker so it persists to D1 (audit trail).
    Silent on transport error — the receipt also exists in our logs."""
    try:
        requests.post(
            f"{WORKER_API_URL}/jobs/{job_id}/deletion_receipt",
            json={
                "r2_key":     receipt.r2_key,
                "deleted_at": receipt.deleted_at,
                "verified":   receipt.verified,
                "attempts":   receipt.attempts,
                "reason":     receipt.reason,
                "requestor":  receipt.requestor,
            },
            headers=AUTH_HEADERS,
            timeout=10,
        )
    except Exception as e:
        log.warning("[%s] could not POST deletion receipt: %s", job_id, e)


def delete_raw_upload(job_id: str, *, reason: str) -> None:
    """Always called after analysis (success OR failure). Raw modern DNA
    must never linger in R2 after the job lifecycle ends."""
    key = upload_key_for(job_id)
    try:
        receipt = r2_delete.delete_and_verify(
            key, reason=reason, requestor="daemon",
        )
        log.info("[%s] raw upload deleted: %s (attempts=%d, verified=%s)",
                 job_id, key, receipt.attempts, receipt.verified)
        post_deletion_receipt(job_id, receipt)
    except RuntimeError as e:
        # Delete-and-verify exhausted retries — escalate.
        log.error("[%s] FAILED to delete raw upload %s: %s. Escalate.",
                  job_id, key, e)


def poll_once() -> None:
    r = requests.get(
        f"{WORKER_API_URL}/jobs",
        params={"status": "queued"},
        headers=AUTH_HEADERS,
        timeout=10,
    )
    r.raise_for_status()
    jobs = r.json()

    for job in jobs:
        job_id = job["id"]
        log.info("[%s] picked up queued job", job_id)

        if not claim_job(job_id):
            log.warning("[%s] could not claim job, skipping", job_id)
            continue

        analysis_error: str | None = None
        try:
            run_analysis(job_id)
        except Exception as exc:
            analysis_error = str(exc)
            log.error("[%s] analysis failed: %s", job_id, exc)
        finally:
            # ALWAYS delete the raw upload, regardless of analysis outcome.
            delete_raw_upload(
                job_id,
                reason="post_analysis" if analysis_error is None else "failed_analysis",
            )

        complete_job(job_id, error=analysis_error)
        log.info("[%s] %s", job_id, "complete" if analysis_error is None else "failed")


def main() -> None:
    log.info("Kinection daemon started (polling every %ds)", POLL_INTERVAL)
    while True:
        try:
            poll_once()
        except Exception as exc:
            log.error("poll error: %s", exc)
        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    main()
