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

SCRIPTS_DIR = Path(__file__).parent


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
    steps = [
        "step1_parse_harmonise.py",
        "step2_haplogroup.py",
        "step3_similarity_pca.py",
    ]
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

        try:
            run_analysis(job_id)
            complete_job(job_id)
            log.info("[%s] complete", job_id)
        except Exception as exc:
            log.error("[%s] failed: %s", job_id, exc)
            complete_job(job_id, error=str(exc))


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
