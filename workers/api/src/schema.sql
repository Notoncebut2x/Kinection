-- Kinection D1 schema (v1)
-- Apply: wrangler d1 execute kinection --file=src/schema.sql
--        wrangler d1 execute kinection --file=src/schema.sql --remote
--
-- All tables are CREATE TABLE IF NOT EXISTS — the schema is idempotent and
-- safe to re-apply. Once we have production data, switch to numbered
-- migration files under workers/api/migrations/ (D1 supports `wrangler d1
-- migrations` natively). See ADR 0015.
--
-- Design principles:
--   - The raw modern DNA file MUST be deletable independently of analysis
--     outputs. Hence uploads (raw) and results (derived) are separate tables.
--   - Every R2 deletion of a raw upload must leave a verifiable receipt.
--   - User accounts are optional in v1 (anonymous jobs allowed); job rows
--     carry user_id when present, NULL otherwise.

------------------------------------------------------------------------
-- users
--   Minimal account model. Auth provider details live in auth_identities
--   so we can support email+password, Google, Apple, etc. without schema
--   churn. Consent timestamps are MANDATORY before any DNA upload.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS users (
  id              TEXT    PRIMARY KEY,            -- UUID
  email           TEXT    UNIQUE,                 -- nullable for anon flows
  created_at      INTEGER NOT NULL,               -- unix seconds
  consent_v       INTEGER,                        -- version of consent text accepted
  consent_at      INTEGER,                        -- unix seconds when consent was given
  deleted_at      INTEGER                         -- soft delete (GDPR Art 17)
);

CREATE INDEX IF NOT EXISTS idx_users_email ON users(email);

------------------------------------------------------------------------
-- auth_identities
--   One row per (provider, provider_user_id) tuple a user has linked.
--   Allows future OAuth providers without altering the users table.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS auth_identities (
  id                TEXT    PRIMARY KEY,
  user_id           TEXT    NOT NULL REFERENCES users(id) ON DELETE CASCADE,
  provider          TEXT    NOT NULL,             -- 'password' | 'google' | 'apple' | ...
  provider_user_id  TEXT    NOT NULL,             -- email for password; OAuth sub claim otherwise
  created_at        INTEGER NOT NULL,
  UNIQUE(provider, provider_user_id)
);

CREATE INDEX IF NOT EXISTS idx_auth_identities_user ON auth_identities(user_id);

------------------------------------------------------------------------
-- jobs
--   One analysis run. Lives until the user deletes their account or
--   manually deletes the job — outputs in R2 are deleted together.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS jobs (
  id              TEXT    PRIMARY KEY,            -- UUID, also used as R2 key prefix
  user_id         TEXT    REFERENCES users(id) ON DELETE SET NULL,
  status          TEXT    NOT NULL DEFAULT 'queued',
                  -- queued|uploading|processing|complete|failed|aborted|deleted
  label           TEXT,                           -- user-facing label, e.g. "rn"
  aadr_version    TEXT,                           -- e.g. "v66" — pinned at job start
  created_at      INTEGER NOT NULL,
  updated_at      INTEGER NOT NULL,
  started_at      INTEGER,                        -- when daemon picked it up
  completed_at    INTEGER,                        -- when results were written
  error           TEXT
);

CREATE INDEX IF NOT EXISTS idx_jobs_status      ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_user        ON jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_jobs_created_at  ON jobs(created_at);

------------------------------------------------------------------------
-- uploads
--   Tracks the raw modern DNA file for a job. Separate from `jobs` so a
--   job row can outlive its raw upload — once analysis completes, the
--   raw file is deleted but the results remain.
--
--   r2_key is always under uploads/<job_id>/raw.txt (or similar);
--   the daily reaper greps this table for r2_key references and asserts
--   no orphan R2 objects exist under uploads/.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS uploads (
  id              TEXT    PRIMARY KEY,            -- UUID
  job_id          TEXT    NOT NULL UNIQUE REFERENCES jobs(id) ON DELETE CASCADE,
  r2_key          TEXT    NOT NULL UNIQUE,        -- uploads/<job_id>/raw.txt
  size_bytes      INTEGER NOT NULL,
  sha256          TEXT    NOT NULL,               -- of the raw bytes (server-side computed)
  content_type    TEXT,
  format          TEXT,                           -- 'ancestrydna' | '23andme' | 'ftdna' | ...
  uploaded_at     INTEGER NOT NULL,
  deleted_at      INTEGER,                        -- NULL = still present in R2
  deletion_receipt_id TEXT REFERENCES deletion_receipts(id)
);

CREATE INDEX IF NOT EXISTS idx_uploads_deleted_at ON uploads(deleted_at);

------------------------------------------------------------------------
-- deletion_receipts
--   Proof that an R2 object was deleted. Created by the post-analysis
--   delete helper after HeadObject returns 404. The daily reaper also
--   writes receipts for objects it deletes (with reason='reaper').
--
--   Receipts are append-only — never updated, never deleted. They are
--   the audit trail for GDPR Article 17 compliance.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS deletion_receipts (
  id              TEXT    PRIMARY KEY,            -- UUID
  r2_key          TEXT    NOT NULL,
  requested_at    INTEGER NOT NULL,
  deleted_at      INTEGER NOT NULL,               -- when HeadObject confirmed 404
  reason          TEXT    NOT NULL,
                  -- 'post_analysis' | 'user_request' | 'reaper' | 'failed_analysis'
  requestor       TEXT,                           -- user_id | 'system' | 'reaper'
  job_id          TEXT    REFERENCES jobs(id) ON DELETE SET NULL,
  verified        INTEGER NOT NULL DEFAULT 0      -- 1 = HeadObject 404 confirmed
);

CREATE INDEX IF NOT EXISTS idx_receipts_r2_key   ON deletion_receipts(r2_key);
CREATE INDEX IF NOT EXISTS idx_receipts_job      ON deletion_receipts(job_id);
CREATE INDEX IF NOT EXISTS idx_receipts_reason   ON deletion_receipts(reason);

------------------------------------------------------------------------
-- results
--   Pointers to per-step output objects stored in R2. The step1.6
--   synthesis output (report.json + map_data.geojson) is the canonical
--   "the answer" — listed with step='synthesis'.
--
--   Storing pointers rather than the JSON itself keeps D1 lean and lets
--   the worker stream large geojsons directly from R2 to the browser
--   via signed URLs.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS results (
  id              TEXT    PRIMARY KEY,            -- UUID
  job_id          TEXT    NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
  step            TEXT    NOT NULL,
                  -- '1.1' | '1.2' | '1.3' | '1.4' | '1.5' | 'synthesis'
  artefact        TEXT    NOT NULL,
                  -- 'report.json' | 'map_data.geojson' | 'pairwise_distances.tsv' | ...
  r2_key          TEXT    NOT NULL UNIQUE,        -- outputs/<job_id>/<step>/<artefact>
  content_type    TEXT,
  size_bytes      INTEGER,
  created_at      INTEGER NOT NULL,
  UNIQUE(job_id, step, artefact)
);

CREATE INDEX IF NOT EXISTS idx_results_job  ON results(job_id);
CREATE INDEX IF NOT EXISTS idx_results_step ON results(job_id, step);

------------------------------------------------------------------------
-- audit_log
--   Security-relevant events: uploads accepted, presigned URLs minted,
--   deletions requested, failed delete retries, reaper actions. Used
--   for incident response and the daily reaper's "no orphans" check.
--
--   Not for application-level logging — those go to Worker logs.
------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS audit_log (
  id              TEXT    PRIMARY KEY,            -- UUID
  at              INTEGER NOT NULL,               -- unix seconds
  actor           TEXT,                           -- user_id | 'system' | 'reaper' | 'cron'
  event           TEXT    NOT NULL,
                  -- 'upload_url_minted' | 'upload_completed' | 'delete_requested'
                  -- | 'delete_verified' | 'delete_failed' | 'reaper_swept' | ...
  job_id          TEXT,
  r2_key          TEXT,
  detail          TEXT                            -- JSON blob, freeform
);

CREATE INDEX IF NOT EXISTS idx_audit_at        ON audit_log(at);
CREATE INDEX IF NOT EXISTS idx_audit_job       ON audit_log(job_id);
CREATE INDEX IF NOT EXISTS idx_audit_event_at  ON audit_log(event, at);
