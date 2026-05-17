-- Kinection D1 schema
-- Run: wrangler d1 execute kinection --file=src/schema.sql

CREATE TABLE IF NOT EXISTS jobs (
  id            TEXT    PRIMARY KEY,
  status        TEXT    NOT NULL DEFAULT 'queued',  -- queued|processing|complete|failed
  upload_key    TEXT    NOT NULL,
  result_prefix TEXT    NOT NULL,
  created_at    INTEGER NOT NULL,
  updated_at    INTEGER NOT NULL,
  error         TEXT
);

CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
