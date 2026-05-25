# ADR-0015: D1 schema v1 (Phase 2.3)

* Status: Accepted
* Date: 2026-05-25
* Supersedes: —
* Phase: 2 (Web platform — Step 2.3)

## Context

Phase 2.3 of the workplan calls for designing the application database schema. Per ADRs 0011 and 0012 the stack is Cloudflare Workers + R2 + D1 (SQLite-on-the-edge), superseding the original FastAPI + Postgres choice (ADRs 0007/0009).

The schema must support:
- **Anonymous and authenticated jobs.** v1 lets users analyse without an account; v2 will require email + consent for retention.
- **Independent lifecycle for raw modern DNA.** The raw upload must be deletable (with verifiable receipt) while keeping analysis outputs. This is the Step 5.1.1 requirement.
- **GDPR Article 17 (right to deletion).** Receipts are append-only proof of compliance.
- **Audit trail.** Security-relevant events need persistent storage for incident response.
- **Multi-step pipeline outputs.** Step 1.1 → 1.6 each produce artefacts; the worker needs to enumerate them per job.

## Decision

Seven tables, all defined in `workers/api/src/schema.sql`:

| Table | Purpose | Key invariant |
|---|---|---|
| `users` | Account identity + consent | `consent_at` required before any upload |
| `auth_identities` | Pluggable auth providers (password, OAuth) | One row per `(provider, provider_user_id)` |
| `jobs` | One analysis run | Status state machine: `queued → uploading → processing → complete/failed/aborted/deleted` |
| `uploads` | Raw modern DNA file pointer | `r2_key` is `uploads/<job_id>/raw.txt`; `deleted_at` set when R2 object is removed |
| `deletion_receipts` | Append-only proof of R2 object deletion | Never updated, never deleted; `verified=1` only after HeadObject 404 |
| `results` | Per-step output object pointers | Pointers only — JSON/TSV/geojson bytes live in R2, streamed to browser via signed URLs |
| `audit_log` | Security-relevant events | Append-only; written by Worker on every privileged action |

### Why pointers in `results`, not blobs

The Step 1.6 synthesis output (`map_data.geojson`) for a typical job is 50–500 KB, and `pairwise_distances.tsv` is ~2 MB. D1 has a 100 KB row size limit and pays per-byte-read on every query. Storing R2 keys in D1 and streaming the actual bytes via signed GETs from R2 is cheaper, faster, and cleaner.

### Why separate `uploads` from `jobs`

If `uploads` was a column on `jobs`, the post-analysis raw-file deletion would either need to mutate the job row (confusing — the job is still "complete") or use a NULL sentinel (also confusing). A separate table with its own `deleted_at` makes the lifecycle explicit and lets `uploads.deletion_receipt_id` reference the receipt row directly.

### Why a single `audit_log` table

For v1, "everything security-relevant in one append-only log" is simpler than fragmenting by event type. The `idx_audit_event_at` index makes per-event-type queries cheap. If volume becomes a problem we can move high-volume events to R2 with a daily JSONL roll-up.

### Migration strategy

For v1 (pre-production), the schema is idempotent (`CREATE TABLE IF NOT EXISTS`) and re-applied via `wrangler d1 execute --file=src/schema.sql`. Acceptable because:
- No production data yet — destructive changes are cheap.
- The schema is small enough to read top-to-bottom and understand.

Once we have real user data, switch to numbered migration files under `workers/api/migrations/` and use `wrangler d1 migrations apply`. This ADR will be amended at that point.

## Consequences

**Positive:**
- The schema supports the Step 5.1.1 lifecycle requirements out of the gate.
- Anonymous jobs work (`jobs.user_id` is nullable) for early/dev usage.
- Auth provider expansion (adding Google/Apple) requires only inserts into `auth_identities`, no schema change.
- The reaper has a definitive source of truth (`uploads.deleted_at IS NULL` + `jobs.created_at`) for finding orphans.

**Negative:**
- Anonymous jobs can't be reclaimed if a user later signs up — `user_id` is only set at job creation.
  - Mitigated: in v2, add a "claim job" flow that updates `user_id` when the user authenticates within N minutes of job creation.
- The idempotent `CREATE TABLE IF NOT EXISTS` pattern silently no-ops if a column is added later. Switch to numbered migrations before production data lands.
- No row-level encryption. R2 SSE-C (per ADR 0014... or rather the lifecycle ADR when written) protects the raw bytes; D1 rows themselves rely on Cloudflare's at-rest encryption.

## Alternatives considered

**A. Postgres / Supabase.**
Original Phase 2 plan (ADR 0009 — now Proposed/superseded). Postgres has richer types and migrations tooling, but adds a non-Cloudflare hosting dependency and per-request latency (edge → origin). For a workload that's mostly key-value lookups and small joins, D1 wins.

**B. KV-only persistence.**
Workers KV is eventually consistent and doesn't support the kind of relational joins (e.g. "list all jobs for user X with their uploads' deletion status") that the security and admin flows need. D1 is the right primitive.

**C. Inline result JSON in D1.**
Would simplify the worker (one query for "everything about job X") but blows up D1 storage, query latency, and the 100 KB row limit. Pointers-to-R2 is the right tradeoff.

## References

- D1 documentation: https://developers.cloudflare.com/d1/
- D1 migrations: https://developers.cloudflare.com/d1/reference/migrations/
- ADR 0011 — Cloudflare R2 for genotype storage
- ADR 0012 — Cloudflare Workers as the API layer
