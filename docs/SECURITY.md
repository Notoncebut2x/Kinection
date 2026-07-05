# Kinection — Security Review

*Last updated: 2026-05-17*

This document records the security review of the Kinection platform, distinguishing what has been **fixed**, what is **acceptable for current single-user use**, and what is **open work blocking public deployment**.

For the architecture this review applies to, see [`ARCHITECTURE.md`](ARCHITECTURE.md).

---

## Threat model

The data being protected is **personal genetic data** — fundamentally non-rotatable. A leak is permanent. The threats we worry about, in priority order:

1. Accidental commit of raw DNA files or derived genetic outputs to a public repository.
2. Leakage of cloud credentials (R2, Worker secrets) granting attacker write access to the dataset or read access to other users' outputs.
3. Unauthorized access to *other users'* result files via guessable job IDs or open endpoints.
4. Pipeline scripts being tricked into uploading raw DNA to the cloud.
5. Denial-of-service or quota exhaustion of cloud resources.

---

## Fixed — completed during the review

### 1. Raw DNA never leaves the user's machine

**Status:** Enforced in code and config.

- `r2_client.py` exposes no `upload_key()` helper for raw DNA. The only R2 keys defined are dataset paths and per-job *output* paths.
- The Worker API has no route that *receives* raw file bytes. The browser uploads directly to per-job R2 storage via a short-lived presigned PUT URL (minted by `POST /uploads/url`), so the raw file never traverses the Worker; the daemon downloads it, analyses it, and deletes it with a verified receipt.
- The pipeline scripts read the modern file only from the `MODERN_DNA` path (set by `run_local.py --dna` locally, or by the daemon after downloading the per-job upload to a generic `modern_individual.txt`). There is **no default/bundled DNA file** — the scripts error out if `MODERN_DNA` is unset, so they never silently analyse anyone else's DNA.
- `.gitignore` covers `data/input_data/` (raw inputs) and `output/` (derived data).

### 2. Personal genetic data purged from git history

**Status:** Removed.

A historical accident had committed the `output/` directory — containing encoded genotype arrays, haplogroup reports, and pairwise genetic distances for two real individuals — to the public GitHub repository. The directory was removed from current state, **and from all historical commits** via `git filter-branch`, and force-pushed. The `output/` path is now in `.gitignore`.

### 3. Exposed credentials rotated

**Status:** Rotated.

R2 access keys and `COMPUTE_API_KEY` were exposed in a shared chat session during initial setup. All affected credentials were destroyed and reissued:
- Old R2 API token deleted, new one generated with `Object Read & Write` scoped to the `kinection` bucket.
- `COMPUTE_API_KEY` re-generated via `openssl rand -hex 32` and updated on both the Worker (via `wrangler secret put`) and in `.env`.

### 4. `.env` files cannot be committed

**Status:** Enforced in `.gitignore`.

`.env`, `.env.local`, and `.dev.vars` are all listed. The template `.env.example` is the only env-pattern file in version control.

### 5. Wrangler local state excluded

**Status:** Enforced.

`.wrangler/` directories — which contain local D1 mirrors and could include test data — are gitignored.

---

## Acceptable for current state — single-user, pre-public

These are documented as conscious trade-offs:

### CORS is permissive (`*`)

`Access-Control-Allow-Origin: *` on the Worker is fine while there is no production frontend. Before public deployment, lock this to the deployed frontend origin.

### Result fetch is unauthenticated

`GET /jobs/:id/results/:filename` requires only the job UUID. Anyone with the UUID can fetch the results. For a single-user deployment where job IDs aren't shared this is acceptable; the UUID is unguessable (122 bits of entropy). For multi-user / public deployment, see "Open work" below.

### `POST /jobs` is unauthenticated

Anyone can create a job in D1 (which the daemon will then try to run). At single-user scale, with a daemon that only reads *your* DNA file, an attacker creating spurious jobs achieves nothing except making the daemon do empty work. At multi-user scale this becomes a DoS vector.

### No retention policy

Job rows in D1 and output files in R2 persist forever. Acceptable while you control all the data; not acceptable for users you don't know.

### Daemon stores credentials in `.env`

Standard practice. For a hosted deployment you'd want a secrets manager. For local development this is fine.

---

## Open work — required before public deployment

### O1. Authenticate result downloads

`GET /jobs/:id/results/:filename` must be tied to the user who created the job. Options:

- Issue a session token at job creation, require it for result fetch.
- Use Cloudflare Access in front of the Worker.
- Issue short-lived R2 signed URLs at fetch time, never streaming through the Worker.

### O2. Rate-limit `POST /jobs`

Either a per-IP token bucket via the Worker's `caches` API, or front the endpoint with Cloudflare Turnstile / Rate Limiting Rules.

### O3. Tighten CORS

Replace `Access-Control-Allow-Origin: *` with the production frontend origin. Restrict allowed methods to those actually used.

### O4. Add security response headers

The Worker currently sends no security headers. Add at least:
- `Strict-Transport-Security: max-age=31536000`
- `X-Content-Type-Options: nosniff`
- `Referrer-Policy: no-referrer`
- A reasonable `Content-Security-Policy` (when the frontend exists)

### O5. Input validation on PATCH body

`handleUpdateStatus` reads JSON without size limits or shape validation. Add a content-length cap and a Zod-style schema check on `{status, error}`.

### O6. Retention and deletion ADR

Document the policy: how long jobs stay in D1, how long outputs stay in R2, whether users can request deletion. Then implement it (a scheduled Worker cron job, or a TTL on R2 keys).

### O7. PII handling documentation

A user-facing privacy notice covering:
- What data is collected (modern DNA stays local; derived analysis outputs only).
- Where that data lives (R2, R2 region, encryption at rest).
- How long it's retained.
- How to delete it.

This is GDPR-grade and is mentioned in [ADR-0009](decisions/0009-database-postgresql.md) as Phase 5 work.

### O8. Dependency hygiene

There is no automated dependency scanning. Add at minimum:
- `pip-audit` in CI for Python deps
- `npm audit` (or Dependabot) for the Worker
- Pin versions in `requirements.txt` (currently dependencies are installed ad-hoc)

### O9. Pre-commit hook

A versioned hook at `.githooks/pre-commit` that blocks commits containing:
- Files under `output/`, `data/input_data/`, or matching `.env*` (except `.env.example`)
- Plain-text patterns resembling API keys (`R2_SECRET_ACCESS_KEY=`, `aws_secret_access_key=`, hex strings ≥32 chars in known credential contexts)
- Files larger than 10 MB

Enable per-clone with `git config core.hooksPath .githooks`.

---

## Operational practices

These aren't code changes — they are process habits to maintain.

| Practice | Detail |
|---|---|
| **Credentials in chat** | Never paste real keys into a Claude session, shared terminal, or PR. The moment a secret is visible to anyone other than the system that needs it, it must be rotated. |
| **R2 token scope** | Always create R2 API tokens with **Object Read & Write** scoped to the `kinection` bucket. Never use account-wide tokens. |
| **Secret rotation cadence** | Rotate `COMPUTE_API_KEY` annually as a matter of hygiene, even without a known compromise. |
| **Local backups** | The user's DNA file should have a local backup. R2 has no copy and never will. |
| **Read access to D1** | `wrangler d1 execute --remote` is the only way to query the prod database. Don't grant additional collaborators access unless they need it. |

---

## Incident-response checklist

If a credential is suspected leaked:

1. **R2 token compromise** — Cloudflare dashboard → R2 → Manage R2 API Tokens → delete the token. Create a new one. Update `.env` on the daemon machine.
2. **`COMPUTE_API_KEY` compromise** — `openssl rand -hex 32` → `wrangler secret put COMPUTE_API_KEY --name kinection-api` → update `.env`.
3. **Cloudflare account compromise** — Reset password, revoke all API tokens, audit Worker deploy history (`wrangler deployments list`), rotate D1 (export → recreate database → reimport).
4. **DNA file leak** — There is no remediation. Acknowledge to the affected individual.
