# Kinection — Technical Architecture

*Last updated: 2026-05-17*

This document describes how Kinection is built end-to-end: the components, what runs where, how data flows, and why the system is shaped this way.

For the rationale behind individual design choices, see [`docs/decisions/`](decisions/).
For an introduction to the *science* the system performs, see [`SCIENCE.md`](SCIENCE.md).

---

## 1. System Overview

Kinection compares a user's personal AncestryDNA file against the Allen Ancient DNA Resource (AADR) — a curated dataset of ~19,000 ancient human genomes — and produces a report describing how the user's DNA relates to those ancient populations.

```
                                ┌─────────────────────────┐
                                │   Cloudflare (cloud)    │
                                │                         │
   ┌────────────┐    POST /jobs │   ┌─────────────────┐   │
   │ Web client ├──────────────►│   │ Workers API     │   │
   │ (browser)  │   GET /jobs   │   │ (TypeScript)    │   │
   └────────────┘◄──────────────┤   └────┬────────────┘   │
                                │        │                │
                                │   ┌────▼────┐  ┌─────┐  │
                                │   │   D1    │  │ R2  │  │
                                │   │ (jobs)  │  │ AADR│  │
                                │   └────┬────┘  │ outs│  │
                                │        │       └──▲──┘  │
                                └────────┼──────────┼─────┘
                                         │ poll     │ read AADR / write results
   ┌─────────────────────────────────────┼──────────┼──┐
   │ User's machine                      │          │  │
   │                                     ▼          │  │
   │   ┌──────────────────┐  ┌───────────────────┐  │  │
   │   │ DNA file (local) ├─►│ Python daemon +   ├──┘  │
   │   │ AncestryDNA.txt  │  │ analysis pipeline │     │
   │   └──────────────────┘  └───────────────────┘     │
   │                                                    │
   └────────────────────────────────────────────────────┘
```

The split between **cloud** (Cloudflare) and **local** (the user's machine) is deliberate: personal DNA never leaves the user's machine. Only derived analysis outputs are stored in the cloud, and only under an opaque job ID.

---

## 2. Components

### 2.1 Worker API — Cloudflare Workers (TypeScript)

**Path:** `workers/api/src/index.ts`
**Deployment:** `wrangler deploy` → `https://kinection-api.<subdomain>.workers.dev`

A stateless HTTP API at the edge. It is the only public-facing component.

| Method | Path | Auth | Caller | Purpose |
|---|---|---|---|---|
| `GET` | `/dataset/version` | — | Web client | Current AADR version (KV-cached) |
| `POST` | `/jobs` | — | Web client | Create a new analysis job |
| `GET` | `/jobs?status=queued` | Bearer | Daemon | Poll for pending jobs |
| `GET` | `/jobs/:id` | — | Web client | Read job status |
| `PATCH` | `/jobs/:id/status` | Bearer | Daemon | Report progress / completion |
| `GET` | `/jobs/:id/results/:filename` | — | Web client | Download a result file |

Bindings (declared in `wrangler.toml`):

- `env.R2` — `kinection` R2 bucket
- `env.DB` — `kinection` D1 database
- `env.MARKER_CACHE` — KV namespace (reserved; not yet used)
- `env.COMPUTE_API_KEY` — shared secret for daemon authentication (set via `wrangler secret put`)

### 2.2 D1 — Job state (SQLite at the edge)

**Schema:** `workers/api/src/schema.sql`

```sql
CREATE TABLE jobs (
  id         TEXT PRIMARY KEY,                    -- crypto.randomUUID()
  status     TEXT NOT NULL DEFAULT 'queued',     -- queued | processing | complete | failed
  created_at INTEGER NOT NULL,                    -- unix ms
  updated_at INTEGER NOT NULL,
  error      TEXT
);
CREATE INDEX idx_jobs_status ON jobs(status);
```

D1 only stores job *state*, never genetic content. Job IDs are unguessable UUIDv4 strings.

### 2.3 R2 — Object storage

**Bucket:** `kinection`
**Key layout:**

```
dataset/v66/v66.1240K.aadr.PUB.geno        ← AADR reference (~7 GB)
dataset/v66/v66.1240K.aadr.PUB.ind         ← individual list
dataset/v66/v66.1240K.aadr.PUB.snp         ← SNP positions
dataset/v66/v66.1240K.aadr.PUB.anno        ← annotation
dataset/current_version.json               ← which AADR version is current

outputs/<job-id>/snp_overlap.tsv           ← step 1 outputs
outputs/<job-id>/step1_summary.json
outputs/<job-id>/haplogroup_report.md      ← step 2 outputs
outputs/<job-id>/ydna_haplogroup.json
outputs/<job-id>/mtdna_haplogroup.json
outputs/<job-id>/top_matches_report.md     ← step 3 outputs
outputs/<job-id>/pairwise_distances.tsv
outputs/<job-id>/pca_coordinates.tsv
```

R2 was chosen for two reasons: zero egress fees (the AADR `.geno` file is read repeatedly via HTTP range requests) and S3-compatible boto3 access from Python. See [ADR-0011](decisions/0011-cloudflare-r2-geno-storage.md).

**What R2 never stores:** the user's raw DNA file. That is read from the user's local disk at analysis time and stays there.

### 2.4 Compute daemon — Python (`scripts/daemon.py`)

A long-running process on the user's machine. Polls the Worker API every `POLL_INTERVAL` seconds, claims queued jobs, runs the three-step pipeline, and reports back.

```
loop {
  jobs ← GET /jobs?status=queued
  for job in jobs:
    PATCH /jobs/{job.id}/status {status: processing}
    try:
      run step1, step2, step3 (env: JOB_ID, USE_R2=1)
      PATCH /jobs/{job.id}/status {status: complete}
    except:
      PATCH /jobs/{job.id}/status {status: failed, error: ...}
  sleep POLL_INTERVAL
}
```

The daemon is the *only* component with access to the user's DNA file.

### 2.5 Analysis pipeline — Python (`scripts/step{1,2,3}_*.py`)

Three sequential steps. Each reads from local disk (modern DNA) and R2 (AADR reference), and writes outputs to R2 under `outputs/<JOB_ID>/`.

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 1 | `step1_parse_harmonise.py` | AncestryDNA + AADR refs | SNP overlap, encoded modern genotypes |
| 2 | `step2_haplogroup.py` | step 1 outputs + Y/mtDNA marker DBs | Y-DNA + mtDNA haplogroup, haplogroup matches |
| 3 | `step3_similarity_pca.py` | step 1 outputs + AADR GENO | Pairwise ASD, population ranking, PCA |

All three steps share `scripts/utils/`:
- `parsers.py` — AncestryDNA, EIGENSTRAT `.ind`, AADR `.anno`, palindromic SNP filter
- `r2_client.py` — boto3-wrapped R2 access; defines the `dataset/v<N>/...` key paths
- `r2_geno.py` — `R2GenoFile` class: reads EIGENSTRAT PACKGENO rows via HTTP range requests, no full download

### 2.6 AADR updater — Python (`scripts/update_aadr.py`)

Standalone tool. Queries the Harvard Dataverse API for the latest AADR release, streams new files directly from Dataverse into R2 (no local disk), writes a manifest, and patches `r2_client.py` constants. Idempotent — safe to run regularly.

```bash
python scripts/update_aadr.py --check   # report status, no upload
python scripts/update_aadr.py           # upload if newer
python scripts/update_aadr.py --force   # re-upload regardless
```

---

## 3. Data Flow — End-to-End

### 3.1 The happy path (one analysis run)

```
1. Web client → POST /jobs
   Worker → INSERT INTO jobs (id=uuid, status='queued')
   Worker → 201 {id, status: 'queued', created_at}

2. Web client polls GET /jobs/{id} every few seconds

3. Daemon (local) polls GET /jobs?status=queued (with Bearer auth)
   Worker → SELECT * FROM jobs WHERE status='queued' LIMIT 10
   Daemon receives [{id, status, ...}]

4. Daemon → PATCH /jobs/{id}/status {status: 'processing'}
   Worker → UPDATE jobs SET status='processing'

5. Daemon runs:
   - step1 reads local AncestryDNA file
   - step1 reads AADR .ind/.snp/.anno from R2 (small)
   - step1 streams AADR .geno from R2 via HTTP range requests
   - step1 writes outputs/<id>/snp_overlap.tsv etc. → R2
   - step2 reads step1 outputs from R2, writes step2 outputs → R2
   - step3 reads step1 outputs + AADR .geno (range reads) → writes step3 outputs → R2

6. Daemon → PATCH /jobs/{id}/status {status: 'complete'}
   Worker → UPDATE jobs SET status='complete'

7. Web client GET /jobs/{id} → sees status='complete'
   Web client GET /jobs/{id}/results/haplogroup_report.md
   Worker reads R2 outputs/<id>/haplogroup_report.md → streams to client
```

### 3.2 What never crosses the cloud boundary

| Asset | Location | Notes |
|---|---|---|
| User's AncestryDNA file | Local disk only | Read by step1, never uploaded |
| Raw genotype dosage array | Local memory during step1 | Encoded version written *only* to R2 under opaque job ID |
| Worker → daemon callback | None | Daemon polls. Worker never reaches out. |

This is enforced in code (`r2_client.py` has no `upload_key()` helper for raw DNA), in the Worker (no upload route exists), and in `.gitignore` (`data/input_data/` and `output/` are never tracked).

---

## 4. Why the system is shaped this way

**Why a daemon instead of running compute in Workers?**
The pipeline uses NumPy/SciPy and reads ~7 GB of binary data via random-access HTTP range requests. Workers' Python runtime (Pyodide) can't run NumPy and has CPU/memory limits unsuitable for this workload. More importantly, the daemon must run *where the DNA file is*, which is on the user's machine.

**Why poll instead of a queue?**
A Cloudflare Queue with a consumer Worker would need to call back into the daemon — but the daemon sits behind NAT on a laptop, with no inbound address. Poll-based pull is the cleanest match.

**Why R2 over S3?**
Zero egress fees. The `.geno` file is read repeatedly via byte-range requests during PCA and ASD computation; on S3 that would be expensive. R2 is also S3-compatible so boto3 works unchanged.

**Why pseudo-haploidisation in ASD?**
The AADR data is pseudo-haploid (one allele per ancient individual due to ancient DNA degradation), while AncestryDNA is diploid. Comparing them directly would introduce a 0.5 bias at heterozygous sites. See [ADR-0003](decisions/0003-pseudo-haploidisation.md).

**Why exclude palindromic SNPs?**
A/T and C/G SNPs are ambiguous under strand flips, which are common between SNP arrays. Excluding them prevents subtle systematic errors. See [ADR-0005](decisions/0005-palindromic-snp-exclusion.md).

---

## 5. Security boundaries

| Boundary | Mechanism |
|---|---|
| Public web → Worker | HTTPS, CORS-permissive while in dev; will be locked down to frontend origin |
| Daemon → Worker (privileged routes) | `Authorization: Bearer <COMPUTE_API_KEY>` |
| Daemon → R2 | boto3 with R2 API token credentials (Object Read & Write on `kinection`) |
| Worker → R2/D1/KV | Cloudflare bindings (no over-the-wire credentials) |
| Job IDs | UUIDv4, unguessable; act as a capability for result fetch |
| User's DNA file | Never leaves local disk; not represented in any cloud key |

Outstanding security work is tracked in [`SECURITY.md`](SECURITY.md).

---

## 6. Deployment

A first-time setup looks like this:

```bash
# Cloudflare resources (one-time)
wrangler login
wrangler d1 create kinection
wrangler kv namespace create MARKER_CACHE
wrangler queues create analysis-jobs           # currently unused, will deprecate
wrangler r2 bucket create kinection

# Apply schema (remote)
cd workers/api && wrangler d1 execute kinection --file=src/schema.sql --remote

# Set the shared secret
wrangler secret put COMPUTE_API_KEY --name kinection-api

# Deploy the Worker
wrangler deploy --name kinection-api

# Local setup
cd ../..
python3 -m venv venv && source venv/bin/activate
pip install boto3 requests python-dotenv

# Fill out .env (see .env.example)
cp .env.example .env

# Upload the AADR dataset (one-time, ~7 GB)
python scripts/update_aadr.py

# Run the daemon
python scripts/daemon.py
```

See [README.md](../README.md) for the user-facing version.

---

## 7. Outstanding architecture work

These are known to be incomplete or under-built:

- **Result URL signing** — `GET /jobs/:id/results/:filename` currently has no rate limit and serves anyone with the job UUID. For a public deployment, this should issue short-lived R2 signed URLs instead of streaming through the Worker.
- **Daemon parallelism** — the daemon processes one job at a time. Fine for personal use; needs revisiting for multi-user beta.
- **PII retention** — there is no automatic deletion of job rows in D1 or output files in R2. A retention/cleanup ADR is needed before any public launch.
- **Extending KV use** — currently only the version manifest is KV-cached. Y-DNA and mtDNA marker JSONs could also be cached if/when the daemon stops bundling them locally.
