# ADR-0012: Use Cloudflare Workers as the Web API Layer

* Status: Accepted
* Date: 2026-05-17
* Supersedes: ADR-0007 (FastAPI web backend)

## Context and Problem Statement

The web platform needs an API layer that handles file uploads, triggers analysis jobs, tracks job status, and serves results to the frontend. ADR-0007 proposed FastAPI (Python) for this role. With the shift to Cloudflare R2 for storage and Cloudflare Queues for job dispatch, the API layer should be reconsidered alongside those choices.

## Decision Drivers

* Low cold-start latency for upload and status endpoints (users waiting on responses)
* Native integration with R2 (stream uploads directly to R2 without buffering through a server)
* Native integration with Cloudflare Queues (publish jobs from the same runtime that handles uploads)
* Eliminate the need to run a persistent API server — reduce operational overhead
* The API logic is lightweight (upload dispatch, status query, result proxy) — no heavy computation happens in the API layer

## Considered Options

* Cloudflare Workers (TypeScript, edge compute)
* FastAPI on a persistent server (Python, as originally planned in ADR-0007)
* FastAPI on a serverless platform (AWS Lambda + API Gateway)

## Decision Outcome

Chosen option: **Cloudflare Workers**, because the API responsibilities are purely I/O — stream a file to R2, write a row to D1, publish to a Queue, read a row from D1, proxy a file from R2. Workers handle all of these natively with zero server management and sub-millisecond cold starts. The heavy Python computation is intentionally kept out of Workers (in external Python compute workers) since Workers' 128 MB memory limit makes NumPy-based analysis impossible there.

Services used:
| Cloudflare service | Role |
|-------------------|------|
| Workers | API routes (upload, job status, result delivery) |
| R2 | User DNA uploads, analysis result files, AADR dataset |
| Queues | Analysis job dispatch (producer in Worker, consumed by Python compute service) |
| D1 (SQLite) | Job status tracking (id, status, upload_key, result_prefix, timestamps) |
| KV | Haplogroup marker databases (globally cached JSON — fast reads from any Worker) |

### Positive Consequences

* Upload handler streams directly from the HTTP request body into R2 — no buffering on a server, no disk I/O
* D1 job status writes and reads are sub-millisecond from Workers
* Queue publish is synchronous from Workers — no risk of a job being lost between API response and job dispatch
* Zero infrastructure to manage for the API layer — Workers scale automatically
* PATCH /jobs/:id/status endpoint lets the Python compute worker report completion with a simple authenticated HTTP call

### Negative Consequences

* Workers use TypeScript — separate language from the Python analysis engine; two codebases to maintain
* Workers have 128 MB memory limit — the analysis pipeline cannot run inside a Worker (by design: Python compute workers handle that)
* D1 is SQLite — limited to ~100k writes/day on the free tier; the paid tier handles production load
* Cloudflare Queues consumer is Worker-only — Python workers cannot directly consume from the Queue; the Worker consumer calls the Python compute service via HTTP instead

## Pros and Cons of the Options

### Cloudflare Workers (chosen)

* Good, because zero server management, automatic scaling
* Good, because native R2, Queue, and D1 bindings — no credentials or SDK needed
* Good, because upload streaming directly into R2 — no intermediate buffering
* Good, because consistent with Cloudflare-first infrastructure (ADR-0011)
* Bad, because TypeScript adds a second language to the project
* Bad, because Workers cannot run the Python analysis pipeline

### FastAPI on a persistent server (original ADR-0007 plan)

* Good, because single Python codebase shared with analysis engine
* Good, because no language boundary between API and analysis code
* Bad, because persistent server to manage, scale, and keep available
* Bad, because file uploads buffer through the server before reaching storage
* Bad, because separate cloud provider from Cloudflare (mixed infrastructure)

### FastAPI on AWS Lambda

* Good, because serverless — no persistent server
* Good, because Python runtime shared with analysis engine
* Bad, because Lambda cold starts can be 1–3 seconds for Python + NumPy packages
* Bad, because AWS egress costs for large file uploads and result downloads
* Bad, because second cloud provider alongside Cloudflare
