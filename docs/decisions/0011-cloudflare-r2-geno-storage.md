# ADR-0011: Store the AADR GENO File in Cloudflare R2 and Access via Range Requests

* Status: Accepted
* Date: 2026-05-17
* Supersedes: ADR-0001 (storage location only — dataset choice unchanged)

## Context and Problem Statement

The AADR v62.0 GENO file is 46 GB. For the web platform each user analysis job needs to read ~394,000 SNP rows from it. Storing and serving this file at scale raises three questions: where to host it, how to avoid downloading 46 GB per job, and how to minimise egress costs given the file will be read constantly in production.

## Decision Drivers

* Zero egress fees — the GENO file will be read for every analysis job; egress charges on AWS S3 would accumulate rapidly at scale
* HTTP range request support — needed to read individual SNP rows without downloading the full file
* S3-compatible API — the existing boto3 codebase can point at R2 with only endpoint/credential changes
* Consistent with the Cloudflare-first infrastructure direction (Workers, Queues, D1 all on Cloudflare)
* Simple pricing model — predictable costs for a dataset that is read-heavy but written only once

## Considered Options

* Cloudflare R2 with HTTP range requests (per-chunk reads)
* AWS S3 with HTTP range requests
* Self-hosted on a VPS alongside the compute workers
* Pre-compute a population frequency matrix and discard per-individual access

## Decision Outcome

Chosen option: **Cloudflare R2 with HTTP range requests**, because R2 has zero egress fees to the internet (unlike S3), is S3-compatible (no code changes beyond endpoint URL), and keeps the full infrastructure on one provider alongside the Workers API and Queues.

The GENO file is accessed via the `R2GenoFile` class (`scripts/utils/r2_geno.py`), which reads one byte-range per chunk covering all SNP rows in that chunk's index span. For 5,000 SNPs per chunk at 4,408 bytes/row, each range request fetches ~22 MB. The full analysis reads ~79 chunks = ~1.7 GB per job.

The `USE_R2=1` environment variable switches all three analysis scripts from local mmap to R2 range reads. Local mode continues to work unchanged for development.

### Positive Consequences

* Zero R2 egress fees — production read costs are limited to R2 operation charges (~$0.36 per million reads), not per-GB transfer
* One range request per 5,000-SNP chunk instead of one read_snp_row call per SNP — reduces HTTP overhead by ~5,000× vs a naïve row-by-row approach
* GENO file uploaded once, read by any number of parallel workers without coordination
* All user uploads and result files also in R2 — single storage service for the entire pipeline

### Negative Consequences

* Each chunk read is a network round-trip (vs local mmap which is near-zero latency); adds ~50–100ms per chunk
* Full analysis involves ~79 range requests totalling ~1.7 GB transfer per job — at R2 operation rates this is negligible but adds wall-clock time vs local
* Compute workers must be co-located in a region close to the R2 bucket to minimise chunk read latency (Cloudflare R2 is globally replicated, but Workers co-location still matters)
* If the AADR dataset is updated (new AADR version), the GENO file must be re-uploaded to R2

## Pros and Cons of the Options

### Cloudflare R2 (chosen)

* Good, because zero egress fees
* Good, because S3-compatible — boto3 with endpoint change only
* Good, because consistent with Cloudflare-first infrastructure
* Bad, because adds network latency per chunk vs local mmap
* Bad, because compute workers must be in a compatible region

### AWS S3

* Good, because proven at extreme scale, mature tooling
* Good, because range requests well-supported
* Bad, because egress pricing (~$0.09/GB) — 1.7 GB per analysis at scale is significant
* Bad, because second cloud provider alongside Cloudflare

### Self-hosted on VPS

* Good, because zero egress (reads are local)
* Good, because fastest possible chunk reads
* Bad, because single point of failure — VPS outage stops all analysis jobs
* Bad, because requires manual scaling when job volume increases
* Bad, because operator responsibility for disk, backups, and replication

### Pre-compute population frequency matrix

* Good, because reduces storage requirements from 46 GB to ~700 MB
* Good, because eliminates per-job GENO reads entirely
* Bad, because loses individual-level matching — results show only population means, not the closest ancient person specifically
* Bad, because pre-computation itself requires a full GENO read (one-time but expensive)
* Note: this optimisation is still planned as a fast-path complement to full individual-level analysis (Phase 3)
