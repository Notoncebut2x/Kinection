# Architecture Decision Records

This directory contains Architecture Decision Records (ADRs) in [MADR format](https://adr.github.io/madr/).

Each record documents a significant technical or architectural choice made in the project:
what the problem was, what options were considered, what was chosen, and why.

For a system-level overview, see [`../ARCHITECTURE.md`](../ARCHITECTURE.md).
For the science behind the pipeline, see [`../SCIENCE.md`](../SCIENCE.md).
For the security posture, see [`../SECURITY.md`](../SECURITY.md).

## Status legend

| Status    | Meaning |
|-----------|---------|
| Accepted  | Decision is in effect and implemented |
| Proposed  | Decision is planned but not yet implemented (Phase 2+ work) |
| Deprecated | No longer applies |
| Superseded | Replaced by a later ADR |

## Index

| # | Title | Status | Phase |
|---|-------|--------|-------|
| [0001](0001-reference-dataset-aadr-v62.md) | Use AADR v62.0 as reference dataset | Accepted | Phase 1 |
| [0002](0002-similarity-metric-allele-sharing-distance.md) | Use allele-sharing distance (ASD) as similarity metric | Accepted | Phase 1 |
| [0003](0003-pseudo-haploidisation-strategy.md) | Pseudo-haploidise modern individual for ASD computation | Accepted | Phase 1 |
| [0004](0004-analysis-language-python.md) | Use Python for the analysis engine | Accepted | Phase 1 |
| [0005](0005-palindromic-snp-exclusion.md) | Exclude palindromic SNPs during strand alignment | Accepted | Phase 1 |
| [0006](0006-haplogroup-reference-databases.md) | Use ISOGG Y-DNA tree and PhyloTree B17 for haplogroup assignment | Accepted | Phase 1 |
| [0007](0007-web-backend-fastapi.md) | Use FastAPI as the web backend framework | Superseded by 0012 | Phase 2 |
| [0008](0008-async-job-queue-celery-redis.md) | Use Celery with Redis for async job processing | Superseded by 0012 | Phase 2 |
| [0009](0009-database-postgresql.md) | Use PostgreSQL as the application database | Proposed | Phase 2 |
| [0010](0010-chunked-geno-processing.md) | Process GENO file in SNP-row chunks | Accepted | Phase 1 |
| [0011](0011-cloudflare-r2-geno-storage.md) | Store AADR GENO file in Cloudflare R2, access via range requests | Accepted | Phase 2 |
| [0012](0012-cloudflare-workers-api.md) | Use Cloudflare Workers as the web API layer | Accepted | Phase 2 |
| [0013](0013-admixture-nnls.md) | Admixture decomposition via constrained NNLS | Accepted | Phase 1 |
| [0014](0014-encoding-and-tmrca-calibration.md) | AADR PACKGENO encoding convention and Y-DNA TMRCA calibration | Accepted | Phase 1 |
