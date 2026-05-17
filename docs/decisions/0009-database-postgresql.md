# ADR-0009: Use PostgreSQL as the Application Database

* Status: Proposed
* Date: 2026-04-21

## Context and Problem Statement

The web platform needs a persistent database for user accounts, upload records, job status, and structured analysis results (haplogroup assignments, top match summaries). Raw genotype files and full distance tables will be stored separately in object storage (S3). Which relational database should be used for application state?

## Decision Drivers

* ACID guarantees — genetic data and consent records must not be lost or corrupted
* JSON column support — analysis results are naturally semi-structured (haplogroup evidence,
  top match lists) and benefit from a JSON-capable relational database
* Strong ecosystem for Python ORMs (SQLAlchemy, psycopg)
* Handles moderate scale — not a high-QPS system; hundreds of concurrent users at most
* Widely supported on cloud platforms (AWS RDS, Supabase, Neon, Fly Postgres)
* GDPR compliance features: row-level deletion, audit logging support

## Considered Options

* PostgreSQL
* MySQL / MariaDB
* SQLite (development only)
* MongoDB (document store)
* DynamoDB (managed NoSQL)

## Decision Outcome

Chosen option: **PostgreSQL**, because it provides ACID guarantees, native JSONB columns for semi-structured result data, excellent Python library support (SQLAlchemy + psycopg3), and is available as a managed service on all major cloud platforms. JSONB in PostgreSQL enables querying inside result documents (e.g., filtering users by Y haplogroup) without sacrificing relational integrity for user and consent records.

### Positive Consequences

* JSONB columns store haplogroup results and top match arrays natively, queryable with
  standard SQL operators — avoids a separate document store
* Full ACID compliance — consent timestamps, deletion requests, and job state transitions
  are atomic and durable
* Row-level security (RLS) enables future multi-tenant isolation at the database layer
* SQLAlchemy ORM provides a clean abstraction; Alembic handles schema migrations
* Managed PostgreSQL (RDS, Supabase, Neon) removes operational burden of backups and
  failover

### Negative Consequences

* More operationally complex than SQLite for local development — requires a running
  Postgres instance (mitigated with Docker Compose)
* For very high read loads, Postgres requires a read replica or connection pooler
  (PgBouncer) — not a concern at early scale but worth planning
* Schema migrations require care when the dataset grows and result schemas evolve

## Pros and Cons of the Options

### PostgreSQL (chosen)

* Good, because ACID, JSONB, mature Python support, managed cloud options
* Good, because widely understood — easy to hire for, extensive documentation
* Good, because GDPR-friendly (row deletion, pgaudit extension)
* Bad, because more setup than SQLite for local development
* Bad, because connection pooling required at scale

### MySQL / MariaDB

* Good, because widely available and understood
* Bad, because JSON support is less capable than PostgreSQL JSONB
* Bad, because stricter ANSI SQL compliance issues with some query patterns
* Bad, because weaker full-text search if ever needed for population/culture filtering

### SQLite

* Good, because zero setup — single file, ships with Python
* Good, because ideal for local development and testing
* Bad, because not suitable for production with concurrent writes (write locking)
* Bad, because no native JSON query operators
* Decision: use SQLite for local development only; deploy with PostgreSQL

### MongoDB

* Good, because schema-flexible for evolving result structures
* Good, because native document storage eliminates JSONB workaround
* Bad, because no joins — user/upload/job/result relationships require application-level
  joins, which is error-prone
* Bad, because weaker ACID guarantees in multi-document transactions
* Bad, because GDPR right-to-erasure across nested documents is operationally complex

### DynamoDB

* Good, because fully managed, scales to any load automatically
* Bad, because query patterns must be known at schema design time (partition key design)
* Bad, because complex relational queries (user → uploads → jobs → results) are expensive
* Bad, because vendor lock-in to AWS
