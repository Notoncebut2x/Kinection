# ADR-0007: Use FastAPI as the Web Backend Framework

* Status: Superseded by [ADR-0012](0012-cloudflare-workers-api.md)
* Date: 2026-04-21

## Context and Problem Statement

The web platform (Phase 2) needs a backend API that handles file uploads, triggers analysis jobs, serves results to the frontend, and manages user accounts. The backend must integrate naturally with the existing Python analysis codebase. Which web framework should serve as the backend?

## Decision Drivers

* Native Python — the analysis engine is Python and sharing code between
  the API layer and the analysis pipeline reduces duplication
* Async support — file upload handling and job status polling benefit from
  async I/O; analysis jobs themselves are async (Celery workers)
* Performance adequate for the expected load (not a high-QPS API)
* Built-in request validation and OpenAPI documentation generation
* Easy to deploy on standard cloud infrastructure (AWS Lambda, ECS, Fly.io)

## Considered Options

* FastAPI (Python async framework with Pydantic validation)
* Django REST Framework (DRF)
* Flask
* Node.js / Express (separate language from analysis engine)

## Decision Outcome

Chosen option: **FastAPI**, because it is the modern Python async API framework with first-class Pydantic validation, automatic OpenAPI/Swagger docs, and native async support for upload handling. It shares the Python runtime with the analysis engine, allowing direct import of parsing utilities and result models without serialisation overhead.

### Positive Consequences

* Pydantic models enforce and document the API contract — reduces a class of input bugs
* Automatic OpenAPI spec generation at `/docs` — aids frontend development and testing
* Async route handlers fit naturally with async file upload and S3 operations
* Shares the Python venv with the analysis engine — no separate runtime to manage
* FastAPI + Uvicorn is lightweight to deploy on containerised infrastructure

### Negative Consequences

* FastAPI does not include an ORM or admin panel out of the box — SQLAlchemy
  and a separate admin tool (e.g., SQLAdmin) must be added
* Django has a more mature ecosystem for authentication, admin, and permissions —
  these will need to be built or sourced from third-party FastAPI packages
* Smaller community than Django for finding pre-built solutions to common patterns

## Pros and Cons of the Options

### FastAPI (chosen)

* Good, because native Python async — shares codebase with analysis engine
* Good, because Pydantic validation and automatic OpenAPI docs
* Good, because lightweight and fast to deploy
* Bad, because no built-in ORM, admin, or auth — requires additional packages
* Bad, because smaller community than Django

### Django REST Framework

* Good, because mature ecosystem (auth, admin, ORM, migrations all built in)
* Good, because battle-tested at scale
* Bad, because synchronous by default — async support is bolted on and less ergonomic
* Bad, because heavier than needed for an API-only backend
* Bad, because Django's ORM is less flexible than SQLAlchemy for complex queries

### Flask

* Good, because minimal and flexible
* Bad, because no async support without significant additional setup
* Bad, because no built-in validation — requires adding marshmallow or Pydantic manually
* Bad, because less ergonomic than FastAPI for building typed APIs

### Node.js / Express

* Good, because excellent async I/O performance
* Good, because large ecosystem for web APIs
* Bad, because separate language from the analysis engine — no shared code
* Bad, because the analysis results must be serialised and passed between runtimes
* Bad, because introduces a second language to maintain
