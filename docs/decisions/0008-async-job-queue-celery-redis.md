# ADR-0008: Use Celery with Redis for Async Analysis Job Processing

* Status: Superseded by [ADR-0012](0012-cloudflare-workers-api.md) (Cloudflare Queues replaces Celery+Redis as the broker; Python compute workers remain but are triggered via HTTP by a Queue consumer Worker rather than consuming from Redis)
* Date: 2026-04-21

## Context and Problem Statement

Each user's DNA analysis takes several minutes of CPU-intensive computation (ASD across 394k SNPs × 17k individuals, PCA projection, haplogroup assignment). This cannot run synchronously within an HTTP request. A job queue is needed to accept uploads, enqueue analysis jobs, and serve results asynchronously. Which queueing architecture should be used?

## Decision Drivers

* Jobs are long-running (target ≤5 minutes per analysis) — must run outside the HTTP request cycle
* Workers must be scalable — auto-scale based on queue depth during traffic spikes
* Retry logic needed — if a worker crashes mid-analysis, the job should be retried
* Must integrate with Python (the analysis engine is Python)
* Results need to be stored and retrievable after the job completes
* Infrastructure should be straightforward to run locally for development

## Considered Options

* Celery + Redis (task queue + broker)
* Celery + RabbitMQ (task queue + broker)
* AWS SQS + Lambda (managed serverless queue)
* RQ (Redis Queue) — simpler Python-native alternative to Celery
* BullMQ (Node.js) — not applicable given Python stack

## Decision Outcome

Chosen option: **Celery with Redis as both the broker and result backend**, because Celery is the most mature Python task queue, Redis is simple to operate locally and in the cloud, and the combined Celery + Redis stack covers the requirements (retries, task status, result storage) without the overhead of managing a full message broker like RabbitMQ.

### Positive Consequences

* Native Python integration — analysis code runs directly inside a Celery worker without
  subprocess overhead or serialisation
* Redis serves as both the broker (task routing) and the result backend (job status/output)
  — single infrastructure component to operate
* Celery provides built-in retry logic with exponential backoff for failed jobs
* Worker count can be scaled horizontally — add more Celery workers behind the same Redis queue
* Excellent developer experience locally: `redis-server` + `celery worker` is all that is needed

### Negative Consequences

* Redis is in-memory — large result payloads (full pairwise distance tables) should not
  be stored in Redis; results should be written to S3/database and only a job-status
  flag stored in Redis
* Celery has significant configuration surface area and subtle edge cases
  (task acknowledgement, visibility timeout) that require care
* Redis persistence must be configured (AOF or RDB snapshots) to survive restarts —
  otherwise in-flight job state is lost on a Redis crash
* Long-running tasks (> Redis visibility timeout) need `acks_late=True` and
  a carefully configured `visibility_timeout`

## Pros and Cons of the Options

### Celery + Redis (chosen)

* Good, because native Python — analysis code runs directly in worker process
* Good, because single additional infrastructure component (Redis)
* Good, because mature retry, scheduling, and monitoring support (Flower)
* Bad, because Celery configuration complexity
* Bad, because Redis must be persisted and sized appropriately

### Celery + RabbitMQ

* Good, because RabbitMQ is purpose-built as a message broker — more reliable
  delivery guarantees than Redis pub/sub
* Good, because better support for complex routing patterns
* Bad, because operating RabbitMQ is more complex than Redis
* Bad, because two separate infrastructure components to manage (RabbitMQ + Redis or
  database for result storage)

### AWS SQS + Lambda

* Good, because fully managed — no broker or worker infrastructure to operate
* Good, because scales to zero when idle
* Bad, because Lambda has a 15-minute maximum execution time — tight for analysis
  jobs that may approach this limit with large files or high-SNP overlap
* Bad, because analysis code must be packaged as a Lambda deployment package
  including all Python dependencies (~hundreds of MB with NumPy/scikit-learn)
* Bad, because cold start latency is unacceptable for interactive user jobs
* Bad, because vendor lock-in

### RQ (Redis Queue)

* Good, because simpler than Celery — fewer configuration footguns
* Good, because Redis already chosen — no additional broker
* Bad, because less mature than Celery — fewer features (scheduling, canvas workflows,
  soft time limits are not supported)
* Bad, because smaller community — harder to find solutions to edge cases
* Bad, because monitoring tooling (RQ Dashboard) is less capable than Celery Flower
