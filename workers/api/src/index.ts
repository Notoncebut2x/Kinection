import { AwsClient } from "aws4fetch";

export interface Env {
  R2: R2Bucket;
  DB: D1Database;
  MARKER_CACHE: KVNamespace;
  COMPUTE_API_KEY: string;

  // R2 S3-compatible credentials, used to mint presigned PUT URLs so the
  // browser uploads directly to R2 (the raw DNA file never traverses the
  // Worker). See ADR 0015 / Step 5.1.1.
  R2_ACCOUNT_ID: string;
  R2_ACCESS_KEY_ID: string;
  R2_SECRET_ACCESS_KEY: string;
  R2_BUCKET_NAME: string;
}

interface JobRow {
  id: string;
  user_id: string | null;
  status: string;
  label: string | null;
  aadr_version: string | null;
  created_at: number;
  updated_at: number;
  started_at: number | null;
  completed_at: number | null;
  error: string | null;
}

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, POST, PATCH, DELETE, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type, Authorization",
};

const UPLOAD_URL_TTL_SECONDS = 15 * 60; // 15 min — ADR 0015 / Step 5.1.1
const MAX_UPLOAD_BYTES = 50 * 1024 * 1024; // 50 MB cap for consumer DNA files

function json(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { ...CORS_HEADERS, "Content-Type": "application/json" },
  });
}

function err(message: string, status: number): Response {
  return json({ error: message }, status);
}

function isAuthorized(request: Request, env: Env): boolean {
  return request.headers.get("Authorization") === `Bearer ${env.COMPUTE_API_KEY}`;
}

function uploadKeyFor(jobId: string): string {
  return `uploads/${jobId}/raw.txt`;
}

function nowSeconds(): number {
  return Math.floor(Date.now() / 1000);
}

async function logAudit(
  env: Env,
  event: string,
  detail: {
    actor?: string;
    job_id?: string;
    r2_key?: string;
    extra?: Record<string, unknown>;
  }
): Promise<void> {
  // Best-effort — audit log writes never block the user response.
  try {
    await env.DB.prepare(
      `INSERT INTO audit_log (id, at, actor, event, job_id, r2_key, detail)
       VALUES (?, ?, ?, ?, ?, ?, ?)`
    )
      .bind(
        crypto.randomUUID(),
        nowSeconds(),
        detail.actor ?? null,
        event,
        detail.job_id ?? null,
        detail.r2_key ?? null,
        detail.extra ? JSON.stringify(detail.extra) : null
      )
      .run();
  } catch (e) {
    console.error("audit_log write failed", e);
  }
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: CORS_HEADERS });
    }

    const url = new URL(request.url);
    const { pathname } = url;

    try {
      // ── Public endpoints ──────────────────────────────────────────────
      if (request.method === "GET" && pathname === "/dataset/version") {
        return handleDatasetVersion(env);
      }

      // POST /uploads/url — mint a presigned PUT URL for raw modern DNA
      if (request.method === "POST" && pathname === "/uploads/url") {
        return handleCreateUploadUrl(request, env);
      }

      if (request.method === "POST" && pathname === "/jobs") {
        return handleCreateJob(env);
      }

      // POST /jobs/:id/upload-complete — browser signals the presigned PUT
      // finished; flip 'uploading' -> 'queued' so the daemon picks it up.
      // Public: the job_id UUID is the capability token (see handleUserDelete).
      const uploadCompleteMatch = pathname.match(
        /^\/jobs\/([^/]+)\/upload-complete$/
      );
      if (request.method === "POST" && uploadCompleteMatch) {
        return handleUploadComplete(uploadCompleteMatch[1], env);
      }

      // ── Daemon endpoints (require auth) ───────────────────────────────
      if (request.method === "GET" && pathname === "/jobs") {
        return handleListJobs(request, url, env);
      }

      const jobMatch = pathname.match(/^\/jobs\/([^/]+)$/);
      if (request.method === "GET" && jobMatch) {
        return handleGetJob(jobMatch[1], env);
      }

      const statusMatch = pathname.match(/^\/jobs\/([^/]+)\/status$/);
      if (request.method === "PATCH" && statusMatch) {
        return handleUpdateStatus(request, statusMatch[1], env);
      }

      const receiptMatch = pathname.match(/^\/jobs\/([^/]+)\/deletion_receipt$/);
      if (request.method === "POST" && receiptMatch) {
        return handleDeletionReceipt(request, receiptMatch[1], env);
      }

      // DELETE /jobs/:id/upload — user-triggered immediate deletion of raw upload
      const userDeleteMatch = pathname.match(/^\/jobs\/([^/]+)\/upload$/);
      if (request.method === "DELETE" && userDeleteMatch) {
        return handleUserDelete(userDeleteMatch[1], env);
      }

      const resultMatch = pathname.match(/^\/jobs\/([^/]+)\/results\/([^/]+)$/);
      if (request.method === "GET" && resultMatch) {
        return handleGetResult(resultMatch[1], resultMatch[2], env);
      }

      return err("Not found", 404);
    } catch (e) {
      console.error(e);
      return err("Internal server error", 500);
    }
  },
};

// ---------------------------------------------------------------------------
// Dataset version manifest (KV-cached)
// ---------------------------------------------------------------------------
const DATASET_VERSION_KEY = "dataset/current_version.json";
const DATASET_VERSION_KV_KEY = "aadr_version_manifest";
const DATASET_VERSION_TTL_SECONDS = 300;

async function handleDatasetVersion(env: Env): Promise<Response> {
  const cached = await env.MARKER_CACHE.get(DATASET_VERSION_KV_KEY, "text");
  if (cached) {
    return new Response(cached, {
      headers: {
        ...CORS_HEADERS,
        "Content-Type": "application/json",
        "X-Cache": "HIT",
      },
    });
  }

  const object = await env.R2.get(DATASET_VERSION_KEY);
  if (!object) {
    return err("Dataset version manifest not found in R2", 404);
  }

  const text = await object.text();
  await env.MARKER_CACHE.put(DATASET_VERSION_KV_KEY, text, {
    expirationTtl: DATASET_VERSION_TTL_SECONDS,
  });

  return new Response(text, {
    headers: {
      ...CORS_HEADERS,
      "Content-Type": "application/json",
      "X-Cache": "MISS",
    },
  });
}

// ---------------------------------------------------------------------------
// Upload flow: create job + mint presigned PUT URL
// ---------------------------------------------------------------------------

async function handleCreateUploadUrl(
  request: Request,
  env: Env
): Promise<Response> {
  let body: { label?: string; format?: string; sha256?: string; size_bytes?: number };
  try {
    body = await request.json();
  } catch {
    return err("Invalid JSON body", 400);
  }

  if (body.size_bytes !== undefined && body.size_bytes > MAX_UPLOAD_BYTES) {
    return err(
      `Upload too large: ${body.size_bytes} bytes (max ${MAX_UPLOAD_BYTES})`,
      413
    );
  }

  const jobId = crypto.randomUUID();
  const uploadId = crypto.randomUUID();
  const now = nowSeconds();
  const r2Key = uploadKeyFor(jobId);

  // Job is 'uploading' until the daemon picks it up. Browser must call
  // POST /jobs/:id/upload-complete (TODO) after the PUT succeeds — for v1
  // we let the daemon discover the file via a `queued` status transition.
  await env.DB.batch([
    env.DB.prepare(
      `INSERT INTO jobs (id, status, label, created_at, updated_at)
       VALUES (?, 'uploading', ?, ?, ?)`
    ).bind(jobId, body.label ?? null, now, now),
    env.DB.prepare(
      `INSERT INTO uploads (id, job_id, r2_key, size_bytes, sha256,
                            content_type, format, uploaded_at)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)`
    ).bind(
      uploadId,
      jobId,
      r2Key,
      body.size_bytes ?? 0,
      body.sha256 ?? "",
      "text/plain",
      body.format ?? "unknown",
      now
    ),
  ]);

  const uploadUrl = await mintPresignedPut(env, r2Key, UPLOAD_URL_TTL_SECONDS);
  const expiresAt = now + UPLOAD_URL_TTL_SECONDS;

  await logAudit(env, "upload_url_minted", {
    actor: "system",
    job_id: jobId,
    r2_key: r2Key,
    extra: { ttl_seconds: UPLOAD_URL_TTL_SECONDS },
  });

  return json(
    {
      job_id: jobId,
      upload_id: uploadId,
      upload_key: r2Key,
      upload_url: uploadUrl,
      expires_at: expiresAt,
    },
    201
  );
}

async function mintPresignedPut(
  env: Env,
  key: string,
  ttlSeconds: number
): Promise<string> {
  const aws = new AwsClient({
    accessKeyId: env.R2_ACCESS_KEY_ID,
    secretAccessKey: env.R2_SECRET_ACCESS_KEY,
    service: "s3",
    region: "auto",
  });

  const endpoint =
    `https://${env.R2_ACCOUNT_ID}.r2.cloudflarestorage.com/` +
    `${env.R2_BUCKET_NAME}/${encodeURI(key)}`;

  // aws4fetch presign: pass the URL with X-Amz-Expires; the signer fills
  // in the rest of the canonical query params.
  const url = new URL(endpoint);
  url.searchParams.set("X-Amz-Expires", String(ttlSeconds));

  const signed = await aws.sign(
    new Request(url.toString(), { method: "PUT" }),
    { aws: { signQuery: true } }
  );
  return signed.url;
}

// ---------------------------------------------------------------------------
// Original jobs endpoints
// ---------------------------------------------------------------------------

async function handleCreateJob(env: Env): Promise<Response> {
  const id = crypto.randomUUID();
  const now = nowSeconds();

  await env.DB.prepare(
    `INSERT INTO jobs (id, status, created_at, updated_at) VALUES (?, 'queued', ?, ?)`
  )
    .bind(id, now, now)
    .run();

  return json({ id, status: "queued", created_at: now }, 201);
}

async function handleUploadComplete(
  jobId: string,
  env: Env
): Promise<Response> {
  // Confirm the job exists and is still awaiting its upload. Anything other
  // than 'uploading' is a no-op success (idempotent — safe under retries and
  // double-clicks; also tolerates the daemon having already advanced it).
  const job = await env.DB.prepare(
    "SELECT id, status FROM jobs WHERE id = ?"
  )
    .bind(jobId)
    .first<{ id: string; status: string }>();
  if (!job) return err("Job not found", 404);
  if (job.status !== "uploading") {
    return json({ ok: true, status: job.status, note: "already advanced" });
  }

  // Defence-in-depth: verify the raw object actually landed in R2 before we
  // queue the job. A presigned PUT that never completed leaves nothing here,
  // and we must not hand the daemon a job with no file to read.
  const upload = await env.DB.prepare(
    "SELECT r2_key FROM uploads WHERE job_id = ?"
  )
    .bind(jobId)
    .first<{ r2_key: string }>();
  if (!upload) return err("Upload record not found", 404);

  const head = await env.R2.head(upload.r2_key);
  if (head === null) {
    return err("Upload not found in storage — did the PUT complete?", 409);
  }

  const now = nowSeconds();
  const result = await env.DB.prepare(
    `UPDATE jobs SET status = 'queued', updated_at = ?
     WHERE id = ? AND status = 'uploading'`
  )
    .bind(now, jobId)
    .run();
  // Lost the race to another writer (e.g. the reaper aborting a stale job).
  if (result.meta.changes === 0) {
    const current = await env.DB.prepare("SELECT status FROM jobs WHERE id = ?")
      .bind(jobId)
      .first<{ status: string }>();
    return json({ ok: true, status: current?.status ?? "unknown" });
  }

  await logAudit(env, "upload_completed", {
    actor: "user",
    job_id: jobId,
    r2_key: upload.r2_key,
    extra: { size_bytes: head.size },
  });

  return json({ ok: true, status: "queued" });
}

async function handleListJobs(
  request: Request,
  url: URL,
  env: Env
): Promise<Response> {
  if (!isAuthorized(request, env)) return err("Unauthorized", 401);

  const status = url.searchParams.get("status") ?? "queued";
  const result = await env.DB.prepare(
    `SELECT * FROM jobs WHERE status = ? ORDER BY created_at ASC LIMIT 10`
  )
    .bind(status)
    .all<JobRow>();

  return json(result.results);
}

async function handleGetJob(id: string, env: Env): Promise<Response> {
  const result = await env.DB.prepare("SELECT * FROM jobs WHERE id = ?")
    .bind(id)
    .first<JobRow>();
  if (!result) return err("Job not found", 404);
  return json(result);
}

async function handleUpdateStatus(
  request: Request,
  id: string,
  env: Env
): Promise<Response> {
  if (!isAuthorized(request, env)) return err("Unauthorized", 401);

  const body = await request.json<{ status: string; error?: string }>();
  if (!body.status) return err("Missing status field", 400);

  const now = nowSeconds();
  // Set started_at on first transition to 'processing'; completed_at on
  // terminal states. Both columns are nullable so re-writes are safe.
  const startedAtClause =
    body.status === "processing" ? ", started_at = COALESCE(started_at, ?)" : "";
  const completedAtClause =
    body.status === "complete" || body.status === "failed" || body.status === "aborted"
      ? ", completed_at = COALESCE(completed_at, ?)"
      : "";

  const sql =
    `UPDATE jobs SET status = ?, error = ?, updated_at = ?` +
    startedAtClause +
    completedAtClause +
    ` WHERE id = ?`;

  const bindings: unknown[] = [body.status, body.error ?? null, now];
  if (startedAtClause) bindings.push(now);
  if (completedAtClause) bindings.push(now);
  bindings.push(id);

  const result = await env.DB.prepare(sql).bind(...bindings).run();
  if (result.meta.changes === 0) return err("Job not found", 404);

  await logAudit(env, `job_status_${body.status}`, {
    actor: "daemon",
    job_id: id,
    extra: body.error ? { error: body.error } : undefined,
  });
  return json({ ok: true });
}

// ---------------------------------------------------------------------------
// Deletion: daemon-driven receipt + user-triggered immediate delete
// ---------------------------------------------------------------------------

interface DeletionReceiptBody {
  r2_key: string;
  deleted_at: number;
  verified: boolean;
  attempts: number;
  reason: string;
  requestor: string;
}

async function handleDeletionReceipt(
  request: Request,
  jobId: string,
  env: Env
): Promise<Response> {
  if (!isAuthorized(request, env)) return err("Unauthorized", 401);

  const body = await request.json<DeletionReceiptBody>();
  if (!body.r2_key || !body.reason) {
    return err("Missing r2_key or reason", 400);
  }

  const receiptId = crypto.randomUUID();

  // Insert receipt + update uploads.deleted_at in one batch. Both are
  // idempotent under the schema's UNIQUE constraints — repeating the
  // call (e.g. retry storm) won't corrupt state.
  await env.DB.batch([
    env.DB.prepare(
      `INSERT INTO deletion_receipts
         (id, r2_key, requested_at, deleted_at, reason, requestor, job_id, verified)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)`
    ).bind(
      receiptId,
      body.r2_key,
      body.deleted_at,
      body.deleted_at,
      body.reason,
      body.requestor,
      jobId,
      body.verified ? 1 : 0
    ),
    env.DB.prepare(
      `UPDATE uploads SET deleted_at = ?, deletion_receipt_id = ?
       WHERE job_id = ? AND deleted_at IS NULL`
    ).bind(body.deleted_at, receiptId, jobId),
  ]);

  await logAudit(env, body.verified ? "delete_verified" : "delete_unverified", {
    actor: body.requestor,
    job_id: jobId,
    r2_key: body.r2_key,
    extra: { attempts: body.attempts, reason: body.reason },
  });

  return json({ ok: true, receipt_id: receiptId });
}

async function handleUserDelete(jobId: string, env: Env): Promise<Response> {
  // User-triggered immediate deletion. v1: no auth — the job_id IS the
  // capability token (UUID, unguessable). When user accounts land, gate
  // this on the session being the job's owner.
  const upload = await env.DB.prepare(
    `SELECT id, r2_key, deleted_at FROM uploads WHERE job_id = ?`
  )
    .bind(jobId)
    .first<{ id: string; r2_key: string; deleted_at: number | null }>();

  if (!upload) return err("Upload not found", 404);
  if (upload.deleted_at !== null) {
    return json({ ok: true, note: "already deleted", deleted_at: upload.deleted_at });
  }

  // Delete from R2 via the binding (faster than presigned DELETE).
  await env.R2.delete(upload.r2_key);

  // Verify via HeadObject (the binding's .head() returns null on 404).
  const head = await env.R2.head(upload.r2_key);
  const verified = head === null;

  const now = nowSeconds();
  const receiptId = crypto.randomUUID();

  await env.DB.batch([
    env.DB.prepare(
      `INSERT INTO deletion_receipts
         (id, r2_key, requested_at, deleted_at, reason, requestor, job_id, verified)
       VALUES (?, ?, ?, ?, 'user_request', 'user', ?, ?)`
    ).bind(receiptId, upload.r2_key, now, now, jobId, verified ? 1 : 0),
    env.DB.prepare(
      `UPDATE uploads SET deleted_at = ?, deletion_receipt_id = ? WHERE id = ?`
    ).bind(now, receiptId, upload.id),
    // Also abort the job if it's still running.
    env.DB.prepare(
      `UPDATE jobs SET status = 'aborted', updated_at = ?
       WHERE id = ? AND status IN ('queued', 'uploading', 'processing')`
    ).bind(now, jobId),
  ]);

  await logAudit(env, "delete_user_requested", {
    actor: "user",
    job_id: jobId,
    r2_key: upload.r2_key,
    extra: { verified },
  });

  if (!verified) {
    return json(
      { ok: false, receipt_id: receiptId, error: "verification failed" },
      500
    );
  }
  return json({ ok: true, receipt_id: receiptId });
}

// ---------------------------------------------------------------------------
// Results download
// ---------------------------------------------------------------------------

async function handleGetResult(
  id: string,
  filename: string,
  env: Env
): Promise<Response> {
  const job = await env.DB.prepare("SELECT * FROM jobs WHERE id = ?")
    .bind(id)
    .first<JobRow>();
  if (!job) return err("Job not found", 404);
  if (job.status !== "complete") {
    return err(`Job is not complete (status: ${job.status})`, 409);
  }

  const objectKey = `outputs/${id}/${filename}`;
  const object = await env.R2.get(objectKey);
  if (!object) return err("Result file not found", 404);

  const ext = filename.split(".").pop()?.toLowerCase();
  const contentTypeMap: Record<string, string> = {
    tsv: "text/tab-separated-values",
    json: "application/json",
    md: "text/markdown",
    geojson: "application/geo+json",
  };
  const contentType = contentTypeMap[ext ?? ""] ?? "application/octet-stream";

  return new Response(object.body, {
    headers: {
      ...CORS_HEADERS,
      "Content-Type": contentType,
      "Content-Disposition": `attachment; filename="${filename}"`,
    },
  });
}
