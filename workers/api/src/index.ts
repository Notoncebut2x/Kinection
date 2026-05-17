/*
 * D1 Schema — run via `wrangler d1 execute kinection --file=schema.sql`
 *
 * CREATE TABLE IF NOT EXISTS jobs (
 *   id TEXT PRIMARY KEY,
 *   status TEXT NOT NULL DEFAULT 'queued',
 *   upload_key TEXT NOT NULL,
 *   result_prefix TEXT NOT NULL,
 *   created_at INTEGER NOT NULL,
 *   updated_at INTEGER NOT NULL,
 *   error TEXT
 * );
 */

export interface Env {
  R2: R2Bucket;
  DB: D1Database;
  ANALYSIS_QUEUE: Queue;
  MARKER_CACHE: KVNamespace;
  COMPUTE_SERVICE_URL: string;
  COMPUTE_API_KEY: string;
}

interface AnalysisJob {
  jobId: string;
  uploadKey: string;
  resultPrefix: string;
  createdAt: number;
}

interface JobRow {
  id: string;
  status: string;
  upload_key: string;
  result_prefix: string;
  created_at: number;
  updated_at: number;
  error: string | null;
}

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, POST, PATCH, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type, Authorization",
};

function json(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { ...CORS_HEADERS, "Content-Type": "application/json" },
  });
}

function err(message: string, status: number): Response {
  return json({ error: message }, status);
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: CORS_HEADERS });
    }

    const url = new URL(request.url);
    const { pathname } = url;

    try {
      // POST /upload
      if (request.method === "POST" && pathname === "/upload") {
        return handleUpload(request, env);
      }

      // GET /jobs/:id
      const jobMatch = pathname.match(/^\/jobs\/([^/]+)$/);
      if (request.method === "GET" && jobMatch) {
        return handleGetJob(jobMatch[1], env);
      }

      // PATCH /jobs/:id/status
      const statusMatch = pathname.match(/^\/jobs\/([^/]+)\/status$/);
      if (request.method === "PATCH" && statusMatch) {
        return handleUpdateStatus(request, statusMatch[1], env);
      }

      // GET /jobs/:id/results/:filename
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

  async queue(batch: MessageBatch<AnalysisJob>, env: Env): Promise<void> {
    for (const message of batch.messages) {
      const job = message.body;

      const response = await fetch(`${env.COMPUTE_SERVICE_URL}/process`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${env.COMPUTE_API_KEY}`,
        },
        body: JSON.stringify(job),
      });

      if (!response.ok) {
        const text = await response.text().catch(() => "(no body)");
        throw new Error(
          `Compute service returned ${response.status}: ${text}`
        );
      }
    }
  },
};

async function handleUpload(request: Request, env: Env): Promise<Response> {
  const contentType = request.headers.get("Content-Type") ?? "";
  if (!contentType.includes("multipart/form-data")) {
    return err("Expected multipart/form-data", 400);
  }

  const formData = await request.formData();
  const file = formData.get("dna_file");
  if (!file || !(file instanceof File)) {
    return err("Missing dna_file field", 400);
  }

  const jobId = crypto.randomUUID();
  const uploadKey = `uploads/${jobId}/raw.txt`;
  const resultPrefix = `outputs/${jobId}`;
  const now = Date.now();

  await env.R2.put(uploadKey, file.stream(), {
    httpMetadata: { contentType: file.type || "text/plain" },
  });

  await env.DB.prepare(
    `INSERT INTO jobs (id, status, upload_key, result_prefix, created_at, updated_at)
     VALUES (?, 'queued', ?, ?, ?, ?)`
  )
    .bind(jobId, uploadKey, resultPrefix, now, now)
    .run();

  const analysisJob: AnalysisJob = {
    jobId,
    uploadKey,
    resultPrefix,
    createdAt: now,
  };

  await env.ANALYSIS_QUEUE.send(analysisJob);

  return json({ jobId }, 202);
}

async function handleGetJob(id: string, env: Env): Promise<Response> {
  const result = await env.DB.prepare("SELECT * FROM jobs WHERE id = ?")
    .bind(id)
    .first<JobRow>();

  if (!result) {
    return err("Job not found", 404);
  }

  return json(result);
}

async function handleUpdateStatus(
  request: Request,
  id: string,
  env: Env
): Promise<Response> {
  const authHeader = request.headers.get("Authorization");
  if (authHeader !== `Bearer ${env.COMPUTE_API_KEY}`) {
    return err("Unauthorized", 401);
  }

  const body = await request.json<{ status: string; error?: string }>();
  if (!body.status) {
    return err("Missing status field", 400);
  }

  const now = Date.now();

  const result = await env.DB.prepare(
    `UPDATE jobs SET status = ?, error = ?, updated_at = ? WHERE id = ?`
  )
    .bind(body.status, body.error ?? null, now, id)
    .run();

  if (result.meta.changes === 0) {
    return err("Job not found", 404);
  }

  return json({ ok: true });
}

async function handleGetResult(
  id: string,
  filename: string,
  env: Env
): Promise<Response> {
  const job = await env.DB.prepare("SELECT * FROM jobs WHERE id = ?")
    .bind(id)
    .first<JobRow>();

  if (!job) {
    return err("Job not found", 404);
  }

  if (job.status !== "complete") {
    return err(`Job is not complete (status: ${job.status})`, 409);
  }

  const objectKey = `outputs/${id}/${filename}`;
  const object = await env.R2.get(objectKey);

  if (!object) {
    return err("Result file not found", 404);
  }

  const ext = filename.split(".").pop()?.toLowerCase();
  const contentTypeMap: Record<string, string> = {
    tsv: "text/tab-separated-values",
    json: "application/json",
    md: "text/markdown",
  };
  const contentType =
    contentTypeMap[ext ?? ""] ?? "application/octet-stream";

  return new Response(object.body, {
    headers: {
      ...CORS_HEADERS,
      "Content-Type": contentType,
      "Content-Disposition": `attachment; filename="${filename}"`,
    },
  });
}
