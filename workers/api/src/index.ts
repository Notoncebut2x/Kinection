export interface Env {
  R2: R2Bucket;
  DB: D1Database;
  MARKER_CACHE: KVNamespace;
  COMPUTE_API_KEY: string;
}

interface JobRow {
  id: string;
  status: string;
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

function isAuthorized(request: Request, env: Env): boolean {
  return request.headers.get("Authorization") === `Bearer ${env.COMPUTE_API_KEY}`;
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: CORS_HEADERS });
    }

    const url = new URL(request.url);
    const { pathname } = url;

    try {
      // POST /jobs — web app creates a new job
      if (request.method === "POST" && pathname === "/jobs") {
        return handleCreateJob(env);
      }

      // GET /jobs?status=queued — daemon polls for pending jobs (requires auth)
      if (request.method === "GET" && pathname === "/jobs") {
        return handleListJobs(request, url, env);
      }

      // GET /jobs/:id — web app polls job status
      const jobMatch = pathname.match(/^\/jobs\/([^/]+)$/);
      if (request.method === "GET" && jobMatch) {
        return handleGetJob(jobMatch[1], env);
      }

      // PATCH /jobs/:id/status — daemon reports job outcome (requires auth)
      const statusMatch = pathname.match(/^\/jobs\/([^/]+)\/status$/);
      if (request.method === "PATCH" && statusMatch) {
        return handleUpdateStatus(request, statusMatch[1], env);
      }

      // GET /jobs/:id/results/:filename — web app fetches result file from R2
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

async function handleCreateJob(env: Env): Promise<Response> {
  const id = crypto.randomUUID();
  const now = Date.now();

  await env.DB.prepare(
    `INSERT INTO jobs (id, status, created_at, updated_at) VALUES (?, 'queued', ?, ?)`
  )
    .bind(id, now, now)
    .run();

  return json({ id, status: "queued", created_at: now }, 201);
}

async function handleListJobs(
  request: Request,
  url: URL,
  env: Env
): Promise<Response> {
  if (!isAuthorized(request, env)) {
    return err("Unauthorized", 401);
  }

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
  if (!isAuthorized(request, env)) {
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
  const contentType = contentTypeMap[ext ?? ""] ?? "application/octet-stream";

  return new Response(object.body, {
    headers: {
      ...CORS_HEADERS,
      "Content-Type": contentType,
      "Content-Disposition": `attachment; filename="${filename}"`,
    },
  });
}
