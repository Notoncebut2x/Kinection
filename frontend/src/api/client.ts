import type {
  CreateUploadRequest,
  CreateUploadResponse,
  DatasetVersion,
  DeleteUploadResponse,
  JobRow,
} from '@/types/api'
import type { GeoJsonFeatureCollection, Report } from '@/types/report'

const API_BASE_URL: string =
  import.meta.env.VITE_API_BASE_URL ?? 'http://localhost:8787'

class ApiError extends Error {
  status: number
  constructor(status: number, message: string) {
    super(message)
    this.name = 'ApiError'
    this.status = status
  }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  // NOTE: We deliberately do NOT send Authorization headers from the browser.
  // The daemon endpoints are auth-gated and not meant for the frontend.
  const res = await fetch(`${API_BASE_URL}${path}`, {
    ...init,
    headers: {
      'Content-Type': 'application/json',
      ...(init?.headers ?? {}),
    },
  })
  if (!res.ok) {
    const body = await res.text().catch(() => '')
    throw new ApiError(res.status, body || `${res.status} ${res.statusText}`)
  }
  return (await res.json()) as T
}

export const api = {
  baseUrl: API_BASE_URL,

  /** POST /uploads/url — create a job + presigned R2 PUT URL. */
  createUpload(body: CreateUploadRequest): Promise<CreateUploadResponse> {
    return request<CreateUploadResponse>('/uploads/url', {
      method: 'POST',
      body: JSON.stringify(body),
    })
  },

  /**
   * PUT the raw file bytes directly to the presigned R2 URL.
   * This does NOT go through the Worker and carries no Authorization header.
   */
  async putFileToUploadUrl(uploadUrl: string, file: File): Promise<void> {
    const res = await fetch(uploadUrl, { method: 'PUT', body: file })
    if (!res.ok) {
      throw new ApiError(res.status, `Upload PUT failed: ${res.status}`)
    }
  },

  /**
   * Flip the job from 'uploading' -> 'queued' after the PUT succeeds.
   *
   * TODO(backend): POST /jobs/:id/upload-complete route needed to flip
   * 'uploading' -> 'queued'. This endpoint does NOT exist yet, so we catch a
   * 404 and treat it as known-pending-backend-work rather than a hard failure.
   */
  async markUploadComplete(jobId: string): Promise<void> {
    try {
      await request<unknown>(`/jobs/${jobId}/upload-complete`, {
        method: 'POST',
      })
    } catch (err) {
      if (err instanceof ApiError && err.status === 404) {
        console.warn(
          `[kinection] POST /jobs/${jobId}/upload-complete returned 404 — ` +
            'backend route not implemented yet. Job will remain in ' +
            "'uploading' until the backend lands this route.",
        )
        return
      }
      throw err
    }
  },

  /** GET /jobs/:id — current job row. */
  getJob(jobId: string): Promise<JobRow> {
    return request<JobRow>(`/jobs/${jobId}`)
  },

  /** GET /jobs/:id/results/report.json (only when status === 'complete'). */
  getReport(jobId: string): Promise<Report> {
    return request<Report>(`/jobs/${jobId}/results/report.json`)
  },

  /** GET /jobs/:id/results/map_data.geojson (only when status === 'complete'). */
  getMapData(jobId: string): Promise<GeoJsonFeatureCollection> {
    return request<GeoJsonFeatureCollection>(
      `/jobs/${jobId}/results/map_data.geojson`,
    )
  },

  /** DELETE /jobs/:id/upload — user-triggered immediate deletion. */
  deleteUpload(jobId: string): Promise<DeleteUploadResponse> {
    return request<DeleteUploadResponse>(`/jobs/${jobId}/upload`, {
      method: 'DELETE',
    })
  },

  /** GET /dataset/version — AADR dataset version manifest. */
  getDatasetVersion(): Promise<DatasetVersion> {
    return request<DatasetVersion>('/dataset/version')
  },
}

export { ApiError }
