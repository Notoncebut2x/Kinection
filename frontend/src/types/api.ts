// Types matching the Kinection Cloudflare Worker API contract.

export type JobStatus =
  | 'uploading'
  | 'queued'
  | 'processing'
  | 'complete'
  | 'failed'
  | 'aborted'

/** Row returned by GET /jobs/:id */
export interface JobRow {
  id: string
  user_id: string | null
  status: JobStatus
  label: string | null
  aadr_version: string | null
  // Unix seconds; nullable until reached.
  created_at: number
  updated_at: number
  started_at: number | null
  completed_at: number | null
  error: string | null
}

/** Body for POST /uploads/url */
export interface CreateUploadRequest {
  label?: string
  format?: string
  sha256?: string
  size_bytes?: number
}

/** 201 response from POST /uploads/url */
export interface CreateUploadResponse {
  job_id: string
  upload_id: string
  upload_key: string
  /** Presigned R2 PUT URL — the raw file is PUT directly here. */
  upload_url: string
  /** Unix seconds at which the presigned URL expires. */
  expires_at: number
}

/** Response from DELETE /jobs/:id/upload */
export interface DeleteUploadResponse {
  ok: boolean
  receipt_id: string
}

/** GET /dataset/version manifest (shape is backend-owned; kept loose). */
export interface DatasetVersion {
  version?: string
  [key: string]: unknown
}
