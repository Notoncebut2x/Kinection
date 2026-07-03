import { useCallback, useEffect, useRef, useState } from 'react'
import { Link, useNavigate, useParams } from 'react-router-dom'
import { api, ApiError } from '@/api/client'
import { Button } from '@/components/ui/button'
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
} from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import type { JobRow, JobStatus } from '@/types/api'

const POLL_MS = 3000

const STATUS_LABEL: Record<JobStatus, string> = {
  uploading: 'Uploading',
  queued: 'Queued',
  processing: 'Processing',
  complete: 'Complete',
  failed: 'Failed',
  aborted: 'Aborted',
}

const IN_PROGRESS: JobStatus[] = ['uploading', 'queued', 'processing']

function statusVariant(
  s: JobStatus,
): 'default' | 'secondary' | 'destructive' | 'success' {
  if (s === 'complete') return 'success'
  if (s === 'failed' || s === 'aborted') return 'destructive'
  if (s === 'processing' || s === 'queued') return 'default'
  return 'secondary'
}

export function StatusPage() {
  const { id } = useParams<{ id: string }>()
  const navigate = useNavigate()
  const [job, setJob] = useState<JobRow | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [deleting, setDeleting] = useState(false)
  const [deleted, setDeleted] = useState(false)
  const timer = useRef<number | null>(null)

  const poll = useCallback(async () => {
    if (!id) return
    try {
      const j = await api.getJob(id)
      setJob(j)
      setError(null)
      if (IN_PROGRESS.includes(j.status)) {
        timer.current = window.setTimeout(poll, POLL_MS)
      }
    } catch (err) {
      const msg =
        err instanceof ApiError
          ? `Could not load job (${err.status}): ${err.message}`
          : 'Could not load job'
      setError(msg)
      timer.current = window.setTimeout(poll, POLL_MS)
    }
  }, [id])

  useEffect(() => {
    poll()
    return () => {
      if (timer.current) window.clearTimeout(timer.current)
    }
  }, [poll])

  async function handleDelete() {
    if (!id) return
    if (!window.confirm('Permanently delete your uploaded file?')) return
    setDeleting(true)
    try {
      await api.deleteUpload(id)
      setDeleted(true)
      if (timer.current) window.clearTimeout(timer.current)
    } catch (err) {
      const msg =
        err instanceof ApiError
          ? `Delete failed (${err.status}): ${err.message}`
          : 'Delete failed'
      setError(msg)
    } finally {
      setDeleting(false)
    }
  }

  const status = job?.status
  const inProgress = status ? IN_PROGRESS.includes(status) : true

  return (
    <div className="mx-auto max-w-xl space-y-4">
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <CardTitle>Job status</CardTitle>
            {status ? (
              <Badge variant={statusVariant(status)}>
                {STATUS_LABEL[status]}
              </Badge>
            ) : null}
          </div>
          <p className="font-mono text-xs text-muted-foreground">{id}</p>
        </CardHeader>
        <CardContent className="space-y-4">
          {deleted ? (
            <p className="text-sm">Your upload was deleted.</p>
          ) : (
            <>
              {inProgress ? (
                <div className="flex items-center gap-3">
                  <span
                    className="inline-block h-5 w-5 animate-spin rounded-full border-2 border-muted border-t-primary"
                    aria-hidden
                  />
                  <span className="text-sm text-muted-foreground">
                    {status === 'uploading'
                      ? 'Waiting for upload to be finalized…'
                      : 'Working on your analysis…'}
                  </span>
                </div>
              ) : null}

              {status === 'uploading' ? (
                <p className="rounded-md bg-muted p-3 text-xs text-muted-foreground">
                  Note: the job may sit in “uploading” until the backend
                  upload-complete route is implemented (see README).
                </p>
              ) : null}

              {status === 'complete' ? (
                <div className="space-y-3">
                  <p className="text-sm">Your report is ready.</p>
                  <Button onClick={() => navigate(`/jobs/${id}/report`)}>
                    View report
                  </Button>
                </div>
              ) : null}

              {(status === 'failed' || status === 'aborted') && job?.error ? (
                <p className="text-sm text-destructive">{job.error}</p>
              ) : null}

              {error ? (
                <p className="text-sm text-destructive">{error}</p>
              ) : null}

              <div className="pt-2">
                <Button
                  variant="destructive"
                  onClick={handleDelete}
                  disabled={deleting}
                >
                  {deleting ? 'Deleting…' : 'Delete my upload'}
                </Button>
              </div>
            </>
          )}
        </CardContent>
      </Card>

      <Link
        to="/"
        className="text-sm text-muted-foreground underline-offset-4 hover:underline"
      >
        ← Upload another file
      </Link>
    </div>
  )
}
