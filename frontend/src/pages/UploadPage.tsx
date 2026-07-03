import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { api, ApiError } from '@/api/client'
import { Button } from '@/components/ui/button'
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from '@/components/ui/card'

async function sha256Hex(file: File): Promise<string | undefined> {
  if (!globalThis.crypto?.subtle) return undefined
  const buf = await file.arrayBuffer()
  const digest = await crypto.subtle.digest('SHA-256', buf)
  return Array.from(new Uint8Array(digest))
    .map((b) => b.toString(16).padStart(2, '0'))
    .join('')
}

export function UploadPage() {
  const navigate = useNavigate()
  const [file, setFile] = useState<File | null>(null)
  const [label, setLabel] = useState('')
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState<string | null>(null)

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault()
    if (!file) return
    setBusy(true)
    setError(null)
    try {
      const sha256 = await sha256Hex(file)
      const created = await api.createUpload({
        label: label.trim() || undefined,
        format: 'txt',
        sha256,
        size_bytes: file.size,
      })

      // PUT the raw bytes straight to the presigned R2 URL (not via Worker).
      await api.putFileToUploadUrl(created.upload_url, file)

      // Ask the backend to flip 'uploading' -> 'queued'. This route does not
      // exist yet; markUploadComplete swallows a 404 as a known gap.
      await api.markUploadComplete(created.job_id)

      navigate(`/jobs/${created.job_id}`)
    } catch (err) {
      const msg =
        err instanceof ApiError
          ? `Upload failed (${err.status}): ${err.message}`
          : err instanceof Error
            ? err.message
            : 'Upload failed'
      setError(msg)
      setBusy(false)
    }
  }

  return (
    <div className="mx-auto max-w-xl">
      <Card>
        <CardHeader>
          <CardTitle>Upload your DNA file</CardTitle>
          <CardDescription>
            Upload a raw autosomal DNA file (.txt) to compare against the
            ancient-DNA dataset.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleSubmit} className="space-y-4">
            <div className="space-y-1.5">
              <label htmlFor="label" className="text-sm font-medium">
                Label (optional)
              </label>
              <input
                id="label"
                type="text"
                value={label}
                onChange={(e) => setLabel(e.target.value)}
                placeholder="e.g. My sample"
                className="flex h-10 w-full rounded-md border border-input bg-background px-3 py-2 text-sm focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
              />
            </div>

            <div className="space-y-1.5">
              <label htmlFor="file" className="text-sm font-medium">
                Raw DNA file (.txt)
              </label>
              <input
                id="file"
                type="file"
                accept=".txt,text/plain"
                onChange={(e) => setFile(e.target.files?.[0] ?? null)}
                className="flex w-full rounded-md border border-input bg-background px-3 py-2 text-sm file:mr-4 file:rounded file:border-0 file:bg-secondary file:px-3 file:py-1 file:text-sm"
              />
            </div>

            {error ? (
              <p className="text-sm text-destructive">{error}</p>
            ) : null}

            <Button type="submit" disabled={!file || busy} className="w-full">
              {busy ? 'Uploading…' : 'Upload & analyze'}
            </Button>
          </form>
        </CardContent>
      </Card>
    </div>
  )
}
