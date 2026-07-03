import { useEffect, useState } from 'react'
import { BrowserRouter, Link, Route, Routes } from 'react-router-dom'
import { api } from '@/api/client'
import { UploadPage } from '@/pages/UploadPage'
import { StatusPage } from '@/pages/StatusPage'
import { ReportPage } from '@/pages/ReportPage'
import type { DatasetVersion } from '@/types/api'

function useDatasetVersion() {
  const [version, setVersion] = useState<string | null>(null)
  useEffect(() => {
    let cancelled = false
    api
      .getDatasetVersion()
      .then((m: DatasetVersion) => {
        if (cancelled) return
        const v =
          typeof m.version === 'string'
            ? m.version
            : typeof m.aadr_version === 'string'
              ? (m.aadr_version as string)
              : null
        setVersion(v)
      })
      .catch(() => {
        /* dataset version is best-effort */
      })
    return () => {
      cancelled = true
    }
  }, [])
  return version
}

function Layout({ children }: { children: React.ReactNode }) {
  const datasetVersion = useDatasetVersion()
  return (
    <div className="flex min-h-screen flex-col">
      <header className="border-b">
        <div className="container flex h-14 items-center justify-between">
          <Link to="/" className="text-lg font-semibold tracking-tight">
            Kinection
          </Link>
          <span className="text-xs text-muted-foreground">
            AADR dataset: {datasetVersion ?? 'unknown'}
          </span>
        </div>
      </header>

      <main className="container flex-1 py-8">{children}</main>

      <footer className="border-t">
        <div className="container flex h-12 items-center justify-between text-xs text-muted-foreground">
          <span>Kinection — personal ancient-DNA comparison</span>
          <span>AADR dataset: {datasetVersion ?? 'unknown'}</span>
        </div>
      </footer>
    </div>
  )
}

export default function App() {
  return (
    <BrowserRouter>
      <Layout>
        <Routes>
          <Route path="/" element={<UploadPage />} />
          <Route path="/jobs/:id" element={<StatusPage />} />
          <Route path="/jobs/:id/report" element={<ReportPage />} />
        </Routes>
      </Layout>
    </BrowserRouter>
  )
}
