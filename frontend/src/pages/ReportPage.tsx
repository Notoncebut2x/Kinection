import { useEffect, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { api, ApiError } from '@/api/client'
import { AdmixtureChart } from '@/components/report/AdmixtureChart'
import { HaplogroupMatches } from '@/components/report/HaplogroupMatches'
import { IndividualMatches } from '@/components/report/IndividualMatches'
import { MapPlaceholder } from '@/components/report/MapPlaceholder'
import { PcaPlaceholder } from '@/components/report/PcaPlaceholder'
import { PopulationMatches } from '@/components/report/PopulationMatches'
import { ReportHeader } from '@/components/report/ReportHeader'
import { TmrcaSection } from '@/components/report/TmrcaSection'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { GeoJsonFeatureCollection, Report } from '@/types/report'

export function ReportPage() {
  const { id } = useParams<{ id: string }>()
  const [report, setReport] = useState<Report | null>(null)
  const [geojson, setGeojson] = useState<GeoJsonFeatureCollection | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    if (!id) return
    let cancelled = false
    setLoading(true)
    ;(async () => {
      try {
        const rep = await api.getReport(id)
        if (cancelled) return
        setReport(rep)
        // Load the map data too; failure here is non-fatal for the report.
        try {
          const geo = await api.getMapData(id)
          if (!cancelled) setGeojson(geo)
        } catch (geoErr) {
          console.warn('[kinection] map_data.geojson unavailable:', geoErr)
        }
      } catch (err) {
        if (cancelled) return
        const msg =
          err instanceof ApiError
            ? err.status === 409
              ? 'The report is not ready yet (job is not complete).'
              : `Could not load report (${err.status}): ${err.message}`
            : 'Could not load report'
        setError(msg)
      } finally {
        if (!cancelled) setLoading(false)
      }
    })()
    return () => {
      cancelled = true
    }
  }, [id])

  if (loading) {
    return (
      <div className="mx-auto max-w-xl">
        <Card>
          <CardContent className="flex items-center gap-3 p-6">
            <span
              className="inline-block h-5 w-5 animate-spin rounded-full border-2 border-muted border-t-primary"
              aria-hidden
            />
            <span className="text-sm text-muted-foreground">
              Loading report…
            </span>
          </CardContent>
        </Card>
      </div>
    )
  }

  if (error || !report) {
    return (
      <div className="mx-auto max-w-xl space-y-4">
        <Card>
          <CardHeader>
            <CardTitle>Report unavailable</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-destructive">
              {error ?? 'No report found.'}
            </p>
          </CardContent>
        </Card>
        <Link
          to={`/jobs/${id}`}
          className="text-sm text-muted-foreground underline-offset-4 hover:underline"
        >
          ← Back to status
        </Link>
      </div>
    )
  }

  return (
    <div className="mx-auto max-w-5xl space-y-6">
      <ReportHeader
        label={report.label}
        generatedAt={report.generated_at}
        modern={report.modern_individual}
        ancient={report.ancient_dataset}
      />

      {report.anomalies.length > 0 ? (
        <Card>
          <CardHeader>
            <CardTitle>Anomalies</CardTitle>
          </CardHeader>
          <CardContent>
            <ul className="list-disc space-y-1 pl-5 text-sm text-muted-foreground">
              {report.anomalies.map((a, i) => (
                <li key={i}>{a}</li>
              ))}
            </ul>
          </CardContent>
        </Card>
      ) : null}

      <AdmixtureChart admixture={report.admixture} />
      <TmrcaSection tmrca={report.y_tmrca} />
      <HaplogroupMatches matches={report.haplogroup_matches} />
      <PopulationMatches matches={report.top_population_matches} />
      <IndividualMatches matches={report.top_individual_matches} />
      <MapPlaceholder geojson={geojson} />
      <PcaPlaceholder pca={report.pca} />

      <Link
        to={`/jobs/${id}`}
        className="text-sm text-muted-foreground underline-offset-4 hover:underline"
      >
        ← Back to status
      </Link>
    </div>
  )
}
