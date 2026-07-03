import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'
import type { TmrcaMatch, YTmrca } from '@/types/report'

function haplo(m: TmrcaMatch): string {
  return m.ancient_y_haplogroup ?? m.ancient_mt_haplogroup ?? '—'
}

function tmrcaCell(m: TmrcaMatch): string {
  if (m.tmrca_yr == null) return '—'
  const yr = Math.round(m.tmrca_yr).toLocaleString()
  if (m.tmrca_lo_95 == null || m.tmrca_hi_95 == null) return `${yr} yr`
  return `${yr} yr (${Math.round(m.tmrca_lo_95).toLocaleString()}–${Math.round(
    m.tmrca_hi_95,
  ).toLocaleString()})`
}

function MatchTable({ matches }: { matches: TmrcaMatch[] }) {
  if (matches.length === 0) {
    return <p className="text-sm text-muted-foreground">No matches.</p>
  }
  return (
    <Table>
      <TableHeader>
        <TableRow>
          <TableHead>Population</TableHead>
          <TableHead>Date</TableHead>
          <TableHead>Haplogroup</TableHead>
          <TableHead className="text-right">n_diff</TableHead>
          <TableHead className="text-right">TMRCA (95% CI)</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {matches.map((m) => (
          <TableRow key={m.genetic_id}>
            <TableCell>{m.population}</TableCell>
            <TableCell>{m.date_display}</TableCell>
            <TableCell>{haplo(m)}</TableCell>
            <TableCell className="text-right">{m.n_diff}</TableCell>
            <TableCell className="text-right">{tmrcaCell(m)}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  )
}

export function TmrcaSection({ tmrca }: { tmrca: YTmrca }) {
  if (tmrca.skipped) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Y-DNA TMRCA</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-sm text-muted-foreground">
            TMRCA analysis was skipped.
          </p>
        </CardContent>
      </Card>
    )
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle>Y-DNA TMRCA</CardTitle>
        <p className="text-sm text-muted-foreground">
          Method: {tmrca.y_method} · Modern Y-haplogroup{' '}
          {tmrca.modern_y_haplogroup} ({tmrca.modern_y_confidence})
        </p>
      </CardHeader>
      <CardContent className="space-y-6">
        <div>
          <h4 className="mb-2 text-sm font-medium">Y matches</h4>
          <MatchTable matches={tmrca.matches} />
        </div>

        <div>
          <h4 className="mb-2 text-sm font-medium">mt matches</h4>
          {tmrca.mt_skipped ? (
            <p className="text-sm text-muted-foreground">
              mtDNA TMRCA skipped
              {tmrca.mt_skip_reason ? `: ${tmrca.mt_skip_reason}` : '.'}
            </p>
          ) : (
            <MatchTable matches={tmrca.mt_matches} />
          )}
        </div>
      </CardContent>
    </Card>
  )
}
