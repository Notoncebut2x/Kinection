import { Badge } from '@/components/ui/badge'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'
import type { Assessment, HaplogroupMatch } from '@/types/report'

function assessmentVariant(
  a: Assessment,
): 'success' | 'warning' | 'destructive' {
  switch (a) {
    case 'PASS':
      return 'success'
    case 'QUESTIONABLE':
      return 'warning'
    case 'CRITICAL':
      return 'destructive'
  }
}

// The AADR's own data-quality verdict for the ANCIENT sample (contamination,
// coverage, damage, duplicate/relatedness checks) — not a rating of you or
// the match. Shown as a hover tooltip on each badge.
const ASSESSMENT_HELP: Record<Assessment, string> = {
  PASS: 'AADR data quality: PASS — the ancient sample cleared quality checks.',
  QUESTIONABLE:
    'AADR data quality: QUESTIONABLE — borderline (e.g. contamination or low coverage); weigh this match less.',
  CRITICAL:
    'AADR data quality: CRITICAL — a serious quality problem (e.g. contamination); the shared haplogroup may be an artifact.',
}

export function HaplogroupMatches({ matches }: { matches: HaplogroupMatch[] }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Haplogroup Matches</CardTitle>
        <p className="text-sm text-muted-foreground">
          Ancient individuals who share your Y-DNA and/or mtDNA lineage. The{' '}
          <span className="font-medium">Assessment</span> badge is the AADR's
          data-quality verdict for that ancient sample (hover for details) — it
          reflects how much to trust the sample's data, not the strength of the
          match.
        </p>
      </CardHeader>
      <CardContent>
        {matches.length === 0 ? (
          <div className="space-y-2 text-sm text-muted-foreground">
            <p>No reliable haplogroup matches.</p>
            <p>
              This usually means your Y-DNA and/or mtDNA haplogroup couldn't be
              confidently resolved from this file — consumer arrays (AncestryDNA,
              23andMe) cover only a small fraction of the markers that define the
              Y and mtDNA trees. Rather than list matches to a deep
              macro-haplogroup (which would be misleading), they're omitted here.
              Your <span className="font-medium">admixture</span> and{' '}
              <span className="font-medium">population matches</span> are the more
              reliable read of your ancestry from array data; for precise
              haplogroups, a dedicated Y or full-mtDNA test is needed.
            </p>
          </div>
        ) : (
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Genetic ID</TableHead>
                <TableHead>Population</TableHead>
                <TableHead>Country</TableHead>
                <TableHead>Date</TableHead>
                <TableHead>Y / mt</TableHead>
                <TableHead className="text-right">Score</TableHead>
                <TableHead>Assessment</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {matches.map((m) => (
                <TableRow key={m.genetic_id}>
                  <TableCell className="font-mono text-xs">
                    {m.genetic_id}
                  </TableCell>
                  <TableCell>{m.population}</TableCell>
                  <TableCell>{m.political_entity}</TableCell>
                  <TableCell>{m.date_display}</TableCell>
                  <TableCell>
                    {m.ancient_y_haplogroup} / {m.ancient_mt_haplogroup}
                  </TableCell>
                  <TableCell className="text-right">
                    {m.combined_score.toFixed(2)}
                  </TableCell>
                  <TableCell>
                    <Badge
                      variant={assessmentVariant(m.assessment)}
                      title={ASSESSMENT_HELP[m.assessment]}
                    >
                      {m.assessment}
                    </Badge>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
        {matches.length > 0 ? (
          <div className="mt-3 flex flex-wrap gap-x-4 gap-y-1 text-xs text-muted-foreground">
            <span>Assessment (AADR sample quality):</span>
            <span><span className="font-medium text-foreground">PASS</span> — cleared QC</span>
            <span><span className="font-medium text-foreground">QUESTIONABLE</span> — borderline</span>
            <span><span className="font-medium text-foreground">CRITICAL</span> — serious quality issue</span>
          </div>
        ) : null}
      </CardContent>
    </Card>
  )
}
