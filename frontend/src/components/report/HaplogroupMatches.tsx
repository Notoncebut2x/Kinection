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

export function HaplogroupMatches({ matches }: { matches: HaplogroupMatch[] }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Haplogroup Matches</CardTitle>
      </CardHeader>
      <CardContent>
        {matches.length === 0 ? (
          <p className="text-sm text-muted-foreground">No matches.</p>
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
                    <Badge variant={assessmentVariant(m.assessment)}>
                      {m.assessment}
                    </Badge>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
      </CardContent>
    </Card>
  )
}
