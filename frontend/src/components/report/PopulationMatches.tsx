import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'
import type { PopulationMatch } from '@/types/report'

export function PopulationMatches({
  matches,
}: {
  matches: PopulationMatch[]
}) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Top Population Matches</CardTitle>
      </CardHeader>
      <CardContent>
        {matches.length === 0 ? (
          <p className="text-sm text-muted-foreground">No matches.</p>
        ) : (
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead className="text-right">Rank</TableHead>
                <TableHead>Population</TableHead>
                <TableHead className="text-right">n</TableHead>
                <TableHead className="text-right">Mean dist</TableHead>
                <TableHead className="text-right">Min dist</TableHead>
                <TableHead>Date</TableHead>
                <TableHead>Locality / Country</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {matches.map((m) => (
                <TableRow key={m.rank}>
                  <TableCell className="text-right">{m.rank}</TableCell>
                  <TableCell>{m.population}</TableCell>
                  <TableCell className="text-right">{m.n_individuals}</TableCell>
                  <TableCell className="text-right">
                    {m.mean_distance.toFixed(4)}
                  </TableCell>
                  <TableCell className="text-right">
                    {m.min_distance.toFixed(4)}
                  </TableCell>
                  <TableCell>{m.date_display}</TableCell>
                  <TableCell>
                    {[m.locality, m.political_entity]
                      .filter(Boolean)
                      .join(', ')}
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
