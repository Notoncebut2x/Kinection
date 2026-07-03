import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'
import type { IndividualMatch } from '@/types/report'

export function IndividualMatches({
  matches,
}: {
  matches: IndividualMatch[]
}) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Top Individual Matches</CardTitle>
      </CardHeader>
      <CardContent>
        {matches.length === 0 ? (
          <p className="text-sm text-muted-foreground">No matches.</p>
        ) : (
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead className="text-right">Rank</TableHead>
                <TableHead>Genetic ID</TableHead>
                <TableHead>Population</TableHead>
                <TableHead className="text-right">ASD</TableHead>
                <TableHead className="text-right">SNPs</TableHead>
                <TableHead>Date</TableHead>
                <TableHead>Country</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {matches.map((m) => (
                <TableRow key={m.genetic_id}>
                  <TableCell className="text-right">{m.rank}</TableCell>
                  <TableCell className="font-mono text-xs">
                    {m.genetic_id}
                  </TableCell>
                  <TableCell>{m.population}</TableCell>
                  <TableCell className="text-right">
                    {m.asd.toFixed(4)}
                  </TableCell>
                  <TableCell className="text-right">
                    {m.snps_compared.toLocaleString()}
                  </TableCell>
                  <TableCell>{m.date_display}</TableCell>
                  <TableCell>{m.political_entity}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
      </CardContent>
    </Card>
  )
}
