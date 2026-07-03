import { Badge } from '@/components/ui/badge'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { AncientDataset, ModernIndividual } from '@/types/report'

interface Props {
  label: string
  generatedAt: string
  modern: ModernIndividual
  ancient: AncientDataset
}

export function ReportHeader({ label, generatedAt, modern, ancient }: Props) {
  return (
    <Card>
      <CardHeader>
        <CardTitle className="text-2xl">{label || 'Report'}</CardTitle>
        <p className="text-sm text-muted-foreground">
          Generated {generatedAt}
        </p>
      </CardHeader>
      <CardContent className="flex flex-wrap gap-6">
        <div className="flex flex-col gap-1">
          <span className="text-xs uppercase text-muted-foreground">
            Y-haplogroup
          </span>
          <div className="flex items-center gap-2">
            <Badge>{modern.y_haplogroup.value}</Badge>
            <span className="text-xs text-muted-foreground">
              {modern.y_haplogroup.confidence}
            </span>
          </div>
          {modern.y_haplogroup.notes ? (
            <span className="text-xs text-muted-foreground">
              {modern.y_haplogroup.notes}
            </span>
          ) : null}
        </div>

        <div className="flex flex-col gap-1">
          <span className="text-xs uppercase text-muted-foreground">
            mt-haplogroup
          </span>
          <div className="flex items-center gap-2">
            <Badge variant="secondary">{modern.mt_haplogroup.value}</Badge>
            <span className="text-xs text-muted-foreground">
              {modern.mt_haplogroup.confidence}
            </span>
          </div>
        </div>

        <div className="flex flex-col gap-1">
          <span className="text-xs uppercase text-muted-foreground">
            SNPs called
          </span>
          <span className="text-lg font-semibold">
            {modern.snps_called.toLocaleString()}
          </span>
        </div>

        <div className="flex flex-col gap-1">
          <span className="text-xs uppercase text-muted-foreground">
            SNP overlap
          </span>
          <span className="text-lg font-semibold">
            {ancient.snp_overlap.toLocaleString()}
          </span>
          <span className="text-xs text-muted-foreground">
            vs {ancient.individuals.toLocaleString()} ancient individuals
          </span>
        </div>
      </CardContent>
    </Card>
  )
}
