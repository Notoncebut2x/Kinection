import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { Pca } from '@/types/report'

export function PcaPlaceholder({ pca }: { pca: Pca }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>PCA Projection</CardTitle>
      </CardHeader>
      <CardContent>
        <div className="flex h-48 flex-col items-center justify-center gap-1 rounded-md border border-dashed text-center">
          <p className="text-sm font-medium">PCA projection — coming soon</p>
          <p className="text-xs text-muted-foreground">
            {pca.n_components} component(s) configured; projection output is
            currently empty.
          </p>
        </div>
      </CardContent>
    </Card>
  )
}
