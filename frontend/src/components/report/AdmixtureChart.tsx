import {
  Bar,
  BarChart,
  CartesianGrid,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { Admixture } from '@/types/report'

interface Props {
  admixture: Admixture
}

export function AdmixtureChart({ admixture }: Props) {
  const data = Object.entries(admixture.proportions).map(([component, frac]) => {
    const ci = admixture.ci95[component]
    return {
      component,
      percent: +(frac * 100).toFixed(1),
      ciLow: ci ? +(ci[0] * 100).toFixed(1) : null,
      ciHigh: ci ? +(ci[1] * 100).toFixed(1) : null,
      description: admixture.sources[component]?.description ?? component,
    }
  })

  return (
    <Card>
      <CardHeader>
        <CardTitle>Admixture</CardTitle>
        <p className="text-sm text-muted-foreground">
          Ancestry proportions ({admixture.n_snps_used.toLocaleString()} SNPs,{' '}
          {admixture.n_bootstrap} bootstraps)
        </p>
      </CardHeader>
      <CardContent>
        <div className="h-72 w-full">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={data} margin={{ top: 8, right: 8, bottom: 8, left: 8 }}>
              <CartesianGrid strokeDasharray="3 3" vertical={false} />
              <XAxis dataKey="component" fontSize={12} />
              <YAxis unit="%" fontSize={12} />
              <Tooltip
                formatter={(value: number, _name, item) => {
                  const p = item.payload as (typeof data)[number]
                  const ci =
                    p.ciLow != null && p.ciHigh != null
                      ? ` (95% CI ${p.ciLow}–${p.ciHigh}%)`
                      : ''
                  return [`${value}%${ci}`, p.description]
                }}
              />
              <Bar dataKey="percent" fill="hsl(222.2 47.4% 11.2%)" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="mt-4 flex flex-wrap gap-x-6 gap-y-2 text-xs text-muted-foreground">
          {data.map((d) => (
            <div key={d.component}>
              <span className="font-medium text-foreground">{d.component}</span>{' '}
              — {d.description}
            </div>
          ))}
        </div>

        <p className="mt-3 text-xs text-muted-foreground">
          Caveat: model residual = {admixture.residual.toFixed(4)}. Higher
          residuals indicate a poorer fit; interpret proportions with care.
        </p>
      </CardContent>
    </Card>
  )
}
