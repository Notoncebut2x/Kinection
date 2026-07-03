import { useEffect } from 'react'
import { MapContainer, TileLayer } from 'react-leaflet'
import 'leaflet/dist/leaflet.css'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { GeoJsonFeatureCollection } from '@/types/report'

interface Props {
  geojson: GeoJsonFeatureCollection | null
}

export function MapPlaceholder({ geojson }: Props) {
  useEffect(() => {
    if (geojson) {
      // TODO: plot map_data.geojson features as markers
      console.log('[kinection] map_data.geojson loaded:', geojson)
    }
  }, [geojson])

  return (
    <Card>
      <CardHeader>
        <CardTitle>Map</CardTitle>
        <p className="text-sm text-muted-foreground">
          Top autosomal + Y-TMRCA matches (marker plotting coming soon).
        </p>
      </CardHeader>
      <CardContent>
        <div className="h-96 w-full overflow-hidden rounded-md border">
          <MapContainer
            center={[30, 15]}
            zoom={2}
            className="h-full w-full"
            scrollWheelZoom={false}
          >
            <TileLayer
              attribution='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
              url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
            />
            {/* TODO: plot map_data.geojson features as markers */}
          </MapContainer>
        </div>
      </CardContent>
    </Card>
  )
}
