// Types derived from the real report.json sample structure.

export type Assessment = 'PASS' | 'QUESTIONABLE' | 'CRITICAL'

export interface HaplogroupCall {
  value: string
  confidence: string
  notes?: string | null
}

export interface ModernIndividual {
  snps_called: number
  y_haplogroup: HaplogroupCall
  mt_haplogroup: HaplogroupCall
}

export interface AncientDataset {
  individuals: number
  snps: number
  snp_overlap: number
}

export interface AdmixtureSource {
  description: string
  n_individuals: number
}

export interface Admixture {
  proportions: Record<string, number>
  ci95: Record<string, [number, number]>
  residual: number
  n_snps_used: number
  n_bootstrap: number
  sources: Record<string, AdmixtureSource>
}

export interface TmrcaMatch {
  genetic_id: string
  population: string
  ancient_y_haplogroup?: string
  ancient_mt_haplogroup?: string
  locality: string
  lat: number
  lon: number
  date_bp: number
  date_display: string
  n_y_sites?: number
  n_mt_sites?: number
  n_diff: number
  diff_rate: number | null
  tmrca_yr: number | null
  tmrca_lo_95: number | null
  tmrca_hi_95: number | null
  below_sample_age: boolean
}

export interface YTmrca {
  skipped: boolean
  y_method: string
  mt_method: string
  modern_y_haplogroup: string
  modern_y_confidence: string
  matches: TmrcaMatch[]
  mt_matches: TmrcaMatch[]
  mt_skipped: boolean
  mt_skip_reason: string | null
}

export interface HaplogroupMatch {
  genetic_id: string
  population: string
  locality: string
  political_entity: string
  lat: number
  lon: number
  date_bp: number
  date_display: string
  molecular_sex: string
  ancient_y_haplogroup: string
  ancient_mt_haplogroup: string
  match_type: string
  combined_score: number
  assessment: Assessment
}

export interface PopulationMatch {
  rank: number
  population: string
  n_individuals: number
  mean_distance: number
  min_distance: number
  median_distance: number
  date_display: string
  locality: string
  political_entity: string
  lat: number
  lon: number
}

export interface IndividualMatch {
  genetic_id: string
  population: string
  locality: string
  political_entity: string
  lat: number
  lon: number
  date_bp: number
  date_display: string
  y_haplogroup: string
  mt_haplogroup: string
  snps_compared: number
  asd: number
  rank: number
}

export interface Pca {
  modern_coords: number[] | null
  variance_explained: number[] | null
  n_components: number
}

export interface Report {
  label: string
  generated_at: string
  schema_version: string
  modern_individual: ModernIndividual
  ancient_dataset: AncientDataset
  admixture: Admixture
  y_tmrca: YTmrca
  haplogroup_matches: HaplogroupMatch[]
  top_population_matches: PopulationMatch[]
  top_individual_matches: IndividualMatch[]
  pca: Pca
  anomalies: string[]
}

// map_data.geojson is a standard GeoJSON FeatureCollection — typed loosely.
export interface GeoJsonFeatureCollection {
  type: 'FeatureCollection'
  features: Array<{
    type: 'Feature'
    geometry: { type: string; coordinates: unknown }
    properties: Record<string, unknown>
  }>
}
