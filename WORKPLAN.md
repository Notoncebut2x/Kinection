# DNA Lineage Platform — Workplan

**Project:** Web platform for personal DNA comparison against ancient population datasets
**Reference dataset:** Allen Ancient DNA Resource v62.0 (17,629 individuals, 406,570 SNPs)
**Methodology reference:** [Uniparental analysis of Deep Maniot Greeks, *Communications Biology* 2026](https://www.nature.com/articles/s42003-026-09597-9)
**Started:** April 2026

---

## Project Vision

A web application where users upload their raw AncestryDNA, 23andMe, or similar genotype files and receive a report showing:
- Their closest ancient population matches by haplogroup and SNP similarity
- Paternal lineage (Y-DNA haplogroup) and maternal lineage (mtDNA haplogroup) placed on the ancient phylogenetic tree
- TMRCA (time to most recent common ancestor) estimates against ancient individuals
- Geographic and temporal distribution of their closest ancient matches
- A narrative explanation of what the genetic data suggests about their deep ancestry

---

## Input Data

| File | Description |
|------|-------------|
| `data/input_data/modern_indvidual.txt` | Individual 1 — AncestryDNA V2.0 array, April 2025 |
| `data/input_data/modern_indv_2.txt` | Individual 2 — AncestryDNA V1.0 array, July 2025 |
| `data/input_data/v62.0_1240k_public.geno` | Ancient genotypes — binary EIGENSTRAT format |
| `data/input_data/v62.0_1240k_public.ind` | Ancient individual metadata (17,629 individuals) |
| `data/input_data/v62.0_1240k_public.snp` | SNP manifest (406,570 positions) |
| `data/input_data/v62.0_1240k_public.anno` | Extended annotations (culture, date, region) |

---

## Phase 1 — Core Analysis Engine (Individual 1 vs Ancient Dataset)

**Goal:** Implement lineage analysis methods matching the Deep Maniot study to compare Individual 1 against the full v62.0 ancient dataset and produce interpretable results. This phase is purely computational — no web interface.

### Step 1.1 — Data Parsing and Harmonisation

**Objective:** Convert Individual 1's AncestryDNA file and the EIGENSTRAT ancient dataset into a common representation for comparison.

**Tasks:**
- Parse `modern_indvidual.txt` (rsID, chromosome, position, allele1, allele2; build GRCh37)
- Parse `v62.0_1240k_public.snp` to extract SNP positions (EIGENSTRAT format: SNP ID, chromosome, genetic position, physical position, ref allele, alt allele)
- Parse `v62.0_1240k_public.ind` to build individual index (ID, sex, population label)
- Parse `v62.0_1240k_public.geno` (binary: 2 bits per genotype, row = SNP, column = individual) into a usable array or memory-mapped structure
- Strand-align modern alleles to ancient reference strand (flip complement where needed; handle ambiguous C/G, A/T SNPs by exclusion)
- Find overlapping SNPs between modern individual and ancient dataset by physical position (GRCh37)
- Record the intersection size — this is the working SNP set

**Expected outputs:**
- `output/step1/snp_overlap.tsv` — list of overlapping SNPs with modern + ancient allele coding
- `output/step1/modern_indv1_encoded.npy` — Individual 1 genotypes numerically encoded at overlap SNPs (0/1/2 dosage)

**Key considerations:**
- AncestryDNA V2.0 array has ~700k SNPs; expect ~80–150k overlap with 1240k panel
- Ancient samples are pseudo-haploid (one allele drawn randomly); modern samples are diploid — account for this in similarity calculations
- Exclude palindromic SNPs (A/T, C/G) that cannot be reliably strand-aligned without additional reference

---

### Step 1.2 — Haplogroup Assignment

**Objective:** Determine Individual 1's Y-DNA (paternal) and mtDNA (maternal) haplogroup, placing them on the ancient phylogenetic tree exactly as done in the Deep Maniot study.

**Tasks:**

**Y-DNA (paternal lineage):**
- Extract chromosome Y SNPs from Individual 1's AncestryDNA file
- Map Y SNPs to the ISOGG/YCC haplogroup tree using branch-defining variants
- Assign terminal haplogroup (e.g., R-M269, J-L26, E-V13)
- Search `v62.0_1240k_public.ind` and `.anno` for ancient individuals sharing haplogroup at same or upstream branches
- Compute phylogenetic distance from Individual 1 to each ancient Y-haplogroup match

**mtDNA (maternal lineage):**
- Extract mitochondrial SNPs from Individual 1's AncestryDNA file (chrMT)
- Assign mtDNA haplogroup using established mtDNA phylogeny (PhyloTree B17)
- Search ancient dataset annotations for matching mtDNA haplogroups
- Identify upstream and downstream branches in ancient samples

**Expected outputs:**
- `output/step2/ydna_haplogroup.json` — assigned Y haplogroup + confidence + branch-defining variants found
- `output/step2/mtdna_haplogroup.json` — assigned mtDNA haplogroup
- `output/step2/ancient_haplogroup_matches.tsv` — ancient individuals sharing Y or mt haplogroup with Individual 1, sorted by phylogenetic proximity
- `output/step2/haplogroup_report.md` — human-readable interpretation

**Reference tools:**
- ISOGG haplogroup trees (freely available, embed as JSON reference)
- PhyloTree Build 17 for mtDNA
- The Deep Maniot study used FamilyTreeDNA haplotree + Yfull for TMRCA — replicate TMRCA logic using mutation rate estimates

---

### Step 1.3 — Genome-wide SNP Similarity and Population Affinity

**Objective:** Beyond haplogroups, compute autosomal (whole-genome) similarity between Individual 1 and each ancient individual at the overlapping SNP set, identifying closest population-level matches.

**Method:** This mirrors the IBD and population affinity analysis used in the Deep Maniot study.

**Tasks:**

**Pairwise genetic distance:**
- For each of the 17,629 ancient individuals, compute a genetic distance score to Individual 1
- Use simple allele-sharing distance (ASD): proportion of overlapping SNPs where alleles differ, accounting for pseudo-haploidy of ancient samples
- Filter ancient individuals by minimum SNP overlap threshold (e.g., ≥10,000 SNPs) to ensure reliable comparisons

**Population-level summary:**
- Aggregate individual distances to population means (using population labels from `.ind` file)
- Rank populations by mean distance to Individual 1
- Identify top 20 closest ancient populations

**PCA projection:**
- Project Individual 1 into the PCA space of ancient individuals using the overlapping SNP set
- Use EIGENSOFT SmartPCA or equivalent Python implementation (scikit-allel / bed-reader)
- Visualise Individual 1's position relative to ancient population clusters on PC1–PC2 and PC1–PC3

**Expected outputs:**
- `output/step3/pairwise_distances.tsv` — distances to all 17,629 ancient individuals
- `output/step3/population_distances.tsv` — mean/median distance per population, ranked
- `output/step3/pca_coordinates.tsv` — PC loadings for ancient individuals + Individual 1 projected in
- `output/step3/top_matches_report.md` — narrative of top ancient population matches with dates and geographies

---

### Step 1.4 — TMRCA Estimation for Closest Haplogroup Matches

**Objective:** Estimate the time to most recent common ancestor between Individual 1 and their closest ancient haplogroup matches, as done in the Deep Maniot study using FTDNATiP-style mutation rate models.

**Tasks:**
- For Y-DNA: count branch-defining SNP differences between Individual 1 and each ancient Y-haplogroup match
- Apply Y-chromosome mutation rate (~0.74 × 10⁻⁹ substitutions/bp/year, or calibrate to known haplogroup ages)
- Estimate TMRCA with confidence intervals for top 10 Y-haplogroup matches
- For mtDNA: apply mtDNA mutation rate (~1.26 × 10⁻⁸ substitutions/bp/year) for maternal lineage TMRCA
- Cross-validate TMRCA estimates against known archaeological dates of ancient samples where available

**Expected outputs:**
- `output/step4/ydna_tmrca.tsv` — TMRCA estimates (point estimate + 95% CI) for top ancient Y matches
- `output/step4/mtdna_tmrca.tsv` — same for mtDNA
- `output/step4/tmrca_timeline.json` — structured data for timeline visualisation

---

### Step 1.5 — AMOVA and Population Structure

**Objective:** Replicate the AMOVA (Analysis of Molecular Variance) approach from the Deep Maniot study to quantify how much of Individual 1's genetic variance is explained by different ancient populations and regions.

**Tasks:**
- Select a panel of reference ancient populations (e.g., Anatolian Neolithic, European Bronze Age, Steppe, Caucasus, Levant, North Africa) from the `.anno` metadata
- Run AMOVA partitioning genetic variance: within-individual, among-individuals-within-population, among-populations
- Compute FST analogues between Individual 1 and each reference population group
- Identify which ancient population groupings minimise residual variance for Individual 1

**Expected outputs:**
- `output/step5/amova_results.tsv`
- `output/step5/fst_table.tsv`

---

### Step 1.6 — Interpretation and Report Generation

**Objective:** Synthesise outputs from steps 1.1–1.5 into a structured, human-readable report for Individual 1.

**Tasks:**
- Write a narrative interpretation: paternal lineage origin and migration story, maternal lineage origin, autosomal closest ancient populations, estimated timeframes
- Map ancient matches geographically and temporally (data for map + timeline visualisation)
- Validate results against known population genetics literature for consistency
- Flag any anomalies (e.g., very low SNP overlap with best matches, conflicting haplogroup vs autosomal signals)

**Expected outputs:**
- `output/step6/individual1_report.json` — full structured result (used as template for web platform)
- `output/step6/individual1_report.md` — human-readable version
- `output/step6/individual1_map_data.geojson` — geographic coordinates of ancient matches for map rendering

---

## Phase 2 — Web Platform Architecture

**Goal:** Design and build the web application that will serve the analysis to end users.

### Step 2.1 — Technology Stack Decision

**Recommended stack:**
- **Backend:** Python (FastAPI) — consistent with existing analysis codebase
- **Frontend:** React + TypeScript with Tailwind CSS
- **Database:** PostgreSQL (user accounts, upload metadata, cached results)
- **File storage:** S3-compatible (uploads + results)
- **Job queue:** Celery + Redis (async analysis jobs)
- **Hosting:** AWS or Azure (review existing Azure infrastructure already in repo)

### Step 2.2 — System Architecture Design

**Components to design:**
- File upload and validation service (accept AncestryDNA, 23andMe, FTDNA, MyHeritage raw formats)
- Async job runner that executes the Phase 1 analysis pipeline per user
- Result caching layer (results are deterministic; cache by file hash)
- User account system (email + password; optionally OAuth via Google/Apple)
- Privacy and consent model — critical for genetic data (see Phase 5)

### Step 2.3 — Database Schema

**Core tables:**
- `users` — account info, consent timestamps
- `uploads` — raw file reference, format, upload date, hash
- `analysis_jobs` — job ID, status, created/completed timestamps
- `results` — structured JSON results linked to job
- `haplogroup_assignments` — Y-DNA and mtDNA haplogroup per result

---

## Phase 3 — Analysis Pipeline Integration

**Goal:** Wrap the Phase 1 Python analysis into a production-ready pipeline that runs reliably on user uploads.

### Step 3.1 — Input Format Support

Extend the parser from Step 1.1 to handle all common consumer DNA formats:
- AncestryDNA (V1.0, V2.0 arrays) — already have examples
- 23andMe (v3, v4, v5 chips)
- FamilyTreeDNA (Family Finder)
- MyHeritage
- Living DNA

All formats use rsID-based SNP identification; strand-alignment and overlap logic is the same.

### Step 3.2 — Pipeline Optimisation

The v62.0 dataset is large (46 GB). For production:
- Pre-compute a population-level SNP matrix (mean allele frequencies per population per SNP) to reduce per-user computation from O(17,629 individuals) to O(~450 populations)
- Pre-compute PCA eigenvectors on the ancient dataset once; project each new user in without recomputing
- Cache haplogroup reference trees in memory
- Target analysis time: ≤5 minutes per user on a standard compute instance

### Step 3.3 — Quality Control Gates

Before running analysis, validate user input:
- Minimum SNP count (reject files with <100k genotyped SNPs)
- Expected chromosomes present (1–22, X; optionally Y and MT)
- SNP overlap with 1240k panel ≥ 5,000 (warn user if lower; reject if <1,000)
- Detect and reject non-human or clearly contaminated files

---

## Phase 4 — Visualisations and User Report

**Goal:** Build the interactive report UI that presents results to users.

### Step 4.1 — Report Components

| Component | Description |
|-----------|-------------|
| Haplogroup badge | Y-DNA and mtDNA haplogroup displayed prominently with phylogenetic context |
| Lineage timeline | Horizontal timeline showing when ancestral haplogroup branches emerged, tied to archaeological cultures |
| Ancient population bar chart | Top 20 closest ancient populations ranked by genetic distance |
| Geographic map | Interactive map with pins on ancient match locations (using geojson from Step 1.6) |
| PCA scatter plot | Interactive PCA showing user projected into ancient population space |
| TMRCA estimates | For top Y and mt haplogroup matches: date ranges with CI visualised |
| Narrative text | Plain-English interpretation of all results |

### Step 4.2 — Charting Libraries

- **Maps:** Mapbox GL JS or Leaflet with custom ancient-world base layer
- **PCA / scatter:** D3.js or Plotly.js
- **Timeline:** D3.js custom or vis-timeline
- **Charts:** Recharts (React-native, easy integration)

---

## Phase 5 — Privacy, Ethics, and Legal

**This phase is non-negotiable and must be completed before any user data is collected.**

### Step 5.1 — Privacy Framework

- Define data retention policy (how long raw uploads are stored)
- Allow users to delete their data at any time (GDPR Article 17)
- Encrypt raw genotype files at rest and in transit
- Do not share or sell genetic data under any circumstances
- Separate personal identity from genetic data in the database (pseudonymisation)

### Step 5.2 — Consent Flow

- Explicit informed consent before upload explaining: what data is collected, how it is used, where it is stored, who can access it, how to request deletion
- Minimum age requirement (18+, or parental consent for minors)
- Separate consent checkbox for any future research use (optional, off by default)

### Step 5.3 — Legal Review

- GDPR compliance (EU users)
- CCPA compliance (California users)
- HIPAA does not apply (consumer ancestry, not medical/diagnostic), but review guidance
- Terms of service drafted by legal counsel before launch

---

## Phase 6 — Launch and Scale

### Step 6.1 — Beta Testing
- Recruit 20–50 beta users with diverse ancestry to validate report quality
- Cross-check haplogroup assignments against users who already know their haplogroups (e.g., from 23andMe reports)
- Collect qualitative feedback on report clarity and accuracy

### Step 6.2 — Infrastructure Scaling
- Auto-scaling job workers based on queue depth
- CDN for static assets
- Rate limiting on uploads (prevent abuse)
- Cost estimate per analysis (compute + storage) to inform pricing model

### Step 6.3 — Dataset Updates
- v62.0_1240k_public is versioned; Allen Ancient DNA Resource releases new versions periodically
- Build a dataset versioning system so results are tied to dataset version
- Plan for reanalysis when major new datasets are released

---

## Immediate Next Action

**Start here:** `Phase 1 — Step 1.1` — parsing Individual 1's AncestryDNA file and the EIGENSTRAT ancient dataset, finding the SNP overlap, and encoding both into a common numerical format.

The output of Step 1.1 is the foundation for all subsequent analysis steps. Nothing else in Phase 1 can proceed without a clean, strand-aligned SNP overlap.

**Working files:**
- Modern: `data/input_data/modern_indvidual.txt` (Individual 1, AncestryDNA V2.0)
- Ancient genotypes: `data/input_data/v62.0_1240k_public.geno` (binary EIGENSTRAT)
- Ancient individuals: `data/input_data/v62.0_1240k_public.ind` (17,629 individuals)
- Ancient SNPs: `data/input_data/v62.0_1240k_public.snp` (406,570 SNPs)
- Ancient annotations: `data/input_data/v62.0_1240k_public.anno`

---

## Timeline Overview

| Phase | Description | Estimated Scope |
|-------|-------------|-----------------|
| Phase 1 | Core analysis engine — Individual 1 vs ancient dataset | ~3–4 weeks |
| Phase 2 | Web platform architecture design | ~1–2 weeks |
| Phase 3 | Production pipeline integration | ~2–3 weeks |
| Phase 4 | Visualisations and user report UI | ~3–4 weeks |
| Phase 5 | Privacy, ethics, legal | Parallel with Phases 2–4 |
| Phase 6 | Beta launch and scale | ~2–3 weeks |

---

*Last updated: April 2026*
