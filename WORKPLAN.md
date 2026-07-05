# DNA Lineage Platform — Workplan

**Project:** Web platform for personal DNA comparison against ancient population datasets
**Reference datasets:** Allen Ancient DNA Resource
  - v62.0 — cached locally (`data/input_data/`), used for offline runs
  - v66 — uploaded to Cloudflare R2 (`dataset/v66/`, 7.17 GB .geno), used for cloud runs
**Pipeline version source:** R2 manifest at `dataset/current_version.json` (read at startup; restart to pick up new versions)
**Methodology reference:** [Uniparental analysis of Deep Maniot Greeks, *Communications Biology* 2026](https://www.nature.com/articles/s42003-026-09597-9)
**Started:** April 2026
**Last updated:** May 2026

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

**Modern individuals (gitignored, local-only, never uploaded to cloud):**

| File | Description |
|------|-------------|
| `data/input_data/AncestryDNA_rn.txt` | Individual "rn" — AncestryDNA raw |
| `data/input_data/AncestryDNA_jn.txt` | Individual "jn" — AncestryDNA raw |

**Ancient reference dataset (v62 on disk, v66 in R2):**

| File | Source | Description |
|------|--------|-------------|
| `v62.0_1240k_public.{geno,ind,snp,anno}` | `data/input_data/` | v62 — 17,629 individuals, ~5.5 GB |
| `v66.1240K.aadr.PUB.{geno,ind,snp,anno}` | R2 `dataset/v66/` | v66 — 7.17 GB .geno, current cloud version |

`run_local.py` auto-detects local AADR files and uses them offline; falls back to R2 when missing. `update_aadr.py` manages R2.

---

## Phase 1 — Core Analysis Engine (Individual 1 vs Ancient Dataset)

**Goal:** Implement lineage analysis methods matching the Deep Maniot study to compare Individual 1 against the full v62.0 ancient dataset and produce interpretable results. This phase is purely computational — no web interface.

**Status overview:**

| Step | Title | Status |
|------|-------|--------|
| 1.1 | Data Parsing & Harmonisation     | ✅ done — `scripts/step1_parse_harmonise.py` |
| 1.2 | Haplogroup Assignment            | ✅ done — `scripts/step2_haplogroup.py` |
| 1.3 | Genome-wide Similarity & PCA     | ✅ done — `scripts/step3_similarity_pca.py` (encoding fix: ADR 0014) |
| 1.4 | TMRCA Estimation                 | ✅ done — `scripts/step1_4_tmrca.py` (ADR 0014; Y-only, mt deferred) |
| 1.5 | Admixture Decomposition (NNLS)   | ✅ done — `scripts/step1_5_admixture.py` (ADR 0013; replaces the original AMOVA scope) |
| 1.6 | Interpretation & Report          | ✅ done — Markdown via `run_local.py`; structured `report.json` + `map_data.geojson` via `scripts/step1_6_synthesis.py` |

**Phase 1 complete.** Next: Phase 2 (web platform), and the production-readiness items below before any external user touches the system.

### Step 1.1 — Data Parsing and Harmonisation

**Status:** ✅ Implemented in `scripts/step1_parse_harmonise.py`. Outputs land in `output/step1_<label>/` (per-individual).

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

**Status:** ✅ Implemented in `scripts/step2_haplogroup.py`. Outputs land in `output/step2_<label>/`.

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

**Status:** ✅ Implemented in `scripts/step3_similarity_pca.py`. Outputs land in `output/step3_<label>/`.

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

### Step 1.5 — Admixture Decomposition (NNLS)  ✅ DONE

**Objective:** Decompose the modern individual's autosomal ancestry into proportions of six ancient source populations (WHG, EHG, EEF, Steppe, Levant_N, Iran_N) using constrained non-negative least squares against population mean allele frequencies. See ADR 0013 for rationale (chose NNLS over qpAdm to avoid f-statistic complexity for v1).

**Implementation:** `scripts/step1_5_admixture.py`
**Outputs:** `output/step1_5_<label>/admixture_decomposition.json`, `admixture_report.md`, `source_coverage.tsv`

**Future (not committed):** Add AMOVA / FST as Step 1.7 if needed for the report narrative.

---

### Step 1.6 — Interpretation and Report Generation

**Status:** ✅ Implemented in `scripts/step1_6_synthesis.py`. Outputs land in `output/step1_6_<label>/`: `report.json` (consolidated structured doc, schema_version 1.0) and `map_data.geojson` (FeatureCollection of top autosomal + Y-TMRCA matches). Markdown synthesis stays in `run_local.py`'s `build_report`.

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

**Decided stack (see ADRs 0011, 0012; supersedes 0007–0009):**
- **API:** Cloudflare Workers (TypeScript) — `workers/api/`
- **Analysis runner:** Local polling daemon (`scripts/daemon.py`) executing the Python pipeline. Replaced Celery+Redis (commit 7671584).
- **Object storage:** Cloudflare R2 — AADR reference + per-job outputs (modern raw files NEVER stored here, per commit 9adabaf)
- **Database:** Cloudflare D1 (SQLite) — job state, deletion receipts, user accounts (when added)
- **KV:** Workers KV — caches the AADR version manifest
- **Frontend:** React + TypeScript + Tailwind — not yet built
- **ADRs to retire / mark superseded:** 0007 (FastAPI), 0008 (Celery+Redis), 0009 (PostgreSQL)

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
- AncestryDNA (V1.0, V2.0 arrays) — ✅ `parse_ancestry_dna`
- 23andMe (v3, v4, v5 chips) — ✅ `parse_23andme` (single genotype column; letter chroms; auto-detected)
- FamilyTreeDNA (Family Finder) — pending
- MyHeritage — pending
- Living DNA — pending

Format is auto-detected by `parse_modern_dna(path, fmt="auto")` (vendor header signature, else 4-vs-5 column count); override via the `MODERN_DNA_FORMAT` env var. All formats collapse to the shared `dict[str, SNP]` representation, so strand-alignment and overlap logic downstream is identical. Covered by `tests/test_parsers_23andme.py`.

### Step 3.2 — Pipeline Optimisation

The AADR .geno is ~5–7 GB (v62: 5.4 GB, v66: 7.2 GB). For production:
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

### Step 5.1.1 — Modern DNA Lifecycle: Upload, Analyse, Permanently Delete

Treat the raw modern-individual file as the most sensitive object in the system. It must exist in the cloud only for the duration of analysis, and never appear in any log, transcript, or version-control artefact.

**Upload (client → R2):**
- Browser uploads directly to a per-job R2 key (e.g. `uploads/<job_id>/raw.txt`) via a short-lived presigned PUT URL minted by the Worker — the file never traverses the Worker, so Worker logs cannot capture it.
- Presigned URL expires in ≤15 minutes; one-shot use enforced by job-state check.
- Object is written with SSE-C or server-side encryption with a per-job key derived from the user's session (so a bucket-wide credential leak does not yield plaintext).
- Set `Cache-Control: no-store` and a short R2 object-lock TTL as a defence-in-depth backstop in case post-analysis deletion fails.

**During analysis:**
- Worker / analysis runner reads the object via short-lived signed GETs scoped to that single key. No copies are made to other buckets or persistent local disk; if a runner needs disk, it uses an ephemeral tmpfs that is wiped on container teardown.
- Derived artefacts (encoded genotypes, overlap tables) are stored under `outputs/<job_id>/` separately, so the raw file can be deleted independently of results.
- Analysis logs must redact any line containing rsID-level genotype calls. Add a structured-logging filter that drops fields named `genotype`, `raw_line`, `allele1`, `allele2`, and similar; unit-test the filter against a sample raw line.

**Post-analysis deletion (must run, must be verified):**
- On successful analysis completion, issue `DeleteObject` against the raw key, then `HeadObject` to confirm 404. Persist a deletion receipt (timestamp, key, requestor) to D1.
- On failure, retry deletion with exponential backoff; if it still fails after N retries, page on-call — never leave the file behind silently.
- Run a daily reaper job that lists everything under `uploads/` older than the max-analysis-window (e.g. 24h) and force-deletes it. This is the backstop for crashed jobs that skipped the inline delete.
- Allow users to trigger immediate deletion from the UI at any time, even mid-analysis (aborts the job).
- Document and test that R2 does not retain object versions for the `uploads/` prefix (versioning OFF, or lifecycle rule purges noncurrent versions within 24h). A delete that leaves a recoverable version is not a delete.

**Prevent leaks into git, chat, and logs:**
- `.gitignore` already excludes `data/input_data/` and `output/` — keep it that way, and add a pre-commit hook that greps staged diffs for AncestryDNA-header signatures (`# AncestryDNA raw data`, `rsid\tchromosome\tposition\tallele1\tallele2`) and `rs\d+\t\d+\t\d+\t[ACGT0]\t[ACGT0]` patterns. Block the commit on match.  <!-- allow-raw-dna -->
- CI runs the same scanner on every PR.
- Worker / runner code must never log raw request bodies. Add a lint rule or code-review checklist item.
- When users (or developers) paste raw DNA into a chat with an AI assistant, the assistant should refuse to echo it back, and the file-handling code must never include raw genotype lines in error messages, exception traces, or telemetry payloads. Sanitise exceptions at the boundary.
- Periodic audit: run `git log --all -p -S "AncestryDNA"` and a pickaxe search for `rs\d+\t\d+\t\d+` against history; if anything shows up, rewrite history with `git filter-repo` and force-push (coordinate with team).

**Acceptance test for this step:**
- End-to-end test: upload a synthetic raw file, run analysis to completion, assert `HeadObject` on the raw key returns 404, assert no log line in the job's log stream contains the synthetic file's known rsIDs, assert git history is clean.

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

**Shipped:** full pipeline (steps 1.1–1.6), Cloudflare backend (Worker + R2 + D1 + KV), presigned-PUT upload + verified deletion, React frontend on Pages, AADR **v66** (tgeno→packed + geno reader fix, ADR 0017), AncestryDNA + 23andMe input, and validated end-to-end (a real upload → correct person-specific report).

What's left, in priority order:

1. **Y-DNA SNP→tree resolution (ISOGG map).** Ancient Y recorded in SNP notation (`R-M269`) is crudely stripped to `R`, so those individuals are *missed* by the haplogroup matcher (it only reads branch notation like `R1b1a1b`). Load an ISOGG SNP→haplogroup map so `M269 → R1b1a1b`, giving richer/correct paternal matches — and improving the modern Y caller's resolution too. (Larger task: thousands of SNPs.)
2. **Wire PCA into the report.** Step 3 computes PCA, but `report.json` (`pca` field) and the frontend scatter plot aren't hooked up — currently a placeholder.
3. **Auto-convert AADR `tgeno` on ingest.** Teach `update_aadr.py` to run `convert_tgeno_to_packed.py` so future AADR releases don't reintroduce the transposed-format problem (ADR 0017).
4. **Reaper Cron Trigger.** `scripts/reaper.py` exists; wire a Cloudflare Cron Trigger (or a cron on the daemon host) to sweep orphaned uploads daily.
5. **More input formats.** Extend the parser to FamilyTreeDNA, MyHeritage, Living DNA (AncestryDNA + 23andMe done).
6. **Move compute off the local daemon.** Containerise and run on Cloudflare Containers or a scale-to-zero runner (Fly Machines / Cloud Run) so no local machine is required.

**Working files:**
- Modern: `data/input_data/AncestryDNA_{rn,jn}.txt` (gitignored)
- Ancient (local): `data/input_data/v62.0_1240k_public.{geno,ind,snp,anno}`
- Ancient (R2):   `dataset/v66/v66.1240K.aadr.PUB.{geno,ind,snp,anno}`
- Entry point:    `python scripts/run_local.py --dna <file> --label <name>`

---

## Timeline Overview

| Phase | Description | Estimated Scope |
|-------|-------------|-----------------|
| Phase 1 | Core analysis engine — Individual 1 vs ancient dataset | ✅ complete (6/6 steps done) |
| Phase 2 | Web platform architecture design | ~1–2 weeks |
| Phase 3 | Production pipeline integration | ~2–3 weeks |
| Phase 4 | Visualisations and user report UI | ~3–4 weeks |
| Phase 5 | Privacy, ethics, legal | Parallel with Phases 2–4 |
| Phase 6 | Beta launch and scale | ~2–3 weeks |

---

*Last updated: May 2026*
