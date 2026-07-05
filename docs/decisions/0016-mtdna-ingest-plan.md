# ADR-0016: mtDNA capture-data ingest

* Status: Accepted (implemented 2026-05-25 in commit 6295598)
* Date: 2026-05-25
* Supersedes: —
* Phase: 1 (extends Steps 1.1, 1.2, 1.4)

## Context

mtDNA TMRCA is the natural sibling of the Y-DNA TMRCA in Step 1.4 — it gives the maternal-line analogue ("your mother's mother's … mother's MRCA with this ancient woman lived ~X years ago"). Step 1.4 currently skips mtDNA because the 1240k Y panel — the only AADR file the pipeline ingests — contains **zero mtDNA positions**.

Survey of the AADR v66 release confirms:

| File | Size | Purpose |
|---|---|---|
| `v66.1240K.aadr.PUB.{geno,ind,snp,anno}` | ~7.2 GB | Autosomal + Y + (no mt). What we already use. |
| `v66.2M.aadr.PUB.*` | ~12.5 GB | Denser 2M-SNP panel. Same composition. |
| `v66.HO.*`, `v66.compatibility_HO.*` | ~5 GB combined | Affymetrix Human Origins genotype panel. |
| **`v66.MT.repo.fa.gz`** | **5.6 MB** | **Per-individual mt-genome consensus FASTA — what we want.** |
| `aadr_v66.0__README_MT.docx` | <1 MB | README for the mt repo |
| `mtdna_uncompress_v66.py` | <1 MB | Helper for parsing the FASTA |

Per the AADR paper, the mt repository has held ~4,122 ancient mt genomes since v52.2 (Mallick et al., *Scientific Data* 2024). Format: **FASTA** of consensus sequences aligned to rCRS (revised Cambridge Reference Sequence, 16,569 bp).

The modern AncestryDNA file covers **~195 mt positions** (encoded as chromosome `26` in the AncestryDNA convention). Empirically confirmed on a test AncestryDNA file.

## Decision (Proposed)

Implement mtDNA TMRCA as a follow-on step. Concrete plan in 5 sub-steps:

### 1. mtDNA data ingest

- Add `v66.MT.repo.fa.gz` to the R2 manifest at `dataset/v66/mt.fa.gz`.
- Extend `scripts/update_aadr.py` to upload the mt FASTA alongside the four PACKGENO files.
- Local mode: extend `scripts/download_v62_local.py` (or new `download_aadr_local.py`) to fetch the mt file too — it's tiny (5.6 MB).

### 2. mtDNA FASTA parser

- New module `scripts/utils/mt_fasta.py` exposing:
  - `parse_mt_repo(path: Path) -> dict[str, str]` returning `{genetic_id: 16569-char consensus}`.
  - Optional `MtRepo.lookup(gid, position) -> str | None` for sparse lookups.
- Handles gzip transparently; tolerates `N` / gap characters; validates length == 16569.

### 3. Modern mt extraction

- Extend `scripts/step1_parse_harmonise.py` to write a separate `modern_mt_positions.tsv` with `rsid \t position \t allele1 \t allele2` for chromosome=26 entries from the AncestryDNA file.
- Outputs to `output/step1_<label>/modern_mt_positions.tsv`.
- Mapping AncestryDNA mt positions (build 37 / rCRS) to rCRS indices is 1:1 — but verify against a known mt position table because some consumer arrays use different mt reference conventions.

### 4. Extend Step 1.4 to compute mt TMRCA

- Add `--lineage {y,mt,both}` flag or default to both when data available.
- For mt, iterate over top mt-haplogroup matches from step 2's `ancient_haplogroup_matches.tsv` (filter `match_type` containing `MT`).
- For each match: extract the ancient's mt consensus from the FASTA, intersect with the modern's ~195 sampled positions, count differences `k` on `L` covered positions.
- TMRCA = k / (2 * μ_mt * L) with **μ_mt = 1.665 × 10⁻⁸ per bp per year** (Soares et al. 2009 calibration of the mt control region).
- Same Poisson-exact CI on k as the Y path.
- **No ascertainment bias issue** — unlike the 1240k Y panel, mt array positions are sparsely chosen but not enriched for polymorphism in the same systematic way, AND we have *full* mt genomes for the ancient side. The standard per-bp formula works here. This is a much cleaner TMRCA estimate than the Y one.

### 5. Tests + report integration

- Unit test for `parse_mt_repo` (synthetic FASTA fixture).
- Unit test for the modern mt extraction (using a small synthetic AncestryDNA slice with the rsIDs allow-listed via `# allow-raw-dna`).
- Update `run_local.py` to render an mt-TMRCA section in the combined report.
- Update `scripts/step1_6_synthesis.py` to include mt-TMRCA matches in `report.json` and `map_data.geojson` (with `match_type: "mt_tmrca"`).

## Consequences

**Positive:**
- Closes the only known gap in Step 1.4 (mtDNA was deferred with caveat).
- Adds a much **more reliable** TMRCA estimate to the report — mt mutation rate is well-calibrated and we have full mt genomes on the ancient side, so the rate-calibration uncertainty that plagues Y TMRCA (ADR 0014, Decision 2) does not apply.
- The mt FASTA is tiny (5.6 MB), so no storage/bandwidth concerns.

**Negative:**
- Adds a new parser and data format that the rest of the pipeline must learn about.
- Modern coverage is sparse (~195 positions out of 16,569) so per-comparison CIs will still be wide — but the formula is honest, not orders-of-magnitude wrong.
- Requires validating that the AncestryDNA mt position numbering matches rCRS exactly. Empirical check needed before trusting outputs.

## Alternatives considered

**A. Skip mt entirely and rely on Y-TMRCA only.**
Current state. Loses the maternal-line narrative which is half the report's appeal.

**B. Use a different mt reference dataset (e.g. MITOMAP, custom Phylotree integration).**
Larger and less curated than AADR's mt repo, and would force a two-source comparison the rest of the pipeline doesn't do. AADR's mt data is already date-annotated and population-labelled — it's the right primitive.

**C. Ingest 23andMe / FTDNA mt data for modern instead of AncestryDNA.**
Other consumer formats sometimes have denser mt coverage (e.g. FTDNA mt-Full sequences cover all 16,569 bp). Not needed for v1 — AncestryDNA's 195 positions gives at least *some* signal, and dense-modern-mt is a Phase 3 input-format expansion (workplan Step 3.1).

## Status

**Accepted, implemented.** Shipped in commit `6295598` (2026-05-25). Actual scope landed at ~700 LOC across `scripts/utils/mt_fasta.py` (new, ~110 LOC), `scripts/data/mt_rcrs.txt` (data), `tests/test_mt_fasta.py` (24 tests), and edits to `step1_parse_harmonise.py`, `step1_4_tmrca.py`, `step1_6_synthesis.py`.

One implementation correction worth noting: the original plan said "no ascertainment correction needed" for mt because we have full ancient mt genomes. That was wrong on the modern side — AncestryDNA's ~190 mt positions ARE ascertained for polymorphism. The implementation uses the full mt-genome length (16,569 bp) as the per-bp formula's denominator, with the assumption that unsampled positions agree between the two samples — a reasonable assumption given mt's high conservation outside the few hundred well-known polymorphic sites. Without this correction, TMRCAs came out ~100× too high (millions of years for k=3 over 170 sampled sites).

Validated on a test individual: 5 mt-haplogroup matches, TMRCAs in the 5,400-year range with Poisson 95% CIs of 1,100–16,000 — biologically sensible for related R-haplogroup lineages.

## References

- Mallick et al., *The Allen Ancient DNA Resource (AADR), a curated compendium of ancient human genomes*, Scientific Data 2024 — https://www.nature.com/articles/s41597-024-03031-7
- Soares et al. 2009 — mt mutation rate calibration: https://doi.org/10.1016/j.ajhg.2009.05.001
- AADR Dataverse listing (queried 2026-05-25): 25 files, including `v66.MT.repo.fa.gz` and the README/uncompress helper.
