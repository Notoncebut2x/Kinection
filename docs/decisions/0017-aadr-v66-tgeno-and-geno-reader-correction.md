# ADR-0017: AADR v66 `tgeno` format and geno reader byte-level correction

* Status: Accepted
* Date: 2026-07-04
* Amends: 0014 (encoding convention — value inversion; this ADR adds the byte-level decoding)
* Phase: 1/2 (Analysis engine + cloud pipeline)

## Context

Moving the cloud pipeline onto AADR **v66** (the current release) surfaced two independent problems, both discovered during the first real end-to-end upload runs.

### 1. v66 ships the 1240K genotypes in the new `tgeno` (transpose_packed) format

AADR v66+ distributes the 1240K `.geno` in a **transposed** packed layout (magic `TGENO`, one record **per individual**), documented as `transpose_packed` in the AdmixTools `convertf` README. The Kinection pipeline reads **SNP rows** (all individuals at one SNP); in the transposed layout a single SNP row is scattered across the entire 7 GB file, which is unworkable over R2 range requests (and pathological even on a local mmap).

Confirmed empirically: the R2 object header was `TGENO 23250 1233013 …` and the file size matched `48 + n_indiv × ⌈n_snp/4⌉` exactly (individual-major), not the SNP-major layout.

### 2. The geno reader was mis-decoding at the byte level

Separately — and orthogonal to ADR 0014's *value* inversion — the reader (`GenoFile`, `R2GenoFile`) decoded PACKEDANCESTRYMAP incorrectly on two counts:

* **Header size.** EIGENSOFT writes the header as one full record padded to `max(48, ⌈n_indiv/4⌉)` bytes — **4408** bytes for v62, **5813** for the converted v66 — not a fixed **64**. The reader's 64-byte assumption misaligned every SNP row.
* **Bit packing.** Four 2-bit genotypes per byte are packed **MSB-first** (first individual in bits 6–7: `(byte >> (2*(3-(i%4)))) & 3`), not LSB-first (`byte & 3`).

Verified against the AdmixTools `mcio.c` source and empirically: with the correct params, v62 and v66 decode to **identical** genotypes for shared individuals (100% on diploids); with the old params, ~50% (barely above chance).

ADR 0014 fixed the *allele-count* semantics (`geno=0` → hom allele2) and reported "sensible European results", but that was on a still-byte-misaligned reader — so autosomal outputs (ASD/PCA/admixture) before this fix were unreliable. (Y/mtDNA haplogroups were unaffected: they derive from the modern file + `.anno` text, never the packed geno.)

### 3. `.anno` header drift

v66 renamed many `.anno` columns (`Lat.`→`Latitude`, the Genetic-ID header text, Y-haplogroup columns gained version numbers) and changed assessment values (`PASS`→`Pass`, plus `PROVISIONAL_PASS`/`MERGE_PASS`). The exact-header parser returned **zero** records on v66, silently disabling PCA, admixture source selection, and all match metadata.

## Decisions

1. **Store the AADR geno in standard SNP-major PACKEDANCESTRYMAP in R2.** `scripts/convert_tgeno_to_packed.py` converts `tgeno` → standard packed (MSB, header `max(48, ⌈n_indiv/4⌉)`), with a round-trip self-test. The v66 R2 object was converted and re-uploaded (verified 100% round-trip). `update_aadr.py` should run this on future ingests (open item).

2. **Correct the geno readers.** `GenoFile`/`R2GenoFile` now compute `header_size = max(48, bytes_per_snp)` and unpack MSB-first. This is the byte-level counterpart to ADR 0014's value inversion; both are required for correct decoding. The value inversion (`geno_to_freq = [1.0, 0.5, 0.0, nan]`) stays downstream in steps 1.3/1.5 unchanged.

3. **Parse `.anno` by distinctive substrings, not exact headers.** `parse_anno_file` matches columns by robust substrings (and prefix for lat/lon), and normalises assessment to canonical `PASS`/`QUESTIONABLE`/`CRITICAL`. Backward-compatible with v62.

4. **Admixture source `group_id` patterns updated for v66** (EHG, EEF renamed), keeping v62 patterns for local runs.

## Consequences

* v66 works end-to-end in the cloud pipeline; a real upload now yields a correct, person-specific report.
* Any autosomal result produced before this fix must be re-run.
* Reader output is now consistent between v62 (local) and v66 (R2).
* The converted R2 geno is standard EIGENSOFT packed and re-readable by external tools.
* Open: teach `update_aadr.py` to auto-convert `tgeno` on ingest so future releases don't reintroduce the transposed-format problem.
