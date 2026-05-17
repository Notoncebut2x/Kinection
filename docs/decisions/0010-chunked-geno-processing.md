# ADR-0010: Process EIGENSTRAT GENO File in SNP-row Chunks

* Status: Accepted
* Date: 2026-04-22

## Context and Problem Statement

The AADR v62.0 GENO binary file encodes genotypes for 17,629 individuals across 1,233,013 SNP positions (406,570 on the 1240k capture panel subset). Loading the full dataset into memory at once would require approximately 17,629 × 406,570 × 1 byte ≈ 7.2 GB, which exceeds available RAM on a typical development machine and would be a bottleneck for production workers processing one job at a time. How should the GENO file be read during ASD computation?

## Decision Drivers

* Fit within typical machine RAM (16 GB development, production workers may be 4–8 GB)
* Maintain acceptable runtime (target: full ASD computation in < 10 minutes locally)
* Support the pseudo-haploid draw loop without rereading the file repeatedly
* Allow progress logging (know what proportion of computation is done)

## Considered Options

* Load full GENO matrix into RAM — read once, compute in full matrix operations
* Memory-map the GENO file — OS manages paging, access as if in RAM
* Chunk processing — read N SNP rows at a time, accumulate ASD statistics
* Pre-filter GENO to overlap SNPs only — write a reduced file before analysis

## Decision Outcome

Chosen option: **Chunk processing with chunks of 5,000 SNPs**, combined with a memory-mapped file reader for random-access SNP row reads. Within each chunk, the 5,000 × 17,629 submatrix is loaded into RAM (< 90 MB), pseudo-haploid draws are run across the chunk, and partial ASD sums are accumulated into a single (17,629,) accumulator array. This keeps peak memory usage well under 1 GB while enabling progress logging every 10 chunks.

### Positive Consequences

* Peak RAM usage is O(chunk_size × n_indiv) ≈ 5,000 × 17,629 × 4 bytes ≈ 350 MB —
  well within a production worker's budget
* Progress is loggable at chunk boundaries — the user-facing job status can report
  percentage completion
* Pseudo-haploid draws run per-chunk against the same loaded submatrix —
  no rereading for additional draws within a chunk
* Chunk size is tunable without code changes (CHUNK_SIZE constant)

### Negative Consequences

* Sequential I/O — chunks are read in GENO row order; random access would require
  seeking across a large file, which is slower than sequential reads
* Chunks must be sorted by GENO row index — the overlap SNP list must be pre-sorted
  (enforced by Step 1.1 which sorts by geno_index)
* Total runtime scales linearly with n_chunks × N_PSEUDO_DRAWS — 79 chunks × 10 draws
  = 790 chunk iterations; each is fast but adds up to several minutes total

## Pros and Cons of the Options

### Chunk processing with 5,000-SNP chunks (chosen)

* Good, because bounded memory — safe on production workers
* Good, because progress logging possible
* Good, because pseudo-haploid draws reuse chunk data (no extra reads per draw)
* Bad, because sequential — cannot parallelise across chunks without shared accumulators
* Bad, because total runtime depends on chunk count × draws

### Full matrix load into RAM

* Good, because single read, then all operations are in-memory matrix math
* Good, because potentially faster for large-RAM machines (vectorised over all SNPs at once)
* Bad, because 7+ GB RAM requirement — not feasible for standard production instances
* Bad, because fails entirely on low-RAM machines with no graceful degradation

### Memory-mapped file

* Good, because the OS handles paging — no explicit chunk management
* Good, because random access without seek overhead (OS manages the cache)
* Bad, because on macOS/Linux, mmap of a 7 GB file still touches all pages eventually —
  equivalent to loading fully if all SNPs are accessed
* Bad, because memory pressure causes thrashing when working set exceeds physical RAM
* Note: GenoFile uses mmap for random row access, but chunk iteration is still used
  to control how many rows are active in memory at once

### Pre-filter GENO to overlap SNPs

* Good, because reduces the file to ~32% of its original size (394k of 1.23M SNPs)
* Good, because subsequent analysis would be faster with fewer SNPs to read
* Bad, because pre-filtering requires reading the full file once anyway — no net savings
  for a one-time analysis
* Bad, because adds a preprocessing step that complicates the pipeline for production
  (per-user SNP sets differ; the pre-filtered file cannot be shared across users)
