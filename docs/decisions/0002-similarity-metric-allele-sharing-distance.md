# ADR-0002: Use Allele-Sharing Distance (ASD) as the Genome-wide Similarity Metric

* Status: Accepted
* Date: 2026-04-21

## Context and Problem Statement

To rank ancient individuals and populations by genetic closeness to a modern individual, we need a pairwise similarity metric. Several metrics exist in population genetics. Which best handles the specific constraints of this project: pseudo-haploid ancient samples, diploid modern samples, variable SNP overlap across comparisons, and the need to run against 17,629 individuals at ~394k SNPs in reasonable compute time?

## Decision Drivers

* Must handle pseudo-haploid ancient samples correctly (no diploid assumption)
* Must be robust to variable SNP overlap (different ancient individuals have different coverage)
* Must be interpretable and explainable to a non-technical user
* Must be computationally tractable at the scale of 17,629 individuals × 394k SNPs
* Must align with methodology used in the Deep Maniot reference study

## Considered Options

* Allele-sharing distance (ASD) — mean absolute allele frequency difference per SNP
* Identity-by-state (IBS) proportion — proportion of shared alleles
* FST (fixation index) — population-level differentiation statistic
* IBD (identity-by-descent) — shared haplotype segments, infers recent common ancestry
* D-statistics / ADMIXTURE — admixture modelling, not a distance metric

## Decision Outcome

Chosen option: **Allele-sharing distance (ASD)** with pseudo-haploidisation of the modern individual, because it is simple, interpretable, handles pseudo-haploid ancient samples directly, scales linearly with SNP count, and is the metric used in the Deep Maniot reference study. Pseudo-haploidisation (10 random draws averaged) resolves the diploid/pseudo-haploid mismatch without introducing systematic bias.

### Positive Consequences

* Scales to the full 17,629 × 394k comparison in a single pass with chunked processing
* No assumptions about population structure, allele frequency distribution, or diploid genotypes
* Naturally handles variable SNP overlap via per-comparison denominators
* Interpretable: value of 0 = identical, 1 = completely different; typical same-population
  distances are 0.12–0.20
* Population-level means are reliable aggregates even when individual-level estimates are noisy
  (due to pseudo-haploid sampling noise in ancient samples)

### Negative Consequences

* Does not detect IBD segments (shared haplotype blocks) — insensitive to very recent shared ancestry
* Population-level means can be dominated by outlier ancient individuals within a population
* Pseudo-haploidisation introduces a small random component; averaged across 10 draws at 394k SNPs,
  variance is negligible but technically not deterministic (seeded for reproducibility)
* Does not directly produce admixture proportions — cannot say "X% Steppe ancestry"

## Pros and Cons of the Options

### ASD with pseudo-haploidisation (chosen)

* Good, because directly comparable to the Deep Maniot reference study
* Good, because handles pseudo-haploid ancient data without requiring imputation
* Good, because O(n_indiv × n_snps) compute — tractable in pure NumPy with chunking
* Good, because interpretable distance range (0–1)
* Bad, because no haplotype phasing — misses IBD segment signals
* Bad, because not directly interpretable as admixture proportions

### IBS proportion

* Good, because intuitive (proportion of matching alleles)
* Bad, because same systematic bias as naive ASD at heterozygous diploid sites
  unless pseudo-haploidisation is applied anyway
* Bad, because not the standard in ancient DNA literature

### FST

* Good, because standard population differentiation measure
* Bad, because requires defining population groups in advance — cannot compute individual-level
  comparisons to a modern individual directly
* Bad, because assumes Hardy-Weinberg equilibrium, which is violated in ancient datasets

### IBD detection

* Good, because detects recent genealogical relatedness (segments > 7 cM ≈ < 10 generations)
* Bad, because ancient samples are pseudo-haploid and not phased — IBD segment detection
  requires phased haplotypes, making it inapplicable here
* Bad, because computationally intensive (requires pairwise haplotype comparison)

### ADMIXTURE / qpAdm

* Good, because produces interpretable ancestry proportion estimates
* Bad, because not a distance metric — cannot rank 17,629 individuals directly
* Bad, because requires a reference panel design and is sensitive to panel choice
* Bad, because runtime is orders of magnitude higher
