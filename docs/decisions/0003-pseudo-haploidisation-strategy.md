# ADR-0003: Pseudo-haploidise the Modern Individual for ASD Computation

* Status: Accepted
* Date: 2026-04-22

## Context and Problem Statement

Modern individuals from consumer DNA tests are diploid — each autosomal SNP has two alleles (one from each parent), giving dosage values of 0 (hom ref), 1 (het), or 2 (hom alt). Ancient samples in the AADR are pseudo-haploid — only one allele is drawn per position during sequencing. When computing allele-sharing distance between a diploid modern sample and a pseudo-haploid ancient sample, heterozygous sites in the modern sample produce a frequency of 0.5, which is never a possible value on the ancient side (0 or 1 only). This creates a systematic upward bias in distances at all heterozygous sites regardless of true genetic closeness.

How should the diploid modern genotype be treated when comparing against pseudo-haploid ancient samples?

## Decision Drivers

* Eliminate systematic distance bias from the diploid/pseudo-haploid mismatch
* Maintain comparability with the Deep Maniot methodology
* Keep the computation deterministic enough to reproduce
* Avoid introducing more noise than necessary — the ancient pseudo-haploid sampling
  is already a source of noise; we should not compound it

## Considered Options

* Treat het sites as 0.5 frequency (no adjustment — naive approach)
* Pseudo-haploidise: randomly pick one allele at each het site for each draw,
  average across N draws
* Impute het sites as the population mean allele frequency
* Discard all heterozygous modern sites

## Decision Outcome

Chosen option: **Pseudo-haploidise with 10 random draws**, seeded for reproducibility (numpy RNG seed=42). At each heterozygous site, randomly pick 0 or 1. Run 10 independent draws and average the resulting ASD scores. This mirrors the standard approach in ancient DNA literature and is statistically stable at 394k autosomal SNPs (standard error of the mean across draws is negligible).

### Positive Consequences

* Removes the systematic 0.5 bias — het sites are treated as 0 or 1 on both sides,
  matching the ancient pseudo-haploid encoding
* 10 draws at 394k SNPs converges: SD of draw-to-draw ASD variation is < 0.001,
  which is below the noise floor of biological interest
* Reproducible: fixed RNG seed (42) produces the same result on every run
* Matches published ancient DNA methodology

### Negative Consequences

* Slight additional computational cost (10 passes per chunk vs 1) — acceptable since
  chunk-level GENO data is read once and reused across draws
* Results are technically seed-dependent; a different seed produces slightly different
  individual-level distances (though population-level means are stable)
* Does not recover the information in the heterozygous call — we lose the signal
  that a site is het, treating it as randomly 0 or 1

## Pros and Cons of the Options

### Pseudo-haploidise with N draws (chosen)

* Good, because eliminates systematic bias at het sites
* Good, because standard approach in ancient DNA pairwise comparison
* Good, because stable with N=10 at 394k SNPs — no need for more draws
* Bad, because adds slight non-determinism (mitigated by fixed seed)
* Bad, because discards the information that a site is heterozygous

### Treat het as 0.5 (naive, no adjustment)

* Good, because no additional computation
* Bad, because systematic inflation of distances at all het sites — a modern individual
  who is 30% het will have all those sites at 0.5, never matching the ancient 0 or 1
* Bad, because distances become dominated by het rate rather than ancestry signal
* Bad, because not used in academic ancient DNA literature for this reason

### Impute het as population mean frequency

* Good, because more statistically principled (expectation under HWE)
* Bad, because requires knowing the reference population allele frequencies in advance,
  creating a circular dependency (we're trying to find which population matches)
* Bad, because not the standard in the literature

### Discard het sites

* Good, because avoids the mismatch entirely
* Bad, because modern individuals are ~30% heterozygous at autosomal SNPs —
  discarding 30% of sites severely reduces power and introduces ascertainment bias
  (remaining sites skew toward low-frequency variants where the individual is hom)
