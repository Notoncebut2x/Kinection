# ADR-0005: Exclude Palindromic SNPs During Strand Alignment

* Status: Accepted
* Date: 2026-04-22

## Context and Problem Statement

Consumer DNA arrays (AncestryDNA, 23andMe) and the AADR ancient dataset may report alleles relative to different strands of the DNA double helix. When aligning the two datasets, most SNPs can be strand-corrected by taking the complement (A↔T, C↔G). However, some SNPs are palindromic — the two possible alleles are complements of each other (A/T or C/G SNPs). For these SNPs, flipping the strand looks the same as not flipping it, making it impossible to determine the correct orientation without additional frequency information.

Should palindromic SNPs be excluded, or should they be retained using allele frequency comparison to determine strand orientation?

## Decision Drivers

* Accuracy — incorrect strand orientation at a SNP produces a dosage error (0↔2) that
  will inflate distance estimates for that individual
* Scale — there are tens of thousands of palindromic SNPs in a typical 700k array overlap
* Complexity — frequency-based strand resolution adds implementation complexity and
  requires a population frequency reference
* Alignment with published methodology — the Deep Maniot study excludes palindromic SNPs

## Considered Options

* Exclude all palindromic (A/T and C/G) SNPs unconditionally
* Retain palindromic SNPs using allele frequency comparison to determine strand orientation
  (if both datasets agree that the minor allele frequency is similar, infer orientation)
* Retain palindromic SNPs only where one allele is very rare (MAF < 0.05 rule) —
  low-frequency alleles are less likely to create ambiguity

## Decision Outcome

Chosen option: **Exclude all palindromic SNPs unconditionally**, because incorrect strand orientation at a palindromic site produces a hard error (dosage 0 miscalled as 2 or vice versa) that silently inflates distances, because the frequency-based approach adds significant complexity and requires a population reference we do not have at parse time, and because the Deep Maniot methodology uses the same exclusion. The loss in SNP count (~10–15%) is acceptable given the 394k overlap still provides excellent statistical power.

### Positive Consequences

* Zero risk of systematic strand-flip errors corrupting distance estimates
* Simple, deterministic logic — no external frequency reference required
* Consistent with the published methodology this project mirrors
* Remaining 394k SNPs provide more than sufficient power for ASD and PCA

### Negative Consequences

* ~10–15% of overlapping SNPs are discarded — small reduction in statistical power
* A/T and C/G SNPs are not evenly distributed across the genome; exclusion may
  slightly skew representation at some loci
* Future users with very low overlap counts (< 10k SNPs) may be more affected
  by this exclusion — worth monitoring in production QC

## Pros and Cons of the Options

### Exclude palindromic SNPs (chosen)

* Good, because eliminates all risk of silent strand-flip errors
* Good, because simple to implement — single frozenset check
* Good, because matches published ancient DNA pipeline standards
* Bad, because discards ~10–15% of potential overlap SNPs

### Frequency-based strand resolution

* Good, because retains more SNPs, increasing statistical power
* Good, because widely used in GWAS harmonisation pipelines (e.g., LDSC)
* Bad, because requires a population allele frequency reference that is not
  available at parse time without a significant additional data dependency
* Bad, because incorrect resolution (when frequencies are very similar) produces
  silent systematic errors — worse than exclusion
* Bad, because adds substantial implementation complexity

### MAF < 0.05 rule

* Good, because a subset of palindromic SNPs (those with very skewed frequencies)
  can be reliably resolved
* Bad, because still requires a frequency reference
* Bad, because the safe-to-retain fraction is small — most palindromic SNPs have
  intermediate frequencies and remain ambiguous
* Bad, because partial retention adds complexity for marginal SNP gain
