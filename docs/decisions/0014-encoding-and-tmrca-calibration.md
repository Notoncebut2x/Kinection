# ADR-0014: AADR PACKGENO encoding convention and Y-DNA TMRCA calibration

* Status: Accepted
* Date: 2026-05-19
* Supersedes: —
* Phase: 1 (Analysis engine — affects Steps 1.1, 1.3, 1.4, 1.5)

## Context

Two related decisions, both forced by surprises discovered while building Step 1.4 (Y-DNA TMRCA).

### 1. The AADR PACKGENO encoding inversion

The EIGENSOFT/EIGENSTRAT documentation describes PACKGENO genotype values as "count of the reference allele" (0 = 0 ref, 2 = 2 ref). The AADR-specific documentation, however, says:

> The .geno file values: 0 = homozygous for allele2 (column 6 of .snp); 2 = homozygous for allele1 (column 5); 1 = heterozygous; 3 = missing.

Step 1.1 (`step1_parse_harmonise.py`) computes `dosage = count of allele2 in the modern individual`. So the encoding alignment between modern and ancient is:

| Modern | Meaning | Ancient (AADR) | Meaning |
|---|---|---|---|
| `dosage = 0` | hom allele1 | `geno = 2` | hom allele1 |
| `dosage = 2` | hom allele2 | `geno = 0` | hom allele2 |

The two encodings are **inverted relative to each other**. Naively comparing `dosage` to `geno` (e.g. `|dosage/2 − geno/2|`) gives maximum distance when the samples are identical and minimum distance when they are opposite — the rankings are fully wrong.

This was discovered when Step 1.4 produced TMRCAs of ~640 million years between an R1b modern and R1b ancients (k/L ≈ 95% Y-SNP differences, which is biologically impossible). Investigation traced the same bug to **Step 1.3** (`step3_similarity_pca.py`): its top "closest" matches for a Northern European individual were Mayan, Bantu, and Brahui samples (the actual most-distant populations), and its mean ASD was 0.57 instead of the expected ~0.25–0.30 for European-vs-European pseudo-haploid comparison.

**Step 1.5** (`step1_5_admixture.py`) had already discovered this independently and applied the inversion locally (line ~233: `source_alt_freq = 1 - mean(geno)/2`), but the discovery was not documented and the fix was not propagated to the rest of the pipeline.

### 2. Y-DNA TMRCA calibration: why per-bp rates don't work for ascertained panels

The workplan (and the methodology paper) specifies a Y mutation rate of 0.74 × 10⁻⁹ substitutions/bp/year (Karmin et al. 2015). This rate is calibrated for *random* Y chromosome bp and is the right number for whole-Y-sequencing studies (Big-Y, YFull, full-Y BAMs).

It is the **wrong** rate for the 1240k Y panel. The panel contains 32,670 Y SNPs selected from the ~10 Mb callable Y because they were observed to be polymorphic. Using the per-bp rate with `L = number of called panel SNPs` gives TMRCAs in the tens of millions of years — three to four orders of magnitude too high — because the panel SNPs each represent variation accumulated over hundreds of bp of evolutionary history.

The cleanest alternative would be to use `L = effective callable bp` (~10 Mb scaled by sample coverage), but this requires assumptions about which non-panel positions would have agreed, which we cannot verify from array data alone.

## Decisions

### Decision 1: Normalise to a single "allele2 frequency" space across the pipeline.

All cross-modern/ancient comparisons normalise both sides to the same convention: the value represents allele2 frequency (where allele2 is column 6 of the AADR `.snp` file). Concretely:

- Modern: `modern_allele2_freq = dosage / 2` (already correct under step 1's encoding).
- Ancient: invert AADR PACKGENO via `ancient_allele2_freq = 1 − geno/2` for geno ∈ {0, 2}, with geno=1 → 0.5 and geno=3 → missing.

Equivalently, the lookup table used in Step 1.3's ASD and PCA is now `geno_to_freq = [1.0, 0.5, 0.0, nan]` (the inverse of the previous `[0.0, 0.5, 1.0, nan]`), and the PCA dosage lookup is `geno_to_dosage = [2, 1, 0, -1]` (inverse of `[0, 1, 2, -1]`).

Step 1.4 uses a normalised `{0, 1}` haploid space derived the same way for Y SNPs.

Step 1.5 already does the equivalent thing locally and is left as-is, with a cross-reference comment to this ADR.

### Decision 2: Y-DNA TMRCA uses a calibrated per-panel-SNP rate, not the per-bp rate.

Step 1.4 uses `μ_panel = 7 × 10⁻⁶ per panel-SNP per year`, anchored to the known R1b clade coalescence age of ~20 ky. The reported TMRCAs are explicitly described as **order-of-magnitude estimates**, accurate to within a factor of ~2; for tighter ages, users need whole-Y sequencing (Big-Y, YFull).

Step 1.4's primary reported metric is the raw `k/L` Y-SNP difference rate (rate-calibration-free, directly comparable across matches). TMRCA in years is reported as a secondary metric with prominent caveats in the markdown report.

## Consequences

**Positive:**
- Step 1.3 ASD rankings now reflect actual genetic similarity. For a Northern European modern individual the closest matches become European/Mediterranean ancients, as expected.
- Step 1.4 produces plausible Y-DNA TMRCAs (3,500–4,000 y for R1b modern vs R1b Bronze Age/Medieval ancients).
- The encoding convention is documented in one place; future steps can rely on it.

**Negative:**
- Any previously-saved Step 1.3 outputs (`pairwise_distances.tsv`, `population_distances.tsv`, `top_matches_report.md`, PCA coordinates) are invalid and must be regenerated. The fix invalidates the existing rn and jn step 3 outputs.
- The reported TMRCA confidence intervals reflect Poisson sampling noise only, not the ±factor-of-2 rate uncertainty. Quoted absolute ages should be widened accordingly in user-facing reports.
- The calibrated rate (7 × 10⁻⁶ /panel-SNP/year) is empirical and ties results to a single anchor (R1b ~20 ky). A future iteration could recalibrate against multiple known coalescence dates (R1a, I1, J1 …) to reduce single-anchor bias.

## Alternatives considered

**A. Rename step 1's `dosage` column to `allele2_dosage` and leave step 3's `geno_to_freq` as `[0.0, 0.5, 1.0, nan]` interpreted as "allele1 frequency".**

Equivalent math, but requires renaming a column written by an upstream script and re-reading all downstream consumers. The chosen fix is minimal and local.

**B. Use whole-Y callable bp as TMRCA denominator with the per-bp rate.**

Requires knowing the effective callable region of the 1240k Y panel and assuming all non-panel positions in that region match between the two samples. The first is poorly documented; the second is empirically unverifiable from array data. The calibrated per-panel-SNP rate sidesteps both problems but ties the estimate to one anchor.

**C. Report `k/L` only and skip TMRCA-in-years entirely.**

Cleanest, but loses the headline narrative payoff ("your common ancestor with this person lived ~3,500 years ago"). Chosen compromise: report both, flag TMRCA as approximate.

## References

- Karmin et al. 2015 — Y mutation rate calibration: https://doi.org/10.1101/gr.186684.114
- AADR documentation (Harvard Dataverse): https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
- Bug discovery: see commit and step1_4 development session, 2026-05-19.
