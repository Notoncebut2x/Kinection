# ADR-0013: Admixture decomposition via constrained NNLS

* Status: Accepted
* Date: 2026-05-17
* Supersedes: —
* Phase: 1 (Analysis engine — Step 1.5)

## Context

Step 1.5 of the pipeline produces the headline "you are X% Steppe, Y% Anatolian Farmer, Z% Levantine, …" result. This is what most ancient-DNA papers report as their primary finding for modern individuals.

The gold-standard method is **qpAdm** (Reich Lab), which uses *outgroup f₄-statistics* to estimate admixture proportions while controlling for shared drift between sources via a panel of "right populations." qpAdm is robust to confounding and produces well-calibrated confidence intervals — but it requires:

- A carefully chosen set of left/right populations
- Computation of f-statistics, ideally with the `admixtools2` R package or `qpAdm` C binary from EIGENSOFT
- Decisions about whether a given source model is *rejected* (P > 0.05) before reading off proportions

For a single-individual personal-genomics workflow, qpAdm is heavy machinery and adds a non-Python dependency.

## Decision

Step 1.5 uses **constrained non-negative least squares** (NNLS) over per-source allele frequencies as a lightweight alternative.

For the target individual's allele frequency vector **b** ∈ [0, 1]^M (where M is the number of overlapping autosomal SNPs) and a matrix **A** ∈ ℝ^(M × K) of source-population allele frequencies, we solve:

```
minimize  ‖A·α − b‖²
subject to   αᵢ ≥ 0,   Σ αᵢ = 1
```

implemented via `scipy.optimize.minimize` with the SLSQP solver, equality constraint on the sum, and box bounds on each component.

95 % confidence intervals come from **block-bootstrap by chromosome** (200 iterations by default, configurable via `ADMIX_BOOTSTRAP=N`). Bootstrapping at the chromosome level rather than the SNP level approximately respects linkage disequilibrium — adjacent SNPs are not independent observations, so SNP-level resampling underestimates variance.

## Source populations

The default 6-source manifest (in `scripts/step1_5_admixture.py`):

| Source | AADR groups (substring patterns) | n typical |
|---|---|---|
| WHG | Germany/France/Spain/Sweden/Belgium Mesolithic | ~90 |
| EHG | Russia_YuzhniyOleniyOstrov, Russia_Minino, Karelia_HG | ~44 |
| EEF | Turkey_Marmara_Barcin_N, Çatalhöyük, Menteşe | ~57 |
| Steppe | Yamnaya (various) + Afanasievo | ~59 |
| Levant_N | Jordan_PPNB + Israel_Natufian + Cyprus_PPNB | ~25 |
| Iran_N | Iran_GanjDareh_N + Georgia_Kotias/Satsurblia | ~14 |

These choices follow the canonical Western Eurasian model from Lazaridis et al. 2014, 2016, 2022 and Allentoft et al. 2024.

## Considered alternatives

**qpAdm / admixtools2.** Rejected for v1 because it adds an R or compiled-C dependency. The Phase 1 goal is a Python-only pipeline. May be added behind an optional flag in a later phase.

**ADMIXTURE / fastNGSadmix.** Rejected because these are *unsupervised* or *semi-supervised* — they discover their own ancestry components from a reference panel rather than fitting against named ancient sources. The user's question is "how much of my genome is Yamnaya?", not "what clusters does my genome contain?"

**Pure NNLS without sum-to-1.** Rejected because the proportions then don't add up to 100 %, complicating interpretation and report wording.

## Consequences

**Positive:**
- Pure Python; no R, no compiled binaries.
- Runs in ~30 seconds for ~380,000 SNPs × 6 sources × 200 bootstrap iterations.
- Headline results are directly comparable in magnitude to published qpAdm decompositions for similar individuals.
- Easy to extend: edit the `SOURCES` dict in the script.

**Negative:**
- NNLS treats source populations as homogeneous and ignores shared drift between sources. Highly correlated sources (e.g. EEF and Levant_N) can lead to weight swaps the model can't distinguish — both are Neolithic farmers and have similar allele frequencies.
- No statistical test of model fit. qpAdm reports a P-value for whether the chosen source set adequately explains the target; NNLS just returns the best linear combination available.
- Confidence intervals reflect statistical (sampling) noise only — they do not capture model misspecification.

**Operational lessons:**
- EIGENSTRAT PACKGENO encodes the *count of the reference allele*. Step 1 of this pipeline computes the modern individual's `dosage` as the *count of the alt allele*. The source allele frequencies must be inverted (`1 − mean(geno)/2`) so target and sources count the same allele. Getting this wrong produces an inverted decomposition that gives Levant_N ~64 % even for non-Levantine targets.
- Caveats around interpretation (especially the Levant_N component as deep Levantine ancestry vs Ashkenazi-specific) are surfaced in the report itself.

## References

- Lazaridis, Patterson, et al. 2014. "Ancient human genomes suggest three ancestral populations for present-day Europeans." *Nature* 513.
- Lazaridis et al. 2016. "Genomic insights into the origin of farming in the ancient Near East." *Nature* 536.
- Lazaridis et al. 2022. "The genetic history of the Southern Arc." *Science* 377.
- Allentoft et al. 2024. "Population genomics of post-glacial western Eurasia." *Nature* 625.
- Haak et al. 2015. Reach Lab. Original NNLS-with-sum-to-1 approach for Yamnaya-related ancestry in modern Europeans.
