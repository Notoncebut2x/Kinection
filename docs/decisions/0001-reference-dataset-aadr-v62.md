# ADR-0001: Use Allen Ancient DNA Resource v62.0 as Reference Dataset

* Status: Accepted
* Date: 2026-04-21

## Context and Problem Statement

The project requires a large, curated reference dataset of ancient human genotypes to compare against modern individual DNA. Several public ancient DNA compendiums exist. Which dataset provides the best combination of scale, data quality, and methodological credibility for a consumer-facing ancestry product?

## Decision Drivers

* Maximum number of ancient individuals for comparison coverage
* High SNP density at consistent genomic positions (required for ASD and PCA)
* Published quality assessments per individual (to filter low-quality samples)
* Active maintenance and versioning (dataset will need to be updated over time)
* Used as the reference in published peer-reviewed population genetics studies
* Free, publicly available without licensing restrictions

## Considered Options

* AADR v62.0 (Allen Ancient DNA Resource, Harvard Dataverse)
* Reich Lab 1240k dataset (subset of AADR, earlier vintage)
* SGDP (Simons Genome Diversity Project) — modern populations only
* Published per-study VCFs (e.g., individual papers from Science/Nature)

## Decision Outcome

Chosen option: **AADR v62.0**, because it is the largest curated compendium of ancient human genotypes available, is actively maintained with versioned releases, includes per-individual quality assessments (PASS/QUESTIONABLE/IGNORE/FAIL), and is the direct reference used in the methodology paper this project mirrors (Deep Maniot Greeks, Communications Biology 2026).

### Positive Consequences

* 17,629 ancient individuals — maximum global coverage by a wide margin
* 406,570 SNPs on the 1240k capture panel — consistent positions across all samples
* Per-individual `.anno` file includes haplogroups, dates, geographic coordinates,
  culture labels, and quality assessments — eliminates need for external annotation
* Version number (v62.0) provides a stable reference for result reproducibility
* Directly cited in the methodology paper this project replicates

### Negative Consequences

* Dataset is large (~46 GB); hosting a full copy for production requires significant storage
* Some ancient individuals have very low SNP coverage, requiring overlap thresholds to filter
* Ancient samples are pseudo-haploid (one allele drawn randomly per position), which
  complicates comparison with diploid modern individuals and requires methodological adjustments
* Dataset versioning means results will shift when AADR releases a new version

## Pros and Cons of the Options

### AADR v62.0

* Good, because largest available dataset (17,629 individuals)
* Good, because includes comprehensive metadata (.anno) used by academic studies
* Good, because versioned and actively updated by Harvard lab
* Good, because directly used in the Deep Maniot reference methodology
* Bad, because 46 GB storage requirement for production
* Bad, because pseudo-haploid encoding requires extra methodological handling

### Reich Lab 1240k (earlier vintage)

* Good, because well-established, heavily cited
* Bad, because superseded by AADR; fewer individuals and less metadata
* Bad, because no active maintenance path

### Per-study VCFs

* Good, because some studies include full genome data (not just 1240k capture)
* Bad, because inconsistent SNP panels across studies — cannot be combined without
  significant harmonisation work
* Bad, because no single unified quality assessment framework
* Bad, because geographically and temporally uneven coverage
