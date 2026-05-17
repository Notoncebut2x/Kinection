# ADR-0006: Use ISOGG Y-DNA Tree and PhyloTree B17 for Haplogroup Assignment

* Status: Accepted
* Date: 2026-04-22

## Context and Problem Statement

Haplogroup assignment requires a reference database of branch-defining markers and their expected allele states. For Y-DNA and mtDNA, several community-maintained and commercial databases exist. Which reference databases should be used for this project, and how should they be embedded?

## Decision Drivers

* Free and openly licensed — no commercial dependency for a product that may charge users
* Stable enough to be embedded as static JSON files in the codebase
* Sufficient marker coverage for reliable haplogroup assignment from AncestryDNA array data
  (which covers a limited subset of Y and mtDNA markers compared to dedicated sequencing)
* Compatible with the haplogroup naming scheme used in the AADR .anno file
  (so that modern assignment can be matched to ancient haplogroup strings)
* Widely cited and accepted in the population genetics literature

## Considered Options

* ISOGG Y-DNA haplogroup tree (isogg.org) + PhyloTree B17 (mtDNA) — embedded as JSON
* Yfull haplogroup tree — more detailed, mutation-rate calibrated
* FamilyTreeDNA haplotree — commercial, used in the Deep Maniot study for TMRCA
* Automated tools: Haplogrep 3 (mtDNA), YHRD tools (Y-DNA)

## Decision Outcome

Chosen option: **ISOGG Y-DNA + PhyloTree B17**, both embedded as static JSON databases in `scripts/data/`. ISOGG provides the standard haplogroup naming scheme that matches the AADR .anno file (e.g., R-M269, I-M423). PhyloTree B17 is the last full revision of the consensus mtDNA phylogenetic tree and covers all major haplogroups detectable on the AncestryDNA array. Both are freely available and do not require API calls or external dependencies at runtime.

### Positive Consequences

* Zero external runtime dependencies — haplogroup assignment runs offline
* ISOGG haplogroup names align directly with AADR .anno strings, enabling prefix-match
  comparison between modern assignment and ancient haplogroup fields
* Both databases are stable (PhyloTree has not been revised since 2016; ISOGG tree
  is versioned) — embedded JSON will not change unexpectedly
* Freely licensed — no legal barrier to commercial use

### Negative Consequences

* AncestryDNA array covers only ~15–34 Y-DNA markers (of ~40,000+ known SNPs on
  the Y chromosome) — assignment will typically resolve to broad haplogroups
  (e.g., R1b1a1b, not R-L21 or R-DF27 sub-branches)
* mtDNA coverage on the AncestryDNA array (~190 positions) resolves major haplogroups
  (H, U, J, T, K) but not sub-haplogroups (e.g., H1, H3, H5) — PhyloTree sub-branch
  markers are present in the database but often not genotyped on the array
* PhyloTree B17 (2016) does not reflect mtDNA phylogenetic revisions from 2017 onward —
  some reclassifications in H and U subclades may differ
* ISOGG tree does not include mutation rate estimates — TMRCA (Step 1.4) will require
  a separate mutation rate calibration (Yfull or published literature)

## Pros and Cons of the Options

### ISOGG + PhyloTree B17 embedded JSON (chosen)

* Good, because free and openly licensed
* Good, because ISOGG naming matches AADR .anno haplogroup strings directly
* Good, because offline — no API call required at analysis time
* Good, because stable — no surprise updates
* Bad, because Y-DNA resolution limited by AncestryDNA array coverage
* Bad, because PhyloTree B17 is slightly dated (2016)

### Yfull haplogroup tree

* Good, because more granular — includes thousands of sub-branches with age estimates
* Good, because better TMRCA calibration (mutation rates per branch)
* Bad, because licensing is unclear for commercial use
* Bad, because haplogroup naming scheme differs from ISOGG (used in AADR) —
  requires mapping layer
* Bad, because requires a more complex tree traversal implementation

### FamilyTreeDNA haplotree

* Good, because used in the Deep Maniot study for TMRCA estimation
* Bad, because proprietary — embedding in a product requires a commercial agreement
* Bad, because API-access only; not available as a static download

### Haplogrep 3 / automated tools

* Good, because well-validated, handles VCF/FASTA input natively
* Good, because maintained by the academic community
* Bad, because requires a subprocess call — adds external dependency
* Bad, because designed for full mtDNA sequences, not microarray data
* Bad, because Y-DNA tools (YHRD) are oriented toward forensic STR data, not SNP arrays
