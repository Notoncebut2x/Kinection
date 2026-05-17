# Kinection — The Science, in Plain English

*Last updated: 2026-05-17*

This document explains what each step of the analysis pipeline actually does, why it works, and what the results mean — without assuming a background in population genetics.

For the system architecture, see [`ARCHITECTURE.md`](ARCHITECTURE.md). For the technical details of each step, see the code in `scripts/`.

---

## What we're trying to do

You have a DNA test from a consumer service (AncestryDNA). It tells you about *modern* populations — "23% Irish, 11% Scandinavian," and so on.

What it doesn't tell you is which *ancient* peoples your DNA most resembles. Were your ancestors among the first farmers who spread from Anatolia 8,000 years ago? The Yamnaya horse-riders who swept across Europe from the steppes 5,000 years ago? The Iron Age Celts? The Vikings?

Researchers have spent the last 15 years sequencing DNA from skeletons across the world — about 19,000 ancient genomes are now public, dated and labelled by archaeologists. This pipeline takes your modern DNA and compares it against all of them.

The output answers three questions:

1. **Which paternal and maternal lineage do you belong to?** (haplogroups)
2. **Which specific ancient individuals look most similar to you across your whole genome?**
3. **Which broader ancient *populations* do you resemble most?**

---

## A 60-second crash course

Your DNA is a 3-billion-letter string. Out of those 3 billion, only about a million places vary between people — these are the **SNPs** (single nucleotide polymorphisms) that DNA tests look at.

At each SNP location, you have two letters — one from each parent. Most are nearly identical between humans, but a small number are informative about geography and history because they shifted in frequency as populations moved, mixed, and split apart.

Two people with **similar SNPs** at many positions are likely descended from a similar ancestral population. That's all this is doing — measuring similarity, very carefully, between you and ~19,000 ancient skeletons.

Three things make it more interesting than that:

- The Y chromosome and mitochondrial DNA pass down *un-recombined* through fathers and mothers respectively, so they trace specific lineages with great precision.
- Ancient DNA is degraded — you usually only recover one of the two alleles at each site. That has to be accounted for or it biases everything.
- Around 10% of SNPs are technically "ambiguous" between DNA strands and have to be thrown out, or they introduce nonsense matches.

The pipeline handles all three.

---

## The Reference Dataset — AADR

**What it is:** The *Allen Ancient DNA Resource* (AADR), curated by David Reich's lab at Harvard. It is a merged, quality-controlled collection of ~19,000 ancient human genomes — every individual reported in scientific publications since ~2010, dated and labelled by archaeologists.

**The format:** A "genotype matrix." Rows are SNP positions (about 1.24 million of them). Columns are individuals. Each cell is one of: 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate), or "missing." For ancient samples, most cells are missing or "pseudo-haploid" — only one allele was recoverable.

**Why we use it:** It's the only comprehensive, professionally-curated resource of ancient genomes. Crucially, the labels (date, location, archaeological culture) are reviewed and standardized. You can't reproduce this work yourself by trawling Genbank.

**Versioning:** New releases come out a few times a year as more skeletons are sequenced. Current: v66 (released April 2026). `scripts/update_aadr.py` checks Harvard Dataverse for new releases and updates the cloud copy.

See [ADR-0001](decisions/0001-reference-dataset-aadr-v62.md) for the formal decision record.

---

## Step 1 — Parse and Harmonise

**Script:** `scripts/step1_parse_harmonise.py`

### What it does

Reads your AncestryDNA file (about 660,000 SNPs) and the AADR reference (about 1.24 million SNPs), finds the SNPs they have in common, and produces an "encoded" version of your genome that the rest of the pipeline can work with.

### Why this is harder than it sounds

Three problems have to be solved before any comparison is possible:

**Problem 1: Different SNP sets.** AncestryDNA and AADR don't measure the same SNPs. Only the overlapping ~395,000 SNPs are usable.

**Problem 2: Strand ambiguity.** DNA has two strands (think of the double helix). The same SNP can be reported as "A vs G" by one lab and "T vs C" by another — both are correct because A↔T and G↔C are complementary. Most of the time the asymmetric letters (A vs G) make this obvious, but for **palindromic SNPs** (A vs T, or C vs G) there's no way to know which strand was reported. These ~10% of SNPs have to be **excluded** or they will produce false matches.

**Problem 3: Allele encoding.** AncestryDNA reports your two alleles as letters ("AG" at a position). AADR reports them as 0/1/2 (the count of "alternate" alleles). They have to be put into the same units.

### What you get

- `snp_overlap.tsv` — one row per common SNP, with its position, reference/alternate alleles, and your genotype encoded as 0/1/2/missing
- `step1_summary.json` — quick stats: how many SNPs overlapped, how many were palindromic and excluded, how many you had genotypes for
- `modern_indv_encoded.npy` — your genotypes as a compact NumPy array (input to steps 2 and 3)

For Rory's actual data: 394,671 SNPs survived the overlap and the palindromic filter. That's a good number — population studies routinely work with fewer.

See [ADR-0005](decisions/0005-palindromic-snp-exclusion.md) for the palindromic exclusion decision.

---

## Step 2 — Haplogroup Assignment

**Script:** `scripts/step2_haplogroup.py`

### What it does

Determines your **paternal lineage** (Y-DNA haplogroup, men only — women inherit no Y chromosome) and **maternal lineage** (mitochondrial DNA haplogroup, everyone). Then finds ancient individuals who share these lineages.

### The science

The Y chromosome and mitochondrial DNA (mtDNA) are special: they pass down without recombining. Your Y chromosome is essentially identical to your father's, his father's, his father's father's, and so on, except for occasional mutations. Same with mtDNA through the maternal line.

Researchers have built **phylogenetic trees** of every Y-DNA and mtDNA lineage on Earth. Specific mutations ("markers") define branching points. For example:

- The mutation **M9** defines haplogroup K (~45,000 years old, originated in Asia)
- **M207** defines R (a descendant of K, ~30,000 years old)
- **M269** defines R1b1a1b (~6,500 years old, the dominant male lineage in Western Europe today, associated with the spread of Bronze Age steppe populations)

By checking which of these marker mutations you have, the script walks down the tree from the root to the deepest branch where your DNA still matches.

### The catch

AncestryDNA tests don't cover the Y chromosome densely — only ~34 ISOGG markers are usually present. So you'll typically get a confident assignment down to a major branch (e.g., R1b1a1b), but not all the way to the leaves (the full ISOGG tree has ~25,000 named sub-haplogroups). For finer resolution you'd need Big-Y or YFull sequencing.

mtDNA is similar: the AncestryDNA array covers ~190 mtDNA positions, enough for the major haplogroup (R, H, U, J, etc.) but not the full HVR1/HVR2 resolution that dedicated mtDNA sequencing gives.

### Matching to ancients

Once your haplogroup is known, the script searches the AADR annotation file for ancient individuals with the same lineage, ranked by how *deep* the shared branching point is (a match at "R1b1a1b" is closer than a match at "R"). The best matches — individuals who share both Y and mtDNA with you — are flagged separately.

### What you get

- `ydna_haplogroup.json` — your Y haplogroup with confidence and supporting marker list
- `mtdna_haplogroup.json` — your mtDNA haplogroup with confidence
- `ancient_haplogroup_matches.tsv` — top ancient individuals sharing your lineages, with archaeological context
- `haplogroup_report.md` — a human-readable narrative of the above

For Rory: Y-DNA **R1b1a1b**, mtDNA **R**. The Y assignment places his paternal line firmly in the post-Bronze Age Western European population that derives from Yamnaya / Steppe ancestry. The mtDNA at the macro-R level is consistent with a West Eurasian maternal lineage but not more specific without dedicated mtDNA sequencing.

See [ADR-0006](decisions/0006-haplogroup-reference-databases.md) for the marker database choices (ISOGG for Y, PhyloTree B17 for mtDNA).

---

## Step 3 — Genome-wide Similarity and PCA

**Script:** `scripts/step3_similarity_pca.py`

### What it does

This is the big one. It computes a single number — the **Allele-Sharing Distance** (ASD) — between you and *every* ancient individual in the AADR. Then it ranks them, groups them by archaeological population, and produces a PCA plot showing where you fall in the space of ancient human genetic variation.

### The science: Allele-Sharing Distance

For each SNP that both you and an ancient individual have genotyped:

- If both have 0 alternate alleles → identical → distance contribution: 0
- If both have 2 alternate alleles → identical → distance contribution: 0
- If you have 0 and they have 2 (or vice versa) → totally different → distance contribution: 1
- If one of you is heterozygous (1) → distance contribution: 0.5

Average that across all shared SNPs. The result is the **mean allele-sharing distance**: 0.0 = identical twins, 0.5 = unrelated humans, 0.7+ = non-human.

### The complication: pseudo-haploidy

Ancient DNA is degraded. From a 5,000-year-old skeleton you usually only recover *one* of the two alleles at any given site. So ancient genotypes in AADR are encoded as pseudo-haploid: they look like 0 or 2 (one allele drawn at random), never 1.

You, on the other hand, are sequenced from saliva. Your genotypes are honestly diploid (0, 1, or 2).

If we compare diploid-vs-haploid directly, an ancient "0" against your "1" looks like distance 0.5 — but it might really have been a 0 against your 0 (distance 0) or a 0 against your 2 (distance 1). The mid-value of 0.5 introduces a systematic bias.

The fix: **pseudo-haploidise yourself**. At each heterozygous site, randomly draw one of your two alleles. Do this 10 times with different random seeds, average the distances. The result is a fair comparison between two pseudo-haploid representations.

See [ADR-0003](decisions/0003-pseudo-haploidisation.md) for the formal rationale.

### Why per-individual then per-population

The per-individual ranking is interesting — you can identify which specific ancient skeleton you most resemble. But individual ancient DNA quality varies enormously: some samples have only 10,000 covered SNPs, others have 1.2 million. A single noisy sample can rank high or low for the wrong reasons.

So the script also averages distances by **population** (e.g., "Britain_IronAge_Roman.SG" or "Spain_Bronze_Age.AG"). Per-population averages are far more stable because they aggregate over multiple individuals.

### PCA — putting it visually

Principal Component Analysis is a way of distilling many-dimensional data down to two (or ten) dimensions that capture the most variation. Apply PCA to the matrix of ancient individuals × SNPs and the first two components typically separate Europe vs. Asia vs. Africa, etc. — the classic "PCA of human variation" plots you've seen in news articles.

The script projects *you* into that ancient-only space (your genome is held out from the computation, then projected onto the axes the ancients defined). Your position on PC1 vs PC2 tells you which broad ancestral cluster you fall into, and which ancient populations you sit closest to in this distilled space.

### Why this takes a while

The AADR `.geno` file is ~7 GB. We read it via HTTP byte-range requests directly from R2 — fetching only the SNP rows we need, in 5,000-row chunks. Total runtime is around 10–20 minutes on a normal connection. See [ADR-0010](decisions/0010-chunked-geno-processing.md) and [ADR-0011](decisions/0011-cloudflare-r2-geno-storage.md).

### What you get

- `pairwise_distances.tsv` — your ASD vs all ~19,000 ancient individuals, ranked best-first
- `population_distances.tsv` — same data averaged by population (~3,400 populations)
- `pca_coordinates.tsv` — your PC1–PC10 scores alongside all ancient individuals
- `pca_variance_explained.json` — how much of human variation each PC captures
- `top_matches_report.md` — narrative of your closest populations

---

## What the results actually *mean*

A few honest caveats:

**"Closest ancient population" ≠ "your ancestors"**

Genetic similarity reflects shared ancestry, but human populations have been mixing continuously for thousands of years. If you match Iron Age Britons closely, it's because the gene pool that produced them and the gene pool you came from share a lot of ancestral components — not because those specific Britons are your direct ancestors.

**Haplogroup ≠ ethnicity**

Y-DNA R1b1a1b is found at high frequency in Western Europe, but also among populations as far apart as Bashkirs (Russia) and parts of Cameroon. A haplogroup is one specific line of descent through your father's-father's-father's line; it represents about 1/2^N of your ancestry going back N generations. Your *autosomal* DNA (everything else) is far more representative of "who you are descended from."

**Confidence depends on coverage**

If your matches are based on 400,000 SNPs and the ancient individual has 50,000 SNPs sequenced, only ~50,000 of those overlap and contribute to the comparison. Low-coverage ancients in the top matches list should be taken with a grain of salt. The population-level results aggregate over many individuals and are more reliable.

**The reference is biased**

AADR contains far more European and Near Eastern samples than African or Pacific samples, simply because more research has been done in those regions. If your ancestry is from an under-represented region, the "closest population" results will be less informative.

---

## What's coming

Currently steps 1.1, 1.2, and 1.3 (parse-harmonise, haplogroup, similarity+PCA) are implemented. Planned next:

- **Step 1.4 — TMRCA estimation** — How long ago did you and a given ancient individual last share a common ancestor? This adds a *time* dimension to the matches.
- **Step 1.5 — AMOVA / admixture decomposition** — Decompose your ancestry into ancient-population components: e.g., "55% Steppe Bronze Age, 30% Anatolian Farmer, 15% Western Hunter-Gatherer." This is the headline number for most ancient-DNA analyses of modern Europeans.
- **Step 1.6 — Synthesis report** — A single human-readable PDF combining steps 1.1–1.5 with maps and a timeline.

Then Phase 2 (web frontend) — a browser-based way for users to upload, kick off, watch, and read their report.

---

## Further reading

- David Reich, *Who We Are and How We Got Here* (2018) — the popular-science book by the lab that produces AADR. Best one-volume introduction to ancient DNA.
- The AADR documentation: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
- ISOGG Y-tree (the canonical Y-DNA reference): https://isogg.org/tree/
- PhyloTree (the canonical mtDNA reference): https://www.phylotree.org/
