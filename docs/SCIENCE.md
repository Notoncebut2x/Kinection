# Kinection — The Science, in Plain English

*Last updated: 2026-05-18*

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

**Versioning:** New releases come out a few times a year as more skeletons are sequenced. Current cloud version: **v66** (~7.2 GB .geno, in Cloudflare R2 at `dataset/v66/`). A v62 copy is also kept on local disk for fully-offline runs. `scripts/update_aadr.py` checks Harvard Dataverse for new releases and updates the cloud copy; the runtime resolves which version to use from the manifest at `dataset/current_version.json`.

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

Typical surviving overlap after the palindromic filter is in the range of 400,000–410,000 SNPs for an AncestryDNA V2 array against AADR v62 — well above what most population genetics studies work with.

See [ADR-0005](decisions/0005-palindromic-snp-exclusion.md) for the palindromic exclusion decision.

---

## Step 2 — Haplogroup Assignment

**Script:** `scripts/step2_haplogroup.py`

### What it does

Determines your **paternal lineage** (Y-DNA haplogroup, men only — women inherit no Y chromosome) and **maternal lineage** (mitochondrial DNA haplogroup, everyone). Then finds ancient individuals who share these lineages.

### The science

The Y chromosome and mitochondrial DNA (mtDNA) are special: they pass down without recombining. Your Y chromosome is essentially identical to your father's, his father's, his father's father's, and so on, except for occasional mutations. Same with mtDNA through the maternal line.

Researchers have built **phylogenetic trees** of every Y-DNA and mtDNA lineage on Earth. Specific mutations ("markers") define branching points. For example, on the Y-DNA tree:

- The mutation **M9** defines haplogroup K (~45,000 years old, originated in Asia)
- **M207** defines R (a descendant of K, ~30,000 years old)
- And so on through dozens more markers down to terminal sub-clades, each tagged with a name like "R-Z2103" or "I-M253" that encodes its position in the tree

By checking which of these marker mutations you have, the script walks down the tree from the root to the deepest branch where your DNA still matches. The current ISOGG tree has roughly 25,000 named branches.

### The catch

AncestryDNA tests don't cover the Y chromosome densely — only ~34 ISOGG markers are usually present. So you'll typically get a confident assignment down to a major branch, but not all the way to the leaves of the tree. For finer resolution you'd need Big-Y or YFull sequencing.

mtDNA is similar: the AncestryDNA array covers ~190 mtDNA positions, enough for the major haplogroup (R, H, U, J, etc.) but not the full HVR1/HVR2 resolution that dedicated mtDNA sequencing gives.

### Matching to ancients

Once your haplogroup is known, the script searches the AADR annotation file for ancient individuals with the same lineage, ranked by how *deep* the shared branching point is (a match at a fine-grained sub-haplogroup is closer than a match only at the macro level). The best matches — individuals who share both Y and mtDNA with you — are flagged separately.

### What you get

- `ydna_haplogroup.json` — your Y haplogroup with confidence and supporting marker list
- `mtdna_haplogroup.json` — your mtDNA haplogroup with confidence
- `ancient_haplogroup_matches.tsv` — top ancient individuals sharing your lineages, with archaeological context
- `haplogroup_report.md` — a human-readable narrative of the above

The depth of your Y-DNA assignment depends on which markers your array happens to cover; a macro-haplogroup will be assigned with high confidence, but a precise sub-clade requires denser Y sequencing. mtDNA from an array test typically lands at the macro-haplogroup level (R, H, U, J, etc.) and won't resolve sub-haplogroups without dedicated mtDNA sequencing.

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

The AADR `.geno` file is ~5–7 GB (v62: 5.4 GB, v66: 7.2 GB). The pipeline can read it two ways:

- **From R2:** HTTP byte-range requests fetch only the SNP rows we need, in 5,000-row chunks. Runtime ~10–20 minutes depending on connection.
- **From local disk:** if `data/input_data/v*_1240k_public.{geno,ind,anno}` are present, `run_local.py` auto-detects them and skips R2 entirely. Full pipeline (steps 1, 2, 3, 1.5) runs in ~4 minutes end-to-end on a normal laptop.

See [ADR-0010](decisions/0010-chunked-geno-processing.md) and [ADR-0011](decisions/0011-cloudflare-r2-geno-storage.md).

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

Even very common haplogroups can be found at high frequency across populations separated by thousands of miles, because Y and mtDNA lineages have moved with men and women throughout prehistory. A haplogroup is one specific line of descent through your father's-father's-father's line (or mother's-mother's-mother's line); it represents about 1/2^N of your ancestry going back N generations. Your *autosomal* DNA (everything else) is far more representative of "who you are descended from."

**Confidence depends on coverage**

If your matches are based on 400,000 SNPs and the ancient individual has 50,000 SNPs sequenced, only ~50,000 of those overlap and contribute to the comparison. Low-coverage ancients in the top matches list should be taken with a grain of salt. The population-level results aggregate over many individuals and are more reliable.

**The reference is biased**

AADR contains far more European and Near Eastern samples than African or Pacific samples, simply because more research has been done in those regions. If your ancestry is from an under-represented region, the "closest population" results will be less informative.

---

## Step 1.5 — Admixture Decomposition

**Script:** `scripts/step1_5_admixture.py`

### What it does

This is the headline result of the whole pipeline. It decomposes your genome into proportions of ancient source populations — the "you are X% Steppe, Y% Anatolian Farmer, Z% Levantine" answer.

For the modern Western Eurasian model, we fit your allele frequencies as a non-negative weighted combination of six well-defined ancient sources:

| Source | What it is |
|---|---|
| WHG | Western European Hunter-Gatherer (Mesolithic, ~9000–5000 BCE) |
| EHG | Eastern European Hunter-Gatherer (around the Volga–Urals) |
| EEF | Anatolian / Early European Farmer (the wave that brought agriculture into Europe ~6500 BCE) |
| Steppe | Yamnaya / Afanasievo herders from the Pontic-Caspian steppe (~3300 BCE) |
| Levant_N | Levantine Neolithic farmer (Jordan PPNB, Israel Natufian) |
| Iran_N | Iranian Neolithic / Caucasus Hunter-Gatherer |

### The science

After the Last Glacial Maximum, human populations across Western Eurasia descended from a small number of distinct ancestral groups. Almost everyone today is a mixture. Modern Europeans, for example, are roughly:
- Some fraction Mesolithic European Hunter-Gatherer
- Some fraction Anatolian Neolithic Farmer (came in 8,000 years ago)
- Some fraction Steppe / Yamnaya (came in 5,000 years ago)
- Varying smaller fractions of other components

If you can pin down the allele frequencies of each ancestral group from preserved skeletons, you can write your own allele frequencies as a weighted sum of those — and solve for the weights.

### The method (NNLS)

We solve:

```
minimize  ‖A · α − b‖²
subject to   αᵢ ≥ 0,   Σ αᵢ = 1
```

where `A` is the matrix of source allele frequencies, `b` is your allele frequency vector, and `α` is the proportions to estimate. This is *constrained non-negative least squares*, implemented with `scipy.optimize.minimize` using SLSQP.

95 % confidence intervals come from **chromosome-block bootstrap**: resample which chromosomes are in the fit (with replacement) and re-fit 200 times, then take the 2.5 / 97.5 percentiles. Chromosome blocks rather than individual SNPs because adjacent SNPs are in linkage disequilibrium and aren't independent observations.

### The caveats — important ones

**The Levant_N component captures deep Levantine farmer ancestry shared by many populations, not Ashkenazi-specific ancestry.** Anyone with significant Greek, Italian, Spanish, or Near Eastern ancestry will show Levant_N. For specifically detecting Ashkenazi Jewish ancestry you'd need an Iron Age / Bronze Age Levantine source population (which AADR mostly lacks) and a careful three-population test.

**NNLS is the lightweight alternative to qpAdm.** Published ancient-DNA papers use *qpAdm*, which is more rigorous — it uses outgroup f₄-statistics and produces a P-value for model fit. NNLS will get you in the same ballpark, but the exact decimal points and especially the split between correlated sources (like EEF vs Levant_N) will differ. See [ADR-0013](decisions/0013-admixture-nnls.md) for why we chose NNLS over qpAdm.

**Correlated sources can swap weights.** EEF and Levant_N are both Neolithic farmer populations and have similar allele frequencies. NNLS sometimes concentrates farmer-ancestry on one or the other arbitrarily. If your EEF or Levant_N looks oddly extreme, the *sum* of the two is usually well-estimated even when the split isn't.

**CIs reflect statistical noise only.** They tell you how the model would jitter if you re-ran with bootstrapped SNPs. They do not capture model misspecification, source-population miscuration, or the fact that the source pops aren't truly ancestral to you.

### What you get

- `admixture_decomposition.json` — proportions, CIs, residual, source metadata
- `admixture_report.md` — narrative summary with the headline table and caveats
- `source_coverage.tsv` — diagnostic: which AADR individuals were used for each source

For someone with mixed Northern European + Ashkenazi Jewish ancestry, expect roughly:
- Levant_N: 20–35 % (driven by the Levantine farmer component of Ashkenazi ancestry plus deep Mediterranean farmer admixture)
- Steppe: 15–25 %
- WHG: 15–25 %
- EEF: 10–20 %
- EHG: 5–15 %
- Iran_N: 0–10 %

For someone of purely Northern European ancestry, expect Steppe + WHG + EEF to dominate (each ~25–35 %) and Levant_N to be smaller (~5–15 %).

---

## Reading your report

The pipeline emits a single `report.json` (schema_version 1.0) that the web app renders. Here is what each section means, with a worked example from a real 23andMe upload.

**Header — who you are, genetically.**
- **SNPs called** — how many usable genotypes were in your file (e.g. 628,854 for a 23andMe v5 array; ~660k for AncestryDNA).
- **SNP overlap** — how many of those fall on the AADR 1240K panel and survive strand-alignment (e.g. 144,613). This is the working set for every autosomal comparison; a 23andMe array overlaps the 1240K panel less than AncestryDNA, so this number is lower.
- **Y-haplogroup / mt-haplogroup** — your paternal and maternal deep lineages, with a confidence flag. Resolution depends on how many Y/mt markers your chip carried; a consumer array often resolves to a broad level (e.g. `R`, `X`) rather than a deep sub-branch.

**Admixture** — your autosomal ancestry decomposed into six ancient sources (WHG, EHG, EEF, Steppe, Levant_N, Iran_N) via constrained NNLS (ADR 0013), with 95% bootstrap intervals. A typical Northern/Western European looks Steppe + EEF heavy with a WHG minority. **Read the residual**: it's the fraction of your ancestry the six sources *cannot* explain (~0.29 is normal); the higher it is, the more your ancestry sits outside this six-source model, so interpret the proportions loosely.

**Y-DNA / mtDNA TMRCA** — for your closest ancient haplogroup matches, an estimate of how long ago you shared a common paternal/maternal ancestor, from counting branch-defining differences. These are **order-of-magnitude** estimates (see ADR 0014): the raw difference count is the reliable signal; the year figure has wide CIs and is a rough guide, not a date.

**Haplogroup matches** — ancient individuals who share your Y and/or mt haplogroup. The **assessment** flag (PASS / QUESTIONABLE / CRITICAL) is the AADR's own data-quality verdict for *that ancient sample* — derived from contamination estimates (X-chromosome and mtDNA), coverage, ancient-DNA damage, and duplicate/relatedness checks. It says nothing about you or the strength of the match; it tells you how much to trust the ancient individual's data. **PASS** = cleared QC; **QUESTIONABLE** = borderline (weigh it less); **CRITICAL** = a serious quality problem (e.g. contamination) — the shared haplogroup could be an artifact, so treat with caution. (v66 reports these as `Pass`/`PROVISIONAL_PASS`/`Questionable`/etc.; the pipeline normalises them to the three canonical levels.)

**Top population / individual matches** — the ancient populations and individuals closest to you by **allele-sharing distance (ASD)** — the fraction of overlapping SNPs where your genotype differs from theirs, corrected for pseudo-haploidy. Lower is closer. For a European, expect ASD in the ~0.30–0.34 range; the *ranking* is what matters, not the absolute value. Each match carries date, locality, and coordinates (used for the map).

**PCA** — where you fall relative to ancient population clusters. *(Computed in step 3, but not yet wired into `report.json`; the frontend shows a placeholder — a known gap.)*

**Anomalies** — automatic flags for things worth caution (very low overlap with a top match, conflicting haplogroup vs. autosomal signals, etc.).

> One thing to internalise: the closest ancient *individuals* are often not from where your recent ancestors lived. ASD measures overall genetic similarity, and many ancient samples are low-coverage or from admixed periods, so a Northern European can match a Sarmatian or a Viking-era individual simply because those genomes sit near the middle of the relevant genetic cloud. Populations and admixture are the more robust read of "where you come from"; individual matches are more of a curiosity.

## What's coming

All six analysis steps (1.1–1.6) and the web frontend are implemented. Planned next:

- **Wire PCA into the report** — step 3 computes it; `report.json` and the frontend scatter plot need hooking up.
- **More input formats** — FamilyTreeDNA, MyHeritage, Living DNA (AncestryDNA + 23andMe are done).
- **Move compute off the local daemon** — containerise and run on Cloudflare Containers or a scale-to-zero runner, so no local machine is required.

---

## Further reading

- David Reich, *Who We Are and How We Got Here* (2018) — the popular-science book by the lab that produces AADR. Best one-volume introduction to ancient DNA.
- The AADR documentation: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
- ISOGG Y-tree (the canonical Y-DNA reference): https://isogg.org/tree/
- PhyloTree (the canonical mtDNA reference): https://www.phylotree.org/
