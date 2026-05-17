# ADR-0004: Use Python for the Analysis Engine

* Status: Accepted
* Date: 2026-04-21

## Context and Problem Statement

Population genetics analysis has historically been dominated by R (via packages like adegenet, hierfstat, PLINK wrappers) and dedicated command-line tools (EIGENSOFT, PLINK, ADMIXTOOLS). The analysis pipeline here involves parsing custom binary formats, large matrix operations, and will eventually feed results directly into a web application backend. Which language should be used for the core analysis engine?

## Decision Drivers

* Reusable codebase between analysis engine and web backend (reduce total language footprint)
* Strong numerical computing support for large matrix operations (394k SNPs × 17k individuals)
* Ability to write a custom binary parser for the EIGENSTRAT .geno format
* Active ecosystem for the specific libraries needed (numpy, scikit-learn)
* Developer familiarity and long-term maintainability

## Considered Options

* Python (NumPy, scikit-learn, custom parsers)
* R (adegenet, PopGenome, hierfstat, vcfR)
* Dedicated tools: EIGENSOFT + PLINK + ADMIXTOOLS called as subprocesses
* Julia (high-performance numerical computing)

## Decision Outcome

Chosen option: **Python**, because it enables a single language across the analysis engine and the planned FastAPI web backend, has mature numerical computing (NumPy, scikit-learn) that can handle the scale requirements, and allows direct low-level binary parsing of EIGENSTRAT .geno files without external tool dependencies.

### Positive Consequences

* Single language across analysis pipeline and web backend — no context switching, shared utilities
* NumPy vectorised operations handle the 394k × 17k matrix comparisons efficiently in pure Python
* Custom GenoFile class reads EIGENSTRAT binary format directly — no EIGENSOFT installation required
* scikit-learn PCA is well-tested and handles projection of new samples naturally
* Easier to package for production deployment than compiled C tools

### Negative Consequences

* R has a richer ecosystem for population genetics-specific analyses (FST, AMOVA, qpAdm wrappers)
  — Steps 1.5 (AMOVA) may require either implementing AMOVA from scratch or calling R via subprocess
* Some advanced analyses (qpAdm, D-statistics) are only natively available in ADMIXTOOLS (C++) or
  its R wrapper (admixtools R package) — will require subprocess calls or reimplementation
* Python's GIL limits true parallelism for CPU-bound tasks; chunked sequential processing is used
  as a workaround

## Pros and Cons of the Options

### Python (chosen)

* Good, because single language for analysis and web backend
* Good, because NumPy + scikit-learn cover ASD, PCA, and distance ranking natively
* Good, because custom binary parsers straightforward to implement
* Good, because excellent package management (pip/venv) for reproducible environments
* Bad, because less idiomatic for population genetics than R
* Bad, because some advanced stats (qpAdm, formal AMOVA) not available as Python packages

### R

* Good, because dominant language in academic population genetics — extensive package ecosystem
* Good, because admixtools, hierfstat, PopGenome cover most planned analyses
* Bad, because separate language from the planned FastAPI web backend — two codebases to maintain
* Bad, because packaging R for production web deployment is significantly harder
* Bad, because R's memory model less efficient for very large matrix operations

### Dedicated tools (EIGENSOFT / PLINK / ADMIXTOOLS)

* Good, because gold standard implementations used in all academic papers
* Good, because highly optimised C/C++ — fast for large datasets
* Bad, because require installation and PATH management on the production server
* Bad, because subprocess integration adds fragility and limits error handling
* Bad, because inputs/outputs are file-based — adds I/O overhead per user analysis
* Bad, because no shared code with the web backend

### Julia

* Good, because high-performance numerical computing, close to C speed
* Bad, because small ecosystem for population genetics specifically
* Bad, because high barrier to entry for web deployment
* Bad, because no shared ecosystem with the web backend
