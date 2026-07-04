# Kinection

Compare your personal DNA (**AncestryDNA** or **23andMe**, auto-detected) against ~23,000 ancient human genomes from the Allen Ancient DNA Resource (AADR) to produce haplogroup assignments, closest ancient-population and individual matches, admixture proportions, and TMRCA estimates.

Two ways to run it:
- **Local-only** (`scripts/run_local.py`) — your raw DNA never leaves your machine.
- **Web app** — [kinection.pages.dev](https://kinection.pages.dev): the file is uploaded to transient per-job storage, analysed by a local compute daemon, and then **deleted with a verified receipt** (the raw file never traverses the API server). See [Privacy](#privacy).

## Documentation

- **[ARCHITECTURE.md](docs/ARCHITECTURE.md)** — How the system is built. Components, data flow, where things run.
- **[SCIENCE.md](docs/SCIENCE.md)** — What the pipeline actually does, in plain English. The science behind each step.
- **[SECURITY.md](docs/SECURITY.md)** — Security review: what's locked down, what's an acceptable trade-off, what's still open.
- **[docs/decisions/](docs/decisions/)** — Architecture Decision Records (ADRs) explaining individual technical choices.

## Status

- Analysis pipeline (parse → haplogroups → similarity+PCA → TMRCA → admixture → report synthesis): **complete**
- Reference dataset: **AADR v66** (23,250 individuals, 1,233,013 SNPs) in R2
- Input formats: **AncestryDNA + 23andMe** (auto-detected)
- Cloudflare backend (Workers API, R2, D1, KV): **deployed live**
- Web frontend (React SPA on Cloudflare Pages): **deployed** — [kinection.pages.dev](https://kinection.pages.dev)
- Local compute daemon (runs the analysis; fetches the uploaded file from R2): **running**
- **End-to-end validated**: a real upload produces a correct, person-specific report

Known gaps: PCA is computed but not yet wired into `report.json`/the frontend; `update_aadr.py` doesn't yet auto-convert the AADR `tgeno` format on ingest (see [ADR 0017](docs/decisions/0017-aadr-v66-tgeno-and-geno-reader-correction.md)); analysis compute still runs on a local daemon (not hosted).

## Run an analysis locally — no Cloudflare writes

Once your `.env` is set up and the AADR is uploaded to R2 (see Quick Start below), you can analyse an individual entirely on your own machine. The AADR is *read* from R2 (so you don't need 7 GB locally) but every output stays on local disk, and the Worker API is never contacted:

```bash
venv/bin/python scripts/run_local.py --dna /path/to/your-AncestryDNA.txt --label myname
```

This runs all three pipeline steps (parse → haplogroups → similarity+PCA) and writes a combined human-readable report to `output/report_myname.md`. Per-step outputs land in `output/step1_rn/`, `output/step2_rn/`, `output/step3_rn/` — none of which are tracked by git.

## Quick start

This assumes a Cloudflare account with R2 enabled and Wrangler installed (`npm install -g wrangler`).

```bash
# Clone and set up Python env
git clone <repo>
cd dna
python3 -m venv venv
source venv/bin/activate
pip install boto3 requests python-dotenv

# Configure
cp .env.example .env       # then edit .env with your R2 credentials and Worker URL

# One-time: provision Cloudflare resources
wrangler login
wrangler d1 create kinection
wrangler kv namespace create MARKER_CACHE
wrangler r2 bucket create kinection
cd workers/api
wrangler d1 execute kinection --file=src/schema.sql --remote
wrangler secret put COMPUTE_API_KEY --name kinection-api      # paste a random 32-byte hex string
wrangler deploy --name kinection-api
cd ../..

# One-time: upload the AADR reference dataset to R2 (~7 GB)
python scripts/update_aadr.py

# Start the daemon
python scripts/daemon.py
```

## Project layout

```
.
├── scripts/                     Python analysis pipeline
│   ├── step1_parse_harmonise.py   parse + harmonise SNPs (AncestryDNA/23andMe)
│   ├── step2_haplogroup.py        Y-DNA + mtDNA haplogroup assignment
│   ├── step3_similarity_pca.py    genome-wide similarity + PCA
│   ├── step1_4_tmrca.py           Y + mtDNA TMRCA estimates
│   ├── step1_5_admixture.py       6-source admixture decomposition (NNLS)
│   ├── step1_6_synthesis.py       consolidated report.json + map_data.geojson
│   ├── daemon.py                  polls Worker API; fetches upload; runs pipeline
│   ├── convert_tgeno_to_packed.py AADR tgeno → standard packed geno (ADR 0017)
│   ├── update_aadr.py             keeps R2 in sync with Harvard Dataverse
│   ├── utils/                     parsers, R2 client, GENO range reader
│   └── data/                      Y-DNA + mtDNA marker reference DBs
├── frontend/                    React + TS SPA (Cloudflare Pages)
├── workers/api/                 Cloudflare Worker (TypeScript)
├── docs/                        Architecture, science, security, ADRs
├── data/input_data/             (gitignored) your DNA file lives here
├── output/                      (gitignored) analysis results
├── .env.example                 template — copy to .env
└── README.md                    you are here
```

## Privacy

There are two paths, with different handling of your raw file:

- **Local-only run** (`scripts/run_local.py`): the raw DNA file is read from `data/input_data/` on your machine and **never uploaded anywhere**. Only *derived* results leave the machine (if at all).
- **Web upload** (kinection.pages.dev): the browser uploads the raw file **directly to per-job R2 storage via a short-lived presigned URL** — so it **never traverses the Worker/API server**, and no plaintext genotype is ever logged. The local daemon downloads it, runs the analysis, and then **deletes the raw file from R2 and records a verified deletion receipt** in D1. Only derived results (haplogroups, ranked matches, admixture, PCA coordinates) persist, under an unguessable job UUID.

In both paths the API/edge never sees a persisted copy of your raw genotypes. See [`SECURITY.md`](docs/SECURITY.md) for the full model (Step 5.1.1 — the modern-DNA lifecycle).

## License

Personal research project. No license granted for redistribution; ask first.
