# Kinection

Compare your personal AncestryDNA results against ~19,000 ancient human genomes from the Allen Ancient DNA Resource (AADR), running locally so that your raw DNA never leaves your machine.

## Documentation

- **[ARCHITECTURE.md](docs/ARCHITECTURE.md)** — How the system is built. Components, data flow, where things run.
- **[SCIENCE.md](docs/SCIENCE.md)** — What the pipeline actually does, in plain English. The science behind each step.
- **[SECURITY.md](docs/SECURITY.md)** — Security review: what's locked down, what's an acceptable trade-off, what's still open.
- **[docs/decisions/](docs/decisions/)** — Architecture Decision Records (ADRs) explaining individual technical choices.

## Status

- Phase 1 analysis pipeline (parse, haplogroup, similarity+PCA): **complete**
- Cloudflare backend (Workers API, R2, D1, KV): **deployed**
- Local compute daemon: **ready**
- Web frontend: **not started** (Phase 2)

For a full status snapshot see [`PROJECT_UPDATE_2026-05-17.txt`](PROJECT_UPDATE_2026-05-17.txt).

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
│   ├── step1_parse_harmonise.py   parse + harmonise SNPs
│   ├── step2_haplogroup.py        Y-DNA + mtDNA haplogroup assignment
│   ├── step3_similarity_pca.py    genome-wide similarity + PCA
│   ├── daemon.py                  polls Worker API for jobs
│   ├── update_aadr.py             keeps R2 in sync with Harvard Dataverse
│   ├── utils/                     parsers, R2 client, GENO range reader
│   └── data/                      Y-DNA + mtDNA marker reference DBs
├── workers/api/                 Cloudflare Worker (TypeScript)
├── docs/                        Architecture, science, security, ADRs
├── data/input_data/             (gitignored) your DNA file lives here
├── output/                      (gitignored) analysis results
├── .env.example                 template — copy to .env
└── README.md                    you are here
```

## Privacy

Your raw DNA file is **never** uploaded to any server. The local compute daemon reads it from `data/input_data/` on your machine, runs the analysis, and uploads only the *derived* results (haplogroup reports, ranked matches, PCA coordinates) to R2 under an unguessable job UUID. The cloud side of the system has no way to ever see the file.

See [`SECURITY.md`](docs/SECURITY.md) for the full security model.

## License

Personal research project. No license granted for redistribution; ask first.
