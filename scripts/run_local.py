"""
Run the Kinection analysis pipeline locally for a single individual.

Reads AADR reference data from R2 (so you don't need 7 GB locally) but
keeps every analysis output on local disk. Nothing is written to
Cloudflare R2 or D1, and the Worker API is not contacted.

Usage:
  python scripts/run_local.py                           # uses default DNA path
  python scripts/run_local.py --dna PATH                # specify DNA file
  python scripts/run_local.py --dna PATH --label NAME   # name the report

Requirements:
  - .env populated with R2 credentials (R2 is used as a read-only source for AADR)
  - AADR uploaded to R2 (run scripts/update_aadr.py first if not already done)

Output:
  output/step1_rn/  — step 1 outputs (overlap, encoded genotypes, summary)
  output/step2_rn/  — step 2 outputs (haplogroup assignments, matches)
  output/step3_rn/  — step 3 outputs (distances, PCA, population matches)
  output/report_<label>.md  — combined human-readable report
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
OUTPUT = ROOT / "output"


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------

def run_step(script: str, env: dict) -> None:
    """Run one pipeline step as a subprocess. Raises if it fails."""
    print(f"\n{'=' * 70}")
    print(f"  Running {script}")
    print(f"{'=' * 70}")
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / script)],
        env=env,
        check=False,
    )
    if result.returncode != 0:
        sys.exit(f"\n{script} exited with code {result.returncode}. Aborting.")


# ---------------------------------------------------------------------------
# Report synthesis
# ---------------------------------------------------------------------------

def read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def read_text(path: Path) -> str | None:
    try:
        return path.read_text()
    except Exception:
        return None


def build_report(label: str, run_started: datetime) -> str:
    """Combine step 1/2/3 outputs into a single Markdown report."""
    step1_summary = read_json(OUTPUT / "step1_rn" / "step1_summary.json")
    step2_y       = read_json(OUTPUT / "step2_rn" / "ydna_haplogroup.json")
    step2_mt      = read_json(OUTPUT / "step2_rn" / "mtdna_haplogroup.json")
    step2_report  = read_text(OUTPUT / "step2_rn" / "haplogroup_report.md")
    step3_report  = read_text(OUTPUT / "step3_rn" / "top_matches_report.md")
    pca_variance  = read_json(OUTPUT / "step3_rn" / "pca_variance_explained.json")
    admixture     = read_json(OUTPUT / "step1_5_rn" / "admixture_decomposition.json")
    admix_report  = read_text(OUTPUT / "step1_5_rn" / "admixture_report.md")

    lines: list[str] = []
    lines.append(f"# Kinection Analysis Report — {label}")
    lines.append("")
    lines.append(f"*Generated: {run_started:%Y-%m-%d %H:%M:%S}*")
    lines.append("")
    lines.append("This report was produced locally. None of the data below was uploaded to any cloud service.")
    lines.append("")
    lines.append("---")
    lines.append("")

    # ── Headline summary ──────────────────────────────────────────────
    lines.append("## Headline Results")
    lines.append("")
    if step2_y or step2_mt:
        lines.append("| Lineage | Haplogroup | Confidence |")
        lines.append("|---|---|---|")
        if step2_y:
            lines.append(
                f"| **Paternal (Y-DNA)** | {step2_y.get('haplogroup', 'n/a')} | "
                f"{step2_y.get('confidence', 'n/a')} |"
            )
        if step2_mt:
            lines.append(
                f"| **Maternal (mtDNA)** | {step2_mt.get('haplogroup', 'n/a')} | "
                f"{step2_mt.get('confidence', 'n/a')} |"
            )
        lines.append("")

    if step1_summary:
        lines.append("**Data overlap:**")
        modern = step1_summary.get("modern_individual", {})
        ancient = step1_summary.get("ancient_dataset", {})
        overlap = step1_summary.get("overlap", {})
        if "total_snps" in modern:
            lines.append(f"- Modern individual SNPs: **{modern['total_snps']:,}**")
        if "n_snps_geno_header" in ancient:
            lines.append(f"- Ancient dataset SNPs: **{ancient['n_snps_geno_header']:,}** "
                         f"across **{ancient.get('n_individuals', 0):,}** individuals")
        if "n_overlap_snps" in overlap:
            lines.append(f"- Overlapping SNPs (post-palindromic-filter): **{overlap['n_overlap_snps']:,}**")
        lines.append("")

    # ── Admixture decomposition headline ──────────────────────────────
    if admixture:
        lines.append("**Ancient ancestry composition:**")
        lines.append("")
        lines.append("| Source | % | 95% CI |")
        lines.append("|---|---:|---:|")
        props = admixture.get("proportions", {})
        cis = admixture.get("ci95", {})
        for name in sorted(props, key=lambda n: -props[n]):
            p = props[name] * 100
            lo, hi = (c * 100 for c in cis.get(name, (0, 0)))
            lines.append(f"| {name} | {p:.1f}% | {lo:.1f}–{hi:.1f}% |")
        lines.append("")

    lines.append("---")
    lines.append("")

    # ── Step 1 ────────────────────────────────────────────────────────
    lines.append("## Step 1 — Data Parsing and Harmonisation")
    lines.append("")
    if step1_summary:
        lines.append("Full statistics:")
        lines.append("```json")
        lines.append(json.dumps(step1_summary, indent=2))
        lines.append("```")
    else:
        lines.append("_Summary not available._")
    lines.append("")
    lines.append("---")
    lines.append("")

    # ── Step 2 ────────────────────────────────────────────────────────
    lines.append("## Step 2 — Haplogroup Assignment")
    lines.append("")
    if step2_report:
        # Drop the first H1 from the embedded report so it nests cleanly
        body = "\n".join(
            line for line in step2_report.splitlines()
            if not line.startswith("# ")
        )
        lines.append(body.strip())
    else:
        lines.append("_Step 2 report not available._")
    lines.append("")
    lines.append("---")
    lines.append("")

    # ── Step 3 ────────────────────────────────────────────────────────
    lines.append("## Step 3 — Genome-wide Similarity and PCA")
    lines.append("")
    if step3_report:
        body = "\n".join(
            line for line in step3_report.splitlines()
            if not line.startswith("# ")
        )
        lines.append(body.strip())
    else:
        lines.append("_Step 3 report not available._")
    lines.append("")

    if pca_variance:
        lines.append("**PCA variance explained:**")
        evals = pca_variance.get("eigenvalues", [])
        for i, ev in enumerate(evals[:10], start=1):
            lines.append(f"- PC{i}: {ev:.2%}")
        lines.append("")

    lines.append("---")
    lines.append("")

    # ── Step 1.5 ──────────────────────────────────────────────────────
    lines.append("## Step 1.5 — Admixture Decomposition")
    lines.append("")
    if admix_report:
        body = "\n".join(
            line for line in admix_report.splitlines()
            if not line.startswith("# ")
        )
        lines.append(body.strip())
    else:
        lines.append("_Step 1.5 report not available._")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## How to read this")
    lines.append("")
    lines.append("See [`docs/SCIENCE.md`](../docs/SCIENCE.md) for plain-English explanations "
                 "of haplogroups, allele-sharing distance, and PCA — and the caveats that come with them.")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the Kinection analysis pipeline locally for one individual.",
    )
    parser.add_argument(
        "--dna",
        help="Path to the AncestryDNA raw file. "
             "Defaults to data/input_data/AncestryDNA_rn.txt.",
    )
    parser.add_argument(
        "--label",
        default="run",
        help="Label appended to the report filename (default: 'run').",
    )
    args = parser.parse_args()

    # Verify .env has R2 credentials before starting a 10+ minute pipeline
    missing = [v for v in (
        "R2_ENDPOINT_URL", "R2_ACCESS_KEY_ID", "R2_SECRET_ACCESS_KEY"
    ) if not os.environ.get(v)]
    if missing:
        # python-dotenv autoloads .env when imported by r2_client; mimic that here
        try:
            from dotenv import load_dotenv
            load_dotenv(ROOT / ".env")
            missing = [v for v in missing if not os.environ.get(v)]
        except ImportError:
            pass
    if missing:
        sys.exit(f"Missing R2 credentials in environment: {missing}. "
                 f"Populate .env from .env.example.")

    run_started = datetime.now()
    job_label = f"local-{run_started:%Y%m%d-%H%M%S}"

    # Auto-detect whether AADR is available on local disk. If so, use local
    # files (faster, fully offline, no R2 access at all). Otherwise read
    # from R2 via HTTP range requests.
    aadr_local = ROOT / "data" / "input_data"
    have_local_aadr = all(
        any(aadr_local.glob(f"v*_1240k_public.{ext}")) or
        any(aadr_local.glob(f"v*.1240K.aadr.PUB.{ext}"))
        for ext in ("geno", "ind", "anno")
    )

    env = {
        **os.environ,
        "LOCAL_OUTPUTS": "1",
        "JOB_ID": job_label,
        "USE_R2": "0" if have_local_aadr else "1",
    }
    if args.dna:
        env["MODERN_DNA"] = str(Path(args.dna).expanduser().resolve())

    print(f"Job label   : {job_label}")
    print(f"DNA file    : {env.get('MODERN_DNA', '(default: data/input_data/AncestryDNA_rn.txt)')}")
    print(f"AADR source : {'local disk (no R2 access)' if have_local_aadr else 'R2 (read-only)'}")
    print(f"Outputs     : local only (no Cloudflare writes)")

    t0 = time.time()
    run_step("step1_parse_harmonise.py", env)
    run_step("step2_haplogroup.py", env)
    run_step("step3_similarity_pca.py", env)
    run_step("step1_5_admixture.py", env)
    elapsed = time.time() - t0

    # Build combined report
    report_path = OUTPUT / f"report_{args.label}.md"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report_path.write_text(build_report(args.label, run_started))

    print(f"\n{'=' * 70}")
    print(f"  Pipeline complete in {elapsed:.0f}s")
    print(f"  Combined report: {report_path}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
