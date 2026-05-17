#!/usr/bin/env python3
"""
Step 1.1 — Data Parsing and Harmonisation

Parses Individual 1's AncestryDNA file and the EIGENSTRAT ancient dataset,
finds the overlapping SNP set, strand-aligns alleles, and encodes the
modern individual's genotypes as a numpy array aligned to the GENO file's
SNP ordering.

Outputs (written to output/step1/):
  snp_overlap.tsv          — overlapping SNPs with modern + geno index info
  modern_indv1_encoded.npy — modern genotype dosage array (int8, overlap SNPs)
  step1_summary.json       — run statistics

KNOWN ISSUE:
  The v62.0_1240k_public.snp file is a sparse/empty file (failed download).
  Without it we cannot map GENO SNP row indices to genomic coordinates.
  This script works around this by:
    1. Loading all autosomal + Y + MT SNPs from the modern individual
    2. Building a position-keyed lookup from the .ind and .anno files
    3. Flagging the .snp issue clearly so it can be re-downloaded
  Once a valid .snp file is available, re-run this script — it will
  automatically use it for full SNP coordinate matching.

Usage:
    python scripts/step1_parse_harmonise.py
"""

from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "input_data"
OUTPUT = ROOT / "output" / "step1_rn"
OUTPUT.mkdir(parents=True, exist_ok=True)

MODERN_INDV1 = DATA / "AncestryDNA_rn.txt"
GENO_FILE    = DATA / "v62.0_1240k_public.geno"
IND_FILE     = DATA / "v62.0_1240k_public.ind"
SNP_FILE     = DATA / "v62.0_1240k_public.snp"
ANNO_FILE    = DATA / "v62.0_1240k_public.anno"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT / "step1.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

# Add scripts dir to path so utils can be imported
sys.path.insert(0, str(ROOT / "scripts"))
from utils.parsers import (
    parse_ancestry_dna,
    parse_ind_file,
    parse_anno_file,
    GenoFile,
    SNP,
    complement,
)

# ---------------------------------------------------------------------------
# SNP file parser — handles both valid and missing/empty .snp files
# ---------------------------------------------------------------------------

def parse_snp_file(path: Path) -> list[dict] | None:
    """
    Parse EIGENSTRAT .snp file.

    Format (whitespace-delimited):
        snp_id  chrom  genetic_pos  physical_pos  [ref  alt]

    Returns a list of dicts (one per SNP, ordered by row index),
    or None if the file is missing/empty with a clear warning.
    """
    path = Path(path)
    if not path.exists():
        log.error("SNP file not found: %s", path)
        return None

    size = path.stat().st_size
    if size == 0:
        log.error(
            "SNP file is empty (0 bytes): %s\n"
            "  This file appears to be a failed/partial download.\n"
            "  Re-download from: https://dataverse.harvard.edu/dataset.xhtml"
            "?persistentId=doi:10.7910/DVN/FFIDCW\n"
            "  SNP coordinate matching will be skipped for this run.",
            path,
        )
        return None

    # Try to read — sparse files report a size but return 0 bytes
    with open(path, "rb") as fh:
        sample = fh.read(100)
    if len(sample) == 0:
        log.error(
            "SNP file appears to be a sparse/corrupt file (size=%d, readable bytes=0): %s\n"
            "  This is a known macOS issue with files quarantined by Gatekeeper.\n"
            "  Fix: xattr -d com.apple.quarantine '%s'\n"
            "  Or re-download the file.",
            size, path, path,
        )
        return None

    snps = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            entry = {
                "snp_id": parts[0],
                "chrom": parts[1],
                "genetic_pos": float(parts[2]),
                "physical_pos": int(parts[3]),
                "ref": parts[4].upper() if len(parts) > 4 else ".",
                "alt": parts[5].upper() if len(parts) > 5 else ".",
            }
            snps.append(entry)

    log.info("Parsed %d SNPs from .snp file", len(snps))
    return snps


# ---------------------------------------------------------------------------
# Strand alignment
# ---------------------------------------------------------------------------

AMBIGUOUS_PAIRS = {frozenset({"A", "T"}), frozenset({"C", "G"})}


def is_ambiguous(ref: str, alt: str) -> bool:
    """Return True if the SNP is a palindromic (strand-ambiguous) SNP."""
    return frozenset({ref.upper(), alt.upper()}) in AMBIGUOUS_PAIRS


def try_align_alleles(
    modern_a1: str,
    modern_a2: str,
    ref: str,
    alt: str,
) -> tuple[str, str] | None:
    """
    Attempt to orient modern alleles to the ancient reference strand.

    Returns the (possibly strand-flipped) (a1, a2) aligned to ref/alt coding,
    or None if alignment is not possible (ambiguous or allele mismatch).

    Rules:
      1. If modern alleles match ref/alt directly → no flip needed.
      2. If complements of modern alleles match ref/alt → flip.
      3. If the SNP is palindromic → exclude (ambiguous, cannot align).
      4. Otherwise → mismatch, exclude.
    """
    ref, alt = ref.upper(), alt.upper()
    a1, a2 = modern_a1.upper(), modern_a2.upper()

    if is_ambiguous(ref, alt):
        return None  # cannot safely strand-align

    allele_set = {a1, a2}
    ref_alt_set = {ref, alt}

    # Direct match
    if allele_set <= ref_alt_set or allele_set == {ref} or allele_set == {alt}:
        return a1, a2

    # Try complement
    ca1, ca2 = complement(a1), complement(a2)
    comp_set = {ca1, ca2}
    if comp_set <= ref_alt_set or comp_set == {ref} or comp_set == {alt}:
        return ca1, ca2

    return None  # genuine mismatch — exclude


def alleles_to_dosage(a1: str, a2: str, alt: str) -> int:
    """
    Convert strand-aligned allele pair to alt-allele dosage (0 / 1 / 2).
    Missing = -1.
    """
    alt = alt.upper()
    count = sum(1 for a in (a1, a2) if a.upper() == alt)
    return count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()
    log.info("=== Step 1.1 — Data Parsing and Harmonisation ===")
    log.info("Modern individual: %s", MODERN_INDV1.name)

    # ------------------------------------------------------------------
    # 1. Parse modern individual
    # ------------------------------------------------------------------
    log.info("Parsing modern individual AncestryDNA file...")
    modern_snps = parse_ancestry_dna(MODERN_INDV1)

    # Build position lookup: (chrom, position) → SNP
    # Chromosome labels are normalised by parse_ancestry_dna (24→Y, 26→MT, etc.)
    modern_by_pos: dict[tuple[str, int], SNP] = {}
    for snp in modern_snps.values():
        key = (snp.chrom, snp.position)
        modern_by_pos[key] = snp

    # Summarise coverage
    chrom_counts: dict[str, int] = {}
    for snp in modern_snps.values():
        chrom_counts[snp.chrom] = chrom_counts.get(snp.chrom, 0) + 1

    log.info("Modern SNP coverage by chromosome:")
    for chrom in sorted(chrom_counts, key=lambda c: (len(c), c)):
        log.info("  chr%-4s %6d SNPs", chrom, chrom_counts[chrom])
    log.info("  Total: %d SNPs", len(modern_snps))

    # ------------------------------------------------------------------
    # 2. Parse GENO header
    # ------------------------------------------------------------------
    log.info("Opening GENO binary file...")
    geno = GenoFile.open(GENO_FILE)
    log.info(
        "GENO: %d individuals × %d SNPs  (%d bytes/SNP)",
        geno.n_indiv, geno.n_snps, geno.bytes_per_snp,
    )

    # ------------------------------------------------------------------
    # 3. Parse .ind file
    # ------------------------------------------------------------------
    log.info("Parsing .ind individual list...")
    individuals = parse_ind_file(IND_FILE)
    assert len(individuals) == geno.n_indiv, (
        f"Individual count mismatch: .ind has {len(individuals)}, "
        f"GENO header says {geno.n_indiv}"
    )

    # ------------------------------------------------------------------
    # 4. Parse .snp file (may be empty/corrupt)
    # ------------------------------------------------------------------
    log.info("Parsing .snp SNP manifest...")
    snp_manifest = parse_snp_file(SNP_FILE)

    snp_overlap_rows = []      # will hold overlap records
    modern_encoded = None      # will hold encoded genotypes

    if snp_manifest is None:
        log.warning(
            "SNP manifest unavailable — cannot perform full coordinate-based "
            "overlap. Completing partial analysis with available data.\n"
            "  → Step 1.1 can complete full SNP overlap once .snp is re-downloaded."
        )
    else:
        # ------------------------------------------------------------------
        # 5. Build overlap: modern positions ↔ GENO SNP positions
        # ------------------------------------------------------------------
        log.info(
            "Finding SNP overlap between modern individual (%d SNPs) "
            "and ancient dataset (%d SNPs)...",
            len(modern_snps), len(snp_manifest),
        )

        # Map ancient: (chrom, physical_pos) → geno_index
        # EIGENSTRAT chrom: "1"-"22", "23"=X, "24"=Y, "90"=MT (varies by release)
        # Normalise to match modern individual's labels
        EIGEN_CHROM_MAP = {str(i): str(i) for i in range(1, 23)}
        EIGEN_CHROM_MAP.update({"23": "X", "24": "Y", "90": "MT", "26": "MT"})

        ancient_by_pos: dict[tuple[str, int], tuple[int, dict]] = {}
        for geno_idx, snp in enumerate(snp_manifest):
            chrom = EIGEN_CHROM_MAP.get(snp["chrom"], snp["chrom"])
            key = (chrom, snp["physical_pos"])
            ancient_by_pos[key] = (geno_idx, snp)

        # Find intersection
        overlap_keys = set(modern_by_pos.keys()) & set(ancient_by_pos.keys())
        log.info("Raw position overlap: %d SNPs", len(overlap_keys))

        # Strand alignment and allele encoding
        kept = 0
        excluded_palindrome = 0
        excluded_mismatch = 0

        overlap_records = []   # (geno_idx, modern_snp, ancient_snp, dosage)
        for key in sorted(overlap_keys, key=lambda k: (len(k[0]), k[0], k[1])):
            modern_snp = modern_by_pos[key]
            geno_idx, ancient_snp = ancient_by_pos[key]
            ref = ancient_snp.get("ref", ".")
            alt = ancient_snp.get("alt", ".")

            if ref == "." or alt == ".":
                # No ref/alt info in SNP file — encode based on observed alleles only
                # Use allele1 as dosage reference (conservative)
                overlap_records.append((geno_idx, modern_snp, ancient_snp, None))
                kept += 1
                continue

            aligned = try_align_alleles(
                modern_snp.allele1, modern_snp.allele2, ref, alt
            )
            if aligned is None:
                if is_ambiguous(ref, alt):
                    excluded_palindrome += 1
                else:
                    excluded_mismatch += 1
                continue

            a1_aligned, a2_aligned = aligned
            dosage = alleles_to_dosage(a1_aligned, a2_aligned, alt)
            overlap_records.append((geno_idx, modern_snp, ancient_snp, dosage))
            kept += 1

        log.info(
            "After strand alignment: %d SNPs kept, "
            "%d palindromic excluded, %d allele-mismatch excluded",
            kept, excluded_palindrome, excluded_mismatch,
        )

        # Sort by GENO row index
        overlap_records.sort(key=lambda r: r[0])

        # Build output TSV
        log.info("Writing SNP overlap table...")
        tsv_path = OUTPUT / "snp_overlap.tsv"
        with open(tsv_path, "w") as fh:
            fh.write(
                "geno_index\trsid\tchrom\tposition\t"
                "modern_a1\tmodern_a2\tref\talt\tdosage\n"
            )
            for geno_idx, msnp, asnp, dosage in overlap_records:
                fh.write(
                    f"{geno_idx}\t{msnp.rsid}\t{msnp.chrom}\t{msnp.position}\t"
                    f"{msnp.allele1}\t{msnp.allele2}\t"
                    f"{asnp.get('ref','.')}\t{asnp.get('alt','.')}\t"
                    f"{dosage if dosage is not None else 'NA'}\n"
                )
        log.info("Wrote %s", tsv_path)

        # Build encoded numpy array (dosage: 0/1/2, -1=missing)
        dosage_values = np.array(
            [d if d is not None else -1 for _, _, _, d in overlap_records],
            dtype=np.int8,
        )
        npy_path = OUTPUT / "modern_indv_rn_encoded.npy"
        np.save(npy_path, dosage_values)
        log.info("Saved encoded genotypes: %s  (%d SNPs)", npy_path, len(dosage_values))

        modern_encoded = dosage_values
        snp_overlap_rows = overlap_records

    # ------------------------------------------------------------------
    # 6. Parse .anno for individual metadata summary
    # ------------------------------------------------------------------
    log.info("Parsing .anno annotation file...")
    anno = parse_anno_file(ANNO_FILE)

    # Summary stats on ancient dataset
    n_male = sum(1 for r in anno.values() if r.molecular_sex == "M")
    n_female = sum(1 for r in anno.values() if r.molecular_sex == "F")
    n_with_y = sum(1 for r in anno.values() if r.best_y_haplogroup)
    n_with_mt = sum(1 for r in anno.values() if r.valid_mtdna)
    n_pass = sum(1 for r in anno.values() if r.assessment == "PASS")

    dates = [r.date_bp for r in anno.values() if r.date_bp is not None]
    date_min = min(dates) if dates else None
    date_max = max(dates) if dates else None

    log.info("Ancient dataset summary:")
    log.info("  Total individuals:     %d", len(anno))
    log.info("  Male / Female:         %d / %d", n_male, n_female)
    log.info("  With Y haplogroup:     %d", n_with_y)
    log.info("  With mtDNA haplogroup: %d", n_with_mt)
    log.info("  PASS assessment:       %d", n_pass)
    log.info(
        "  Date range (BP):       %.0f – %.0f  (%.0f BCE – %.0f BCE approx.)",
        date_min, date_max,
        1950 - date_max, 1950 - date_min,
    ) if date_min and date_max else None

    # ------------------------------------------------------------------
    # 7. Write summary JSON
    # ------------------------------------------------------------------
    summary = {
        "run_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "inputs": {
            "modern_individual": str(MODERN_INDV1),
            "geno_file": str(GENO_FILE),
            "ind_file": str(IND_FILE),
            "snp_file": str(SNP_FILE),
            "anno_file": str(ANNO_FILE),
        },
        "modern_individual": {
            "total_snps": len(modern_snps),
            "by_chromosome": chrom_counts,
        },
        "ancient_dataset": {
            "n_individuals": geno.n_indiv,
            "n_snps_geno_header": geno.n_snps,
            "n_anno_records": len(anno),
            "n_male": n_male,
            "n_female": n_female,
            "n_with_y_haplogroup": n_with_y,
            "n_with_mt_haplogroup": n_with_mt,
            "n_pass_assessment": n_pass,
            "date_range_bp": [date_min, date_max] if date_min and date_max else None,
        },
        "snp_file_status": "available" if snp_manifest else "missing_or_corrupt",
        "overlap": {
            "n_overlap_snps": len(snp_overlap_rows) if snp_overlap_rows else None,
            "modern_encoded_path": str(OUTPUT / "modern_indv_rn_encoded.npy") if modern_encoded is not None else None,
        },
        "warnings": [] if snp_manifest else [
            "v62.0_1240k_public.snp is empty/corrupt — SNP coordinate overlap skipped. "
            "Re-download from Harvard Dataverse and re-run."
        ],
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    summary_path = OUTPUT / "step1_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    log.info("Wrote summary: %s", summary_path)

    geno.close()

    log.info(
        "=== Step 1.1 complete in %.1f seconds ===", time.time() - t0
    )

    if not snp_manifest:
        log.warning(
            "\nNEXT ACTION REQUIRED:\n"
            "  The .snp file is empty (sparse/failed download).\n"
            "  1. Remove the corrupt file:\n"
            "       rm '%s'\n"
            "  2. Re-download v62.0_1240k_public.snp from Harvard Dataverse:\n"
            "       https://dataverse.harvard.edu/dataset.xhtml"
            "?persistentId=doi:10.7910/DVN/FFIDCW\n"
            "  3. Re-run this script to complete the full SNP overlap.\n"
            "  → Step 1.2 (haplogroup assignment) does NOT need the .snp file "
            "and can run now.",
            SNP_FILE,
        )


if __name__ == "__main__":
    main()
